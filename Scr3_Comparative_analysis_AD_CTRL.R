# ==============================================================================
# Pipeline: CHIP Mutation Burden Analysis (AD vs. Control) and Publication Figures
# Description: Takes the merged SNV + indel calls from Scr2, attaches sample-level
#              clinical metadata, builds an age/sex/CVD-matched CTRL vs. AD subset
#              (MatchIt), and compares somatic / CHIP mutation burden between groups:
#                - Coverage QC across all samples                      ED Fig 2
#                - Age distribution before vs. after matching          Fig 1A
#                - Age, sex and APOE composition by group              ED Fig 1
#                - Mutation rate per bp by variant class               Fig 1B, 1C
#                - Mutation landscape heatmap (matched CHIP / all)     Fig 1D, ED Fig 4
#                - Logistic regression of AD status on CHIP burden     Fig 1E
#                - Mutation rate stratified by APOE genotype           Fig 1F
#                - Odds ratio as a function of minimum VAF             Fig 1G
#                - CHIP rate per VAF bin (+ meta-analysis, e3/e3)      Fig 1H
#
# Inputs:
#   AD_panel.merged.addage.tsv          Merged SNV + indel calls (Scr2 output; must
#                                       include the Location and Type columns)
#   SupplTable1_Sample_Information.xlsx Sample-level clinical metadata for the 353
#                                       QC-passed samples (sheet "Sample_Information")
#   AD_panel.coverage.summary           Per-sample mean depth and fraction of
#                                       bases >=200X (Scr1 output)
#   3277911_Covered.bed                 Capture panel design (2 header lines)
#   whitelist_Genename.tsv              CHIP whitelist gene names (Scr2 output)
#
# Group definition (column "Group" of SupplTable1):
#   AD   : clinically diagnosed AD
#   C    : non-AD samples older than 90 y (labelled "Centenarian" in figures)
#   CTRL : all other non-AD samples
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
  library(tidyr)
  library(readxl)
  library(MatchIt)
  library(meta)
})

# ==============================================================================
# Configuration
# ==============================================================================
merged_file       <- "AD_panel.merged.addage.tsv"
sample_info_file  <- "SupplTable1_Sample_Information.xlsx"
sample_info_sheet <- "Sample_Information"
coverage_file     <- "AD_panel.coverage.summary"
panel_bed_file    <- "3277911_Covered.bed"
chip_gene_file    <- "whitelist_Genename.tsv"
output_header     <- "AD_panel.merged.addage.pub.clinFilt"

MIN_AGE        <- 60                    # minimum age for the matched analysis
MATCH_RATIO    <- 2                     # controls matched per AD case
MATCH_CALIPER  <- 0.1                   # propensity-score caliper (SD units)
MIN_VAF_STEPS  <- seq(0, 0.1, by = 0.005)
VAF_BIN_BREAKS <- c(-Inf, 0.02, 0.05, 0.1, Inf)
VAF_BIN_LABELS <- c("VAF<=0.02", "0.02<VAF<=0.05", "0.05<VAF<=0.1", "0.1<VAF")

GROUP_LEVELS <- c("CTRL", "C", "AD")    # all samples
MATCH_LEVELS <- c("CTRL", "AD")         # matched subset (C excluded)
GROUP_LABELS <- c(CTRL = "Control", C = "Centenarian", AD = "AD")
GROUP_COLORS <- c(CTRL = "royalblue1", C = "darkblue", AD = "indianred3")
GROUP_COLORS_DARK <- c(CTRL = "royalblue4", C = "black", AD = "indianred4")
# Fill order of interaction(Group, Class): CTRL.indel, AD.indel, CTRL.snp, AD.snp
CLASS_FILL <- c("royalblue4", "indianred4", "royalblue1", "indianred3")

# Mutation-type labels for the heatmaps (Fig 1D, ED Fig 4)
TYPE_LABELS <- c("nonsynonymous"           = "Missense mutation",
                 "stopgain"                = "Nonsense mutation",
                 "nonframeshift deletion"  = "In-frame deletion",
                 "nonframeshift insertion" = "In-frame insertion",
                 "frameshift deletion"     = "Frameshift deletion",
                 "frameshift insertion"    = "Frameshift insertion")
MUT_TYPE_LEVELS <- c(TYPE_LABELS, "Splice-site mutation", "Multiple hits")
# Tile colours per mutation type, as in the published figures. The CHIP heatmap uses
# its own palette (the whitelist yields no in-frame indels; colours are kept for them
# in case any appear).
HEATMAP_COLORS_ALL  <- c("Missense mutation"    = "gold1",      "Nonsense mutation"    = "red2",
                         "In-frame deletion"    = "forestgreen", "In-frame insertion"   = "royalblue1",
                         "Frameshift deletion"  = "purple",     "Frameshift insertion" = "darkorange1",
                         "Splice-site mutation" = "black",      "Multiple hits"        = "white")
HEATMAP_COLORS_CHIP <- c("Missense mutation"    = "gold1",      "Nonsense mutation"    = "red2",
                         "In-frame deletion"    = "forestgreen", "In-frame insertion"   = "turquoise3",
                         "Frameshift deletion"  = "royalblue1", "Frameshift insertion" = "purple",
                         "Splice-site mutation" = "darkorange1", "Multiple hits"       = "black")

# Axis labels for logistic-regression coefficients (Fig 1E)
TERM_LABELS <- c("scale(rate_chip)"     = "CHIP burden",
                 "factor(CVD.history)1" = "CVD history",
                 "SexM"                 = "Sex: Male",
                 "scale(Age)"           = "Age",
                 "factor(ApoE)34"       = "APOE: e3/e4")

# ==============================================================================
# Helper Functions
# ==============================================================================

# Load sample-level clinical metadata (SupplTable1) and rename columns for analysis
load_sample_info <- function(path, sheet) {
  raw <- data.frame(read_excel(path, sheet = sheet), check.names = FALSE)
  info <- data.frame(
    ID                 = raw[["ID"]],
    ApoE               = as.character(raw[["APOE genotype"]]),        # "33" (e3/e3) or "34" (e3/e4)
    Age                = raw[["Age (at blood collection)"]],
    Age_at_diagnosis   = raw[["Age (at diagnosis)"]],                 # AD cases only
    Sex                = raw[["Sex"]],                                # "F" / "M"
    MMSE               = raw[["MMSE"]],
    Group              = factor(raw[["Group"]], levels = GROUP_LEVELS),
    Batch              = raw[["Batch"]],
    Dementia.type      = raw[["Dementia type"]],
    CVD.history        = raw[["CVD history"]],                        # 0 / 1
    Selected_published = as.logical(raw[["Selected for subset"]]),
    stringsAsFactors   = FALSE
  )
  info[order(info$ID), ]   # fixed row order keeps nearest-neighbour matching reproducible
}

# Named per-group counts / summed coverage (bp covered >=200X)
count_by_group <- function(df, lv = MATCH_LEVELS) {
  setNames(as.vector(table(factor(df$Group, levels = lv))), lv)
}
coverage_by_group <- function(samples, lv = MATCH_LEVELS) {
  cov <- merge(samples[, c("ID", "Group")], id_coverage[, c("ID", "Coverage")], by = "ID")
  setNames(as.vector(tapply(cov$Coverage, factor(cov$Group, levels = lv), sum)), lv)
}

# Add rate and 95% binomial CI columns: rate = rate_count / n, CI from ci_count / n
add_rate <- function(tab, prefix, rate_count, ci_count, n) {
  ci <- t(mapply(function(k, m) suppressWarnings(prop.test(k, m))$conf.int, tab[[ci_count]], tab[[n]]))
  tab[[prefix]] <- tab[[rate_count]] / tab[[n]]
  tab[[paste0(prefix, "_LB")]] <- ci[, 1]
  tab[[paste0(prefix, "_UB")]] <- ci[, 2]
  tab
}

# One-sided test (AD > CTRL) of mutation rate, plus the AD/CTRL rate ratio
test_ad_vs_ctrl <- function(tab, count_col, n_col, label) {
  ad   <- tab[tab$Group == "AD", ][1, ]
  ctrl <- tab[tab$Group == "CTRL", ][1, ]
  cat(sprintf("\n--- %s ---\n", label))
  print(prop.test(c(ad[[count_col]], ctrl[[count_col]]), c(ad[[n_col]], ctrl[[n_col]]), alternative = "greater"))
  cat(sprintf("Rate ratio (AD/CTRL): %.3f\n",
              (ad[[count_col]] / ad[[n_col]]) / (ctrl[[count_col]] / ctrl[[n_col]])))
}

# Mutation rates per group split by variant class (snp / indel).
#   Rate1: all somatic mutations per bp of panel
#   Rate3: CHIP (whitelist) mutations per bp of CHIP genes
# Bar heights use class-specific counts (stacked); CIs use the class-combined total.
build_class_rate_table <- function(v, s) {
  base <- data.frame(Group    = factor(MATCH_LEVELS, levels = MATCH_LEVELS),
                     Count1   = count_by_group(v),
                     Count3   = count_by_group(v[which(v$CHIP), ]),
                     Coverage = coverage_by_group(s))
  tab <- do.call(rbind, lapply(c("snp", "indel"), function(cl) {
    vc <- v[v$Class == cl, ]
    cbind(base, Class = cl,
          Count_sep1 = count_by_group(vc),
          Count_sep3 = count_by_group(vc[which(vc$CHIP), ]))
  }))
  tab$Class     <- factor(tab$Class, levels = c("indel", "snp"))
  tab$Coverage3 <- tab$Coverage * chip_fraction
  tab <- add_rate(tab, "Rate1", "Count_sep1", "Count1", "Coverage")
  tab <- add_rate(tab, "Rate3", "Count_sep3", "Count3", "Coverage3")
  rownames(tab) <- NULL
  tab
}

# CHIP mutation rate per group within each VAF bin
build_vaf_rate_table <- function(v, s) {
  v$VAF_bin <- cut(v$MAF, VAF_BIN_BREAKS, labels = VAF_BIN_LABELS, right = TRUE)
  cov <- coverage_by_group(s)
  tab <- do.call(rbind, lapply(VAF_BIN_LABELS, function(b) {
    vb <- v[which(v$VAF_bin == b), ]
    data.frame(Group    = factor(MATCH_LEVELS, levels = MATCH_LEVELS),
               Type     = b,
               Count3   = count_by_group(vb[which(vb$CHIP), ]),
               Coverage = cov)
  }))
  tab$Type      <- factor(tab$Type, levels = VAF_BIN_LABELS)
  tab$Coverage3 <- tab$Coverage * chip_fraction
  tab <- add_rate(tab, "Rate3", "Count3", "Count3", "Coverage3")
  rownames(tab) <- NULL
  tab
}

# Fisher odds ratio (AD vs CTRL, mutations vs covered bp) for mutations with VAF >= threshold
#   OR1: all somatic mutations; OR2: CHIP mutations (CHIP-gene coverage)
build_or_table <- function(v, s) {
  cov      <- coverage_by_group(s)
  cov_chip <- round(cov * chip_fraction)   # fisher.test requires integer counts
  do.call(rbind, lapply(MIN_VAF_STEPS, function(t) {
    vt     <- v[which(v$MAF >= t), ]
    n_all  <- count_by_group(vt)
    n_chip <- count_by_group(vt[which(vt$CHIP), ])
    f1 <- fisher.test(matrix(c(n_all["AD"],  cov["AD"],      n_all["CTRL"],  cov["CTRL"]),      ncol = 2))
    f2 <- fisher.test(matrix(c(n_chip["AD"], cov_chip["AD"], n_chip["CTRL"], cov_chip["CTRL"]), ncol = 2))
    data.frame(MAF = t,
               OR1 = f1$estimate, OR1_LB = f1$conf.int[1], OR1_UB = f1$conf.int[2],
               OR2 = f2$estimate, OR2_LB = f2$conf.int[1], OR2_UB = f2$conf.int[2],
               row.names = NULL)
  }))
}

# Prepare one row per sample x gene for the mutation heatmaps:
#   - recode ANNOVAR consequences to display labels (splicing overrides the exonic type)
#   - drop synonymous / unannotated calls
#   - order genes by number of calls (most frequently mutated gene at the top)
#   - collapse samples with >1 call in the same gene into "Multiple hits"
prepare_heatmap_calls <- function(v) {
  type   <- v$Type
  mapped <- which(type %in% names(TYPE_LABELS))
  type[mapped] <- TYPE_LABELS[type[mapped]]
  type[which(v$Location == "splicing")] <- "Splice-site mutation"
  v$Type <- type
  v <- v[!is.na(v$Type) & v$Type != "synonymous", ]

  v$Symbol <- factor(v$Symbol, levels = names(sort(table(v$Symbol))))
  key <- v[, c("ID", "Symbol")]
  v$Type[duplicated(key) | duplicated(key, fromLast = TRUE)] <- "Multiple hits"
  v$Type <- factor(v$Type, levels = MUT_TYPE_LEVELS)
  v[!duplicated(key), c("ID", "Symbol", "Type", "CHIP")]
}

# Sample x gene heatmap faceted by group. Samples without calls are shown as empty
# columns; within each group, samples are sorted by their most frequently mutated genes.
plot_mutation_heatmap <- function(calls, samples, group_names, colors, y_title = waiver()) {
  lv  <- names(group_names)
  lab <- setNames(sprintf("%s (N=%d)", group_names, count_by_group(samples, lv)), lv)

  carriers <- calls[calls$ID %in% samples$ID, c("ID", "Symbol", "Type")]
  carriers$score <- 2^as.integer(carriers$Symbol)   # rank-weighted score for column ordering
  empty_ids <- setdiff(samples$ID, carriers$ID)
  noncarriers <- data.frame(ID     = empty_ids,
                            Symbol = rep(calls$Symbol[1], length(empty_ids)),
                            Type   = factor(rep(NA, length(empty_ids)), levels = MUT_TYPE_LEVELS),
                            score  = rep(0, length(empty_ids)))
  df <- rbind(carriers, noncarriers)
  df$Group <- factor(lab[as.character(samples$Group[match(df$ID, samples$ID)])], levels = lab)
  df <- df %>% group_by(Group) %>% complete(ID, Symbol, fill = list(score = 0)) %>% ungroup()

  ggplot(df, aes(x = reorder(ID, score, sum, decreasing = TRUE), y = Symbol, fill = Type)) +
    geom_tile(colour = "grey", linewidth = 0.2) +
    facet_grid(~ Group, scales = "free", space = "free_x") +
    scale_fill_manual(values = colors, na.value = "white") + ylab(y_title) + theme_bw() +
    theme(axis.text.x = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank())
}

# Tidy logistic-regression coefficients (intercept dropped, CI = estimate +/- 2 SE)
tidy_glm <- function(fit) {
  co <- data.frame(coef(summary(fit)))[-1, ]
  colnames(co) <- c("Estimate", "SE", "Zvalue", "Pvalue")
  co$Coefficient <- factor(rownames(co), levels = rev(rownames(co)))
  co$LowCI  <- co$Estimate - 2 * co$SE
  co$HighCI <- co$Estimate + 2 * co$SE
  co
}

# ------------------------------------------------------------------------------
# Plot helpers
# ------------------------------------------------------------------------------
save_pdf <- function(file, plots, width, height) {
  pdf(file, width = width, height = height)
  on.exit(dev.off())
  for (p in plots) print(p)
}

plot_class_rate <- function(tab, rate, ylab, text_size = 15) {
  ggplot(tab, aes(x = Group, y = .data[[rate]])) +
    geom_bar(aes(fill = interaction(Group, Class)), stat = "identity", color = "black") +
    geom_errorbar(aes(ymin = .data[[paste0(rate, "_LB")]], ymax = .data[[paste0(rate, "_UB")]]), width = 0.2) +
    scale_fill_manual(values = CLASS_FILL, name = "Mut.type") +
    scale_x_discrete(labels = GROUP_LABELS) +
    theme_classic() + theme(text = element_text(size = text_size)) +
    xlab("") + ylab(ylab)
}

plot_or <- function(tab, or, title, ymax = 8) {
  ggplot(tab, aes(x = MAF, y = .data[[or]])) +
    geom_ribbon(aes(ymin = .data[[paste0(or, "_LB")]], ymax = .data[[paste0(or, "_UB")]]), fill = "grey90") +
    geom_line() + geom_point() + geom_hline(yintercept = 1, linetype = "dashed") +
    coord_cartesian(ylim = c(0, ymax)) + xlab("Min VAF") + ylab("Odds ratio") +
    theme_classic() + theme(text = element_text(size = 12)) + ggtitle(title)
}

plot_forest <- function(co) {
  ggplot(co, aes(x = Coefficient, y = Estimate)) +
    geom_errorbar(aes(ymin = LowCI, ymax = HighCI), width = 0.2) +
    geom_point(shape = 16, size = 3) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_x_discrete(labels = TERM_LABELS) + coord_flip() +
    xlab("") + ylab("Estimate") + theme_classic() + theme(text = element_text(size = 12)) +
    ggtitle("Contribution to AD risk")
}

plot_vaf_rate <- function(tab, title) {
  ggplot(tab, aes(x = Type, y = Rate3, color = Group)) +
    geom_pointrange(aes(ymin = Rate3_LB, ymax = Rate3_UB), position = position_dodge(width = 0.5)) +
    scale_color_manual(values = GROUP_COLORS, name = "Group") +
    theme_classic() + theme(text = element_text(size = 10)) +
    xlab("VAF") + ylab("Somatic CHIP mutations per bp") + ggtitle(title)
}

# Stacked proportion bar of a categorical variable (Sex / APOE) per group
plot_composition <- function(df, var, title, lv) {
  df$Group <- factor(as.character(df$Group), levels = lv, labels = GROUP_LABELS[lv])
  var_levels <- levels(df[[var]])
  ggplot(df, aes(x = Group, fill = interaction(Group, .data[[var]]))) +
    geom_bar(position = "fill") +
    geom_text(aes(label = after_stat(count)), stat = "count",
              position = position_fill(vjust = 0.5), color = "white", size = 5) +
    coord_flip() + xlab("") + ylab("") + ggtitle(title) + theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, face = "italic"),
          legend.position = "top", legend.title = element_blank()) +
    scale_fill_manual(values = unname(c(GROUP_COLORS_DARK[lv], GROUP_COLORS[lv])),
                      labels = rep(var_levels, each = length(lv))) +
    guides(fill = guide_legend(ncol = 2, byrow = TRUE))
}

# Split violin (left half = first group, right half = second group)
GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, draw_group = function(self, data, ...) {
  data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
  left <- data[1, "group"] %% 2 == 1
  newdata <- transform(data, x = if (left) xminv else xmaxv)
  newdata <- newdata[order(if (left) newdata$y else -newdata$y), ]
  newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
  newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
  ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
})
geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ...,
                              trim = TRUE, scale = "area", na.rm = FALSE, show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, position = position,
        show.legend = show.legend, inherit.aes = inherit.aes,
        params = list(trim = trim, scale = scale, na.rm = na.rm, ...))
}


print("=== Step 1: Load Samples, Coverage, Panel Design and Variant Calls ===")
id_all <- load_sample_info(sample_info_file, sample_info_sheet)
cat(sprintf("QC-passed samples: %d\n", nrow(id_all))); print(table(id_all$Group))

# Panel design: total panel length and the fraction covered by CHIP whitelist genes
panel_design <- read.delim(panel_bed_file, header = FALSE, skip = 2, stringsAsFactors = FALSE)[, 1:4]
colnames(panel_design) <- c("Chr", "Start", "End", "Symbol")
panel_design$Length <- panel_design$End - panel_design$Start
chip_genes    <- read.delim(chip_gene_file, header = FALSE, stringsAsFactors = FALSE)[, 1]
chip_fraction <- sum(panel_design$Length[panel_design$Symbol %in% chip_genes]) / sum(panel_design$Length)

# Coverage: number of bases covered by >=200 deduplicated reads per sample
id_coverage <- read.delim(coverage_file, header = FALSE, stringsAsFactors = FALSE)
colnames(id_coverage) <- c("ID", "Depth", "Proportion")
id_coverage$Coverage <- round(id_coverage$Proportion * sum(panel_design$Length))

# Variant calls: keep variant-level columns and attach clinical data from the sample table
variant_cols <- c("Chr", "Start", "End", "Ref", "Alt", "Symbol", "Location", "Type", "Impact",
                  "Ref_num", "Alt_num", "Gnomad", "CHIP", "MAF")
variants_raw <- read.delim(merged_file, stringsAsFactors = FALSE,
                           colClasses = c(Ref = "character", Alt = "character"))
variants_raw <- variants_raw[, c("ID", intersect(variant_cols, colnames(variants_raw)))]
variants_raw$CHIP  <- as.logical(variants_raw$CHIP)
variants_raw$Class <- ifelse(variants_raw$Ref %in% c("A", "C", "G", "T") &
                             variants_raw$Alt %in% c("A", "C", "G", "T"), "snp", "indel")
clin_cols    <- c("ID", "ApoE", "Age", "Sex", "MMSE", "Group", "Batch")
variants_all <- merge(variants_raw, id_all[, clin_cols], by = "ID")   # QC-passed samples only

print("=== Step 2: Coverage QC Across All Samples (ED Fig 2) ===")
df_cov <- merge(id_all, id_coverage, by = "ID")
save_pdf(sprintf("%s.pub.EDFig2.pdf", output_header), width = 6, height = 6, plots = list(
  ggplot(df_cov, aes(Depth + 1, color = Group)) + stat_ecdf() +
    scale_color_manual(values = GROUP_COLORS, labels = GROUP_LABELS) +
    xlab("Depth") + ylab("Empirical cumulative density") +
    theme_classic() + theme(text = element_text(size = 15)),
  ggplot(df_cov, aes(x = Group, y = Proportion, color = Group, fill = Group)) +
    geom_boxplot(outlier.shape = NA, width = 0.7, color = "black", alpha = 0.4) +
    geom_jitter(position = position_jitterdodge(), size = 2, alpha = 0.7) +
    scale_color_manual(values = GROUP_COLORS, labels = GROUP_LABELS) +
    scale_fill_manual(values = GROUP_COLORS, labels = GROUP_LABELS) +
    scale_x_discrete(labels = GROUP_LABELS) + ylim(0, 1) +
    xlab("Cohort") + ylab("Proportion of bases covered by >200X") +
    theme_classic() + theme(text = element_text(size = 15))
))


print("=== Step 3: Propensity-Score Matching of CTRL vs. AD ===")
# Eligible: CTRL and AD samples aged >= MIN_AGE; matched on age, sex and CVD history
id_eligible <- droplevels(id_all[id_all$Group != "C" & id_all$Age >= MIN_AGE, ])
match_obj <- matchit(Group ~ Age + Sex + CVD.history, data = id_eligible,
                     distance = "glm", method = "nearest", caliper = MATCH_CALIPER, ratio = MATCH_RATIO)
print(summary(match_obj))
id_matched <- match.data(match_obj)
id_matched$Group <- factor(as.character(id_matched$Group), levels = MATCH_LEVELS)
print(id_matched %>% group_by(Group) %>% summarize(n = n(), mean_age = mean(Age, na.rm = TRUE)))
print(wilcox.test(Age ~ Group, data = id_matched))

# The recomputed subset should reproduce the published "Selected for subset" flag
print(table(Published = id_all$Selected_published, Recomputed = id_all$ID %in% id_matched$ID))
if (!identical(id_all$Selected_published, id_all$ID %in% id_matched$ID)) {
  warning("Recomputed matched subset differs from 'Selected for subset' in SupplTable1")
}

# Restrict variant calls to the matched subset; all downstream analyses use these
variants <- variants_all[variants_all$ID %in% id_matched$ID, ]
variants$Group <- factor(as.character(variants$Group), levels = MATCH_LEVELS)


print("=== Step 4: Cohort Description Figures (Fig 1A, ED Fig 1) ===")
# Fig 1A: age distribution of CTRL vs. AD before (all samples) and after matching
common_cols <- intersect(names(id_all), names(id_matched))
combined_age <- rbind(cbind(id_all[id_all$Group != "C", common_cols], Dataset = "All"),
                      cbind(id_matched[, common_cols], Dataset = "Matched"))
combined_age$Dataset <- factor(combined_age$Dataset, levels = c("All", "Matched"))
combined_age$Group   <- factor(as.character(combined_age$Group), levels = MATCH_LEVELS)
save_pdf(sprintf("%s.pub.Fig1A.pdf", output_header), width = 4, height = 3, plots = list(
  ggplot(combined_age, aes(x = Dataset, y = Age, fill = Group)) +
    geom_split_violin(trim = FALSE, scale = "count") +
    geom_boxplot(aes(fill = NULL), width = 0.05, outlier.shape = NA) +
    scale_fill_manual(values = GROUP_COLORS[MATCH_LEVELS], name = "Group") +
    labs(x = "", y = "Age", title = "Age Distribution: Before vs After Sub-sampling") + theme_classic()
))

# ED Fig 1: age distribution per group (matched subset, all samples)
violin_age <- function(df, lv, title) {
  ggplot(df, aes(x = Group, y = Age, fill = Group)) + geom_violin(trim = FALSE) +
    stat_summary(fun = median, geom = "point", shape = 95, size = 8, color = "black") +
    scale_x_discrete(limits = rev(lv)) +
    scale_fill_manual(values = GROUP_COLORS[lv], name = "Group") + theme_classic() + ggtitle(title)
}
save_pdf(sprintf("%s.pub.EDFig1.pdf", output_header), width = 4, height = 3, plots = list(
  violin_age(id_matched, MATCH_LEVELS, "Matched Subset"),
  violin_age(id_all, c("C", "CTRL", "AD"), "All Samples")
))

# ED Fig 1: sex and APOE composition per group (matched subset, then all samples)
df_comp <- id_all
df_comp$Sex  <- factor(df_comp$Sex,  levels = c("M", "F"),   labels = c("Male", "Female"))
df_comp$ApoE <- factor(df_comp$ApoE, levels = c("34", "33"), labels = c("e3/e4", "e3/e3"))
df_comp_matched <- df_comp[df_comp$ID %in% id_matched$ID, ]
save_pdf(sprintf("%s.pub.EDFig1.add.pdf", output_header), width = 10, height = 4, plots = list(
  ggarrange(plot_composition(df_comp_matched, "Sex",  "Sex",  MATCH_LEVELS),
            plot_composition(df_comp_matched, "ApoE", "APOE", MATCH_LEVELS), ncol = 2, nrow = 1),
  ggarrange(plot_composition(df_comp, "Sex",  "Sex",  c("C", "CTRL", "AD")),
            plot_composition(df_comp, "ApoE", "APOE", c("C", "CTRL", "AD")), ncol = 2, nrow = 1)
))


print("=== Step 5: Mutation Rate per bp by Variant Class (Fig 1B, 1C) ===")
rate_table <- build_class_rate_table(variants, id_matched)
test_ad_vs_ctrl(rate_table, "Count1", "Coverage",  "All somatic mutations")
test_ad_vs_ctrl(rate_table, "Count3", "Coverage3", "CHIP mutations")

save_pdf(sprintf("%s.pub.Fig1B.pdf", output_header), width = 4, height = 3,
         plots = list(plot_class_rate(rate_table, "Rate1", "Somatic mutations per bp")))
save_pdf(sprintf("%s.pub.Fig1C.pdf", output_header), width = 4, height = 3,
         plots = list(plot_class_rate(rate_table, "Rate3", "Somatic CHIP mutations per bp")))


print("=== Step 6: Mutation Landscape Heatmaps (Fig 1D, ED Fig 4) ===")
if (!all(c("Location", "Type") %in% colnames(variants_all))) {
  stop("Heatmaps need the 'Location' and 'Type' columns in ", merged_file)
}
heatmap_calls <- prepare_heatmap_calls(variants_all)

# ED Fig 4: all somatic mutations in all QC-passed samples
save_pdf(sprintf("%s.pub.EDFig4.pdf", output_header), width = 15, height = 10, plots = list(
  plot_mutation_heatmap(heatmap_calls, id_all,
                        c(AD = "Alzheimer's disease", CTRL = "Control (Age-matched)", C = "Centenarians"),
                        HEATMAP_COLORS_ALL, y_title = "")
))

# Fig 1D: CHIP mutations in the matched subset (genes re-ordered by CHIP call frequency)
chip_calls <- heatmap_calls[heatmap_calls$ID %in% id_matched$ID & heatmap_calls$CHIP %in% TRUE, ]
chip_symbols <- droplevels(chip_calls$Symbol)
chip_calls$Symbol <- factor(chip_symbols, levels = names(sort(table(chip_symbols))))
save_pdf(sprintf("%s.pub.Fig1D.pdf", output_header), width = 10, height = 2, plots = list(
  plot_mutation_heatmap(chip_calls, id_matched,
                        c(AD = "Alzheimer's disease (Matched)", CTRL = "Control (Matched)"),
                        HEATMAP_COLORS_CHIP)
))


print("=== Step 7: Logistic Regression of AD Status on CHIP Burden (Fig 1E) ===")
# Per-sample CHIP mutation rate (CHIP mutations per bp of CHIP genes covered >=200X)
sample_burden <- id_matched
sample_burden$Freq_chip <- as.vector(table(factor(variants$ID[which(variants$CHIP)], levels = sample_burden$ID)))
sample_burden <- merge(sample_burden, id_coverage, by = "ID")
sample_burden$rate_chip <- sample_burden$Freq_chip / (sample_burden$Coverage * chip_fraction)

glm_formulas <- list(
  "CHIP burden"       = Group ~ scale(rate_chip) + Sex + scale(Age) + factor(ApoE),
  "CHIP burden + CVD" = Group ~ scale(rate_chip) + factor(CVD.history) + Sex + scale(Age) + factor(ApoE)
)
glm_coefs <- lapply(names(glm_formulas), function(nm) {
  fit <- glm(glm_formulas[[nm]], data = sample_burden, family = "binomial")
  cat(sprintf("\n--- Model: %s ---\n", nm)); print(summary(fit))
  tidy_glm(fit)
})
save_pdf(sprintf("%s.pub.Fig1E.pdf", output_header), width = 3, height = 4, plots = lapply(glm_coefs, plot_forest))


print("=== Step 8: Mutation Rate Stratified by APOE Genotype (Fig 1F) ===")
rate_table_apoe <- do.call(rbind, lapply(c("34", "33"), function(g) {
  tab <- build_class_rate_table(variants[variants$ApoE == g, ], id_matched[id_matched$ApoE == g, ])
  tab$Type <- paste0("ApoE", g)
  test_ad_vs_ctrl(tab, "Count3", "Coverage3", sprintf("CHIP mutations, APOE %s", g))
  tab
}))

save_pdf(sprintf("%s.pub.Fig1F.pdf", output_header), width = 4, height = 3, plots = list(
  plot_class_rate(rate_table_apoe, "Rate1", "Somatic mutations per bp", text_size = 10) + facet_grid(. ~ Type),
  plot_class_rate(rate_table_apoe, "Rate3", "Somatic CHIP mutations per bp", text_size = 10) + facet_grid(. ~ Type)
))


print("=== Step 9: Odds Ratio by Minimum VAF (Fig 1G) ===")
id_matched_33 <- id_matched[id_matched$ApoE == "33", ]
variants_33   <- variants[variants$ApoE == "33", ]
or_table    <- build_or_table(variants, id_matched)
or_table_33 <- build_or_table(variants_33, id_matched_33)

save_pdf(sprintf("%s.pub.Fig1G.pdf", output_header), width = 6, height = 3, plots = list(
  plot_or(or_table,    "OR1", "All somatic mutations"),
  plot_or(or_table,    "OR2", "Somatic CHIP mutations"),
  plot_or(or_table_33, "OR1", "All somatic mutations - e3/e3 only"),
  plot_or(or_table_33, "OR2", "Somatic CHIP mutations - e3/e3 only", ymax = 10)
))


print("=== Step 10: CHIP Mutation Rate by VAF Bin (Fig 1H) ===")
vaf_table    <- build_vaf_rate_table(variants, id_matched)
vaf_table_33 <- build_vaf_rate_table(variants_33, id_matched_33)
for (b in VAF_BIN_LABELS) {
  test_ad_vs_ctrl(vaf_table[vaf_table$Type == b, ],       "Count3", "Coverage3", sprintf("CHIP, %s", b))
  test_ad_vs_ctrl(vaf_table_33[vaf_table_33$Type == b, ], "Count3", "Coverage3", sprintf("CHIP, %s, e3/e3 only", b))
}

save_pdf(sprintf("%s.pub.Fig1H.pdf", output_header), width = 5, height = 2, plots = list(
  plot_vaf_rate(vaf_table,    "CHIP mutation burden by VAF"),
  plot_vaf_rate(vaf_table_33, "CHIP mutation burden by VAF - e3/e3 only")
))

# Heterogeneity of the AD/CTRL rate ratio across VAF bins (e3/e3 only), Mantel-Haenszel
meta_data <- data.frame(
  vaf_bins = VAF_BIN_LABELS,
  event.e  = vaf_table_33$Count3[vaf_table_33$Group == "AD"],
  n.e      = vaf_table_33$Coverage3[vaf_table_33$Group == "AD"],
  event.c  = vaf_table_33$Count3[vaf_table_33$Group == "CTRL"],
  n.c      = vaf_table_33$Coverage3[vaf_table_33$Group == "CTRL"]
)
meta_q_test <- metabin(event.e = event.e, n.e = n.e, event.c = event.c, n.c = n.c, studlab = vaf_bins,
                       data = meta_data, sm = "RR", method = "MH", incr = 0.5)
print(meta_q_test)
# Subgroup comparison: low (<=5%) vs. high (>5%) VAF bins
print(update(meta_q_test, subgroup = c("Lower_VAF", "Lower_VAF", "Higher_VAF", "Higher_VAF"),
             print.subgroup.name = FALSE))

print("=== CHIP burden analysis finished ===")