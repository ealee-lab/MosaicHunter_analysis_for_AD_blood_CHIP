# ==============================================================================
# Pipeline: Age Trends of Somatic / CHIP Mutations and Association with MMSE
# Description: Uses all QC-passed samples (no matching) to describe how somatic and
#              CHIP mutations change with age in CTRL, C and AD, and how CHIP status
#              relates to cognitive score (MMSE) within AD:
#                - Age distribution by group                          Fig 1Aa
#                - CHIP carrier proportion by age (+ trend lm)        Fig 2A
#                - CHIP mutation rate per bp by age group             Fig 2B
#                - VAF of somatic / CHIP mutations vs. age            ED Fig 8
#                - CHIP carriers / burden by MMSE and APOE (AD only)  ED Fig 10
#                  and logistic regression of CHIP status on MMSE
#                - Global dN/dS, e3/e3 SNVs (dNdScv)                  Fig 2C
#                - Per-gene dN/dS of selected CHIP genes, e3/e3 SNVs  Fig 2D
#                - Per-gene dN/dS significance, e3/e3 SNVs            ED Fig 9
#
# Inputs (same as Scr3):
#   AD_panel.merged.addage.tsv          Merged SNV + indel calls (Scr2 output)
#   SupplTable1_Sample_Information.xlsx Sample-level clinical metadata for the 353
#                                       QC-passed samples (sheet "Sample_Information")
#   AD_panel.coverage.summary           Per-sample mean depth and fraction of
#                                       bases >=200X (Scr1 output)
#   3277911_Covered.bed                 Capture panel design (2 header lines)
#   whitelist_Genename.tsv              CHIP whitelist gene names (Scr2 output)
#   CancerGeneCensus_01212022.tsv       Optional: COSMIC Cancer Gene Census, used only to
#                                       label TSGs in the dNdScv gene table
#
# Group definition (column "Group" of SupplTable1):
#   AD   : clinically diagnosed AD
#   C    : non-AD samples older than 90 y (labelled "Centenarian" in figures)
#   CTRL : all other non-AD samples
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggbreak)
  library(readxl)
  library(dndscv)
  library(ggrepel)
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
cgc_file          <- "CancerGeneCensus_01212022.tsv"   # optional
output_header     <- "AD_panel.merged.addage.pub.clinFilt"

# Age bins: (0,60], (60,70], (70,80], (80,90], (90,110]; plotted at the bin centres
AGE_BREAKS     <- c(0, 60, 70, 80, 90, 110)
AGE_LABELS     <- c("50-60", "60-70", "70-80", "80-90", "C(90-110)")
AGE_CENTRES    <- c(55, 65, 75, 85, 100)
AGE_AXIS_BREAK <- c(90, 95)             # gap drawn on the age axis before the oldest bin
MIN_CELL_SIZE  <- 5                     # group x age cells need more samples than this

# MMSE bins: [0,5), [5,10), ..., [25,30]
MMSE_BREAKS  <- seq(0, 30, by = 5)
MMSE_LABELS  <- c("[0,5]", "[5,10]", "[10,15]", "[15,20]", "[20,25]", "[25,30]")
MMSE_CENTRES <- seq(2.5, 27.5, by = 5)

# The only AD case older than 90 y; excluded from the VAF-vs-age plots (ED Fig 8)
VAF_AGE_EXCLUDE <- "RESEARCH_03_S86"

GROUP_LEVELS <- c("CTRL", "C", "AD")
GROUP_COLORS <- c(CTRL = "royalblue1", C = "darkblue", AD = "indianred3")
APOE_LABELS  <- c("33" = "e3/e3", "34" = "e3/e4")
APOE_COLORS  <- c("e3/e3" = "royalblue1", "e3/e4" = "indianred3")
# Stacked-bar fills: group (or APOE) x variant class
AGE_CLASS_FILL  <- c("CTRL.indel" = "royalblue4", "C.indel" = "black",    "AD.indel" = "indianred4",
                     "CTRL.snp"   = "royalblue1", "C.snp"   = "darkblue", "AD.snp"   = "indianred3")
APOE_CLASS_FILL <- c("e3/e3.indel" = "royalblue1", "e3/e4.indel" = "indianred3",
                     "e3/e3.snp"   = "royalblue4", "e3/e4.snp"   = "indianred4")
CLASS_SHAPES    <- c(indel = 24, snp = 21)

# dNdScv settings
DNDS_GENE_COLS <- c(Gene = "gene_name", N_synonymous = "n_syn", N_missense = "n_mis",
                    N_nonsense = "n_non", N_splicing = "n_spl",
                    Ratio_missense = "wmis_cv", Ratio_nonsense = "wnon_cv", Ratio_splicing = "wspl_cv",
                    Pvalue = "pallsubs_cv", Qvalue = "qallsubs_cv")
DNDS_GLOBAL_TYPES <- c(wmis = "Missense", wnon = "Nonsense", wspl = "Splicing", wtru = "Truncating", wall = "All")
DIRECTION_SHAPES  <- c("+" = 1, "-" = 4)

# Axis labels for the MMSE regression coefficients (ED Fig 10)
TERM_LABELS <- c("scale(MMSE)"    = "MMSE",
                 "factor(ApoE)34" = "APOE: e3/e4",
                 "scale(Age)"     = "Age",
                 "SexM"           = "Sex: Male")

dodge_val <- position_dodge(width = 4, preserve = "single")

# ==============================================================================
# Helper Functions
# ==============================================================================

# Load sample-level clinical metadata (SupplTable1) and rename columns for analysis
load_sample_info <- function(path, sheet) {
  raw <- data.frame(read_excel(path, sheet = sheet), check.names = FALSE)
  info <- data.frame(
    ID          = raw[["ID"]],
    ApoE        = as.character(raw[["APOE genotype"]]),                 # "33" / "34"
    Age         = raw[["Age (at blood collection)"]],
    Sex         = raw[["Sex"]],                                         # "F" / "M"
    MMSE        = suppressWarnings(as.numeric(raw[["MMSE"]])),          # AD cases only
    Group       = factor(raw[["Group"]], levels = GROUP_LEVELS),
    CVD.history = raw[["CVD history"]],
    stringsAsFactors = FALSE
  )
  info[order(info$ID), ]
}

# 95% binomial CI of k / n (lower, upper)
prop_ci <- function(k, n) suppressWarnings(prop.test(k, n))$conf.int

# CHIP carrier proportions with 95% CIs, per stratum
summarise_carriers <- function(df, by) {
  df %>% group_by(across(all_of(by))) %>%
    summarise(Prop_chip    = mean(Carrier_chip),
              Prop_chip_LB = prop_ci(sum(Carrier_chip), n())[1],
              Prop_chip_UB = prop_ci(sum(Carrier_chip), n())[2],
              Total        = n(), .groups = "drop") %>%
    as.data.frame()
}

# CHIP mutation rate per bp in every combination of two strata (e.g. Group x Age_group),
# split by variant class. Rate3 bar heights use class-specific counts (stacked);
# CIs use the class-combined CHIP count.
build_strata_rate_table <- function(v, s, strata) {
  lv    <- lapply(strata, function(x) levels(s[[x]]))
  cells <- expand.grid(setNames(lapply(lv, function(l) factor(l, levels = l)), strata))
  count_cells <- function(df) {
    as.vector(table(factor(df[[strata[1]]], levels = lv[[1]]), factor(df[[strata[2]]], levels = lv[[2]])))
  }
  cov  <- merge(s[, c("ID", strata)], id_coverage[, c("ID", "Coverage")], by = "ID")
  base <- cbind(cells,
                Count3   = count_cells(v[which(v$CHIP), ]),
                Total    = count_cells(s),
                Coverage = as.vector(tapply(cov$Coverage, list(cov[[strata[1]]], cov[[strata[2]]]), sum)))
  tab <- do.call(rbind, lapply(c("snp", "indel"), function(cl) {
    vc <- v[which(v$Class == cl & v$CHIP), ]
    cbind(base, Class = cl, Count_sep3 = count_cells(vc))
  }))
  tab <- tab[!is.na(tab$Coverage), ]
  tab$Coverage3 <- tab$Coverage * chip_fraction
  ci <- t(mapply(prop_ci, tab$Count3, tab$Coverage3))
  tab$Rate3    <- tab$Count_sep3 / tab$Coverage3
  tab$Rate3_LB <- ci[, 1]
  tab$Rate3_UB <- ci[, 2]
  rownames(tab) <- NULL
  tab
}

# One-sided test (AD > CTRL) of CHIP mutation rate, plus the AD/CTRL rate ratio
test_ad_vs_ctrl <- function(tab, label) {
  ad   <- tab[tab$Group == "AD", ][1, ]
  ctrl <- tab[tab$Group == "CTRL", ][1, ]
  cat(sprintf("\n--- %s ---\n", label))
  print(prop.test(c(ad$Count3, ctrl$Count3), c(ad$Coverage3, ctrl$Coverage3), alternative = "greater"))
  cat(sprintf("Rate ratio (AD/CTRL): %.3f\n", (ad$Count3 / ad$Coverage3) / (ctrl$Count3 / ctrl$Coverage3)))
}

# Tidy regression coefficients (intercept dropped, 95% CI = estimate +/- 1.96 SE)
tidy_model <- function(fit) {
  co <- data.frame(coef(summary(fit)))[-1, ]
  colnames(co) <- c("Estimate", "SE", "Stat", "Pvalue")
  co$Coefficient <- factor(rownames(co), levels = rev(rownames(co)))
  co$LowCI  <- co$Estimate - 1.96 * co$SE
  co$HighCI <- co$Estimate + 1.96 * co$SE
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

# Carrier proportion per age bin and group, with a per-group linear trend.
#   trend = "smooth": dashed lm line (geom_smooth)
#   trend = "line"  : semi-transparent dashed lm line, with an extra "NA" tick at 110
plot_carrier_by_age <- function(tab, y, ylab, title, trend = c("smooth", "line")) {
  trend <- match.arg(trend)
  x_breaks <- if (trend == "line") c(AGE_CENTRES, 110) else AGE_CENTRES
  x_labels <- if (trend == "line") c(AGE_LABELS, "NA") else AGE_LABELS
  trend_layer <- if (trend == "line") {
    geom_line(stat = "smooth", method = "lm", formula = y ~ x, se = FALSE, linetype = "dashed", alpha = 0.6)
  } else {
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linetype = "dashed")
  }
  ggplot(tab, aes(x = Age_group_val, y = .data[[y]], group = Group, color = Group, fill = Group)) +
    geom_bar(stat = "identity", position = dodge_val, width = 4) +
    geom_errorbar(aes(ymin = .data[[paste0(y, "_LB")]], ymax = .data[[paste0(y, "_UB")]]),
                  color = "black", width = 2, position = dodge_val) +
    scale_x_continuous(breaks = x_breaks, labels = x_labels) + scale_x_break(AGE_AXIS_BREAK) +
    scale_color_manual(values = GROUP_COLORS) + scale_fill_manual(values = GROUP_COLORS) +
    trend_layer + theme_classic() + ylim(0, 1) +
    xlab("Age") + ylab(ylab) + ggtitle(title)
}

# Stacked CHIP rate per group, faceted by age bin (C and the youngest bin drawn narrower)
plot_rate_by_age <- function(tab, title) {
  tab$Mut.type <- factor(paste(tab$Group, tab$Class, sep = "."), levels = names(AGE_CLASS_FILL))
  tab$w <- ifelse(tab$Group == "C" | tab$Age_group == AGE_LABELS[1], 0.45, 0.9)
  ggplot(tab, aes(x = Group, y = Rate3, fill = Mut.type)) +
    suppressWarnings(geom_bar(aes(width = w), stat = "identity", position = "stack")) +  # per-row bar width
    geom_errorbar(aes(ymin = Rate3_LB, ymax = Rate3_UB), width = 0.2) +
    facet_grid(. ~ Age_group, scales = "free") +
    scale_fill_manual(values = AGE_CLASS_FILL, name = "Mut.type") +
    theme_classic() + xlab("Age") + ylab("Somatic CHIP mutations per bp") + ggtitle(title)
}

# Per-gene dN/dS significance for each group (genes along x, sorted by AD P-value)
plot_dnds_genes <- function(df, y, sig_col, ylab) {
  ggplot(df, aes(x = Gene, y = .data[[y]], color = Group)) +
    geom_point(aes(shape = Direction), size = 2) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    geom_text_repel(data = df[which(df[[sig_col]] < 0.05), ], aes(label = Gene), max.overlaps = Inf) +
    scale_shape_manual(values = DIRECTION_SHAPES) + scale_color_manual(values = GROUP_COLORS) +
    theme_classic() + theme(text = element_text(size = 15), axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    ylab(ylab)
}

plot_dnds_ratio <- function(df) {
  ggplot(df, aes(x = Type, y = MLE, fill = Group, color = Group)) +
    geom_pointrange(aes(ymin = LowCI, ymax = HighCI), position = position_dodge(width = 0.5)) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    scale_fill_manual(values = GROUP_COLORS) + scale_color_manual(values = GROUP_COLORS) +
    theme_classic() + theme(text = element_text(size = 10)) + xlab("Mutation type") + ylab("dN/dS ratio")
}

plot_forest <- function(co, ylab, title) {
  ggplot(co, aes(x = Coefficient, y = Estimate)) +
    geom_errorbar(aes(ymin = LowCI, ymax = HighCI), width = 0.2) +
    geom_point(shape = 16, size = 3) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_x_discrete(labels = TERM_LABELS) + coord_flip() +
    xlab("") + ylab(ylab) + theme_classic() + theme(text = element_text(size = 12)) + ggtitle(title)
}


print("=== Step 1: Load Samples, Coverage, Panel Design and Variant Calls ===")
id_all <- load_sample_info(sample_info_file, sample_info_sheet)
id_all$Age_group <- cut(id_all$Age, breaks = AGE_BREAKS, labels = AGE_LABELS, right = TRUE)
cat(sprintf("QC-passed samples: %d\n", nrow(id_all))); print(table(id_all$Group, id_all$Age_group))

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

# Variant calls with clinical data from the sample table (QC-passed samples only)
variants_raw <- read.delim(merged_file, stringsAsFactors = FALSE,
                           colClasses = c(Ref = "character", Alt = "character"))
variants_raw <- variants_raw[, c("ID", "Chr", "Start", "Symbol", "Ref", "Alt", "CHIP", "MAF")]
variants_raw$CHIP  <- as.logical(variants_raw$CHIP)
variants_raw$Class <- ifelse(variants_raw$Ref %in% c("A", "C", "G", "T") &
                             variants_raw$Alt %in% c("A", "C", "G", "T"), "snp", "indel")
variants_all <- merge(variants_raw, id_all[, c("ID", "ApoE", "Age", "Age_group", "Group")], by = "ID")

# Per-sample CHIP mutation counts and carrier status (samples without coverage are dropped)
sample_burden <- id_all
sample_burden$Freq_chip    <- as.vector(table(factor(variants_all$ID[which(variants_all$CHIP)], levels = sample_burden$ID)))
sample_burden$Carrier_chip <- as.integer(sample_burden$Freq_chip > 0)
sample_burden$Age_group_val <- AGE_CENTRES[as.integer(sample_burden$Age_group)]
sample_burden <- merge(sample_burden, id_coverage, by = "ID")


print("=== Step 2: Age Distribution by Group (Fig 1Aa) ===")
save_pdf(sprintf("%s.pub.Fig1Aa.pdf", output_header), width = 4, height = 3, plots = list(
  ggplot(id_all, aes(x = Age, fill = Group)) + geom_histogram(bins = 10, position = "dodge") +
    scale_fill_manual(values = GROUP_COLORS, name = "Group") + theme_classic()
))


print("=== Step 3: CHIP Carrier Proportion by Age (Fig 2A) ===")
# Carrier proportions per group x age bin; cells with <= MIN_CELL_SIZE samples are dropped
carrier_age    <- summarise_carriers(sample_burden, c("Group", "Age_group_val"))
carrier_age    <- carrier_age[carrier_age$Total > MIN_CELL_SIZE, ]
carrier_age_33 <- summarise_carriers(sample_burden[sample_burden$ApoE == "33", ], c("Group", "Age_group_val"))
carrier_age_33 <- carrier_age_33[carrier_age_33$Total > MIN_CELL_SIZE, ]

# Does the age trend of CHIP carrier proportion differ between AD and CTRL? (bin-level lm)
cat("\n--- CHIP carrier proportion ~ age bin x group (CTRL vs AD) ---\n")
print(summary(lm(Prop_chip ~ Age_group_val * Group, data = droplevels(carrier_age[carrier_age$Group != "C", ]))))
cat("\n--- CHIP carrier proportion ~ age bin x group (CTRL vs AD), e3/e3 only ---\n")
print(summary(lm(Prop_chip ~ Age_group_val * Group, data = droplevels(carrier_age_33[carrier_age_33$Group != "C", ]))))

save_pdf(sprintf("%s.pub.Fig2A.pdf", output_header), width = 6, height = 3, plots = list(
  plot_carrier_by_age(carrier_age,    "Prop_chip", "Proportion of CHIP mutation carriers",
                      "CHIP carrier proportion by age", trend = "line"),
  plot_carrier_by_age(carrier_age,    "Prop_chip", "Proportion of CHIP mutation carriers",
                      "CHIP carrier proportion by age", trend = "smooth"),
  plot_carrier_by_age(carrier_age_33, "Prop_chip", "Proportion of CHIP mutation carriers",
                      "CHIP carrier proportion by age - e3/e3 only", trend = "line"),
  plot_carrier_by_age(carrier_age_33, "Prop_chip", "Proportion of CHIP mutation carriers",
                      "CHIP carrier proportion by age - e3/e3 only", trend = "smooth")
))


print("=== Step 4: CHIP Mutation Rate per bp by Age Group (Fig 2B) ===")
# Group x age cells with <= MIN_CELL_SIZE samples are dropped
rate_age <- build_strata_rate_table(variants_all, id_all, c("Group", "Age_group"))
rate_age <- rate_age[rate_age$Total > MIN_CELL_SIZE, ]
rate_age_33 <- build_strata_rate_table(variants_all[variants_all$ApoE == "33", ],
                                       id_all[id_all$ApoE == "33", ], c("Group", "Age_group"))
rate_age_33 <- rate_age_33[rate_age_33$Total > MIN_CELL_SIZE, ]

save_pdf(sprintf("%s.pub.Fig2B.pdf", output_header), width = 6, height = 3, plots = list(
  plot_rate_by_age(rate_age,    "CHIP mutation burden by Age"),
  plot_rate_by_age(rate_age_33, "CHIP mutation burden by Age - e3/e3 only")
))

# AD vs CTRL within every age bin that retains both groups
for (tab_name in c("rate_age", "rate_age_33")) {
  tab <- get(tab_name)
  for (b in AGE_LABELS) {
    tb <- tab[tab$Age_group == b, ]
    if (all(c("AD", "CTRL") %in% tb$Group)) {
      test_ad_vs_ctrl(tb, sprintf("CHIP mutations, age %s%s", b, if (tab_name == "rate_age_33") ", e3/e3 only" else ""))
    }
  }
}


print("=== Step 5: VAF vs. Age (ED Fig 8) ===")
vaf_age      <- variants_all[variants_all$ID != VAF_AGE_EXCLUDE, ]
vaf_age_chip <- vaf_age[which(vaf_age$CHIP), ]   # for the VAF model below

save_pdf(sprintf("%s.pub.EDFig8.pdf", output_header), width = 6, height = 3, plots = list(
  # All somatic mutations, with linear trend per group
  ggplot(vaf_age, aes(x = Age, y = MAF, fill = Group)) +
    geom_point(aes(shape = Class, col = Group), size = 1) +
    geom_smooth(method = "lm", formula = y ~ x, color = "black", se = TRUE, linetype = "dashed") +
    scale_color_manual(values = GROUP_COLORS) + scale_fill_manual(values = GROUP_COLORS) +
    scale_shape_manual(values = CLASS_SHAPES) + ylab("VAF") + theme_classic()
))

cat("\n--- CHIP mutation VAF ~ group + age ---\n")
print(summary(lm(MAF ~ Group + Age, data = vaf_age_chip)))


print("=== Step 6: CHIP Mutations and MMSE in AD (ED Fig 10) ===")
# All AD cases, binned into 5-point MMSE groups (cases without MMSE stay in the data so
# that scale(Age) in the models below is standardised over all AD cases)
mmse_bin <- function(x) cut(x, breaks = MMSE_BREAKS, labels = MMSE_LABELS, include.lowest = TRUE, right = FALSE)
burden_AD <- sample_burden[sample_burden$Group == "AD", ]
burden_AD$MMSE_group     <- mmse_bin(burden_AD$MMSE)
burden_AD$MMSE_group_val <- MMSE_CENTRES[as.integer(burden_AD$MMSE_group)]
burden_AD$ApoE_label     <- factor(APOE_LABELS[burden_AD$ApoE], levels = APOE_LABELS)
burden_AD$CHIP_burden    <- burden_AD$Freq_chip / (burden_AD$Coverage * chip_fraction)

carrier_mmse_apoe <- summarise_carriers(burden_AD[!is.na(burden_AD$MMSE), ], c("ApoE_label", "MMSE_group_val"))

# CHIP rate per bp by APOE x MMSE group (AD only)
samples_AD <- id_all[id_all$Group == "AD" & !is.na(id_all$MMSE), ]
samples_AD$MMSE_group <- mmse_bin(samples_AD$MMSE)
samples_AD$ApoE_label <- factor(APOE_LABELS[samples_AD$ApoE], levels = APOE_LABELS)
variants_AD <- merge(variants_all[, c("ID", "CHIP", "Class")], samples_AD[, c("ID", "ApoE_label", "MMSE_group")], by = "ID")
rate_mmse <- build_strata_rate_table(variants_AD, samples_AD, c("ApoE_label", "MMSE_group"))
rate_mmse$Mut.type <- factor(paste(rate_mmse$ApoE_label, rate_mmse$Class, sep = "."), levels = names(APOE_CLASS_FILL))

# Is MMSE associated with CHIP carrier status / CHIP burden in AD, adjusting for APOE, age and sex?
model_carrier <- glm(Carrier_chip ~ scale(MMSE) + factor(ApoE) + scale(Age) + Sex, data = burden_AD, family = binomial)
model_burden  <- lm(scale(CHIP_burden) ~ scale(MMSE) + factor(ApoE) + scale(Age) + Sex, data = burden_AD)
cat("\n--- CHIP carrier status ~ MMSE (AD) ---\n"); print(summary(model_carrier))
cat("\n--- CHIP burden ~ MMSE (AD) ---\n");         print(summary(model_burden))

mmse_axis <- scale_x_continuous(breaks = MMSE_CENTRES, labels = MMSE_LABELS)
save_pdf(sprintf("%s.pub.EDFig10.pdf", output_header), width = 6, height = 5, plots = list(
  # CHIP carrier proportion by MMSE and APOE genotype
  ggplot(carrier_mmse_apoe, aes(x = MMSE_group_val, y = Prop_chip, group = ApoE_label, color = ApoE_label, fill = ApoE_label)) +
    geom_bar(stat = "identity", position = dodge_val, width = 4) +
    geom_errorbar(aes(ymin = Prop_chip_LB, ymax = Prop_chip_UB), color = "black", width = 2, position = dodge_val) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linetype = "dashed") +
    scale_color_manual(values = APOE_COLORS, name = "ApoE") + scale_fill_manual(values = APOE_COLORS, name = "ApoE") +
    mmse_axis + theme_classic() + ylim(0, 1) +
    xlab("MMSE") + ylab("Proportion of CHIP mutation carriers"),
  # CHIP mutation rate per bp by APOE x MMSE group
  ggplot(rate_mmse, aes(x = ApoE_label, y = Rate3, fill = Mut.type)) +
    geom_bar(stat = "identity", position = "stack") +
    geom_errorbar(aes(ymin = Rate3_LB, ymax = Rate3_UB), width = 0.2) +
    facet_grid(. ~ MMSE_group, scales = "free") +
    scale_fill_manual(values = APOE_CLASS_FILL, name = "Mut.type") +
    theme_classic() + xlab("ApoE x MMSE") + ylab("Somatic CHIP mutations per bp"),
  # Logistic regression of CHIP carrier status
  plot_forest(tidy_model(model_carrier), "Log-Odds Estimate (CHIP Carrier)", "Predictors of CHIP Mutation Presence")
))

print("=== Step 7: dN/dS Selection Analysis on e3/e3 SNVs (dNdScv) ===")
# Genes targeted by the panel (non-gene backbone targets are named "chr...")
panel_genes <- unique(panel_design$Symbol)
panel_genes <- panel_genes[!grepl("chr", panel_genes)]

# dNdScv per group on SNVs from APOE e3/e3 carriers (all QC-passed samples, no matching).
# Every table and figure below (gene table, ED Fig 9, Fig 2C, Fig 2D) is derived from dnds_33.
snv_33 <- variants_all[variants_all$Class == "snp" & variants_all$ApoE == "33", ]
stopifnot(all(snv_33$ID %in% id_all$ID[id_all$ApoE == "33"]))
dnds_33 <- lapply(setNames(nm = c("AD", "CTRL", "C")), function(g) {
  muts <- snv_33[which(snv_33$Group == g), c("ID", "Chr", "Start", "Ref", "Alt")]
  dndscv(muts, gene_list = panel_genes, outmats = TRUE,
         max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf)
})

# Gene-level selection (e3/e3): counts, dN/dS ratios and P/Q-values (all substitution types)
dnds_gene <- do.call(rbind, lapply(names(dnds_33), function(g) {
  x <- dnds_33[[g]]$sel_cv[, DNDS_GENE_COLS]
  colnames(x) <- names(DNDS_GENE_COLS)
  cbind(x, Group = g)
}))
dnds_gene$Gene    <- factor(dnds_gene$Gene, levels = rev(dnds_33$AD$sel_cv$gene_name))
dnds_gene$Group   <- factor(dnds_gene$Group, levels = GROUP_LEVELS)
dnds_gene$Pconvert <- -log10(dnds_gene$Pvalue)
dnds_gene$Qconvert <- -log10(dnds_gene$Qvalue)
dnds_gene$Direction <- factor(ifelse(dnds_gene$Ratio_missense > 1 | dnds_gene$Ratio_nonsense > 1, "+", "-"),
                              levels = names(DIRECTION_SHAPES))

# Gene type: tumour suppressors (COSMIC CGC, if available) or other
tsg_genes <- character(0)
if (file.exists(cgc_file)) {
  cgc <- read.delim(cgc_file, stringsAsFactors = FALSE)
  tsg_genes <- cgc$Gene.Symbol[grepl("TSG", cgc$Role.in.Cancer) & !grepl("oncogene", cgc$Role.in.Cancer)]
} else {
  message("Cancer Gene Census not found (", cgc_file, "); no genes are labelled TSG")
}
dnds_gene$Genetype <- factor(ifelse(dnds_gene$Gene %in% tsg_genes, "TSG", "Other"), levels = c("TSG", "Other"))

# Table of genes with P < 0.05 in any group
write.table(dnds_gene[dnds_gene$Gene %in% dnds_gene$Gene[which(dnds_gene$Pvalue < 0.05)], ],
            file = sprintf("%s.dNdScv.gene.tsv", output_header), quote = FALSE, sep = "\t", row.names = FALSE)

# ED Fig 9: per-gene dN/dS significance by group (e3/e3)
save_pdf(sprintf("%s.pub.EDFig9.pdf", output_header), width = 8, height = 3, plots = list(
  plot_dnds_genes(dnds_gene, "Pconvert", "Pvalue", "-log10(P-value)"),
  plot_dnds_genes(dnds_gene, "Qconvert", "Qvalue", "-log10(adjusted P-value)")
))

# Fig 2C: global dN/dS across all panel genes (e3/e3)
dnds_global <- do.call(rbind, lapply(names(dnds_33), function(g) {
  x <- dnds_33[[g]]$globaldnds
  data.frame(Type = DNDS_GLOBAL_TYPES[as.character(x$name)], MLE = x$mle,
             LowCI = x$cilow, HighCI = x$cihigh, Group = g, row.names = NULL)
}))
dnds_global$Type  <- factor(dnds_global$Type, levels = c("All", "Missense", "Nonsense", "Splicing", "Truncating"))
dnds_global$Group <- factor(dnds_global$Group, levels = GROUP_LEVELS)
save_pdf(sprintf("%s.pub.Fig2C.pdf", output_header), width = 6, height = 3, plots = list(
  plot_dnds_ratio(dnds_global[dnds_global$Type %in% c("All", "Missense", "Truncating"), ])
))

# Does global truncating dN/dS (e3/e3) differ between AD and CTRL? (heterogeneity of log ratios)
dnds_trunc <- dnds_global[dnds_global$Group != "C" & dnds_global$Type == "Truncating", ]
print(metagen(TE = log(dnds_trunc$MLE), seTE = (log(dnds_trunc$HighCI) - log(dnds_trunc$LowCI)) / (2 * 1.96),
              studlab = as.character(dnds_trunc$Group), sm = "ROM"))

# Fig 2D: per-gene missense / truncating dN/dS (e3/e3) for CHIP whitelist genes with Q < 0.05 in any group
sig_genes <- intersect(unique(as.character(dnds_gene$Gene[which(dnds_gene$Qvalue < 0.05)])), chip_genes)
if (length(sig_genes) > 0) {
  dnds_sig <- do.call(rbind, lapply(names(dnds_33), function(g) {
    ci <- geneci(dnds_33[[g]], gene_list = sig_genes)
    rbind(data.frame(Gene = ci$gene, MLE = ci$mis_mle, LowCI = ci$mis_low, HighCI = ci$mis_high, Group = g, Type = "Missense"),
          data.frame(Gene = ci$gene, MLE = ci$tru_mle, LowCI = ci$tru_low, HighCI = ci$tru_high, Group = g, Type = "Truncating"))
  }))
  dnds_sig$Type  <- factor(dnds_sig$Type, levels = c("Missense", "Truncating"))
  dnds_sig$Group <- factor(dnds_sig$Group, levels = GROUP_LEVELS)
  save_pdf(sprintf("%s.pub.Fig2D.pdf", output_header), width = 8, height = 3, plots = list(
    plot_dnds_ratio(dnds_sig) + scale_y_log10() + facet_grid(. ~ Gene)
  ))
} else {
  message("No CHIP whitelist gene has dNdScv Q < 0.05 in any group; Fig2D.pdf not written")
}

print("=== CHIP age-trend analysis finished ===")
