# ==============================================================================
# Pipeline: CHIP Variant Whitelist Filtering, Gene Annotation, and Merging
# Description: Merged R script that performs SNV and indel whitelist filtering,
#              gene annotation processing, artifact removal, clinical metadata
#              integration (via SupplTable1 clinical information), and final combined dataset generation.
# ==============================================================================

# Load necessary libraries for Excel reading, table joining, and bioinformatics
library(readxl)
library(limma)
library(dplyr)
library(tidyr)
library(org.Hs.eg.db)

# ==============================================================================
# Helper Function: Load Clinical Metadata from SupplTable1
# ==============================================================================
load_clinical_data <- function() {
  suppl_path <- "SupplTable1_Sample_Information.xlsx"
  # Read the sample information sheet which contains clinical attributes
  raw_data <- read_excel(suppl_path, sheet = "Sample_Information")
  
  # Extract and standardize clinical attributes per sample ID
  raw_data %>% 
    select(ID, `APOE genotype`, `Age (at blood collection)`, Sex, MMSE, Group, Batch) %>% 
    distinct() %>%
    rename(
      ApoE = `APOE genotype`,
      Age = `Age (at blood collection)`
    ) %>%
    mutate(
      Sex = factor(Sex, levels = c("Female", "Male", "F", "M")),
      ApoE = factor(ApoE, levels = sort(unique(ApoE))),
      Cogdx = factor(Group, levels = c("CTRL", "C", "AD"))
    ) %>%
    arrange(ID)
}

print("=== Step 1: SNV Whitelist Filtering ===")
snv_header <- "AD_panel.MH.loose"

# Load literature-curated CHIP mutation whitelist from Excel
df_snv_wl <- data.frame(read_excel("whitelist_NatMed.xlsx", sheet = 1, skip = 3))
write.table(df_snv_wl$Gene.name, "whitelist_Genename.tsv", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
df_snv_wl$nonsense <- grepl('nonsense', df_snv_wl[, 2], fixed = TRUE)
df_snv_wl$splice <- grepl('splice', df_snv_wl[, 2], fixed = TRUE)
df_snv_wl$missense <- grepl('missense', df_snv_wl[, 2], fixed = TRUE)

# Parse protein coordinate boundaries for whitelisted SNVs
df_snv_wl$start <- NA
df_snv_wl$start[grepl("p\\.[0-9]+", df_snv_wl[, 2])] <- as.numeric(gsub("-[0-9]+.*", "", gsub(".*p\\.", "", df_snv_wl[grepl("p\\.[0-9]+", df_snv_wl[, 2]), 2])))
df_snv_wl$end <- NA
df_snv_wl$end[grepl("p\\.[0-9]+", df_snv_wl[, 2])] <- as.numeric(gsub(")\\..*", "", gsub(".*p\\..+-", "", df_snv_wl[grepl("p\\.[0-9]+", df_snv_wl[, 2]), 2])))
df_snv_wl$ex_start <- NA
df_snv_wl$ex_start[df_snv_wl$Gene.name == "TET2"] <- 1482
df_snv_wl$ex_end <- NA
df_snv_wl$ex_end[df_snv_wl$Gene.name == "TET2"] <- 1842

# Handle gene aliases and formatting
df_snv_wl[, 2] <- gsub("\\([0-9]+\\)", "", df_snv_wl[, 2])
df_snv_wl[, 1] <- alias2SymbolTable(df_snv_wl[, 1], species = "Hs")

# Load raw SNV TSV file and filter for exonic/splicing non-synonymous variants
df_snv_in <- read.table(paste(snv_header, ".tsv", sep = ""), fill = TRUE, sep = "\t")
df_snv_ori <- df_snv_in
df_snv_in <- df_snv_in[df_snv_in[, 11] == "exonic" | df_snv_in[, 11] == "splicing", ]
df_snv_in <- df_snv_in[df_snv_in[, 13] != "synonymous" & df_snv_in[, 13] != "unknown", ]

# Iteratively evaluate each SNV against the whitelist criteria
df_snv_out <- NULL
for (i in 1:nrow(df_snv_in)) {
  df_tmp <- df_snv_wl[df_snv_wl[, 1] == gsub("\\(.*", "", df_snv_in[i, 12]), ]
  if (nrow(df_tmp) == 0) next
  if (df_snv_in[i, 11] == "splicing") {
    entry_list <- strsplit(as.character(df_snv_in[i, 12]), split = ',')[[1]]
    entry <- entry_list[grep(df_tmp$Accession, entry_list, fixed = TRUE)]
    if (length(entry) == 0) next
    exon <- gsub(":.*", "", gsub(".*exon", "", entry))
    pos <- as.numeric(gsub("[^0-9].*", "", gsub(".*:c\\.", "", entry))) / 3
  } else {
    entry_list <- strsplit(as.character(df_snv_in[i, 14]), split = ',')[[1]]
    entry <- entry_list[grep(df_tmp$Accession, entry_list, fixed = TRUE)]
    if (length(entry) == 0) next
    exon <- gsub(":.*", "", gsub(".*exon", "", entry))
    pos <- as.numeric(gsub("[^0-9]", "", gsub(".*:p\\..", "", entry)))
    var <- gsub(".*:p\\.", "", entry)
  }
  
  if (!is.na(df_tmp$ex_start)) {
    if (df_snv_in[i, 13] == "nonsynonymous") {
      if (pos >= df_tmp$start & pos <= df_tmp$ex_start) { df_snv_out <- rbind(df_snv_out, df_snv_in[i, ]); next }
      if (pos >= df_tmp$ex_end & pos <= df_tmp$end) { df_snv_out <- rbind(df_snv_out, df_snv_in[i, ]); next }
    }
  } else if (!is.na(df_tmp$start)) {
    if (df_tmp$missense == TRUE && df_snv_in[i, 13] == "nonsynonymous") {
      if (pos >= df_tmp$start & pos <= df_tmp$end) { df_snv_out <- rbind(df_snv_out, df_snv_in[i, ]); next }
    } else if ((df_tmp$nonsense == TRUE | df_tmp$splice == TRUE) && (df_snv_in[i, 11] == "splicing" | df_snv_in[i, 13] == "stopgain")) {
      if (pos >= df_tmp$start & pos <= df_tmp$end) { df_snv_out <- rbind(df_snv_out, df_snv_in[i, ]); next }
    }
    next
  }
  if (df_snv_in[i, 11] == "splicing" && df_tmp$splice == TRUE) { df_snv_out <- rbind(df_snv_out, df_snv_in[i, ]); next }
  if (df_snv_in[i, 13] == "stopgain" && df_tmp$nonsense == TRUE) {
    df_snv_out <- rbind(df_snv_out, df_snv_in[i, ])
    next
  } else if (df_snv_in[i, 13] == "nonsynonymous") {
    if (grepl(var, df_tmp[2])) { df_snv_out <- rbind(df_snv_out, df_snv_in[i, ]); next }
    if (df_snv_in[i, 12] == "PPM1D" && (exon == 5 | exon == 6)) { df_snv_out <- rbind(df_snv_out, df_snv_in[i, ]); next }
    else if (df_snv_in[i, 12] == "BRAF" && grepl("R603", var)) { df_snv_out <- rbind(df_snv_out, df_snv_in[i, ]); next }
  }
}

# Write out filtered and annotated SNV whitelist tables
write.table(df_snv_out, paste(snv_header, ".WLfiltered.tsv", sep = ""), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
df_snv_out$whitelist <- "whitelist"
df_snv_out2 <- merge(df_snv_ori, merge(df_snv_ori, df_snv_out, all.x = TRUE, sort = FALSE), sort = FALSE)
write.table(df_snv_out2, paste(snv_header, ".WLannotated.tsv", sep = ""), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)


print("=== Step 2: Indel Whitelist Filtering ===")
indel_header <- "AD_panel.indel.loose"

# Load whitelist specifically for indels (frameshift / splice)
df_indel_wl <- data.frame(read_excel("whitelist_NatMed.xlsx", sheet = 1, skip = 3))
df_indel_wl$splice <- grepl('splice', df_indel_wl[, 2], fixed = TRUE)
df_indel_wl$frameshift <- grepl('Frameshift', df_indel_wl[, 2], fixed = TRUE)

df_indel_wl$start <- NA
df_indel_wl$start[grepl("p\\.[0-9]+", df_indel_wl[, 2])] <- as.numeric(gsub("-[0-9]+.*", "", gsub(".*p\\.", "", df_indel_wl[grepl("p\\.[0-9]+", df_indel_wl[, 2]), 2])))
df_indel_wl$end <- NA
df_indel_wl$end[grepl("p\\.[0-9]+", df_indel_wl[, 2])] <- as.numeric(gsub(")\\..*", "", gsub(".*p\\..+-", "", df_indel_wl[grepl("p\\.[0-9]+", df_indel_wl[, 2]), 2])))

df_indel_wl[, 2] <- gsub("\\([0-9]+\\)", "", df_indel_wl[, 2])
df_indel_wl[df_indel_wl$Gene.name == "TET2", "start"] <- NA
df_indel_wl[df_indel_wl$Gene.name == "TET2", "end"] <- NA
df_indel_wl[, 1] <- alias2SymbolTable(df_indel_wl[, 1], species = "Hs")

# Load raw indel TSV file
df_indel_in <- read.table(paste(indel_header, ".tsv", sep = ""), fill = TRUE, sep = "\t")
df_indel_ori <- df_indel_in
df_indel_in <- df_indel_in[df_indel_in[, 9] == "exonic" | df_indel_in[, 9] == "splicing", ]

# Iteratively evaluate each indel against the whitelist criteria
df_indel_out <- NULL
for (i in 1:nrow(df_indel_in)) {
  df_tmp <- df_indel_wl[df_indel_wl[, 1] == gsub("\\(.*", "", df_indel_in[i, 10]), ]
  if (nrow(df_tmp) == 0) next
  if (df_indel_in[i, 9] == "splicing") {
    entry_list <- strsplit(as.character(df_indel_in[i, 10]), split = ',')[[1]]
    entry <- entry_list[grep(df_tmp$Accession, entry_list, fixed = TRUE)]
    if (length(entry) == 0) next
    exon <- gsub(":.*", "", gsub(".*exon", "", entry))
    pos_start <- as.numeric(gsub("[^0-9].*", "", gsub(".*:c\\.", "", entry))) / 3
    pos_end <- pos_start + 1
  } else {
    entry_list <- strsplit(as.character(df_indel_in[i, 12]), split = ',')[[1]]
    entry <- entry_list[grep(df_tmp$Accession, entry_list, fixed = TRUE)]
    if (length(entry) == 0) next
    exon <- gsub(":.*", "", gsub(".*exon", "", entry))
    pos_start <- as.numeric(gsub("[^0-9].*", "", gsub(".*:p\\.[A-Z]*", "", entry)))
    pos_end <- as.numeric(gsub("[^0-9]*", "", gsub(".*_", "", gsub(".*:p\\.[A-Z]*", "", entry))))
    var <- gsub(".*:p\\.", "", entry)
  }
  
  if (!is.na(df_tmp$start)) {
    if (df_tmp$frameshift == TRUE && (pos_start >= df_tmp$start & pos_end <= df_tmp$end)) {
      df_indel_out <- rbind(df_indel_out, df_indel_in[i, ])
    }
    next
  }
  if (df_indel_in[i, 9] == "splicing" && df_tmp$splice == TRUE) {
    df_indel_out <- rbind(df_indel_out, df_indel_in[i, ])
  } else if ((df_indel_in[i, 11] == "frameshift insertion" | df_indel_in[i, 11] == "frameshift deletion") && df_tmp$frameshift == TRUE) {
    df_indel_out <- rbind(df_indel_out, df_indel_in[i, ])
  }
}

# Write out filtered and annotated indel whitelist tables
write.table(df_indel_out, paste(indel_header, ".WLfiltered.tsv", sep = ""), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
df_indel_out$whitelist <- "whitelist"
df_indel_out2 <- merge(df_indel_ori, merge(df_indel_ori, df_indel_out, all.x = TRUE, sort = FALSE), sort = FALSE)
write.table(df_indel_out2, paste(indel_header, ".WLannotated.tsv", sep = ""), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)


print("=== Step 3: SNV Annotation & Filtering ===")
snv_out_header <- "AD_panel.MH.rev.addage.loose"
input_snv <- read.delim(sprintf("%s.WLannotated.tsv", snv_header), header = FALSE, stringsAsFactors = FALSE)[, c(1:4, 6:7, 9, 10, 15)]
colnames(input_snv) <- c("Chr", "Pos", "Ref", "Alt", "Ref_num", "Alt_num", "ID", "Gnomad", "CHIP")
input_snv$MAF <- input_snv$Alt_num / (input_snv$Ref_num + input_snv$Alt_num)
input_snv$CHIP <- sapply(input_snv$CHIP, function(x) isTRUE(x == "whitelist"))

gene_anno_snv <- read.delim(sprintf("%s.gene_anno.tsv", snv_header), header = FALSE, stringsAsFactors = FALSE)
colnames(gene_anno_snv) <- c("Chr", "Pos", "Ref", "Alt", "Symbol", "Location", "Type", "Impact")

# Load clinical metadata from SupplTable1 and merge with SNVs
id_clinical <- load_clinical_data()
merged_snv <- unique(merge(merge(gene_anno_snv, input_snv, by = c("Chr", "Pos", "Ref", "Alt")), id_clinical, by = "ID"))
write.table(merged_snv[, which(colnames(merged_snv) != "Location")], file = sprintf("%s.merged.tsv", snv_out_header), quote = FALSE, sep = "\t", row.names = FALSE)

# Filter SNVs by read depth and gnomAD population frequency thresholds
merged_snv_filtered <- merged_snv[merged_snv$Ref_num + merged_snv$Alt_num >= 100 & merged_snv$Gnomad < 1e-4, ]
write.table(merged_snv_filtered[, which(colnames(merged_snv_filtered) != "Location")], file = sprintf("%s.merged_filtered.tsv", snv_out_header), quote = FALSE, sep = "\t", row.names = FALSE)
write.table(id_clinical, file = sprintf("%s.samples_filtered.tsv", snv_out_header), quote = FALSE, sep = "\t", row.names = FALSE)


print("=== Step 4: Indel Annotation & Filtering ===")
indel_out_header <- "AD_panel.indel.addage.loose"
input_indel <- read.delim(sprintf("%s.WLannotated.tsv", indel_header), header = FALSE, stringsAsFactors = FALSE)[, c(1:8, 13)]
colnames(input_indel) <- c("Chr", "Start", "End", "Ref", "Alt", "Info", "ID", "Gnomad", "CHIP")
input_indel$Ref_num <- sapply(input_indel$Info, function(x) as.numeric(strsplit(strsplit(x, split = ':')[[1]][3], split = ',')[[1]][1]))
input_indel$Alt_num <- sapply(input_indel$Info, function(x) as.numeric(strsplit(strsplit(x, split = ':')[[1]][3], split = ',')[[1]][2]))
input_indel$MAF <- input_indel$Alt_num / (input_indel$Ref_num + input_indel$Alt_num)
input_indel$CHIP <- sapply(input_indel$CHIP, function(x) isTRUE(x == "whitelist"))

gene_anno_indel <- read.delim(sprintf("%s.gene_anno.tsv", indel_header), header = FALSE, stringsAsFactors = FALSE)
colnames(gene_anno_indel) <- c("Chr", "Start", "End", "Ref", "Alt", "Symbol", "Location", "Type", "Impact")

# Merge indel annotations with clinical metadata from SupplTable1
merged_indel <- unique(merge(merge(gene_anno_indel, input_indel, by = c("Chr", "Start", "End", "Ref", "Alt")), id_clinical, by = "ID"))
write.table(merged_indel, file = sprintf("%s.merged.tsv", indel_out_header), quote = FALSE, sep = "\t", row.names = FALSE)

# Filter indels by depth and gnomAD frequency
merged_indel_filtered <- merged_indel[merged_indel$Ref_num + merged_indel$Alt_num >= 100 & merged_indel$Gnomad < 1e-4, ]
write.table(merged_indel_filtered, file = sprintf("%s.merged_filtered.tsv", indel_out_header), quote = FALSE, sep = "\t", row.names = FALSE)


print("=== Step 5: Merge Filtered SNVs and Indels ===")
df_snv_filt <- read.table(sprintf("%s.merged_filtered.tsv", snv_out_header), sep = "\t", header = TRUE)
df_indel_filt <- read.table(sprintf("%s.merged_filtered.tsv", indel_out_header), sep = "\t", header = TRUE)

# Align columns and combine SNVs and indels into a unified final dataset
df_tmp1 <- df_snv_filt[, c(1, 2, 3, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 21, 20)]
df_tmp2 <- df_indel_filt[, c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 14, 15, 12, 13, 16, 17, 18, 19, 20, 21, 23, 22)]
colnames(df_tmp1) <- colnames(df_tmp2)

df_final_merged <- rbind(df_tmp1, df_tmp2)
write.table(df_final_merged, "AD_panel.merged.addage.tsv", row.names = FALSE, sep = "\t", quote = FALSE)
