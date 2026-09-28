# MosaicHunter_analysis_for_AD_blood_CHIP

Analysis code for the study:

> **Clonal Hematopoiesis Mutations are Associated with Alzheimer's Disease Stratified by *APOE* ε3/ε3 Genotype**
> Jaejoon Choi†, Kyung Sun Park†, Yann Le Guen†, Jong-Ho Park, Zinan Zhou, Liz Enyenihi, Ila Rosen, Yoon-Ho Choi, Christopher A. Walsh, Jong-Won Kim\*, August Yue Huang\*, Eunjung Alice Lee\*
> († equal contribution; \* corresponding authors)

An earlier version of this work is available as a preprint:
*Clonal Hematopoiesis Mutations Increase Risk of Alzheimer's Disease with APOE ε3/ε3 Genotype*, bioRxiv (2025), https://www.biorxiv.org/content/10.1101/2025.05.19.654981v1.
The code in this repository corresponds to the **revised manuscript**, which differs from the preprint in cohort size, analyses and figure numbering.

## Overview

We used molecular-barcoded deep targeted panel sequencing (~400X) of 149 cancer driver genes to profile somatic mutations in the blood of 353 individuals aged 50 years and older: 108 Alzheimer's disease (AD) patients, 196 controls and 49 centenarians without AD. The analysis focuses on mutations that drive clonal hematopoiesis of indeterminate potential (CHIP). Burden comparisons between AD and controls use a subset of 83 AD patients and 123 controls matched for age, sex and cardiovascular disease history.

This repository contains the scripts used to go from per-sample variant calls to the figures and statistics reported in the paper:

1. Filtering somatic SNVs (MosaicHunter) and indels (Pisces) and annotating them with ANNOVAR
2. Restricting to curated CHIP mutations and merging SNVs and indels into one call set
3. Comparing somatic and CHIP mutation burden between AD patients and a propensity-score-matched control set
4. Describing age trends, the relationship with cognitive score (MMSE), and selection (dN/dS) of somatic mutations

## Repository contents

| Script | Language | Purpose |
|---|---|---|
| `Scr1_variant_processing.sh` | Bash | Sample coverage filtering, cohort-level recurrence checks and Fisher's-exact filtering of MosaicHunter SNVs and Pisces indels, ANNOVAR annotation |
| `Scr2_CHIP_variant_annotation.R` | R | CHIP whitelist filtering, gene annotation, depth/gnomAD filtering, merging of SNVs and indels |
| `Scr3_CHIP_burden_analysis.R` | R | Propensity-score matching of AD and controls; mutation burden, mutation landscape, logistic regression, VAF- and *APOE*-stratified analyses (Fig. 1) |
| `Scr4_CHIP_age_trend_analysis.R` | R | Age trends of CHIP carriers and burden, VAF vs. age, MMSE association in AD, dN/dS selection analysis (Fig. 2) |

The scripts are meant to be run in order, from a single project directory:

```
MosaicHunter / Pisces calls
        │
        ▼
Scr1_variant_processing.sh      →  AD_panel.MH.loose.*, AD_panel.indel.loose.*, AD_panel.coverage.summary
        │
        ▼
Scr2_CHIP_variant_annotation.R  →  AD_panel.merged.addage.tsv, whitelist_Genename.tsv
        │
        ├──▶ Scr3_CHIP_burden_analysis.R     →  Fig. 1, ED Figs. 1, 2, 4
        └──▶ Scr4_CHIP_age_trend_analysis.R  →  Fig. 2, ED Figs. 8, 9, 10
```

## Requirements

### Command-line tools (Scr1)

- [ANNOVAR](https://annovar.openbioinformatics.org/) with the hg19 `refGene`, `gnomad_genome`, `ljb26_all` and `cosmic70` databases
- Python 3, for the Fisher's-exact filtering helpers `fisher_filter_v2_2sided_0.1proportion.py` and `fisher_filter_indel_v2_2sided_0.1proportion.py`
- The table-joining helpers `my.grep` and `myjoin`, on your `PATH`

Scr1 reads two environment variables:

```bash
export REFERENCE_DIR=/path/to/reference
export ANNOVAR_DIR=/path/to/annovar/humandb
```

### R packages (R ≥ 4.1)

| Package | Source | Used in |
|---|---|---|
| readxl, dplyr, tidyr | CRAN | Scr2, Scr3, Scr4 |
| limma, org.Hs.eg.db | Bioconductor | Scr2 |
| ggplot2, ggpubr, MatchIt, meta | CRAN | Scr3 |
| ggplot2, ggbreak, ggrepel, meta | CRAN | Scr4 |
| [dndscv](https://github.com/im3sanger/dndscv) | GitHub | Scr4 |

```r
install.packages(c("readxl", "dplyr", "tidyr", "ggplot2", "ggpubr", "MatchIt",
                   "meta", "ggbreak", "ggrepel", "remotes", "BiocManager"))
BiocManager::install(c("limma", "org.Hs.eg.db"))
remotes::install_github("im3sanger/dndscv")
```

## Input files

| File | Used by | Description |
|---|---|---|
| `MosaicHunter/<sample>.MH.tsv` | Scr1 | Per-sample MosaicHunter SNV calls |
| `Pisces/<sample>.final.vcf` | Scr1 | Per-sample Pisces indel calls |
| `finalBam/<sample>.final.coverageBed` | Scr1 | Per-sample coverage histograms over the panel |
| `summary_on_target_depth_Batch*.tsv` | Scr1 | Mean on-target depth per sample and batch (samples below 100X are excluded) |
| `whitelist_NatMed.xlsx` | Scr2 | CHIP mutations previously reported by Bouzid et al. (*Nat Med* 2023), used as the whitelist |
| `SupplTable1_Sample_Information.xlsx` | Scr3, Scr4 | Clinical metadata of the 353 QC-passed samples (sheet `Sample_Information`) |
| `SupplTable2_Somatic_Variant_Calls.xlsx` | Scr2 | Somatic variant calls with sample attributes (sheet `Somatic_Variant_Calls`) |
| `3277911_Covered.bed` | Scr3, Scr4 | Capture panel target regions |
| `CancerGeneCensus_01212022.tsv` | Scr4 (optional) | COSMIC Cancer Gene Census; only used to label tumour suppressor genes in the dN/dS gene table. Not distributed here because of COSMIC licensing. |

Sample groups, as defined in `SupplTable1_Sample_Information.xlsx`:

- **AD**: clinically diagnosed Alzheimer's disease
- **C**: individuals without AD who are older than 90 years (49; labelled "Centenarian" in figures)
- **CTRL**: all other individuals without AD (196)

AD has 108 individuals.

Only individuals with *APOE* ε3/ε3 or ε3/ε4 genotypes are included.

## Usage

```bash
bash    Scr1_variant_processing.sh
Rscript Scr2_CHIP_variant_annotation.R
Rscript Scr3_CHIP_burden_analysis.R
Rscript Scr4_CHIP_age_trend_analysis.R
```

Thresholds and file paths are set in the **Configuration** block at the top of each R script. Statistical test results (proportion tests, regression summaries, meta-analyses) are printed to standard output, so keep a log of each run:

```bash
Rscript Scr3_CHIP_burden_analysis.R > Scr3.log 2>&1
```

## Outputs

All figure files from Scr3 and Scr4 are prefixed with `AD_panel.merged.addage.pub.clinFilt`.

### Scr3: CHIP burden in AD vs. matched controls

AD patients and controls aged ≥ 60 are matched on age, sex and cardiovascular disease history (MatchIt, nearest neighbour, up to 2 controls per case, caliper 0.1), giving 83 AD patients and 123 controls. The script checks that the recomputed matched set reproduces the `Selected for subset` column of Supplementary Table 1.

| Output | Content |
|---|---|
| `.pub.EDFig2.pdf` | Sequencing depth and coverage across all samples |
| `.pub.Fig1A.pdf` | Age distribution before and after matching (study design) |
| `.pub.EDFig1.pdf`, `.pub.EDFig1.add.pdf` | Age distribution per group; sex and *APOE* composition |
| `.pub.Fig1B.pdf`, `.pub.Fig1C.pdf` | Rate of all somatic mutations and of CHIP mutations per bp |
| `.pub.Fig1D.pdf` | CHIP mutation landscape in the matched set, including carriers of multiple CHIP genes |
| `.pub.EDFig4.pdf` | Somatic mutation landscape in all samples |
| `.pub.Fig1E.pdf` | Logistic regression of AD status on CHIP burden, adjusted for sex, age and *APOE* genotype (page 2 also adjusts for cardiovascular disease history) |
| `.pub.Fig1F.pdf` | Mutation rates stratified by *APOE* genotype (ε3/ε3 vs. ε3/ε4) |
| `.pub.Fig1G.pdf` | Odds ratio as a function of minimum VAF (matched set, and its *APOE* ε3/ε3 subset) |
| `.pub.Fig1H.pdf` | CHIP mutation rate per VAF bin (with Mantel–Haenszel meta-analysis across bins) |

### Scr4: Age trends, MMSE and selection

These analyses use all QC-passed samples (no matching).

| Output | Content |
|---|---|
| `.pub.Fig1Aa.pdf` | Age histogram per group |
| `.pub.Fig2A.pdf` | Proportion of CHIP mutation carriers by age: all samples (ED Fig. 7A) and *APOE* ε3/ε3 only (Fig. 2A) |
| `.pub.Fig2B.pdf` | CHIP mutation rate per bp by age group: all samples (ED Fig. 7B) and *APOE* ε3/ε3 only (Fig. 2B) |
| `.pub.EDFig8.pdf` | VAF (mutant clone size) of somatic mutations vs. age (ED Fig. 8A) |
| `.pub.EDFig10.pdf` | CHIP carriers and burden by MMSE and *APOE* in AD; regression of CHIP carrier status on MMSE |
| `.pub.EDFig9.pdf` | Per-gene dN/dS significance by group (*APOE* ε3/ε3 SNVs) |
| `.pub.Fig2C.pdf` | Global dN/dS by group (*APOE* ε3/ε3 SNVs) |
| `.pub.Fig2D.pdf` | Per-gene dN/dS for significantly selected CHIP genes (*APOE* ε3/ε3 SNVs) |
| `.dNdScv.gene.tsv` | Gene-level dN/dS results for genes with P < 0.05 in any group |

Figure numbers follow the revised manuscript and differ from those in the preprint. The following panels are not produced by these scripts: amplicon validation (ED Fig. 3), the ADSP replication analysis (Fig. 1I, ED Fig. 6) and the comparison with standard-depth sequencing (ED Fig. 8B).

## Citation

If you use this code, please cite the paper. Until the revised version is published, please cite the preprint:

```
Choi J, Park KS, Le Guen Y, Park J-H, Zhou Z, Enyenihi L, Rosen I, Choi Y-H, Walsh CA,
Kim J-W, Huang AY, Lee EA. Clonal Hematopoiesis Mutations Increase Risk of Alzheimer's
Disease with APOE ε3/ε3 Genotype. bioRxiv (2025). doi:10.1101/2025.05.19.654981
```
