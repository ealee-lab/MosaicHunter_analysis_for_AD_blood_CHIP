# MosaicHunter_analysis_for_AD_blood_CHIP

Analysis code for the study:

> **Clonal Hematopoiesis Mutation Burden is Increased in Alzheimer's Disease Patients with *APOE* ε3/ε3 Genotype**
> Jaejoon Choi†, Kyung Sun Park†, Yann Le Guen†, Jong-Ho Park, Zinan Zhou, Liz Enyenihi, Ila Rosen, Yoon-Ho Choi, Christopher A. Walsh, Jong-Won Kim\*, August Yue Huang\*, Eunjung Alice Lee\*
> († equal contribution; \* corresponding authors)

An earlier version of this work is available as a preprint:
*Clonal Hematopoiesis Mutations Increase Risk of Alzheimer's Disease with APOE ε3/ε3 Genotype*, bioRxiv (2025), https://www.biorxiv.org/content/10.1101/2025.05.19.654981v1.
The code in this repository corresponds to the **revised manuscript**, which differs from the preprint in cohort size, analyses and figure numbering.

Repository: https://github.com/ealee-lab/MosaicHunter_analysis_for_AD_blood_CHIP

## Contents

- [Overview](#overview)
- [Repository contents](#repository-contents)
- [1. System requirements](#1-system-requirements)
- [2. Installation guide](#2-installation-guide)
- [3. Demo](#3-demo)
- [4. Instructions for use](#4-instructions-for-use)
- [Reproducing the manuscript results](#reproducing-the-manuscript-results)
- [License](#license)
- [Citation](#citation)

## Overview

We used molecular-barcoded deep targeted panel sequencing (~400X) of 149 cancer driver genes to profile somatic mutations in the blood of 353 individuals aged 50 years and older: 108 Alzheimer's disease (AD) patients, 196 controls and 49 centenarians without AD. The analysis focuses on mutations that drive clonal hematopoiesis of indeterminate potential (CHIP). Burden comparisons between AD and controls use a subset of 83 AD patients and 123 controls matched for age, sex and cardiovascular disease history.

The pipeline has four scripts, run in order:

1. **Scr1** filters somatic SNVs (MosaicHunter) and indels (Pisces) across the cohort and annotates them with ANNOVAR.
2. **Scr2** restricts calls to curated CHIP mutations, adds gene annotation, applies depth and gnomAD filters, and merges SNVs and indels.
3. **Scr3** compares somatic and CHIP mutation burden between AD patients and a propensity-score-matched control set (Fig. 1).
4. **Scr4** describes age trends, the relationship with cognitive score (MMSE), and selection (dN/dS) of somatic mutations (Fig. 2).

```
MosaicHunter / Pisces raw calls   (not distributed)
        │
        ▼
Scr1_variant_processing.sh              →  AD_panel.MH.loose.tsv, AD_panel.MH.loose.gene_anno.tsv,
        │                                   AD_panel.indel.loose.tsv, AD_panel.indel.loose.gene_anno.tsv,
        │                                   AD_panel.coverage.summary        (all provided in this repository)
        ▼
Scr2_CHIP_variant_annotation.R          →  AD_panel.merged.addage.tsv, whitelist_Genename.tsv
        │
        ├──▶ Scr3_Comparative_analysis_AD_CTRL.R  →  Fig. 1, ED Figs. 1, 2, 4
        └──▶ Scr4_Age_Trend_Analysis.R            →  Fig. 2, ED Figs. 7, 8, 9, 10
```

**Scr1 cannot be run from this repository.** It needs the per-sample MosaicHunter and Pisces call sets, which are not distributed. Its outputs are provided instead, so **Scr2, Scr3 and Scr4 can be run directly** on the files in this repository.

## Repository contents

### Scripts

| File | Language | Purpose |
|---|---|---|
| `Scr1_variant_processing.sh` | Bash | Sample coverage filtering, cohort-level recurrence checks, Fisher's-exact filtering of MosaicHunter SNVs and Pisces indels, ANNOVAR annotation |
| `fisher_filter_v2_2sided_0.1proportion.py` | Python | Fisher's-exact filter for SNVs (called by Scr1) |
| `fisher_filter_indel_v2_2sided_0.1proportion.py` | Python | Fisher's-exact filter for indels (called by Scr1) |
| `Scr2_CHIP_variant_annotation.R` | R | CHIP whitelist filtering, gene annotation, depth/gnomAD filtering, merging of SNVs and indels |
| `Scr3_Comparative_analysis_AD_CTRL.R` | R | Propensity-score matching of AD and controls; mutation burden, mutation landscape, logistic regression, VAF- and *APOE*-stratified analyses |
| `Scr4_Age_Trend_Analysis.R` | R | Age trends of CHIP carriers and burden, VAF vs. age, MMSE association in AD, dN/dS selection analysis |

### Data files

| File | Used by | Description |
|---|---|---|
| `AD_panel.MH.loose.tsv` | Scr2 | Filtered, ANNOVAR-annotated somatic SNV calls for all samples (Scr1 output) |
| `AD_panel.MH.loose.gene_anno.tsv` | Scr2 | Gene symbol, location, consequence and predicted impact for each SNV (Scr1 output) |
| `AD_panel.indel.loose.tsv` | Scr2 | Filtered, ANNOVAR-annotated somatic indel calls for all samples (Scr1 output) |
| `AD_panel.indel.loose.gene_anno.tsv` | Scr2 | Gene symbol, location, consequence and predicted impact for each indel (Scr1 output) |
| `AD_panel.coverage.summary` | Scr3, Scr4 | Per-sample mean depth and fraction of panel bases covered ≥ 200X (Scr1 output) |
| `whitelist_NatMed.xlsx` | Scr2 | CHIP mutations previously reported by Bouzid et al. (*Nat Med* 2023), used as the whitelist |
| `SupplTable1_Sample_Information.xlsx` | Scr2, Scr3, Scr4 | Clinical metadata of the 353 QC-passed samples (sheet `Sample_Information`) |
| `SupplTable2_Somatic_Variant_Calls.xlsx` | – | Final somatic variant calls (Supplementary Table 2) |
| `3277911_Covered.bed` | Scr3, Scr4 | Capture panel target regions |
| `CancerGeneCensus_01212022.tsv` | Scr4 | COSMIC Cancer Gene Census (21 Jan 2022); only used to label tumour suppressor genes in the dN/dS gene table |
| `AD_panel.merged.addage.tsv` | Scr3, Scr4 | Merged SNV + indel calls (Scr2 output), provided for reference; running Scr2 overwrites it |
| `AD_panel.sample.list`, `AD_panel.sample.clinical.tsv`, `summary_on_target_depth_Batch*.tsv` | Scr1 | Sample lists, batch-level clinical tables and per-sample on-target depth used for sample QC (samples below 100X are excluded) |

Sample groups, as defined in `SupplTable1_Sample_Information.xlsx`:

- **AD**: clinically diagnosed Alzheimer's disease (108)
- **CTRL**: individuals without AD aged 90 years or younger (196)
- **C**: individuals without AD older than 90 years (49; labelled "Centenarian" in figures)

Only individuals with *APOE* ε3/ε3 or ε3/ε4 genotypes are included.

## 1. System requirements

### Operating system

The R scripts (Scr2–Scr4) run on any operating system supported by R (Linux, macOS, Windows). Scr1 requires a Unix-like shell (Linux or macOS).

The demo (Scr2–Scr4) was tested on **Ubuntu 24.04 LTS** with **R 4.3.3**.

### Software dependencies

R ≥ 4.1 is required. The demo was tested with these package versions:

| Package | Tested version | Source | Used in |
|---|---|---|---|
| readxl | 1.4.3 | CRAN | Scr2, Scr3, Scr4 |
| dplyr | 1.1.4 | CRAN | Scr2, Scr3, Scr4 |
| tidyr | 1.3.1 | CRAN | Scr2, Scr3 |
| limma | 3.58.1 | Bioconductor | Scr2 |
| org.Hs.eg.db | 3.18.0 | Bioconductor | Scr2 |
| ggplot2 | 3.4.4 | CRAN | Scr3, Scr4 |
| ggpubr | 0.6.0 | CRAN | Scr3 |
| MatchIt | 4.5.5 | CRAN | Scr3 |
| meta | 7.0-0 | CRAN | Scr3, Scr4 |
| ggbreak | 0.1.2 | CRAN | Scr4 |
| ggrepel | 0.9.5 | CRAN | Scr4 |
| dndscv | 0.0.1.0 | [GitHub](https://github.com/im3sanger/dndscv) | Scr4 |

Scr1 additionally requires:

- Bash with standard Unix tools (`awk`, `grep`, `cut`, `paste`)
- Python 3 with SciPy, for the Fisher's-exact filters
- [ANNOVAR](https://annovar.openbioinformatics.org/) (version 2015-03-22) with the hg19 `refGene`, `gnomad_genome`, `ljb26_all` and `cosmic70` databases
- The table-joining utilities `my.grep` and `myjoin` on your `PATH`

### Hardware

No non-standard hardware is required. The demo runs on a standard desktop or laptop computer; it was tested with 1 CPU core and 4 GB of RAM.

## 2. Installation guide

1. Install [R](https://cran.r-project.org/) (≥ 4.1).

2. Clone the repository:

   ```bash
   git clone https://github.com/ealee-lab/MosaicHunter_analysis_for_AD_blood_CHIP.git
   cd MosaicHunter_analysis_for_AD_blood_CHIP
   ```

3. Install the R packages from within R:

   ```r
   install.packages(c("readxl", "dplyr", "tidyr", "ggplot2", "ggpubr", "MatchIt",
                      "meta", "ggbreak", "ggrepel", "remotes", "BiocManager"))
   BiocManager::install(c("limma", "org.Hs.eg.db"))
   remotes::install_github("im3sanger/dndscv")
   ```

The scripts themselves need no installation or compilation.

**Typical install time:** about 10 minutes on a normal desktop computer where CRAN provides precompiled packages (Windows, macOS). On Linux, where packages are compiled from source, allow 20–30 minutes.

## 3. Demo

The demo runs Scr2, Scr3 and Scr4 on the data files in this repository. It reproduces the panel-sequencing figures and statistics of the manuscript.

### Instructions

From the repository directory:

```bash
Rscript Scr2_CHIP_variant_annotation.R      > Scr2.log 2>&1
Rscript Scr3_Comparative_analysis_AD_CTRL.R > Scr3.log 2>&1
Rscript Scr4_Age_Trend_Analysis.R           > Scr4.log 2>&1
```

Run Scr2 before Scr3 and Scr4. It regenerates `AD_panel.merged.addage.tsv` from the Scr1 outputs and also writes `whitelist_Genename.tsv`, which Scr3 and Scr4 need.

Statistical results (proportion tests, regression summaries, meta-analyses) are printed to the log files. Figures are written as PDF files.

### Expected output

**Scr2** writes `AD_panel.merged.addage.tsv` (563 somatic variants in the 353 QC-passed samples), `whitelist_Genename.tsv`, and intermediate tables (`*.WLfiltered.tsv`, `*.WLannotated.tsv`, `*.merged.tsv`, `*.merged_filtered.tsv`).

**Scr3 and Scr4** write 20 PDF figures and one table, all prefixed with `AD_panel.merged.addage.pub.clinFilt` (listed under [Reproducing the manuscript results](#reproducing-the-manuscript-results)).

Key statistics in `Scr3.log` should match the manuscript:

| Result (AD vs. matched controls) | Expected value |
|---|---|
| Matched set | 83 AD patients, 123 controls (`Published` / `Recomputed` table: 147 FALSE, 206 TRUE) |
| All somatic mutations | p = 0.0154, rate ratio 1.27 |
| CHIP mutations | p = 0.00259, rate ratio 2.06 |
| CHIP mutations, *APOE* ε3/ε3 | p = 3.96 × 10⁻⁵, rate ratio 3.17 |
| CHIP mutations, *APOE* ε3/ε4 | p = 0.627, rate ratio 0.76 |
| CHIP mutations, VAF ≤ 0.02, *APOE* ε3/ε3 | p = 0.0113, rate ratio 2.48 |
| CHIP mutations, 0.02 < VAF ≤ 0.05, *APOE* ε3/ε3 | p = 0.0255, rate ratio 3.77 |
| Logistic regression, CHIP burden | p = 0.0218, OR 1.44 |
| Logistic regression, *APOE* ε3/ε4 | p = 0.0132, OR 2.18 |

The proportion tests are one-sided (AD > control). The rate ratio is the AD/control ratio of mutations per sequenced base pair, reported as the odds ratio in the manuscript.

### Expected run time

Measured on the test system (1 CPU core, 4 GB RAM): about 9 s for Scr2, 11 s for Scr3 and 55 s for Scr4. Most of Scr4's time is spent in the dN/dS analysis. The full demo takes **under 2 minutes** on a normal desktop computer.

## 4. Instructions for use

### Running on your own data

The scripts read fixed file names from the working directory. Input paths and analysis thresholds are set in the **Configuration** block at the top of Scr3 and Scr4.

To analyse a new cohort, prepare files in the same formats as the provided ones (tab-separated, no header, hg19 coordinates):

- **`AD_panel.MH.loose.tsv`** (SNVs, 14 columns): chromosome, position, reference allele and alternative allele; the MosaicHunter read-count fields; sample ID; then the ANNOVAR location, gene, exonic consequence and amino-acid change. Use the provided file as a template.
- **`AD_panel.indel.loose.tsv`** (indels, 12 columns): chromosome, start, end, reference allele and alternative allele; the Pisces genotype field (`GT:GQ:AD:DP:VF:NL:SB`); sample ID; then the ANNOVAR location, gene, exonic consequence and amino-acid change. Use the provided file as a template.
- **`*.gene_anno.tsv`**: variant coordinates and alleles, then gene, location, consequence and `deleterious`/`neutral`.
- **`AD_panel.coverage.summary`**: sample ID, mean depth, fraction of panel bases covered ≥ 200X.
- **`SupplTable1_Sample_Information.xlsx`**: one row per sample, with the columns `ID`, `APOE genotype`, `Age (at blood collection)`, `Age (at diagnosis)`, `Sex`, `MMSE`, `Group` (`AD`, `CTRL` or `C`), `Batch`, `Dementia type`, `CVD history` and `Selected for subset`.
- **`3277911_Covered.bed`**: your capture panel's target regions, with gene symbols in the fourth column.

Then run Scr2, Scr3 and Scr4 as in the demo. If your matched subset differs from the `Selected for subset` column, Scr3 prints a warning; this is expected for a new cohort.

To start from raw MosaicHunter and Pisces calls, run Scr1 first. It expects:

- `MosaicHunter/<sample>.MH.tsv` and `Pisces/<sample>.final.vcf` for each sample
- `finalBam/<sample>.final.coverageBed`
- the sample lists and depth summaries listed under [Data files](#data-files)

It also reads two environment variables:

```bash
export REFERENCE_DIR=/path/to/reference
export ANNOVAR_DIR=/path/to/annovar/humandb
bash Scr1_variant_processing.sh
```

## Reproducing the manuscript results

Running the [demo](#3-demo) reproduces the following figures and statistics of the revised manuscript.

### Scr3: CHIP burden in AD vs. matched controls

AD patients and controls aged ≥ 60 are matched on age, sex and cardiovascular disease history (MatchIt, nearest neighbour, up to 2 controls per case, caliper 0.1). The script checks that the recomputed matched set reproduces the `Selected for subset` column of Supplementary Table 1.

| Output | Content |
|---|---|
| `.pub.EDFig2.pdf` | Sequencing depth and coverage across all samples |
| `.pub.Fig1A.pdf` | Age distribution before and after matching |
| `.pub.EDFig1.pdf`, `.pub.EDFig1.add.pdf` | Age distribution per group; sex and *APOE* composition |
| `.pub.Fig1B.pdf`, `.pub.Fig1C.pdf` | Rate of all somatic mutations and of CHIP mutations per bp |
| `.pub.Fig1D.pdf` | CHIP mutation landscape in the matched set, including carriers of multiple CHIP genes |
| `.pub.EDFig4.pdf` | Somatic mutation landscape in all samples |
| `.pub.Fig1E.pdf` | Logistic regression of AD status on CHIP burden, adjusted for sex, age and *APOE* genotype (page 2 also adjusts for cardiovascular disease history) |
| `.pub.Fig1F.pdf` | Mutation rates stratified by *APOE* genotype (ε3/ε3 vs. ε3/ε4) |
| `.pub.Fig1G.pdf` | Odds ratio as a function of minimum VAF (matched set, and its *APOE* ε3/ε3 subset) |
| `.pub.Fig1H.pdf` | CHIP mutation rate per VAF bin (with Mantel–Haenszel meta-analysis across bins) |

### Scr4: Age trends, MMSE and selection

These analyses use all 353 QC-passed samples (no matching).

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

## License

This code is released under the [MIT License](LICENSE).

The COSMIC Cancer Gene Census (`CancerGeneCensus_01212022.tsv`) and the CHIP whitelist (`whitelist_NatMed.xlsx`) are third-party data and remain subject to their original terms of use.

## Citation

If you use this code, please cite the paper. Until the revised version is published, please cite the preprint:

```
Choi J, Park KS, Le Guen Y, Park J-H, Zhou Z, Enyenihi L, Rosen I, Choi Y-H, Walsh CA,
Kim J-W, Huang AY, Lee EA. Clonal Hematopoiesis Mutations Increase Risk of Alzheimer's
Disease with APOE ε3/ε3 Genotype. bioRxiv (2025). doi:10.1101/2025.05.19.654981
```
