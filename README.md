# Henry et al. — Anole cloacal microbiome analysis scripts

Scripts and metadata for:

**Henry CJ, Stelling EG, Lattin CR** (in review). Effects of captivity and host species on the cloacal microbiomes of a native and an invasive anole. *Animal Microbiome*.

## Overview

This repository contains the QIIME2 shell script, R scripts, and sample metadata used to analyze 16S rRNA amplicon sequencing data from cloacal lavage samples of green anoles (*Anolis carolinensis*) and brown anoles (*Anolis sagrei*) collected at Louisiana State University, Baton Rouge, Louisiana, USA. Samples were collected at capture (pre-captivity) and again after two months in captivity (post-captivity).

**Target region:** 16S rRNA V4 (515F–806R)  
**Pipeline:** QIIME2 amplicon 2026.1  
**Sequencing platform:** AVITI (Element Biosciences), paired-end, at MSU Research Technology Support Facility Genomics Core

Raw sequencing reads are deposited in NCBI SRA under BioProject PRJNA1478682.

---

## Repository contents

```
qiime2_pipeline.sh                  Full QIIME2 pipeline from import to export
R/
    alpha_models_0p1.R              Alpha diversity linear mixed models (lmerTest)
    beta_permanova_full.R           PERMANOVA for all four beta diversity metrics
    check_sequencing_depth.R        Tests whether read depth differs by species/timepoint
    clr_mixed_model_taxa.R          CLR-transformed mixed models for Helicobacteraceae
                                    and Mycoplasmataceae
    genus_within_target_families.R  Genus-level breakdown of target families
    plot_alpha_paired_0p1.R         Alpha diversity paired plots (Fig. 1)
    plot_pcoa_all_metrics.R         PCoA plots for all four beta diversity metrics
    taxonomy_multilevel.R           Phylum/order/family bar charts and summary tables
metadata/
    sample-metadata-with-subject.tsv  Sample metadata including Subject column for
                                      paired mixed models (R analyses)
```

---

## Dependencies

### QIIME2
- QIIME2 amplicon 2026.1 (`conda activate qiime2-amplicon-2026.1`)
- SILVA 138.2 SSURef NR99 classifier (10.82364/138.2/SSU/Ref-NR99/QIIME2/2025.7/V4-515f-806r/uniform) trimmed to V4 region (515F–806R), available at https://www.arb-silva.de/

### R packages
| Package | Version used | Purpose |
|---|---|---|
| lme4 | 1.1-23 | Linear mixed models |
| lmerTest | 3.2-0 | p-values for LMMs (Satterthwaite's method) |
| vegan | 2.5-6 | PERMANOVA (adonis2), betadisper |
| ggplot2 | 3.3.2 | Figures |
| dplyr | 1.0.2 | Data manipulation |
| tidyr | 1.1.2 | Data reshaping |
| readr | 1.4.0 | File I/O |

---

## Usage

### QIIME2 pipeline

Run from the project root directory with QIIME2 activated:

```bash
conda activate qiime2-amplicon-2026.1
bash qiime2_pipeline.sh
```

**Note:** The pipeline pauses at Step 5 for decontamination in R. Run the decontam R script (`04_decontam/thr0p1/run_decontam_thr0.1.R`) before proceeding to Step 6. See the comments in `qiime2_pipeline.sh` for details.

### R scripts

All R scripts find the project root automatically and can be run from any subdirectory:

```bash
Rscript R/beta_permanova_full.R
Rscript R/taxonomy_multilevel.R
# etc.
```

Scripts expect the folder structure produced by `qiime2_pipeline.sh` to be in place (i.e., `10_exports/` populated with exported distance matrices and taxonomy tables).

---

## Citation

If you use these scripts, please cite:

Henry CJ, Stelling EG, Lattin CR (in review). Effects of captivity and host species on the cloacal microbiomes of a native and an invasive anole. *Animal Microbiome*.

Scripts archived at: https://doi.org/10.5281/zenodo.20818038

---

## Contact

Christine R. Lattin  
Department of Biological Sciences, Louisiana State University  
[christinelattin@lsu.edu](mailto:christinelattin@lsu.edu)
