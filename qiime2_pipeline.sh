#!/bin/bash
# =============================================================================
# qiime2_pipeline.sh
#
# Title:       Effects of captivity and species on the cloacal microbiomes
#              of an invasive and native anole
# Description: 16S rRNA amplicon sequencing analysis of green vs brown anole
#              gut microbiomes from cloacal lavage water
# Authors:     Caroline J. Henry, Emily G. Stelling, Christine R. Lattin
# Analysis:    2026
# Pipeline:    QIIME2 amplicon 2026.1
# Target:      16S rRNA V4 region (515F-806R)
#
# Run from the project root directory with QIIME2 activated:
#   conda activate qiime2-amplicon-2026.1
#   bash 08_scripts/qiime2_pipeline.sh
#
# Note: decontam (step 5) is run in R between the QIIME2 steps.
# The R script is 04_decontam/thr0p1/run_decontam_thr0.1.R
# Run that script after step 5 and before step 6.
# =============================================================================

# =============================================================================
# STEP 1: Import paired-end reads using manifest
# =============================================================================
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path 02_metadata/manifest.tsv \
  --output-path 03_qiime2_artifacts/paired-end-demux.qza \
  --input-format PairedEndFastqManifestPhred33

# =============================================================================
# STEP 2: DADA2 denoising
# Truncation lengths: forward 240 bp, reverse 220 bp
# =============================================================================
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs 03_qiime2_artifacts/paired-end-demux.qza \
  --p-trunc-len-f 240 \
  --p-trunc-len-r 220 \
  --o-table 03_qiime2_artifacts/table-240-220.qza \
  --o-representative-sequences 03_qiime2_artifacts/rep-seqs-240-220.qza \
  --o-denoising-stats 03_qiime2_artifacts/denoising-stats-240-220.qza

# =============================================================================
# STEP 3: Taxonomic assignment
# Classifier: SILVA 138.2 SSURef NR99, trimmed to V4 region (515F-806R)
# =============================================================================
qiime feature-classifier classify-sklearn \
  --i-classifier 09_reference/SILVA/SILVA138.2_SSURef_NR99_uniform_classifier_V4-515f-806r.qza \
  --i-reads 03_qiime2_artifacts/rep-seqs-240-220.qza \
  --o-classification 03_qiime2_artifacts/taxonomy-240-220.qza \
  --p-n-jobs 0

qiime metadata tabulate \
  --m-input-file 03_qiime2_artifacts/taxonomy-240-220.qza \
  --o-visualization 03_qiime2_artifacts/taxonomy-240-220.qzv

# =============================================================================
# STEP 4: Remove negative and positive controls; summarize to assess depth
# =============================================================================

# Filter to biological samples only
qiime feature-table filter-samples \
  --i-table 03_qiime2_artifacts/table-240-220.qza \
  --m-metadata-file 02_metadata/sample-metadata.tsv \
  --p-where "[SampleType]='Biological'" \
  --o-filtered-table 03_qiime2_artifacts/table-no-controls.qza

# Summarize to inspect read depth distribution and determine rarefaction cutoff
qiime feature-table summarize \
  --i-table 03_qiime2_artifacts/table-no-controls.qza \
  --m-metadata-file 02_metadata/sample-metadata.tsv \
  --o-feature-frequencies 03_qiime2_artifacts/feature-frequencies-no-controls.qza \
  --o-sample-frequencies 03_qiime2_artifacts/sample-frequencies-no-controls.qza \
  --o-summary 03_qiime2_artifacts/table-no-controls.qzv

# =============================================================================
# STEP 5: Prepare table for decontam in R
#
# decontam requires controls to be present (to identify contaminants by
# prevalence), so the export uses the original table with controls retained,
# but with the 6 low-depth paired samples removed first.
# paired-failures.tsv lists the 3 low-depth samples and their 3 paired samples.
#
# After this step, run the R decontam script:
#   04_decontam/thr0p1/run_decontam_thr0.1.R
# That script identifies contaminant ASVs and writes
#   04_decontam/thr0p1/contaminant-asvs_prevalence_thr-0.10.tsv
# which is used in Step 6 below.
# =============================================================================

# Remove the 6 low-depth paired samples (controls still included)
qiime feature-table filter-samples \
  --i-table 03_qiime2_artifacts/table-240-220.qza \
  --m-metadata-file 02_metadata/samples-to-remove.tsv \
  --p-exclude-ids \
  --o-filtered-table 03_qiime2_artifacts/table-decontam-input.qza

# Export feature table for decontam
mkdir -p 04_decontam/decontam_export/table
mkdir -p 04_decontam/decontam_export/tax

qiime tools export \
  --input-path 03_qiime2_artifacts/table-decontam-input.qza \
  --output-path 04_decontam/decontam_export/table

biom convert \
  -i 04_decontam/decontam_export/table/feature-table.biom \
  -o 04_decontam/decontam_export/table/feature-table.tsv \
  --to-tsv

# Export taxonomy for decontam
qiime tools export \
  --input-path 03_qiime2_artifacts/taxonomy-240-220.qza \
  --output-path 04_decontam/decontam_export/tax

# --- RUN R DECONTAM SCRIPT HERE before proceeding to Step 6 ---

# =============================================================================
# STEP 6: Remove contaminant ASVs identified by decontam (threshold 0.1)
# =============================================================================

# The table-paired-clean.qza is the decontam-input table filtered to
# biological samples only (i.e. controls already removed)
qiime feature-table filter-features \
  --i-table 04_decontam/thr0p1/table-paired-clean.qza \
  --m-metadata-file 04_decontam/thr0p1/contaminant-asvs_prevalence_thr-0.10.tsv \
  --p-exclude-ids \
  --o-filtered-table 04_decontam/thr0p1/table-paired-clean-decontam-0p1.qza

# Filter representative sequences to match the cleaned feature table
qiime feature-table filter-seqs \
  --i-data 03_qiime2_artifacts/rep-seqs-240-220.qza \
  --i-table 04_decontam/thr0p1/table-paired-clean-decontam-0p1.qza \
  --o-filtered-data 03_qiime2_artifacts/rep-seqs-paired-clean-decontam-0p1.qza

# =============================================================================
# STEP 7: Build phylogenetic tree (MAFFT + FastTree)
# Used for Faith's PD and UniFrac metrics
# =============================================================================
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences 03_qiime2_artifacts/rep-seqs-paired-clean-decontam-0p1.qza \
  --o-alignment 05_diversity/thr0p1/aligned-rep-seqs-clean-0p1.qza \
  --o-masked-alignment 05_diversity/thr0p1/masked-aligned-rep-seqs-clean-0p1.qza \
  --o-tree 05_diversity/thr0p1/unrooted-tree-clean-0p1.qza \
  --o-rooted-tree 05_diversity/thr0p1/rooted-tree-clean-0p1.qza

# =============================================================================
# STEP 8: Core diversity metrics (rarefaction depth = 20,000 reads)
# Computes Shannon, Faith's PD, Observed ASVs, Bray-Curtis, wUF, uUF, Jaccard
# =============================================================================
qiime diversity core-metrics-phylogenetic \
  --i-table 04_decontam/thr0p1/table-paired-clean-decontam-0p1.qza \
  --i-phylogeny 05_diversity/thr0p1/rooted-tree-clean-0p1.qza \
  --p-sampling-depth 20000 \
  --m-metadata-file 02_metadata/sample-metadata.tsv \
  --output-dir 05_diversity/thr0p1/core-metrics-20k-0p1

# =============================================================================
# STEP 9: Export alpha and beta diversity results for R
# =============================================================================

# Alpha diversity
qiime tools export \
  --input-path 05_diversity/thr0p1/core-metrics-20k-0p1/shannon_vector.qza \
  --output-path 10_exports/exported_shannon_0p1

qiime tools export \
  --input-path 05_diversity/thr0p1/core-metrics-20k-0p1/faith_pd_vector.qza \
  --output-path 10_exports/exported_faith_pd_0p1

qiime tools export \
  --input-path 05_diversity/thr0p1/core-metrics-20k-0p1/observed_features_vector.qza \
  --output-path 10_exports/exported_observed_0p1

# Beta diversity distance matrices
qiime tools export \
  --input-path 05_diversity/thr0p1/core-metrics-20k-0p1/bray_curtis_distance_matrix.qza \
  --output-path 10_exports/exported_bray_full_0p1

qiime tools export \
  --input-path 05_diversity/thr0p1/core-metrics-20k-0p1/weighted_unifrac_distance_matrix.qza \
  --output-path 10_exports/exported_wuf_full_0p1

qiime tools export \
  --input-path 05_diversity/thr0p1/core-metrics-20k-0p1/unweighted_unifrac_distance_matrix.qza \
  --output-path 10_exports/exported_uuf_full_0p1

qiime tools export \
  --input-path 05_diversity/thr0p1/core-metrics-20k-0p1/jaccard_distance_matrix.qza \
  --output-path 10_exports/exported_jaccard_full_0p1

# =============================================================================
# STEP 10: Collapse feature table to taxonomic levels for taxonomy figures
# SILVA 138.2 levels: 2=phylum, 4=order, 5=family, 6=genus
# =============================================================================

# Phylum (level 2)
qiime taxa collapse \
  --i-table 04_decontam/thr0p1/table-paired-clean-decontam-0p1.qza \
  --i-taxonomy 03_qiime2_artifacts/taxonomy-240-220.qza \
  --p-level 2 \
  --o-collapsed-table 03_qiime2_artifacts/table-phylum-level.qza

qiime tools export \
  --input-path 03_qiime2_artifacts/table-phylum-level.qza \
  --output-path 10_exports/exported_phylum_level

biom convert \
  -i 10_exports/exported_phylum_level/feature-table.biom \
  -o 10_exports/exported_phylum_level/feature-table.tsv \
  --to-tsv

# Order (level 4)
qiime taxa collapse \
  --i-table 04_decontam/thr0p1/table-paired-clean-decontam-0p1.qza \
  --i-taxonomy 03_qiime2_artifacts/taxonomy-240-220.qza \
  --p-level 4 \
  --o-collapsed-table 03_qiime2_artifacts/table-order-level.qza

qiime tools export \
  --input-path 03_qiime2_artifacts/table-order-level.qza \
  --output-path 10_exports/exported_order_level

biom convert \
  -i 10_exports/exported_order_level/feature-table.biom \
  -o 10_exports/exported_order_level/feature-table.tsv \
  --to-tsv

# Family (level 5)
qiime taxa collapse \
  --i-table 04_decontam/thr0p1/table-paired-clean-decontam-0p1.qza \
  --i-taxonomy 03_qiime2_artifacts/taxonomy-240-220.qza \
  --p-level 5 \
  --o-collapsed-table 03_qiime2_artifacts/table-family-level.qza

qiime tools export \
  --input-path 03_qiime2_artifacts/table-family-level.qza \
  --output-path 10_exports/exported_family_level

biom convert \
  -i 10_exports/exported_family_level/feature-table.biom \
  -o 10_exports/exported_family_level/feature-table.tsv \
  --to-tsv

# Genus (level 6)
qiime taxa collapse \
  --i-table 04_decontam/thr0p1/table-paired-clean-decontam-0p1.qza \
  --i-taxonomy 03_qiime2_artifacts/taxonomy-240-220.qza \
  --p-level 6 \
  --o-collapsed-table 03_qiime2_artifacts/table-genus-level.qza

qiime tools export \
  --input-path 03_qiime2_artifacts/table-genus-level.qza \
  --output-path 10_exports/exported_genus_level

biom convert \
  -i 10_exports/exported_genus_level/feature-table.biom \
  -o 10_exports/exported_genus_level/feature-table.tsv \
  --to-tsv

echo "QIIME2 pipeline complete."
echo "Next steps: run R scripts in 08_scripts/ for statistical analysis and figures."
