#!/usr/bin/env Rscript
# =============================================================================
# clr_mixed_model_taxa.R
# Lattin Lab | LSU
#
# CLR-transformed linear mixed model testing whether two pre-specified
# families -- Helicobacteraceae and Mycoplasmataceae -- differ by host
# species, timepoint, or their interaction.
#
# WHY THIS APPROACH (vs. a full DA package like ANCOM-BC/ALDEx2/LinDA):
# these two taxa were flagged by inspecting the descriptive taxonomy
# tables, not from a screen across hundreds of features needing heavy
# multiple-testing correction. A targeted CLR mixed model mirrors the same
# fixed/random-effects structure already used for alpha diversity
# (Group * Timepoint + (1|Subject)) and is easier to justify in Methods
# than introducing a new package for two taxa.
#
# IMPORTANT METHODS NOTE: CLR's geometric-mean reference frame must be
# computed across the FULL composition in each sample, not just the taxa
# being tested -- filtering down to only the 2 taxa of interest first would
# destroy the compositional structure CLR is meant to preserve. This script:
#   1. Loads the FULL family-level table (every family, not the >=1%
#      filtered table used for the manuscript's descriptive tables)
#   2. Applies a modest prevalence filter (default: present in >=10% of
#      samples) to drop very sparse/noisy features before computing CLR
#   3. Computes CLR per sample across that full filtered composition
#   4. THEN extracts just the two columns of interest for modeling
#
# DEPTH COVARIATE: sequencing depth was shown not to differ by Species or
# Timepoint (check_sequencing_depth.R), so it is not a confound here. It is
# still tested as an optional covariate because these two taxa are
# low-abundance/near-zero in much of the dataset, and detection of a rare
# taxon is sensitive to depth independent of any confounding with the
# predictors of interest. Both models are reported. The without-depth model
# parallels the alpha-diversity model and is treated as primary; the
# depth-adjusted model is a sensitivity check.
#
# A NOTE ON ZEROS: Helicobacteraceae is undetected (0 reads) in every green
# anole sample in this dataset. With one whole group at zero, this is closer
# to a presence/absence comparison than a continuous differential-abundance
# test, and the CLR estimate for that taxon is sensitive to the PSEUDOCOUNT
# choice below; this should be noted when interpreting that result.
# Mycoplasmataceae is present (at low level) in both species and is a more
# conventional differential-abundance comparison.
#
# Outputs:
#   07_figures/clr_helicobacteraceae_mycoplasmataceae.pdf/.png
#   Console: model summaries for both taxa, with and without depth covariate
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(lmerTest)
})

# ---- Settings (adjust here if needed) ----
PREVALENCE_MIN <- 0.10   # keep families present (count > 0) in >=10% of samples
PSEUDOCOUNT    <- 1      # added to all counts before log transform
TARGET_TAXA    <- c("Helicobacteraceae", "Mycoplasmataceae")

# ---- Locate project root regardless of where this script is launched from ----
find_project_root <- function() {
  for (up in 0:4) {
    candidate <- if (up == 0) "." else paste(rep("..", up), collapse = "/")
    if (dir.exists(file.path(candidate, "02_metadata")) &&
        dir.exists(file.path(candidate, "10_exports"))) {
      return(normalizePath(candidate))
    }
  }
  stop("Could not locate project root (a folder containing both 02_metadata/ ",
       "and 10_exports/) within 4 levels of the current directory. Either cd ",
       "to the project root before running this script, or check that those ",
       "folders exist.")
}

setwd(find_project_root())
cat("Working directory set to project root:", getwd(), "\n")
dir.create("07_figures", showWarnings = FALSE, recursive = TRUE)

# ---- Load metadata (requires Subject, for the paired model) ----
meta_path <- "02_metadata/sample-metadata-with-subject.tsv"
if (!file.exists(meta_path)) {
  stop("'", meta_path, "' not found -- this script requires the Subject ",
       "column for the paired mixed model.")
}
meta <- read.delim(meta_path, sep = "\t", comment.char = "", check.names = FALSE)
meta$SampleID  <- trimws(meta$`#SampleID`)
meta$Group     <- trimws(meta$Group)
meta$Timepoint <- trimws(meta$Timepoint)
meta$Subject   <- factor(trimws(meta$Subject))

if ("SampleType" %in% colnames(meta)) {
  meta <- meta[meta$SampleType == "Biological", ]
}

meta$Species   <- ifelse(meta$Group == "G", "Green anole", "Brown anole")
meta$Species   <- factor(meta$Species,   levels = c("Green anole", "Brown anole"))
meta$Timepoint <- factor(meta$Timepoint, levels = c("PRE", "POST"))

cat("Metadata loaded:", nrow(meta), "biological samples\n")

# ---- Load the FULL family-level table (not the >=1%-filtered summary table) ----
tsv_path <- "10_exports/exported_family_level/feature-table.tsv"
if (!file.exists(tsv_path)) {
  stop("'", tsv_path, "' not found. Run qiime2_pipeline.sh first.")
}

raw_lines <- readLines(tsv_path)
header_idx <- which(grepl("^#\\s*(OTU ID|Feature ID|FeatureID)", raw_lines,
                          ignore.case = TRUE))
skip_n <- if (length(header_idx) > 0) header_idx[1] - 1 else 1
df <- read.delim(tsv_path, skip = skip_n,
                 header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
colnames(df)[1] <- "Taxon"

# ---- Parse SILVA taxonomy string down to family name ----
extract_family <- function(tax_string) {
  parts <- trimws(strsplit(trimws(tax_string), ";")[[1]])
  match_part <- parts[startsWith(parts, "f__")]
  if (length(match_part) == 0) return("Unclassified")
  name <- sub("f__", "", match_part[1])
  if (name == "" || name == "uncultured" || name == "unidentified") return("Unclassified")
  name
}
df$Family <- sapply(df$Taxon, extract_family)

# Collapse any rows sharing the same family name (summing), in case the
# original table wasn't already fully collapsed to this level
count_cols <- setdiff(colnames(df), c("Taxon", "Family"))
df_fam <- df %>%
  group_by(Family) %>%
  summarise(across(all_of(count_cols), sum), .groups = "drop")

sample_cols <- intersect(count_cols, meta$SampleID)
cat("Samples matched between feature table and metadata:", length(sample_cols),
    "of", length(count_cols), "\n")

count_mat <- as.matrix(df_fam[, sample_cols])
rownames(count_mat) <- df_fam$Family

# Total depth per sample, from the FULL (unfiltered) table -- used as the
# optional covariate below
total_depth <- colSums(count_mat)

# ---- Prevalence filter (applied to the full composition, before CLR) ----
prevalence <- rowMeans(count_mat > 0)
keep <- prevalence >= PREVALENCE_MIN
cat("\nFamilies retained at >=", PREVALENCE_MIN * 100, "% prevalence: ",
    sum(keep), " of ", nrow(count_mat), " total\n", sep = "")

missing_targets <- setdiff(TARGET_TAXA, rownames(count_mat)[keep])
if (length(missing_targets) > 0) {
  warning("Target taxa not retained after prevalence filtering: ",
          paste(missing_targets, collapse = ", "),
          "; lower PREVALENCE_MIN to retain them.")
}

count_mat <- count_mat[keep, , drop = FALSE]

# ---- CLR transform (per sample, across the full retained composition) ----
count_mat_pc <- count_mat + PSEUDOCOUNT
log_mat <- log(count_mat_pc)
sample_log_mean <- colMeans(log_mat)                 # mean of logs per sample
clr_mat <- sweep(log_mat, 2, sample_log_mean, "-")    # subtract per-sample mean

cat("CLR transform complete:", nrow(clr_mat), "families x", ncol(clr_mat), "samples\n")

# ---- Build modeling data frame for the two target taxa ----
present_targets <- TARGET_TAXA[TARGET_TAXA %in% rownames(clr_mat)]
clr_long <- as.data.frame(t(clr_mat[present_targets, , drop = FALSE]))
clr_long$SampleID  <- rownames(clr_long)
clr_long$Depth     <- total_depth[clr_long$SampleID]
clr_long$LogDepth  <- log(clr_long$Depth)

model_df <- clr_long %>%
  left_join(meta[, c("SampleID", "Species", "Timepoint", "Subject")],
            by = "SampleID")

cat("\nFinal modeling dataset:", nrow(model_df), "samples\n")

# ---- Run models for each target taxon ----
for (taxon in present_targets) {

  cat("\n\n============================================================\n")
  cat(taxon, "\n")
  cat("============================================================\n")

  model_df$CLR_y <- model_df[[taxon]]

  cat("\n--- Model WITHOUT depth covariate (primary) ---\n")
  cat("CLR_", taxon, " ~ Species * Timepoint + (1|Subject)\n", sep = "")
  m1 <- lmerTest::lmer(CLR_y ~ Species * Timepoint + (1 | Subject), data = model_df)
  print(anova(m1))

  cat("\n--- Model WITH log(depth) covariate (sensitivity check) ---\n")
  cat("CLR_", taxon, " ~ Species * Timepoint + LogDepth + (1|Subject)\n", sep = "")
  m2 <- lmerTest::lmer(CLR_y ~ Species * Timepoint + LogDepth + (1 | Subject), data = model_df)
  print(anova(m2))

  cat("\nFixed effect estimates (without depth covariate):\n")
  print(summary(m1)$coefficients)
}

# ---- Multiple comparisons note ----
cat("\n\n============================================================\n")
cat("MULTIPLE COMPARISONS\n")
cat("============================================================\n")
cat("Two pre-specified taxa were tested, not a screen across all families,\n")
cat("so heavy correction is not strictly required. For BH-adjusted Species-\n")
cat("effect p-values across both taxa, take the relevant p-values from the\n")
cat("anova() tables above and run p.adjust(c(p1, p2), method = 'BH').\n")

# ---- Visualization ----
plot_df <- model_df %>%
  select(SampleID, Species, Timepoint, all_of(present_targets)) %>%
  pivot_longer(cols = all_of(present_targets), names_to = "Taxon", values_to = "CLR")

p <- ggplot(plot_df, aes(x = Timepoint, y = CLR, fill = Species)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6, position = position_dodge(0.75)) +
  geom_jitter(aes(color = Species),
             position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.75),
             size = 1.5, alpha = 0.7) +
  facet_wrap(~ Taxon, scales = "free_y") +
  scale_fill_manual(values = c("Green anole" = "#2E8B57", "Brown anole" = "#8B4513")) +
  scale_color_manual(values = c("Green anole" = "#2E8B57", "Brown anole" = "#8B4513")) +
  theme_classic(base_size = 12) +
  theme(strip.background = element_blank(),
       strip.text = element_text(face = "bold.italic")) +
  labs(y = "CLR-transformed abundance", x = "Timepoint",
       title = "Candidate taxa: CLR abundance by species and timepoint")

ggsave("07_figures/clr_helicobacteraceae_mycoplasmataceae.pdf", p, width = 8, height = 5)
ggsave("07_figures/clr_helicobacteraceae_mycoplasmataceae.png", p, width = 8, height = 5, dpi = 300)
cat("\nSaved: 07_figures/clr_helicobacteraceae_mycoplasmataceae.pdf/.png\n")

cat("\nDone.\n")
