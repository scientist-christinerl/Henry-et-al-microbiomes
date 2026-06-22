#!/usr/bin/env Rscript
# =============================================================================
# check_sequencing_depth.R
# Lattin Lab | LSU
#
# Checks whether post-filtering sequencing depth (total reads per sample in
# the final cleaned/decontaminated feature table) varies systematically by
# Species, Timepoint, or their interaction. This matters because depth
# varies substantially in this dataset (min ~22,440, median ~505,180 reads,
# a ~22-fold spread) -- if that variation tracks the comparison of interest
# rather than being scattered randomly, depth is a real confound that needs
# to be modeled explicitly (e.g. as a covariate), not just noise that
# averages out.
#
# Reuses the family-level collapsed feature table already exported by
# qiime2_pipeline.sh. Taxa-collapsing doesn't lose or duplicate
# any reads, so summing all taxa within a sample gives the exact same total
# depth as the original (uncollapsed) cleaned/decontaminated table used
# throughout this analysis -- no separate export needed.
#
# Outputs:
#   07_figures/sequencing_depth_by_group.pdf/.png
#   Console: summary stats and statistical tests for depth ~ Species*Timepoint
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

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

# ---- Load metadata (prefer the version with Subject, for the paired model) ----
meta_path <- "02_metadata/sample-metadata-with-subject.tsv"
if (!file.exists(meta_path)) {
  warning("'", meta_path, "' not found -- falling back to sample-metadata.tsv. ",
          "The paired mixed model with (1|Subject) will not be available; ",
          "only a fixed-effects model will run.")
  meta_path <- "02_metadata/sample-metadata.tsv"
}
meta <- read.delim(meta_path, sep = "\t", comment.char = "", check.names = FALSE)
meta$SampleID  <- trimws(meta$`#SampleID`)
meta$Group     <- trimws(meta$Group)
meta$Timepoint <- trimws(meta$Timepoint)
has_subject <- "Subject" %in% colnames(meta)
if (has_subject) meta$Subject <- factor(trimws(meta$Subject))

if ("SampleType" %in% colnames(meta)) {
  meta <- meta[meta$SampleType == "Biological", ]
}

meta$Species   <- ifelse(meta$Group == "G", "Green anole", "Brown anole")
meta$Species   <- factor(meta$Species,   levels = c("Green anole", "Brown anole"))
meta$Timepoint <- factor(meta$Timepoint, levels = c("PRE", "POST"))

cat("Metadata loaded:", nrow(meta), "biological samples\n")
cat("Subject column available for paired model:", has_subject, "\n")

# ---- Load family-level table and sum to total depth per sample ----
# Reuses the same header-detection fix as taxonomy_multilevel.R: a standard
# biom-convert export has a header row that itself starts with "#", which a
# naive "skip every line starting with #" approach would mis-skip.
tsv_path <- "10_exports/exported_family_level/feature-table.tsv"
if (!file.exists(tsv_path)) {
  stop("'", tsv_path, "' not found. Run qiime2_pipeline.sh first ",
       "to generate the family-level export.")
}

raw_lines <- readLines(tsv_path)
header_idx <- which(grepl("^#\\s*(OTU ID|Feature ID|FeatureID)", raw_lines,
                          ignore.case = TRUE))
skip_n <- if (length(header_idx) > 0) header_idx[1] - 1 else 1
df <- read.delim(tsv_path, skip = skip_n,
                 header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
colnames(df)[1] <- "Taxon"

join_cols <- c("SampleID", "Species", "Timepoint")
if (has_subject) join_cols <- c(join_cols, "Subject")

depth <- df %>%
  pivot_longer(-Taxon, names_to = "SampleID", values_to = "Count") %>%
  group_by(SampleID) %>%
  summarise(Depth = sum(Count, na.rm = TRUE), .groups = "drop") %>%
  inner_join(meta[, join_cols], by = "SampleID")

cat("\nSamples with depth + metadata matched:", nrow(depth), "\n")

# Sanity check: expected to be near min ~22,440 / median ~505,180 reads if
# this is reading the same final filtered table used throughout the analysis.
cat("\nOverall depth summary:\n")
cat("  Min:    ", min(depth$Depth), "\n")
cat("  Median: ", median(depth$Depth), "\n")
cat("  Max:    ", max(depth$Depth), "\n")
cat("  Mean:   ", round(mean(depth$Depth), 0), "\n")

# ---- Summary by group ----
cat("\n--- Depth summary by Species x Timepoint ---\n")
group_summary <- depth %>%
  group_by(Species, Timepoint) %>%
  summarise(N = n(),
            Mean   = round(mean(Depth), 0),
            Median = round(median(Depth), 0),
            Min    = min(Depth),
            Max    = max(Depth),
            SD     = round(sd(Depth), 0),
            .groups = "drop")
print(as.data.frame(group_summary), row.names = FALSE)

# ---- Boxplot (log scale, given the skew) ----
p <- ggplot(depth, aes(x = Timepoint, y = Depth, fill = Species)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6, position = position_dodge(0.75)) +
  geom_jitter(aes(color = Species),
             position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75),
             size = 1.5, alpha = 0.7) +
  scale_y_log10(labels = scales::comma) +
  scale_fill_manual(values = c("Green anole" = "#2E8B57", "Brown anole" = "#8B4513")) +
  scale_color_manual(values = c("Green anole" = "#2E8B57", "Brown anole" = "#8B4513")) +
  theme_classic(base_size = 12) +
  labs(y = "Post-filtering sequencing depth (log scale)",
       x = "Timepoint",
       title = "Sequencing depth by species and timepoint")

ggsave("07_figures/sequencing_depth_by_group.pdf", p, width = 6, height = 5)
ggsave("07_figures/sequencing_depth_by_group.png", p, width = 6, height = 5, dpi = 300)
cat("\nSaved: 07_figures/sequencing_depth_by_group.pdf/.png\n")

# ---- Statistical test: does depth differ by Species / Timepoint / interaction? ----
depth$LogDepth <- log(depth$Depth)

cat("\n--- Testing whether log(depth) differs by Species x Timepoint ---\n")

if (has_subject) {
  suppressPackageStartupMessages(library(lmerTest))
  cat("\nLinear mixed model: log(Depth) ~ Species * Timepoint + (1|Subject)\n")
  m <- lmerTest::lmer(LogDepth ~ Species * Timepoint + (1 | Subject), data = depth)
  print(summary(m))
  cat("\n")
  print(anova(m))
} else {
  cat("\nNo Subject column available -- running a fixed-effects model only.\n")
  m <- lm(LogDepth ~ Species * Timepoint, data = depth)
  print(summary(m))
  print(anova(m))
}

# Non-parametric cross-check, since depth distributions are often heavily
# right-skewed even after log transformation
cat("\n--- Kruskal-Wallis (non-parametric cross-check) ---\n")
cat("\nBy Species:\n")
print(kruskal.test(Depth ~ Species, data = depth))
cat("\nBy Timepoint:\n")
print(kruskal.test(Depth ~ Timepoint, data = depth))

# Interpretation: a significant Species and/or Timepoint effect on depth
# indicates depth covaries with the comparison of interest and should be
# modeled as a covariate rather than assumed to average out. A significant
# Timepoint effect specifically warrants checking whether PRE and POST were
# run on different sequencing batches (a technical rather than biological
# source of the pattern).
cat("\nDone.\n")
