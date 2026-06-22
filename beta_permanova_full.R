#!/usr/bin/env Rscript
# =============================================================================
# beta_permanova_full.R
# Lattin Lab | LSU
#
# Beta diversity PERMANOVA and dispersion tests for all four distance metrics:
#   Bray-Curtis, Weighted UniFrac, Unweighted UniFrac, Jaccard
#
# Note on paired PERMANOVA: adonis2() does not accept lme4-style random
# effects (e.g. (1|Subject)) in its formula. For the paired/repeated-measures
# design, permutations are constrained within individuals using strata =.
#
# Three tests are run for each metric:
#   1. Species effect:        PRE-only samples, no strata (between-subjects)
#   2. Timepoint effect:      all samples, strata = Subject (within-individual)
#   3. Species x Timepoint:   marginal term, by = "margin", strata = Subject
#
# Dispersion tests (betadisper) are run for Group and Timepoint.
# =============================================================================

suppressPackageStartupMessages(library(vegan))

# ---- Locate project root regardless of where this script is launched from ----
# Works when run from the project root, from 06_stats/thr0p1,
# or anywhere else inside the project tree -- searches up to 4 levels up
# for the folder that contains both 02_metadata/ and 10_exports/.
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

project_root <- find_project_root()
setwd(project_root)
cat("Working directory set to project root:", getwd(), "\n")

# ---- Set seed for reproducibility ----
set.seed(42)

# ---- Load metadata ----
# Uses sample-metadata-with-subject.tsv, which includes the Subject column
# needed for the strata argument in paired PERMANOVA tests.
meta_path <- "02_metadata/sample-metadata-with-subject.tsv"
meta <- read.delim(meta_path, sep = "\t", comment.char = "", check.names = FALSE)

required_cols <- c("#SampleID", "Group", "Timepoint", "Subject")
missing_cols <- setdiff(required_cols, colnames(meta))
if (length(missing_cols) > 0) {
  stop("Metadata file '", meta_path, "' is missing required column(s): ",
       paste(missing_cols, collapse = ", "), ".\n",
       "Columns found: ", paste(colnames(meta), collapse = ", "))
}

meta$SampleID  <- trimws(meta$`#SampleID`)
meta$Group     <- factor(trimws(meta$Group))
meta$Timepoint <- factor(trimws(meta$Timepoint))
meta$Subject   <- factor(trimws(meta$Subject))

cat("Metadata loaded:", nrow(meta), "rows\n")
cat("Groups:", levels(meta$Group), "\n")
cat("Timepoints:", levels(meta$Timepoint), "\n")
cat("Unique subjects:", nlevels(meta$Subject), "\n")

# ---- Helper: subset a square distance matrix to a set of IDs ----
subset_dist <- function(dm_matrix, ids) {
  ids_present <- ids[ids %in% rownames(dm_matrix)]
  as.dist(dm_matrix[ids_present, ids_present])
}

# ---- Main PERMANOVA function ----
run_permanova <- function(dist_path, label, out_dir = "06_stats/thr0p1") {

  cat("\n\n============================================================\n")
  cat(label, "\n")
  cat("============================================================\n")

  # Load and validate distance matrix
  if (!file.exists(dist_path)) {
    cat("  !! File not found:", dist_path, "-- skipping\n")
    return(invisible(NULL))
  }

  dm <- read.delim(dist_path, header = TRUE, row.names = 1, check.names = FALSE)
  dm <- as.matrix(dm)
  mode(dm) <- "numeric"

  # Align to metadata (keep only samples present in both)
  ids    <- intersect(rownames(dm), meta$SampleID)
  dm     <- dm[ids, ids]
  meta_m <- meta[match(ids, meta$SampleID), ]
  dist_full <- as.dist(dm)

  cat("\nSamples in analysis:", length(ids), "\n")
  cat("PRE samples:", sum(meta_m$Timepoint == "PRE"),
      "  POST samples:", sum(meta_m$Timepoint == "POST"), "\n")

  # ------------------------------------------------------------------
  # 1. Species effect: PRE-only samples, no permutation strata
  # ------------------------------------------------------------------
  pre_ids  <- meta_m$SampleID[meta_m$Timepoint == "PRE"]
  meta_pre <- meta_m[match(pre_ids, meta_m$SampleID), ]
  dist_pre <- subset_dist(dm, pre_ids)

  cat("\n--- 1. Species effect (PRE only, n =", length(pre_ids), ") ---\n")
  res_sp <- adonis2(dist_pre ~ Group, data = meta_pre, permutations = 999)
  print(res_sp)

  # ------------------------------------------------------------------
  # 2. Timepoint effect: all samples, strata = Subject
  # ------------------------------------------------------------------
  cat("\n--- 2. Timepoint effect (all samples, strata = Subject) ---\n")
  res_tp <- adonis2(dist_full ~ Timepoint, data = meta_m,
                    permutations = 999, strata = meta_m$Subject)
  print(res_tp)

  # ------------------------------------------------------------------
  # 3. Species x Timepoint interaction: marginal, strata = Subject
  # ------------------------------------------------------------------
  cat("\n--- 3. Species x Timepoint interaction (marginal, strata = Subject) ---\n")
  res_int <- adonis2(dist_full ~ Group * Timepoint, data = meta_m,
                     by = "margin", permutations = 999,
                     strata = meta_m$Subject)
  print(res_int)

  # ------------------------------------------------------------------
  # 4. Dispersion tests
  # ------------------------------------------------------------------
  cat("\n--- 4a. Dispersion: Group (PRE only) ---\n")
  bd_sp <- betadisper(dist_pre, meta_pre$Group)
  perm_sp <- permutest(bd_sp)
  print(perm_sp)

  cat("\n--- 4b. Dispersion: Timepoint (all samples) ---\n")
  bd_tp <- betadisper(dist_full, meta_m$Timepoint)
  perm_tp <- permutest(bd_tp)
  print(perm_tp)

  # ------------------------------------------------------------------
  # Save results to file
  # ------------------------------------------------------------------
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  safe_label <- gsub("[^A-Za-z0-9]", "_", tolower(label))
  out_file <- file.path(out_dir, paste0("permanova_", safe_label, ".txt"))

  sink(out_file)
  cat("PERMANOVA results:", label, "\n")
  cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("Distance matrix:", dist_path, "\n")
  cat("Samples used:", length(ids), "\n\n")
  cat("--- 1. Species effect (PRE only) ---\n"); print(res_sp)
  cat("\n--- 2. Timepoint effect (strata = Subject) ---\n"); print(res_tp)
  cat("\n--- 3. Interaction (marginal, strata = Subject) ---\n"); print(res_int)
  cat("\n--- 4a. Dispersion: Group (PRE only) ---\n"); print(perm_sp)
  cat("\n--- 4b. Dispersion: Timepoint ---\n"); print(perm_tp)
  sink()

  cat("\nResults saved to:", out_file, "\n")
  invisible(list(species = res_sp, timepoint = res_tp, interaction = res_int))
}

# ---- Run for all four metrics ----
# Bray-Curtis (existing export)
run_permanova(
  "10_exports/exported_bray_full_0p1/distance-matrix.tsv",
  "Bray-Curtis"
)

# Weighted UniFrac (existing export)
run_permanova(
  "10_exports/exported_wuf_full_0p1/distance-matrix.tsv",
  "Weighted UniFrac"
)

# Unweighted UniFrac (existing export)
run_permanova(
  "10_exports/exported_uuf_full_0p1/distance-matrix.tsv",
  "Unweighted UniFrac"
)

# Jaccard
run_permanova(
  "10_exports/exported_jaccard_full_0p1/distance-matrix.tsv",
  "Jaccard"
)

cat("\n\nAll PERMANOVA analyses complete.\n")
