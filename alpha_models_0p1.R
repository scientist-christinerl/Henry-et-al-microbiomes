#!/usr/bin/env Rscript
# =============================================================================
# alpha_models_0p1.R
# Lattin Lab | LSU
#
# Linear mixed models for the three alpha diversity metrics (Shannon,
# Faith's PD, observed ASVs) at the 0.1 decontam threshold. Each metric is
# modeled as:
#   metric ~ Group * Timepoint + (1 | Subject)
# i.e. the paired PRE/POST design is handled with a per-individual random
# intercept. Type III ANOVA p-values are written to a summary CSV.
#
# Outputs:
#   06_stats/thr0p1/alpha_mixed_models_0p1_results.csv
#   Console: model summary and Type III ANOVA for each metric
# =============================================================================

suppressPackageStartupMessages({
  library(lme4)
  library(lmerTest)
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
dir.create("06_stats/thr0p1", showWarnings = FALSE, recursive = TRUE)

# ---- Inputs ----
meta_file <- "02_metadata/sample-metadata-with-subject.tsv"

alpha_files <- list(
  Shannon  = "10_exports/exported_shannon_0p1/alpha-diversity.tsv",
  FaithPD  = "10_exports/exported_faith_pd_0p1/alpha-diversity.tsv",
  Observed = "10_exports/exported_observed_0p1/alpha-diversity.tsv"
)

out_csv <- "06_stats/thr0p1/alpha_mixed_models_0p1_results.csv"

# ---- Helper: run model for one metric ----
run_alpha_model <- function(metric_name, alpha_path, meta_path) {
  if (!file.exists(alpha_path)) {
    stop("Alpha diversity file not found: ", alpha_path)
  }
  if (!file.exists(meta_path)) {
    stop("Metadata file not found: ", meta_path)
  }

  # Read alpha vector
  alpha <- read.delim(alpha_path,
                      header = TRUE, sep = "\t",
                      comment.char = "", check.names = FALSE)

  # Standardize column names:
  # exported alpha-diversity.tsv is: sample-id <tab> alpha-diversity
  colnames(alpha)[1] <- "SampleID"
  colnames(alpha)[2] <- metric_name

  # Read metadata
  meta <- read.delim(meta_path,
                     header = TRUE, sep = "\t",
                     comment.char = "", check.names = FALSE)

  # Merge
  dat <- merge(alpha, meta, by.x = "SampleID", by.y = "#SampleID", all = FALSE)

  # Trim and set factors (same coding used in the PERMANOVA script)
  dat$Group <- factor(trimws(dat$Group))
  dat$Timepoint <- factor(trimws(dat$Timepoint))
  dat$Subject <- factor(trimws(dat$Subject))

  # Drop rows with missing values
  dat <- dat[complete.cases(dat[, c(metric_name, "Group", "Timepoint", "Subject")]), ]

  if (nrow(dat) < 4) {
    stop("Too few samples after merging for metric ", metric_name,
         " (n=", nrow(dat), "). Check SampleIDs match.")
  }

  cat("\n==============================\n")
  cat("Metric:", metric_name, "\n")
  cat("Alpha file:", alpha_path, "\n")
  cat("Samples used:", nrow(dat), "\n")
  cat("==============================\n")

  # Mixed model (paired via random intercept)
  form <- as.formula(paste0(metric_name, " ~ Group * Timepoint + (1 | Subject)"))
  model <- lmer(form, data = dat, REML = TRUE)

  cat("\nModel summary:\n")
  print(summary(model))

  cat("\nANOVA (Type III fixed effects):\n")
  a <- anova(model, type = 3)
  print(a)

  # Pull p-values (Type III)
  # rownames include: (Intercept), Group, Timepoint, Group:Timepoint
  out <- data.frame(
    Metric = metric_name,
    Term = rownames(a),
    NumDF = a$NumDF,
    DenDF = a$DenDF,
    F_value = a$`F value`,
    P_value = a$`Pr(>F)`,
    stringsAsFactors = FALSE
  )

  return(out)
}

# ---- Run all metrics ----
all_results <- do.call(rbind, lapply(names(alpha_files), function(m) {
  run_alpha_model(m, alpha_files[[m]], meta_file)
}))

# Write summary CSV
write.csv(all_results, out_csv, row.names = FALSE)
cat("\nWrote summary table:", out_csv, "\n")
