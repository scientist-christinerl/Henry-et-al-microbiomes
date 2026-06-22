#!/usr/bin/env Rscript
# =============================================================================
# genus_within_target_families.R
# Lattin Lab | LSU
#
# Identifies which genus (or genera) make up Helicobacteraceae and
# Mycoplasmataceae specifically. This is NOT a full systematic genus table
# (a full genus table is deliberately avoided here: V4 amplicon data resolves
# genus less reliably than family, and most genera would just be noise).
# Instead this is a targeted pull restricted to the two families already
# shown to differ significantly by species.
#
# For each genus found within each target family, reports:
#   1. Relative abundance as % of the WHOLE gut microbiome -- directly
#      comparable to the family-level numbers in taxonomy_summary_family.csv
#   2. Relative abundance as % WITHIN the family -- i.e. is this family
#      essentially one genus, or split across several
#   3. The fraction of family-level reads that could NOT be resolved to a
#      named genus. V4 doesn't always resolve genus as cleanly as it
#      resolves family -- this is reported explicitly rather than silently
#      dropped, so it is clear whether the genus call is clean or partial
#      before citing it.
#
# Requires the genus-level export from qiime2_pipeline.sh
# (Level 6 block) to have been run first.
#
# Outputs:
#   07_figures/genus_helicobacteraceae.csv
#   07_figures/genus_mycoplasmataceae.csv
#   Console: summary tables and genus-resolution diagnostics
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
})

TARGET_FAMILIES <- c("Helicobacteraceae", "Mycoplasmataceae")

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

# ---- Load metadata ----
meta_path <- "02_metadata/sample-metadata-with-subject.tsv"
if (!file.exists(meta_path)) meta_path <- "02_metadata/sample-metadata.tsv"
meta <- read.delim(meta_path, sep = "\t", comment.char = "", check.names = FALSE)
meta$SampleID  <- trimws(meta$`#SampleID`)
meta$Group     <- trimws(meta$Group)
meta$Timepoint <- trimws(meta$Timepoint)

if ("SampleType" %in% colnames(meta)) {
  meta <- meta[meta$SampleType == "Biological", ]
}

meta$Species   <- ifelse(meta$Group == "G", "Green anole", "Brown anole")
meta$Species   <- factor(meta$Species,   levels = c("Green anole", "Brown anole"))
meta$Timepoint <- factor(meta$Timepoint, levels = c("PRE", "POST"))

cat("Metadata loaded:", nrow(meta), "biological samples\n")

# ---- Load the genus-level table ----
tsv_path <- "10_exports/exported_genus_level/feature-table.tsv"
if (!file.exists(tsv_path)) {
  stop("'", tsv_path, "' not found. Run the genus-level (Level 6) block in ",
       "qiime2_pipeline.sh first.")
}

raw_lines <- readLines(tsv_path)
header_idx <- which(grepl("^#\\s*(OTU ID|Feature ID|FeatureID)", raw_lines,
                          ignore.case = TRUE))
skip_n <- if (length(header_idx) > 0) header_idx[1] - 1 else 1
df <- read.delim(tsv_path, skip = skip_n,
                 header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
colnames(df)[1] <- "Taxon"

sample_cols <- intersect(setdiff(colnames(df), "Taxon"), meta$SampleID)
cat("Samples matched between feature table and metadata:", length(sample_cols), "\n")

if (length(sample_cols) == 0) {
  stop("No sample IDs in '", tsv_path, "' match meta$SampleID. Check for ",
       "naming mismatches.")
}

# ---- Parse family and genus from each taxonomy string ----
extract_level <- function(tax_string, prefix) {
  parts <- trimws(strsplit(trimws(tax_string), ";")[[1]])
  match_part <- parts[startsWith(parts, prefix)]
  if (length(match_part) == 0) return(NA_character_)   # level not present at all
  name <- trimws(sub(prefix, "", match_part[1]))
  if (name == "" || grepl("^uncultured", name, ignore.case = TRUE) ||
      grepl("^unidentified", name, ignore.case = TRUE)) return("Unclassified")
  name
}

df$Family <- sapply(df$Taxon, extract_level, prefix = "f__")
df$Genus  <- sapply(df$Taxon, extract_level, prefix = "g__")
df$Genus[is.na(df$Genus)] <- "Unclassified"   # no g__ field at all in the string

# ---- Total depth per sample, from the FULL genus-level table ----
count_mat_full <- as.matrix(df[, sample_cols])
total_depth <- colSums(count_mat_full)
depth_lookup <- data.frame(SampleID = names(total_depth), TotalDepth = total_depth)

# ---- Process each target family ----
all_results <- list()

for (fam in TARGET_FAMILIES) {

  cat("\n\n============================================================\n")
  cat(fam, "\n")
  cat("============================================================\n")

  fam_rows <- df[df$Family == fam, , drop = FALSE]

  if (nrow(fam_rows) == 0) {
    cat("No rows found for '", fam, "' in the genus-level table -- check ",
        "spelling, or that the genus-level collapse ran correctly.\n", sep = "")
    next
  }

  fam_long <- fam_rows %>%
    select(Genus, all_of(sample_cols)) %>%
    pivot_longer(-Genus, names_to = "SampleID", values_to = "Count") %>%
    group_by(Genus, SampleID) %>%
    summarise(Count = sum(Count, na.rm = TRUE), .groups = "drop") %>%
    filter(!is.na(Genus), Genus != "")   # drop any residual NA/empty rows

  fam_total_per_sample <- fam_long %>%
    group_by(SampleID) %>%
    summarise(FamilyTotal = sum(Count, na.rm = TRUE), .groups = "drop")

  # ---- Genus-resolution diagnostic, pooled across all samples ----
  pooled <- fam_long %>%
    group_by(Genus) %>%
    summarise(TotalReads = sum(Count, na.rm = TRUE), .groups = "drop") %>%
    mutate(PctOfFamily = round(100 * TotalReads / sum(TotalReads, na.rm = TRUE), 1)) %>%
    arrange(desc(TotalReads))

  unresolved_pct <- pooled$PctOfFamily[pooled$Genus == "Unclassified"]
  if (length(unresolved_pct) == 0) unresolved_pct <- 0
  cat("\nGenus resolution: ", round(100 - unresolved_pct, 1), "% of ", fam,
      " reads resolved to a named genus; ", unresolved_pct,
      "% unresolved (uncultured/unidentified/no genus-level call)\n", sep = "")
  cat("\nGenera found within ", fam, " (pooled across all samples):\n", sep = "")
  print(as.data.frame(pooled), row.names = FALSE)

  # ---- Per-sample relative abundance: both denominators ----
  fam_long <- fam_long %>%
    left_join(depth_lookup, by = "SampleID") %>%
    left_join(fam_total_per_sample, by = "SampleID") %>%
    mutate(PctOfMicrobiome = 100 * Count / TotalDepth,
           PctOfFamily     = 100 * Count / FamilyTotal) %>%
    left_join(meta[, c("SampleID", "Species", "Timepoint")], by = "SampleID")

  # ---- Summarize by group: % of whole microbiome ----
  group_summary_microbiome <- fam_long %>%
    group_by(Genus, Species, Timepoint) %>%
    summarise(MeanPct = round(mean(PctOfMicrobiome, na.rm = TRUE), 2), .groups = "drop") %>%
    mutate(GroupLabel = paste(Species, Timepoint)) %>%
    select(Genus, GroupLabel, MeanPct) %>%
    pivot_wider(names_from = GroupLabel, values_from = MeanPct, values_fill = 0)

  cat("\n--- Genus relative abundance, % of WHOLE gut microbiome ---\n")
  print(as.data.frame(group_summary_microbiome), row.names = FALSE)

  # ---- Summarize by group: % within the family ----
  group_summary_within <- fam_long %>%
    group_by(Genus, Species, Timepoint) %>%
    summarise(MeanPct = round(mean(PctOfFamily, na.rm = TRUE), 1), .groups = "drop") %>%
    mutate(GroupLabel = paste(Species, Timepoint)) %>%
    select(Genus, GroupLabel, MeanPct) %>%
    pivot_wider(names_from = GroupLabel, values_from = MeanPct, values_fill = 0)

  cat("\n--- Genus relative abundance, % WITHIN ", fam, " ---\n", sep = "")
  print(as.data.frame(group_summary_within), row.names = FALSE)

  # ---- Save combined CSV (both denominators side by side) ----
  out_path <- paste0("07_figures/genus_", tolower(fam), ".csv")
  combined <- group_summary_microbiome %>%
    rename_with(~ paste0(.x, "_pctMicrobiome"), -Genus) %>%
    left_join(
      group_summary_within %>% rename_with(~ paste0(.x, "_pctWithinFamily"), -Genus),
      by = "Genus"
    )
  write_csv(combined, out_path)
  cat("\nSaved:", out_path, "\n")

  all_results[[fam]] <- combined
}

cat("\n\n============================================================\n")
cat("Done. See the 'Genus resolution' line for each family above:\n")
cat("where a large fraction is unresolved, the named-genus result is\n")
cat("partial and should be reported as such.\n")
cat("============================================================\n")
