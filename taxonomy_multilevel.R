#!/usr/bin/env Rscript
# =============================================================================
# taxonomy_multilevel.R
# Lattin Lab | LSU
#
# Produces stacked bar charts and summary tables at three taxonomic levels:
#   phylum, order, and family
# Panels are arranged as: Green PRE | Green POST | Brown PRE | Brown POST
#
# Bar charts: top-N taxa shown per panel (standard for readability), with
#   remaining taxa collapsed to "Other".
# Summary tables: THRESHOLD-based, not top-N. A taxon is included if its
#   mean relative abundance reaches min_pct (default 1%) in at least one of
#   the four Species x Timepoint groups. Its real value is then reported in
#   ALL FOUR groups (including values below the threshold), so there are no
#   gaps from a taxon being "below the cutoff in this group, true value
#   unknown" the way a per-group top-N table would produce.
#
# Requires: ggplot2, dplyr, tidyr, readr, forcats
# Run after qiime2_pipeline.sh has been executed.
#
# Outputs (in 07_figures/):
#   taxa_phylum_by_species_timepoint.pdf/.png
#   taxa_order_by_species_timepoint.pdf/.png
#   taxa_family_by_species_timepoint.pdf/.png
#   taxonomy_summary_phylum.csv   (taxa >=1% in at least one group, wide format)
#   taxonomy_summary_order.csv
#   taxonomy_summary_family.csv
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(forcats)
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

# ---- Load metadata ----
meta <- read.delim("02_metadata/sample-metadata.tsv",
                   sep = "\t", comment.char = "", check.names = FALSE)
meta$SampleID  <- trimws(meta$`#SampleID`)
meta$Group     <- trimws(meta$Group)
meta$Timepoint <- trimws(meta$Timepoint)

cat("Metadata loaded:", nrow(meta), "rows\n")
cat("Unique Group values found:", paste(unique(meta$Group), collapse = ", "), "\n")
cat("Unique Timepoint values found:", paste(unique(meta$Timepoint), collapse = ", "), "\n")

# Map Group codes to readable labels (adjust if Group codes differ).
# Flags rather than silently mislabeling if Group isn't exactly "G"/"B".
if (!all(unique(meta$Group) %in% c("G", "B"))) {
  warning("meta$Group contains values other than 'G'/'B': ",
          paste(unique(meta$Group), collapse = ", "),
          ". Check the Species mapping below is still correct.")
}
meta$Species <- ifelse(meta$Group == "G", "Green anole", "Brown anole")

# Panel order: PRE before POST, Green before Brown
meta$Species   <- factor(meta$Species,   levels = c("Green anole", "Brown anole"))
meta$Timepoint <- factor(meta$Timepoint, levels = c("PRE", "POST"))

# Keep only biological samples, IF the metadata distinguishes them this way.
# (Many pipelines already exclude controls upstream -- e.g. this project's
# table-no-controls.qza -- in which case sample-metadata.tsv may only ever
# contain biological samples and this column won't exist. Only filter if
# the column is actually present, so this doesn't silently empty the table.)
if ("SampleType" %in% colnames(meta)) {
  cat("SampleType column found -- filtering to Biological samples only.\n")
  meta <- meta[meta$SampleType == "Biological", ]
} else {
  cat("No SampleType column found in metadata -- assuming all",
      nrow(meta), "rows are biological samples (controls already removed upstream).\n")
}

# ---- Parse SILVA taxonomy string ----
# Collapsed tables use full taxonomy strings as feature IDs, e.g.:
#   "d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Pseudomonadales"
# This function extracts the name at the target level.
extract_taxon_name <- function(tax_string, level_prefix) {
  # level_prefix: "p__", "o__", "f__", etc.
  parts <- strsplit(trimws(tax_string), ";")[[1]]
  parts <- trimws(parts)
  match_part <- parts[startsWith(parts, level_prefix)]
  if (length(match_part) == 0) return("Unclassified")
  name <- sub(level_prefix, "", match_part[1])
  if (name == "" || name == "uncultured" || name == "unidentified" ||
      grepl("^metagenome$", name, ignore.case = TRUE)) return("Unclassified")
  name
}

# ---- Load and process a collapsed feature table ----
load_collapsed_table <- function(tsv_path, level_prefix, level_label) {
  raw_lines <- readLines(tsv_path)

  # A standard `biom convert --to-tsv` export has a true comment line
  # ("# Constructed from biom file") followed by the real header row,
  # which ALSO starts with "#" (e.g. "#OTU ID\tSample1\tSample2...").
  # Skipping every line that starts with "#" would skip the header too,
  # turning the first row of counts into column names. Instead, find the
  # header row specifically and skip only the line(s) before it.
  header_idx <- which(grepl("^#\\s*(OTU ID|Feature ID|FeatureID)", raw_lines,
                            ignore.case = TRUE))
  if (length(header_idx) == 0) {
    # Fallback for non-standard files: assume the standard single
    # leading comment line used by biom convert.
    skip_n <- 1
    warning("Could not find a '#OTU ID' / '#Feature ID' header row in '",
            tsv_path, "'. Assuming 1 leading comment line. If sample IDs ",
            "look wrong below, check this file's format manually.")
  } else {
    skip_n <- header_idx[1] - 1
  }

  df <- read.delim(tsv_path, skip = skip_n,
                   header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  colnames(df)[1] <- "Taxon"

  # Diagnostic check: confirm sample IDs in the feature table actually
  # overlap with metadata before proceeding, so a mismatch fails clearly
  # here rather than as a cryptic error several steps downstream.
  table_samples <- setdiff(colnames(df), "Taxon")
  n_overlap <- sum(table_samples %in% meta$SampleID)
  if (n_overlap == 0) {
    stop("No sample IDs in '", tsv_path, "' match meta$SampleID.\n",
         "First few column names in the table: ",
         paste(head(table_samples, 5), collapse = ", "), "\n",
         "First few SampleIDs in metadata: ",
         paste(head(meta$SampleID, 5), collapse = ", "), "\n",
         "Check for naming mismatches (e.g. punctuation, case, or extra characters).")
  }
  cat("  ", level_label, ": ", n_overlap, " of ", length(table_samples),
      " table samples matched in metadata\n", sep = "")

  # Melt to long format
  long <- df %>%
    pivot_longer(-Taxon, names_to = "SampleID", values_to = "Count") %>%
    filter(SampleID %in% meta$SampleID)

  # Add readable taxon name
  long$TaxonName <- sapply(long$Taxon, extract_taxon_name,
                           level_prefix = level_prefix)

  # Calculate relative abundance within each sample
  long <- long %>%
    group_by(SampleID) %>%
    mutate(RelAbund = Count / sum(Count, na.rm = TRUE)) %>%
    ungroup()

  # Add metadata
  long <- long %>%
    left_join(meta[, c("SampleID", "Species", "Timepoint")],
              by = "SampleID")

  long$Level <- level_label
  long
}

# ---- Build a stacked bar chart ----
make_bar_chart <- function(long_df, level_label, top_n = 9,
                           out_base, w = 12, h = 5) {
  if (nrow(long_df) == 0) {
    stop("make_bar_chart received an empty data frame for '", level_label,
         "'. This usually means no samples matched between the feature ",
         "table and metadata -- check the diagnostic message printed by ",
         "load_collapsed_table() above.")
  }

  # Identify top N taxa by mean relative abundance across all samples
  top_taxa <- long_df %>%
    group_by(TaxonName) %>%
    summarise(MeanAbund = mean(RelAbund, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(MeanAbund)) %>%
    slice_head(n = top_n) %>%
    pull(TaxonName)

  # Collapse everything else to "Other"
  long_df$TaxonPlot <- ifelse(long_df$TaxonName %in% top_taxa,
                               long_df$TaxonName, "Other")

  # Fixed order: top taxa by abundance, then Other, then Unclassified at bottom
  taxon_order <- c(top_taxa[top_taxa != "Unclassified"],
                   "Other", "Unclassified")
  taxon_order <- taxon_order[taxon_order %in% unique(long_df$TaxonPlot)]

  long_df$TaxonPlot <- factor(long_df$TaxonPlot, levels = rev(taxon_order))

  # Color palette: distinguishable colors for up to 11 categories
  pal_base <- c(
    "#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F",
    "#EDC948", "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC",
    "#A0CBE8"
  )
  n_cats <- length(taxon_order)
  pal <- setNames(
    c(pal_base[seq_len(min(n_cats - 1, length(pal_base)))], "#CCCCCC"),
    c(taxon_order[taxon_order != "Unclassified" & taxon_order != "Other"],
      "Other")
  )
  if ("Unclassified" %in% taxon_order) pal["Unclassified"] <- "#EEEEEE"

  p <- ggplot(long_df,
              aes(x = SampleID, y = RelAbund, fill = TaxonPlot)) +
    geom_col(width = 0.9) +
    facet_grid(. ~ Species + Timepoint, scales = "free_x", space = "free_x",
               labeller = label_wrap_gen(multi_line = TRUE)) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                       expand = c(0, 0)) +
    scale_fill_manual(values = pal, name = level_label,
                      guide = guide_legend(reverse = TRUE,
                                           ncol = 1)) +
    theme_classic(base_size = 10) +
    theme(
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank(),
      strip.background = element_blank(),
      strip.text  = element_text(face = "bold", size = 9),
      legend.text = element_text(face = "italic"),
      legend.key.size = unit(0.4, "cm"),
      panel.spacing = unit(0.3, "lines")
    ) +
    labs(x = "Individual sample", y = "Relative abundance",
         title = paste("Bacterial", tolower(level_label), "composition"))

  ggsave(paste0(out_base, ".pdf"), p, width = w, height = h)
  ggsave(paste0(out_base, ".png"), p, width = w, height = h, dpi = 300)
  cat("Saved:", out_base, ".pdf/.png\n")
  invisible(p)
}

# ---- Build a summary table (threshold-based, not top-N) ----
# Keeps every taxon whose mean relative abundance reaches `min_pct` in AT
# LEAST ONE of the four Species x Timepoint groups, then reports that
# taxon's real value in ALL FOUR groups (including small values below the
# threshold) -- so there are no "below the cutoff, true value unknown" gaps
# the way there would be with a per-group top-N table. Output is wide
# format (one row per taxon, one column per group) for direct use as a
# manuscript table.
make_summary_table <- function(long_df, level_label, min_pct = 1.0) {
  group_means <- long_df %>%
    group_by(Species, Timepoint, TaxonName) %>%
    summarise(MeanPct = 100 * mean(RelAbund, na.rm = TRUE), .groups = "drop")

  # Taxa qualifying in at least one group
  qualifying_taxa <- group_means %>%
    filter(MeanPct >= min_pct) %>%
    pull(TaxonName) %>%
    unique()

  # Report the real value for every qualifying taxon in every group,
  # filling in 0 only for true absence (taxon never observed in that
  # group), not for "didn't make the cutoff."
  wide <- group_means %>%
    filter(TaxonName %in% qualifying_taxa) %>%
    mutate(GroupLabel = paste(Species, Timepoint)) %>%
    select(TaxonName, GroupLabel, MeanPct) %>%
    pivot_wider(names_from = GroupLabel, values_from = MeanPct, values_fill = 0)

  # Sort rows by overall mean abundance across groups, descending
  group_cols <- setdiff(colnames(wide), "TaxonName")
  wide$OverallMean <- rowMeans(wide[, group_cols], na.rm = TRUE)
  wide <- wide %>% arrange(desc(OverallMean))

  # Round for display, keep OverallMean for sorting/reference, add Level
  wide <- wide %>%
    mutate(across(all_of(group_cols), ~ round(.x, 1)),
           OverallMean = round(OverallMean, 1),
           Level = level_label) %>%
    relocate(Level, .before = TaxonName)

  cat("  ", level_label, ": ", nrow(wide), " taxa at >=", min_pct,
      "% in at least one group\n", sep = "")

  wide
}

# =============================================================================
# MAIN ANALYSIS
# =============================================================================

all_tables <- list()

# ---- Phylum ----
cat("\n--- Phylum level ---\n")
phylum_long <- load_collapsed_table(
  "10_exports/exported_phylum_level/feature-table.tsv",
  level_prefix = "p__",
  level_label  = "Phylum"
)
make_bar_chart(phylum_long, "Phylum", top_n = 8,
               out_base = "07_figures/taxa_phylum_by_species_timepoint")
all_tables[["Phylum"]] <- make_summary_table(phylum_long, "Phylum", min_pct = 1.0)

# ---- Order ----
cat("\n--- Order level ---\n")
order_long <- load_collapsed_table(
  "10_exports/exported_order_level/feature-table.tsv",
  level_prefix = "o__",
  level_label  = "Order"
)
make_bar_chart(order_long, "Order", top_n = 9,
               out_base = "07_figures/taxa_order_by_species_timepoint")
all_tables[["Order"]] <- make_summary_table(order_long, "Order", min_pct = 1.0)

# ---- Family ----
cat("\n--- Family level ---\n")
family_long <- load_collapsed_table(
  "10_exports/exported_family_level/feature-table.tsv",
  level_prefix = "f__",
  level_label  = "Family"
)
make_bar_chart(family_long, "Family", top_n = 9,
               out_base = "07_figures/taxa_family_by_species_timepoint")
all_tables[["Family"]] <- make_summary_table(family_long, "Family", min_pct = 1.0)

# ---- Write one summary CSV per taxonomic level ----
# Each level becomes its own table since that's how they'll be used in
# the manuscript (e.g. "Order-level table" next to "Family-level Table 1"),
# rather than one combined file that would need re-splitting later.
dir.create("07_figures", showWarnings = FALSE, recursive = TRUE)
for (lev in names(all_tables)) {
  out_path <- paste0("07_figures/taxonomy_summary_", tolower(lev), ".csv")
  write_csv(all_tables[[lev]], out_path)
  cat("Saved:", out_path, "\n")
}

# Print tables to console for quick review
for (lev in c("Phylum", "Order", "Family")) {
  cat("\n===", lev, "level: taxa >=1% in at least one group ===\n")
  print(as.data.frame(all_tables[[lev]]), row.names = FALSE)
}

cat("\nAll taxonomy analyses complete.\n")
