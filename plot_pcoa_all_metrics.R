#!/usr/bin/env Rscript
# =============================================================================
# plot_pcoa_all_metrics.R
# Lattin Lab | LSU
#
# PCoA plots for all four beta diversity metrics:
#   Bray-Curtis, Weighted UniFrac, Unweighted UniFrac, Jaccard
#
# All plots use PRE-only samples (between-species comparison at capture).
#
# Outputs (in 07_figures/):
#   bray_pcoa_0p1.pdf/.png
#   wuf_pcoa_0p1.pdf/.png
#   uuf_pcoa_0p1.pdf/.png
#   jaccard_pcoa_0p1.pdf/.png
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
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

# ---- Colors (green anole = green, brown anole = brown) ----
COLORS <- c(G = "#2E8B57", B = "#8B4513")
LABELS <- c(G = "Green anole", B = "Brown anole")

# ---- Load metadata ----
meta <- read.delim("02_metadata/sample-metadata.tsv",
                   sep = "\t", comment.char = "", check.names = FALSE)
meta$SampleID  <- trimws(meta$`#SampleID`)
meta$Group     <- trimws(meta$Group)
meta$Timepoint <- trimws(meta$Timepoint)

cat("Unique Group values found:", paste(unique(meta$Group), collapse = ", "), "\n")

# Keep PRE samples, restricting to biological samples only IF the metadata
# distinguishes them this way (some pipelines already exclude controls
# upstream, in which case this column won't exist -- only filter if present
# so this doesn't silently produce an empty data frame).
if ("SampleType" %in% colnames(meta)) {
  meta_pre <- meta[meta$Timepoint == "PRE" & meta$SampleType == "Biological", ]
} else {
  meta_pre <- meta[meta$Timepoint == "PRE", ]
}
cat("PRE samples for PCoA plots:", nrow(meta_pre), "\n")

# ---- PCoA function ----
make_pcoa_plot <- function(dist_path, metric_label, out_base) {

  if (!file.exists(dist_path)) {
    cat("File not found, skipping:", dist_path, "\n")
    return(invisible(NULL))
  }

  # Load full distance matrix, then subset to PRE samples
  dm <- read.delim(dist_path, header = TRUE, row.names = 1, check.names = FALSE)
  dm <- as.matrix(dm)
  mode(dm) <- "numeric"

  pre_ids <- intersect(meta_pre$SampleID, rownames(dm))
  dm_pre  <- dm[pre_ids, pre_ids]
  dist_obj <- as.dist(dm_pre)

  # PCoA via cmdscale
  pcoa <- cmdscale(dist_obj, k = 2, eig = TRUE)
  scores <- as.data.frame(pcoa$points)
  colnames(scores) <- c("PC1", "PC2")
  scores$SampleID <- rownames(scores)

  # Variance explained (only positive eigenvalues count)
  eigvals  <- pcoa$eig
  pos_eig  <- eigvals[eigvals > 0]
  var_expl <- eigvals / sum(pos_eig)
  pc1_lab  <- sprintf("PC1 (%.1f%%)", 100 * var_expl[1])
  pc2_lab  <- sprintf("PC2 (%.1f%%)", 100 * var_expl[2])

  # Merge with metadata
  dat <- left_join(scores, meta_pre[, c("SampleID", "Group")], by = "SampleID")
  dat$Group <- factor(dat$Group, levels = names(COLORS))

  n_per_group <- table(dat$Group)
  cat(metric_label, ": n per group:", paste(names(n_per_group),
      n_per_group, sep = "=", collapse = ", "), "\n")

  p <- ggplot(dat, aes(PC1, PC2, color = Group, fill = Group)) +
    geom_point(size = 3, alpha = 0.9, shape = 21, color = "white", stroke = 0.3) +
    stat_ellipse(level = 0.95, linewidth = 0.8, geom = "polygon",
                 alpha = 0.12, linetype = "solid") +
    stat_ellipse(level = 0.95, linewidth = 0.8) +
    scale_color_manual(values = COLORS, labels = LABELS, name = NULL) +
    scale_fill_manual(values  = COLORS, labels = LABELS, name = NULL) +
    theme_classic(base_size = 12) +
    theme(
      legend.position  = "bottom",
      legend.text      = element_text(face = "italic"),
      panel.border     = element_rect(fill = NA, color = "grey70", linewidth = 0.5)
    ) +
    labs(x = pc1_lab, y = pc2_lab,
         title = paste0("Beta diversity: ", metric_label,
                        " (pre-captivity samples)"))

  ggsave(paste0(out_base, ".pdf"), p, width = 5.5, height = 5)
  ggsave(paste0(out_base, ".png"), p, width = 5.5, height = 5, dpi = 300)
  cat("  Saved:", out_base, ".pdf/.png\n")
  invisible(p)
}

# ---- Generate all four plots ----
make_pcoa_plot(
  "10_exports/exported_bray_full_0p1/distance-matrix.tsv",
  "Bray-Curtis dissimilarity",
  "07_figures/bray_pcoa_0p1"
)

make_pcoa_plot(
  "10_exports/exported_uuf_full_0p1/distance-matrix.tsv",
  "Unweighted UniFrac distance",
  "07_figures/uuf_pcoa_0p1"
)

make_pcoa_plot(
  "10_exports/exported_wuf_full_0p1/distance-matrix.tsv",
  "Weighted UniFrac distance",
  "07_figures/wuf_pcoa_0p1"
)

# Jaccard
make_pcoa_plot(
  "10_exports/exported_jaccard_full_0p1/distance-matrix.tsv",
  "Jaccard dissimilarity",
  "07_figures/jaccard_pcoa_0p1"
)

cat("\nAll PCoA plots complete.\n")
