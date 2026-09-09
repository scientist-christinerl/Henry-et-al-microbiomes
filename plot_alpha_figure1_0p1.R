#!/usr/bin/env Rscript
# =============================================================================
# plot_alpha_figure1_0p1.R
# Lattin Lab | LSU
#
# Reproduces manuscript Figure 1: alpha diversity (observed ASV richness,
# Shannon, Faith's PD) by species, with pre- and post-captivity samples shown
# separately. Species is on the x-axis; timepoint is encoded by point shape
# (PRE = triangle, POST = circle) and separated by dodging within each species.
# Colored diamonds with error bars show group means +/- standard error. There
# are no connecting lines between individuals.
#
# This is the grouped-by-species design used in the manuscript, and is a
# separate figure from plot_alpha_paired_0p1.R (which plots PRE vs POST on the
# x-axis with per-individual connecting lines).
#
# Reads the same three alpha diversity exports as the other alpha scripts.
#
# Outputs (in 07_figures/):
#   figure1_alpha_diversity.pdf   (vector, for submission)
#   figure1_alpha_diversity.png   (300 dpi raster)
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(purrr)
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

# ---- Inputs ----
meta_file <- "02_metadata/sample-metadata-with-subject.tsv"

alpha_files <- tibble::tribble(
  ~Metric,    ~Path,
  "Observed", "10_exports/exported_observed_0p1/alpha-diversity.tsv",
  "Shannon",  "10_exports/exported_shannon_0p1/alpha-diversity.tsv",
  "FaithPD",  "10_exports/exported_faith_pd_0p1/alpha-diversity.tsv"
)

out_pdf <- "07_figures/figure1_alpha_diversity.pdf"
out_png <- "07_figures/figure1_alpha_diversity.png"

# Species colors (keyed to the short Group labels used on the axis)
group_colors <- c(
  "Brown" = "#A0522D",
  "Green" = "#2E8B57"
)

# Timepoint shapes: PRE = filled triangle, POST = filled circle
timepoint_shapes <- c("PRE" = 17, "POST" = 16)

# -----------------------------
# 1) Read metadata
# -----------------------------
meta <- read_tsv(meta_file, show_col_types = FALSE) %>%
  mutate(
    SampleID = trimws(`#SampleID`),
    Group = trimws(Group),
    Timepoint = trimws(Timepoint),
    Subject = trimws(Subject)
  )

# Map coded values to the short labels used in the figure. Values already in
# label form are left unchanged.
meta <- meta %>%
  mutate(
    Group = case_when(
      Group %in% c("B", "Brown", "brown", "Brown anole") ~ "Brown",
      Group %in% c("G", "Green", "green", "Green anole") ~ "Green",
      TRUE ~ Group
    ),
    Timepoint = case_when(
      Timepoint %in% c("PRE", "Pre", "pre") ~ "PRE",
      Timepoint %in% c("POST", "Post", "post") ~ "POST",
      TRUE ~ Timepoint
    )
  )

meta$Group <- factor(meta$Group, levels = c("Brown", "Green"))
meta$Timepoint <- factor(meta$Timepoint, levels = c("PRE", "POST"))

# -----------------------------
# 2) Read + stack alpha diversity vectors
# -----------------------------
read_alpha <- function(metric, path) {
  if (!file.exists(path)) stop("Missing alpha file: ", path)
  a <- read_tsv(path, show_col_types = FALSE, comment = "")
  colnames(a)[1] <- "SampleID"
  colnames(a)[2] <- "Value"
  a %>% mutate(SampleID = trimws(SampleID), Metric = metric)
}

alpha_long <- map2_dfr(alpha_files$Metric, alpha_files$Path, read_alpha)

dat <- alpha_long %>%
  left_join(meta, by = "SampleID") %>%
  filter(!is.na(Group), !is.na(Timepoint)) %>%
  mutate(Metric = factor(Metric, levels = c("Observed", "Shannon", "FaithPD")))

# -----------------------------
# 3) Build the figure
# -----------------------------
dodge_w <- 0.7
pd <- position_dodge(width = dodge_w)

p <- ggplot(dat, aes(x = Group, y = Value, color = Group, shape = Timepoint)) +
  facet_wrap(~ Metric, scales = "free_y") +
  # Individual samples, jittered within each species x timepoint cluster
  geom_point(
    position = position_jitterdodge(jitter.width = 0.15,
                                    dodge.width = dodge_w,
                                    seed = 1),
    size = 2, alpha = 0.9
  ) +
  # Group means (diamonds) with mean +/- 1 SE error bars
  stat_summary(fun.data = mean_se, geom = "errorbar",
               width = 0.15, linewidth = 0.6, position = pd,
               show.legend = FALSE) +
  stat_summary(fun = mean, geom = "point",
               shape = 18, size = 3.4, position = pd,
               show.legend = FALSE) +
  scale_color_manual(values = group_colors) +
  scale_shape_manual(values = timepoint_shapes) +
  labs(x = "Group", y = "Diversity", color = "Group", shape = "Timepoint") +
  guides(
    color = guide_legend(order = 1,
                         override.aes = list(shape = 18, size = 3, linetype = 0)),
    shape = guide_legend(order = 2)
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

# -----------------------------
# 4) Save at submission resolution
# -----------------------------
ggsave(out_pdf, p, width = 9.5, height = 3.8)
ggsave(out_png, p, width = 9.5, height = 3.8, dpi = 300)

cat("Wrote:\n- ", out_pdf, "\n- ", out_png, "\n", sep = "")
