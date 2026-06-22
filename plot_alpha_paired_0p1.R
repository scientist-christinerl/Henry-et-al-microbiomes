#!/usr/bin/env Rscript
# =============================================================================
# plot_alpha_paired_0p1.R
# Lattin Lab | LSU
#
# Paired alpha diversity figure (PRE vs POST) for the three metrics at the
# 0.1 decontam threshold: observed ASVs, Shannon, Faith's PD. Within-lizard
# repeated measures are connected by thin lines; group means with 95% CI are
# overlaid in black.
#
# Outputs (in 07_figures/):
#   alpha_diversity_paired_0p1.pdf/.png
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(tidyr)
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
  "Shannon",  "10_exports/exported_shannon_0p1/alpha-diversity.tsv",
  "FaithPD",  "10_exports/exported_faith_pd_0p1/alpha-diversity.tsv",
  "Observed", "10_exports/exported_observed_0p1/alpha-diversity.tsv"
)

out_pdf <- "07_figures/alpha_diversity_paired_0p1.pdf"
out_png <- "07_figures/alpha_diversity_paired_0p1.png"

# Species colors
species_colors <- c(
  "Brown anole" = "#8B4513",
  "Green anole" = "#2E8B57"
)

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

# Map coded values (B/G and PRE/POST) to publication labels. Values that are
# already in label form are left unchanged.
meta <- meta %>%
  mutate(
    Group = case_when(
      Group %in% c("B", "Brown", "brown") ~ "Brown anole",
      Group %in% c("G", "Green", "green") ~ "Green anole",
      Group %in% c("Brown anole", "Green anole") ~ Group,
      TRUE ~ Group
    ),
    Timepoint = case_when(
      Timepoint %in% c("PRE", "Pre", "pre") ~ "PRE",
      Timepoint %in% c("POST", "Post", "post") ~ "POST",
      TRUE ~ Timepoint
    )
  )

meta$Group <- factor(meta$Group, levels = c("Brown anole", "Green anole"))
meta$Timepoint <- factor(meta$Timepoint, levels = c("PRE", "POST"))
meta$Subject <- factor(meta$Subject)

# -----------------------------
# 2) Read + stack alpha diversity vectors
# -----------------------------
read_alpha <- function(metric, path) {
  if (!file.exists(path)) stop("Missing alpha file: ", path)

  a <- read_tsv(path, show_col_types = FALSE, comment = "")
  colnames(a)[1] <- "SampleID"
  colnames(a)[2] <- "Value"

  a %>%
    mutate(
      SampleID = trimws(SampleID),
      Metric = metric
    )
}

alpha_long <- purrr::map2_dfr(alpha_files$Metric, alpha_files$Path, read_alpha)

dat <- alpha_long %>%
  left_join(meta, by = "SampleID") %>%
  filter(!is.na(Group), !is.na(Timepoint), !is.na(Subject)) %>%
  mutate(Metric = factor(Metric, levels = c("Observed", "Shannon", "FaithPD")))

# -----------------------------
# 3) Publication-quality paired plot
# -----------------------------
# Summary points (mean +/- 95% CI) per Group x Timepoint x Metric
sum_df <- dat %>%
  group_by(Metric, Group, Timepoint) %>%
  summarise(
    mean = mean(Value),
    se = sd(Value) / sqrt(n()),
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(
    ci = 1.96 * se,
    ymin = mean - ci,
    ymax = mean + ci
  )

p <- ggplot(dat, aes(x = Timepoint, y = Value, color = Group)) +
  facet_wrap(~Metric, scales = "free_y", nrow = 1) +
  theme_classic(base_size = 12) +
  labs(
    x = "",
    y = "Alpha diversity",
    color = "Species",
    title = "Alpha diversity (paired PRE vs POST)",
    subtitle = "Lines connect repeated measures within each lizard (Subject)."
  ) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  ) +
  scale_color_manual(values = species_colors)

# Paired lines (thin, unobtrusive)
p <- p +
  geom_line(
    aes(group = Subject),
    linewidth = 0.35,
    alpha = 0.35,
    show.legend = FALSE
  )

# Points (jittered a bit)
p <- p +
  geom_point(
    position = position_jitter(width = 0.08, height = 0, seed = 1),
    size = 2.2,
    alpha = 0.85
  )

# Overlay mean +/- 95% CI (black so it's visible on both colors)
p <- p +
  geom_errorbar(
    data = sum_df,
    aes(x = Timepoint, ymin = ymin, ymax = ymax, group = Group),
    inherit.aes = FALSE,
    width = 0.12,
    linewidth = 0.6,
    color = "black"
  ) +
  geom_point(
    data = sum_df,
    aes(x = Timepoint, y = mean),
    inherit.aes = FALSE,
    size = 2.4,
    shape = 21,
    fill = "white",
    color = "black"
  )

# -----------------------------
# 4) Save
# -----------------------------
ggsave(out_pdf, p, width = 10.5, height = 4.2)
ggsave(out_png, p, width = 10.5, height = 4.2, dpi = 300)

cat("Wrote:\n- ", out_pdf, "\n- ", out_png, "\n", sep = "")
