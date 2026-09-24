# sensitivity_analysis.R
# Lattin Lab | LSU
#
# Sensitivity analysis for the microbiome comparisons: the smallest effect
# detectable at 80% power given the achieved sample sizes. Reports Cohen's d
# for the between-species comparison and dz for the within-individual
# (pre- vs post-captivity) comparison. alpha = 0.05, two-sided.

library(pwr)

# Between-species comparison: 15 green vs 13 brown (two independent groups)
between_species <- pwr.t2n.test(n1 = 15, n2 = 13,
                                sig.level = 0.05, power = 0.80)$d

# Within-individual captivity comparison: 28 individuals measured pre and post
within_captivity <- pwr.t.test(n = 28, type = "paired",
                               sig.level = 0.05, power = 0.80)$d

cat("Minimum detectable between-species effect (Cohen's d):",
    round(between_species, 2), "\n")
cat("Minimum detectable within-individual effect (dz):     ",
    round(within_captivity, 2), "\n")