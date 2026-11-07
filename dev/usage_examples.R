# Complete workflow of a user of bwebs (run in a clean session with
# Rscript --vanilla in an empty folder; it is the basis of the README, the
# vignette and tests/testthat/test-workflow.R).
library(bwebs)

# 1. Nodes of the community and their matching ladder
nodes <- prepare_nodes(example_taxa)
nodes

# 2. Deterministic metaweb with the evidence of each link
mw <- build_metaweb(nodes, gray2015)
mw
table(mw$edges$evidence_class)

# 3. Bayesian calibration of the Gray registry (bundled; computed again with
#    calibrate_links(validate_rule(gray2015, sources = ...)))
gray2015_calibration

# 3b. How the default rule was chosen (cross-validation by system and test on
#     Tuesday Lake)
head(gray2015_rule_comparison[, c("rule", "brier", "brier_diff", "brier_diff_se", "tl_brier")])

# 4. Plausible metawebs
draws <- sample_metawebs(mw, gray2015_calibration, n_draws = 50, seed = 2026)
draws

# 5. Indicators of each sample in each draw (detritus is not sampled: it is a
#    permanent resource)
ind <- calc_indicators(draws, example_community, permanent = c("cpom", "fpom"))
ind

# 6. Summaries by sample and by treatment
summarize_indicators(ind, indicators = c("S", "C", "TL_mean"))
by_treatment <- summarize_indicators(ind, indicators = c("S", "C", "TL_mean"),
                                     samples = example_samples, by = "treatment")
by_treatment

# 7. Deterministic reference (the whole metaweb)
det <- calc_indicators(mw, example_community, permanent = c("cpom", "fpom"))
summarize_indicators(det, indicators = "C", samples = example_samples, by = "treatment")

# 8. Export
out <- tempfile(fileext = ".csv")
utils::write.csv(by_treatment, out, row.names = FALSE)
file.exists(out)
