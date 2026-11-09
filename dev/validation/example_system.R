# Example of the README with the bundled calibration: unknown system
# (system = "random", the default) against a typical system ("average"),
# 1000 draws. Run from the root of the package.
pkgload::load_all(".", quiet = TRUE)
# Without the taxonomic decay (matching_rule(decay = FALSE)), as when these
# results were computed (before version 0.3.0).
nodes <- prepare_nodes(example_taxa, rule = matching_rule(decay = FALSE))
mw <- build_metaweb(nodes, gray2015)
for (system in c("random", "average")) {
  dr <- sample_metawebs(mw, gray2015_calibration, n_draws = 1000, system = system, seed = 2026)
  ind <- calc_indicators(dr, example_community, permanent = c("cpom", "fpom"))
  s <- summarize_indicators(ind, indicators = c("S", "C", "TL_mean"), samples = example_samples,
                            by = "treatment", na_rm = TRUE)
  k <- contrast_indicators(ind, example_samples, by = "treatment", reference = "control",
                           indicators = c("S", "C", "TL_mean"), na_rm = TRUE)
  cat("\nsystem =", system, "; median links per draw:", stats::median(colSums(dr$z)), "\n")
  print(s[, c("treatment", "indicator", "estimate", "lower", "upper")], digits = 3, row.names = FALSE)
  print(k[, c("indicator", "estimate", "lower", "upper", "prob_positive")], digits = 3,
        row.names = FALSE)
}
