# Final test on Tuesday Lake (TL84 and TL86, from cheddar), independent of the
# selection and of the calibration: the registry and the 12 validation webs do
# not contain study 12 of Gray (which is Tuesday Lake).
#
# For each rule of the comparison (dev/validation/compare_rules.R):
#   1. Bayesian calibration with its 12 validation webs;
#   2. metaweb of the 82 taxa of Tuesday Lake with the registry without study 12;
#   3. Brier score, TSS, precision and recall against the observed web of each
#      year (only the taxa present that year).
# For the chosen rule and the rule of the reanalysis: indicators of 1000
# sampled webs against those of the observed web.
#
# Run from the root of the package: Rscript dev/validation/tuesday_lake.R
pkgload::load_all(".", quiet = TRUE)
cmp <- readRDS(file.path("dev", "validation", "rule_comparison_gray.rds"))
reg <- gray2015[gray2015$study_id != "12", ]
years <- c("TL84", "TL86")
present <- split(tuesday_lake_community$node_id, tuesday_lake_community$sample_id)

cals <- mws <- list()
scores <- list()
for (nm in cmp$table$rule) {
  v <- cmp$validations[[nm]]
  cals[[nm]] <- calibrate_links(v, seed = 2026, quiet = TRUE)
  mws[[nm]] <- build_metaweb(prepare_nodes(tuesday_lake_taxa, rule = v$settings$rule), reg)
  for (y in years) {
    obs <- tuesday_lake_links[tuesday_lake_links$sample_id == y, ]
    s <- score_metaweb(mws[[nm]], cals[[nm]], obs, nodes = present[[y]])$summary
    scores[[length(scores) + 1L]] <- cbind(rule = nm, year = y, s, stringsAsFactors = FALSE)
  }
  message(nm)
}
tl <- do.call(rbind, scores)
tl_mean <- stats::aggregate(cbind(brier, skill, tss, precision, recall) ~ rule, tl, mean,
                            na.rm = TRUE, na.action = stats::na.pass)
both <- merge(cmp$table[c("rule", "brier", "brier_diff", "brier_diff_se", "skill", "tss")],
              tl_mean, by = "rule", suffixes = c("_cv", "_tl"))
both <- both[order(both$brier_cv), ]
both$rank_cv <- seq_len(nrow(both))
both$rank_tl <- rank(both$brier_tl)
print(both, digits = 3, row.names = FALSE)
cat("Spearman (Brier CV vs Tuesday Lake):",
    round(stats::cor(both$brier_cv, both$brier_tl, method = "spearman"), 2), "\n")

# Indicators: chosen rule (best cross-validated Brier score) and rule of the
# reanalysis
best <- cmp$table$rule[1]
reanalysis <- "climb/family/class/1"
observed_ind <- do.call(rbind, lapply(years, function(y) {
  ids <- present[[y]]
  A <- matrix(0, length(ids), length(ids), dimnames = list(ids, ids))
  l <- tuesday_lake_links[tuesday_lake_links$sample_id == y, ]
  A[cbind(l$resource, l$consumer)] <- 1
  cbind(sample_id = y, web_indicators(A))
}))
ind_cols <- c("S", "L", "LS", "C", "prop_basal", "prop_top", "TL_mean", "prop_omnivores",
              "MxSim_trophic_species")
compare_ind <- do.call(rbind, lapply(unique(c(best, reanalysis)), function(nm) {
  draws <- sample_metawebs(mws[[nm]], cals[[nm]], n_draws = 1000, seed = 2026)
  ind <- calc_indicators(draws, tuesday_lake_community)
  s <- summarize_indicators(ind, indicators = ind_cols, na_rm = TRUE)
  s$observed <- vapply(seq_len(nrow(s)), function(i) {
    observed_ind[observed_ind$sample_id == s$sample_id[i], s$indicator[i]]
  }, numeric(1))
  s$covered <- s$observed >= s$lower & s$observed <= s$upper
  cbind(rule = nm, s[c("sample_id", "indicator", "estimate", "lower", "upper", "observed",
                       "covered", "n_omitted")])
}))
print(compare_ind, digits = 3, row.names = FALSE)
saveRDS(list(scores = tl, summary = both, indicators = compare_ind, observed = observed_ind,
             calibrations = cals[unique(c(best, reanalysis))]),
        file.path("dev", "validation", "tuesday_lake.rds"))
