# Comparison of matching rules -> data/gray2015_rule_comparison.rda
#
# Built from the results of dev/validation/compare_rules.R (cross-validation
# by system on the 12 validation webs of gray2015 without study 12, Tuesday
# Lake; about 30 minutes in two processes) and dev/validation/tuesday_lake.R
# (final test on Tuesday Lake 1984 and 1986). Packages used only here:
# pkgload, usethis.
#
# Run from the package root after those two scripts:
#   Rscript data-raw/gray2015_rule_comparison.R

pkgload::load_all(".", quiet = TRUE)
cmp <- readRDS(file.path("dev", "validation", "rule_comparison_gray.rds"))
tl <- readRDS(file.path("dev", "validation", "tuesday_lake.rds"))$scores
cv <- cmp$table[c("rule", "mode", "base_level", "climb_max", "min_partners", "downward",
                  "brier", "brier_se", "brier_diff", "brier_diff_se", "skill", "tss",
                  "precision", "recall", "n_predicted")]
tl_cols <- c("brier", "skill", "tss", "precision", "recall")
tl_mean <- as.data.frame(lapply(stats::setNames(tl_cols, paste0("tl_", tl_cols)), function(col) {
  as.numeric(tapply(tl[[col]], tl$rule, mean, na.rm = TRUE)[cv$rule])
}))
gray2015_rule_comparison <- cbind(cv, tl_mean)
gray2015_rule_comparison$tl_precision[is.nan(gray2015_rule_comparison$tl_precision)] <- NA_real_
rownames(gray2015_rule_comparison) <- NULL
print(gray2015_rule_comparison, digits = 3)
usethis::use_data(gray2015_rule_comparison, overwrite = TRUE, compress = "xz")
