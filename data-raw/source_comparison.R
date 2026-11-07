# Comparison of sources of interactions -> data/source_comparison.rda
#
# Built from the results of dev/validation/compare_sources.R: cross-validation
# by system on the 12 validation webs of gray2015 without study 12 (Tuesday
# Lake) with the Gray registry, GloBI and both, and the final test on Tuesday
# Lake 1984 and 1986. The GloBI download (dev/globi/download_globi.R,
# processed by dev/globi/process_globi.R) is not bundled; only the scores
# are. Packages used only here: pkgload, usethis.
#
# Run from the package root after dev/validation/compare_sources.R:
#   Rscript data-raw/source_comparison.R

pkgload::load_all(".", quiet = TRUE)
cmp <- readRDS(file.path("dev", "validation", "source_comparison.rds"))
tl <- readRDS(file.path("dev", "validation", "tuesday_lake_sources.rds"))$scores
# how each source was used: registries, primary registry (NA: every registry
# at every level) and highest level of the supplementary registry
setup <- data.frame(
  key = c("gray", "globi", "both_family", "both_class", "both_all"),
  source = c("Gray", "GloBI", "Gray + GloBI", "Gray + GloBI", "Gray + GloBI"),
  primary_registry = c(NA, NA, "gray2015", "gray2015", NA),
  supplementary_max_level = c(NA, NA, "family", "class", NA),
  stringsAsFactors = FALSE
)
cv <- cmp$table[c("source", "rule", "brier", "brier_se", "brier_diff", "brier_diff_se",
                  "diff_vs_gray", "diff_vs_gray_se", "skill", "tss", "precision", "recall",
                  "n_predicted")]
key <- paste(cv$source, cv$rule)
tl_cols <- c("brier", "skill", "tss", "precision", "recall")
tl_key <- paste(tl$source, tl$rule)
tl_mean <- as.data.frame(lapply(stats::setNames(tl_cols, paste0("tl_", tl_cols)), function(col) {
  as.numeric(tapply(tl[[col]], tl_key, mean, na.rm = TRUE)[key])
}))
s <- setup[match(cv$source, setup$key), ]
source_comparison <- cbind(source = s$source, primary_registry = s$primary_registry,
                           supplementary_max_level = s$supplementary_max_level,
                           cv[-1], tl_mean, stringsAsFactors = FALSE)
source_comparison$tl_precision[is.nan(source_comparison$tl_precision)] <- NA_real_
rownames(source_comparison) <- NULL
print(source_comparison, digits = 3)
usethis::use_data(source_comparison, overwrite = TRUE, compress = "xz")
