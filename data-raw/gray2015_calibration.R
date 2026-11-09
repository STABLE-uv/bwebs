# Calibration of link probabilities for the bundled registry ->
# data/gray2015_calibration.rda
#
# Source: data/gray2015.rda (see data-raw/gray2015.R) and the provisional
# validation set of data/gray2015_sources.rda (column `validation`), all 13
# webs: Tuesday Lake (study 12) was kept apart only while the matching rule
# was chosen and tested (dev/validation/compare_rules.R, tuesday_lake.R).
# Rule: matching_rule() with its defaults (the combination with the lowest
# cross-validated Brier score, dev/validation/compare_sources.R; with a single
# registry, using every registry at every level is the same as naming it the
# primary one), with the taxonomic decay of the pairs not inferred. Bayesian
# fit, 4 chains x 4000 draws, seed 2026 (the decay model keeps 4000 draws).
# Takes about 15 minutes. Packages used only here: pkgload, usethis.
#
# Run from the package root after data-raw/gray2015.R:
#   Rscript data-raw/gray2015_calibration.R

pkgload::load_all(".", quiet = TRUE)
sources <- gray2015_sources$study_id[gray2015_sources$validation]
# Studies left out together because they describe the same system: 62% of the
# links of study 13 (Ledger, Edwards & Woodward, unpublished, in Brose et al.
# 2005) are also recorded in study 50 (Ledger et al. 2012); no other study
# holds more than 9% of the links of a validation web (shared_studies()).
# Provisional: 46 and 49 (Woodward, Hildrew et al.) and 15 and 59 (P. H.
# Warren) may also describe the same system; see dev/validation/results.md.
groups <- list(c("13", "50"))
validation <- validate_rule(gray2015, sources = sources, groups = groups, rule = matching_rule())
gray2015_calibration <- calibrate_links(validation, method = "bayes", n_chains = 4L,
                                        n_iter = 4000L, n_warmup = 2000L, quadrature = 20L,
                                        seed = 2026)
print(gray2015_calibration)
usethis::use_data(gray2015_calibration, overwrite = TRUE, compress = "xz")
