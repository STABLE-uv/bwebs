# Calibration of link probabilities for the Gray et al. (2015) registry
# combined with GloBI -> data/gray2015_globi_calibration.rda
#
# Source: data/gray2015.rda and the GloBI download for the consumers of the
# validation webs, of Tuesday Lake and of the azoxystrobin experiment
# (dev/globi/download_globi.R, 832,003 records downloaded on 2026-11-03,
# processed by dev/globi/process_globi.R). The GloBI records are not bundled
# (their licences depend on each source); only the counts of the validation
# and the posterior draws are.
#
# All 13 validation webs, 13 and 50 left out together; when each web is
# rebuilt, the GloBI studies holding at least 10% of its links are left out
# too (copies of the web in other compilations), and for Tuesday Lake (study
# 12) also the GloBI studies that copy it (dev/globi/globi_tuesday_lake.csv).
# Rule: matching_rule() with its defaults (GloBI used at every level). Bayesian
# fit, 4 chains x 4000 draws, seed 2026. Takes about 15 minutes. Packages used
# only here: pkgload, usethis.
#
# Run from the package root after dev/globi/process_globi.R:
#   Rscript data-raw/gray2015_globi_calibration.R <path>/globi_registry.rds

pkgload::load_all(".", quiet = TRUE)
g <- readRDS(commandArgs(trailingOnly = TRUE)[1])
both <- combine_registries(gray2015, g$globi)
sources <- gray2015_sources$study_id[gray2015_sources$validation]
groups <- list(c("13", "50"))
shared <- shared_studies(both, sources, threshold = 0.1)
exclusions <- rbind(
  shared[shared$study_id %in% g$globi$study_id, c("web", "study_id")],
  data.frame(web = "12", study_id = g$tuesday_lake$study_id, stringsAsFactors = FALSE)
)
exclusions <- exclusions[!duplicated(exclusions), , drop = FALSE]
rownames(exclusions) <- NULL
print(table(exclusions$web))
validation <- validate_rule(both, sources = sources, groups = groups, rule = matching_rule(),
                            exclusions = exclusions)
gray2015_globi_calibration <- calibrate_links(validation, method = "bayes", n_chains = 4L,
                                              n_iter = 4000L, n_warmup = 2000L,
                                              quadrature = 20L, seed = 2026)
print(gray2015_globi_calibration)
usethis::use_data(gray2015_globi_calibration, overwrite = TRUE, compress = "xz")
