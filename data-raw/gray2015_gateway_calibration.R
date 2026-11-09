# Calibration of link probabilities for the Gray et al. (2015) registry
# combined with GATEWAy -> data/gray2015_gateway_calibration.rda
#
# Source: data/gray2015.rda and data/gateway.rda (the freshwater links of
# GATEWAy version 1.0, Brose et al. 2019, CC BY 4.0; data-raw/gateway.R).
#
# All 13 validation webs, 13 and 50 left out together; when each web is
# rebuilt, the webs of GATEWAy that copy it are left out too
# (dev/gateway/gateway_exclusions.csv: webs holding at least 10% of its
# links, alone or pooled by publication; GATEWAy has no copy of Tuesday
# Lake). Rule: matching_rule() with its defaults (every registry at every
# level, with the taxonomic decay of the pairs not inferred). Bayesian fit,
# 4 chains x 4000 draws, seed 2026 (the decay model keeps 4000 draws). Takes
# about 15 minutes. Packages used only here: pkgload, usethis.
#
# Run from the package root after data-raw/gateway.R:
#   Rscript data-raw/gray2015_gateway_calibration.R

pkgload::load_all(".", quiet = TRUE)
both <- combine_registries(gray2015, gateway)
sources <- gray2015_sources$study_id[gray2015_sources$validation]
groups <- list(c("13", "50"))
exclusions <- utils::read.csv(file.path("dev", "gateway", "gateway_exclusions.csv"),
                              colClasses = "character")
print(table(exclusions$web))
validation <- validate_rule(both, sources = sources, groups = groups, rule = matching_rule(),
                            exclusions = exclusions)
gray2015_gateway_calibration <- calibrate_links(validation, method = "bayes", n_chains = 4L,
                                                n_iter = 4000L, n_warmup = 2000L,
                                                quadrature = 20L, seed = 2026)
print(gray2015_gateway_calibration)
usethis::use_data(gray2015_gateway_calibration, overwrite = TRUE, compress = "xz")
