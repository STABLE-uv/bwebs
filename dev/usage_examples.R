# Complete workflow of a user of bwebs (run in a clean session with
# Rscript --vanilla in an empty folder; it is the basis of the README, the
# vignette and tests/testthat/test-workflow.R).
library(bwebs)

# 1. Nodes of the community and their matching ladder
nodes <- prepare_nodes(example_taxa)
nodes

# 2. Deterministic metaweb with the evidence of each link (the bundled
#    registries of Gray et al. 2015 and of the freshwater webs of GATEWAy)
registry <- combine_registries(gray2015, gateway)
mw <- build_metaweb(nodes, registry)
mw
table(mw$edges$evidence_class)

# 3. Bayesian calibration of those registries (bundled; computed again with
#    data-raw/gray2015_gateway_calibration.R)
gray2015_gateway_calibration

# 3b. How the default rule was chosen (cross-validation by system and test on
#     Tuesday Lake)
head(gray2015_rule_comparison[, c("rule", "brier", "brier_diff", "brier_diff_se", "tl_brier")])

# 3c. Probability of every resource -> consumer pair (inferred links by their
#     class of evidence, the other pairs by the taxonomic decay)
P <- link_probabilities(mw, gray2015_gateway_calibration)
round(P[c("cocconeis_placentula", "cpom", "baetis_rhodani", "simulium"),
        c("gammarus_pulex", "cordulegaster_boltonii")], 2)

# 4. Plausible metawebs drawn from those probabilities
draws <- sample_metawebs(mw, gray2015_gateway_calibration, n_draws = 50, seed = 2026)
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

# 7b. Energy fluxes of every sample in every realization: the size term and
#     the fluxes with the masses and biomasses of the populations (in an
#     analysis, the draws of bmassa: realization d uses draw d of both);
#     engine = "fluxweb", the default, needs the package fluxweb
masses <- node_masses(mw, example_populations)
draws_m <- sample_metawebs(mw, gray2015_gateway_calibration, n_draws = 50, masses = masses,
                           seed = 2026)
fl <- calc_fluxes(draws_m, example_populations, temperature = 15, engine = "internal")
fl
summarize_indicators(fl, indicators = c("total_flux", "herbivory", "detritivory", "carnivory"),
                     samples = example_samples, by = "treatment")
contrast_indicators(fl, example_samples, by = "treatment", reference = "control",
                    indicators = "total_flux")

# 8. Export
out <- tempfile(fileext = ".csv")
utils::write.csv(by_treatment, out, row.names = FALSE)
file.exists(out)
