test_that("the documented workflow runs from nodes to contrasts", {
  nodes <- prepare_nodes(example_taxa)
  mw <- build_metaweb(nodes, gray2015)
  expect_s3_class(mw, "bwebs_metaweb")
  expect_true(all(mw$edges$evidence_class %in% c(evidence_classes, "manual")))
  # the shrimp takes the records of Asellus
  expect_true(any(mw$edges$consumer == "dugastella_valentina"))
  draws <- sample_metawebs(mw, gray2015_calibration, n_draws = 5, seed = 2026)
  expect_identical(dim(draws$z), c(nrow(draws$links), 5L))
  ind <- calc_indicators(draws, example_community, permanent = c("cpom", "fpom"))
  expect_identical(nrow(ind$values), 5L * 6L)
  s <- summarize_indicators(ind, indicators = c("S", "C"), samples = example_samples, by = "treatment")
  expect_identical(s$treatment, c("control", "control", "exposed", "exposed"))
  k <- contrast_indicators(ind, example_samples, by = "treatment", reference = "control",
                           indicators = "S")
  # the exposed samples lose taxa in every draw
  expect_lt(k$upper, 0)
  det <- calc_indicators(mw, example_community, permanent = c("cpom", "fpom"))
  expect_identical(unique(det$values$draw), 0L)
})

test_that("the workflow runs with the registries and calibration of Gray + GATEWAy", {
  nodes <- prepare_nodes(example_taxa)
  mw <- build_metaweb(nodes, combine_registries(gray2015, gateway))
  draws <- sample_metawebs(mw, gray2015_gateway_calibration, n_draws = 5, seed = 2026)
  expect_identical(ncol(draws$z), 5L)
  expect_false(any(is.na(draws$links$p_mean[draws$links$status == "calibrated"])))
  ind <- calc_indicators(draws, example_community, permanent = c("cpom", "fpom"))
  expect_identical(nrow(ind$values), 5L * 6L)
  # the matrix the draws come from
  P <- link_probabilities(mw, gray2015_gateway_calibration)
  expect_identical(dim(P), rep(nrow(nodes), 2L))
  expect_true(all(P >= 0 & P <= 1))
  l <- draws$links[draws$links$status == "calibrated", ]
  expect_true(all(P[cbind(l$resource, l$consumer)] > 0))
})

test_that("the workflow runs from the populations to the contrasts of the fluxes", {
  nodes <- prepare_nodes(example_taxa)
  mw <- build_metaweb(nodes, combine_registries(gray2015, gateway))
  masses <- node_masses(mw, example_populations)
  expect_identical(sum(!is.na(masses)), 17L)
  draws <- sample_metawebs(mw, gray2015_gateway_calibration, n_draws = 10, masses = masses,
                           seed = 2026)
  fl <- calc_fluxes(draws, example_populations, temperature = 15, engine = "internal")
  expect_identical(nrow(fl$values), 10L * 6L)
  # draw 10 keeps only 5 inferred links: its consumers eat mostly one another
  # through the links of the decay, and the web has no steady state
  expect_true(all(fl$values$status %in% c("ok", "unsolved")))
  ok <- fl$values[fl$values$status == "ok", ]
  expect_gt(nrow(ok), 50L)
  expect_equal(ok$total_flux, ok$herbivory + ok$detritivory + ok$carnivory)
  k <- contrast_indicators(fl, example_samples, by = "treatment", reference = "control",
                           indicators = c("total_flux", "metabolism"), na_rm = TRUE)
  # the exposed samples lose biomass of consumers: less metabolism in every draw
  expect_lt(k$upper[k$indicator == "metabolism"], 0)
})
