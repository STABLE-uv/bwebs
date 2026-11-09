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
