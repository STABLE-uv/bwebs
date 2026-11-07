test_that("links are kept with the probability of their class", {
  mw <- build_metaweb(prepare_nodes(tiny_taxa(), rule = tiny_rule()), tiny_registry())
  cal <- fake_calibration(c(genus = 2, family = 0, order = -1.5))
  d <- sample_metawebs(mw, cal, n_draws = 4000, system = "average",
                       parameter_uncertainty = FALSE, seed = 1)
  expect_s3_class(d, "bwebs_draws")
  p <- stats::plogis(c(genus = 2, family = 0, order = -1.5))[d$links$class_used]
  mcse <- sqrt(p * (1 - p) / 4000)
  expect_true(all(abs(d$links$inclusion - p) < 4 * mcse))
  expect_equal(d$links$p_mean, unname(p))
})

test_that("sampling is reproducible and leaves the global stream alone", {
  mw <- build_metaweb(prepare_nodes(tiny_taxa(), rule = tiny_rule()), tiny_registry())
  cal <- fake_calibration(c(genus = 0.5, family = 0, order = -1))
  a <- sample_metawebs(mw, cal, n_draws = 5, seed = 3)
  b <- sample_metawebs(mw, cal, n_draws = 5, seed = 3)
  expect_identical(a$z, b$z)
  set.seed(10); x <- stats::runif(1)
  set.seed(10); sample_metawebs(mw, cal, n_draws = 5, seed = 3); y <- stats::runif(1)
  expect_identical(x, y)
})

test_that("a random system effect makes draws vary together", {
  mw <- build_metaweb(prepare_nodes(tiny_taxa(), rule = tiny_rule()), tiny_registry())
  cal <- fake_calibration(c(genus = 0, family = 0, order = 0), sigma = 2)
  d <- sample_metawebs(mw, cal, n_draws = 300, system = "random", seed = 2)
  expect_gt(stats::sd(d$system_effect), 1)
  fixed <- sample_metawebs(mw, cal, n_draws = 10, system = 1.5, parameter_uncertainty = FALSE, seed = 2)
  expect_true(all(fixed$system_effect == 1.5))
  expect_equal(unname(fixed$p[, 1]), rep(stats::plogis(1.5), 3))
})

test_that("manual links are always present and uncalibrated links need a decision", {
  nodes <- prepare_nodes(tiny_taxa(), rule = tiny_rule())
  ml <- data.frame(resource = "cocc", consumer = "asel", action = "add")
  mw <- build_metaweb(nodes, tiny_registry(), manual_links = ml)
  cal <- fake_calibration(c(genus = 0, family = 0))
  expect_error(sample_metawebs(mw, cal, n_draws = 3, seed = 1), "no calibrated probability")
  d <- sample_metawebs(mw, cal, n_draws = 3, uncalibrated = "exclude", seed = 1)
  expect_true(all(d$z[d$links$status == "manual", ]))
  expect_false(any(d$z[d$links$status == "uncalibrated_excluded", ]))
  k <- sample_metawebs(mw, cal, n_draws = 3, uncalibrated = "keep", seed = 1)
  expect_true(all(k$z[k$links$status == "uncalibrated_kept", ]))
})

test_that("links from registries outside the calibration are not calibrated", {
  mw <- build_metaweb(prepare_nodes(tiny_taxa(), rule = tiny_rule()), tiny_registry())
  cal <- fake_calibration(c(genus = 0, family = 0, order = 0), registries = "other")
  expect_warning(expect_error(sample_metawebs(mw, cal, n_draws = 2), "registry not calibrated"),
                 class = "bwebs_warning_registries")
  # a calibration with more registries than the metaweb covers its links, but
  # its probabilities were estimated for another registry
  both <- fake_calibration(c(genus = 0, family = 0, order = 0), registries = c("gray2015", "globi"))
  expect_warning(d <- sample_metawebs(mw, both, n_draws = 2, seed = 1), "registries")
  expect_true(all(d$links$status == "calibrated"))
  expect_no_warning(sample_metawebs(mw, fake_calibration(c(genus = 0, family = 0, order = 0)),
                                    n_draws = 2, seed = 1))
})

test_that("a calibration of another rule is rejected", {
  mw <- build_metaweb(prepare_nodes(tiny_taxa(), rule = tiny_rule()), tiny_registry())
  cal <- fake_calibration(c(genus = 0, family = 0, order = 0), rule = matching_rule(mode = "fixed"))
  expect_error(sample_metawebs(mw, cal, n_draws = 2), "different matching rule")
  expect_error(sample_metawebs(mw, cal, n_draws = 0), class = "bwebs_error_input")
})

test_that("parameters are taken from the posterior draws when the calibration has them", {
  mw <- build_metaweb(prepare_nodes(tiny_taxa(), rule = tiny_rule()), tiny_registry())
  cal <- fake_calibration(c(genus = 0, family = 0, order = 0))
  cal$draws <- cbind(genus = rep(c(-3, 3), 500), family = 0, order = 0, log_sigma = log(0.01))
  d <- sample_metawebs(mw, cal, n_draws = 2000, system = "average", seed = 4)
  p_genus <- round(d$p["genus", ], 8)
  expect_true(all(p_genus %in% round(stats::plogis(c(-3, 3)), 8)))
  expect_equal(mean(p_genus > 0.5), 0.5, tolerance = 0.05)
  fixed <- sample_metawebs(mw, cal, n_draws = 10, system = "average",
                           parameter_uncertainty = FALSE, seed = 4)
  expect_equal(unname(fixed$p["genus", 1]), stats::plogis(0))
})
