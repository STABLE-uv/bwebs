test_that("link_probabilities() gives each link the probability of its class", {
  nodes <- prepare_nodes(tiny_taxa(), rule = tiny_rule())
  ml <- data.frame(resource = c("cocc", "cpom"), consumer = c("asel", "asel"),
                   action = c("add", "remove"))
  mw <- build_metaweb(nodes, tiny_registry(), manual_links = ml)
  alpha <- c(genus = 1, family = 0, order = -1)
  cal <- fake_calibration(alpha, sigma = 1.2)
  P <- link_probabilities(mw, cal)
  expect_identical(dim(P), c(nrow(nodes), nrow(nodes)))
  expect_identical(dimnames(P), list(resource = nodes$node_id, consumer = nodes$node_id))
  e <- mw$edges[mw$edges$evidence_class %in% names(alpha), ]
  expect_equal(P[cbind(e$resource, e$consumer)],
               p_new_system(unname(alpha[e$evidence_class]), 1.2))
  expect_identical(P["cocc", "asel"], 1)      # manual link
  expect_identical(P["cpom", "asel"], 0)      # removed
  expect_identical(sum(P > 0), nrow(mw$edges) - sum(mw$edges$cannibal))
  expect_true(all(diag(P) == 0))
  # other system effects
  Pa <- link_probabilities(mw, cal, system = "average")
  expect_equal(Pa[cbind(e$resource, e$consumer)], stats::plogis(unname(alpha[e$evidence_class])))
  P1 <- link_probabilities(mw, cal, system = 1)
  expect_equal(P1[cbind(e$resource, e$consumer)], stats::plogis(unname(alpha[e$evidence_class]) + 1))
  expect_error(link_probabilities(mw, cal, system = "new"), "must be one of")
})

test_that("link_probabilities() agrees with the inclusion of the sampled webs", {
  st <- decay_setup()
  P <- link_probabilities(st$mw, st$cal)
  d <- st$mw$decay
  expect_equal(P[cbind(d$resource, d$consumer)],
               decay_p_new_of(st$cal$decay, d, n_nodes = nrow(st$nodes)))
  dr <- sample_metawebs(st$mw, st$cal, n_draws = 4000, parameter_uncertainty = FALSE, seed = 4)
  l <- dr$links
  expect_lt(max(abs(l$inclusion - P[cbind(l$resource, l$consumer)])), 0.04)
  bg <- dr$background$pairs
  expect_lt(max(abs(bg$inclusion - P[cbind(bg$resource, bg$consumer)])), 0.04)
  # without the decay, the pairs not inferred have probability 0
  mw0 <- st$mw
  mw0$decay <- NULL
  P0 <- link_probabilities(mw0, st$cal)
  expect_true(all(P0[cbind(d$resource, d$consumer)] == 0))
  cal0 <- st$cal
  cal0$decay <- NULL
  expect_error(link_probabilities(st$mw, cal0), "no decay model")
})

test_that("body masses lower the probability of predation pairs only", {
  st <- decay_setup()
  ids <- st$nodes$node_id
  masses <- stats::setNames(rep(1, length(ids)), ids)   # equal masses: far from the optimum
  P <- link_probabilities(st$mw, st$cal)
  Pm <- link_probabilities(st$mw, st$cal, masses = masses)
  animal <- ids[st$nodes$role == "consumer"]
  pred <- outer(ids %in% animal, rep(TRUE, length(ids))) & P > 0 & P < 1
  expect_true(all(Pm[pred] < P[pred]))
  expect_equal(Pm[!pred], P[!pred])
  expect_error(link_probabilities(st$mw, st$cal, masses = c(x = 1)), "not in the metaweb")
})

test_that("uncalibrated links need a decision, as when sampling", {
  nodes <- prepare_nodes(tiny_taxa(), rule = tiny_rule())
  mw <- build_metaweb(nodes, tiny_registry())
  cal <- fake_calibration(c(genus = 0, family = 0))
  expect_error(link_probabilities(mw, cal), "no calibrated probability")
  Pe <- link_probabilities(mw, cal, uncalibrated = "exclude")
  Pk <- link_probabilities(mw, cal, uncalibrated = "keep")
  e <- mw$edges[mw$edges$evidence_class == "order", ]
  expect_true(all(Pe[cbind(e$resource, e$consumer)] == 0))
  expect_true(all(Pk[cbind(e$resource, e$consumer)] == 1))
})
