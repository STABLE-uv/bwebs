test_that("the size effect is 0 at the optimum ratio and negative elsewhere", {
  m <- size_model
  expect_s3_class(m, "bwebs_size_model")
  expect_lt(m$coefficients[["x2"]], 0)
  g <- size_effect(m, c(m$optimum, m$optimum - 1, m$optimum + 2, NA))
  expect_equal(g[1], 0)
  expect_true(all(g[2:3] < 0))
  expect_identical(g[4], 0)
  expect_equal(g[2], m$coefficients[["x2"]])
  expect_output(print(m), "optimum")
})

test_that("check_masses() accepts vectors and matrices of masses by node", {
  ids <- c("a", "b", "c")
  lm <- check_masses(c(a = 10, c = 0.1), ids, n_draws = 2L)
  expect_identical(dim(lm), c(3L, 2L))
  expect_equal(lm["a", ], c(1, 1))
  expect_true(all(is.na(lm["b", ])))
  mat <- matrix(c(1, 10, 100, 1000), 2, 2, dimnames = list(c("a", "b"), NULL))
  expect_equal(check_masses(mat, ids, n_draws = 2L)["b", ], c(1, 3))
  # a one-dimensional array, such as the result of tapply()
  arr <- tapply(c(10, 30, 0.1), c("a", "a", "c"), mean)
  expect_equal(check_masses(arr, ids)[, 1], c(a = log10(20), b = NA, c = -1))
  expect_error(check_masses(c(a = 1, z = 1), ids), "not in the metaweb")
  expect_error(check_masses(c(a = -1), ids), "positive")
  expect_error(check_masses(mat, ids, n_draws = 3L), "one per draw")
  expect_error(check_masses(c(1, 2), ids), "named")
})

# Metaweb of the tiny registry with the decay and a calibration with known
# classes and decay (as in test-decay.R).
size_setup <- function() {
  nodes <- prepare_nodes(tiny_taxa(), rule = tiny_rule(decay = TRUE))
  mw <- build_metaweb(nodes, tiny_registry())
  cal <- fake_calibration(c(exact = 1, genus = 0.5, family = 0, order = -0.5, class = -1),
                          rule = tiny_rule(decay = TRUE), decay = fake_decay())
  animals <- nodes$node_id[nodes$role == "consumer"]
  list(nodes = nodes, mw = mw, cal = cal, animals = animals)
}

test_that("the size term only lowers the probability of predation pairs", {
  st <- size_setup()
  # every node of the same mass: x = 0 for every pair, below the optimum
  masses <- stats::setNames(rep(1, nrow(st$nodes)), st$nodes$node_id)
  d0 <- sample_metawebs(st$mw, st$cal, n_draws = 500, seed = 3)
  d1 <- sample_metawebs(st$mw, st$cal, n_draws = 500, masses = masses, seed = 3)
  pred <- d0$links$resource %in% st$animals
  # the same random numbers: a lower probability can only remove links
  expect_true(all(d1$z <= d0$z))
  expect_identical(d1$z[!pred, ], d0$z[!pred, ])
  expect_true(any(d1$z[pred, ] < d0$z[pred, ]))
  expect_true(all(vapply(seq_len(500), function(k) {
    all(d1$background$present[[k]] %in% d0$background$present[[k]])
  }, logical(1))))
  bpred <- d0$background$pairs$resource %in% st$animals
  expect_equal(d1$background$pairs$p_mean[!bpred], d0$background$pairs$p_mean[!bpred])
  expect_true(all(d1$background$pairs$p_mean[bpred] < d0$background$pairs$p_mean[bpred]))
  expect_s3_class(d1$settings$size_model, "bwebs_size_model")
  # weights of 0 change nothing
  m0 <- size_model
  m0$weights[] <- 0
  d2 <- sample_metawebs(st$mw, st$cal, n_draws = 500, masses = masses, size_model = m0, seed = 3)
  expect_identical(d2$z, d0$z)
})

test_that("each draw uses its own masses and keeps them", {
  st <- size_setup()
  ids <- st$nodes$node_id
  # draw 1: every node of the same mass; draw 2: masses spread over four orders
  m <- matrix(1, length(ids), 2, dimnames = list(ids, NULL))
  m[, 2] <- 10^seq(-2, 2, length.out = length(ids))
  d <- sample_metawebs(st$mw, st$cal, n_draws = 2, masses = m, seed = 1)
  expect_identical(ncol(d$z), 2L)
  expect_error(sample_metawebs(st$mw, st$cal, n_draws = 3, masses = m), "one per draw")
  # the same draws as with the masses of draw 1 or of draw 2 in every draw
  d1 <- sample_metawebs(st$mw, st$cal, n_draws = 2, masses = m[, 1], seed = 1)
  d2 <- sample_metawebs(st$mw, st$cal, n_draws = 2, masses = m[, 2], seed = 1)
  expect_identical(d$z[, 1], d1$z[, 1])
  expect_identical(d$z[, 2], d2$z[, 2])
  expect_identical(d$background$present[[1]], d1$background$present[[1]])
  expect_identical(d$background$present[[2]], d2$background$present[[2]])
  # the draws keep the masses of each draw
  expect_equal(d$masses, m[st$mw$nodes$node_id, ])
  expect_output(print(d), "drawn with each draw")
  expect_output(print(d1), "the same in every draw")
  expect_null(sample_metawebs(st$mw, st$cal, n_draws = 2, seed = 1)$masses)
})

test_that("score_metaweb() applies the size term with the weights of the model", {
  st <- size_setup()
  e <- st$mw$edges[!st$mw$edges$cannibal, ]
  observed <- data.frame(resource = e$resource[1:2], consumer = e$consumer[1:2])
  masses <- stats::setNames(rep(1, nrow(st$nodes)), st$nodes$node_id)
  s0 <- score_metaweb(st$mw, st$cal, observed)
  s1 <- score_metaweb(st$mw, st$cal, observed, masses = masses)
  k0 <- paste(s0$pairs$resource, s0$pairs$consumer)
  k1 <- paste(s1$pairs$resource, s1$pairs$consumer)
  p1 <- s1$pairs$p[match(k0, k1)]
  pred <- s0$pairs$resource %in% st$animals & s0$pairs$p > 0 & s0$pairs$p < 1
  g <- size_effect(size_model, 0)
  w <- ifelse(s0$pairs$source == "rule", size_model$weights[["rule"]],
              size_model$weights[["decay"]])
  expect_equal(p1[pred], stats::plogis(stats::qlogis(s0$pairs$p[pred]) + w[pred] * g))
  expect_equal(p1[!pred], s0$pairs$p[!pred])
  expect_output(print(s1), "size term")
})

test_that("size_curve_fit() gives the estimates of glm() with one intercept per consumer", {
  set.seed(1)
  d <- do.call(rbind, lapply(seq_len(60), function(k) {
    x <- stats::rnorm(15, 2, 2)
    data.frame(cons = paste0("c", k), x = x,
               y = stats::rbinom(15, 1, stats::plogis(stats::rnorm(1) + 0.7 * x - 0.16 * x^2)))
  }))
  f <- size_curve_fit(d$y, d$x, d$cons)
  ok <- tapply(d$y, d$cons, function(v) any(v > 0) && any(v < 1))
  dd <- d[d$cons %in% names(ok)[ok], ]
  g <- stats::glm(y ~ 0 + factor(cons) + x + I(x^2), stats::binomial, data = dd)
  expect_true(f$converged)
  expect_equal(unname(f$coefficients), unname(stats::coef(g)[c("x", "I(x^2)")]), tolerance = 1e-6)
  expect_equal(unname(f$std_error), unname(sqrt(diag(stats::vcov(g)))[c("x", "I(x^2)")]),
               tolerance = 1e-5)
  expect_identical(f$pairs, nrow(dd))
})
