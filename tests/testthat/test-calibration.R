test_that("Gauss-Hermite nodes and weights integrate polynomials exactly", {
  gh <- gauss_hermite(5)
  # integral of exp(-x^2) = sqrt(pi); of x^2 exp(-x^2) = sqrt(pi) / 2;
  # of x^4 exp(-x^2) = 3 sqrt(pi) / 4; odd moments vanish
  expect_equal(sum(gh$w), sqrt(pi))
  expect_equal(sum(gh$w * gh$x^2), sqrt(pi) / 2)
  expect_equal(sum(gh$w * gh$x^4), 3 * sqrt(pi) / 4)
  expect_equal(sum(gh$w * gh$x^3), 0)
  expect_identical(gauss_hermite(1), list(x = 0, w = sqrt(pi)))
})

test_that("the marginal likelihood matches numerical integration", {
  # one web, two classes: integrate the conditional likelihood by brute force
  Y <- matrix(c(7, 2), 1, 2)
  M <- matrix(c(10, 9), 1, 2)
  theta <- c(a = 0.3, b = -0.8, log_sigma = log(0.9))
  f <- function(u) {
    p <- stats::plogis(theta[1:2] + u)
    exp(sum(Y * log(p) + (M - Y) * log1p(-p))) * stats::dnorm(u, 0, exp(theta[3]))
  }
  ref <- log(stats::integrate(Vectorize(f), -Inf, Inf, rel.tol = 1e-12)$value)
  expect_equal(marginal_loglik(theta, Y, M, gauss_hermite(20)), ref, tolerance = 1e-8)
})

test_that("fit_calibration() recovers simulated coefficients and agrees with lme4", {
  withr::local_seed(42)
  alpha <- c(exact = 1, genus = 0, family = -1)
  webs <- paste0("w", 1:12)
  u <- stats::rnorm(length(webs), 0, 0.7)
  counts <- do.call(rbind, lapply(seq_along(webs), function(w) {
    m <- c(80, 120, 200)
    y <- stats::rbinom(3, m, stats::plogis(alpha + u[w]))
    data.frame(study_id = webs[w], evidence_class = names(alpha), predicted = m, observed = y)
  }))
  fit <- fit_calibration(counts, quadrature = 20)
  expect_true(fit$converged)
  expect_true(fit$information_pd)
  est <- stats::setNames(fit$coefficients$estimate, fit$coefficients$evidence_class)
  se <- stats::setNames(fit$coefficients$std_error, fit$coefficients$evidence_class)
  expect_true(all(abs(est[names(alpha)] - alpha) < 4 * se[names(alpha)]))
  skip_if_not_installed("lme4")
  counts$evidence_class <- factor(counts$evidence_class, levels = fit$coefficients$evidence_class)
  g <- lme4::glmer(cbind(observed, predicted - observed) ~ 0 + evidence_class + (1 | study_id),
                   data = counts, family = stats::binomial, nAGQ = 20)
  expect_equal(unname(est), unname(lme4::fixef(g)), tolerance = 1e-4)
  expect_equal(fit$sigma, attr(lme4::VarCorr(g)$study_id, "stddev")[[1]], tolerance = 1e-3)
})

test_that("validate_web() compares inferred and recorded links of a held-out web", {
  # Web A: Baetis rhodani eats Cocconeis placentula; Perla bipunctata eats
  # Baetis rhodani. Web B (training): the same genera, plus Perla eating
  # Cocconeis. Inferred for A at genus level: 3 links, 2 recorded in A.
  rec <- two_web_records()
  reg <- as_registry(rec, registry = "gray2015")
  vn <- registry_nodes(reg, "A", tiny_rule())
  expect_identical(stats::setNames(vn$nodes$rank, vn$nodes$node_id),
                   c(`Cocconeis placentula` = "species", `Baetis rhodani` = "species",
                     `Perla bipunctata` = "species"))
  expect_identical(vn$nodes$role, c("producer", "consumer", "consumer"))
  res <- validate_web(reg, "A", tiny_rule())
  expect_identical(res$web$n_predicted, 3L)
  expect_identical(res$web$n_predicted_observed, 2L)
  expect_identical(res$counts$evidence_class, "genus")
  expect_identical(res$counts$observed, 2L)
  res_b <- validate_web(reg, "B", tiny_rule())
  expect_identical(res_b$web$precision, 1)
  expect_identical(res$web$left_out, "A")

  # Study C repeats the records of A (same system): left in the training
  # registry, it predicts the two links of A exactly
  rec_c <- rbind(rec, transform(rec[rec$study_id == "A", ], study_id = "C"))
  reg_c <- as_registry(rec_c, registry = "gray2015")
  alone <- validate_web(reg_c, "A", tiny_rule())
  expect_identical(alone$counts$evidence_class, c("exact", "genus"))
  expect_identical(alone$counts$predicted, c(2L, 1L))
  expect_identical(alone$counts$observed, c(2L, 0L))
  together <- validate_web(reg_c, "A", tiny_rule(), exclude = "C")
  expect_identical(together$counts, res$counts)
  expect_identical(together$web$left_out, "A, C")
  expect_identical(left_out_sets(list(c("C", "A")), c("A", "B"), reg_c), list(c("A", "C"), "B"))
  expect_identical(left_out_sets(NULL, c("A", "B"), reg_c), list("A", "B"))
})

test_that("validation keeps the same links when only records matched on both sides are kept", {
  full <- function(registry, nodes, rule, sides) rep(TRUE, nrow(registry))
  filtered <- validate_web(gray2015, "49", tiny_rule())
  local_mocked_bindings(relevant_rows = full)
  expect_identical(validate_web(gray2015, "49", tiny_rule()), filtered)
})

test_that("validate_rule() leaves out the studies in `exclusions` for their web only", {
  # Study C copies web A; excluded when A is rebuilt, A is predicted from B at
  # genus level only, and the rebuilding of B does not change
  rec <- two_web_records()
  rec_c <- rbind(rec, transform(rec[rec$study_id == "A", ], study_id = "C"))
  reg_c <- as_registry(rec_c, registry = "gray2015")
  ex <- data.frame(web = "A", study_id = "C")
  v0 <- validate_rule(reg_c, sources = c("A", "B"), rule = tiny_rule(), quiet = TRUE)
  v1 <- validate_rule(reg_c, sources = c("A", "B"), rule = tiny_rule(), exclusions = ex,
                      quiet = TRUE)
  expect_identical(v0$webs$n_excluded, c(0L, 0L))
  expect_identical(v1$webs$n_excluded, c(1L, 0L))
  expect_identical(v1$webs$left_out, c("A", "B"))
  a1 <- v1$counts[v1$counts$study_id == "A", ]
  expect_identical(a1$evidence_class, "genus")
  expect_identical(c(a1$predicted, a1$observed), c(3L, 2L))
  expect_identical(v1$webs[v1$webs$study_id == "B", ], v0$webs[v0$webs$study_id == "B", ])
  expect_identical(v1$settings$exclusions, ex)
  expect_error(validate_rule(reg_c, "A", rule = tiny_rule(), exclusions = data.frame(web = "A")),
               "not found")
  # only the records of the training registries rebuild the web
  other <- reg_c[reg_c$study_id == "C", ]
  other$registry <- "other"
  mixed <- rbind(reg_c[reg_c$study_id != "C", ], other)
  v2 <- validate_rule(mixed, sources = "A", rule = tiny_rule(), train_registries = "gray2015",
                      quiet = TRUE)
  expect_identical(v2$counts$evidence_class, "genus")
  expect_identical(v2$settings$registries, "gray2015")
  expect_error(validate_rule(mixed, "A", rule = tiny_rule(), train_registries = "nope"),
               "train_registries")
})

test_that("calibrate_links() validates its input", {
  expect_error(calibrate_links(tiny_registry()), "sources")
  expect_error(calibrate_links(tiny_registry(), sources = "zz"), "not in `registry`")
  expect_error(calibrate_links(tiny_registry(), sources = c("s1", "s1")), "repeated")
  expect_error(calibrate_links(tiny_registry(), sources = "s1", groups = c("s1", "s2")),
               "list of character vectors")
  expect_error(calibrate_links(tiny_registry(), sources = "s1", groups = list("s1")),
               "at least two studies")
  expect_error(calibrate_links(tiny_registry(), sources = "s1", groups = list(c("s1", "zz"))),
               "not in `registry`")
  expect_error(calibrate_links(tiny_registry(), sources = "s1",
                               groups = list(c("s1", "s2"), c("s2", "s1"))),
               "more than once")
})

test_that("the bundled calibration has the documented structure", {
  cal <- gray2015_calibration
  expect_s3_class(cal, "bwebs_calibration")
  expect_identical(nrow(cal$webs), 13L)
  expect_true(cal$diagnostics$converged)
  expect_true(cal$diagnostics$information_pd)
  expect_true(all(c("exact", "genus", "family", "order") %in% cal$coefficients$evidence_class))
  expect_identical(sum(cal$counts$predicted), sum(cal$webs$n_predicted))
  expect_identical(cal$settings$groups, list(c("13", "50")))
  expect_identical(cal$webs$left_out[cal$webs$study_id %in% c("13", "50", "46")],
                   c("13, 50", "46", "50, 13"))
  expect_snapshot(print(cal))
})

test_that("the bundled calibration with GloBI has the documented structure", {
  cal <- gray2015_globi_calibration
  expect_s3_class(cal, "bwebs_calibration")
  expect_identical(cal$settings$registries, c("globi", "gray2015"))
  expect_true(same_rule(cal$settings$rule, matching_rule()))
  expect_true(same_rule(gray2015_calibration$settings$rule, matching_rule()))
  expect_identical(nrow(cal$webs), 13L)
  expect_true(cal$diagnostics$converged)
  expect_identical(cal$settings$groups, list(c("13", "50")))
  ex <- cal$settings$exclusions
  expect_true(all(c("web", "study_id") %in% names(ex)))
  expect_true(all(ex$web %in% cal$webs$study_id))
  expect_identical(cal$webs$n_excluded[cal$webs$study_id == "12"] > 0L, TRUE)
  expect_identical(sum(cal$counts$predicted), sum(cal$webs$n_predicted))
})

test_that("the bundled comparison of sources has the documented structure", {
  sc <- source_comparison
  expect_true(all(c("source", "primary_registry", "supplementary_max_level", "rule", "brier",
                    "brier_diff", "diff_vs_gray", "skill", "tss", "precision", "recall",
                    "tl_brier", "tl_skill", "tl_tss") %in% names(sc)))
  expect_true(all(sc$source %in% c("Gray", "GloBI", "Gray + GloBI")))
  expect_identical(sc$brier, sort(sc$brier))
  # the defaults of matching_rule() are the best combination
  best <- sc[1, ]
  expect_identical(best$rule, rule_label(matching_rule()))
  expect_identical(best$source, "Gray + GloBI")
  expect_true(is.na(best$primary_registry))
  expect_identical(sc$diff_vs_gray[sc$source == "Gray"], rep(0, sum(sc$source == "Gray")))
})

test_that("calibrate_links() runs end to end on small validation webs", {
  skip_on_cran()
  v <- validate_rule(gray2015, sources = c("46", "49", "58"), quiet = TRUE)
  expect_s3_class(v, "bwebs_validation")
  expect_identical(v$webs$study_id, c("46", "49", "58"))
  cal <- calibrate_links(gray2015, sources = c("46", "49", "58"), method = "ml", quiet = TRUE)
  expect_s3_class(cal, "bwebs_calibration")
  expect_identical(cal$counts, v$counts)
  expect_identical(cal$webs$n_predicted, v$webs$n_predicted)
  expect_null(cal$draws)
  expect_error(calibrate_links(v, sources = "46"), "taken from the validation")
})

test_that("the Bayesian fit samples the posterior of a small problem", {
  # One class, three webs: the posterior of (alpha, log sigma) on a grid.
  counts <- data.frame(study_id = c("a", "b", "c"), evidence_class = "genus",
                       predicted = c(10, 15, 20), observed = c(3, 7, 12), stringsAsFactors = FALSE)
  cm <- count_matrices(counts)
  gh <- gauss_hermite(20)
  lp <- function(a, ls) {
    marginal_loglik(c(a, ls), cm$Y, cm$M, gh) + stats::dnorm(a, 0, 1.5, log = TRUE) +
      stats::dnorm(exp(ls), 0, 2.5, log = TRUE) + ls
  }
  ga <- seq(-3, 3, length.out = 121)
  gs <- seq(-6, 2.5, length.out = 171)
  L <- outer(ga, gs, Vectorize(lp))
  W <- exp(L - max(L))
  W <- W / sum(W)
  mean_a <- sum(W * ga)
  mean_s <- sum(W * outer(ga * 0, exp(gs), "+"))
  fit <- fit_calibration_bayes(counts, prior_sd = 1.5, prior_sigma = 2.5, n_chains = 4L,
                               n_iter = 3000L, n_warmup = 1000L, quadrature = 20L, seed = 11)
  expect_lt(max(fit$rhat), 1.01)
  mcse_a <- stats::sd(fit$draws[, 1]) / sqrt(fit$ess[1])
  mcse_s <- stats::sd(exp(fit$draws[, 2])) / sqrt(fit$ess[2])
  expect_lt(abs(mean(fit$draws[, 1]) - mean_a), 4 * mcse_a)
  expect_lt(abs(mean(exp(fit$draws[, 2])) - mean_s), 4 * mcse_s)
})

test_that("R-hat and effective sample size detect poor mixing", {
  withr::local_seed(5)
  iid <- matrix(stats::rnorm(4000), 1000, 4)
  expect_lt(abs(rhat_split(iid) - 1), 0.01)
  expect_gt(ess_bulk(iid), 3000)
  shifted <- iid + matrix(c(0, 0, 0, 3), 1000, 4, byrow = TRUE)
  expect_gt(rhat_split(shifted), 1.1)
  ar <- apply(iid, 2, function(z) as.numeric(stats::filter(z, 0.9, method = "recursive")))
  expect_lt(ess_bulk(ar), 600)
})

test_that("the probability in a new system integrates over the web effect", {
  f <- function(u) stats::plogis(0.4 + u) * stats::dnorm(u, 0, 1.3)
  expect_equal(p_new_system(0.4, 1.3), stats::integrate(f, -Inf, Inf)$value, tolerance = 1e-8)
  expect_equal(p_new_system(c(0, 0), c(0.5, 2)), c(0.5, 0.5))
})
