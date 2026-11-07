test_that("the Brier score of a web matches a hand calculation", {
  test <- data.frame(evidence_class = c("genus", "family", "order"), predicted = c(4, 6, 2),
                     observed = c(3, 1, 1), stringsAsFactors = FALSE)
  p <- c(genus = 0.7, family = 0.2)   # no probability for "order": p = 0
  s <- web_brier(test, p, n_pairs = 50, n_recorded = 7, n_predicted = 12, n_predicted_observed = 5,
                 p_ref = 0.1)
  # inferred: 3 * 0.3^2 + 1 * 0.7^2 + 1 * 0.8^2 + 5 * 0.2^2 + 1 (order, observed, p = 0)
  # missed: 7 - 5 = 2
  brier <- (3 * 0.09 + 0.49 + 0.64 + 5 * 0.04 + 1 + 2) / 50
  ref <- (7 * 0.9^2 + 43 * 0.1^2) / 50
  expect_equal(s$brier, brier)
  expect_equal(s$brier_reference, ref)
  expect_equal(s$skill, 1 - brier / ref)
  expect_equal(s$tss, 5 / 7 - 7 / 43)
  expect_equal(s$precision, 5 / 12)
  expect_equal(s$recall, 5 / 7)
  expect_identical(s$n_uncalibrated, 2L)
})

test_that("score_rule() refits without the system of each web", {
  # Four webs, two classes; webs a and b are the same system.
  counts <- data.frame(
    study_id = rep(c("a", "b", "c", "d"), each = 2), evidence_class = rep(c("genus", "family"), 4),
    predicted = c(20, 30, 25, 35, 18, 40, 22, 28), observed = c(14, 9, 15, 12, 9, 10, 15, 6),
    stringsAsFactors = FALSE
  )
  tot <- stats::aggregate(cbind(predicted, observed) ~ study_id, counts, sum)
  webs <- data.frame(study_id = tot$study_id, left_out = c("a, b", "b, a", "c", "d"),
                     n_pairs = 400, n_recorded = tot$observed + 5, n_predicted = tot$predicted,
                     n_predicted_observed = tot$observed, stringsAsFactors = FALSE)
  v <- structure(list(counts = counts, webs = webs,
                      settings = list(rule = matching_rule(), sources = webs$study_id)),
                 class = "bwebs_validation")
  s <- score_rule(v, method = "ml")
  expect_s3_class(s, "bwebs_score")
  fit_a <- fit_calibration(counts[counts$study_id %in% c("c", "d"), ], 20)
  p_a <- stats::setNames(fit_a$coefficients$p_new, fit_a$coefficients$evidence_class)
  p_ref_a <- mean(webs$n_recorded[3:4] / webs$n_pairs[3:4])
  ref_a <- web_brier(counts[counts$study_id == "a", ], p_a, 400, webs$n_recorded[1],
                     webs$n_predicted[1], webs$n_predicted_observed[1], p_ref_a)
  expect_equal(s$webs$brier[1], ref_a$brier)
  fit_c <- fit_calibration(counts[counts$study_id != "c", ], 20)
  p_c <- stats::setNames(fit_c$coefficients$p_new, fit_c$coefficients$evidence_class)
  p_ref_c <- mean(webs$n_recorded[-3] / webs$n_pairs[-3])
  ref_c <- web_brier(counts[counts$study_id == "c", ], p_c, 400, webs$n_recorded[3],
                     webs$n_predicted[3], webs$n_predicted_observed[3], p_ref_c)
  expect_equal(s$webs$brier[3], ref_c$brier)
  expect_equal(s$webs$skill[1], ref_a$skill)
  expect_identical(sort(unique(s$probabilities$study_id)), c("a", "b", "c", "d"))
  expect_equal(s$summary$brier, mean(s$webs$brier))
  expect_error(score_rule(list()), "validate_rule")
})

test_that("candidate rules form the documented grid", {
  r <- candidate_rules()
  expect_length(r, 21L)
  expect_true(all(vapply(r, inherits, logical(1), "bwebs_rule")))
  expect_identical(names(r), unname(vapply(r, rule_label, character(1))))
  expect_true(all(c("fixed/exact", "climb/family/class/1", "ladder/class/5",
                    "climb/family/class/1 no-downward") %in% names(r)))
  expect_identical(candidate_rules(manual_diet = "x")[[1]]$manual_diet, "x")
})

test_that("score_metaweb() scores the probabilities of a metaweb against observed links", {
  nodes <- prepare_nodes(tiny_taxa(), rule = matching_rule(mode = "fixed", base_level = "family", climb_max = "family"))
  mw <- build_metaweb(nodes, tiny_registry())
  cal <- fake_calibration(c(genus = 1, family = -1, order = 0), sigma = 1,
                          rule = matching_rule(mode = "fixed", base_level = "family", climb_max = "family"))
  e <- mw$edges[!mw$edges$cannibal, ]
  observed <- data.frame(resource = c(e$resource[1], "cocc"), consumer = c(e$consumer[1], "gerr"))
  s <- score_metaweb(mw, cal, observed)
  expect_s3_class(s, "bwebs_web_score")
  consumers <- unique(observed$consumer)
  n_pairs <- nrow(nodes) * length(consumers) - length(consumers)
  ev <- e[e$consumer %in% consumers, ]
  p <- cal$coefficients$p_new[match(ev$evidence_class, cal$coefficients$evidence_class)]
  y <- paste(ev$resource, ev$consumer) %in% paste(observed$resource, observed$consumer)
  missed <- sum(!paste(observed$resource, observed$consumer) %in% paste(ev$resource, ev$consumer))
  expect_equal(s$summary$brier, (sum((p - y)^2) + missed) / n_pairs)
  expect_identical(s$summary$n_pairs, as.integer(n_pairs))
  expect_error(score_metaweb(mw, cal, data.frame(resource = "zz", consumer = "asel")), "not in the metaweb")
  other <- fake_calibration(c(genus = 1), rule = matching_rule(mode = "ladder"))
  expect_error(score_metaweb(mw, other, observed), "different matching rule")
})

test_that("compare_rules() ranks rules by their cross-validated Brier score", {
  skip_on_cran()
  rules <- candidate_rules()[c("fixed/genus", "climb/family/class/1")]
  cmp <- compare_rules(gray2015, sources = c("46", "49", "58", "59"), rules = rules,
                       method = "ml", quiet = TRUE)
  expect_s3_class(cmp, "bwebs_comparison")
  expect_identical(sort(cmp$table$rule), sort(names(rules)))
  expect_equal(min(cmp$table$brier_diff), 0)
  expect_true(all(diff(cmp$table$brier) >= 0))
  expect_s3_class(cmp$validations[[1]], "bwebs_validation")
  expect_error(compare_rules(gray2015, sources = "46", rules = list(matching_rule())), "named list")
})
