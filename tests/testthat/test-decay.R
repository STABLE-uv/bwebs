test_that("evidence profiles give the closest relation at each consumer level", {
  nodes <- prepare_nodes(tiny_taxa(), rule = tiny_rule())
  pr <- decay_profiles(nodes, tiny_registry(), consumers = c("baet", "simu", "perl", "gerr"))
  lev <- function(r, c) decay_levels[pr$level[r, c, ]]
  # Baetis muticus: Baetis eats Cocconeis (same class as Navicula); the
  # Ephemeroptera also eat Navicula gregaria (same genus)
  expect_identical(lev("navi", "baet"),
                   c(NA, "class", "class", "class", "genus", "genus", "genus", "genus"))
  # CPOM: only by name; the Insecta (Simulium) eat it
  expect_identical(lev("cpom", "baet"), c(NA, rep("any", 4), rep("exact", 3)))
  expect_identical(lev("cpom", "simu"), c(NA, rep("exact", 7)))
  # a coarse record (Gerris eats "Diptera") relates at the order
  expect_identical(lev("simu", "gerr"), c(NA, rep("order", 7)))
  expect_identical(lev("baet", "perl"), c(NA, rep("genus", 7)))
  expect_false(any(pr$generic))
  # distinct resources known at each level for Baetis
  expect_identical(unname(pr$known["baet", ]), c(0L, 1L, 1L, 1L, 2L, 5L, 5L, NA))
})

test_that("records of generic groups and life stages are taken into account", {
  reg <- tiny_registry()
  extra <- reg[reg$consumer == "Ecdyonurus venosus", ]
  extra$record_id <- "x1"
  extra$resource <- extra$resource_std <- "Diatom"
  extra[c("res_genus", "res_family", "res_order", "res_class")] <- NA_character_
  reg <- rbind(reg, extra)
  nodes <- prepare_nodes(tiny_taxa(), rule = tiny_rule())
  pr <- decay_profiles(nodes, reg, consumers = c("baet", "gerr"))
  # the Ephemeroptera eat "Diatom", which applies to the diatoms (class)
  expect_identical(unname(pr$generic["navi", "baet", ]), rep(c(FALSE, TRUE), c(4L, 4L)))
  expect_false(any(pr$generic["cpom", "baet", ]))
  # a record of adults does not count for a larva
  reg$con_stage[reg$consumer == "Baetis rhodani"] <- "adult"
  taxa <- tiny_taxa()
  taxa$life_stage <- ifelse(taxa$node_id == "baet", "larva", NA)
  pr2 <- decay_profiles(prepare_nodes(taxa, rule = tiny_rule()), reg, consumers = "baet")
  expect_true(is.na(pr2$level["cocc", "baet", "genus"]))
  expect_identical(decay_levels[pr2$level["cocc", "baet", "order"]], "class")
})

# Two consumers more: Caenis (a family without records, in the order of
# Baetis and Ecdyonurus) and Hydra (a phylum without records).
decay_taxa <- function() {
  extra <- data.frame(
    node_id = c("caen", "hydr"), taxon = c("Caenis", "Hydra"), rank = "genus", role = "consumer",
    node_type = "taxon", genus = c("Caenis", "Hydra"), subfamily = NA_character_,
    family = c("Caenidae", "Hydridae"), order = c("Ephemeroptera", "Anthoathecata"),
    class = c("Insecta", "Hydrozoa"), stringsAsFactors = FALSE
  )
  rbind(tiny_taxa(), extra)
}

test_that("decay_evidence() resolves the consumer at its closest level with records", {
  nodes <- prepare_nodes(decay_taxa(), rule = tiny_rule())
  pr <- decay_profiles(nodes, tiny_registry(), consumers = c("baet", "caen", "hydr"))
  ev <- decay_evidence(pr, c("navi", "cpom", "navi", "cpom", "navi", "perl"),
                       c("baet", "baet", "caen", "caen", "hydr", "hydr"))
  # Baetis: its family eats Cocconeis (same class as Navicula), not CPOM
  # Caenis: the Ephemeroptera eat Navicula gregaria (same genus)
  # Hydra: no cnidarian has records; Navicula gregaria is eaten (genus), and
  # insects (Baetis, Diptera), of the class of Perla
  expect_identical(ev$consumer_level, c("family", "family", "order", "order", "any", "any"))
  expect_identical(ev$resource_level, c("class", "any", "genus", "any", "genus", "class"))
  expect_false(any(ev$generic))
  expect_identical(nrow(decay_evidence(pr, character(0), character(0))), 0L)
})

test_that("decay counts tabulate pairs and recorded links by evidence profile", {
  ev <- data.frame(consumer_level = c("family", "family", "order"),
                   resource_level = c("any", "any", "class"), generic = c(FALSE, FALSE, TRUE))
  tab <- decay_count_table("w", ev, c(TRUE, FALSE, TRUE))
  expect_identical(tab$pairs, c(2L, 1L))
  expect_identical(tab$observed, c(1L, 1L))
  expect_identical(tab$resource_level, c("any", "class"))
  all <- decay_all_profiles()
  expect_identical(nrow(all), 80L)
  expect_identical(decay_profile_id(all$consumer_level, all$resource_level, all$generic), 1:80)
})

test_that("validate_rule() with the decay accounts for every pair not inferred", {
  v <- validate_rule(tiny_registry(), c("s1", "s2", "s3"), rule = tiny_rule(decay = TRUE),
                     quiet = TRUE)
  d <- v$decay
  expect_true(all(c("study_id", "consumer_level", "resource_level", "generic", "pairs",
                    "observed") %in% names(d)))
  w <- v$webs
  pairs <- vapply(w$study_id, function(s) sum(d$pairs[d$study_id == s]), numeric(1))
  observed <- vapply(w$study_id, function(s) sum(d$observed[d$study_id == s]), numeric(1))
  expect_equal(unname(pairs), w$n_pairs - w$n_predicted)
  expect_equal(unname(observed), w$n_recorded - w$n_predicted_observed)
  expect_null(validate_rule(tiny_registry(), "s1", rule = tiny_rule(), quiet = TRUE)$decay)
})

test_that("the likelihood of the decay model has the right gradient", {
  counts <- data.frame(study_id = rep(c("a", "b", "c"), each = 4),
                       consumer_level = rep(c("family", "order", "phylum", "family"), 3),
                       resource_level = rep(c("any", "class", "genus", "order"), 3),
                       generic = rep(c(FALSE, FALSE, FALSE, TRUE), 3),
                       pairs = c(50, 20, 10, 5, 40, 30, 5, 6, 60, 10, 8, 4),
                       observed = c(10, 6, 1, 3, 4, 10, 2, 2, 15, 1, 3, 3))
  gh <- gauss_hermite(10)
  check <- function(counts, th) {
    webs <- decay_webs(counts)
    f <- decay_marginal(th, webs, gh, gradient = TRUE)
    num <- vapply(seq_along(th), function(k) {
      e <- replace(numeric(length(th)), k, 1e-5)
      (decay_marginal(th + e, webs, gh) - decay_marginal(th - e, webs, gh)) / 2e-5
    }, numeric(1))
    expect_equal(attr(f, "gradient"), num, tolerance = 1e-5)
  }
  check(counts, c(-1, 0.3, 0.2, 0.5, 0.1, 0.2, 0.3, 0.1, 0.4, 0.5, -0.2))
  # with the size of the web
  counts$web_nodes <- rep(c(20, 60, 150), each = 4)
  check(counts, c(-1, 0.3, 0.2, 0.5, 0.1, 0.2, 0.3, 0.1, 0.4, 0.5, -0.6, -0.2))
})

test_that("closer evidence never gives a lower probability and none is 0", {
  th <- c(-2, 0.4, 0, 1, 0.5, 0.2, 0, 0.7, 0.3, 0.8, log(0.9))
  tab <- decay_probability_table(th)
  p <- matrix(tab$p_new[!tab$generic], 5, 8, byrow = TRUE,
              dimnames = list(decay_consumer_levels, decay_levels))
  expect_true(all(p > 0 & p < 1))
  expect_true(all(diff(p) <= 1e-12))                # down the consumer levels
  expect_true(all(apply(p, 1, diff) <= 1e-12))      # down the resource levels
  # the levels within the family share one effect
  expect_equal(p[, "exact"], p[, "family"])
  # a generic record adds its effect
  ef <- decay_level_effects(th)
  expect_equal(decay_linear(th, 1, 5, TRUE) - decay_linear(th, 1, 5, FALSE), ef$generic)
})

test_that("fit_decay() recovers the parameters of simulated webs", {
  theta <- c(-1.5, 0.2, 0.6, 0.8, 0.5, 0.3, 0.4, 0.5, 0.6, 0.7, log(0.5))
  prof <- decay_all_profiles()
  use <- prof[prof$resource_level %in% c("family", "order", "class", "phylum", "any") &
                prof$consumer_level != "any", ]
  counts <- withr::with_seed(1, do.call(rbind, lapply(sprintf("w%02d", 1:15), function(w) {
    u <- stats::rnorm(1, 0, 0.5)
    m <- rep(300L, nrow(use))
    cl <- match(use$consumer_level, decay_consumer_levels)
    rg <- decay_resource_group(match(use$resource_level, decay_levels))
    y <- stats::rbinom(nrow(use), m, stats::plogis(decay_linear(theta, cl, rg, use$generic) + u))
    cbind(data.frame(study_id = w, stringsAsFactors = FALSE), use, pairs = m, observed = y)
  })))
  fit <- fit_decay(counts, method = "ml")
  ef_true <- unlist(decay_level_effects(theta))
  ef_fit <- fit$effects$estimate
  ok <- c(1:4, 6:11)          # no data for consumers resolved as "any"
  expect_lt(max(abs(ef_fit[ok] - ef_true[ok])), 0.25)
  expect_lt(abs(fit$sigma - 0.5), 0.3)
  expect_true(all(fit$probabilities$p_new > 0))
})

test_that("the Bayesian decay fit keeps evenly spaced draws within the bounds", {
  skip_on_cran()
  counts <- data.frame(study_id = rep(c("a", "b", "c", "d"), each = 3),
                       consumer_level = rep(c("family", "order", "phylum"), 4),
                       resource_level = rep(c("any", "class", "order"), 4), generic = FALSE,
                       pairs = rep(c(60, 30, 20), 4),
                       observed = c(12, 9, 2, 8, 5, 1, 15, 8, 3, 10, 4, 1))
  fit <- fit_decay(counts, n_chains = 2L, n_iter = 300L, n_warmup = 300L, seed = 1,
                   max_draws = 200L, thin = 2L)
  expect_identical(dim(fit$draws), c(200L, 11L))
  expect_true(all(fit$draws[, 2:9] >= 0))
  expect_identical(names(fit$rhat), decay_names())
  fit2 <- fit_decay(counts, n_chains = 2L, n_iter = 300L, n_warmup = 300L, seed = 1,
                    max_draws = 200L, thin = 2L)
  expect_identical(fit$draws, fit2$draws)
})

test_that("calibrate_links() fits the decay of the pairs not inferred", {
  v <- validate_rule(tiny_registry(), c("s1", "s2", "s3"), rule = tiny_rule(decay = TRUE),
                     quiet = TRUE)
  cal <- calibrate_links(v, method = "ml", quiet = TRUE)
  # 5 consumer levels, 5 resource levels, generic records and the size of the web
  expect_identical(nrow(cal$decay$effects), 12L)
  expect_identical(cal$decay$effects$side[12], "web")
  expect_true(all(v$decay$web_nodes == v$webs$n_nodes[match(v$decay$study_id, v$webs$study_id)]))
  expect_identical(nrow(cal$decay$probabilities), 80L)
  expect_identical(cal$decay$counts, v$decay)
  expect_identical(cal$decay$webs$study_id, cal$webs$study_id)
  expect_null(calibrate_links(validate_rule(tiny_registry(), c("s1", "s2", "s3"),
                                            rule = tiny_rule(), quiet = TRUE),
                              method = "ml", quiet = TRUE)$decay)
})

test_that("build_metaweb() lists the pairs without a link with their evidence", {
  st <- decay_setup()
  mw <- st$mw
  d <- mw$decay
  key <- function(x) paste(x$resource, x$consumer)
  consumers <- st$nodes$node_id[st$nodes$role == "consumer"]
  expect_true(all(d$consumer %in% consumers))
  expect_false(any(d$resource == d$consumer))
  expect_false(any(key(d) %in% key(mw$edges)))
  # every pair of a consumer is a link, cannibalism or a candidate
  n_all <- length(consumers) * nrow(st$nodes)
  n_links <- sum(mw$edges$consumer %in% consumers)
  n_cannibal_free <- sum(!consumers %in% mw$edges$consumer[mw$edges$cannibal])
  expect_identical(nrow(d) + n_links, as.integer(n_all - n_cannibal_free))
  # removed manual links and manual diets stay at 0
  ml <- data.frame(resource = mw$edges$resource[1], consumer = mw$edges$consumer[1],
                   action = "remove")
  mw2 <- build_metaweb(st$nodes, tiny_registry(), manual_links = ml)
  expect_false(key(ml) %in% key(mw2$decay))
  nodes3 <- prepare_nodes(tiny_taxa(), rule = tiny_rule(decay = TRUE, manual_diet = "gerr"))
  expect_false("gerr" %in% build_metaweb(nodes3, tiny_registry())$decay$consumer)
  expect_null(build_metaweb(prepare_nodes(tiny_taxa(), rule = tiny_rule()), tiny_registry())$decay)
})

test_that("sample_metawebs() draws the pairs without a link from the decay", {
  st <- decay_setup()
  dr <- sample_metawebs(st$mw, st$cal, n_draws = 2000, system = "average",
                        parameter_uncertainty = FALSE, seed = 1)
  bg <- dr$background
  expect_identical(nrow(bg$pairs), nrow(st$mw$decay))
  expect_identical(length(bg$present), 2000L)
  p <- decay_p_new_of(st$cal$decay, bg$pairs)
  # with the system effect at 0, the probability of each pair is p_average
  p_avg <- st$cal$decay$probabilities$p_average[
    decay_profile_id(bg$pairs$consumer_level, bg$pairs$resource_level, bg$pairs$generic)]
  expect_equal(bg$pairs$p_mean, p_avg)
  expect_lt(max(abs(bg$pairs$inclusion - p_avg)), 0.05)
  expect_true(all(p > 0))
  # the inferred links do not change with the decay
  mw0 <- st$mw
  mw0$decay <- NULL
  dr0 <- sample_metawebs(mw0, st$cal, n_draws = 2000, system = "average",
                         parameter_uncertainty = FALSE, seed = 1)
  expect_identical(dr0$z, dr$z)
  expect_null(dr0$background)
  # the indicators include the links drawn from the decay
  comm <- data.frame(sample_id = "s", node_id = st$nodes$node_id, value = 1)
  i1 <- calc_indicators(dr, comm)$values
  i0 <- calc_indicators(dr0, comm)$values
  expect_true(all(i1$L >= i0$L))
  expect_gt(mean(i1$L), mean(i0$L))
  # a calibration without the decay cannot sample it
  cal0 <- st$cal
  cal0$decay <- NULL
  expect_error(sample_metawebs(st$mw, cal0, n_draws = 2), "no decay model")
  expect_s3_class(sample_metawebs(mw0, cal0, n_draws = 2), "bwebs_draws")
})

test_that("score_metaweb() and web_brier() score the pairs without a link", {
  st <- decay_setup()
  e <- st$mw$edges[!st$mw$edges$cannibal, ]
  d <- st$mw$decay
  observed <- data.frame(resource = c(e$resource[1], d$resource[1]),
                         consumer = c(e$consumer[1], d$consumer[1]))
  s <- score_metaweb(st$mw, st$cal, observed)
  consumers <- unique(observed$consumer)
  n_pairs <- nrow(st$nodes) * length(consumers) - length(consumers)
  pr <- s$pairs
  expect_identical(sort(unique(pr$source)), c("decay", "rule"))
  bg <- d[d$consumer %in% consumers, ]
  expect_equal(sum(pr$source == "decay"), nrow(bg))
  expect_equal(s$summary$expected_background, sum(decay_p_new_of(st$cal$decay, bg)))
  expect_equal(s$summary$brier, sum((pr$p - pr$observed)^2) / n_pairs)
  # web_brier(): a recorded link not inferred costs (1 - p)^2 instead of 1
  test <- data.frame(evidence_class = "genus", predicted = 4, observed = 2)
  b0 <- web_brier(test, c(genus = 0.5), 20, 3, 4, 2, 0.1)
  bg2 <- data.frame(pairs = c(10, 6), observed = c(1, 0), p = c(0.2, 0.1))
  b1 <- web_brier(test, c(genus = 0.5), 20, 3, 4, 2, 0.1, bg2)
  inferred <- 2 * 0.5^2 + 2 * 0.5^2
  expect_equal(b0$brier, (inferred + 1) / 20)
  expect_equal(b1$brier, (inferred + 0.8^2 + 9 * 0.2^2 + 6 * 0.1^2) / 20)
  expect_equal(b1$expected_background, 10 * 0.2 + 6 * 0.1)
  expect_identical(b1$n_missed, 1L)
})
