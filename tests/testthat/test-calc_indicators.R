tiny_community <- function() {
  data.frame(
    sample_id = c("s1", "s1", "s1", "s2", "s2", "s2"),
    node_id = c("cocc", "baet", "perl", "cocc", "baet", "perl"),
    value = c(10, 3, 1, 8, 2, 0)
  )
}

test_that("deterministic webs contain the present nodes and the permanent resources", {
  mw <- build_metaweb(prepare_nodes(tiny_taxa(), rule = tiny_rule()), tiny_registry())
  ind <- calc_indicators(mw, tiny_community(), permanent = "cpom")
  v <- ind$values
  expect_identical(v$draw, c(0L, 0L))
  expect_identical(v$sample_id, c("s1", "s2"))
  # s1: cocc -> baet -> perl plus isolated cpom; s2: cocc -> baet plus cpom
  expect_identical(v$S_taxa, c(4L, 3L))
  expect_identical(v$n_isolated, c(1L, 1L))
  expect_identical(v$L, c(2L, 1L))
  expect_equal(v$TL_max, c(3, 2))
  expect_identical(v$n_present, c(3L, 2L))
})

test_that("missing values block a sample unless na_rm = TRUE", {
  mw <- build_metaweb(prepare_nodes(tiny_taxa(), rule = tiny_rule()), tiny_registry())
  com <- tiny_community()
  com$value[3] <- NA
  blocked <- calc_indicators(mw, com)
  expect_identical(blocked$values$status, c("missing_values", "ok"))
  expect_true(is.na(blocked$values$C[1]))
  kept <- calc_indicators(mw, com, na_rm = TRUE)
  expect_identical(kept$values$n_missing, c(1L, 0L))
  expect_identical(kept$values$S_taxa[1], 2L)
})

test_that("results do not depend on row order and zero rows give zero rows", {
  mw <- build_metaweb(prepare_nodes(tiny_taxa(), rule = tiny_rule()), tiny_registry())
  com <- tiny_community()
  a <- calc_indicators(mw, com)$values
  b <- calc_indicators(mw, com[rev(seq_len(nrow(com))), ])$values
  expect_identical(a, b)
  z <- calc_indicators(mw, com[0, ])$values
  expect_identical(nrow(z), 0L)
  expect_identical(names(z), names(a))
})

test_that("draws give one row per draw and sample", {
  mw <- build_metaweb(prepare_nodes(tiny_taxa(), rule = tiny_rule()), tiny_registry())
  d <- sample_metawebs(mw, fake_calibration(c(genus = 3, family = 3, order = 3)), n_draws = 4,
                       system = "average", parameter_uncertainty = FALSE, seed = 1)
  ind <- calc_indicators(d, tiny_community(), permanent = "cpom")
  expect_identical(nrow(ind$values), 8L)
  expect_identical(ind$values$draw, rep(1:4, each = 2))
})

test_that("calc_indicators() rejects invalid input with classed errors", {
  mw <- build_metaweb(prepare_nodes(tiny_taxa(), rule = tiny_rule()), tiny_registry())
  com <- tiny_community()
  expect_error(calc_indicators(mw, rbind(com, com)), "one row per combination")
  bad <- com; bad$node_id[1] <- "zzz"
  expect_error(calc_indicators(mw, bad), "not in the metaweb")
  bad <- com; bad$value[1] <- -1
  expect_error(calc_indicators(mw, bad), class = "bwebs_error_input")
  expect_error(calc_indicators(mw, com, permanent = "zzz"), "permanent")
  expect_error(calc_indicators(mw, com, value_col = "abundance"), "value_col")
  expect_error(calc_indicators(list(), com), class = "bwebs_error_input")
})

test_that("print method is stable", {
  mw <- build_metaweb(prepare_nodes(tiny_taxa(), rule = tiny_rule()), tiny_registry())
  expect_snapshot(print(calc_indicators(mw, tiny_community())))
})
