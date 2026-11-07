hand_indicators <- function(source = "draws") {
  structure(
    list(values = data.frame(draw = rep(1:3, each = 2), sample_id = rep(c("a", "b"), 3),
                             n_present = 2L, n_missing = 0L,
                             C = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6), status = "ok"),
         settings = list(source = source, system = "average", parameter_uncertainty = FALSE)),
    class = "bwebs_indicators"
  )
}

test_that("summaries over draws match a hand calculation", {
  s <- summarize_indicators(hand_indicators(), level = 0.5)
  a <- s[s$sample_id == "a", ]
  # values of sample a: 0.1, 0.3, 0.5
  expect_equal(c(a$estimate, a$median, a$sd, a$lower, a$upper), c(0.3, 0.3, 0.2, 0.2, 0.4))
  expect_identical(a$interval, "credible")
  expect_identical(a$uncertainty, "link_inference")
  grp <- data.frame(sample_id = c("a", "b"), habitat = c("x", "x"))
  g <- summarize_indicators(hand_indicators(), samples = grp, by = "habitat")
  # per-draw means 0.15, 0.35, 0.55
  expect_equal(g$estimate, 0.35)
  expect_identical(g$habitat, "x")
})

test_that("contrasts are computed within draws", {
  grp <- data.frame(sample_id = c("a", "b"), trt = c("ctrl", "exp"))
  k <- contrast_indicators(hand_indicators(), samples = grp, by = "trt", reference = "ctrl")
  # b - a is 0.1 in every draw
  expect_equal(c(k$estimate, k$sd, k$prob_positive), c(0.1, 0, 1))
  expect_identical(c(k$level_compared, k$reference), c("exp", "ctrl"))
  expect_error(contrast_indicators(hand_indicators(), samples = grp, by = "trt", reference = "zz"),
               "reference")
})

test_that("missing values are reported, not dropped, by default", {
  x <- hand_indicators()
  x$values$C[1] <- NA
  s <- summarize_indicators(x)
  expect_identical(s$status[s$sample_id == "a"], "missing_values")
  s2 <- summarize_indicators(x, na_rm = TRUE)
  expect_equal(s2$estimate[s2$sample_id == "a"], 0.4)
  expect_identical(s2$n_missing[s2$sample_id == "a"], 1L)
  # in a group, the missing sample is left out of the draw mean and counted
  grp <- data.frame(sample_id = c("a", "b"), habitat = c("x", "x"))
  g <- summarize_indicators(x, samples = grp, by = "habitat", na_rm = TRUE)
  expect_identical(c(g$n_missing, g$n_omitted), c(0L, 1L))
  expect_equal(g$estimate, mean(c(0.2, 0.35, 0.55)))
})

test_that("a deterministic metaweb has no interval", {
  x <- hand_indicators(source = "metaweb")
  x$values <- x$values[x$values$draw == 1, ]
  s <- summarize_indicators(x)
  expect_identical(unique(s$interval), "none")
  expect_true(all(is.na(s$lower)))
})

test_that("summaries reject unknown indicators and incomplete sample tables", {
  expect_error(summarize_indicators(hand_indicators(), indicators = "nope"), "Unknown")
  expect_error(summarize_indicators(hand_indicators(), samples = data.frame(sample_id = "a", g = 1),
                                    by = "g"), "lacks samples")
})
