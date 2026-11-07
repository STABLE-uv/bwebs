# Tests for the condition helpers (template from the desarrollo-paquetes-r skill).

test_that("abort_input() signals a classed error that names the caller", {
  f <- function(x) abort_input("`x` must be positive.", argument = "x")
  err <- expect_error(f(-1), class = "bwebs_error_input")
  expect_s3_class(err, "bwebs_error")
  expect_identical(conditionCall(err), quote(f(-1)))
  expect_identical(conditionMessage(err), "`x` must be positive.")
  expect_identical(err$argument, "x")
})

test_that("check_installed() explains how to install a missing package", {
  expect_true(check_installed("stats", "fit models"))
  err <- expect_error(
    check_installed("bwebsnotapackage", "draw the figure"),
    class = "bwebs_error_missing_package"
  )
  expect_match(conditionMessage(err), "required to draw the figure", fixed = TRUE)
  expect_match(conditionMessage(err), 'install.packages("bwebsnotapackage")', fixed = TRUE)
  expect_identical(err$package, "bwebsnotapackage")
})

test_that("format_items() abbreviates long lists", {
  expect_identical(format_items(character()), "")
  expect_identical(format_items(letters[1:3]), "a, b, c")
  expect_identical(format_items(letters[1:7], max = 5), "a, b, c, d, e and 2 more")
})

test_that("format_rows() names one row or several", {
  expect_identical(format_rows(3L), "row 3")
  expect_identical(format_rows(c(3L, 8L)), "rows 3, 8")
})
