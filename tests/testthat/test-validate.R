# Tests for the input validators (template from the desarrollo-paquetes-r skill).
# Delete the blocks of validators the package does not use.

test_that("check_data_frame() accepts data frames and rejects other objects", {
  df <- data.frame(a = 1:2)
  expect_identical(check_data_frame(df), df)
  expect_error(check_data_frame(list(a = 1)), class = "bwebs_error_input")
  expect_error(check_data_frame(df[0, , drop = FALSE], allow_empty = FALSE), "has no rows")
  err <- expect_error(check_data_frame(1:3, arg = "data"), class = "bwebs_error_input")
  expect_match(conditionMessage(err), "`data` must be a data frame", fixed = TRUE)
})

test_that("check_data_frame() rejects grouped data frames instead of ignoring groups", {
  grouped <- structure(
    data.frame(site = c("a", "b"), n = 1:2),
    class = c("grouped_df", "tbl_df", "tbl", "data.frame")
  )
  expect_error(check_data_frame(grouped), "grouped data frame", class = "bwebs_error_input")
})

test_that("as_plain_data_frame() keeps names and resets row names", {
  df <- data.frame(`a b` = 1:3, check.names = FALSE)[c(3, 1), , drop = FALSE]
  out <- as_plain_data_frame(df)
  expect_identical(names(out), "a b")
  expect_identical(rownames(out), c("1", "2"))
})

test_that("check_columns() lists missing and available columns and the argument to use", {
  df <- data.frame(site = "a", count = 1)
  expect_identical(check_columns(df, c("site", "count")), df)
  err <- expect_error(
    check_columns(df, c(taxon_col = "taxon", value_col = "count")),
    class = "bwebs_error_input"
  )
  expect_match(conditionMessage(err), "Column 'taxon' not found in `data`", fixed = TRUE)
  expect_match(conditionMessage(err), "Available columns: site, count", fixed = TRUE)
  expect_match(conditionMessage(err), "`taxon_col`", fixed = TRUE)
  expect_identical(err$missing_columns, "taxon")
})

test_that("check_numeric_column() explains decimal commas and factors", {
  expect_error(
    check_numeric_column(data.frame(x = c("1,5", "2")), "x"),
    "decimal comma",
    class = "bwebs_error_input"
  )
  expect_error(
    check_numeric_column(data.frame(x = factor(c("10", "5"))), "x"),
    "level codes"
  )
})

test_that("check_numeric_column() reports NA, NaN, Inf and their rows", {
  err <- expect_error(check_numeric_column(data.frame(x = c(1, NA, NaN)), "x"))
  expect_match(conditionMessage(err), "NA in row 2; NaN in row 3", fixed = TRUE)
  expect_identical(err$rows, c(2L, 3L))
  expect_identical(check_numeric_column(data.frame(x = c(1, NA)), "x", allow_na = TRUE)$x, c(1, NA))
  expect_error(check_numeric_column(data.frame(x = c(1, Inf)), "x"), "infinite values in row 2")
})

test_that("check_numeric_column() applies inclusive and exclusive bounds", {
  df <- data.frame(volume = c(0, 1, 2))
  expect_identical(check_numeric_column(df, "volume", lower = 0), df)
  err <- expect_error(check_numeric_column(df, "volume", lower = 0, include_lower = FALSE))
  expect_match(conditionMessage(err), "must be in (0, Inf)", fixed = TRUE)
  expect_identical(err$rows, 1L)
  expect_error(check_numeric_column(df, "volume", upper = 1), "row 3")
})

test_that("check_numeric_column() can require whole numbers", {
  expect_error(
    check_numeric_column(data.frame(n = c(1, 2.5)), "n", whole_numbers = TRUE),
    "whole numbers"
  )
  expect_silent(check_numeric_column(data.frame(n = c(1, 2)), "n", whole_numbers = TRUE))
})

test_that("check_unique_key() reports every row of a duplicated key", {
  df <- data.frame(site = c("a", "a", "b"), taxon = c("x", "x", "x"))
  err <- expect_error(check_unique_key(df, c("site", "taxon")), class = "bwebs_error_input")
  expect_identical(err$rows, c(1L, 2L))
  expect_identical(check_unique_key(df[2:3, ], c("site", "taxon")), df[2:3, ])
})

test_that("check_choice() matches exactly and rejects abbreviations and vectors", {
  method <- "linear"
  expect_identical(check_choice(method, c("linear", "log")), "linear")
  method <- "lin"
  err <- expect_error(check_choice(method, c("linear", "log")), class = "bwebs_error_input")
  expect_match(conditionMessage(err), "`method` must be one of 'linear', 'log', not 'lin'", fixed = TRUE)
  expect_error(check_choice(c("linear", "log"), c("linear", "log"), arg = "method"))
  expect_error(check_choice(NA_character_, "a", arg = "method"))
})

test_that("check_flag() accepts only TRUE or FALSE", {
  expect_error(check_flag(NA, arg = "quiet"), "TRUE or FALSE")
  expect_error(check_flag(c(TRUE, FALSE), arg = "quiet"))
  expect_identical(check_flag(TRUE), TRUE)
})

test_that("check_number() rejects wrong types and values outside its bounds", {
  expect_error(check_number("1", arg = "level"), "single finite number")
  expect_error(check_number(1, lower = 0, upper = 1, include_upper = FALSE, arg = "level"), "must be in [0, 1)", fixed = TRUE)
  expect_null(check_number(NULL, allow_null = TRUE))
})

test_that("check_string() rejects empty strings", {
  expect_error(check_string("", arg = "unit"), "non-empty string")
  expect_identical(check_string("mg/L"), "mg/L")
})

