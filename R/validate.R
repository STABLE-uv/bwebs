# Input validators (internal) -------------------------------------------------
#
# From the template of the desarrollo-paquetes-r skill. Call them at the start
# of exported functions; internal helpers assume validated input. Before
# delivery, delete the validators no function uses, together with their tests.
# Needs conditions.R (abort_input(), format_items(), format_rows(), ...).
#
# Every validator returns its input invisibly, so it can be used as a
# statement, and accepts `call` so the error shows the user's call.

# `x` must be a data frame (tibbles are data frames). Grouped data frames are
# rejected because their groups would be silently ignored.
check_data_frame <- function(x, arg = arg_name(substitute(x)), allow_empty = TRUE,
                             call = sys.call(-1)) {
  force(arg)
  if (!is.data.frame(x)) {
    abort_input(
      sprintf("`%s` must be a data frame, not %s.", arg, describe_value(x)),
      call = call
    )
  }
  if (inherits(x, "grouped_df")) {
    abort_input(
      sprintf(
        paste(
          "`%s` is a grouped data frame. Its groups would be ignored here;",
          "ungroup it first (dplyr::ungroup()) and, if the function has a `by`",
          "argument, give the grouping columns there."
        ),
        arg
      ),
      call = call
    )
  }
  if (!allow_empty && nrow(x) == 0L) {
    abort_input(sprintf("`%s` has no rows.", arg), call = call)
  }
  invisible(x)
}

# Tibbles and data.tables become plain data frames with unchanged column names
# and default row names.
as_plain_data_frame <- function(x) {
  x <- as.data.frame(x)
  rownames(x) <- NULL
  x
}

# Required columns. Name the elements of `columns` after the arguments that
# supply them, e.g. c(taxon_col = "taxon", value_col = "count"), so the
# message can say which argument to change.
check_columns <- function(data, columns, arg = "data", call = sys.call(-1)) {
  missing <- columns[!columns %in% names(data)]
  if (length(missing) == 0L) {
    return(invisible(data))
  }
  hint <- ""
  arg_names <- names(missing)
  if (!is.null(arg_names) && any(nzchar(arg_names))) {
    arg_names <- arg_names[nzchar(arg_names)]
    hint <- sprintf(
      " Use %s to give the name of the column to use.",
      format_items(paste0("`", arg_names, "`"), max = Inf)
    )
  }
  abort_input(
    sprintf(
      "%s %s not found in `%s`. Available columns: %s.%s",
      if (length(missing) == 1L) "Column" else "Columns",
      format_items(quote_items(missing), max = Inf),
      arg,
      if (ncol(data) == 0L) "none" else format_items(names(data), max = 20L),
      hint
    ),
    missing_columns = unname(missing),
    call = call
  )
}

# A numeric column with optional range, NA, Inf and whole-number rules.
# The interval is [lower, upper] by default; use include_lower = FALSE for
# strictly positive quantities (volumes, areas, sampling effort).
check_numeric_column <- function(data, column, lower = -Inf, upper = Inf,
                                 include_lower = TRUE, include_upper = TRUE,
                                 allow_na = FALSE, allow_infinite = FALSE,
                                 whole_numbers = FALSE, arg = "data",
                                 call = sys.call(-1)) {
  x <- data[[column]]
  label <- sprintf("Column '%s' of `%s`", column, arg)
  if (!is.numeric(x)) {
    hint <- ""
    if (is.factor(x)) {
      hint <- paste(
        " It is a factor: as.numeric() would return the level codes,",
        "not the values; convert it with as.numeric(as.character(x))."
      )
    } else if (is.character(x) && any(grepl("^\\s*-?[0-9]+,[0-9]+\\s*$", x))) {
      hint <- paste(
        " Some values look like numbers with a decimal comma (such as '1,5').",
        "Check the source file and read it with dec = \",\" (or read.csv2())."
      )
    }
    abort_input(
      sprintf("%s must be numeric, not %s.%s", label, class(x)[1L], hint),
      column = column,
      call = call
    )
  }
  if (!allow_na) {
    na_rows <- which(is.na(x) & !is.nan(x))
    nan_rows <- which(is.nan(x))
    if (length(na_rows) > 0L || length(nan_rows) > 0L) {
      parts <- c(
        if (length(na_rows) > 0L) paste("NA in", format_rows(na_rows)),
        if (length(nan_rows) > 0L) paste("NaN in", format_rows(nan_rows))
      )
      abort_input(
        sprintf(
          paste(
            "%s has missing values (%s). Missing values are not removed",
            "automatically: decide what they mean (not sampled, unknown,",
            "below the detection limit) and handle them before this step."
          ),
          label, paste(parts, collapse = "; ")
        ),
        column = column,
        rows = c(na_rows, nan_rows),
        call = call
      )
    }
  }
  if (!allow_infinite) {
    inf_rows <- which(is.infinite(x))
    if (length(inf_rows) > 0L) {
      abort_input(
        sprintf("%s has infinite values in %s.", label, format_rows(inf_rows)),
        column = column,
        rows = inf_rows,
        call = call
      )
    }
  }
  finite <- !is.na(x) & is.finite(x)
  below <- if (include_lower) x < lower else x <= lower
  above <- if (include_upper) x > upper else x >= upper
  out_rows <- which(finite & (below | above))
  if (length(out_rows) > 0L) {
    abort_input(
      sprintf(
        "%s must be in %s; %s outside this range (%s).",
        label,
        format_interval(lower, upper, include_lower, include_upper),
        if (length(out_rows) == 1L) "1 value is" else paste(length(out_rows), "values are"),
        format_rows(out_rows)
      ),
      column = column,
      rows = out_rows,
      call = call
    )
  }
  if (whole_numbers) {
    frac_rows <- which(finite & abs(x - round(x)) > sqrt(.Machine$double.eps))
    if (length(frac_rows) > 0L) {
      abort_input(
        sprintf(
          "%s must contain whole numbers (counts); found decimals in %s.",
          label, format_rows(frac_rows)
        ),
        column = column,
        rows = frac_rows,
        call = call
      )
    }
  }
  invisible(data)
}

# One row per combination of the key columns.
check_unique_key <- function(data, key, arg = "data", call = sys.call(-1)) {
  dup <- duplicated(data[key]) | duplicated(data[key], fromLast = TRUE)
  if (!any(dup)) {
    return(invisible(data))
  }
  rows <- which(dup)
  abort_input(
    sprintf(
      paste(
        "`%s` must have one row per combination of %s, but %d rows share a",
        "combination (%s). Aggregate or remove the duplicates first."
      ),
      arg, format_items(quote_items(key), max = Inf), length(rows), format_rows(rows)
    ),
    rows = rows,
    call = call
  )
}

# One option from a closed list, matched exactly (no abbreviations, no
# implicit default taken from the first element).
check_choice <- function(x, choices, arg = arg_name(substitute(x)),
                         call = sys.call(-1)) {
  force(arg)
  if (!is.character(x) || length(x) != 1L || is.na(x) || !x %in% choices) {
    abort_input(
      sprintf(
        "`%s` must be one of %s, not %s.",
        arg, format_items(quote_items(choices), max = Inf), describe_value(x)
      ),
      call = call
    )
  }
  invisible(x)
}

# TRUE or FALSE.
check_flag <- function(x, arg = arg_name(substitute(x)), call = sys.call(-1)) {
  force(arg)
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    abort_input(
      sprintf("`%s` must be TRUE or FALSE, not %s.", arg, describe_value(x)),
      call = call
    )
  }
  invisible(x)
}

# A single finite number in [lower, upper] (bounds optionally exclusive).
check_number <- function(x, lower = -Inf, upper = Inf, include_lower = TRUE,
                         include_upper = TRUE, allow_null = FALSE,
                         arg = arg_name(substitute(x)), call = sys.call(-1)) {
  force(arg)
  if (allow_null && is.null(x)) {
    return(invisible(x))
  }
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x)) {
    abort_input(
      sprintf("`%s` must be a single finite number, not %s.", arg, describe_value(x)),
      call = call
    )
  }
  below <- if (include_lower) x < lower else x <= lower
  above <- if (include_upper) x > upper else x >= upper
  if (below || above) {
    abort_input(
      sprintf(
        "`%s` must be in %s, not %s.",
        arg, format_interval(lower, upper, include_lower, include_upper), format(x)
      ),
      call = call
    )
  }
  invisible(x)
}

# A single non-empty string (column names, units, labels).
check_string <- function(x, allow_null = FALSE, arg = arg_name(substitute(x)),
                         call = sys.call(-1)) {
  force(arg)
  if (allow_null && is.null(x)) {
    return(invisible(x))
  }
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    abort_input(
      sprintf("`%s` must be a single non-empty string, not %s.", arg, describe_value(x)),
      call = call
    )
  }
  invisible(x)
}

# "[0, 1]", "(0, Inf)": infinite bounds are always open.
format_interval <- function(lower, upper, include_lower = TRUE, include_upper = TRUE) {
  paste0(
    if (include_lower && is.finite(lower)) "[" else "(",
    format(lower), ", ", format(upper),
    if (include_upper && is.finite(upper)) "]" else ")"
  )
}

make_key <- function(data, by) {
  cols <- lapply(data[by], function(v) if (is.factor(v)) as.character(v) else v)
  do.call(paste, c(unname(cols), sep = "\u001f"))
}
