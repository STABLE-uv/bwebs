# Condition helpers (internal) ---------------------------------------------
#
# From the template of the desarrollo-paquetes-r skill. crear_paquete.R
# replaces bwebs with the package name. Keep only what the package uses.
#
# The package keeps three situations apart:
#   * invalid input          -> abort_input(): classed error whose message says
#                               what is wrong, where (argument, column, rows)
#                               and how to fix it.
#   * estimate not available -> not an error: the result holds NA plus a
#                               `status` column (for example "empty_web").
#   * programming error      -> the natural R error; never turned into NA.
#
# Base R only. errorCondition() and warningCondition() need R >= 3.6.0.

# Signal an input error. `call` defaults to the call of the function that
# called abort_input(), so the user sees their own call in the message.
# Extra named arguments in `...` are stored in the condition object
# (e.g. column = "abundance", rows = c(3L, 8L)) for programmatic handling.
abort_input <- function(message, ..., class = character(), call = sys.call(-1)) {
  stop(errorCondition(
    message,
    ...,
    class = c(class, "bwebs_error_input", "bwebs_error"),
    call = call
  ))
}

# Optional dependencies (Suggests) ------------------------------------------

# Thin wrapper so that tests can simulate a missing package with
# testthat::local_mocked_bindings(is_installed = function(pkg) FALSE)
# (needs testthat >= 3.1.7).
is_installed <- function(pkg) {
  requireNamespace(pkg, quietly = TRUE)
}

# Stop with an informative message when an optional package is missing.
# `purpose` completes the sentence "Package 'x' is required to ...".
check_installed <- function(pkg, purpose, call = sys.call(-1)) {
  missing <- pkg[!vapply(pkg, is_installed, logical(1))]
  if (length(missing) == 0L) {
    return(invisible(TRUE))
  }
  stop(errorCondition(
    sprintf(
      "%s %s required to %s. Install %s with install.packages(%s).",
      if (length(missing) == 1L) "Package" else "Packages",
      paste0(quote_items(missing), collapse = " and "),
      purpose,
      if (length(missing) == 1L) "it" else "them",
      if (length(missing) == 1L) {
        paste0("\"", missing, "\"")
      } else {
        paste0("c(", paste0("\"", missing, "\"", collapse = ", "), ")")
      }
    ),
    package = missing,
    class = c("bwebs_error_missing_package", "bwebs_error"),
    call = call
  ))
}

# Message formatting ----------------------------------------------------------

# "a, b, c, d, e and 2 more"
format_items <- function(x, max = 5L) {
  x <- as.character(x)
  n <- length(x)
  if (n == 0L) {
    return("")
  }
  shown <- x[seq_len(min(n, max))]
  out <- paste(shown, collapse = ", ")
  if (n > max) {
    out <- paste0(out, " and ", n - max, " more")
  }
  out
}

# "row 3" / "rows 3, 8, 12 and 4 more"
format_rows <- function(rows, max = 5L) {
  paste0(if (length(rows) == 1L) "row " else "rows ", format_items(rows, max))
}

quote_items <- function(x) {
  paste0("'", x, "'")
}

# Text used in messages to name an argument from its expression.
arg_name <- function(expr) {
  paste(deparse(expr, width.cutoff = 500L), collapse = " ")
}

# Short description of an unexpected value, for messages.
describe_value <- function(x) {
  if (is.null(x)) {
    return("NULL")
  }
  if (is.character(x) && length(x) == 1L && !is.na(x)) {
    return(quote_items(x))
  }
  if (length(x) == 1L && (is.numeric(x) || is.logical(x))) {
    return(format(x))
  }
  sprintf("an object of class '%s' and length %d", class(x)[1L], length(x))
}
