# Posterior summaries of the indicators ----------------------------------------

#' Summarise indicators over draws
#'
#' `summarize_indicators()` summarises the indicators of [calc_indicators()]
#' over the draws, for every sample or for groups of samples. For a group,
#' the indicator is first averaged over the samples of the group within each
#' draw, and the draws of that mean are then summarised.
#'
#' @param x An object from [calc_indicators()].
#' @param indicators Names of the indicators to summarise; `NULL` (default)
#'   for all the numeric indicators.
#' @param samples Optional data frame with one row per sample: the column
#'   `sample_id` and the grouping columns named in `by`.
#' @param by Names of the columns of `samples` that define the groups;
#'   `NULL` (default) summarises every sample.
#' @param level Probability of the central credible interval (in (0, 1)).
#' @param na_rm If `FALSE` (default), a group or draw with a missing value
#'   (for example a sample left without links in a draw) gives a missing
#'   summary (`status = "missing_values"`). If `TRUE`, missing values are left
#'   out: `n_omitted` counts the values of samples left out of the group
#'   means and `n_missing` the draws without any value; the summary then
#'   describes the draws and samples where the indicator is defined.
#'
#' @return A data frame with one row per group (or sample) and indicator,
#'   ordered by group and indicator, with the grouping columns, `indicator`,
#'   `estimate` (mean over draws), `median`, `sd`, `lower` and `upper` (the
#'   `(1 - level) / 2` and `(1 + level) / 2` quantiles over draws), `level`,
#'   `interval` (`"credible"`, or `"none"` for the deterministic metaweb),
#'   `uncertainty` (components included), `n_draws`, `n_missing`,
#'   `n_omitted` and `status`.
#'
#' @details
#' The intervals describe the uncertainty of the links given the calibration
#' (and, with `system = "random"` in [sample_metawebs()], the variation
#' between systems, and with `parameter_uncertainty = TRUE` the uncertainty
#' of the calibration coefficients). They do not include the sampling error
#' of the community data.
#'
#' @export
#' @examples
#' nodes <- prepare_nodes(example_taxa)
#' mw <- build_metaweb(nodes, combine_registries(gray2015, gateway))
#' draws <- sample_metawebs(mw, gray2015_gateway_calibration, n_draws = 10, seed = 1)
#' ind <- calc_indicators(draws, example_community, permanent = c("cpom", "fpom"))
#' summarize_indicators(ind, indicators = c("C", "TL_mean"))
#' summarize_indicators(ind, indicators = "C", samples = example_samples, by = "treatment")
summarize_indicators <- function(x, indicators = NULL, samples = NULL, by = NULL,
                                 level = 0.95, na_rm = FALSE) {
  if (!inherits(x, "bwebs_indicators")) {
    abort_input("`x` must be created with calc_indicators().")
  }
  check_number(level, lower = 0, upper = 1, include_lower = FALSE, include_upper = FALSE)
  check_flag(na_rm)
  v <- x$values
  id_cols <- c("draw", "sample_id", "n_present", "n_missing")
  numeric_cols <- names(v)[vapply(v, is.numeric, logical(1))]
  available <- setdiff(numeric_cols, id_cols)
  if (is.null(indicators)) {
    indicators <- available
  } else {
    bad <- setdiff(indicators, available)
    if (!is.character(indicators) || length(bad)) {
      abort_input(sprintf("Unknown indicators: %s. Available: %s.",
                          format_items(bad), format_items(available, 30L)))
    }
  }
  if (is.null(by)) {
    if (!is.null(samples)) {
      abort_input("`samples` is only used with `by`; give the grouping columns in `by`.")
    }
    v$.group <- v$sample_id
    groups <- data.frame(sample_id = sort(unique(v$sample_id), method = "radix"),
                         stringsAsFactors = FALSE)
    group_key <- groups$sample_id
  } else {
    if (!is.character(by) || length(by) == 0L) abort_input("`by` must name columns of `samples`.")
    check_data_frame(samples, allow_empty = FALSE)
    check_columns(samples, c("sample_id", by), arg = "samples")
    samples <- as_plain_data_frame(samples[c("sample_id", by)])
    check_unique_key(samples, "sample_id", arg = "samples")
    miss <- setdiff(unique(v$sample_id), as.character(samples$sample_id))
    if (length(miss)) {
      abort_input(sprintf("`samples` lacks samples of `x`: %s.", format_items(miss)))
    }
    gk <- make_key(samples, by)
    v$.group <- gk[match(v$sample_id, as.character(samples$sample_id))]
    groups <- unique(samples[by])
    groups <- groups[do.call(order, c(unname(as.list(groups)), method = "radix")), , drop = FALSE]
    group_key <- make_key(groups, by)
    groups <- groups[group_key %in% v$.group, , drop = FALSE]
    group_key <- group_key[group_key %in% v$.group]
  }
  deterministic <- identical(x$settings$source, "metaweb")
  uncertainty <- if (deterministic) {
    "not_quantified"
  } else {
    paste(c("link_inference",
            if (identical(x$settings$system, "random")) "between_system",
            if (isTRUE(x$settings$parameter_uncertainty)) "parameter"), collapse = "+")
  }
  probs <- c((1 - level) / 2, (1 + level) / 2)
  rows <- vector("list", length(group_key) * length(indicators))
  k <- 0L
  for (g in seq_along(group_key)) {
    vg <- v[v$.group == group_key[g], , drop = FALSE]
    draws <- sort(unique(vg$draw))
    for (ind in indicators) {
      k <- k + 1L
      n_omitted <- if (na_rm) sum(is.na(vg[[ind]])) else 0L
      # mean over the samples of the group within each draw
      per_draw <- vapply(draws, function(d) {
        z <- vg[[ind]][vg$draw == d]
        if (anyNA(z) && !na_rm) return(NA_real_)
        z <- z[!is.na(z)]
        if (length(z) == 0L) NA_real_ else mean(z)
      }, numeric(1))
      n_na <- sum(is.na(per_draw))
      vals <- per_draw[!is.na(per_draw)]
      ok <- length(vals) > 0L && (n_na == 0L || na_rm)
      est <- med <- sdv <- lo <- hi <- NA_real_
      if (ok) {
        est <- mean(vals)
        med <- stats::median(vals)
        if (!deterministic && length(vals) > 1L) {
          sdv <- stats::sd(vals)
          q <- stats::quantile(vals, probs, names = FALSE)
          lo <- q[1]
          hi <- q[2]
        }
      }
      rows[[k]] <- data.frame(
        .g = g, indicator = ind, estimate = est, median = med, sd = sdv, lower = lo,
        upper = hi, level = if (deterministic) NA_real_ else level,
        interval = if (deterministic) "none" else "credible", uncertainty = uncertainty,
        n_draws = length(draws), n_missing = n_na, n_omitted = n_omitted,
        status = if (ok) "ok" else "missing_values", stringsAsFactors = FALSE
      )
    }
  }
  out <- do.call(rbind, rows)
  gcols <- groups[out$.g, , drop = FALSE]
  rownames(gcols) <- NULL
  out$.g <- NULL
  out <- cbind(gcols, out)
  rownames(out) <- NULL
  out
}

#' Contrasts of indicators between groups of samples
#'
#' `contrast_indicators()` compares groups of samples within each draw: for
#' every draw it averages the indicator over the samples of each group and
#' subtracts the mean of the reference group, so that the uncertainty of the
#' links, which is shared by all the samples of a draw, is carried into the
#' contrast. The draws of the difference are then summarised.
#'
#' @inheritParams summarize_indicators
#' @param samples Data frame with one row per sample: the column
#'   `sample_id`, the column named in `by` and those named in `within`.
#' @param by Name of the column of `samples` whose levels are compared.
#' @param reference Level of `by` subtracted from the others.
#' @param within Optional names of columns of `samples` that define strata
#'   (for example the sampling date): the contrast is computed within each
#'   stratum.
#'
#' @return A data frame with one row per stratum, compared level and
#'   indicator, with the columns of `within`, `level_compared`,
#'   `reference`, `indicator`, `estimate` (mean difference over draws),
#'   `median`, `sd`, `lower`, `upper`, `level`, `prob_positive` (proportion of
#'   draws with a positive difference), `n_draws`, `n_missing` (draws without
#'   a difference), `n_omitted` (sample values left out with `na_rm = TRUE`)
#'   and `status`.
#'   For the deterministic metaweb (one draw) only `estimate` and `median`
#'   are given.
#'
#' @export
#' @examples
#' nodes <- prepare_nodes(example_taxa)
#' mw <- build_metaweb(nodes, combine_registries(gray2015, gateway))
#' draws <- sample_metawebs(mw, gray2015_gateway_calibration, n_draws = 10, seed = 1)
#' ind <- calc_indicators(draws, example_community, permanent = c("cpom", "fpom"))
#' contrast_indicators(ind, samples = example_samples, by = "treatment",
#'                     reference = "control", indicators = c("C", "TL_mean"))
contrast_indicators <- function(x, samples, by, reference, within = NULL, indicators = NULL,
                                level = 0.95, na_rm = FALSE) {
  if (!inherits(x, "bwebs_indicators")) {
    abort_input("`x` must be created with calc_indicators().")
  }
  check_string(by)
  check_number(level, lower = 0, upper = 1, include_lower = FALSE, include_upper = FALSE)
  check_flag(na_rm)
  if (!is.null(within) && (!is.character(within) || by %in% within)) {
    abort_input("`within` must name columns of `samples` other than `by`.")
  }
  check_data_frame(samples, allow_empty = FALSE)
  check_columns(samples, c("sample_id", by, within), arg = "samples")
  samples <- as_plain_data_frame(samples[c("sample_id", by, within)])
  check_unique_key(samples, "sample_id", arg = "samples")
  lev_by <- as.character(samples[[by]])
  if (length(reference) != 1L || !as.character(reference) %in% lev_by) {
    abort_input(sprintf("`reference` must be one of the values of '%s': %s.",
                        by, format_items(sort(unique(lev_by), method = "radix"))))
  }
  reference <- as.character(reference)
  v <- x$values
  miss <- setdiff(unique(v$sample_id), as.character(samples$sample_id))
  if (length(miss)) {
    abort_input(sprintf("`samples` lacks samples of `x`: %s.", format_items(miss)))
  }
  id_cols <- c("draw", "sample_id", "n_present", "n_missing")
  available <- setdiff(names(v)[vapply(v, is.numeric, logical(1))], id_cols)
  if (is.null(indicators)) {
    indicators <- available
  } else if (!is.character(indicators) || length(setdiff(indicators, available))) {
    abort_input(sprintf("Unknown indicators: %s.", format_items(setdiff(indicators, available))))
  }
  m <- match(v$sample_id, as.character(samples$sample_id))
  v$.level <- lev_by[m]
  v$.stratum <- if (is.null(within)) rep("", nrow(v)) else make_key(samples[m, within, drop = FALSE], within)
  strata <- NULL
  strata_key <- ""
  if (!is.null(within)) {
    strata <- unique(samples[within])
    strata <- strata[do.call(order, c(unname(as.list(strata)), method = "radix")), , drop = FALSE]
    strata_key <- make_key(strata, within)
  }
  deterministic <- identical(x$settings$source, "metaweb")
  probs <- c((1 - level) / 2, (1 + level) / 2)
  group_mean <- function(z) {
    if (anyNA(z) && !na_rm) return(NA_real_)
    z <- z[!is.na(z)]
    if (length(z) == 0L) NA_real_ else mean(z)
  }
  rows <- list()
  for (s in seq_along(strata_key)) {
    vs <- v[v$.stratum == strata_key[s], , drop = FALSE]
    levs <- setdiff(sort(unique(vs$.level), method = "radix"), reference)
    if (!reference %in% vs$.level) next
    draws <- sort(unique(vs$draw))
    for (lv in levs) {
      for (ind in indicators) {
        used <- vs$.level %in% c(lv, reference)
        n_omitted <- if (na_rm) sum(is.na(vs[[ind]][used])) else 0L
        diff <- vapply(draws, function(d) {
          a <- group_mean(vs[[ind]][vs$draw == d & vs$.level == lv])
          b <- group_mean(vs[[ind]][vs$draw == d & vs$.level == reference])
          a - b
        }, numeric(1))
        n_na <- sum(is.na(diff))
        vals <- diff[!is.na(diff)]
        ok <- length(vals) > 0L && (n_na == 0L || na_rm)
        est <- med <- sdv <- lo <- hi <- pp <- NA_real_
        if (ok) {
          est <- mean(vals)
          med <- stats::median(vals)
          if (!deterministic && length(vals) > 1L) {
            sdv <- stats::sd(vals)
            q <- stats::quantile(vals, probs, names = FALSE)
            lo <- q[1]
            hi <- q[2]
            pp <- mean(vals > 0)
          }
        }
        rows[[length(rows) + 1L]] <- data.frame(
          .s = s, level_compared = lv, reference = reference, indicator = ind,
          estimate = est, median = med, sd = sdv, lower = lo, upper = hi,
          level = if (deterministic) NA_real_ else level, prob_positive = pp,
          n_draws = length(draws), n_missing = n_na, n_omitted = n_omitted,
          status = if (ok) "ok" else "missing_values", stringsAsFactors = FALSE
        )
      }
    }
  }
  if (!length(rows)) {
    abort_input("No stratum has both the reference and another level of `by`.")
  }
  out <- do.call(rbind, rows)
  if (!is.null(within)) {
    sc <- strata[out$.s, , drop = FALSE]
    rownames(sc) <- NULL
    out <- cbind(sc, out)
  }
  out$.s <- NULL
  rownames(out) <- NULL
  out
}
