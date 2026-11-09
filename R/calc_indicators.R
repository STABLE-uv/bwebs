# Indicators of the local food webs of a community ------------------------------

#' Topological indicators of the food web of every sample and draw
#'
#' `calc_indicators()` builds the local food web of every sample of a
#' community table (the nodes present in the sample plus the permanent
#' resources) from each draw of [sample_metawebs()] (its inferred links and,
#' with the taxonomic decay, the pairs without a link drawn in that draw), or
#' from the deterministic metaweb of [build_metaweb()] (inferred links only),
#' and computes its indicators with [web_indicators()].
#'
#' @param x A `bwebs_draws` object from [sample_metawebs()] or a
#'   `bwebs_metaweb` object from [build_metaweb()] (one deterministic web per
#'   sample, reported as draw 0).
#' @param community A data frame with one row per sample and node, in long
#'   format.
#' @param sample_col,node_col,value_col Names of the columns of `community`
#'   with the sample identifier, the node identifier (as in the metaweb) and
#'   the abundance or biomass (non-negative, any unit; only presence is used).
#' @param permanent Identifiers of nodes present in every sample, such as
#'   unsampled detritus or macrophytes. They are also treated as permanent in
#'   the robustness analysis.
#' @param threshold A node is present in a sample when its value is greater
#'   than `threshold` (default 0: any observed abundance or biomass).
#' @param na_rm If `FALSE` (default), a sample with a missing value (`NA`)
#'   gets no indicators (`status = "missing_values"`). If `TRUE`, the nodes
#'   with missing values are left out of that sample's web and counted in
#'   `n_missing`; the indicators then describe the web of the nodes with
#'   known values.
#' @param modularity,robustness,louvain_runs,robustness_reps,seed Passed to
#'   [web_indicators()]; `seed` is applied to every web.
#' @param quiet If `FALSE`, report progress with messages.
#'
#' @return An object of class `bwebs_indicators`, a list with `values`, a
#'   data frame with one row per draw and sample (ordered by draw and by
#'   sample in C-locale order) and columns `draw` (0 for a deterministic
#'   metaweb), `sample_id`, `n_present` (nodes detected), `n_missing` (nodes
#'   with missing values), the indicators of [web_indicators()] and `status`
#'   (`"ok"`, `"empty_web"` or `"missing_values"`); and `settings`. A
#'   community with zero rows gives zero rows.
#'
#' @section Missing values and zeros:
#' A zero is an observed absence. A node that has no row for a sample is
#' absent from it: give every sample its full list of nodes if absence
#' should not be assumed. `NA` is unknown and never treated as absence (see
#' `na_rm`).
#'
#' @section Assumptions:
#' Every node detected in a sample interacts as in the metaweb draw; the
#' indicators are topological (abundances and biomasses only define
#' presence).
#'
#' @seealso [summarize_indicators()].
#' @export
#' @examples
#' nodes <- prepare_nodes(example_taxa)
#' mw <- build_metaweb(nodes, combine_registries(gray2015, gateway))
#' draws <- sample_metawebs(mw, gray2015_gateway_calibration, n_draws = 5, seed = 1)
#' ind <- calc_indicators(draws, example_community, permanent = c("cpom", "fpom"))
#' head(ind$values[, c("draw", "sample_id", "S", "L", "C", "TL_mean")])
calc_indicators <- function(x, community, sample_col = "sample_id", node_col = "node_id",
                            value_col = "value", permanent = character(0), threshold = 0,
                            na_rm = FALSE, modularity = FALSE, robustness = FALSE,
                            louvain_runs = 10L, robustness_reps = 100L, seed = NULL,
                            quiet = TRUE) {
  if (!inherits(x, c("bwebs_draws", "bwebs_metaweb"))) {
    abort_input("`x` must be created with sample_metawebs() or build_metaweb().")
  }
  check_data_frame(community)
  check_string(sample_col)
  check_string(node_col)
  check_string(value_col)
  check_columns(community, c(sample_col = sample_col, node_col = node_col,
                             value_col = value_col), arg = "community")
  check_numeric_column(community, value_col, lower = 0, allow_na = TRUE, arg = "community")
  check_unique_key(community, c(sample_col, node_col), arg = "community")
  check_number(threshold, lower = 0)
  check_flag(na_rm)
  check_flag(modularity)
  check_flag(robustness)
  check_number(louvain_runs, lower = 1)
  check_number(robustness_reps, lower = 1)
  check_number(seed, allow_null = TRUE)
  check_flag(quiet)
  if (modularity) check_installed("igraph", "compute modularity")
  node_ids <- if (inherits(x, "bwebs_draws")) x$nodes$node_id else x$nodes$node_id
  samp <- as.character(community[[sample_col]])
  node <- as.character(community[[node_col]])
  val <- community[[value_col]]
  if (anyNA(samp) || anyNA(node)) {
    abort_input(sprintf("Columns '%s' and '%s' of `community` must not have missing values.",
                        sample_col, node_col))
  }
  unknown <- setdiff(node, node_ids)
  if (length(unknown)) {
    abort_input(sprintf(
      "`community` has nodes that are not in the metaweb: %s. Use the node_id of the nodes.",
      format_items(unknown)
    ))
  }
  if (!is.character(permanent) || anyNA(permanent) || length(setdiff(permanent, node_ids))) {
    abort_input(sprintf("`permanent` must contain node identifiers of the metaweb; unknown: %s.",
                        format_items(setdiff(permanent, node_ids))))
  }
  samples <- sort(unique(samp), method = "radix")
  # nodes of each sample, in the order of the metaweb nodes
  by_sample <- lapply(samples, function(s) {
    rows <- which(samp == s)
    present <- node[rows][!is.na(val[rows]) & val[rows] > threshold]
    missing <- node[rows][is.na(val[rows])]
    ids <- node_ids[node_ids %in% c(present, permanent)]
    list(ids = ids, n_present = length(present), n_missing = length(missing),
         blocked = length(missing) > 0L && !na_rm)
  })
  if (inherits(x, "bwebs_draws")) {
    links <- x$links
    bg <- x$background
    draws <- seq_len(ncol(x$z))
    adj_of <- function(d) {
      a <- matrix(0L, length(node_ids), length(node_ids), dimnames = list(node_ids, node_ids))
      on <- x$z[, d]
      if (any(on)) a[cbind(links$resource[on], links$consumer[on])] <- 1L
      if (!is.null(bg) && length(bg$present[[d]])) {
        b <- bg$present[[d]]
        a[cbind(bg$pairs$resource[b], bg$pairs$consumer[b])] <- 1L
      }
      a
    }
  } else {
    draws <- 0L
    adj_of <- function(d) x$adjacency
  }
  out <- vector("list", length(draws) * length(samples))
  k <- 0L
  for (d in draws) {
    if (!quiet) message(sprintf("Draw %d of %d", max(d, 1L), length(draws)))
    adj <- adj_of(d)
    for (s in seq_along(samples)) {
      k <- k + 1L
      info <- by_sample[[s]]
      ids <- info$ids
      met <- web_metrics_core(adj[ids, ids, drop = FALSE], permanent, modularity, robustness,
                              as.integer(louvain_runs), as.integer(robustness_reps), seed)
      if (info$blocked) {
        num <- vapply(met, is.numeric, logical(1))
        met[num] <- lapply(met[num], function(v) v * NA)
        met$TL_status <- NA_character_
        met$status <- "missing_values"
      }
      out[[k]] <- cbind(data.frame(draw = d, sample_id = samples[s], n_present = info$n_present,
                                   n_missing = info$n_missing, stringsAsFactors = FALSE), met)
    }
  }
  values <- if (k) do.call(rbind, out) else empty_indicator_values(modularity, robustness)
  rownames(values) <- NULL
  structure(
    list(values = values,
         settings = list(source = if (inherits(x, "bwebs_draws")) "draws" else "metaweb",
                         n_draws = length(draws), system = if (inherits(x, "bwebs_draws")) x$settings$system else NA,
                         parameter_uncertainty = if (inherits(x, "bwebs_draws")) x$settings$parameter_uncertainty else NA,
                         permanent = permanent, threshold = threshold, na_rm = na_rm,
                         modularity = modularity, robustness = robustness, seed = seed)),
    class = "bwebs_indicators"
  )
}

empty_indicator_values <- function(modularity, robustness) {
  met <- empty_web_metrics(0L, 0L, modularity, robustness)[0, , drop = FALSE]
  cbind(data.frame(draw = integer(), sample_id = character(), n_present = integer(),
                   n_missing = integer(), stringsAsFactors = FALSE), met)
}

#' @export
print.bwebs_indicators <- function(x, ...) {
  v <- x$values
  s <- x$settings
  cat(sprintf("<bwebs_indicators> %d sample(s) x %d draw(s) from %s\n",
              length(unique(v$sample_id)), s$n_draws,
              if (s$source == "draws") "sampled metawebs" else "the deterministic metaweb"))
  st <- table(v$status)
  cat("  status: ", paste(names(st), st, sep = " = ", collapse = ", "), "\n", sep = "")
  cat("  use summarize_indicators() for posterior summaries; values in $values\n")
  invisible(x)
}

#' @export
as.data.frame.bwebs_indicators <- function(x, ...) {
  x$values
}
