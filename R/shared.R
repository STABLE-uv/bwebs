# Studies that share links with a web ----------------------------------------------

#' Find studies that share links with a web
#'
#' `shared_studies()` measures, for each web, the proportion of its links
#' that are also recorded in each study of a registry. A study holding a
#' large share of the links of a validation web probably describes the same
#' system (or copies its data), and its records would predict that web when
#' it is left out; such studies belong in the `groups` of [validate_rule()],
#' or out of the registry when the web is an external test.
#'
#' @param registry A registry in the format of [as_registry()].
#' @param webs Either identifiers (`study_id`) of studies of `registry`, or a
#'   data frame of links with columns `web`, `resource` and `consumer` (names
#'   comparable with those of the registry, such as the taxa of an observed
#'   web).
#' @param threshold Smallest proportion of the links of a web reported
#'   (between 0 and 1).
#'
#' @return A data frame with one row per web and study above `threshold`:
#'   `web`, `study_id`, `n_links` (links of the web), `n_shared` and `share`,
#'   from the largest share. Its attribute `groups` lists the connected
#'   groups of studies (each web with the studies that share its links), in
#'   the format of `groups` in [validate_rule()] (only when `webs` are study
#'   identifiers).
#'
#' @details
#' Links are compared by their harmonised names in lower case
#' (`resource_std`, `consumer_std`), without cannibalism. The share is
#' asymmetric: a small web can be contained in a large compilation.
#'
#' @seealso [validate_rule()].
#' @export
#' @examples
#' val <- gray2015_sources$study_id[gray2015_sources$validation]
#' shared_studies(gray2015, val, threshold = 0.2)
shared_studies <- function(registry, webs, threshold = 0.2) {
  check_registry(registry)
  check_number(threshold, lower = 0, upper = 1)
  key <- function(r, c) paste(tolower(r), tolower(c), sep = "\u001f")
  keep <- tolower(registry$resource_std) != tolower(registry$consumer_std)
  pairs <- data.frame(study_id = registry$study_id[keep],
                      key = key(registry$resource_std[keep], registry$consumer_std[keep]),
                      stringsAsFactors = FALSE)
  pairs <- pairs[!duplicated(pairs), , drop = FALSE]
  if (is.character(webs)) {
    unknown <- setdiff(webs, registry$study_id)
    if (length(unknown)) {
      abort_input(sprintf("`webs` has studies that are not in `registry`: %s.", format_items(unknown)))
    }
    web_links <- lapply(stats::setNames(webs, webs), function(w) pairs$key[pairs$study_id == w])
  } else {
    check_data_frame(webs)
    check_columns(webs, c("web", "resource", "consumer"), arg = "webs")
    w <- webs[tolower(webs$resource) != tolower(webs$consumer), , drop = FALSE]
    web_links <- lapply(split(key(w$resource, w$consumer), w$web), unique)
  }
  rows <- lapply(names(web_links), function(wn) {
    lk <- web_links[[wn]]
    hit_rows <- pairs$study_id != wn & pairs$key %in% lk
    if (!any(hit_rows) || !length(lk)) return(NULL)
    n <- table(pairs$study_id[hit_rows])
    share <- as.numeric(n) / length(lk)
    hit <- share >= threshold
    if (!any(hit)) return(NULL)
    data.frame(web = wn, study_id = names(n)[hit], n_links = length(lk),
               n_shared = as.integer(n[hit]), share = share[hit], stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  if (is.null(out)) {
    out <- data.frame(web = character(), study_id = character(), n_links = integer(),
                      n_shared = integer(), share = numeric(), stringsAsFactors = FALSE)
  }
  out <- out[order(-out$share, out$web, out$study_id, method = "radix"), , drop = FALSE]
  rownames(out) <- NULL
  if (is.character(webs)) attr(out, "groups") <- connected_groups(out$web, out$study_id)
  out
}

# Connected components of the pairs (a[i], b[i]), as a list of sorted
# character vectors with at least two members.
connected_groups <- function(a, b) {
  if (!length(a)) return(list())
  ids <- sort(unique(c(a, b)), method = "radix")
  comp <- stats::setNames(seq_along(ids), ids)
  repeat {
    changed <- FALSE
    for (i in seq_along(a)) {
      m <- min(comp[[a[i]]], comp[[b[i]]])
      if (comp[[a[i]]] != m || comp[[b[i]]] != m) {
        old <- c(comp[[a[i]]], comp[[b[i]]])
        comp[comp %in% old] <- m
        changed <- TRUE
      }
    }
    if (!changed) break
  }
  g <- split(names(comp), comp)
  g <- lapply(g, sort, method = "radix")
  unname(g[vapply(g, length, integer(1)) >= 2L])
}
