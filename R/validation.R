# Validation of a matching rule on held-out webs ---------------------------------

#' Validate a matching rule on held-out webs of the registry
#'
#' `validate_rule()` rebuilds each validation study of the registry from the
#' rest of the registry and compares the inferred links with the recorded
#' ones: the taxa of the study become the nodes, the metaweb is built with
#' the registry without that study (and without the other studies of its
#' system, see `groups`) and the same [matching_rule()], and every inferred
#' link is classified by its class of evidence and by whether the study
#' recorded it, as in the validation of Gray et al. (2015, section 2.3,
#' p. 15). The result feeds [calibrate_links()] (link probabilities) and
#' [score_rule()] (cross-validated scores of the rule).
#'
#' @param registry A registry in the format of [as_registry()].
#' @param sources Identifiers (`study_id`) of the validation studies: whole
#'   communities in which the diets of the consumers were recorded. No
#'   default: which studies qualify is a scientific decision. For the bundled
#'   registry, `gray2015_sources$study_id[gray2015_sources$validation]`
#'   gives a provisional set.
#' @param groups `NULL` (each validation study is left out alone) or a list of
#'   character vectors of study identifiers, each listing studies of the same
#'   system, which share records. When a validation study is left out, the
#'   other studies of its group are left out with it, so that records of the
#'   same system do not predict it; studies that are not validation webs can
#'   be listed too. Which studies describe the same system is a scientific
#'   decision; see [gray2015_calibration] for the groups used with the bundled
#'   registry.
#' @param rule A [matching_rule()].
#' @param train_registries `NULL` (default: every registry of `registry`) or
#'   names of the registries (column `registry`) used to rebuild the
#'   validation webs, for example `"globi"` to validate GloBI alone on webs of
#'   [gray2015]. The validation webs always come from `registry`.
#' @param exclusions `NULL` or a data frame with columns `web` and `study_id`,
#'   such as the result of [shared_studies()]: studies left out of the
#'   registry when `web` is rebuilt, in addition to the studies of its group
#'   (for example copies of that web in other compilations).
#' @param quiet If `FALSE`, report progress with messages.
#'
#' @return An object of class `bwebs_validation`, a list with
#' * `counts`: one row per web and class of evidence (`study_id`,
#'   `evidence_class`, `predicted` links and how many of them were
#'   `observed`);
#' * `webs`: one row per validation web with the number of nodes, consumers
#'   with recorded diet, pairs evaluated, recorded links, predicted links,
#'   predicted links observed, precision, recall, whether the matching
#'   levels converged, the studies of its group left out to rebuild it
#'   (`left_out`) and the number of other studies left out (`n_excluded`, see
#'   `exclusions`);
#' * `links`: one row per link inferred in each validation web (`study_id`,
#'   `resource`, `consumer`, `evidence_class`, the amount of evidence of the
#'   link as in the `edges` of [build_metaweb()]: `n_records`, `n_sources`,
#'   `n_consumer_taxa`, `n_resource_taxa`, `n_consumer_known`,
#'   `n_resource_known`; and whether the web recorded it, `observed`);
#' * `decay` (with `decay = TRUE` in the rule, otherwise `NULL`): the pairs
#'   without an inferred link of each web by their evidence for the
#'   taxonomic decay (`study_id`, `consumer_level`, `resource_level`,
#'   `generic`; see [matching_rule()]), with the number of `pairs`, how
#'   many of them were `observed` (recorded links that the rule missed) and
#'   the number of nodes of the web (`web_nodes`);
#' * `settings`: rule, sources, groups, registries and number of records.
#'
#' @details
#' For the validation web \eqn{w}, the pairs evaluated are the resource
#' \eqn{\to} consumer pairs among its taxa whose consumer has its diet
#' recorded in \eqn{w}, without cannibalism. The taxa of a validation web are
#' classified with the taxonomy of its records (the first classification
#' found for each name); a name equal to its genus is a genus, a name
#' starting with its genus is a species, and a name without any rank is a
#' non-taxonomic node. Taxa that are only eaten in the web are producers when
#' the registry says so, detritus when they have no classification and
#' consumers otherwise; their diets are not evaluated. Nodes with a manual
#' diet in `rule` are matched like any other node. With `decay = TRUE` in the
#' rule, the evidence of the pairs without an inferred link comes from the
#' same records (the registry without the studies left out).
#'
#' @section Assumptions:
#' No study left in the registry describes the same system as the web left
#' out; otherwise its records predict the web (list such studies in
#' `groups`).
#'
#' @references
#' Gray C, Figueroa DH, Hudson LN, Ma A, Perkins D, Woodward G (2015).
#' Joining the dots: an automated method for constructing food webs from
#' compendia of published interactions. *Food Webs*, 5, 11--20.
#' \doi{10.1016/j.fooweb.2015.09.001}
#'
#' @seealso [calibrate_links()], [score_rule()], [compare_rules()].
#' @export
#' @examplesIf interactive()
#' v <- validate_rule(gray2015, sources = c("46", "49", "58"))
#' v
validate_rule <- function(registry, sources, groups = NULL, rule = matching_rule(),
                          train_registries = NULL, exclusions = NULL, quiet = FALSE) {
  check_registry(registry)
  check_rule(rule)
  if (missing(sources) || !is.character(sources) || length(sources) == 0L || anyNA(sources)) {
    abort_input("`sources` must be a character vector of study identifiers (study_id).")
  }
  unknown <- setdiff(sources, registry$study_id)
  if (length(unknown)) {
    abort_input(sprintf("`sources` has studies that are not in `registry`: %s.",
                        format_items(unknown)))
  }
  if (anyDuplicated(sources)) {
    abort_input("`sources` has repeated study identifiers.")
  }
  left_out <- left_out_sets(groups, sources, registry)
  regs <- sort(unique(registry$registry), method = "radix")
  if (is.null(train_registries)) {
    train_registries <- regs
  } else if (!is.character(train_registries) || !length(train_registries) ||
             !all(train_registries %in% regs)) {
    abort_input(sprintf("`train_registries` must be names of registries in `registry` (%s).",
                        format_items(regs, Inf)))
  }
  if (!is.null(exclusions)) {
    check_data_frame(exclusions)
    check_columns(exclusions, c("web", "study_id"), arg = "exclusions")
  }
  check_flag(quiet)
  rule_cv <- rule
  rule_cv$manual_diet <- character(0)
  webs <- vector("list", length(sources))
  counts <- vector("list", length(sources))
  links <- vector("list", length(sources))
  decay <- vector("list", length(sources))
  for (k in seq_along(sources)) {
    w <- sources[k]
    if (!quiet) message(sprintf("Validation web %s (%d/%d)", w, k, length(sources)))
    extra <- if (is.null(exclusions)) character(0) else exclusions$study_id[exclusions$web == w]
    res <- validate_web(registry, w, rule_cv, exclude = left_out[[k]],
                        train_registries = train_registries, extra = extra)
    webs[[k]] <- res$web
    counts[[k]] <- res$counts
    links[[k]] <- res$links
    decay[[k]] <- res$decay
  }
  webs <- do.call(rbind, webs)
  counts <- do.call(rbind, counts)
  links <- do.call(rbind, links)
  rownames(webs) <- NULL
  rownames(counts) <- NULL
  rownames(links) <- NULL
  if (isTRUE(rule$decay)) {
    decay <- do.call(rbind, decay)
    rownames(decay) <- NULL
  } else {
    decay <- NULL
  }
  structure(
    list(counts = counts, webs = webs, links = links, decay = decay,
         settings = list(rule = rule, sources = sources, groups = groups,
                         exclusions = exclusions,
                         registries = sort(train_registries, method = "radix"),
                         n_records = sum(registry$registry %in% train_registries))),
    class = "bwebs_validation"
  )
}

#' @export
print.bwebs_validation <- function(x, ...) {
  w <- x$webs
  r <- x$settings$rule
  cat(sprintf("<bwebs_validation> %d validation webs, %d links inferred, %.0f%% observed\n",
              nrow(w), sum(w$n_predicted),
              100 * sum(w$n_predicted_observed) / max(sum(w$n_predicted), 1)))
  cat(sprintf("  rule: %s; registries: %s\n", rule_label(r),
              paste(x$settings$registries, collapse = ", ")))
  cat(sprintf("  recall (recorded links inferred) per web from %.2f to %.2f\n",
              min(w$recall, na.rm = TRUE), max(w$recall, na.rm = TRUE)))
  if (!is.null(x$decay)) {
    cat(sprintf("  pairs not inferred: %d, with %d recorded links (evidence of the taxonomic decay)\n",
                as.integer(sum(x$decay$pairs)), as.integer(sum(x$decay$observed))))
  }
  if (length(x$settings$groups)) {
    cat("  studies left out together: ",
        paste(vapply(x$settings$groups, paste, character(1), collapse = " + "), collapse = "; "),
        "\n", sep = "")
  }
  invisible(x)
}

# Short description of a rule: mode/base[/climb_max/min_partners][ no downward].
rule_label <- function(r) {
  lab <- switch(r$mode,
    fixed = sprintf("fixed/%s", r$base_level),
    climb = sprintf("climb/%s/%s/%d", r$base_level, r$climb_max, r$min_partners),
    ladder = sprintf("ladder/%s/%d", r$climb_max, r$min_partners)
  )
  if (!length(r$downward_ranks)) lab <- paste(lab, "no-downward")
  lab
}

# Studies left out with each validation web: the web itself and the other
# studies of its group.
left_out_sets <- function(groups, sources, registry) {
  if (is.null(groups)) {
    return(as.list(sources))
  }
  ok <- is.list(groups) && length(groups) > 0L &&
    all(vapply(groups, function(g) is.character(g) && length(g) >= 2L && !anyNA(g), logical(1)))
  if (!ok) {
    abort_input(paste0(
      "`groups` must be NULL or a list of character vectors of study identifiers ",
      "(study_id), each with at least two studies, e.g. list(c(\"13\", \"50\"))."
    ))
  }
  ids <- unlist(groups, use.names = FALSE)
  unknown <- setdiff(ids, registry$study_id)
  if (length(unknown)) {
    abort_input(sprintf("`groups` has studies that are not in `registry`: %s.",
                        format_items(unknown)))
  }
  repeated <- unique(ids[duplicated(ids)])
  if (length(repeated)) {
    abort_input(sprintf("`groups` lists these studies more than once: %s.",
                        format_items(repeated)))
  }
  lapply(sources, function(w) {
    g <- as.character(unlist(groups[vapply(groups, function(x) w %in% x, logical(1))],
                             use.names = FALSE))
    c(w, sort(setdiff(g, w), method = "radix"))
  })
}

# Nodes of a validation web, classified with the taxonomy of its records.
registry_nodes <- function(registry, study, rule) {
  d <- registry[registry$study_id == study, , drop = FALSE]
  side <- function(pre, nm) {
    data.frame(name = d[[nm]], genus = d[[paste0(pre, "genus")]],
               subfamily = d[[paste0(pre, "subfamily")]], family = d[[paste0(pre, "family")]],
               order = d[[paste0(pre, "order")]], class = d[[paste0(pre, "class")]],
               stringsAsFactors = FALSE)
  }
  t <- rbind(side("res_", "resource"), side("con_", "consumer"))
  t <- t[!is.na(t$name) & !duplicated(t$name), , drop = FALSE]
  lname <- tolower(t$name)
  eq <- function(v) !is.na(v) & lname == tolower(v)
  rank <- rep("informal", nrow(t))
  rank[eq(t$class)] <- "class"
  rank[eq(t$order)] <- "order"
  rank[eq(t$family)] <- "family"
  rank[eq(t$subfamily)] <- "subfamily"
  rank[!is.na(t$genus) & startsWith(lname, paste0(tolower(t$genus), " "))] <- "species"
  rank[eq(t$genus)] <- "genus"
  notax <- is.na(t$genus) & is.na(t$subfamily) & is.na(t$family) & is.na(t$order) & is.na(t$class)
  rank[notax] <- "non_taxonomic"
  consumers <- unique(d$consumer)
  cat_res <- registry$res_category[match(t$name, registry$resource)]
  role <- rep("consumer", nrow(t))
  role[!t$name %in% consumers & !is.na(cat_res) & cat_res == "producer"] <- "producer"
  role[!t$name %in% consumers & notax & (is.na(cat_res) | cat_res != "producer")] <- "detritus"
  taxa <- data.frame(
    node_id = t$name, taxon = t$name, rank = rank, role = role,
    node_type = ifelse(notax, "non_taxonomic", "taxon"),
    genus = t$genus, subfamily = t$subfamily, family = t$family, order = t$order,
    class = t$class, stringsAsFactors = FALSE
  )
  nodes <- prepare_nodes(taxa, rule = rule)
  list(nodes = nodes, studied = t$name %in% consumers, records = d)
}

# Leave one web out: rebuild its metaweb from the rest of the registry and
# compare the inferred links with the recorded ones.
validate_web <- function(registry, study, rule, exclude = study, train_registries = NULL,
                         extra = character(0)) {
  vn <- registry_nodes(registry, study, rule)
  nodes <- vn$nodes
  exclude <- unique(c(study, exclude))
  keep <- !registry$study_id %in% c(exclude, extra)
  if (!is.null(train_registries)) keep <- keep & registry$registry %in% train_registries
  train <- registry[keep, , drop = FALSE]
  mw <- metaweb_core(nodes, train, check_manual_links(NULL, nodes$node_id), prefilter = "both")
  e <- mw$edges[!mw$edges$cannibal, , drop = FALSE]
  studied_ids <- nodes$node_id[vn$studied]
  e <- e[e$consumer %in% studied_ids, , drop = FALSE]
  key <- function(r, c) paste(r, c, sep = "\u001f")
  rec <- vn$records
  rec <- rec[rec$resource != rec$consumer, , drop = FALSE]
  truth <- unique(key(rec$resource, rec$consumer))
  e$observed <- key(e$resource, e$consumer) %in% truth
  n_nodes <- nrow(nodes)
  n_pairs <- n_nodes * length(studied_ids) - length(studied_ids)
  web <- data.frame(
    study_id = study, n_nodes = n_nodes, n_consumers = length(studied_ids),
    n_pairs = n_pairs, n_recorded = length(truth), n_predicted = nrow(e),
    n_predicted_observed = sum(e$observed),
    precision = if (nrow(e)) sum(e$observed) / nrow(e) else NA_real_,
    recall = if (length(truth)) sum(e$observed) / length(truth) else NA_real_,
    levels_converged = mw$diagnostics$converged,
    left_out = paste(exclude, collapse = ", "),
    n_excluded = length(setdiff(unique(extra), exclude)),
    stringsAsFactors = FALSE
  )
  links <- data.frame(study_id = rep(study, nrow(e)), e[c("resource", "consumer", "evidence_class",
                                                          link_evidence_columns, "observed")],
                      stringsAsFactors = FALSE)
  rownames(links) <- NULL
  if (nrow(e) == 0L) {
    counts <- data.frame(study_id = character(), evidence_class = character(),
                         predicted = integer(), observed = integer(), stringsAsFactors = FALSE)
  } else {
    cls <- factor(e$evidence_class, levels = evidence_classes)
    pred <- tabulate(cls, length(evidence_classes))
    obs <- tabulate(cls[e$observed], length(evidence_classes))
    keep <- pred > 0
    counts <- data.frame(study_id = rep(study, sum(keep)), evidence_class = evidence_classes[keep],
                         predicted = pred[keep], observed = obs[keep], stringsAsFactors = FALSE)
  }
  # evidence of the taxonomic decay for the pairs without an inferred link
  decay <- NULL
  if (isTRUE(rule$decay)) {
    pairs <- expand.grid(resource = nodes$node_id, consumer = studied_ids, stringsAsFactors = FALSE)
    pairs <- pairs[pairs$resource != pairs$consumer, , drop = FALSE]
    pairs <- pairs[!key(pairs$resource, pairs$consumer) %in% key(e$resource, e$consumer), ,
                   drop = FALSE]
    prof <- decay_profiles(nodes, train, rule, consumers = studied_ids)
    ev <- decay_evidence(prof, pairs$resource, pairs$consumer)
    decay <- decay_count_table(study, ev, key(pairs$resource, pairs$consumer) %in% truth)
    decay$web_nodes <- rep(nrow(nodes), nrow(decay))
  }
  list(web = web, counts = counts, links = links, decay = decay)
}

# Columns of the edges with the amount of evidence of each link (see
# build_metaweb()).
link_evidence_columns <- c("n_records", "n_sources", "n_consumer_taxa", "n_resource_taxa",
                           "n_consumer_known", "n_resource_known")
