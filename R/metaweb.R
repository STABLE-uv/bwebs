# Deterministic metaweb ------------------------------------------------------

#' Build the metaweb of potential trophic links
#'
#' `build_metaweb()` infers the potential trophic links among the nodes from
#' a registry of published interactions with the WebBuilder logic of Gray et
#' al. (2015) and the extensions set in [matching_rule()]: a link
#' resource → consumer exists if a record of the registry matches the
#' resource in its resource role and the consumer in its consumer role.
#' Manual links are then added or removed. Every link keeps the evidence
#' that supports it, including the class of evidence used by
#' [sample_metawebs()] to assign it a probability.
#'
#' @param nodes Nodes prepared with [prepare_nodes()]; their rule is used.
#' @param registry A registry in the format of [as_registry()], such as
#'   [gray2015].
#' @param manual_links Optional data frame of curated links with columns
#'   `resource`, `consumer` (node identifiers), `action` (`"add"` or
#'   `"remove"`) and, optionally, `source` and `notes`.
#'
#' @return An object of class `bwebs_metaweb`, a list with
#' * `edges`: one row per link, ordered by consumer and resource (C locale),
#'   with columns `resource`, `consumer`, `support` (`"direct"`: a record
#'   matches both nodes at their resolved levels; `"downward"`: only through
#'   downward inheritance; `"manual"`), `registries`, `n_records` and
#'   `n_sources` (records and distinct studies that support the link),
#'   `res_level` and `con_level` (levels at which resource and consumer were
#'   resolved), `evidence` (types of evidence of the records), `source_ids`
#'   (up to 20 studies), `examples` (up to 5 records), `citations` (up to 3),
#'   `evidence_class` (see Details), the amount of evidence of the link
#'   (`n_consumer_taxa` and `n_resource_taxa`: distinct consumer and resource
#'   taxa among its records, that is the consumer and its relatives that eat
#'   the resource or its relatives, and those relatives of the resource;
#'   `n_consumer_known`: distinct consumer taxa among the records of the
#'   consumer node whose resource is a node of the network, its relatives
#'   with a recorded diet in the network at the level at which it was
#'   resolved; `n_resource_known`: distinct resource taxa among the records
#'   of the resource node whose consumer is a node of the network),
#'   `cannibal` (logical) and `extrapolated` (logical);
#' * `nodes`: the nodes with the levels used in each role
#'   (`resource_level_used`, `consumer_level_used`) and whether the evidence
#'   threshold was met;
#' * `adjacency`: an integer matrix resource (rows) x consumer (columns)
#'   with 1 for a link;
#' * `decay` (with `decay = TRUE` in the rule, otherwise `NULL`): the pairs
#'   without a link that [sample_metawebs()] draws with a probability from
#'   the taxonomic decay, with their evidence (`consumer_level`,
#'   `resource_level`, `generic`; see [matching_rule()]): every pair whose
#'   consumer is a consumer node without a manual diet, except cannibalism,
#'   the links and the pairs removed by `manual_links`;
#' * `settings`: the rule and the registries used;
#' * `diagnostics`: nodes left without links (`unlinked`), and whether the
#'   levels converged.
#'
#' @details
#' The class of evidence of a link is the coarsest taxonomic generalisation
#' needed by its best supporting record. For each record, each side gets the
#' finest rank shared by the node and the taxon of the record (`"exact"` when
#' the record names the node or one of its aliases, then `"genus"`,
#' `"subfamily"`, `"family"`, `"order"`, `"class"`); a record of a generic
#' resource group (such as "Algae") gets `"generic"`, and a record that
#' shares no rank up to class with the node (only through surrogates) gets
#' `"beyond"`. The record's class is the coarser of its two sides and the
#' link takes the finest class among its records. Links that rely on a
#' surrogate sharing no rank with the node are marked `extrapolated` and
#' treated as `"class"` by [sample_metawebs()]. Manual links have class
#' `"manual"`.
#'
#' @section Missing values and zeros:
#' A pair of nodes without any matching record has no link in the metaweb:
#' as in WebBuilder, absence of a record is treated as absence of the link.
#' With `decay = TRUE` in the rule, [sample_metawebs()] still draws such
#' pairs with a small probability that depends on how close the records of
#' the relatives of the consumer come to the resource (`decay`).
#'
#' @section Assumptions:
#' The registry records are correct and their taxonomy uses the same names as
#' the node classification.
#'
#' @references
#' Gray C, Figueroa DH, Hudson LN, Ma A, Perkins D, Woodward G (2015).
#' Joining the dots: an automated method for constructing food webs from
#' compendia of published interactions. *Food Webs*, 5, 11--20.
#' \doi{10.1016/j.fooweb.2015.09.001}
#'
#' @seealso [sample_metawebs()] for food webs sampled from link
#'   probabilities; [calc_indicators()].
#' @export
#' @examples
#' nodes <- prepare_nodes(example_taxa)
#' mw <- build_metaweb(nodes, combine_registries(gray2015, gateway))
#' mw
#' head(mw$edges[, c("resource", "consumer", "support", "evidence_class")])
build_metaweb <- function(nodes, registry, manual_links = NULL) {
  if (!inherits(nodes, "bwebs_nodes")) {
    abort_input("`nodes` must be prepared with prepare_nodes().")
  }
  check_registry(registry)
  ml <- check_manual_links(manual_links, nodes$node_id)
  mw <- metaweb_core(nodes, registry, ml, decay = isTRUE(attr(nodes, "rule")$decay))
  if (!mw$diagnostics$converged) {
    warning(sprintf("The matching levels did not stabilise in %d passes; see `diagnostics`.",
                    attr(nodes, "rule")$max_iter), call. = FALSE)
  }
  mw
}

# Metaweb from validated input (used by build_metaweb() and calibrate_links()).
# With decay = TRUE, also the pairs without a link that take a probability
# from the taxonomic decay, with their evidence (from the whole registry).
metaweb_core <- function(nodes, registry, ml, prefilter = "any", decay = FALSE) {
  rule <- attr(nodes, "rule")
  n_records <- nrow(registry)
  regs <- sort(unique(registry$registry), method = "radix")
  n_studies <- length(unique(registry$study_id))
  full <- registry
  registry <- registry[relevant_rows(registry, nodes, rule, sides = prefilter), , drop = FALSE]
  rx <- registry_index(registry)
  m <- compute_matching(nodes, rx, rule)
  edges <- build_edges(nodes, registry, rx, m, rule)
  basal_ids <- nodes$node_id[nodes$role %in% c("producer", "detritus")]
  edges <- curate_edges(edges, ml, basal_ids)
  adj <- edges_to_adjacency(edges, nodes$node_id)
  nodes_out <- as_plain_data_frame(unclass(nodes))
  nodes_out$resource_level_used <- m$res_level
  nodes_out$consumer_level_used <- m$con_level
  nodes_out$resource_threshold_met <- m$res_met
  nodes_out$consumer_threshold_met <- m$con_met
  nodes_out$predator <- m$predator
  nodes_out$n_res_direct <- unname(rowSums(m$Rd))
  nodes_out$n_res_downward <- unname(rowSums(m$Rw))
  nodes_out$n_con_direct <- unname(rowSums(m$Cd))
  nodes_out$n_con_downward <- unname(rowSums(m$Cw))
  candidates <- if (decay) decay_candidates(nodes, full, edges, ml, rule) else NULL
  structure(
    list(
      edges = edges, nodes = nodes_out, adjacency = adj, decay = candidates,
      settings = list(rule = rule, registries = regs, n_records = n_records,
                      n_studies = n_studies),
      diagnostics = list(unlinked = unlinked_nodes(nodes_out, adj, rule$manual_diet),
                         converged = m$converged)
    ),
    class = "bwebs_metaweb"
  )
}

check_manual_links <- function(manual_links, node_ids, call = sys.call(-1)) {
  empty <- data.frame(resource = character(), consumer = character(), action = character(),
                      source = character(), notes = character(), stringsAsFactors = FALSE)
  if (is.null(manual_links)) {
    return(empty)
  }
  check_data_frame(manual_links, arg = "manual_links", call = call)
  check_columns(manual_links, c("resource", "consumer", "action"), arg = "manual_links", call = call)
  ml <- as_plain_data_frame(manual_links)
  for (col in c("source", "notes")) if (!col %in% names(ml)) ml[[col]] <- rep(NA_character_, nrow(ml))
  for (col in c("resource", "consumer", "action", "source", "notes")) {
    ml[[col]] <- trimws(as.character(ml[[col]]))
  }
  bad <- setdiff(c(ml$resource, ml$consumer), node_ids)
  if (length(bad)) {
    abort_input(sprintf("`manual_links` has node identifiers that are not in `nodes`: %s.",
                        format_items(bad)), call = call)
  }
  bad <- setdiff(ml$action, c("add", "remove"))
  if (length(bad)) {
    abort_input(sprintf("`manual_links$action` must be 'add' or 'remove', not %s.",
                        format_items(quote_items(bad))), call = call)
  }
  check_unique_key(ml, c("resource", "consumer"), arg = "manual_links", call = call)
  ml[c("resource", "consumer", "action", "source", "notes")]
}

evidence_classes <- c("exact", "genus", "subfamily", "family", "order", "class",
                      "generic", "beyond")

# Rank shared between a node and the taxon of each record (one side).
side_rank <- function(node, rx, rows, role, generic_names) {
  pre <- role_prefix(role)
  out <- rep(8L, length(rows))
  out[rx[[role]][rows] %in% generic_names] <- 7L
  for (k in rev(seq_along(rank_cols))) {
    v <- node[[rank_cols[k]]]
    if (is.na(v)) next
    out[rx[[paste0(pre, rank_cols[k])]][rows] %in% tolower(v)] <- k + 1L
  }
  names_ok <- tolower(c(node$taxon, split_alias(node$alias)))
  out[rx[[role]][rows] %in% names_ok] <- 1L
  out
}

# Links from a matching. support = "direct" if a direct-direct pair of
# records supports it; otherwise "downward".
build_edges <- function(nodes, reg, rx, m, rule) {
  R_ <- m$Rd
  C_ <- m$Cd
  n_nodes <- nrow(R_)
  # Records of each resource node (direct, and direct or downward): the
  # records of a link are those of its resource that its consumer also
  # matches (in registry order)
  rows_direct <- lapply(seq_len(n_nodes), function(i) which(R_[i, ]))
  rows_all <- lapply(seq_len(n_nodes), function(i) which(R_[i, ] | m$Rw[i, ]))
  # Pairs with shared records, and whether direct records alone link them
  pairs <- lapply(seq_len(n_nodes), function(i) {
    ri <- rows_all[[i]]
    if (!length(ri)) return(NULL)
    j <- which(rowSums(C_[, ri, drop = FALSE] | m$Cw[, ri, drop = FALSE]) > 0)
    if (!length(j)) return(NULL)
    rd <- rows_direct[[i]]
    d <- if (length(rd)) rowSums(C_[j, rd, drop = FALSE]) > 0 else rep(FALSE, length(j))
    cbind(i, j, d)
  })
  idx <- do.call(rbind, pairs)
  n <- if (is.null(idx)) 0L else nrow(idx)
  if (n == 0L) {
    return(empty_edges())
  }
  idx <- idx[order(idx[, 2L], idx[, 1L]), , drop = FALSE]   # as which(arr.ind = TRUE)
  generic_names <- tolower(rule$generic_groups$generic)
  node_rows <- lapply(seq_len(nrow(nodes)), function(k) nodes[k, , drop = FALSE])
  # Relatives of each node with records among the nodes: distinct consumer
  # taxa among the records of a consumer node whose resource is matched by a
  # node (relatives with a recorded diet in the network), and distinct
  # resource taxa among the records of a resource node whose consumer is
  # matched by a node. Records matched on one side only do not count, so
  # that these numbers do not depend on the records left out before the
  # matching (relevant_rows()).
  res_any <- colSums(R_ | m$Rw) > 0
  con_any <- colSums(C_ | m$Cw) > 0
  n_distinct <- function(M, names, other) {
    vapply(seq_len(n_nodes), function(k) length(unique(names[M[k, ] & other])), integer(1))
  }
  known_con <- n_distinct(C_ | m$Cw, rx$consumer_std, res_any)
  known_res <- n_distinct(R_ | m$Rw, rx$resource_std, con_any)
  chr <- function() character(n)
  support <- registries <- evidence <- source_ids <- examples <- citations <- chr()
  evidence_class <- chr()
  n_records <- n_sources <- n_con_taxa <- n_res_taxa <- integer(n)
  for (e in seq_len(n)) {
    i <- idx[e, 1L]
    j <- idx[e, 2L]
    is_direct <- idx[e, 3L] == 1
    rows <- if (is_direct) {
      ri <- rows_direct[[i]]
      ri[C_[j, ri]]
    } else {
      ri <- rows_all[[i]]
      ri[C_[j, ri] | m$Cw[j, ri]]
    }
    sr <- side_rank(node_rows[[i]], rx, rows, "resource", generic_names)
    sc <- side_rank(node_rows[[j]], rx, rows, "consumer", generic_names)
    evidence_class[e] <- evidence_classes[min(pmax(sr, sc))]
    support[e] <- if (is_direct) "direct" else "downward"
    registries[e] <- paste(sort(unique(reg$registry[rows]), method = "radix"), collapse = ";")
    n_records[e] <- length(rows)
    n_sources[e] <- length(unique(reg$study_id[rows]))
    n_con_taxa[e] <- length(unique(rx$consumer_std[rows]))
    n_res_taxa[e] <- length(unique(rx$resource_std[rows]))
    evidence[e] <- paste(sort(unique(stats::na.omit(reg$evidence[rows])), method = "radix"),
                         collapse = " | ")
    source_ids[e] <- paste(utils::head(sort(unique(stats::na.omit(reg$study_id[rows])),
                                            method = "radix"), 20L), collapse = ";")
    examples[e] <- paste(utils::head(unique(paste(reg$consumer[rows], "->", reg$resource[rows])), 5L),
                         collapse = " | ")
    citations[e] <- paste(utils::head(unique(stats::na.omit(reg$citation[rows])), 3L),
                          collapse = " || ")
  }
  data.frame(
    resource = nodes$node_id[idx[, 1L]], consumer = nodes$node_id[idx[, 2L]],
    support = support, registries = registries, n_records = n_records, n_sources = n_sources,
    res_level = m$res_level[idx[, 1L]], con_level = m$con_level[idx[, 2L]],
    evidence = evidence, source_ids = source_ids, examples = examples, citations = citations,
    evidence_class = evidence_class,
    n_consumer_taxa = n_con_taxa, n_resource_taxa = n_res_taxa,
    n_consumer_known = known_con[idx[, 2L]], n_resource_known = known_res[idx[, 1L]],
    stringsAsFactors = FALSE
  )
}

empty_edges <- function() {
  data.frame(resource = character(), consumer = character(), support = character(),
             registries = character(), n_records = integer(), n_sources = integer(),
             res_level = character(), con_level = character(), evidence = character(),
             source_ids = character(), examples = character(), citations = character(),
             evidence_class = character(), n_consumer_taxa = integer(),
             n_resource_taxa = integer(), n_consumer_known = integer(),
             n_resource_known = integer(), stringsAsFactors = FALSE)
}

# Curation: remove/add manual links and drop basal nodes as consumers.
curate_edges <- function(edges, ml, basal_ids) {
  key <- function(r, c) paste(r, c, sep = "\u001f")
  rm_key <- key(ml$resource[ml$action == "remove"], ml$consumer[ml$action == "remove"])
  edges <- edges[!key(edges$resource, edges$consumer) %in% rm_key, , drop = FALSE]
  add <- ml[ml$action == "add", , drop = FALSE]
  add <- add[!key(add$resource, add$consumer) %in% key(edges$resource, edges$consumer), , drop = FALSE]
  if (nrow(add)) {
    na_chr <- rep(NA_character_, nrow(add))
    na_int <- rep(NA_integer_, nrow(add))
    edges <- rbind(edges, data.frame(
      resource = add$resource, consumer = add$consumer, support = rep("manual", nrow(add)),
      registries = rep("manual", nrow(add)), n_records = na_int, n_sources = na_int,
      res_level = na_chr, con_level = na_chr, evidence = add$notes, source_ids = add$source,
      examples = na_chr, citations = na_chr, evidence_class = rep("manual", nrow(add)),
      n_consumer_taxa = na_int, n_resource_taxa = na_int, n_consumer_known = na_int,
      n_resource_known = na_int, stringsAsFactors = FALSE
    ))
  }
  edges <- edges[!edges$consumer %in% basal_ids, , drop = FALSE]
  edges$cannibal <- edges$resource == edges$consumer
  edges$extrapolated <- edges$evidence_class == "beyond"
  edges <- edges[order(edges$consumer, edges$resource, method = "radix"), , drop = FALSE]
  rownames(edges) <- NULL
  edges
}

# Pairs without a link whose probability comes from the taxonomic decay:
# every resource -> consumer pair whose consumer is a consumer node without a
# manual diet, except cannibalism, the links of the metaweb and the pairs
# removed by manual links; with the evidence of decay_evidence().
decay_candidates <- function(nodes, registry, edges, ml, rule) {
  empty <- data.frame(resource = character(), consumer = character(),
                      consumer_level = character(), resource_level = character(),
                      generic = logical(), stringsAsFactors = FALSE)
  consumers <- nodes$node_id[nodes$role == "consumer" & !nodes$node_id %in% rule$manual_diet]
  if (!length(consumers)) return(empty)
  pairs <- expand.grid(resource = nodes$node_id, consumer = consumers, stringsAsFactors = FALSE)
  pairs <- pairs[pairs$resource != pairs$consumer, , drop = FALSE]
  key <- function(r, c) paste(r, c, sep = "\u001f")
  removed <- ml[ml$action == "remove", , drop = FALSE]
  taken <- c(key(edges$resource, edges$consumer), key(removed$resource, removed$consumer))
  pairs <- pairs[!key(pairs$resource, pairs$consumer) %in% taken, , drop = FALSE]
  if (!nrow(pairs)) return(empty)
  prof <- decay_profiles(nodes, registry, rule, consumers = unique(pairs$consumer))
  out <- cbind(pairs, decay_evidence(prof, pairs$resource, pairs$consumer))
  out <- out[order(out$consumer, out$resource, method = "radix"), , drop = FALSE]
  rownames(out) <- NULL
  out
}

edges_to_adjacency <- function(edges, node_ids) {
  a <- matrix(0L, length(node_ids), length(node_ids),
              dimnames = list(resource = node_ids, consumer = node_ids))
  if (nrow(edges)) a[cbind(edges$resource, edges$consumer)] <- 1L
  a
}

# Nodes that need a surrogate or manual links: consumers without prey (they
# would count as basal) and basal nodes without consumers (isolated).
# Consumers with prey but no predators (possible top predators) are reported
# as "none_top" for information.
unlinked_nodes <- function(nodes, adj, manual) {
  A <- adj
  diag(A) <- 0L
  basal <- nodes$role %in% c("producer", "detritus")
  n_prey <- unname(colSums(A)[nodes$node_id])
  n_pred <- unname(rowSums(A)[nodes$node_id])
  need <- rep(NA_character_, nrow(nodes))
  need[!basal & n_pred == 0] <- "none_top"
  need[basal & n_pred == 0] <- "resource"
  need[!basal & n_prey == 0] <- "consumer"
  need[!basal & n_prey == 0 & nodes$node_id %in% manual] <- "manual_links"
  out <- data.frame(node_id = nodes$node_id, taxon = nodes$taxon, node_role = nodes$role,
                    n_prey = n_prey, n_predators = n_pred, need = need,
                    stringsAsFactors = FALSE)
  out <- out[!is.na(out$need), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' @export
print.bwebs_metaweb <- function(x, ...) {
  e <- x$edges
  cat(sprintf("<bwebs_metaweb> %d nodes, %d links (%d cannibalistic, %d manual)\n",
              nrow(x$nodes), nrow(e), sum(e$cannibal), sum(e$support == "manual")))
  cat(sprintf("  registries: %s (%d records, %d studies)\n",
              paste(x$settings$registries, collapse = ", "), x$settings$n_records,
              x$settings$n_studies))
  cls <- table(factor(e$evidence_class[!e$cannibal], levels = c(evidence_classes, "manual")))
  cls <- cls[cls > 0]
  cat("  evidence of non-cannibalistic links: ",
      paste(names(cls), cls, sep = " = ", collapse = ", "), "\n", sep = "")
  if (!is.null(x$decay)) {
    cat(sprintf("  pairs without a link, with a probability from the taxonomic decay: %d\n",
                nrow(x$decay)))
  }
  u <- x$diagnostics$unlinked
  need <- u[u$need %in% c("consumer", "resource", "manual_links"), , drop = FALSE]
  if (nrow(need)) {
    cat(sprintf("  nodes without links: %s\n", format_items(need$node_id, 10L)))
  }
  invisible(x)
}
