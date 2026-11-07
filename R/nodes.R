# Nodes and matching ladder ---------------------------------------------------

#' Prepare the nodes of a metaweb and their matching ladder
#'
#' `prepare_nodes()` validates a table with one row per node (taxa and
#' non-taxonomic resources such as detritus) and adds, for each node and
#' role, the levels at which it is matched with the registry: its own rank
#' (`*_method`), the starting level (`*_base`), the highest level it may
#' climb to (`*_max`), and the rule that set them (`*_rule`, `*_rule_reason`),
#' following [matching_rule()].
#'
#' @param taxa A data frame with one row per node and columns:
#'   * `node_id` (character, unique): identifier used everywhere else;
#'   * `taxon` (character): name matched exactly against the consumer and
#'     resource names of the registry;
#'   * `rank` (character): rank of identification (`"species"`, `"genus"`,
#'     `"tribe"`, `"subfamily"`, `"family"`, `"suborder"`, `"order"`,
#'     `"class"`; any other value, such as `"phylum"` or `"informal"`, matches
#'     by name only);
#'   * `role`: `"consumer"`, `"producer"` or `"detritus"` (producers and
#'     detritus are basal and never consume);
#'   * `node_type`: `"taxon"`, `"non_taxonomic"` or `"functional"`
#'     (non-taxonomic and functional nodes match by name and aliases only);
#'   * `genus`, `subfamily`, `family`, `order`, `class` (character, `NA`
#'     where not resolved): classification aligned with the registry;
#'   * optional `phylum` (used by generic resource groups), `species`,
#'     `alias` (other names, separated by `;`, matched exactly against the
#'     registry names and genera), `life_stage` (`"larva"` or `"adult"`
#'     restrict the records to that stage; other values do not filter),
#'     `resource_surrogate` and `consumer_surrogate` (taxa of any rank whose
#'     records replace the ladder in that role, separated by `;`).
#' @param surrogates Optional data frame with one row per node and role:
#'   columns `node_id`, `role` (`"resource"` or `"consumer"`) and
#'   `surrogate`. It fills `resource_surrogate` and `consumer_surrogate`;
#'   give the surrogates either here or as columns of `taxa`, not both.
#' @param rule A [matching_rule()].
#' @param taxonomy_fixes Named character vector of names replaced in the
#'   classification columns (`genus` to `class`), the same harmonisation that
#'   [as_registry()] applies to the registry; `default_taxonomy_fixes()` by
#'   default, `character(0)` to keep the classification as given.
#'
#' @return A data frame of class `bwebs_nodes` with the rows of `taxa` in the
#'   same order, the columns above (missing optional columns are added with
#'   `NA`; classification after `taxonomy_fixes`) and, for `resource` and
#'   `consumer`, the columns `*_method`,
#'   `*_base`, `*_max`, `*_rule` and `*_rule_reason`. Levels are one of
#'   `"exact"`, `"genus"`, `"subfamily"`, `"family"`, `"order"`, `"class"`,
#'   or `"none"` when the node does not act in that role. The rule is stored
#'   in the attribute `"rule"`.
#'
#' @section Missing values and zeros:
#' A missing rank (`NA`) means the rank was not resolved; it is never filled
#' in, and the node cannot be matched at that level.
#'
#' @section Assumptions:
#' The classification of the nodes uses the same names as the registry at
#' every rank: matching at family level only works if the family is spelt as
#' in the registry. `taxonomy_fixes` aligns the synonyms it lists (such as
#' `"Cladocera"` and `"Diplostraca"`); other differences between
#' classifications (a family placed in another order) are not resolved.
#'
#' @seealso [build_metaweb()] to infer the links.
#' @export
#' @examples
#' nodes <- prepare_nodes(example_taxa)
#' nodes[, c("node_id", "resource_base", "consumer_base", "consumer_rule")]
prepare_nodes <- function(taxa, surrogates = NULL, rule = matching_rule(),
                          taxonomy_fixes = default_taxonomy_fixes()) {
  check_data_frame(taxa, allow_empty = FALSE)
  check_rule(rule)
  if (!is.character(taxonomy_fixes) || (length(taxonomy_fixes) && is.null(names(taxonomy_fixes)))) {
    abort_input("`taxonomy_fixes` must be a named character vector (names: as found; values: replacement).")
  }
  required <- c("node_id", "taxon", "rank", "role", "node_type", rank_cols)
  check_columns(taxa, stats::setNames(required, rep("", length(required))), arg = "taxa")
  nodes <- as_plain_data_frame(taxa)
  optional <- c("phylum", "species", "alias", "life_stage",
                "resource_surrogate", "consumer_surrogate")
  for (col in optional) {
    if (!col %in% names(nodes)) nodes[[col]] <- rep(NA_character_, nrow(nodes))
  }
  text_cols <- c(required, optional)
  for (col in text_cols) {
    v <- nodes[[col]]
    if (is.factor(v)) v <- as.character(v)
    if (is.logical(v) && all(is.na(v))) v <- as.character(v)
    if (!is.character(v)) {
      abort_input(sprintf("Column '%s' of `taxa` must be text, not %s.", col, class(v)[1L]))
    }
    v <- trimws(v)
    v[!is.na(v) & v == ""] <- NA_character_
    nodes[[col]] <- v
  }
  for (col in c("node_id", "taxon", "rank", "role", "node_type")) {
    rows <- which(is.na(nodes[[col]]))
    if (length(rows)) {
      abort_input(sprintf("Column '%s' of `taxa` has missing values in %s.",
                          col, format_rows(rows)))
    }
  }
  check_unique_key(nodes, "node_id", arg = "taxa")
  for (col in rank_cols) nodes[[col]] <- apply_fixes(nodes[[col]], taxonomy_fixes)
  bad <- which(!nodes$role %in% c("consumer", "producer", "detritus"))
  if (length(bad)) {
    abort_input(sprintf(
      "Column 'role' of `taxa` must be 'consumer', 'producer' or 'detritus'; check %s (%s).",
      format_rows(bad), format_items(quote_items(unique(nodes$role[bad])))
    ))
  }
  bad <- which(!nodes$node_type %in% c("taxon", "non_taxonomic", "functional"))
  if (length(bad)) {
    abort_input(sprintf(
      "Column 'node_type' of `taxa` must be 'taxon', 'non_taxonomic' or 'functional'; check %s.",
      format_rows(bad)
    ))
  }
  if (!is.null(surrogates)) {
    if (any(!is.na(nodes$resource_surrogate)) || any(!is.na(nodes$consumer_surrogate))) {
      abort_input(paste(
        "Surrogates were given both in `surrogates` and as columns of `taxa`;",
        "give them in one place only."
      ))
    }
    nodes <- add_surrogates(nodes, surrogates)
  }
  for (role in c("resource", "consumer")) {
    lad <- lapply(seq_len(nrow(nodes)), function(k) ladder_for_node(nodes[k, , drop = FALSE], role, rule))
    for (f in c("method", "base", "max", "rule", "reason")) {
      col <- if (f == "reason") paste0(role, "_rule_reason") else paste0(role, "_", f)
      nodes[[col]] <- vapply(lad, function(x) x[[f]], character(1))
    }
  }
  rownames(nodes) <- NULL
  structure(nodes, rule = rule, class = c("bwebs_nodes", "data.frame"))
}

check_rule <- function(rule, call = sys.call(-1)) {
  if (!inherits(rule, "bwebs_rule")) {
    abort_input("`rule` must be created with matching_rule().", call = call)
  }
  invisible(rule)
}

add_surrogates <- function(nodes, surrogates, call = sys.call(-1)) {
  check_data_frame(surrogates, arg = "surrogates", call = call)
  check_columns(surrogates, c("node_id", "role", "surrogate"), arg = "surrogates", call = call)
  s <- as_plain_data_frame(surrogates[c("node_id", "role", "surrogate")])
  for (col in names(s)) s[[col]] <- trimws(as.character(s[[col]]))
  s <- s[!is.na(s$surrogate) & s$surrogate != "", , drop = FALSE]
  bad <- setdiff(s$role, c("resource", "consumer"))
  if (length(bad)) {
    abort_input(sprintf("`surrogates$role` must be 'resource' or 'consumer', not %s.",
                        format_items(quote_items(bad))), call = call)
  }
  unknown <- setdiff(s$node_id, nodes$node_id)
  if (length(unknown)) {
    abort_input(sprintf("`surrogates` has node_id values that are not in `taxa`: %s.",
                        format_items(unknown)), call = call)
  }
  check_unique_key(s, c("node_id", "role"), arg = "surrogates", call = call)
  for (role in c("resource", "consumer")) {
    sr <- s[s$role == role, , drop = FALSE]
    nodes[[paste0(role, "_surrogate")]][match(sr$node_id, nodes$node_id)] <- sr$surrogate
  }
  nodes
}

# Ranks with their own step in the WebBuilder ladder (tribes and suborders are
# matched at subfamily and order level).
ladder_ranks <- c("species", "genus", "tribe", "subfamily", "family", "suborder",
                  "order", "class")

# Matching level that corresponds to the rank of identification (rule R1).
# Ranks without a step (subclass, phylum, informal groups) and non-taxonomic
# nodes match by name and aliases ("exact").
own_match_level <- function(rank) {
  r <- tolower(rank)
  out <- rep("exact", length(r))
  out[r == "genus"] <- "genus"
  out[r %in% c("tribe", "subfamily")] <- "subfamily"
  out[r == "family"] <- "family"
  out[r %in% c("suborder", "order")] <- "order"
  out[r == "class"] <- "class"
  out
}

split_alias <- function(x) {
  if (length(x) == 0L || is.na(x)) {
    return(character(0))
  }
  a <- trimws(strsplit(x, ";", fixed = TRUE)[[1]])
  a[a != ""]
}

node_surrogates <- function(node, role) {
  split_alias(node[[paste0(role, "_surrogate")]])
}

has_surrogate <- function(node, role) {
  length(node_surrogates(node, role)) > 0L
}

exception_hits <- function(node, role, exceptions) {
  if (nrow(exceptions) == 0L) {
    return(exceptions)
  }
  ok <- vapply(seq_len(nrow(exceptions)), function(e) {
    col <- exceptions$rank[e]
    v <- if (col %in% names(node)) node[[col]] else NA_character_
    exceptions$role[e] %in% c(role, "both") && !is.na(v) &&
      tolower(v) == tolower(exceptions$taxon[e])
  }, logical(1))
  exceptions[ok, , drop = FALSE]
}

# Levels of one node in one role. Order of precedence: the node does not act
# in that role (basal, manual) > explicit surrogate (R5) > non-taxonomic
# node (NT) > starting level (R2, R3) > exception by taxon. Climbing (R4)
# reaches climb_max.
ladder_for_node <- function(node, role, rule) {
  lv <- function(x) match(x, match_levels)
  res <- function(mn, base, mx, rl, why) {
    list(method = mn, base = base, max = mx, rule = rl, reason = why)
  }
  if (role == "consumer" && node$role %in% c("producer", "detritus")) {
    return(res("none", "none", "none", "basal", "basal node: does not consume"))
  }
  if (role == "consumer" && node$node_id %in% rule$manual_diet) {
    return(res("none", "none", "none", "manual", "no diet in the registry: manual links only"))
  }
  if (has_surrogate(node, role)) {
    return(res("exact", "exact", "exact", "R5",
               paste("trophic surrogate:", node[[paste0(role, "_surrogate")]])))
  }
  if (!identical(node$node_type, "taxon")) {
    return(res("exact", "exact", "exact", "NT", "non-taxonomic node: matches by name and alias"))
  }
  lo <- own_match_level(node$rank)
  coarse <- !tolower(node$rank) %in% ladder_ranks || lv(lo) > lv(rule$base_level)
  out <- if (coarse) {
    res(lo, lo, NA_character_, "R3", paste("identified to", node$rank))
  } else {
    res(lo, rule$base_level, NA_character_, "R2", paste("starting level:", rule$base_level))
  }
  ex <- exception_hits(node, role, rule$exceptions)
  if (nrow(ex) > 1L) {
    abort_input(sprintf("More than one exception applies to node '%s' as %s: %s.",
                        node$node_id, role, format_items(ex$taxon)))
  }
  if (nrow(ex) == 1L && lv(ex$base) >= lv(lo)) {
    out$base <- ex$base
    out$rule <- ex$rule
    out$reason <- ex$reason
  }
  out$max <- match_levels[max(lv(out$base), lv(rule$climb_max))]
  out
}

#' @export
print.bwebs_nodes <- function(x, ...) {
  needed <- c("node_id", "rank", "role", "resource_base", "consumer_base", "consumer_max",
              "resource_rule", "consumer_rule")
  if (!all(needed %in% names(x))) {
    # a subset of columns: print it as a plain data frame
    print(as.data.frame(unclass(x), stringsAsFactors = FALSE), ...)
    return(invisible(x))
  }
  cat(sprintf("<bwebs_nodes> %d nodes (%d consumers, %d producers, %d detritus)\n",
              nrow(x), sum(x$role == "consumer"), sum(x$role == "producer"),
              sum(x$role == "detritus")))
  tab <- table(factor(x$consumer_rule, levels = c("R2", "R3", "R5", "NT", "basal", "manual")))
  cat("  rules as consumer: ", paste(names(tab), tab, sep = " = ", collapse = ", "), "\n", sep = "")
  tab <- table(factor(x$resource_rule, levels = c("R2", "R3", "R5", "NT")))
  cat("  rules as resource: ", paste(names(tab), tab, sep = " = ", collapse = ", "), "\n", sep = "")
  cols <- c("node_id", "rank", "role", "resource_base", "consumer_base", "consumer_max")
  df <- as.data.frame(unclass(x), stringsAsFactors = FALSE)[cols]
  n_show <- min(nrow(df), 10L)
  print(df[seq_len(n_show), , drop = FALSE], row.names = FALSE)
  if (nrow(df) > n_show) cat(sprintf("  ... and %d more nodes\n", nrow(df) - n_show))
  invisible(x)
}
