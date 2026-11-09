# Matching rule ---------------------------------------------------------------

#' Matching rule for inferring trophic links from a registry
#'
#' `matching_rule()` collects every decision that turns a table of taxa and a
#' registry of published trophic interactions into a metaweb: the taxonomic
#' level at which each node is matched with the registry, when and how far the
#' match climbs to coarser ranks, which records are extended downwards and
#' which registry counts as primary. The same rule object is used by
#' [prepare_nodes()], [build_metaweb()] and [calibrate_links()], so that link
#' probabilities are always calibrated for the rule that produced the links.
#'
#' The defaults are the combination with the lowest cross-validated Brier
#' score among the rules of [candidate_rules()] and the sources compared
#' (the registry of Gray et al. 2015, GloBI and both), with the validation
#' webs of the registry left out by system and Tuesday Lake kept apart as an
#' independent test (see [compare_rules()], [source_comparison] and the
#' vignette): start at family and climb up to class until a node links with
#' at least five other nodes, using the records of every registry at every
#' level. The rule of the reanalysis of Amador et al. (2024) that this package
#' was written from is the same with `min_partners = 1`, whose score was
#' indistinguishable. The taxonomic decay of the pairs not inferred
#' (`decay`), added afterwards, lowered the cross-validated Brier score of
#' the default rule with the Gray registry, left it unchanged with Gray +
#' GloBI and improved the scores on Tuesday Lake. The ladder of taxonomic
#' generalisation follows the WebBuilder method of Gray et al. (2015); the
#' climbing rule, the evidence threshold, the downward inheritance and the
#' generic resource groups are extensions of that method (see Details).
#'
#' @param base_level Starting level of the match for every node and role
#'   (rule R2): one of `"exact"`, `"genus"`, `"subfamily"`, `"family"`,
#'   `"order"`, `"class"`. Nodes identified at a coarser rank start at their
#'   own rank (rule R3).
#' @param climb_max Highest level the match may climb to (rule R4).
#' @param mode How a node collects registry records (see Details):
#'   `"climb"`, `"fixed"` or `"ladder"`.
#' @param min_partners Minimum number of other nodes of the network a node
#'   must link with at a level for the match to stop climbing (integer >= 1).
#'   `1` is the criterion of Gray et al. (2015) and of the reanalysis of
#'   Amador et al. (2024).
#' @param climb_predators If `FALSE`, predators (nodes for which more than
#'   `predator_majority` of their prey are animals) do not climb as
#'   resources, so that top predators without recorded predators in their
#'   family do not inherit the predators of their whole order or class.
#' @param predator_majority Proportion of animal prey above which a node is
#'   treated as a predator (in (0, 1)).
#' @param primary_registry `NULL` (default: every registry is used at every
#'   level and for downward inheritance) or the name of the registry, in the
#'   `registry` column of the registry table, that is used above
#'   `supplementary_max_level` and for downward inheritance, the other
#'   registries being supplementary.
#' @param supplementary_max_level Highest matching level at which records of
#'   registries other than `primary_registry` are used (ignored when
#'   `primary_registry` is `NULL`).
#' @param downward_ranks Ranks from which a record resolved only to that rank
#'   (for example "Gerris eats Diptera") is extended to every node of that
#'   taxon. `character(0)` switches downward inheritance off.
#' @param downward_roles Roles in which downward inheritance applies
#'   (`"resource"`, `"consumer"` or both).
#' @param generic_groups Generic, non-taxonomic resources of the registry
#'   (such as "Algae") that are extended to every node whose rank `rank` is
#'   one of `values`: a data frame with columns `generic`, `rank` and
#'   `values` (values separated by `;`), or `NULL` for none.
#' @param exceptions Exceptions to `base_level` by taxon: a data frame with
#'   columns `role` (`"resource"`, `"consumer"` or `"both"`), `rank` (column
#'   of the node table compared with `taxon`), `taxon`, `base` (starting
#'   level), `rule` and `reason`, or `NULL` for none.
#' @param manual_diet Node identifiers of consumers whose diet comes only
#'   from manual links (they are not matched as consumers).
#' @param max_iter Maximum number of passes used to stabilise the levels
#'   (the number of partners of a node depends on the levels of its
#'   partners).
#' @param decay If `TRUE` (default), the resource → consumer pairs that the
#'   rule does not infer get a probability from the taxonomic distance decay
#'   instead of 0 (see Details): the closer the recorded diets of the
#'   relatives of the consumer come to the resource, the more probable the
#'   link. The decay is fitted by [calibrate_links()] and used by
#'   [sample_metawebs()], [score_rule()] and [score_metaweb()]; the
#'   deterministic metaweb of [build_metaweb()] keeps only the inferred
#'   links. `FALSE` gives those pairs probability 0, as in WebBuilder.
#'
#' @return An object of class `bwebs_rule`: a named list with the arguments
#'   above, validated.
#'
#' @details
#' Rules applied by [prepare_nodes()] to each node and role:
#' * **R1** the records of the node's own taxon (its name, or its rank for
#'   tribes and suborders) are always used;
#' * **R2** starting level `base_level`, in every group and both roles;
#' * **R3** nodes identified at a rank coarser than `base_level` start at
#'   their own rank;
#' * **R4** the match climbs to `climb_max` (mode `"climb"`);
#' * **R5** an explicit trophic surrogate replaces the ladder;
#' * **NT** non-taxonomic nodes match by name and aliases;
#' * basal nodes (producers and detritus) are not matched as consumers, and
#'   `manual_diet` nodes take their diet only from manual links.
#'
#' Modes. `"climb"`: a node takes all the records of its own taxon and of
#' its starting level; if with them it links with fewer than `min_partners`
#' other nodes, it climbs to the next level, accumulating records. `"fixed"`:
#' the same without climbing. `"ladder"`: the WebBuilder ladder of Gray et
#' al. (2015): from the node's own rank, the first level that links the node
#' with at least `min_partners` nodes, without accumulating records. If no
#' level reaches the threshold, the finest level with at least one partner is
#' used, and failing that the finest level with at least one record.
#'
#' Taxonomic distance decay (`decay = TRUE`). For a pair without an inferred
#' link, the consumer is resolved at the closest level, from its family on
#' (its order, class, phylum, or any consumer), whose consumers have at
#' least one resource recorded in the registry (records of the consumer
#' itself, its genus or its subfamily count at the family), and the resource
#' at the closest level (its own taxon, genus, subfamily, family, order,
#' class, phylum or none) at which it is related to a resource of those
#' records; a record of a generic group that applies to the resource also
#' counts. The probability of the pair decreases with both levels and is
#' never 0 (see [calibrate_links()]). Cannibalism, the pairs removed by
#' manual links and the diets of basal and `manual_diet` nodes stay at 0.
#'
#' @references
#' Gray C, Figueroa DH, Hudson LN, Ma A, Perkins D, Woodward G (2015).
#' Joining the dots: an automated method for constructing food webs from
#' compendia of published interactions. *Food Webs*, 5, 11--20.
#' \doi{10.1016/j.fooweb.2015.09.001}
#'
#' @export
#' @examples
#' matching_rule()
#' matching_rule(mode = "ladder", base_level = "genus")
matching_rule <- function(base_level = "family", climb_max = "class", mode = "climb",
                          min_partners = 5L, climb_predators = FALSE,
                          predator_majority = 0.5, primary_registry = NULL,
                          supplementary_max_level = "family",
                          downward_ranks = c("subfamily", "family", "order"),
                          downward_roles = "resource",
                          generic_groups = default_generic_groups(),
                          exceptions = default_exceptions(),
                          manual_diet = character(0), max_iter = 5L, decay = TRUE) {
  check_choice(base_level, match_levels)
  check_choice(climb_max, match_levels)
  check_choice(mode, c("climb", "fixed", "ladder"))
  check_number(min_partners, lower = 1)
  if (min_partners != round(min_partners)) {
    abort_input("`min_partners` must be a whole number.")
  }
  check_flag(climb_predators)
  check_number(predator_majority, lower = 0, upper = 1,
               include_lower = FALSE, include_upper = FALSE)
  if (!is.null(primary_registry)) check_string(primary_registry)
  check_choice(supplementary_max_level, match_levels)
  if (!is.character(downward_ranks) ||
      !all(downward_ranks %in% c("subfamily", "family", "order"))) {
    abort_input(paste(
      "`downward_ranks` must be a character vector with values among",
      "'subfamily', 'family' and 'order' (downward inheritance from classes",
      "is not supported)."
    ))
  }
  if (!is.character(downward_roles) ||
      !all(downward_roles %in% c("resource", "consumer"))) {
    abort_input("`downward_roles` must contain 'resource' and/or 'consumer'.")
  }
  generic_groups <- check_generic_groups(generic_groups)
  exceptions <- check_exceptions(exceptions)
  if (!is.character(manual_diet) || anyNA(manual_diet)) {
    abort_input("`manual_diet` must be a character vector of node identifiers.")
  }
  check_number(max_iter, lower = 1)
  check_flag(decay)
  rule <- list(
    base_level = base_level, climb_max = climb_max, mode = mode,
    min_partners = as.integer(min_partners), climb_predators = climb_predators,
    predator_majority = predator_majority, primary_registry = primary_registry,
    supplementary_max_level = supplementary_max_level,
    downward_ranks = downward_ranks, downward_roles = downward_roles,
    generic_groups = generic_groups, exceptions = exceptions,
    manual_diet = manual_diet, max_iter = as.integer(max_iter), decay = decay
  )
  structure(rule, class = "bwebs_rule")
}

#' @export
print.bwebs_rule <- function(x, ...) {
  cat("<bwebs_rule>\n")
  cat(sprintf("  mode: %s; start at %s, climb up to %s; min_partners = %d\n",
              x$mode, x$base_level, x$climb_max, x$min_partners))
  cat(sprintf("  predators climb as resources: %s (predator if > %s animal prey)\n",
              if (x$climb_predators) "yes" else "no", format(x$predator_majority)))
  if (is.null(x$primary_registry)) {
    cat("  registries: all, at every level\n")
  } else {
    cat(sprintf("  primary registry: %s; other registries up to %s\n",
                x$primary_registry, x$supplementary_max_level))
  }
  cat(sprintf("  downward inheritance: %s (roles: %s)\n",
              if (length(x$downward_ranks)) paste(x$downward_ranks, collapse = ", ") else "none",
              paste(x$downward_roles, collapse = ", ")))
  cat(sprintf("  generic resource groups: %d; exceptions: %d; manual-diet nodes: %d\n",
              nrow(x$generic_groups), nrow(x$exceptions), length(x$manual_diet)))
  cat(sprintf("  pairs not inferred: %s\n",
              if (isTRUE(x$decay)) "probability from the taxonomic distance decay" else
                "probability 0"))
  invisible(x)
}

# Levels of the matching ladder, from the finest to the coarsest.
match_levels <- c("exact", "genus", "subfamily", "family", "order", "class")
# Rank columns of node and registry tables (registry columns carry the prefix
# res_ or con_).
rank_cols <- c("genus", "subfamily", "family", "order", "class")

# Generic resources of the Gray et al. (2015) registry and the nodes they are
# extended to (rule of the original analysis).
default_generic_groups <- function() {
  data.frame(
    generic = c("Algae", "Diatom", "Pennatae", "Centricae", "Cyanobacteria"),
    rank = c("phylum", "class", "class", "class", "phylum"),
    values = c(
      "Chlorophyta;Charophyta;Ochrophyta;Cryptophyta;Myzozoa;Euglenozoa",
      "Bacillariophyceae;Coscinodiscophyceae", "Bacillariophyceae",
      "Coscinodiscophyceae", "Cyanobacteria"
    ),
    stringsAsFactors = FALSE
  )
}

# Exception of the original analysis: water mites (Hydracarina) are an
# informal group; in the registry, Trombidiformes are aquatic mites.
default_exceptions <- function() {
  data.frame(
    role = "both", rank = "order", taxon = "Trombidiformes", base = "order",
    rule = "R3",
    reason = "informal group; in the registry, Trombidiformes are aquatic mites",
    stringsAsFactors = FALSE
  )
}

check_generic_groups <- function(x, call = sys.call(-1)) {
  cols <- c("generic", "rank", "values")
  if (is.null(x)) {
    return(data.frame(generic = character(), rank = character(), values = character(),
                      stringsAsFactors = FALSE))
  }
  check_data_frame(x, arg = "generic_groups", call = call)
  check_columns(x, cols, arg = "generic_groups", call = call)
  x <- as_plain_data_frame(x[cols])
  for (col in cols) {
    if (!is.character(x[[col]]) || anyNA(x[[col]])) {
      abort_input(sprintf("Column '%s' of `generic_groups` must be text without NA.", col),
                  call = call)
    }
  }
  x
}

check_exceptions <- function(x, call = sys.call(-1)) {
  cols <- c("role", "rank", "taxon", "base", "rule", "reason")
  if (is.null(x)) {
    return(data.frame(role = character(), rank = character(), taxon = character(),
                      base = character(), rule = character(), reason = character(),
                      stringsAsFactors = FALSE))
  }
  check_data_frame(x, arg = "exceptions", call = call)
  check_columns(x, cols, arg = "exceptions", call = call)
  x <- as_plain_data_frame(x[cols])
  for (col in cols) {
    if (!is.character(x[[col]])) {
      abort_input(sprintf("Column '%s' of `exceptions` must be text.", col), call = call)
    }
  }
  bad_role <- setdiff(x$role, c("resource", "consumer", "both"))
  if (length(bad_role)) {
    abort_input(sprintf("`exceptions$role` must be 'resource', 'consumer' or 'both', not %s.",
                        format_items(quote_items(bad_role))), call = call)
  }
  bad_base <- setdiff(x$base, match_levels)
  if (length(bad_base)) {
    abort_input(sprintf("`exceptions$base` must be a matching level (%s), not %s.",
                        format_items(match_levels, Inf), format_items(quote_items(bad_base))),
                call = call)
  }
  x
}

# Two rules are the same if every setting is identical (used to check that a
# calibration was fitted with the rule that produced a metaweb).
same_rule <- function(a, b) {
  identical(unclass(a), unclass(b))
}
