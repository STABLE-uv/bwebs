# Matrix of link probabilities ---------------------------------------------------

#' Link probabilities of a metaweb
#'
#' `link_probabilities()` gives the probability of every resource → consumer
#' pair of a metaweb: the matrix from which [sample_metawebs()] draws its
#' food webs, so that a pair of probability 0.25 is present in about one
#' draw in four. A link inferred by the matching rule takes the calibrated
#' probability of its class of evidence (the taxonomic level of its closest
#' record; see [build_metaweb()] and [calibrate_links()]); with `decay = TRUE`
#' in the rule, a pair without an inferred link takes the probability of its
#' evidence for the taxonomic decay; with `masses`, the probability of every
#' predation pair is lowered by the size term. Manual links have probability
#' 1; cannibalism, the pairs removed by manual links and the diets of basal
#' nodes, 0.
#'
#' @inheritParams sample_metawebs
#' @param system Effect of the system: `"random"` (default) gives the
#'   probability in a new system, averaged over the distribution of system
#'   effects (the probability with which [sample_metawebs()] includes a pair
#'   with `system = "random"`); `"average"` fixes it at 0; a number fixes it
#'   at that value (on the logit scale).
#' @param masses `NULL` (default, no size term) or the body masses of the
#'   nodes, a vector named by `node_id` (see [sample_metawebs()]).
#'
#' @return A numeric matrix resource (rows) x consumer (columns) named by
#'   `node_id`, like the `adjacency` of the metaweb.
#'
#' @details
#' The probabilities are posterior means: with `system = "random"`, the
#' probability of a class (or of an evidence profile of the decay) is
#' \eqn{E[\mathrm{logit}^{-1}(\alpha + u)]}{E[plogis(alpha + u)]} over the
#' posterior draws of the coefficients and \eqn{u \sim N(0, \sigma^2)}{u ~
#' N(0, sigma^2)} (`p_new` of the calibration). The decay is evaluated for a
#' web of as many nodes as the metaweb, as in [sample_metawebs()]. With
#' `masses`, the size term of [sample_metawebs()] shifts the logit of these
#' probabilities, as in [score_metaweb()].
#'
#' In the draws of [sample_metawebs()] the links of a draw share its
#' parameters and system effect, so the number of links varies more among
#' draws than in webs drawn link by link from this matrix; the probability
#' of each pair is the same.
#'
#' @seealso [sample_metawebs()], [calibrate_links()].
#' @export
#' @examples
#' nodes <- prepare_nodes(example_taxa)
#' mw <- build_metaweb(nodes, combine_registries(gray2015, gateway))
#' P <- link_probabilities(mw, gray2015_gateway_calibration)
#' round(P[1:5, 1:5], 2)
#' # without the taxonomic decay of the pairs not inferred
#' nodes0 <- prepare_nodes(example_taxa, rule = matching_rule(decay = FALSE))
#' mw0 <- build_metaweb(nodes0, combine_registries(gray2015, gateway))
#' sum(link_probabilities(mw0, gray2015_gateway_calibration) > 0)
link_probabilities <- function(metaweb, calibration, system = "random", uncalibrated = NULL,
                               masses = NULL, size_model = NULL) {
  if (!inherits(metaweb, "bwebs_metaweb")) {
    abort_input("`metaweb` must be created with build_metaweb().")
  }
  if (!inherits(calibration, "bwebs_calibration")) {
    abort_input("`calibration` must be created with calibrate_links().")
  }
  if (is.character(system)) {
    check_choice(system, c("random", "average"))
  } else {
    check_number(system, arg = "system")
  }
  if (!is.null(uncalibrated)) check_choice(uncalibrated, c("exclude", "keep"))
  if (!rules_compatible(metaweb$settings$rule, calibration$settings$rule)) {
    abort_input(paste(
      "The calibration was fitted with a different matching rule from the one",
      "that built the metaweb; calibrate the links with the same rule",
      "(calibrate_links(..., rule = )) first."
    ))
  }
  warn_registries(metaweb, calibration)
  use_decay <- !is.null(metaweb$decay)
  if (use_decay && is.null(calibration$decay)) {
    abort_input(paste(
      "The metaweb gives the pairs without a link a probability from the taxonomic",
      "decay (matching_rule(decay = TRUE)), but the calibration has no decay model;",
      "calibrate the links with that rule, or build the metaweb with",
      "matching_rule(decay = FALSE)."
    ))
  }
  ids <- metaweb$nodes$node_id
  P <- matrix(0, length(ids), length(ids), dimnames = list(resource = ids, consumer = ids))
  use_size <- !is.null(masses)
  if (use_size) {
    size_model <- check_size_model(size_model)
    log_masses <- check_masses(masses, ids)
    animal <- ids[metaweb$nodes$role == "consumer"]
  }
  e <- metaweb$edges[!metaweb$edges$cannibal, , drop = FALSE]
  if (nrow(e)) {
    st <- link_status(e, calibration, uncalibrated)
    p <- link_class_probability(calibration, st$class_used, system)
    p[st$status == "manual"] <- 1
    p[st$status == "uncalibrated_kept"] <- 1
    p[st$status == "uncalibrated_excluded"] <- 0
    if (use_size) {
      p <- size_shift_probability(p, e$resource, e$consumer, size_model$weights[["rule"]],
                                  size_model, log_masses, animal)
    }
    P[cbind(e$resource, e$consumer)] <- p
  }
  if (use_decay && nrow(metaweb$decay)) {
    cand <- metaweb$decay
    pd <- decay_probability_of(calibration$decay, cand, n_nodes = length(ids), system = system)
    if (use_size) {
      pd <- size_shift_probability(pd, cand$resource, cand$consumer, size_model$weights[["decay"]],
                                   size_model, log_masses, animal)
    }
    P[cbind(cand$resource, cand$consumer)] <- pd
  }
  P
}
