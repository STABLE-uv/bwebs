# Food webs sampled from link probabilities ---------------------------------

#' Sample metawebs from calibrated link probabilities
#'
#' `sample_metawebs()` turns a deterministic metaweb into a set of plausible
#' metawebs. In each draw, the probability of every class of evidence is
#' drawn from the calibration (coefficients and, optionally, the effect of
#' the system), and each inferred link is kept with the probability of its
#' class. With `decay = TRUE` in the rule, every pair without a link of the
#' metaweb is drawn too, with the probability of its evidence for the
#' taxonomic distance decay (see [matching_rule()]). The same realisation is
#' later shared by all the samples of a community in [calc_indicators()], so
#' that differences between samples are computed within each draw.
#'
#' @param metaweb A metaweb from [build_metaweb()].
#' @param calibration A calibration from [calibrate_links()] fitted with the
#'   same matching rule and registries, such as [gray2015_gateway_calibration]
#'   for metawebs built with `combine_registries(gray2015, gateway)` or
#'   [gray2015_calibration] for metawebs built with [gray2015] (a warning is
#'   given when the registries differ).
#' @param n_draws Number of draws (positive whole number).
#' @param system Effect of the system on the logit scale: `"random"` draws it
#'   from the distribution of web effects in each draw (an unknown system,
#'   the more conservative choice, with wide intervals); `"average"` fixes it
#'   at 0 (a system like the average validation web); a number fixes it at
#'   that value. The decay model has its own system effect, drawn
#'   independently from its own distribution, set to 0 or to the same
#'   number.
#' @param parameter_uncertainty If `TRUE`, the class coefficients and the
#'   standard deviation of the web effect of each draw are taken from the
#'   posterior draws of the calibration (for a calibration fitted with
#'   `method = "ml"`, from their approximate normal sampling distribution); if
#'   `FALSE`, the posterior means (or the estimates) are used.
#' @param uncalibrated What to do with links without a calibrated
#'   probability, because their class of evidence has no coefficient or
#'   because they are supported only by registries that were not used in the
#'   calibration (such as GloBI with [gray2015_calibration]): `NULL`
#'   (default) stops and lists them; `"exclude"` gives them probability 0;
#'   `"keep"` probability 1.
#' @param masses `NULL` (default, no size term) or the body masses of the
#'   nodes (mean individual mass, in any unit common to all nodes): a vector
#'   named by `node_id`, the same in every draw, or a matrix with one row per
#'   node (row names `node_id`) and one column per draw, such as the masses
#'   of the draws of a biomass model, so that each draw uses its own masses.
#'   Nodes without a mass (absent or `NA`) get no size term.
#' @param size_model The effect of the predator/prey mass ratio used with
#'   `masses`; `NULL` (default) uses [size_model].
#' @param seed `NULL` (default) to use the current random number stream, or a
#'   number. When given, the result is reproducible and the global random
#'   number stream is restored afterwards.
#'
#' @return An object of class `bwebs_draws`, a list with
#' * `links`: one row per non-cannibalistic link of the metaweb with
#'   `resource`, `consumer`, `evidence_class`, `class_used` (class whose
#'   probability applies; `"class"` for extrapolated surrogate links, `NA`
#'   for manual links), `status` (`"calibrated"`, `"manual"`,
#'   `"uncalibrated_excluded"`, `"uncalibrated_kept"`), `p_mean` (mean
#'   probability over draws) and `inclusion` (proportion of draws in which
#'   the link is present);
#' * `z`: logical matrix links x draws, `TRUE` when the link is present;
#' * `p`: matrix classes x draws of the probabilities used;
#' * `system_effect`: the effect of the system in each draw;
#' * `background` (with `decay = TRUE` in the rule, otherwise `NULL`): the
#'   pairs without a link drawn from the taxonomic decay, a list with `pairs`
#'   (the `decay` table of the metaweb with `p_mean` and `inclusion`),
#'   `present` (one integer vector per draw with the rows of `pairs` present
#'   in it), `p` (probability of each evidence profile in each draw) and
#'   `system_effect` (of the decay model);
#' * `masses` (with `masses`, otherwise `NULL`): the body masses used, a
#'   matrix nodes x draws (`NA` for nodes without a mass), so that draw
#'   \eqn{d} holds the links of a web (column \eqn{d} of `z` and element
#'   \eqn{d} of `background$present`) and the masses of its nodes (column
#'   \eqn{d} of `masses`);
#' * `nodes`: the nodes of the metaweb (identity, role and classification);
#' * `settings`: the options used (with `masses`, the size model).
#'
#' With `masses`, `p_mean` of the links and of the pairs of the decay
#' includes the size term; `p` does not.
#'
#' @details
#' In draw \eqn{d}, with class coefficients \eqn{\alpha_c^{(d)}} and system
#' effect \eqn{u^{(d)}}, link \eqn{l} of class \eqn{c(l)} is present with
#' probability
#' \deqn{p_{c(l)}^{(d)} = \mathrm{logit}^{-1}\left(\alpha_{c(l)}^{(d)} + u^{(d)}\right)}{
#' p = plogis(alpha[c(l)] + u)}
#' independently of the other links. Manual links are always present;
#' cannibalistic links are left out, as in the indicators. With the decay, a
#' pair without a link whose consumer is resolved at level \eqn{c}, whose
#' resource is related at level \eqn{r} and with \eqn{g = 1} when a generic
#' record applies is present with probability
#' \deqn{\mathrm{logit}^{-1}\left(a_c^{(d)} + b_r^{(d)} + \gamma^{(d)} g +
#' \beta^{(d)} \log(n / 50) + v^{(d)}\right)}{plogis(a[c] + b[r] + gamma g +
#' beta log(n / 50) + v)}
#' (see [calibrate_links()]), where \eqn{n} is the number of nodes of the
#' metaweb, with its own system effect \eqn{v^{(d)}}.
#'
#' Body masses. With `masses`, the logit of every predation pair (an animal
#' resource, role `"consumer"`) whose two nodes have a mass in draw \eqn{d}
#' is shifted by \eqn{w\, b\,(x - x^*)^2 \le 0}{w b (x - x*)^2 <= 0}, where
#' \eqn{x = \log_{10}(M_{consumer} / M_{resource})}{x = log10(M_consumer /
#' M_resource)} with the masses of that draw, \eqn{x^*} is the optimum ratio
#' and \eqn{b < 0} the effect of [size_model], and \eqn{w} is the weight of
#' the links inferred by the rule or of the pairs of the decay. Links and
#' pairs of other kinds, manual links and nodes without a mass are not
#' affected. The same random numbers are used with and without `masses`, so
#' that a draw changes only where the size term changes a probability.
#'
#' @section Reproducibility:
#' The same seed can give different draws in other versions of R or with
#' other BLAS/LAPACK libraries.
#'
#' @seealso [calibrate_links()], [calc_indicators()].
#' @export
#' @examples
#' nodes <- prepare_nodes(example_taxa)
#' mw <- build_metaweb(nodes, combine_registries(gray2015, gateway))
#' draws <- sample_metawebs(mw, gray2015_gateway_calibration, n_draws = 20, seed = 1)
#' draws
#' # with body masses (here invented, one per node)
#' masses <- stats::setNames(10^stats::runif(nrow(nodes), -4, 1), nodes$node_id)
#' sample_metawebs(mw, gray2015_gateway_calibration, n_draws = 20, masses = masses, seed = 1)
sample_metawebs <- function(metaweb, calibration, n_draws, system = "random",
                            parameter_uncertainty = TRUE, uncalibrated = NULL, masses = NULL,
                            size_model = NULL, seed = NULL) {
  if (!inherits(metaweb, "bwebs_metaweb")) {
    abort_input("`metaweb` must be created with build_metaweb().")
  }
  if (!inherits(calibration, "bwebs_calibration")) {
    abort_input("`calibration` must be created with calibrate_links().")
  }
  check_number(n_draws, lower = 1)
  if (n_draws != round(n_draws)) abort_input("`n_draws` must be a whole number.")
  if (is.character(system)) {
    check_choice(system, c("random", "average"))
  } else {
    check_number(system, arg = "system")
  }
  check_flag(parameter_uncertainty)
  if (!is.null(uncalibrated)) check_choice(uncalibrated, c("exclude", "keep"))
  check_number(seed, allow_null = TRUE)
  if (!rules_compatible(metaweb$settings$rule, calibration$settings$rule)) {
    abort_input(paste(
      "The calibration was fitted with a different matching rule from the one",
      "that built the metaweb; calibrate the links with the same rule",
      "(calibrate_links(..., rule = )) before sampling."
    ))
  }
  warn_registries(metaweb, calibration)
  if (parameter_uncertainty && is.null(calibration$draws) &&
      !isTRUE(calibration$diagnostics$information_pd)) {
    abort_input(paste(
      "The information matrix of the calibration is not positive definite, so",
      "parameter uncertainty cannot be sampled; use parameter_uncertainty = FALSE."
    ))
  }
  use_decay <- !is.null(metaweb$decay)
  if (use_decay && is.null(calibration$decay)) {
    abort_input(paste(
      "The metaweb gives the pairs without a link a probability from the taxonomic",
      "decay (matching_rule(decay = TRUE)), but the calibration has no decay model;",
      "calibrate the links with that rule, or build the metaweb with",
      "matching_rule(decay = FALSE)."
    ))
  }
  if (use_decay && parameter_uncertainty && is.null(calibration$decay$draws) &&
      !isTRUE(calibration$decay$diagnostics$information_pd)) {
    abort_input(paste(
      "The information matrix of the decay model of the calibration is not positive",
      "definite, so parameter uncertainty cannot be sampled; use",
      "parameter_uncertainty = FALSE or a Bayesian calibration."
    ))
  }
  n_draws <- as.integer(n_draws)
  use_size <- !is.null(masses)
  if (use_size) {
    size_model <- check_size_model(size_model)
    log_masses <- check_masses(masses, metaweb$nodes$node_id, n_draws)
    animal <- metaweb$nodes$node_id[metaweb$nodes$role == "consumer"]
  }
  e <- metaweb$edges[!metaweb$edges$cannibal, , drop = FALSE]
  classes <- calibration$coefficients$evidence_class
  st <- link_status(e, calibration, uncalibrated)
  class_used <- st$class_used
  status <- st$status
  theta_hat <- calibration$theta
  C <- length(classes)
  draw <- function() {
    theta <- if (parameter_uncertainty && !is.null(calibration$draws)) {
      idx <- sample.int(nrow(calibration$draws), n_draws, replace = TRUE)
      t(calibration$draws[idx, , drop = FALSE])
    } else if (parameter_uncertainty) {
      R <- chol(calibration$vcov)
      z <- matrix(stats::rnorm(length(theta_hat) * n_draws), length(theta_hat), n_draws)
      theta_hat + t(R) %*% z
    } else {
      matrix(theta_hat, length(theta_hat), n_draws)
    }
    sigma <- exp(theta[C + 1L, ])
    u <- if (identical(system, "random")) {
      stats::rnorm(n_draws, 0, sigma)
    } else if (identical(system, "average")) {
      rep(0, n_draws)
    } else {
      rep(system, n_draws)
    }
    L <- theta[seq_len(C), , drop = FALSE] + matrix(u, C, n_draws, byrow = TRUE)
    P <- stats::plogis(L)
    dimnames(P) <- list(classes, NULL)
    idx <- match(class_used, classes)
    Z <- matrix(FALSE, nrow(e), n_draws)
    p_sum <- numeric(nrow(e))
    if (use_size) si <- size_index(e$resource, e$consumer, animal, log_masses)
    for (d in seq_len(n_draws)) {
      l <- L[idx, d]
      if (use_size) {
        l <- l + size_shift(size_model, size_model$weights[["rule"]], si, log_masses, d)
      }
      p_l <- stats::plogis(l)
      Z[, d] <- stats::runif(nrow(e)) < p_l
      p_sum <- p_sum + p_l
    }
    bg <- if (use_decay) {
      draw_background(metaweb$decay, calibration$decay, n_draws, system, parameter_uncertainty,
                      n_nodes = nrow(metaweb$nodes),
                      size = if (use_size) list(model = size_model, log_masses = log_masses,
                                                animal = animal))
    }
    list(P = P, u = u, Z = Z, p_mean = p_sum / n_draws, bg = bg)
  }
  res <- if (is.null(seed)) draw() else withr::with_seed(seed, draw())
  Z <- res$Z
  Z[status == "manual", ] <- TRUE
  Z[status == "uncalibrated_kept", ] <- TRUE
  Z[status == "uncalibrated_excluded", ] <- FALSE
  p_mean <- rep(NA_real_, nrow(e))
  ok <- status == "calibrated"
  p_mean[ok] <- res$p_mean[ok]
  p_mean[status %in% c("manual", "uncalibrated_kept")] <- 1
  p_mean[status == "uncalibrated_excluded"] <- 0
  links <- data.frame(resource = e$resource, consumer = e$consumer,
                      evidence_class = e$evidence_class, class_used = class_used,
                      status = status, p_mean = p_mean, inclusion = rowMeans(Z),
                      stringsAsFactors = FALSE)
  background <- NULL
  if (use_decay) {
    bg <- res$bg
    cand <- metaweb$decay
    background <- list(
      pairs = cbind(cand, data.frame(
        p_mean = bg$p_mean,
        inclusion = tabulate(unlist(bg$present, use.names = FALSE), nrow(cand)) / n_draws
      )),
      present = bg$present, p = bg$P, system_effect = bg$v
    )
  }
  structure(
    list(links = links, z = Z, p = res$P, system_effect = res$u, background = background,
         masses = if (use_size) 10^log_masses,
         nodes = metaweb$nodes[intersect(draw_node_columns, names(metaweb$nodes))],
         settings = list(n_draws = n_draws, system = system,
                         parameter_uncertainty = parameter_uncertainty,
                         uncalibrated = uncalibrated, seed = seed,
                         size_model = if (use_size) size_model,
                         rule = metaweb$settings$rule,
                         calibration_sources = calibration$settings$sources)),
    class = "bwebs_draws"
  )
}

# Columns of the nodes kept in the draws: identity, role and classification
# (the classification gives the metabolic types of calc_fluxes()).
draw_node_columns <- c("node_id", "taxon", "rank", "role", "node_type", "genus", "subfamily",
                       "family", "order", "class", "phylum", "life_stage")

# Class whose probability applies to each link of `e` (edges without
# cannibalism; "beyond" uses "class", manual links NA) and its status:
# "calibrated", "manual" or, for links of a class without a coefficient or
# supported only by registries that were not part of the calibration (for
# example GloBI with a calibration of Gray), "uncalibrated_excluded" or
# "uncalibrated_kept" as chosen by `uncalibrated` (NULL stops and lists
# them).
link_status <- function(e, calibration, uncalibrated, call = sys.call(-1)) {
  classes <- calibration$coefficients$evidence_class
  class_used <- e$evidence_class
  class_used[class_used == "beyond"] <- "class"
  class_used[class_used == "manual"] <- NA_character_
  status <- ifelse(is.na(class_used), "manual", "calibrated")
  cal_regs <- calibration$settings$registries
  from_cal <- vapply(strsplit(e$registries, ";", fixed = TRUE),
                     function(r) any(r %in% cal_regs), logical(1))
  missing <- !is.na(class_used) & (!class_used %in% classes | !from_cal)
  if (any(missing)) {
    if (is.null(uncalibrated)) {
      why <- ifelse(!from_cal[missing], "registry not calibrated",
                    paste("class", class_used[missing]))
      tab <- table(why)
      abort_input(sprintf(
        paste("%d link(s) have no calibrated probability (%s). Choose",
              "uncalibrated = \"exclude\" (probability 0) or \"keep\" (probability 1)."),
        sum(missing), paste(names(tab), tab, sep = ": ", collapse = ", ")
      ), call = call)
    }
    status[missing] <- paste0("uncalibrated_", if (uncalibrated == "exclude") "excluded" else "kept")
  }
  list(class_used = class_used, status = status)
}

# Probability of the classes `class_used` (NA for classes without a
# coefficient): in a new system (system = "random", the posterior mean of
# p_new), with a system effect of 0 ("average") or of a given value, averaged
# over the posterior draws (or at the estimates of a calibration without
# draws).
link_class_probability <- function(calibration, class_used, system = "random") {
  co <- calibration$coefficients
  if (identical(system, "random")) {
    p <- if (is.null(co$p_new)) p_new_system(co$estimate, calibration$sigma) else co$p_new
    return(p[match(class_used, co$evidence_class)])
  }
  C <- nrow(co)
  th <- if (is.null(calibration$draws)) {
    matrix(calibration$theta[seq_len(C)], 1L)
  } else {
    calibration$draws[, seq_len(C), drop = FALSE]
  }
  s <- if (identical(system, "average")) 0 else system
  colMeans(stats::plogis(th + s))[match(class_used, co$evidence_class)]
}

# Pairs without a link drawn from the taxonomic decay (`cand`, the `decay`
# table of the metaweb) with the decay model of the calibration (`dk`): in
# each draw, the parameters (from the posterior draws, the approximate normal
# sampling distribution or the estimates) and a system effect of the decay
# model, independent of that of the classes; with `size` (the size model,
# the log masses and the animal nodes), the size term of each pair. Returns
# the probability of each evidence profile in each draw (`P`, without the
# size term), the system effects (`v`), the pairs present in each draw
# (`present`, indices of the rows of `cand`), the profile of each pair
# (`idx`, rows of `P`) and the mean probability of each pair (`p_mean`).
draw_background <- function(cand, dk, n_draws, system, parameter_uncertainty, n_nodes = 50,
                            size = NULL) {
  K <- length(dk$theta)
  th <- if (parameter_uncertainty && !is.null(dk$draws)) {
    t(dk$draws[sample.int(nrow(dk$draws), n_draws, replace = TRUE), , drop = FALSE])
  } else if (parameter_uncertainty) {
    dk$theta + t(chol(dk$vcov)) %*% matrix(stats::rnorm(K * n_draws), K, n_draws)
  } else {
    matrix(dk$theta, K, n_draws)
  }
  th[2:9, ] <- pmax(th[2:9, ], 0)
  sigma <- exp(th[K, ])
  v <- if (identical(system, "random")) {
    stats::rnorm(n_draws, 0, sigma)
  } else if (identical(system, "average")) {
    rep(0, n_draws)
  } else {
    rep(system, n_draws)
  }
  if (!nrow(cand)) {
    return(list(P = matrix(numeric(0), 0L, n_draws), v = v,
                present = replicate(n_draws, integer(0), simplify = FALSE), idx = integer(0),
                p_mean = numeric(0)))
  }
  z <- decay_web_size(n_nodes)
  id <- decay_profile_id(cand$consumer_level, cand$resource_level, cand$generic)
  prof <- sort(unique(id))
  ap <- decay_all_profiles()[prof, , drop = FALSE]
  cl <- match(ap$consumer_level, decay_consumer_levels)
  rg <- decay_resource_group(match(ap$resource_level, decay_levels))
  L <- matrix(0, length(prof), n_draws,
              dimnames = list(paste(ap$consumer_level, ap$resource_level,
                                    ifelse(ap$generic, "generic", "-"), sep = "/"), NULL))
  for (d in seq_len(n_draws)) L[, d] <- decay_linear(th[, d], cl, rg, ap$generic, z) + v[d]
  idx <- match(id, prof)
  if (!is.null(size)) si <- size_index(cand$resource, cand$consumer, size$animal, size$log_masses)
  present <- vector("list", n_draws)
  p_sum <- numeric(length(idx))
  for (d in seq_len(n_draws)) {
    l <- L[idx, d]
    if (!is.null(size)) {
      l <- l + size_shift(size$model, size$model$weights[["decay"]], si, size$log_masses, d)
    }
    p <- stats::plogis(l)
    present[[d]] <- which(stats::runif(length(idx)) < p)
    p_sum <- p_sum + p
  }
  list(P = stats::plogis(L), v = v, present = present, idx = idx, p_mean = p_sum / n_draws)
}

# The probability of a class of evidence depends on the registry that infers
# the links (a larger registry infers more links of each class): warn when the
# metaweb and the calibration were built with different registries.
warn_registries <- function(metaweb, calibration, call = sys.call(-1)) {
  m <- metaweb$settings$registries
  k <- calibration$settings$registries
  if (is.null(m) || is.null(k) || setequal(m, k)) {
    return(invisible(FALSE))
  }
  warning(warningCondition(
    sprintf(paste("The calibration was fitted with the registries %s and the metaweb was built",
                  "with %s; the probability of each class of evidence depends on the",
                  "registries, so calibrate with the registries of the metaweb."),
            format_items(quote_items(sort(k, method = "radix")), Inf),
            format_items(quote_items(sort(m, method = "radix")), Inf)),
    class = c("bwebs_warning_registries", "bwebs_warning"), call = call
  ))
  invisible(TRUE)
}

# Rules that differ only in manual_diet produce the same validation webs;
# the decay does not change the inferred links (whether it is used is
# checked where it applies).
rules_compatible <- function(a, b) {
  a$manual_diet <- character(0)
  b$manual_diet <- character(0)
  a$decay <- NULL
  b$decay <- NULL
  same_rule(a, b)
}

#' @export
print.bwebs_draws <- function(x, ...) {
  s <- x$settings
  bg <- x$background
  n_bg <- if (is.null(bg)) 0L else lengths(bg$present)
  n_links <- colSums(x$z) + n_bg
  cat(sprintf("<bwebs_draws> %d draws of a metaweb with %d nodes and %d candidate links\n",
              s$n_draws, nrow(x$nodes), nrow(x$links)))
  cat(sprintf("  system effect: %s; parameter uncertainty: %s\n",
              if (is.character(s$system)) s$system else format(s$system),
              if (s$parameter_uncertainty) "yes" else "no"))
  q <- stats::quantile(n_links, c(0.025, 0.5, 0.975), names = FALSE)
  cat(sprintf("  links per draw: median %.0f (95%% of draws between %.0f and %.0f)\n",
              q[2], q[1], q[3]))
  if (!is.null(bg)) {
    qb <- stats::quantile(n_bg, c(0.025, 0.5, 0.975), names = FALSE)
    cat(sprintf(paste("  of which from the taxonomic decay: median %.0f (%.0f to %.0f),",
                      "among %d pairs without a link\n"), qb[2], qb[1], qb[3], nrow(bg$pairs)))
  }
  if (!is.null(x$masses)) {
    m <- x$masses
    with_mass <- rowSums(!is.na(m)) > 0
    varies <- any(apply(m[with_mass, , drop = FALSE], 1L, function(r) {
      diff(range(r, na.rm = TRUE)) > 0
    }))
    w <- s$size_model$weights
    cat(sprintf("  body masses: %d of %d nodes, %s\n", sum(with_mass), nrow(m),
                if (varies) "drawn with each draw" else "the same in every draw"))
    cat(sprintf("  size term: weight %.2g on the inferred links, %.2g on the pairs of the decay\n",
                w[["rule"]], w[["decay"]]))
  }
  st <- table(x$links$status)
  cat("  links: ", paste(names(st), st, sep = " = ", collapse = ", "), "\n", sep = "")
  invisible(x)
}
