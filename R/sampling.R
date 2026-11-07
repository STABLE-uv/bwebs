# Food webs sampled from link probabilities ---------------------------------

#' Sample metawebs from calibrated link probabilities
#'
#' `sample_metawebs()` turns a deterministic metaweb into a set of plausible
#' metawebs. In each draw, the probability of every class of evidence is
#' drawn from the calibration (coefficients and, optionally, the effect of
#' the system), and each inferred link is kept with the probability of its
#' class. The same realisation is later shared by all the samples of a
#' community in [calc_indicators()], so that differences between samples are
#' computed within each draw.
#'
#' @param metaweb A metaweb from [build_metaweb()].
#' @param calibration A calibration from [calibrate_links()] fitted with the
#'   same matching rule and registries, such as [gray2015_calibration] for
#'   metawebs built with [gray2015] (a warning is given when the registries
#'   differ).
#' @param n_draws Number of draws (positive whole number).
#' @param system Effect of the system on the logit scale: `"random"` draws it
#'   from the distribution of web effects in each draw (an unknown system,
#'   the more conservative choice, with wide intervals); `"average"` fixes it
#'   at 0 (a system like the average validation web); a number fixes it at
#'   that value.
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
#' * `nodes`: the nodes of the metaweb;
#' * `settings`: the options used.
#'
#' @details
#' In draw \eqn{d}, with class coefficients \eqn{\alpha_c^{(d)}} and system
#' effect \eqn{u^{(d)}}, link \eqn{l} of class \eqn{c(l)} is present with
#' probability
#' \deqn{p_{c(l)}^{(d)} = \mathrm{logit}^{-1}\left(\alpha_{c(l)}^{(d)} + u^{(d)}\right)}{
#' p = plogis(alpha[c(l)] + u)}
#' independently of the other links. Manual links are always present;
#' cannibalistic links are left out, as in the indicators.
#'
#' @section Reproducibility:
#' The same seed can give different draws in other versions of R or with
#' other BLAS/LAPACK libraries.
#'
#' @seealso [calibrate_links()], [calc_indicators()].
#' @export
#' @examples
#' nodes <- prepare_nodes(example_taxa)
#' mw <- build_metaweb(nodes, gray2015)
#' draws <- sample_metawebs(mw, gray2015_calibration, n_draws = 20, seed = 1)
#' draws
sample_metawebs <- function(metaweb, calibration, n_draws, system = "random",
                            parameter_uncertainty = TRUE, uncalibrated = NULL, seed = NULL) {
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
  e <- metaweb$edges[!metaweb$edges$cannibal, , drop = FALSE]
  classes <- calibration$coefficients$evidence_class
  class_used <- e$evidence_class
  class_used[class_used == "beyond"] <- "class"
  class_used[class_used == "manual"] <- NA_character_
  status <- ifelse(is.na(class_used), "manual", "calibrated")
  # Links supported only by registries that were not part of the calibration
  # (for example GloBI) have no calibrated probability either.
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
      ))
    }
    status[missing] <- paste0("uncalibrated_", if (uncalibrated == "exclude") "excluded" else "kept")
  }
  n_draws <- as.integer(n_draws)
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
    P <- stats::plogis(theta[seq_len(C), , drop = FALSE] + matrix(u, C, n_draws, byrow = TRUE))
    dimnames(P) <- list(classes, NULL)
    idx <- match(class_used, classes)
    Z <- matrix(FALSE, nrow(e), n_draws)
    for (d in seq_len(n_draws)) {
      p_l <- P[idx, d]
      Z[, d] <- stats::runif(nrow(e)) < p_l
    }
    list(P = P, u = u, Z = Z)
  }
  res <- if (is.null(seed)) draw() else withr::with_seed(seed, draw())
  Z <- res$Z
  Z[status == "manual", ] <- TRUE
  Z[status == "uncalibrated_kept", ] <- TRUE
  Z[status == "uncalibrated_excluded", ] <- FALSE
  p_mean <- rep(NA_real_, nrow(e))
  ok <- status == "calibrated"
  p_mean[ok] <- rowMeans(res$P[match(class_used[ok], classes), , drop = FALSE])
  p_mean[status %in% c("manual", "uncalibrated_kept")] <- 1
  p_mean[status == "uncalibrated_excluded"] <- 0
  links <- data.frame(resource = e$resource, consumer = e$consumer,
                      evidence_class = e$evidence_class, class_used = class_used,
                      status = status, p_mean = p_mean, inclusion = rowMeans(Z),
                      stringsAsFactors = FALSE)
  structure(
    list(links = links, z = Z, p = res$P, system_effect = res$u,
         nodes = metaweb$nodes[c("node_id", "taxon", "role")],
         settings = list(n_draws = n_draws, system = system,
                         parameter_uncertainty = parameter_uncertainty,
                         uncalibrated = uncalibrated, seed = seed,
                         rule = metaweb$settings$rule,
                         calibration_sources = calibration$settings$sources)),
    class = "bwebs_draws"
  )
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

# Rules that differ only in manual_diet produce the same validation webs.
rules_compatible <- function(a, b) {
  a$manual_diet <- character(0)
  b$manual_diet <- character(0)
  same_rule(a, b)
}

#' @export
print.bwebs_draws <- function(x, ...) {
  s <- x$settings
  n_links <- colSums(x$z)
  cat(sprintf("<bwebs_draws> %d draws of a metaweb with %d nodes and %d candidate links\n",
              s$n_draws, nrow(x$nodes), nrow(x$links)))
  cat(sprintf("  system effect: %s; parameter uncertainty: %s\n",
              if (is.character(s$system)) s$system else format(s$system),
              if (s$parameter_uncertainty) "yes" else "no"))
  q <- stats::quantile(n_links, c(0.025, 0.5, 0.975), names = FALSE)
  cat(sprintf("  links per draw: median %.0f (95%% of draws between %.0f and %.0f)\n",
              q[2], q[1], q[3]))
  st <- table(x$links$status)
  cat("  links: ", paste(names(st), st, sep = " = ", collapse = ", "), "\n", sep = "")
  invisible(x)
}
