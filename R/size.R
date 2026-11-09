# Body-mass ratio of predator and prey -----------------------------------------
#
# The size term (sample_metawebs(masses = ...), score_metaweb(masses = ...)):
# for a predation pair (an animal resource) whose two nodes have a mass, the
# logit of its probability is shifted by w * g(x), with x = log10(M_consumer /
# M_resource), g(x) = b2 * (x - x_opt)^2 <= 0 the within-consumer effect of
# size_model (0 at the optimum ratio x_opt) and w the weight of the links
# inferred by the rule or of the pairs of the taxonomic decay.

# Effect of the mass ratio (log10) on the logit of a link, <= 0.
size_effect <- function(model, x) {
  g <- model$coefficients[["x2"]] * (x - model$optimum)^2
  g[is.na(g)] <- 0
  g
}

# Body masses as a matrix nodes x draws (log10), from a named vector (the
# same in every draw) or a matrix with one column per draw; NA for nodes
# without a mass.
check_masses <- function(masses, node_ids, n_draws = 1L) {
  if (is.data.frame(masses)) masses <- as.matrix(masses)
  if (!is.numeric(masses)) abort_input("`masses` must be numeric.")
  if (length(dim(masses)) < 2L) masses <- stats::setNames(as.vector(masses), names(masses))
  if (is.null(dim(masses))) {
    if (is.null(names(masses))) abort_input("`masses` must be named by node_id.")
    masses <- matrix(masses, length(masses), n_draws, dimnames = list(names(masses), NULL))
  }
  if (is.null(rownames(masses))) abort_input("The rows of `masses` must be named by node_id.")
  if (ncol(masses) != n_draws) {
    abort_input(sprintf("`masses` has %d columns; it needs one per draw (%d), or a vector.",
                        ncol(masses), n_draws))
  }
  unknown <- setdiff(rownames(masses), node_ids)
  if (length(unknown)) {
    abort_input(sprintf("`masses` has nodes that are not in the metaweb: %s.",
                        format_items(unknown, 10L)))
  }
  if (any(duplicated(rownames(masses)))) abort_input("`masses` has repeated nodes.")
  if (any(!is.na(masses) & (!is.finite(masses) | masses <= 0))) {
    abort_input("`masses` must be positive (NA for nodes without a mass).")
  }
  out <- matrix(NA_real_, length(node_ids), n_draws, dimnames = list(node_ids, NULL))
  out[rownames(masses), ] <- log10(masses)
  out
}

check_size_model <- function(size_model) {
  if (is.null(size_model)) size_model <- bwebs::size_model
  if (!inherits(size_model, "bwebs_size_model")) {
    abort_input("`size_model` must be a size model such as `size_model`.")
  }
  w <- size_model$weights
  if (!is.numeric(w) || !all(c("rule", "decay") %in% names(w)) ||
      any(!is.finite(w[c("rule", "decay")])) || any(w[c("rule", "decay")] < 0)) {
    abort_input("The weights of `size_model` must be two numbers >= 0 named rule and decay.")
  }
  size_model
}

# Predation pairs (animal resources) among the pairs (resource, consumer),
# with the rows of their nodes in `log_masses` (from check_masses()).
size_index <- function(resource, consumer, animal, log_masses) {
  pred <- which(resource %in% animal)
  list(pred = pred, ci = match(consumer[pred], rownames(log_masses)),
       ri = match(resource[pred], rownames(log_masses)), n = length(resource))
}

# Shift of the logit of the pairs of `index` by the size term in draw d (0
# for pairs that are not predation or without masses).
size_shift <- function(model, weight, index, log_masses, d) {
  out <- numeric(index$n)
  if (!length(index$pred) || weight == 0) return(out)
  x <- log_masses[index$ci, d] - log_masses[index$ri, d]
  out[index$pred] <- weight * size_effect(model, x)
  out
}

# Probabilities `q` of the pairs (resource, consumer) with the logit shifted
# by the size term (one set of masses, the first column of `log_masses`);
# probabilities 0 and 1 are kept.
size_shift_probability <- function(q, resource, consumer, weight, model, log_masses, animal) {
  s <- size_shift(model, weight, size_index(resource, consumer, animal, log_masses), log_masses, 1L)
  out <- stats::plogis(stats::qlogis(pmin(pmax(q, 1e-12), 1 - 1e-12)) + s)
  out[q %in% c(0, 1)] <- q[q %in% c(0, 1)]
  out
}

# Within-consumer logistic regression of the links of each predator with the
# animals of its web: logit p = a[consumer] + X b, X = (x, x^2) with x =
# log10(M_consumer / M_resource), and one intercept per consumer (fixed
# effects) profiled out. The estimates of b are those of
# glm(y ~ 0 + factor(consumer) + x + I(x^2), binomial), without its design
# matrix (one column per consumer), so that webs with thousands of
# consumers can be fitted. Consumers whose pairs are all links or none carry
# no information on b and are dropped. Used by data-raw/size_model.R and the
# validation scripts.
size_curve_fit <- function(y, x, consumer, max_iter = 100L, tol = 1e-9) {
  ok <- !is.na(y) & !is.na(x) & !is.na(consumer)
  y <- as.numeric(y[ok])
  x <- x[ok]
  consumer <- consumer[ok]
  inf <- tapply(y, consumer, function(v) any(v > 0) && any(v < 1))
  keep <- consumer %in% names(inf)[inf]
  y <- y[keep]
  X <- cbind(x = x[keep], x2 = x[keep]^2)
  cons <- as.integer(factor(consumer[keep]))
  n_c <- max(cons)
  alpha <- stats::qlogis(pmin(pmax(as.vector(rowsum(y, cons)) / tabulate(cons, n_c), 0.01), 0.99))
  solve_alpha <- function(beta, alpha) {
    eta_x <- as.vector(X %*% beta)
    for (k in seq_len(50L)) {
      p <- stats::plogis(alpha[cons] + eta_x)
      step <- as.vector(rowsum(y - p, cons)) / pmax(as.vector(rowsum(p * (1 - p), cons)), 1e-12)
      step <- pmax(pmin(step, 5), -5)
      alpha <- alpha + step
      if (max(abs(step)) < 1e-10) break
    }
    alpha
  }
  loglik <- function(beta, alpha) {
    eta <- alpha[cons] + as.vector(X %*% beta)
    sum(y * eta - log1p(exp(eta)))
  }
  beta <- c(x = 0, x2 = 0)
  alpha <- solve_alpha(beta, alpha)
  ll <- loglik(beta, alpha)
  converged <- FALSE
  for (it in seq_len(max_iter)) {
    p <- stats::plogis(alpha[cons] + as.vector(X %*% beta))
    w <- p * (1 - p)
    U <- as.vector(crossprod(X, y - p))
    WX <- rowsum(w * X, cons)
    W <- as.vector(rowsum(w, cons))
    info <- crossprod(X, w * X) - crossprod(WX / sqrt(pmax(W, 1e-12)))
    delta <- as.vector(solve(info, U))
    step <- 1
    repeat {
      b_new <- beta + step * delta
      a_new <- solve_alpha(b_new, alpha)
      ll_new <- loglik(b_new, a_new)
      if (ll_new >= ll - 1e-12 || step < 1e-4) break
      step <- step / 2
    }
    beta <- stats::setNames(b_new, c("x", "x2"))
    alpha <- a_new
    done <- abs(ll_new - ll) < tol * (abs(ll) + tol) && max(abs(step * delta)) < 1e-7
    ll <- ll_new
    if (done) {
      converged <- TRUE
      break
    }
  }
  p <- stats::plogis(alpha[cons] + as.vector(X %*% beta))
  w <- p * (1 - p)
  WX <- rowsum(w * X, cons)
  info <- crossprod(X, w * X) - crossprod(WX / sqrt(pmax(as.vector(rowsum(w, cons)), 1e-12)))
  vc <- solve(info)
  dimnames(vc) <- list(c("x", "x2"), c("x", "x2"))
  list(coefficients = beta, std_error = sqrt(diag(vc)), vcov = vc, optimum = -beta[["x"]] / (2 * beta[["x2"]]),
       pairs = length(y), links = sum(y), consumers = n_c, loglik = ll, converged = converged)
}

#' @export
print.bwebs_size_model <- function(x, ...) {
  cat("<bwebs_size_model> effect of the predator/prey body-mass ratio on links\n")
  cat(sprintf("  logit shift: %.3f * (log10(M_consumer / M_resource) - %.2f)^2\n",
              x$coefficients[["x2"]], x$optimum))
  cat(sprintf("  optimum: consumer %.0f times heavier than its prey\n", 10^x$optimum))
  cat(sprintf("  weights: links inferred by the rule %.2f; pairs of the decay %.2f\n",
              x$weights[["rule"]], x$weights[["decay"]]))
  n_webs <- if ("webs" %in% names(x$data)) sum(x$data$webs) else nrow(x$data)
  cat(sprintf("  fitted to %d predator x animal pairs (%d links) of %d webs\n",
              sum(x$data$pairs), sum(x$data$links), n_webs))
  invisible(x)
}
