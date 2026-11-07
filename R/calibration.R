# Calibration of link probabilities by class of evidence ------------------------

#' Calibrate link probabilities by class of evidence
#'
#' `calibrate_links()` estimates how often links inferred with each class of
#' evidence (see [build_metaweb()]) are observed in independent, well-studied
#' communities, from the held-out webs of [validate_rule()]: each validation
#' study of the registry is rebuilt from the rest of the registry (without
#' the studies of its system) and its inferred links are compared with the
#' recorded ones, as in the validation of Gray et al. (2015, section 2.3,
#' p. 15). A logistic model with one coefficient per class of evidence and a
#' random effect of the validation web is fitted to the counts, by default as
#' a Bayesian model whose posterior is sampled by Markov chain Monte Carlo.
#'
#' @param x A [validate_rule()] result, or a registry in the format of
#'   [as_registry()]; with a registry, [validate_rule()] is run first with
#'   `sources`, `groups` and `rule`.
#' @param sources,groups,rule Passed to [validate_rule()] when `x` is a
#'   registry; with a validation they are taken from it and must not be given.
#'   The calibration only applies to metawebs built with the same rule.
#' @param method `"bayes"` (default): posterior by MCMC with the priors below;
#'   `"ml"`: maximum marginal likelihood with a normal approximation.
#' @param prior_sd Standard deviation of the normal prior, centred at 0, of
#'   each class coefficient on the logit scale (positive number).
#' @param prior_sigma Scale of the half-normal prior of the standard deviation
#'   of the web effect, on the logit scale (positive number).
#' @param n_chains,n_iter,n_warmup Number of chains, of draws kept per chain
#'   and of warm-up iterations per chain (adaptive random-walk Metropolis).
#' @param quadrature Number of adaptive Gauss-Hermite quadrature points used
#'   to integrate the web effects (integer >= 1; 1 is the Laplace
#'   approximation).
#' @param seed `NULL` (default) to use the current random number stream, or a
#'   number (`method = "bayes"`). When given, the result is reproducible and
#'   the global random number stream is restored afterwards.
#' @param quiet If `FALSE`, report progress with messages.
#'
#' @return An object of class `bwebs_calibration`, a list with
#' * `coefficients`: data frame with one row per calibrated class
#'   (`evidence_class`; `estimate`, `std_error`, `lower` and `upper` of the
#'   coefficient on the logit scale: posterior mean, standard deviation and
#'   95% credible interval, or estimate, standard error and Wald interval with
#'   `method = "ml"`; `p_average`, the probability in a web with effect 0;
#'   `p_new`, the probability in a new system, averaged over web effects; and
#'   the number of `webs`, `predicted` and `observed` links);
#' * `sigma`: standard deviation of the web effect (posterior mean, or
#'   estimate);
#' * `theta` and `vcov`: class coefficients and `log_sigma` (posterior means,
#'   or estimates) and their covariance (posterior, or inverse observed
#'   information);
#' * `draws`: matrix of posterior draws of `theta` (`NULL` with
#'   `method = "ml"`);
#' * `counts` and `webs`: from the validation, with the web effects
#'   (conditional modes at `theta`);
#' * `settings` (rule, sources, groups, exclusions, registries, method, priors,
#'   chains)
#'   and `diagnostics` (convergence, positive-definite information at the
#'   mode, `rhat` and `ess` of every parameter, acceptance rate of each chain,
#'   classes without data).
#'
#' @details
#' For the validation web \eqn{w}, with \eqn{y_{cw}} of the \eqn{m_{cw}}
#' links of class \eqn{c} inferred for \eqn{w} recorded in \eqn{w}:
#' \deqn{y_{cw} \sim \mathrm{Binomial}(m_{cw}, p_{cw}), \quad
#' \mathrm{logit}\, p_{cw} = \alpha_c + u_w, \quad u_w \sim N(0, \sigma^2)}{
#' y[c,w] ~ Binomial(m[c,w], p[c,w]), logit(p[c,w]) = alpha[c] + u[w],
#' u[w] ~ N(0, sigma^2)}
#' with priors \eqn{\alpha_c \sim N(0, s_\alpha^2)}{alpha[c] ~ N(0,
#' prior_sd^2)} and \eqn{\sigma \sim \textrm{half-}N(0, s_\sigma^2)}{sigma ~
#' half-normal(0, prior_sigma^2)}. The web effects are integrated out with
#' adaptive Gauss-Hermite quadrature centred at their conditional modes, and
#' the posterior of \eqn{(\alpha, \log\sigma)} is sampled with an adaptive
#' random-walk Metropolis algorithm started at the posterior mode; the
#' proposal is tuned during the warm-up and fixed afterwards. Check
#' `diagnostics$rhat` (close to 1) and `diagnostics$ess`.
#'
#' Interpretation. The probability is that of observing an inferred link in a
#' comparable, well-studied community; because the validation webs miss rare
#' links, it is a lower bound for the existence of the link. Pairs without an
#' inferred link are not modelled: in a registry web, a taxon is listed
#' because it eats or is eaten, so the proportion of recorded links among the
#' pairs not inferred in those webs does not transfer to communities sampled
#' independently of their diets.
#'
#' @section Missing values and zeros:
#' A class of evidence that never occurs in the validation webs gets no
#' coefficient; [sample_metawebs()] then asks what to do with its links.
#'
#' @section Assumptions:
#' Links are recorded independently given their class and web; the
#' validation webs are a sample of the communities to which the calibration
#' is applied; the links recorded in a validation web are true links; no
#' study left in the registry describes the same system as the web left out
#' (see `groups` in [validate_rule()]).
#'
#' @references
#' Gray C, Figueroa DH, Hudson LN, Ma A, Perkins D, Woodward G (2015).
#' Joining the dots: an automated method for constructing food webs from
#' compendia of published interactions. *Food Webs*, 5, 11--20.
#' \doi{10.1016/j.fooweb.2015.09.001}
#'
#' @seealso [gray2015_calibration] for the calibration of the bundled
#'   registry; [sample_metawebs()]; [score_rule()].
#' @export
#' @examplesIf interactive()
#' v <- validate_rule(gray2015, sources = c("46", "49", "58", "59", "63"))
#' cal <- calibrate_links(v, seed = 1)
#' cal
calibrate_links <- function(x, sources, groups = NULL, rule = matching_rule(),
                            method = "bayes", prior_sd = 1.5, prior_sigma = 2.5,
                            n_chains = 4L, n_iter = 2000L, n_warmup = 1000L,
                            quadrature = 20L, seed = NULL, quiet = FALSE) {
  if (inherits(x, "bwebs_validation")) {
    if (!missing(sources) || !is.null(groups) || !missing(rule)) {
      abort_input(paste("`sources`, `groups` and `rule` are taken from the validation;",
                        "do not give them with a validate_rule() result."))
    }
    validation <- x
  } else {
    if (missing(sources)) {
      abort_input("`sources` must be a character vector of study identifiers (study_id).")
    }
    validation <- validate_rule(x, sources = sources, groups = groups, rule = rule, quiet = quiet)
  }
  check_choice(method, c("bayes", "ml"))
  check_number(prior_sd, lower = 0, include_lower = FALSE)
  check_number(prior_sigma, lower = 0, include_lower = FALSE)
  check_number(n_chains, lower = 1)
  check_number(n_iter, lower = 10)
  check_number(n_warmup, lower = 0)
  check_number(quadrature, lower = 1)
  check_number(seed, allow_null = TRUE)
  check_flag(quiet)
  fit <- if (method == "ml") {
    fit_calibration(validation$counts, quadrature)
  } else {
    if (!quiet) message("Sampling the posterior of the calibration model")
    fit_calibration_bayes(validation$counts, prior_sd = prior_sd, prior_sigma = prior_sigma,
                          n_chains = as.integer(n_chains), n_iter = as.integer(n_iter),
                          n_warmup = as.integer(n_warmup), quadrature = quadrature, seed = seed)
  }
  webs <- validation$webs
  webs$web_effect <- fit$web_effect[match(webs$study_id, names(fit$web_effect))]
  s <- validation$settings
  structure(
    list(
      coefficients = fit$coefficients, sigma = fit$sigma, theta = fit$theta,
      vcov = fit$vcov, draws = fit$draws, counts = validation$counts, webs = webs,
      settings = list(rule = s$rule, sources = s$sources, groups = s$groups,
                      exclusions = s$exclusions,
                      registries = s$registries, n_records = s$n_records,
                      method = method, quadrature = as.integer(quadrature),
                      prior = if (method == "bayes") list(prior_sd = prior_sd, prior_sigma = prior_sigma),
                      mcmc = if (method == "bayes") list(n_chains = as.integer(n_chains),
                                                         n_iter = as.integer(n_iter),
                                                         n_warmup = as.integer(n_warmup),
                                                         seed = seed)),
      diagnostics = list(converged = fit$converged, information_pd = fit$information_pd,
                         message = fit$message, rhat = fit$rhat, ess = fit$ess,
                         acceptance = fit$acceptance,
                         classes_without_data = setdiff(evidence_classes,
                                                        fit$coefficients$evidence_class))
    ),
    class = "bwebs_calibration"
  )
}

# Matrices of observed (Y) and predicted (M) links, webs x classes.
count_matrices <- function(counts) {
  if (nrow(counts) == 0L || sum(counts$predicted) == 0) {
    abort_input("No link was inferred in the validation webs: nothing to calibrate.")
  }
  classes <- evidence_classes[evidence_classes %in% counts$evidence_class]
  webs <- sort(unique(counts$study_id), method = "radix")
  Y <- M <- matrix(0, length(webs), length(classes), dimnames = list(webs, classes))
  idx <- cbind(match(counts$study_id, webs), match(counts$evidence_class, classes))
  Y[idx] <- counts$observed
  M[idx] <- counts$predicted
  list(Y = Y, M = M, classes = classes, webs = webs)
}

# Logistic model with class coefficients and a normal random effect of the
# web, fitted by maximum marginal likelihood (adaptive Gauss-Hermite).
fit_calibration <- function(counts, quadrature) {
  cm <- count_matrices(counts)
  Y <- cm$Y
  M <- cm$M
  classes <- cm$classes
  C <- length(classes)
  gh <- gauss_hermite(quadrature)
  lchoose_sum <- sum(lchoose(M, Y))
  nll <- function(theta) -marginal_loglik(theta, Y, M, gh) - lchoose_sum
  pooled <- (colSums(Y) + 0.5) / (colSums(M) + 1)
  start <- c(stats::qlogis(pooled), log_sigma = 0)
  names(start) <- c(classes, "log_sigma")
  opt <- stats::optim(start, nll, method = "BFGS", control = list(maxit = 1000L, reltol = 1e-12))
  theta <- opt$par
  H <- stats::optimHess(theta, nll)
  info_pd <- tryCatch({ chol(H); TRUE }, error = function(e) FALSE)
  V <- if (info_pd) solve(H) else matrix(NA_real_, length(theta), length(theta))
  dimnames(V) <- list(names(theta), names(theta))
  alpha <- theta[seq_len(C)]
  se <- sqrt(pmax(diag(V)[seq_len(C)], 0))
  sigma <- exp(theta[["log_sigma"]])
  coefs <- data.frame(
    evidence_class = classes, estimate = unname(alpha), std_error = unname(se),
    lower = unname(alpha - 1.959964 * se), upper = unname(alpha + 1.959964 * se),
    p_average = unname(stats::plogis(alpha)),
    p_new = unname(p_new_system(alpha, sigma)),
    webs = unname(colSums(M > 0)), predicted = unname(colSums(M)), observed = unname(colSums(Y)),
    stringsAsFactors = FALSE
  )
  list(coefficients = coefs, sigma = sigma, theta = theta, vcov = V, draws = NULL,
       web_effect = web_effects(theta, Y, M),
       converged = opt$convergence == 0L, information_pd = info_pd,
       message = if (is.null(opt$message)) NA_character_ else opt$message,
       rhat = NULL, ess = NULL, acceptance = NULL, loglik = -opt$value)
}

# Bayesian version: normal priors on the class coefficients, half-normal
# prior on sigma, posterior of (alpha, log sigma) by adaptive random-walk
# Metropolis on the marginal posterior (web effects integrated out).
fit_calibration_bayes <- function(counts, prior_sd, prior_sigma, n_chains, n_iter, n_warmup,
                                  quadrature, seed) {
  cm <- count_matrices(counts)
  Y <- cm$Y
  M <- cm$M
  classes <- cm$classes
  C <- length(classes)
  pnames <- c(classes, "log_sigma")
  gh <- gauss_hermite(quadrature)
  lp <- function(theta) {
    ls <- theta[[C + 1L]]
    if (!is.finite(ls) || ls > 5 || ls < -10) return(-Inf)
    marginal_loglik(theta, Y, M, gh) +
      sum(stats::dnorm(theta[seq_len(C)], 0, prior_sd, log = TRUE)) +
      stats::dnorm(exp(ls), 0, prior_sigma, log = TRUE) + ls
  }
  pooled <- (colSums(Y) + 0.5) / (colSums(M) + 1)
  start <- stats::setNames(c(stats::qlogis(pooled), 0), pnames)
  opt <- stats::optim(start, function(t) -lp(t), method = "BFGS",
                      control = list(maxit = 1000L, reltol = 1e-12))
  mode <- opt$par
  H <- stats::optimHess(mode, function(t) -lp(t))
  info_pd <- tryCatch({ chol(H); TRUE }, error = function(e) FALSE)
  V0 <- if (info_pd) solve(H) else diag(0.05, C + 1L)
  run <- function() {
    lapply(seq_len(n_chains), function(k) {
      init <- mode + drop(crossprod(chol(V0), stats::rnorm(C + 1L)))
      metropolis_chain(init, V0, n_iter, n_warmup, lp)
    })
  }
  chains <- if (is.null(seed)) run() else withr::with_seed(seed, run())
  arr <- array(NA_real_, c(n_iter, n_chains, C + 1L), dimnames = list(NULL, NULL, pnames))
  for (k in seq_len(n_chains)) arr[, k, ] <- chains[[k]]$draws
  draws <- matrix(aperm(arr, c(1, 2, 3)), n_iter * n_chains, C + 1L,
                  dimnames = list(NULL, pnames))
  rhat <- vapply(pnames, function(p) rhat_split(arr[, , p, drop = TRUE]), numeric(1))
  ess <- vapply(pnames, function(p) ess_bulk(arr[, , p, drop = TRUE]), numeric(1))
  theta <- colMeans(draws)
  alpha_d <- draws[, seq_len(C), drop = FALSE]
  sigma_d <- exp(draws[, C + 1L])
  pnew_d <- vapply(seq_len(C), function(c) p_new_system(alpha_d[, c], sigma_d), numeric(nrow(draws)))
  pnew_d <- matrix(pnew_d, nrow(draws), C)
  q <- apply(alpha_d, 2, stats::quantile, probs = c(0.025, 0.975), names = FALSE)
  coefs <- data.frame(
    evidence_class = classes, estimate = unname(colMeans(alpha_d)),
    std_error = unname(apply(alpha_d, 2, stats::sd)),
    lower = unname(q[1, ]), upper = unname(q[2, ]),
    p_average = unname(colMeans(stats::plogis(alpha_d))),
    p_new = unname(colMeans(pnew_d)),
    webs = unname(colSums(M > 0)), predicted = unname(colSums(M)), observed = unname(colSums(Y)),
    stringsAsFactors = FALSE
  )
  V <- stats::cov(draws)
  acceptance <- vapply(chains, function(ch) ch$acceptance, numeric(1))
  list(coefficients = coefs, sigma = mean(sigma_d), theta = theta, vcov = V, draws = draws,
       web_effect = web_effects(theta, Y, M),
       converged = all(rhat < 1.01) && all(ess > 400), information_pd = info_pd,
       message = if (is.null(opt$message)) NA_character_ else opt$message,
       rhat = rhat, ess = ess, acceptance = acceptance)
}

# One chain of Metropolis-Hastings on the marginal posterior. Warm-up:
# random-walk proposals whose covariance is that of the chain so far (from
# iteration 200, every 100 iterations) and whose scale is adapted towards an
# acceptance rate of 0.234. After the warm-up each iteration makes an
# independence step, from a multivariate t (5 degrees of freedom) with the
# mean and 1.5 times the covariance of the second half of the warm-up, and a
# random-walk step with the final warm-up proposal; both proposals are then
# fixed, so the chain is a valid Markov chain for the posterior.
metropolis_chain <- function(init, V0, n_iter, n_warmup, lp) {
  d <- length(init)
  S <- V0
  log_scale <- log(2.38^2 / d)
  R <- chol(exp(log_scale) * S)
  cur <- init
  cur_lp <- lp(cur)
  tries <- 0L
  while (!is.finite(cur_lp) && tries < 100L) {
    cur <- init + drop(crossprod(R, stats::rnorm(d))) * 0.1
    cur_lp <- lp(cur)
    tries <- tries + 1L
  }
  hist <- matrix(NA_real_, max(n_warmup, 1L), d)
  for (t in seq_len(n_warmup)) {
    prop <- cur + drop(crossprod(R, stats::rnorm(d)))
    prop_lp <- lp(prop)
    ok <- is.finite(prop_lp) && log(stats::runif(1)) < prop_lp - cur_lp
    if (ok) {
      cur <- prop
      cur_lp <- prop_lp
    }
    hist[t, ] <- cur
    log_scale <- log_scale + (as.numeric(ok) - 0.234) / sqrt(t)
    if (t >= 200L && t %% 100L == 0L) {
      Sn <- stats::cov(hist[ceiling(t / 2):t, , drop = FALSE]) + diag(1e-6, d)
      if (all(is.finite(Sn))) S <- Sn
    }
    Rn <- tryCatch(chol(exp(log_scale) * S), error = function(e) NULL)
    if (!is.null(Rn)) R <- Rn
  }
  # independence proposal
  nu <- 5
  if (n_warmup >= 200L) {
    half <- hist[ceiling(n_warmup / 2):n_warmup, , drop = FALSE]
    mu <- colMeans(half)
    Sig <- 1.5 * (stats::cov(half) + diag(1e-6, d))
  } else {
    mu <- init
    Sig <- 1.5 * V0
  }
  Ri <- tryCatch(chol(Sig), error = function(e) chol(1.5 * V0))
  Si_inv <- chol2inv(Ri)
  lq <- function(x) {
    z <- x - mu
    -(nu + d) / 2 * log1p(drop(crossprod(z, Si_inv %*% z)) / nu)
  }
  cur_lq <- lq(cur)
  out <- matrix(NA_real_, n_iter, d)
  accepted <- 0L
  for (t in seq_len(n_iter)) {
    prop <- mu + drop(crossprod(Ri, stats::rnorm(d))) * sqrt(nu / stats::rchisq(1, nu))
    prop_lp <- lp(prop)
    prop_lq <- lq(prop)
    if (is.finite(prop_lp) &&
        log(stats::runif(1)) < (prop_lp - prop_lq) - (cur_lp - cur_lq)) {
      cur <- prop
      cur_lp <- prop_lp
      cur_lq <- prop_lq
      accepted <- accepted + 1L
    }
    prop <- cur + drop(crossprod(R, stats::rnorm(d)))
    prop_lp <- lp(prop)
    if (is.finite(prop_lp) && log(stats::runif(1)) < prop_lp - cur_lp) {
      cur <- prop
      cur_lp <- prop_lp
      cur_lq <- lq(cur)
    }
    out[t, ] <- cur
  }
  list(draws = out, acceptance = accepted / n_iter)
}

# Split R-hat of a matrix of draws (iterations x chains).
rhat_split <- function(x) {
  x <- as.matrix(x)
  n <- floor(nrow(x) / 2)
  if (n < 2L) return(NA_real_)
  s <- cbind(x[seq_len(n), , drop = FALSE], x[nrow(x) - n + seq_len(n), , drop = FALSE])
  W <- mean(apply(s, 2, stats::var))
  B <- n * stats::var(colMeans(s))
  if (W <= 0) return(NA_real_)
  sqrt(((n - 1) / n * W + B / n) / W)
}

# Effective sample size of a matrix of draws (iterations x chains), from the
# autocorrelations of the split chains, truncated at the first negative sum
# of two consecutive autocorrelations.
ess_bulk <- function(x) {
  x <- as.matrix(x)
  n <- floor(nrow(x) / 2)
  if (n < 4L) return(NA_real_)
  s <- cbind(x[seq_len(n), , drop = FALSE], x[nrow(x) - n + seq_len(n), , drop = FALSE])
  m <- ncol(s)
  W <- mean(apply(s, 2, stats::var))
  B <- n * stats::var(colMeans(s))
  var_plus <- (n - 1) / n * W + B / n
  if (var_plus <= 0) return(NA_real_)
  rho <- numeric(0)
  for (t in seq_len(n - 1L)) {
    vt <- mean(vapply(seq_len(m), function(j) mean((s[(t + 1):n, j] - s[1:(n - t), j])^2),
                      numeric(1)))
    rho[t] <- 1 - vt / (2 * var_plus)
    if (t %% 2L == 0L && rho[t - 1L] + rho[t] < 0) {
      rho <- rho[seq_len(t - 2L)]
      break
    }
  }
  m * n / (1 + 2 * sum(rho))
}

# Probability of a class in a new system: mean of plogis(alpha + u) over
# u ~ N(0, sigma^2), by Gauss-Hermite quadrature (vectorised over alpha and
# sigma of the same length, or a scalar sigma).
p_new_system <- function(alpha, sigma, n = 30L) {
  gh <- gauss_hermite(n)
  sigma <- rep_len(sigma, length(alpha))
  vapply(seq_along(alpha), function(i) {
    sum(gh$w * stats::plogis(alpha[i] + sqrt(2) * sigma[i] * gh$x)) / sqrt(pi)
  }, numeric(1))
}

# Conditional modes of the web effects at theta.
web_effects <- function(theta, Y, M) {
  C <- ncol(Y)
  alpha <- theta[seq_len(C)]
  sigma <- exp(theta[[C + 1L]])
  out <- vapply(seq_len(nrow(Y)), function(w) web_mode(alpha, sigma, Y[w, ], M[w, ])$u, numeric(1))
  stats::setNames(out, rownames(Y))
}

# Nodes and weights of Gauss-Hermite quadrature (weight exp(-x^2)), by the
# Golub-Welsch algorithm.
gauss_hermite <- function(n) {
  n <- as.integer(n)
  if (n == 1L) {
    return(list(x = 0, w = sqrt(pi)))
  }
  b <- sqrt(seq_len(n - 1L) / 2)
  J <- matrix(0, n, n)
  J[cbind(seq_len(n - 1L), 2:n)] <- b
  J[cbind(2:n, seq_len(n - 1L))] <- b
  ev <- eigen(J, symmetric = TRUE)
  o <- order(ev$values)
  list(x = ev$values[o], w = sqrt(pi) * ev$vectors[1L, o]^2)
}

# log(1 + exp(x)) without overflow.
softplus <- function(x) {
  pmax(x, 0) + log1p(exp(-abs(x)))
}

# Mode of the conditional log density of a web effect, by Newton's method.
web_mode <- function(alpha, sigma, y, m) {
  u <- 0
  for (it in seq_len(100L)) {
    p <- stats::plogis(alpha + u)
    g1 <- sum(y - m * p) - u / sigma^2
    g2 <- -sum(m * p * (1 - p)) - 1 / sigma^2
    step <- g1 / g2
    step <- max(min(step, 5), -5)
    u <- u - step
    if (abs(step) < 1e-10) break
  }
  p <- stats::plogis(alpha + u)
  list(u = u, curvature = -sum(m * p * (1 - p)) - 1 / sigma^2)
}

# Marginal log likelihood (up to the binomial coefficients) with the web
# effects integrated by adaptive Gauss-Hermite quadrature.
marginal_loglik <- function(theta, Y, M, gh) {
  C <- ncol(Y)
  alpha <- theta[seq_len(C)]
  sigma <- exp(theta[[C + 1L]])
  lw <- log(gh$w) + gh$x^2
  total <- 0
  for (w in seq_len(nrow(Y))) {
    y <- Y[w, ]
    m <- M[w, ]
    md <- web_mode(alpha, sigma, y, m)
    s <- 1 / sqrt(-md$curvature)
    u <- md$u + sqrt(2) * s * gh$x
    eta <- outer(alpha, u, "+")
    g <- colSums(y * eta - m * softplus(eta)) - u^2 / (2 * sigma^2)
    a <- lw + g
    lse <- max(a) + log(sum(exp(a - max(a))))
    total <- total + log(sqrt(2) * s) - log(sqrt(2 * pi) * sigma) + lse
  }
  total
}

#' @export
print.bwebs_calibration <- function(x, ...) {
  w <- x$webs
  method <- if (is.null(x$settings$method)) "ml" else x$settings$method
  cat(sprintf("<bwebs_calibration> %d validation webs, %d links inferred, %.0f%% observed\n",
              nrow(w), sum(w$n_predicted), 100 * sum(w$n_predicted_observed) / sum(w$n_predicted)))
  if (identical(method, "bayes")) {
    m <- x$settings$mcmc
    cat(sprintf("  Bayesian: %d chains x %d draws; max R-hat %.3f, min ESS %.0f\n",
                m$n_chains, m$n_iter, max(x$diagnostics$rhat), min(x$diagnostics$ess)))
  } else {
    cat("  maximum marginal likelihood\n")
  }
  cat(sprintf("  web effect: sd = %.2f on the logit scale; precision per web from %.2f to %.2f\n",
              x$sigma, min(w$precision, na.rm = TRUE), max(w$precision, na.rm = TRUE)))
  co <- x$coefficients
  p_new <- if (is.null(co$p_new)) p_new_system(co$estimate, x$sigma) else co$p_new
  out <- data.frame(class = co$evidence_class, p_average = round(co$p_average, 3),
                    p_new = round(p_new, 3), logit = round(co$estimate, 2),
                    sd = round(co$std_error, 2), webs = co$webs, predicted = co$predicted,
                    observed = co$observed)
  print(out, row.names = FALSE)
  if (length(x$settings$groups)) {
    cat("  studies left out together: ",
        paste(vapply(x$settings$groups, paste, character(1), collapse = " + "), collapse = "; "),
        "\n", sep = "")
  }
  if (length(x$diagnostics$classes_without_data)) {
    cat("  classes without validation data: ",
        paste(x$diagnostics$classes_without_data, collapse = ", "), "\n", sep = "")
  }
  if (identical(method, "bayes") && !isTRUE(x$diagnostics$converged)) {
    cat("  warning: R-hat above 1.01 or effective sample size below 400; run longer chains (n_iter)\n")
  } else if (!identical(method, "bayes") &&
             (!isTRUE(x$diagnostics$converged) || !isTRUE(x$diagnostics$information_pd))) {
    cat("  warning: the fit did not converge or its information matrix is not positive definite\n")
  }
  invisible(x)
}
