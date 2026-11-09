# Experimental versions of the taxonomic distance decay (first exploration,
# before the decay was part of the package): models over every consumer
# level ("linear", "steps", "steps_known"), the decay used instead of the
# matching rule and the decay as the background of the rule. Used only by
# dev/validation/decay_background.R; the package uses the model of R/decay.R
# (consumer resolved from the family, see decay_evidence()).
#
# Needs the package loaded with pkgload::load_all() (internal functions).

# Distance of each level on the scale of the "linear" decay (a subfamily lies
# between the genus and the family).
decay_scores <- c(exact = 0, genus = 1, subfamily = 1.5, family = 2, order = 3, class = 4,
                  phylum = 5, any = 6)

# Profiles of a set of pairs as matrices (pairs x consumer levels) and a key.
decay_pair_profiles <- function(profiles, resource, consumer) {
  ri <- match(resource, dimnames(profiles$level)$resource)
  ci <- match(consumer, dimnames(profiles$level)$consumer)
  L <- length(decay_levels)
  P <- matrix(NA_integer_, length(ri), L, dimnames = list(NULL, decay_levels))
  G <- matrix(FALSE, length(ri), L, dimnames = list(NULL, decay_levels))
  for (c in seq_len(L)) {
    P[, c] <- profiles$level[cbind(ri, ci, c)]
    G[, c] <- profiles$generic[cbind(ri, ci, c)]
  }
  # resources known at the levels closer than each level (0 for the consumer itself)
  kn <- profiles$known[ci, , drop = FALSE]
  Q <- log1p(cbind(0, kn[, -L, drop = FALSE]))
  dimnames(Q) <- list(NULL, decay_levels)
  P0 <- P
  P0[is.na(P0)] <- 0L
  key <- paste(apply(P0, 1L, paste, collapse = ""), apply(G * 1L, 1L, paste, collapse = ""),
               apply(cbind(0, kn[, -L, drop = FALSE]), 1L, paste, collapse = ","), sep = "|")
  list(P = P, G = G, Q = Q, key = key)
}

# Parameters of the decay model. "linear": a0 - b_con * d_con - b_res * d_res
# with the distances of decay_scores; "steps": a free, non-increasing effect of
# each level on each side (a0 for the consumer itself, then steps exp(.) down
# to "any"; 0 for the resource itself, then steps down). Both have an effect
# of the records of a generic group (relative to the resource itself) and the
# log standard deviation of the web effect.
decay_parameter_names <- function(model) {
  switch(model,
    linear = c("a0", "log_b_con", "log_b_res", "d_generic", "log_sigma"),
    steps = c("a0", paste0("con_step_", decay_levels[-1L]), paste0("res_step_", decay_levels[-1L]),
              "d_generic", "log_sigma"),
    steps_known = c("a0", paste0("con_step_", decay_levels[-1L]),
                    paste0("res_step_", decay_levels[-1L]), "d_generic", "log_gamma", "log_sigma")
  )
}

# Effects of the consumer levels (a) and of the resource levels (b) on the
# logit scale, from theta.
decay_effects <- function(theta, model) {
  L <- length(decay_levels)
  if (model == "linear") {
    a <- theta[[1L]] - exp(theta[[2L]]) * decay_scores
    b <- -exp(theta[[3L]]) * decay_scores
    g <- -theta[[4L]]
  } else {
    a <- theta[[1L]] - cumsum(c(0, exp(theta[2:L])))
    b <- -cumsum(c(0, exp(theta[(L + 1L):(2L * L - 1L)])))
    g <- -theta[[2L * L]]
  }
  gamma <- if (model == "steps_known") exp(theta[[2L * L + 1L]]) else 0
  list(a = unname(a), b = unname(b), g = unname(g), gamma = gamma)
}

# Linear predictor of each profile: the best evidence among its consumer
# levels, a[consumer level] + b[resource level], or a[consumer level] - d_gen
# for records of a generic group. With `design = TRUE`, also the derivatives
# of the linear predictor with respect to theta (without log_sigma), from the
# evidence that gives the maximum.
decay_eta <- function(theta, P, G, model = "steps", design = FALSE, Q = NULL) {
  ef <- decay_effects(theta, model)
  n <- nrow(P)
  L <- length(decay_levels)
  A <- matrix(ef$a, n, L, byrow = TRUE)
  if (ef$gamma > 0) A <- A - ef$gamma * Q
  E <- A + matrix(ef$b[P], n, L)
  E[is.na(E)] <- -Inf
  Eg <- A + ef$g
  Eg[!G] <- -Inf
  cand <- cbind(E, Eg)
  best <- max.col(cand, ties.method = "first")
  eta <- cand[cbind(seq_len(n), best)]
  if (!design) return(eta)
  generic <- best > L
  cl <- (best - 1L) %% L + 1L
  rl <- P[cbind(seq_len(n), cl)]
  rl[generic | is.na(rl)] <- 1L
  K <- length(theta) - 1L
  D <- matrix(0, n, K)
  D[, 1L] <- 1
  if (model == "linear") {
    D[, 2L] <- -exp(theta[[2L]]) * decay_scores[cl]
    D[, 3L] <- -exp(theta[[3L]]) * decay_scores[rl]
    D[, 4L] <- -generic
  } else {
    for (k in 2:L) {
      D[, k] <- -exp(theta[[k]]) * (cl >= k)
      D[, L - 1L + k] <- -exp(theta[[L - 1L + k]]) * (rl >= k & !generic)
    }
    D[, 2L * L] <- -generic
    if (model == "steps_known") D[, 2L * L + 1L] <- -ef$gamma * Q[cbind(seq_len(n), cl)]
  }
  list(eta = eta, D = D)
}

# Marginal log likelihood of the decay model over webs (web effects
# integrated by adaptive Gauss-Hermite quadrature, as in marginal_loglik()),
# and its gradient with respect to theta (the quadrature nodes held at the
# conditional mode of each web effect).
decay_loglik <- function(theta, webs, gh, model = "steps", gradient = FALSE) {
  K <- length(theta)
  sigma <- exp(theta[[K]])
  lw <- log(gh$w) + gh$x^2
  total <- 0
  grad <- numeric(K)
  for (w in webs) {
    ed <- decay_eta(theta, w$P, w$G, model, design = gradient, Q = w$Q)
    eta <- if (gradient) ed$eta else ed
    eta <- pmax(eta, -50)
    md <- web_mode(eta, sigma, w$y, w$m)
    s <- 1 / sqrt(-md$curvature)
    u <- md$u + sqrt(2) * s * gh$x
    E <- outer(eta, u, "+")
    g <- colSums(w$y * E - w$m * softplus(E)) - u^2 / (2 * sigma^2)
    a <- lw + g
    amax <- max(a)
    total <- total + log(sqrt(2) * s) - log(sqrt(2 * pi) * sigma) + amax + log(sum(exp(a - amax)))
    if (gradient) {
      om <- exp(a - amax)
      om <- om / sum(om)
      d_eta <- drop((w$y - w$m * stats::plogis(E)) %*% om)
      grad[-K] <- grad[-K] + drop(crossprod(ed$D, d_eta))
      grad[K] <- grad[K] + sum(om * (u^2 / sigma^2 - 1))
    }
  }
  if (gradient) attr(total, "gradient") <- grad
  total
}

# Leave-one-system-out evidence of the validation webs (as validate_rule()):
# for each web, the evidence profiles of its pairs from the registry without
# the studies of its system (and its exclusions), grouped by profile, with
# the number of pairs and of recorded links of each profile. `inferred`: NULL
# or a list, named by web, of the pairs (resource, consumer) inferred by a
# matching rule; the groups are then also split by whether the rule infers
# the pair (`inferred`).
decay_validate <- function(registry, sources, groups = NULL, rule = matching_rule(),
                           train_registries = NULL, exclusions = NULL, inferred = NULL,
                           quiet = TRUE) {
  check_registry(registry)
  left_out <- left_out_sets(groups, sources, registry)
  if (is.null(train_registries)) train_registries <- sort(unique(registry$registry), method = "radix")
  webs <- vector("list", length(sources))
  for (k in seq_along(sources)) {
    w <- sources[k]
    if (!quiet) message(sprintf("Validation web %s (%d/%d)", w, k, length(sources)))
    extra <- if (is.null(exclusions)) character(0) else exclusions$study_id[exclusions$web == w]
    vn <- registry_nodes(registry, w, rule)
    nodes <- vn$nodes
    keep <- !registry$study_id %in% c(left_out[[k]], extra) & registry$registry %in% train_registries
    studied <- nodes$node_id[vn$studied]
    prof <- decay_profiles(nodes, registry[keep, , drop = FALSE], rule, consumers = studied)
    pairs <- expand.grid(resource = nodes$node_id, consumer = studied, stringsAsFactors = FALSE)
    pairs <- pairs[pairs$resource != pairs$consumer, , drop = FALSE]
    rec <- vn$records
    rec <- rec[rec$resource != rec$consumer, , drop = FALSE]
    truth <- unique(paste(rec$resource, rec$consumer, sep = "\u001f"))
    y <- paste(pairs$resource, pairs$consumer, sep = "\u001f") %in% truth
    pp <- decay_pair_profiles(prof, pairs$resource, pairs$consumer)
    inf <- rep(FALSE, nrow(pairs))
    if (!is.null(inferred)) {
      e <- inferred[[w]]
      inf <- paste(pairs$resource, pairs$consumer, sep = "\u001f") %in%
        paste(e$resource, e$consumer, sep = "\u001f")
      pp$key <- paste(pp$key, inf, sep = "|")
    }
    first <- !duplicated(pp$key)
    g <- match(pp$key, pp$key[first])
    webs[[k]] <- list(
      study_id = w, left_out = left_out[[k]],
      P = pp$P[first, , drop = FALSE], G = pp$G[first, , drop = FALSE],
      Q = pp$Q[first, , drop = FALSE], inferred = inf[first],
      m = tabulate(g, sum(first)), y = tabulate(g[y], sum(first)),
      n_nodes = nrow(nodes), n_consumers = length(studied), n_pairs = nrow(pairs),
      n_recorded = length(truth)
    )
  }
  names(webs) <- sources
  structure(list(webs = webs, settings = list(sources = sources, groups = groups,
                                              registries = train_registries,
                                              exclusions = exclusions)),
            class = "bwebs_decay_validation")
}

# Fit of the decay model: maximum marginal likelihood ("ml") or posterior by
# adaptive random-walk Metropolis ("bayes"). Priors: a0 ~ N(0, 1.5^2); with
# "linear", log b ~ N(0, 1); with "steps", log step ~ N(-1, 1.5^2);
# log gamma ~ N(-2, 1.5^2); d_generic ~ N(0, 2.5^2); sigma ~ half-N(0, 2.5^2).
decay_fit <- function(webs, model = c("steps", "linear", "steps_known"), method = c("ml", "bayes"),
                      n_chains = 2L, n_iter = 500L, n_warmup = 500L, quadrature = 20L,
                      seed = NULL, start = NULL) {
  model <- match.arg(model)
  method <- match.arg(method)
  pn <- decay_parameter_names(model)
  K <- length(pn)
  slopes <- if (model == "linear") 2:3 else 2:(2L * length(decay_levels) - 1L)
  if (model == "steps_known") slopes <- c(slopes, K - 1L)
  gh <- gauss_hermite(quadrature)
  lower <- rep(-Inf, K)
  upper <- rep(Inf, K)
  lower[slopes] <- -8
  upper[slopes] <- 8
  lower[K] <- -10
  upper[K] <- 5
  # prior: means and standard deviations of the normal priors of theta (the
  # half-normal prior of sigma is added on the log scale), and start values
  z_names <- pn[-K]
  pm <- ifelse(grepl("_step_", z_names), -1, ifelse(z_names == "log_gamma", -2, 0))
  ps <- ifelse(z_names == "a0", 1.5, ifelse(grepl("^log_b", z_names), 1,
               ifelse(z_names == "d_generic", 2.5, 1.5)))
  lp <- function(theta, gradient = FALSE) {
    if (any(!is.finite(theta)) || any(theta < lower) || any(theta > upper)) return(-Inf)
    ll <- decay_loglik(theta, webs, gh, model, gradient = gradient)
    if (method == "ml") return(ll)
    ls <- theta[[K]]
    z <- theta[-K]
    out <- ll + sum(stats::dnorm(z, pm, ps, log = TRUE)) +
      stats::dnorm(exp(ls), 0, 2.5, log = TRUE) + ls
    if (gradient) {
      attr(out, "gradient") <- attr(ll, "gradient") + c(-(z - pm) / ps^2, 1 - exp(2 * ls) / 2.5^2)
    }
    out
  }
  if (is.null(start)) start <- c(ifelse(z_names == "d_generic", 1, pm), 0)
  start <- stats::setNames(start, pn)
  cache <- new.env()
  fg <- function(t) {
    if (!identical(t, cache$t)) {
      v <- lp(t, gradient = TRUE)
      cache$t <- t
      cache$v <- -as.numeric(v)
      cache$g <- -attr(v, "gradient")
    }
    cache
  }
  opt <- stats::optim(start, function(t) fg(t)$v, function(t) fg(t)$g, method = "L-BFGS-B",
                      lower = lower + 1e-6, upper = upper - 1e-6,
                      control = list(maxit = 5000L, factr = 1e5))
  mode <- opt$par
  if (method == "ml") {
    return(list(model = model, theta = mode, draws = NULL, converged = opt$convergence == 0L,
                loglik = -opt$value))
  }
  H <- stats::optimHess(mode, function(t) fg(t)$v, function(t) fg(t)$g)
  V0 <- tryCatch(solve(H), error = function(e) diag(0.05, K))
  if (any(!is.finite(V0)) || any(diag(V0) <= 0) ||
      inherits(tryCatch(chol(V0), error = function(e) e), "error")) {
    V0 <- diag(0.05, K)
  }
  run <- function() {
    lapply(seq_len(n_chains), function(k) {
      init <- mode + drop(crossprod(chol(V0), stats::rnorm(K))) * 0.5
      metropolis_chain(init, V0, n_iter, n_warmup, lp)
    })
  }
  chains <- if (is.null(seed)) run() else withr::with_seed(seed, run())
  draws <- do.call(rbind, lapply(chains, `[[`, "draws"))
  colnames(draws) <- pn
  arr <- array(NA_real_, c(n_iter, n_chains, K))
  for (k in seq_len(n_chains)) arr[, k, ] <- chains[[k]]$draws
  rhat <- vapply(seq_len(K), function(p) rhat_split(arr[, , p, drop = TRUE]), numeric(1))
  list(model = model, theta = colMeans(draws), draws = draws, mode = mode,
       converged = all(rhat < 1.01, na.rm = TRUE), rhat = stats::setNames(rhat, pn),
       acceptance = vapply(chains, `[[`, numeric(1), "acceptance"))
}

# Probability of each profile in a new system (web effect averaged out),
# averaged over the posterior draws when there are draws.
decay_p_new <- function(fit, P, G, Q = NULL) {
  th <- if (is.null(fit$draws)) matrix(fit$theta, 1L) else fit$draws
  K <- ncol(th)
  p <- vapply(seq_len(nrow(th)), function(d) {
    p_new_system(decay_eta(th[d, ], P, G, fit$model, Q = Q), exp(th[d, K]))
  }, numeric(nrow(P)))
  if (is.matrix(p)) rowMeans(p) else mean(p)
}

# Cross-validated scores of the decay model (as score_rule()): each web is
# predicted by the model fitted to the webs of the other systems.
decay_score <- function(validation, model = "steps", method = "bayes", n_chains = 2L,
                        n_iter = 500L, n_warmup = 500L, n_thin = 10L, seed = NULL, quiet = TRUE) {
  webs <- validation$webs
  ids <- names(webs)
  run <- function() {
    lapply(ids, function(w) {
      out_ids <- webs[[w]]$left_out
      train <- webs[!ids %in% out_ids]
      if (!quiet) message(sprintf("Scoring web %s", w))
      fit <- decay_fit(train, model = model, method = method, n_chains = n_chains,
                       n_iter = n_iter, n_warmup = n_warmup)
      if (!is.null(fit$draws)) {
        fit$draws <- fit$draws[seq(1L, nrow(fit$draws), by = n_thin), , drop = FALSE]
      }
      t <- webs[[w]]
      p <- decay_p_new(fit, t$P, t$G, t$Q)
      brier <- sum(t$y * (1 - p)^2 + (t$m - t$y) * p^2) / t$n_pairs
      p_ref <- mean(vapply(train, function(o) o$n_recorded / o$n_pairs, numeric(1)))
      ref <- (t$n_recorded * (1 - p_ref)^2 + (t$n_pairs - t$n_recorded) * p_ref^2) / t$n_pairs
      list(score = data.frame(study_id = w, n_pairs = t$n_pairs, n_recorded = t$n_recorded,
                              expected = sum(t$m * p), brier = brier, p_reference = p_ref,
                              brier_reference = ref, skill = 1 - brier / ref,
                              stringsAsFactors = FALSE),
           theta = fit$theta)
    })
  }
  rows <- if (is.null(seed)) run() else withr::with_seed(seed, run())
  out <- do.call(rbind, lapply(rows, `[[`, "score"))
  list(webs = out,
       theta = do.call(rbind, lapply(rows, `[[`, "theta")),
       summary = data.frame(n_webs = nrow(out), brier = mean(out$brier),
                            brier_se = stats::sd(out$brier) / sqrt(nrow(out)),
                            skill = mean(out$skill)),
       settings = list(model = model, method = method))
}

# Groups of pairs of a web that the matching rule does not infer.
decay_not_inferred <- function(w) {
  k <- !w$inferred
  w$P <- w$P[k, , drop = FALSE]
  w$G <- w$G[k, , drop = FALSE]
  w$Q <- w$Q[k, , drop = FALSE]
  w$m <- w$m[k]
  w$y <- w$y[k]
  w$inferred <- w$inferred[k]
  w
}

# Cross-validated scores of a matching rule whose pairs without an inferred
# link take the probability of the decay model instead of 0. `reference` is
# the score_rule() result of the rule (its webs give the contribution of the
# inferred links); the decay model is fitted, for each web left out, to the
# pairs not inferred of the webs of the other systems.
decay_score_background <- function(validation, reference, model = "steps", method = "ml",
                                   n_chains = 2L, n_iter = 500L, n_warmup = 500L, n_thin = 10L,
                                   seed = NULL) {
  webs <- lapply(validation$webs, decay_not_inferred)
  ids <- names(webs)
  ref <- reference$webs
  run <- function() {
    lapply(ids, function(w) {
      train <- webs[!ids %in% webs[[w]]$left_out]
      fit <- decay_fit(train, model = model, method = method, n_chains = n_chains,
                       n_iter = n_iter, n_warmup = n_warmup)
      if (!is.null(fit$draws)) {
        fit$draws <- fit$draws[seq(1L, nrow(fit$draws), by = n_thin), , drop = FALSE]
      }
      t <- webs[[w]]
      p <- decay_p_new(fit, t$P, t$G, t$Q)
      r <- ref[ref$study_id == w, ]
      missed <- r$n_recorded - r$n_predicted_observed
      inferred_sum <- r$brier * r$n_pairs - missed
      background_sum <- sum(t$y * (1 - p)^2 + (t$m - t$y) * p^2)
      brier <- (inferred_sum + background_sum) / r$n_pairs
      data.frame(study_id = w, n_pairs = r$n_pairs, n_recorded = r$n_recorded,
                 missed = missed, missed_check = sum(t$y), expected_missed = sum(t$m * p),
                 n_background = sum(t$m), brier_rule = r$brier, brier = brier,
                 brier_reference = r$brier_reference,
                 skill = 1 - brier / r$brier_reference, stringsAsFactors = FALSE)
    })
  }
  rows <- if (is.null(seed)) run() else withr::with_seed(seed, run())
  out <- do.call(rbind, rows)
  d <- out$brier - out$brier_rule
  list(webs = out,
       summary = data.frame(n_webs = nrow(out), brier = mean(out$brier),
                            brier_rule = mean(out$brier_rule), diff = mean(d),
                            diff_se = stats::sd(d) / sqrt(nrow(out)), skill = mean(out$skill)),
       settings = list(model = model, method = method))
}

# Keep only the evidence of the closest consumer level, from `from` on, whose
# consumers have at least `min_known` distinct resources recorded (the
# consumer is resolved at that level, as in the climbing of WebBuilder); the
# evidence of the other levels is dropped.
decay_resolve <- function(P, G, Q, min_known = 1, from = "family") {
  L <- length(decay_levels)
  K <- cbind(round(expm1(Q[, -1L, drop = FALSE])), Inf)
  ok <- K >= min_known
  ok[, seq_len(match(from, decay_levels) - 1L)] <- FALSE
  c0 <- max.col(ok * 1, ties.method = "first")
  keep <- outer(c0, seq_len(L), "==")
  P[!keep] <- NA_integer_
  G[!keep] <- FALSE
  list(P = P, G = G)
}
