# Taxonomic distance decay of link probabilities --------------------------------
#
# With matching_rule(decay = TRUE), the resource -> consumer pairs that the
# matching rule does not infer get a probability from the taxonomic distance
# decay instead of 0: the closer the records of the relatives of the consumer
# come to the resource, the more probable the link.
#
# Evidence of a pair. decay_profiles() gives, for each level of relatedness
# on the consumer side (the consumer itself, genus, subfamily, family, order,
# class, phylum or any consumer), the closest relation on the resource side
# among the records of those consumers. decay_evidence() resolves the
# consumer at the closest of those levels, from the family on, whose
# consumers have at least one resource recorded (as the climbing of
# WebBuilder), and keeps the closest relation of the resource at that level
# and whether a record of a generic group (such as "Algae") applies to the
# resource. The model of the probability of each evidence profile is
# fit_decay().

decay_levels <- c("exact", "genus", "subfamily", "family", "order", "class", "phylum", "any")
# Levels at which the consumer is resolved; "family" includes the closer
# relatives (the consumer itself, its genus and subfamily).
decay_consumer_levels <- c("family", "order", "class", "phylum", "any")
# rank_cols (R/rules.R) and the phylum
decay_ranks <- c("genus", "subfamily", "family", "order", "class", "phylum")

decay_low <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x[!is.na(x) & x == ""] <- NA_character_
  x
}

# Phylum of each class (internal table class_phylum, data-raw/class_phylum.R).
phylum_of_class <- function(class) {
  class_phylum$phylum[match(class, class_phylum$class)]
}

# Ranks of the nodes in lower case, with the phylum of their class (or their
# own phylum when the class is unknown), names by role and life stage.
decay_node_info <- function(nodes) {
  n <- nrow(nodes)
  ranks <- lapply(stats::setNames(rank_cols, rank_cols), function(r) decay_low(nodes[[r]]))
  ph <- phylum_of_class(nodes$class)
  if ("phylum" %in% names(nodes)) ph[is.na(ph)] <- nodes$phylum[is.na(ph)]
  ranks$phylum <- decay_low(ph)
  col <- function(nm) if (nm %in% names(nodes)) nodes[[nm]] else rep(NA_character_, n)
  alias <- col("alias")
  names_of <- function(role) {
    sur <- col(paste0(role, "_surrogate"))
    lapply(seq_len(n), function(k) {
      unique(decay_low(c(nodes$taxon[k], split_alias(alias[k]), split_alias(sur[k]))))
    })
  }
  surrogates_of <- function(role) {
    sur <- col(paste0(role, "_surrogate"))
    lapply(seq_len(n), function(k) unique(decay_low(split_alias(sur[k]))))
  }
  stage <- col("life_stage")
  stage[is.na(stage) | !stage %in% c("larva", "adult")] <- ""
  list(n = n, ranks = ranks, stage = stage,
       res_names = names_of("resource"), con_names = names_of("consumer"),
       res_sur = surrogates_of("resource"), con_sur = surrogates_of("consumer"))
}

# Registry columns in lower case (with the phylum of each class) and, for
# the consumer side, the rows of each name and of each taxon of each rank.
decay_index <- function(registry) {
  n <- nrow(registry)
  rx <- list(n = n)
  rows <- seq_len(n)
  hash <- function(v) {
    ok <- !is.na(v)
    list2env(split(rows[ok], v[ok]), hash = TRUE)
  }
  for (side in c("res", "con")) {
    nm <- if (side == "res") "resource" else "consumer"
    rx[[paste0(side, "_name")]] <- decay_low(registry[[nm]])
    rx[[paste0(side, "_std")]] <- decay_low(registry[[paste0(nm, "_std")]])
    for (r in rank_cols) rx[[paste0(side, "_", r)]] <- decay_low(registry[[paste0(side, "_", r)]])
    rx[[paste0(side, "_phylum")]] <- decay_low(phylum_of_class(registry[[paste0(side, "_class")]]))
    st <- registry[[paste0(side, "_stage")]]
    st[is.na(st) | !st %in% c("larva", "adult")] <- ""
    rx[[paste0(side, "_stage")]] <- st
  }
  rx$con_rows <- c(list(name = hash(rx$con_name), std = hash(rx$con_std)),
                   lapply(stats::setNames(decay_ranks, decay_ranks),
                          function(r) hash(rx[[paste0("con_", r)]])))
  rx
}

# Rows of the records whose consumer has `value` at `rank` (a name with
# rank = "name"), as an integer vector.
decay_rows <- function(rx, rank, value) {
  if (is.na(value)) return(integer(0))
  out <- if (rank == "name") {
    c(rx$con_rows$name[[value]], rx$con_rows$std[[value]])
  } else {
    rx$con_rows[[rank]][[value]]
  }
  if (is.null(out)) integer(0) else unique(out)
}

# Generic resource names that apply to each node (groups of the rule).
decay_generic <- function(nodes, info, rule) {
  g <- rule$generic_groups
  out <- replicate(info$n, character(0), simplify = FALSE)
  for (k in seq_len(nrow(g))) {
    rk <- g$rank[k]
    v <- if (rk %in% names(info$ranks)) info$ranks[[rk]] else rep(NA_character_, info$n)
    hit <- !is.na(v) & v %in% decay_low(split_alias(g$values[k]))
    for (i in which(hit)) out[[i]] <- c(out[[i]], decay_low(g$generic[k]))
  }
  out
}

# For the records `rows` (already filtered by the life stage of the
# consumer): the closest level (index of decay_levels) at which each node is
# related to the resource of some record, and whether a record of a generic
# group that applies to the node is among them. NA when there is no record.
decay_resource_levels <- function(rows, rx, info, generic, generic_names) {
  n <- info$n
  level <- rep(NA_integer_, n)
  gen <- rep(FALSE, n)
  if (!length(rows)) return(list(level = level, generic = gen, names = character(0)))
  for (st in unique(info$stage)) {
    nodes_st <- which(info$stage == st)
    r <- if (st == "") rows else rows[rx$res_stage[rows] %in% c("", st)]
    if (!length(r)) next
    names_set <- unique(c(rx$res_name[r], rx$res_std[r]))
    lev <- rep(length(decay_levels), length(nodes_st))
    for (k in rev(seq_along(decay_ranks))) {
      vals <- info$ranks[[decay_ranks[k]]][nodes_st]
      set <- unique(rx[[paste0("res_", decay_ranks[k])]][r])
      lev[!is.na(vals) & vals %in% set] <- k + 1L
    }
    sur_set <- if (any(lengths(info$res_sur[nodes_st]) > 0L)) {
      unique(unlist(lapply(rank_cols, function(rk) rx[[paste0("res_", rk)]][r]), use.names = FALSE))
    } else {
      character(0)
    }
    # hash the large sets once: keep only the names of these nodes found in them
    all_names <- unique(unlist(info$res_names[nodes_st], use.names = FALSE))
    hit_names <- all_names[all_names %in% names_set]
    all_sur <- unique(unlist(info$res_sur[nodes_st], use.names = FALSE))
    hit_sur <- all_sur[all_sur %in% sur_set]
    exact <- vapply(nodes_st, function(i) {
      any(info$res_names[[i]] %in% hit_names) || any(info$res_sur[[i]] %in% hit_sur)
    }, logical(1))
    lev[exact] <- 1L
    level[nodes_st] <- lev
    eaten_generic <- intersect(names_set, generic_names)
    if (length(eaten_generic)) {
      gen[nodes_st] <- vapply(nodes_st, function(i) any(generic[[i]] %in% eaten_generic),
                              logical(1))
    }
  }
  list(level = level, generic = gen, names = unique(rx$res_std[rows]))
}

# Evidence profiles of every resource -> consumer pair.
#
# Returns a list with `level`, an integer array resource x consumer x
# consumer level (the 8 decay_levels): the closest level of relation on the
# resource side among the records whose consumer is related to the consumer
# of the pair at that level or closer (NA: no such record); and `generic`, a
# logical array of the same dimensions (a record of a generic group that
# applies to the resource); and `known`, a matrix consumer x consumer level
# with the number of distinct resources recorded for the consumers related to
# each consumer at that level or closer (how well its diet is known). Only the
# consumers in `consumers` are computed.
decay_profiles <- function(nodes, registry, rule = attr(nodes, "rule"),
                           consumers = nodes$node_id) {
  info <- decay_node_info(nodes)
  rx <- decay_index(registry)
  generic <- decay_generic(nodes, info, rule)
  generic_names <- decay_low(rule$generic_groups$generic)
  J <- match(consumers, nodes$node_id)
  L <- length(decay_levels)
  level <- array(NA_integer_, c(info$n, length(J), L),
                 dimnames = list(resource = nodes$node_id, consumer = consumers, level = decay_levels))
  gen <- array(FALSE, dim(level), dimnames = dimnames(level))
  known <- matrix(0L, length(J), L, dimnames = list(consumer = consumers, level = decay_levels))
  cache <- new.env(parent = emptyenv(), hash = TRUE)
  all_rows <- seq_len(rx$n)
  for (jj in seq_along(J)) {
    j <- J[jj]
    st <- info$stage[j]
    stage_ok <- function(r) if (st == "") r else r[rx$con_stage[r] %in% c("", st)]
    cur_level <- rep(NA_integer_, info$n)
    cur_gen <- rep(FALSE, info$n)
    cur_names <- character(0)
    for (c in seq_len(L)) {
      lev_name <- decay_levels[c]
      key <- if (lev_name == "exact") {
        paste("exact", paste(info$con_names[[j]], collapse = "\u0003"),
              paste(info$con_sur[[j]], collapse = "\u0003"), st, sep = "\u0001")
      } else if (lev_name == "any") {
        paste("any", st, sep = "\u0001")
      } else {
        paste(lev_name, info$ranks[[lev_name]][j], st, sep = "\u0001")
      }
      v <- cache[[key]]
      if (is.null(v)) {
        rows <- if (lev_name == "exact") {
          r <- unlist(lapply(info$con_names[[j]], decay_rows, rx = rx, rank = "name"),
                      use.names = FALSE)
          for (s in info$con_sur[[j]]) {
            for (rk in decay_ranks) r <- c(r, decay_rows(rx, rk, s))
          }
          unique(r)
        } else if (lev_name == "any") {
          all_rows
        } else {
          decay_rows(rx, lev_name, info$ranks[[lev_name]][j])
        }
        v <- decay_resource_levels(stage_ok(rows), rx, info, generic, generic_names)
        assign(key, v, envir = cache)
      }
      cur_level <- pmin(cur_level, v$level, na.rm = TRUE)
      cur_gen <- cur_gen | v$generic
      level[, jj, c] <- cur_level
      gen[, jj, c] <- cur_gen
      if (lev_name != "any") {
        cur_names <- union(cur_names, v$names)
        known[jj, c] <- length(cur_names)
      } else {
        known[jj, c] <- NA_integer_
      }
    }
  }
  list(level = level, generic = gen, known = known)
}

# Evidence of the pairs (resource, consumer) from their profiles: the level
# at which the consumer is resolved (decay_consumer_levels), the closest
# level at which the resource is related to a resource of the records of
# those consumers (decay_levels; "any" when the records are only of other
# life stages of the resource) and whether a record of a generic group
# applies to the resource.
decay_evidence <- function(profiles, resource, consumer) {
  if (!length(resource)) {
    return(data.frame(consumer_level = character(), resource_level = character(),
                      generic = logical(), stringsAsFactors = FALSE))
  }
  L <- length(decay_levels)
  first <- match("family", decay_levels)
  ri <- match(resource, dimnames(profiles$level)$resource)
  ci <- match(consumer, dimnames(profiles$level)$consumer)
  kn <- profiles$known[ci, first:(L - 1L), drop = FALSE]
  ok <- !is.na(kn) & kn >= 1L
  cl <- rep(length(decay_consumer_levels), length(ci))
  found <- rowSums(ok) > 0L
  cl[found] <- max.col(ok[found, , drop = FALSE] * 1, ties.method = "first")
  lev <- first - 1L + cl
  rl <- profiles$level[cbind(ri, ci, lev)]
  rl[is.na(rl)] <- L
  data.frame(consumer_level = decay_consumer_levels[cl], resource_level = decay_levels[rl],
             generic = profiles$generic[cbind(ri, ci, lev)], stringsAsFactors = FALSE)
}

# Every evidence profile (consumer level x resource level x generic record),
# in the order of decay_profile_id().
decay_all_profiles <- function() {
  data.frame(consumer_level = rep(decay_consumer_levels, each = 16L),
             resource_level = rep(rep(decay_levels, each = 2L), 5L),
             generic = rep(c(FALSE, TRUE), 40L), stringsAsFactors = FALSE)
}

decay_profile_id <- function(consumer_level, resource_level, generic) {
  (match(consumer_level, decay_consumer_levels) - 1L) * 16L +
    (match(resource_level, decay_levels) - 1L) * 2L + as.integer(generic) + 1L
}

# Pairs and recorded links of each evidence profile of a web.
decay_count_table <- function(study, evidence, observed) {
  id <- decay_profile_id(evidence$consumer_level, evidence$resource_level, evidence$generic)
  n <- tabulate(id, 80L)
  y <- tabulate(id[observed], 80L)
  k <- which(n > 0L)
  out <- cbind(data.frame(study_id = rep(study, length(k)), stringsAsFactors = FALSE),
               decay_all_profiles()[k, , drop = FALSE],
               data.frame(pairs = n[k], observed = y[k]))
  rownames(out) <- NULL
  out
}

# Model of the decay ----------------------------------------------------------
#
# logit p = a[consumer level] + b[resource level] + gamma [generic record] +
# v[web]: b is relative to a resource not related to the recorded diet
# ("any"). The resource levels within the family of a recorded resource (the
# resource itself, its genus, subfamily or family) share one effect: the
# rule infers almost every such pair.

# Groups of the resource levels in the model.
decay_resource_groups <- c("family", "order", "class", "phylum", "any")
decay_resource_group <- function(rl) c(1L, 1L, 1L, 1L, 2L, 3L, 4L, 5L)[rl]

# Parameters: the logit of a pair whose consumer is resolved at its family
# and whose resource is not related to its recorded diet (in a web of 50
# nodes); the steps down for consumers resolved at their order, class,
# phylum or as any consumer; the steps up for resources related at the
# phylum, class, order and within the family (all steps >= 0); the effect of
# a record of a generic group that applies to the resource (added); with
# web_size = TRUE, the effect of the log of the number of nodes of the web
# relative to 50 (decay_web_size()); and the log standard deviation of the
# web effect.
decay_names <- function(web_size = FALSE) {
  c("intercept", paste0("step_consumer_", decay_consumer_levels[-1L]),
    paste0("step_resource_", rev(decay_resource_groups[-5L])), "generic",
    if (web_size) "web_size", "log_sigma")
}

# Covariate of the size of a web: log(number of nodes / 50).
decay_web_size <- function(n_nodes) log(n_nodes / 50)

# Effects of the consumer levels (with a resource not related, "any") and of
# the resource levels (added; 0 for "any") on the logit scale, and of a
# generic record (added).
decay_level_effects <- function(theta) {
  list(consumer = theta[[1L]] - cumsum(c(0, theta[2:5])),
       resource = rev(cumsum(c(0, theta[6:9]))),
       generic = theta[[10L]])
}

# Linear predictor of evidence profiles given by the index of their consumer
# level (cl, in decay_consumer_levels), of their resource group (rg, in
# decay_resource_groups) and whether a generic record applies (g), in a web
# whose size covariate is z (used when theta has the web_size parameter).
# With design = TRUE, also its derivatives with respect to theta (without
# log_sigma).
decay_linear <- function(theta, cl, rg, g, z = 0, design = FALSE) {
  ef <- decay_level_effects(theta)
  with_size <- length(theta) == 12L
  eta <- ef$consumer[cl] + ef$resource[rg] + ef$generic * g +
    if (with_size) theta[[11L]] * z else 0
  if (!design) return(eta)
  D <- matrix(0, length(eta), length(theta) - 1L)
  D[, 1L] <- 1
  for (k in 2:5) D[, k] <- -as.numeric(cl >= k)
  for (k in 1:4) D[, k + 5L] <- as.numeric(rg <= 5L - k)
  D[, 10L] <- as.numeric(g)
  if (with_size) D[, 11L] <- z
  list(eta = eta, D = D)
}

# Derivatives of the effects (consumer levels, resource groups, generic)
# with respect to theta, for the delta method.
decay_effect_jacobian <- function(theta) {
  with_size <- length(theta) == 12L
  J <- matrix(0, 11L + with_size, length(theta))
  J[1:5, 1L] <- 1
  for (c in 2:5) J[c, 2:c] <- -1
  for (r in 1:4) J[5L + r, 6:(10L - r)] <- 1
  J[11L, 10L] <- 1
  if (with_size) J[12L, 11L] <- 1
  J
}

# Evidence profiles of each web of a table of decay counts, with the size
# covariate of the web (0 when the counts have no `web_nodes`).
decay_webs <- function(counts) {
  ids <- sort(unique(counts$study_id), method = "radix")
  out <- lapply(ids, function(w) {
    d <- counts[counts$study_id == w, , drop = FALSE]
    list(cl = match(d$consumer_level, decay_consumer_levels),
         rg = decay_resource_group(match(d$resource_level, decay_levels)), g = d$generic,
         y = d$observed, m = d$pairs,
         z = if (is.null(d$web_nodes)) 0 else decay_web_size(d$web_nodes[1L]))
  })
  names(out) <- ids
  out
}

# Marginal log likelihood of the decay model (up to the binomial
# coefficients), with the web effects integrated by adaptive Gauss-Hermite
# quadrature as in marginal_loglik(), and its gradient with respect to theta.
decay_marginal <- function(theta, webs, gh, gradient = FALSE) {
  K <- length(theta)
  sigma <- exp(theta[[K]])
  lw <- log(gh$w) + gh$x^2
  total <- 0
  grad <- numeric(K)
  for (w in webs) {
    ed <- decay_linear(theta, w$cl, w$rg, w$g, w$z, design = gradient)
    eta <- if (gradient) ed$eta else ed
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

# Fit of the decay model to the counts of validate_rule() (its `decay`
# table; with a column `web_nodes`, the model has the effect of the size of
# the web): maximum marginal likelihood ("ml") or posterior by the adaptive
# random-walk Metropolis of calibrate_links() ("bayes"). Priors: intercept,
# generic and web size ~ N(0, prior_sd^2), steps ~ half-N(0, prior_sd^2), sigma ~
# half-N(0, prior_sigma^2). At most `max_draws` posterior draws are kept,
# evenly spaced.
fit_decay <- function(counts, method = "bayes", prior_sd = 1.5, prior_sigma = 2.5,
                      n_chains = 4L, n_iter = 2000L, n_warmup = 1000L, quadrature = 20L,
                      seed = NULL, max_draws = 4000L, thin = 20L) {
  if (!nrow(counts) || sum(counts$pairs) == 0) {
    abort_input("No pair without an inferred link in the validation webs: nothing to fit.")
  }
  webs <- decay_webs(counts)
  with_size <- !is.null(counts$web_nodes)
  pn <- decay_names(with_size)
  K <- length(pn)
  steps <- 2:9
  gh <- gauss_hermite(quadrature)
  lower <- c(-30, rep(0, 8L), -30, if (with_size) -30, -10)
  upper <- c(30, rep(20, 8L), 30, if (with_size) 30, 5)
  lp <- function(theta, gradient = FALSE) {
    if (any(!is.finite(theta)) || any(theta < lower) || any(theta > upper)) return(-Inf)
    ll <- decay_marginal(theta, webs, gh, gradient = gradient)
    if (method == "ml") return(ll)
    ls <- theta[[K]]
    z <- theta[-K]
    out <- ll - sum(z^2) / (2 * prior_sd^2) + stats::dnorm(exp(ls), 0, prior_sigma, log = TRUE) + ls
    if (gradient) {
      attr(out, "gradient") <- attr(ll, "gradient") +
        c(-z / prior_sd^2, 1 - exp(2 * ls) / prior_sigma^2)
    }
    out
  }
  rate <- sum(counts$observed) / sum(counts$pairs)
  start <- stats::setNames(c(stats::qlogis(min(max(rate, 1e-3), 0.5)), rep(0.5, 8L), 0,
                             if (with_size) 0, 0), pn)
  cache <- new.env(parent = emptyenv())
  fg <- function(t) {
    t <- pmin(pmax(t, lower), upper)
    if (!identical(t, cache$t)) {
      v <- lp(t, gradient = TRUE)
      cache$t <- t
      cache$v <- -as.numeric(v)
      cache$g <- if (is.finite(cache$v)) -attr(v, "gradient") else rep(0, K)
    }
    cache
  }
  opt <- stats::optim(start, function(t) fg(t)$v, function(t) fg(t)$g, method = "L-BFGS-B",
                      lower = lower, upper = upper, control = list(maxit = 5000L, factr = 1e5))
  mode <- stats::setNames(opt$par, pn)
  at_bound <- mode <= lower + 1e-6 | mode >= upper - 1e-6
  hess <- function(t) {
    H <- tryCatch(stats::optimHess(t, function(x) fg(x)$v, function(x) fg(x)$g),
                  error = function(e) NULL)
    ok <- !is.null(H) && all(is.finite(H)) &&
      !inherits(tryCatch(chol(H), error = function(e) e), "error")
    if (ok) H else NULL
  }
  rhat <- ess <- acceptance <- NULL
  if (method == "ml") {
    H <- if (any(at_bound)) NULL else hess(mode)
    info_pd <- !is.null(H)
    theta <- mode
    V <- if (info_pd) solve(H) else matrix(NA_real_, K, K)
    dimnames(V) <- list(pn, pn)
    draws <- NULL
    converged <- opt$convergence == 0L
  } else {
    inside <- mode
    inside[steps] <- pmax(inside[steps], 0.05)
    H <- hess(inside)
    info_pd <- !is.null(H)
    V0 <- if (info_pd) solve(H) else diag(0.05, K)
    if (any(!is.finite(V0)) || inherits(tryCatch(chol(V0), error = function(e) e), "error")) {
      V0 <- diag(0.05, K)
    }
    run <- function() {
      lapply(seq_len(n_chains), function(k) {
        init <- inside + drop(crossprod(chol(V0), stats::rnorm(K))) * 0.5
        init[steps] <- abs(init[steps])
        init <- pmin(pmax(init, lower), upper)
        metropolis_reflect(init, V0, n_iter, n_warmup, lp, reflect = steps, thin = thin)
      })
    }
    chains <- if (is.null(seed)) run() else withr::with_seed(seed, run())
    arr <- array(NA_real_, c(n_iter, n_chains, K), dimnames = list(NULL, NULL, pn))
    for (k in seq_len(n_chains)) arr[, k, ] <- chains[[k]]$draws
    draws <- matrix(arr, n_iter * n_chains, K, dimnames = list(NULL, pn))
    rhat <- vapply(pn, function(p) rhat_split(arr[, , p, drop = TRUE]), numeric(1))
    ess <- vapply(pn, function(p) ess_bulk(arr[, , p, drop = TRUE]), numeric(1))
    acceptance <- vapply(chains, `[[`, numeric(1), "acceptance")
    converged <- all(rhat < 1.01, na.rm = TRUE) && all(ess > 400, na.rm = TRUE)
    keep <- unique(round(seq(1, nrow(draws), length.out = min(max_draws, nrow(draws)))))
    draws <- draws[keep, , drop = FALSE]
    theta <- colMeans(draws)
    V <- stats::cov(draws)
  }
  th <- if (is.null(draws)) matrix(theta, 1L, K, dimnames = list(NULL, pn)) else draws
  # effects of the levels
  eff <- t(apply(th, 1L, function(t) {
    c(unlist(decay_level_effects(t), use.names = FALSE), if (with_size) t[[11L]])
  }))
  if (is.null(draws)) {
    est <- eff[1L, ]
    J <- decay_effect_jacobian(theta)
    se <- if (info_pd) sqrt(pmax(diag(J %*% V %*% t(J)), 0)) else rep(NA_real_, nrow(J))
    lo <- est - 1.959964 * se
    hi <- est + 1.959964 * se
  } else {
    est <- colMeans(eff)
    se <- apply(eff, 2L, stats::sd)
    q <- apply(eff, 2L, stats::quantile, probs = c(0.025, 0.975), names = FALSE)
    lo <- q[1L, ]
    hi <- q[2L, ]
  }
  effects <- data.frame(side = c(rep("consumer", 5L), rep("resource", 5L), "generic",
                                 if (with_size) "web"),
                        level = c(decay_consumer_levels, decay_resource_groups, "generic",
                                  if (with_size) "log(nodes / 50)"),
                        estimate = unname(est), std_error = unname(se), lower = unname(lo),
                        upper = unname(hi), stringsAsFactors = FALSE)
  prof <- decay_probability_table(th, counts)
  sigma <- exp(theta[[K]])
  web_effect <- vapply(webs, function(w) {
    web_mode(decay_linear(theta, w$cl, w$rg, w$g, w$z), sigma, w$y, w$m)$u
  }, numeric(1))
  list(effects = effects, probabilities = prof,
       sigma = if (is.null(draws)) sigma else mean(exp(draws[, K])), theta = theta, vcov = V,
       draws = draws, web_effect = web_effect, converged = converged, information_pd = info_pd,
       message = if (is.null(opt$message)) NA_character_ else opt$message, rhat = rhat,
       ess = ess, acceptance = acceptance)
}

# Probability of every evidence profile (decay_all_profiles()) with a web
# effect of 0 (`p_average`) and in a new system (`p_new`), in a web of 50
# nodes, averaged over the rows of `th` (draws, or a single row of
# estimates), and the number of webs, pairs and recorded links behind each
# profile in `counts`.
decay_probability_table <- function(th, counts = NULL) {
  if (is.null(dim(th))) th <- matrix(th, nrow = 1L)
  K <- ncol(th)
  prof <- decay_all_profiles()
  cl <- match(prof$consumer_level, decay_consumer_levels)
  rg <- decay_resource_group(match(prof$resource_level, decay_levels))
  E <- vapply(seq_len(nrow(th)), function(d) decay_linear(th[d, ], cl, rg, prof$generic),
              numeric(nrow(prof)))
  E <- matrix(E, nrow(prof))
  prof$p_average <- rowMeans(stats::plogis(E))
  prof$p_new <- rowMeans(p_new_matrix(E, exp(th[, K])))
  prof$webs <- 0L
  prof$pairs <- 0
  prof$observed <- 0
  if (!is.null(counts) && nrow(counts)) {
    id <- decay_profile_id(counts$consumer_level, counts$resource_level, counts$generic)
    prof$webs <- tabulate(id, nrow(prof))
    s <- rowsum(cbind(counts$pairs, counts$observed), id)
    prof$pairs[as.integer(rownames(s))] <- s[, 1L]
    prof$observed[as.integer(rownames(s))] <- s[, 2L]
  }
  prof
}

# Adaptive random-walk Metropolis for the decay model. The proposals of the
# parameters bounded below by 0 (`reflect`, the steps) are reflected at 0,
# which keeps the proposal symmetric, so that parameters at the bound do not
# make every proposal fail. Warm-up as in metropolis_chain() (scale adapted
# towards an acceptance rate of 0.234, covariance of the chain so far from
# iteration 200, every 100 iterations); afterwards the proposal is fixed.
# Each iteration makes `thin` steps.
metropolis_reflect <- function(init, V0, n_iter, n_warmup, lp, reflect, thin = 1L) {
  d <- length(init)
  S <- V0
  log_scale <- log(2.38^2 / d)
  R <- chol(exp(log_scale) * S)
  propose <- function(x) {
    prop <- x + drop(crossprod(R, stats::rnorm(d)))
    prop[reflect] <- abs(prop[reflect])
    prop
  }
  cur <- init
  cur_lp <- lp(cur)
  tries <- 0L
  while (!is.finite(cur_lp) && tries < 100L) {
    cur <- propose(init)
    cur_lp <- lp(cur)
    tries <- tries + 1L
  }
  hist <- matrix(NA_real_, max(n_warmup, 1L), d)
  n_steps <- 0L
  for (t in seq_len(n_warmup)) {
    for (j in seq_len(thin)) {
      prop <- propose(cur)
      prop_lp <- lp(prop)
      ok <- is.finite(prop_lp) && log(stats::runif(1)) < prop_lp - cur_lp
      if (ok) {
        cur <- prop
        cur_lp <- prop_lp
      }
      n_steps <- n_steps + 1L
      log_scale <- log_scale + (as.numeric(ok) - 0.234) / sqrt(n_steps)
    }
    hist[t, ] <- cur
    if (t >= 200L && t %% 100L == 0L) {
      Sn <- stats::cov(hist[ceiling(t / 2):t, , drop = FALSE]) + diag(1e-8, d)
      if (all(is.finite(Sn))) S <- Sn
    }
    Rn <- tryCatch(chol(exp(log_scale) * S), error = function(e) NULL)
    if (!is.null(Rn)) R <- Rn
  }
  out <- matrix(NA_real_, n_iter, d)
  accepted <- 0L
  for (t in seq_len(n_iter)) {
    for (j in seq_len(thin)) {
      prop <- propose(cur)
      prop_lp <- lp(prop)
      if (is.finite(prop_lp) && log(stats::runif(1)) < prop_lp - cur_lp) {
        cur <- prop
        cur_lp <- prop_lp
        accepted <- accepted + 1L
      }
    }
    out[t, ] <- cur
  }
  list(draws = out, acceptance = accepted / (n_iter * thin))
}

# Probability in a new system (web effects averaged out, as p_new_system())
# of a matrix of linear predictors (profiles x draws), with the standard
# deviation of the web effect of each draw.
p_new_matrix <- function(E, sigma, n = 30L) {
  gh <- gauss_hermite(n)
  S <- matrix(sigma, nrow(E), ncol(E), byrow = TRUE)
  P <- 0
  for (k in seq_along(gh$x)) P <- P + gh$w[k] * stats::plogis(E + sqrt(2) * gh$x[k] * S)
  P / sqrt(pi)
}

# Probability of the profiles of a table with columns consumer_level,
# resource_level and generic: in a new system (system = "random", as
# decay_p_new_of()), or with a system effect of 0 ("average") or of a given
# value, averaged over the posterior draws, in a web of `n_nodes` nodes.
decay_probability_of <- function(fit, x, n_nodes = 50, system = "random") {
  if (identical(system, "random")) return(decay_p_new_of(fit, x, n_nodes))
  id <- decay_profile_id(x$consumer_level, x$resource_level, x$generic)
  if (!length(id)) return(numeric(0))
  th <- if (is.null(fit$draws)) matrix(fit$theta, nrow = 1L) else fit$draws
  u <- sort(unique(id))
  ap <- decay_all_profiles()[u, , drop = FALSE]
  cl <- match(ap$consumer_level, decay_consumer_levels)
  rg <- decay_resource_group(match(ap$resource_level, decay_levels))
  z <- decay_web_size(n_nodes)
  s <- if (identical(system, "average")) 0 else system
  E <- matrix(vapply(seq_len(nrow(th)), function(d) {
    decay_linear(th[d, ], cl, rg, ap$generic, z)
  }, numeric(length(u))), length(u))
  rowMeans(stats::plogis(E + s))[match(id, u)]
}

# Probability in a new system of the profiles of a table with columns
# consumer_level, resource_level and generic, from a fit_decay() result (or
# the `decay` element of a calibration), in a web of `n_nodes` nodes (only
# used by a model with the effect of the size of the web).
decay_p_new_of <- function(fit, x, n_nodes = 50) {
  id <- decay_profile_id(x$consumer_level, x$resource_level, x$generic)
  th <- if (is.null(fit$draws)) matrix(fit$theta, nrow = 1L) else fit$draws
  if (ncol(th) < 12L || !length(id)) return(fit$probabilities$p_new[id])
  u <- sort(unique(id))
  ap <- decay_all_profiles()[u, , drop = FALSE]
  cl <- match(ap$consumer_level, decay_consumer_levels)
  rg <- decay_resource_group(match(ap$resource_level, decay_levels))
  z <- decay_web_size(n_nodes)
  E <- matrix(vapply(seq_len(nrow(th)), function(d) {
    decay_linear(th[d, ], cl, rg, ap$generic, z)
  }, numeric(length(u))), length(u))
  rowMeans(p_new_matrix(E, exp(th[, ncol(th)])))[match(id, u)]
}
