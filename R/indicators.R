# Topological indicators of one food web --------------------------------------

#' Topological indicators of a food web
#'
#' `web_indicators()` computes the structural indicators of one binary food
#' web given as an adjacency matrix with resources in rows and consumers in
#' columns. Cannibalistic links are removed, isolated nodes are left out (and
#' counted), and nodes with identical prey and predators are aggregated into
#' trophic species before computing the indicators.
#'
#' @param adjacency A square numeric or logical matrix with identical row and
#'   column names (node identifiers); `adjacency[i, j]` is non-zero when
#'   consumer `j` eats resource `i`.
#' @param permanent Identifiers of permanent resources (for example detritus
#'   and macrophytes) that cannot go extinct in the robustness analysis.
#' @param modularity If `TRUE`, compute the modularity of the best of
#'   `louvain_runs` partitions found with the Louvain algorithm on the
#'   undirected web (requires the igraph package).
#' @param robustness If `TRUE`, compute R50: the proportion of primary
#'   removals needed to lose at least half of the removable trophic species,
#'   with removals in random order (mean of `robustness_reps` sequences) and
#'   from the most to the least connected; a consumer goes extinct when it
#'   loses all its prey.
#' @param louvain_runs,robustness_reps Number of Louvain runs and of random
#'   removal sequences (positive whole numbers).
#' @param seed `NULL` (default) to use the current random number stream for
#'   the stochastic indicators (modularity and random R50), or a number:
#'   each of them is then computed after setting this seed, and the global
#'   random number stream is restored afterwards.
#'
#' @return A data frame with one row and columns:
#' * `S_taxa`, `n_isolated`: nodes given and nodes without links;
#' * `S`, `L`, `LS`, `C`: trophic species, links, links per species
#'   (\eqn{L/S}) and connectance (\eqn{L/S^2});
#' * `prop_basal`, `prop_intermediate`, `prop_top`: proportions of trophic
#'   species without prey, with prey and predators, and with prey but without
#'   predators;
#' * `gen_mean`, `vul_mean`: mean number of prey of consumers and of
#'   predators of resources; `gen_sd`, `vul_sd`: standard deviations of the
#'   number of prey and of predators of all trophic species divided by
#'   \eqn{L/S};
#' * `TL_mean`, `TL_max`: mean and maximum prey-averaged trophic level;
#'   `TL_short_mean`: mean shortest-chain trophic level (of the species
#'   connected to a basal species); `TL_status`: `"solve"`, `"ginv"`
#'   (solved with a generalised inverse) or `"unresolved"`;
#' * `prop_omnivores`, `omnivory_mean`: proportion of consumers whose prey
#'   have trophic levels with standard deviation above 0.01, and mean of that
#'   standard deviation;
#' * `MxSim_taxa`, `MxSim_trophic_species`: mean maximum trophic similarity
#'   among taxa (before aggregation) and among trophic species;
#' * with `modularity = TRUE`: `modularity`, `n_modules`; with
#'   `robustness = TRUE`: `R50_random`, `R50_connected`, `R50_random_living`
#'   and `R50_connected_living` (the last two on the web without the
#'   permanent resources);
#' * `status`: `"ok"` or `"empty_web"` (fewer than two linked nodes); in an
#'   empty web every indicator except `S_taxa` and `n_isolated` is `NA`.
#'
#' @details
#' Trophic species: nodes with the same set of prey and the same set of
#' predators. Prey-averaged trophic level: \eqn{TL_j = 1 + \frac{1}{n_j}
#' \sum_{i \in prey(j)} TL_i}{TL[j] = 1 + mean(TL of the prey of j)}, with
#' \eqn{TL = 1} for species without prey, solved as a linear system; when it
#' is singular (for example a loop without input from a basal species) a
#' generalised inverse is tried and accepted only if it satisfies the
#' equations. Trophic similarity of two species: number of shared prey and
#' predators divided by the number of prey and predators of either.
#'
#' @section Missing values and zeros:
#' The adjacency matrix must not contain `NA`. Indicators that are undefined
#' (for example the mean number of predators when no species has predators)
#' are `NA`.
#'
#' @seealso [calc_indicators()] to compute them for every sample and draw.
#' @export
#' @examples
#' # A three-level chain: alga -> grazer -> predator
#' a <- matrix(0, 3, 3, dimnames = list(c("alga", "grazer", "predator"),
#'                                      c("alga", "grazer", "predator")))
#' a["alga", "grazer"] <- 1
#' a["grazer", "predator"] <- 1
#' web_indicators(a)
web_indicators <- function(adjacency, permanent = character(0), modularity = FALSE,
                           robustness = FALSE, louvain_runs = 10L,
                           robustness_reps = 100L, seed = NULL) {
  check_adjacency(adjacency)
  if (!is.character(permanent) || anyNA(permanent)) {
    abort_input("`permanent` must be a character vector of node identifiers.")
  }
  check_flag(modularity)
  check_flag(robustness)
  check_number(louvain_runs, lower = 1)
  check_number(robustness_reps, lower = 1)
  check_number(seed, allow_null = TRUE)
  if (modularity) check_installed("igraph", "compute modularity")
  A <- adjacency
  storage.mode(A) <- "integer"
  A[A != 0L] <- 1L
  web_metrics_core(A, permanent, modularity, robustness, as.integer(louvain_runs),
                   as.integer(robustness_reps), seed)
}

check_adjacency <- function(a, call = sys.call(-1)) {
  if (!is.matrix(a) || !(is.numeric(a) || is.logical(a))) {
    abort_input("`adjacency` must be a numeric or logical matrix.", call = call)
  }
  if (nrow(a) != ncol(a) || is.null(rownames(a)) || !identical(rownames(a), colnames(a))) {
    abort_input("`adjacency` must be square with identical row and column names (node identifiers).",
                call = call)
  }
  if (anyNA(a)) abort_input("`adjacency` must not contain NA.", call = call)
  invisible(a)
}

# Indicators of one web (A: integer 0/1 matrix, resource x consumer, with
# dimnames). Port of network_metrics() of the original analysis.
web_metrics_core <- function(A, permanent, modularity, robustness, louvain_runs,
                             robustness_reps, seed) {
  S_taxa <- nrow(A)
  diag(A) <- 0L
  isolated <- rowSums(A) == 0 & colSums(A) == 0
  A <- A[!isolated, !isolated, drop = FALSE]
  if (nrow(A) < 2L || sum(A) == 0) {
    return(empty_web_metrics(S_taxa, sum(isolated), modularity, robustness))
  }
  ts <- trophic_species(A)
  B <- aggregate_web(A, ts)
  S <- nrow(B)
  L <- sum(B)
  n_prey <- colSums(B)
  n_pred <- rowSums(B)
  basal <- n_prey == 0
  top <- n_prey > 0 & n_pred == 0
  ls <- L / S
  tl <- trophic_level_prey_avg(B)
  tl_s <- trophic_level_short(B)
  omn <- omnivory_index(B, tl)
  out <- data.frame(
    S_taxa = S_taxa, n_isolated = sum(isolated), S = S, L = L, LS = ls, C = L / S^2,
    prop_basal = mean(basal), prop_intermediate = mean(!basal & !top), prop_top = mean(top),
    gen_mean = mean(n_prey[!basal]), gen_sd = stats::sd(n_prey / ls),
    vul_mean = if (any(n_pred > 0)) mean(n_pred[n_pred > 0]) else NA_real_,
    vul_sd = stats::sd(n_pred / ls),
    TL_mean = mean(tl), TL_max = if (all(is.na(tl))) NA_real_ else max(tl),
    TL_short_mean = if (all(is.na(tl_s))) NA_real_ else mean(tl_s, na.rm = TRUE),
    TL_status = attr(tl, "status"),
    prop_omnivores = mean_or_na(omn[!basal] > 0.01),
    omnivory_mean = mean_or_na(omn[!basal]),
    MxSim_taxa = max_trophic_similarity(A),
    MxSim_trophic_species = max_trophic_similarity(B),
    stringsAsFactors = FALSE
  )
  with_seed_maybe <- function(expr) {
    if (is.null(seed)) force(expr) else withr::with_seed(seed, expr)
  }
  if (modularity) {
    mods <- with_seed_maybe(modules_louvain(B, louvain_runs))
    out$modularity <- mods$Q
    out$n_modules <- length(unique(mods$membership))
  }
  if (robustness) {
    ts_rob <- trophic_species(A, extra = rownames(A) %in% permanent)
    B_rob <- aggregate_web(A, ts_rob)
    ts_perm <- vapply(seq_len(nrow(B_rob)), function(i) {
      all(rownames(A)[ts_rob == i] %in% permanent)
    }, logical(1))
    rob <- with_seed_maybe(robustness_r50(B_rob, removable = !ts_perm, reps = robustness_reps))
    keep <- !rownames(A) %in% permanent
    A_liv <- A[keep, keep, drop = FALSE]
    iso_liv <- rowSums(A_liv) == 0 & colSums(A_liv) == 0
    A_liv <- A_liv[!iso_liv, !iso_liv, drop = FALSE]
    B_liv <- if (nrow(A_liv)) aggregate_web(A_liv, trophic_species(A_liv)) else A_liv
    rob_liv <- with_seed_maybe(robustness_r50(B_liv, removable = rep(TRUE, nrow(B_liv)),
                                              reps = robustness_reps))
    out$R50_random <- rob[["random"]]
    out$R50_connected <- rob[["connected"]]
    out$R50_random_living <- rob_liv[["random"]]
    out$R50_connected_living <- rob_liv[["connected"]]
  }
  out$status <- "ok"
  out
}

mean_or_na <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0L) NA_real_ else mean(x)
}

empty_web_metrics <- function(S_taxa, n_isolated, modularity, robustness) {
  out <- data.frame(
    S_taxa = S_taxa, n_isolated = n_isolated, S = 0L, L = 0L, LS = NA_real_, C = NA_real_,
    prop_basal = NA_real_, prop_intermediate = NA_real_, prop_top = NA_real_,
    gen_mean = NA_real_, gen_sd = NA_real_, vul_mean = NA_real_, vul_sd = NA_real_,
    TL_mean = NA_real_, TL_max = NA_real_, TL_short_mean = NA_real_, TL_status = "empty",
    prop_omnivores = NA_real_, omnivory_mean = NA_real_, MxSim_taxa = NA_real_,
    MxSim_trophic_species = NA_real_, stringsAsFactors = FALSE
  )
  if (modularity) {
    out$modularity <- NA_real_
    out$n_modules <- NA_integer_
  }
  if (robustness) {
    out$R50_random <- NA_real_
    out$R50_connected <- NA_real_
    out$R50_random_living <- NA_real_
    out$R50_connected_living <- NA_real_
  }
  out$status <- "empty_web"
  out
}

# Trophic species: nodes with the same set of prey and of predators. `extra`
# (one value per node) is added to the key, so that only nodes that also
# share it are merged. Returns the trophic species of each node (integer, in
# order of first appearance).
trophic_species <- function(A, extra = NULL) {
  key <- vapply(seq_len(nrow(A)), function(i) {
    paste(paste(which(A[, i] > 0), collapse = ","), "|",
          paste(which(A[i, ] > 0), collapse = ","))
  }, character(1))
  if (!is.null(extra)) key <- paste(key, extra, sep = "#")
  match(key, unique(key))
}

# Web aggregated into trophic species (names join the nodes of each group).
aggregate_web <- function(A, ts) {
  k <- max(ts)
  B <- matrix(0L, k, k)
  idx <- which(A > 0, arr.ind = TRUE)
  if (nrow(idx)) B[cbind(ts[idx[, 1L]], ts[idx[, 2L]])] <- 1L
  diag(B) <- 0L
  lab <- vapply(seq_len(k), function(i) paste(rownames(A)[ts == i], collapse = "+"), character(1))
  dimnames(B) <- list(lab, lab)
  B
}

# Prey-averaged trophic level: TL_j = 1 + mean TL of its prey; species
# without prey have TL 1. Solves (I - Q) TL = 1; if singular, tries the
# generalised inverse and keeps it only if it satisfies the equations.
trophic_level_prey_avg <- function(B, tol = 1e-8) {
  n <- nrow(B)
  n_prey <- colSums(B)
  Q <- t(B) / pmax(n_prey, 1)
  M <- diag(n) - Q
  status <- "solve"
  tl <- tryCatch(solve(M, rep(1, n)), error = function(e) NULL)
  if (is.null(tl)) {
    status <- "ginv"
    tl <- as.numeric(MASS::ginv(M) %*% rep(1, n))
  }
  resid <- max(abs(M %*% tl - 1))
  if (!is.finite(resid) || resid > tol || any(!is.finite(tl))) {
    status <- "unresolved"
    tl <- rep(NA_real_, n)
  }
  structure(stats::setNames(as.numeric(tl), rownames(B)), status = status)
}

# Shortest-chain trophic level: 1 + number of links of the shortest chain
# from any species without prey (breadth-first search); NA when the species
# cannot be reached from one.
trophic_level_short <- function(B) {
  n <- nrow(B)
  dist <- rep(Inf, n)
  frontier <- which(colSums(B) == 0)
  if (!length(frontier)) {
    return(stats::setNames(rep(NA_real_, n), rownames(B)))
  }
  dist[frontier] <- 0
  d <- 0
  while (length(frontier)) {
    reach <- colSums(B[frontier, , drop = FALSE]) > 0
    nxt <- which(reach & is.infinite(dist))
    d <- d + 1
    dist[nxt] <- d
    frontier <- nxt
  }
  tl <- 1 + dist
  tl[is.infinite(tl)] <- NA_real_
  stats::setNames(as.numeric(tl), rownames(B))
}

# Omnivory of each consumer: standard deviation of the trophic levels of its
# prey (0 with one prey, NA without prey).
omnivory_index <- function(B, tl) {
  vapply(seq_len(ncol(B)), function(j) {
    p <- which(B[, j] > 0)
    if (length(p) < 2L) {
      return(if (length(p)) 0 else NA_real_)
    }
    stats::sd(tl[p])
  }, numeric(1))
}

# Mean maximum trophic similarity: shared prey and predators / prey and
# predators of either.
max_trophic_similarity <- function(B) {
  n <- nrow(B)
  if (n < 2L) {
    return(NA_real_)
  }
  X <- cbind(t(B), B) > 0
  inter <- X %*% t(X)
  tot <- rowSums(X)
  uni <- outer(tot, tot, "+") - inter
  sim <- ifelse(uni > 0, inter / uni, 0)
  diag(sim) <- NA
  mean(apply(sim, 1, max, na.rm = TRUE))
}

# Modularity of the best of `runs` Louvain partitions of the undirected web.
modules_louvain <- function(B, runs) {
  if (sum(B) == 0) {
    return(list(Q = NA_real_, membership = rep(NA_integer_, nrow(B))))
  }
  g <- igraph::graph_from_adjacency_matrix(pmax(B, t(B)), mode = "undirected")
  best <- NULL
  for (r in seq_len(runs)) {
    cl <- igraph::cluster_louvain(g)
    q <- igraph::modularity(cl)
    if (is.null(best) || q > best$Q) {
      best <- list(Q = q, membership = as.integer(igraph::membership(cl)))
    }
  }
  best
}

# R50 for one removal order: proportion of primary removals needed to lose at
# least `threshold` of the removable species (secondary extinction when a
# consumer loses all its prey).
r50 <- function(B, removable, removal_order, threshold = 0.5) {
  n <- nrow(B)
  S <- sum(removable)
  if (S == 0) {
    return(NA_real_)
  }
  had_prey <- colSums(B) > 0
  alive <- rep(TRUE, n)
  primary <- 0
  for (v in removal_order) {
    if (!alive[v]) next
    alive[v] <- FALSE
    primary <- primary + 1
    repeat {
      prey_alive <- colSums(B[alive, , drop = FALSE]) > 0
      dead <- alive & had_prey & !prey_alive & removable
      if (!any(dead)) break
      alive[dead] <- FALSE
    }
    if (sum(!alive & removable) >= threshold * S) {
      return(primary / S)
    }
  }
  1
}

robustness_r50 <- function(B, removable, reps) {
  idx <- which(removable)
  if (!length(idx)) {
    return(c(random = NA_real_, connected = NA_real_))
  }
  rnd <- mean(vapply(seq_len(reps), function(r) {
    r50(B, removable, idx[sample.int(length(idx))])
  }, numeric(1)))
  deg <- rowSums(B) + colSums(B)
  con <- r50(B, removable, idx[order(-deg[idx], rownames(B)[idx], method = "radix")])
  c(random = rnd, connected = con)
}
