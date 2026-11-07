# Matching engine (internal) ------------------------------------------------
#
# WebBuilder logic (Gray et al. 2015, Food Webs 5:11-20) with the extensions
# of the original analysis (matching levels by node and role, three modes,
# aliases and surrogates, life-stage filter, downward inheritance from coarse
# records and generic resources). Ported from R/functions/webbuilder.R of the
# original project; dev/reference/ checks that the metaweb is identical.
#
# `nodes` is a bwebs_nodes data frame; `rx` is a registry index made by
# registry_index(). Matrices are nodes x registry rows.

role_prefix <- function(role) {
  if (role == "resource") "res_" else "con_"
}

# Registry rows that can match some node of `nodes`, at any level, by name,
# alias, surrogate, rank or generic group: in one role at least
# (`sides = "any"`) or in both (`"both"`). Rows that match no node in either
# role take no part in the matching, so dropping them gives the same metaweb
# (build_metaweb()). A link needs a node on each side, and partner counts and
# the predator rule only count rows matched on both sides; but rows matched
# on one side only still set the level of a node without partners (last
# resort) and its record counts, so `"both"` keeps the same links only when
# the matching levels converge. validate_web() uses it: only the links of
# the web are scored, and large registries become much smaller.
relevant_rows <- function(registry, nodes, rule, sides = c("any", "both")) {
  sides <- match.arg(sides)
  low <- function(x) unique(tolower(x[!is.na(x) & x != ""]))
  alias <- unlist(lapply(nodes$alias, split_alias), use.names = FALSE)
  sur <- c(if ("resource_surrogate" %in% names(nodes)) unlist(lapply(nodes$resource_surrogate, split_alias)),
           if ("consumer_surrogate" %in% names(nodes)) unlist(lapply(nodes$consumer_surrogate, split_alias)))
  extra <- low(c(alias, sur))
  names_all <- c(low(nodes$taxon), extra, low(rule$generic_groups$generic))
  rank_sets <- lapply(stats::setNames(rank_cols, rank_cols), function(r) c(low(nodes[[r]]), extra))
  side_ok <- function(pre, name_col) {
    ok <- tolower(registry[[name_col]]) %in% names_all |
      tolower(registry[[paste0(name_col, "_std")]]) %in% names_all
    for (r in rank_cols) ok <- ok | tolower(registry[[paste0(pre, r)]]) %in% rank_sets[[r]]
    ok
  }
  res <- side_ok("res_", "resource")
  con <- side_ok("con_", "consumer")
  if (sides == "both") res & con else res | con
}

# Lower-case copies of the registry columns used for matching, computed once.
registry_index <- function(reg) {
  cols <- c("resource", "consumer", paste0("res_", rank_cols), paste0("con_", rank_cols))
  rx <- lapply(reg[cols], tolower)
  rx$resource_std <- tolower(reg$resource_std)
  rx$consumer_std <- tolower(reg$consumer_std)
  rx$res_stage <- reg$res_stage
  rx$con_stage <- reg$con_stage
  rx$registry <- reg$registry
  rx$study_id <- reg$study_id
  rx$n <- nrow(reg)
  rx
}

# Records compatible with the life stage of the node (larva / adult).
stage_ok <- function(node, rx, role) {
  st <- node$life_stage
  if (is.na(st) || !st %in% c("larva", "adult")) {
    return(rep(TRUE, rx$n))
  }
  rs <- rx[[paste0(role_prefix(role), "stage")]]
  is.na(rs) | rs == st
}

# Coarse rank of each record in a role: the name IS the taxon of that rank
# and no finer rank is resolved ("Diptera" with res_order = "Diptera").
coarse_rank <- function(rx, role, ranks) {
  pre <- role_prefix(role)
  nm <- rx[[paste0(role, "_std")]]
  out <- rep(NA_character_, rx$n)
  for (r in ranks) {
    finer <- rank_cols[seq_len(match(r, rank_cols) - 1L)]
    finer_na <- rep(TRUE, rx$n)
    for (f in finer) finer_na <- finer_na & is.na(rx[[paste0(pre, f)]])
    v <- rx[[paste0(pre, r)]]
    hit <- is.na(out) & finer_na & !is.na(nm) & !is.na(v) & nm == v
    out[hit] <- r
  }
  out
}

# Records matched DIRECTLY by a node at level `lev` in a role.
match_rows <- function(node, rx, lev, role) {
  pre <- role_prefix(role)
  if (lev == "exact") {
    alias <- tolower(split_alias(node$alias))
    sur <- tolower(node_surrogates(node, role))
    names_ok <- c(tolower(node$taxon), alias, sur)
    hit <- rx[[role]] %in% names_ok
    if (length(c(alias, sur))) {
      hit <- hit | rx[[paste0(pre, "genus")]] %in% c(alias, sur)
    }
    # A surrogate may be a taxon of any rank (for example Tricladida)
    if (length(sur)) {
      for (r in setdiff(rank_cols, "genus")) hit <- hit | rx[[paste0(pre, r)]] %in% sur
    }
  } else {
    target <- node[[lev]]
    if (is.na(target)) {
      return(rep(FALSE, rx$n))
    }
    hit <- rx[[paste0(pre, lev)]] %in% tolower(target)
  }
  hit & stage_ok(node, rx, role)
}

# Cache of the records matched in compute_matching(). The records a node
# matches at a level depend only on the names, ranks and life stage involved
# (not on the pass of the matching), and nodes of the same genus, family,
# order or class share them; they are stored once as row indices.
row_cache <- function() new.env(parent = emptyenv(), hash = TRUE)

stage_key <- function(node) {
  st <- node$life_stage
  if (is.na(st) || !st %in% c("larva", "adult")) "" else st
}

key_part <- function(x) {
  if (length(x) == 0L) return("")
  x <- ifelse(is.na(x), "\u0002", tolower(x))
  paste(x, collapse = "\u0003")
}

cached_rows <- function(cache, key, n, compute) {
  if (is.null(cache)) return(compute())
  idx <- cache[[key]]
  if (is.null(idx)) {
    idx <- which(compute())
    assign(key, idx, envir = cache)
  }
  out <- logical(n)
  out[idx] <- TRUE
  out
}

# match_rows() through the cache.
match_rows_cached <- function(node, rx, lev, role, cache) {
  key <- if (lev == "exact") {
    paste("exact", role, key_part(node$taxon), key_part(split_alias(node$alias)),
          key_part(node_surrogates(node, role)), stage_key(node), sep = "\u0001")
  } else {
    paste(lev, role, key_part(node[[lev]]), stage_key(node), sep = "\u0001")
  }
  cached_rows(cache, key, rx$n, function() match_rows(node, rx, lev, role))
}

# Records of the primary registry (all records when the rule names none).
is_primary <- function(rx, rule) {
  if (is.null(rule$primary_registry)) rep(TRUE, rx$n) else rx$registry == rule$primary_registry
}

# Records of a node in a role according to its levels (see matching_rule()).
resolve_direct <- function(node, rx, role, rule, partner_mat = NULL, mode = rule$mode,
                           cache = NULL) {
  primary <- is_primary(rx, rule)
  mn <- node[[paste0(role, "_method")]]
  mx <- node[[paste0(role, "_max")]]
  bs <- node[[paste0(role, "_base")]]
  if (is.na(bs)) bs <- mn
  empty <- rep(FALSE, rx$n)
  if (mn == "none") {
    return(list(level = NA_character_, rows = empty, met = NA))
  }
  # An explicit surrogate is a decision of the user: it replaces the ladder
  if (has_surrogate(node, role)) {
    h <- match_rows_cached(node, rx, "exact", role, cache)
    return(list(level = if (any(h)) "surrogate" else "no_match", rows = h, met = NA))
  }
  n_partners <- function(h) {
    if (is.null(partner_mat)) {
      return(Inf)
    }
    others <- rownames(partner_mat) != node$node_id
    sum(rowSums(partner_mat[others, h, drop = FALSE]) > 0)
  }
  i <- match(mn, match_levels)
  b <- max(i, match(bs, match_levels))
  j <- max(b, match(mx, match_levels), na.rm = TRUE)
  supp_max <- match(rule$supplementary_max_level, match_levels)
  rows_at <- function(lev) {
    h <- match_rows_cached(node, rx, lev, role, cache)
    # above supplementary_max_level only records of the primary registry
    if (match(lev, match_levels) > supp_max) h <- h & primary
    h
  }
  per_level <- lapply(match_levels[i:j], rows_at)
  if (mode == "ladder") {
    steps <- Map(function(lev, h) list(level = lev, rows = h), match_levels[i:j], per_level)
  } else {
    pooled <- Reduce(`|`, per_level[seq_len(b - i + 1L)])
    steps <- list(list(level = match_levels[b], rows = pooled))
    if (mode == "climb" && j > b) {
      for (L in (b + 1L):j) {
        pooled <- pooled | per_level[[L - i + 1L]]
        steps <- c(steps, list(list(level = match_levels[L], rows = pooled)))
      }
    }
  }
  fallback <- NULL      # finest level with >= 1 partner node
  fallback_rec <- NULL  # finest level with >= 1 record (last resort)
  for (s in steps) {
    if (!any(s$rows)) next
    if (is.null(fallback_rec)) fallback_rec <- list(level = s$level, rows = s$rows, met = FALSE)
    n_part <- n_partners(s$rows)
    if (n_part >= rule$min_partners) {
      return(list(level = s$level, rows = s$rows, met = TRUE))
    }
    if (is.null(fallback) && n_part >= 1) {
      fallback <- list(level = s$level, rows = s$rows, met = FALSE)
    }
  }
  if (!is.null(fallback)) {
    return(fallback)
  }
  if (!is.null(fallback_rec)) {
    return(fallback_rec)
  }
  list(level = "no_match", rows = empty, met = NA)
}

# Records that reach a node by downward inheritance (excluding direct ones).
resolve_downward <- function(node, rx, role, direct_rows, coarse, rule, cache = NULL) {
  if (node[[paste0(role, "_method")]] == "none") {
    return(rep(FALSE, rx$n))
  }
  groups <- rule$generic_groups
  group_vals <- vapply(unique(groups$rank), function(col) {
    if (col %in% names(node)) key_part(node[[col]]) else key_part(NA_character_)
  }, character(1))
  key <- paste("downward", role, key_part(unlist(node[rule$downward_ranks])),
               paste(group_vals, collapse = "\u0004"), stage_key(node), sep = "\u0001")
  hit <- cached_rows(cache, key, rx$n, function() {
    pre <- role_prefix(role)
    h <- rep(FALSE, rx$n)
    coarse_ok <- role %in% rule$downward_roles & is_primary(rx, rule)
    for (r in rule$downward_ranks) {
      v <- node[[r]]
      if (is.na(v)) next
      h <- h | (coarse_ok & coarse %in% r & rx[[paste0(pre, r)]] %in% tolower(v))
    }
    for (g in seq_len(nrow(groups))) {
      col <- groups$rank[g]
      v <- if (col %in% names(node)) node[[col]] else NA_character_
      if (is.na(v) || !v %in% split_alias(groups$values[g])) next
      h <- h | rx[[role]] %in% tolower(groups$generic[g])
    }
    h & stage_ok(node, rx, role)
  })
  hit & !direct_rows
}

# Nodes x records: TRUE if the record side `role` matches the node at any
# level of its ladder (including aliases and surrogates). Used to count how
# many nodes of the network a set of records brings as partners.
relevant_matrix <- function(nodes, rx, role, cache = NULL) {
  M <- matrix(FALSE, nrow(nodes), rx$n, dimnames = list(nodes$node_id, NULL))
  for (k in seq_len(nrow(nodes))) {
    n <- nodes[k, , drop = FALSE]
    mn <- n[[paste0(role, "_method")]]
    mx <- n[[paste0(role, "_max")]]
    if (mn == "none") next
    if (has_surrogate(n, role)) {
      M[k, ] <- match_rows_cached(n, rx, "exact", role, cache)
      next
    }
    i <- match(mn, match_levels)
    j <- max(i, match(mx, match_levels), na.rm = TRUE)
    for (lev in match_levels[i:j]) M[k, ] <- M[k, ] | match_rows_cached(n, rx, lev, role, cache)
  }
  M
}

# Predators by diet: more than `majority` of their prey (network nodes that
# match as resource in their consumer records) are animals.
predator_nodes <- function(nodes, Cd, partner_mat, majority) {
  animal <- nodes$role == "consumer"
  vapply(seq_len(nrow(nodes)), function(k) {
    h <- Cd[k, ]
    if (!any(h)) {
      return(FALSE)
    }
    prey <- rowSums(partner_mat[, h, drop = FALSE]) > 0
    prey[k] <- FALSE
    if (!any(prey)) {
      return(FALSE)
    }
    sum(prey & animal) / sum(prey) > majority
  }, logical(1))
}

# Full matching of all nodes. Returns the levels used and the node x record
# matrices of direct (Rd, Cd) and downward (Rw, Cw) records.
compute_matching <- function(nodes, rx, rule) {
  n_nodes <- nrow(nodes)
  k_seq <- seq_len(n_nodes)
  # Node x record matrices are filled row by row (no list of rows kept), so
  # that large registries fit in memory.
  new_mat <- function() matrix(FALSE, n_nodes, rx$n, dimnames = list(nodes$node_id, NULL))
  coarse <- list(resource = coarse_rank(rx, "resource", rule$downward_ranks),
                 consumer = coarse_rank(rx, "consumer", rule$downward_ranks))
  cache <- row_cache()
  node_rows <- lapply(k_seq, function(k) nodes[k, , drop = FALSE])
  # Resolve the nodes in one role: the matrix of their records and, for each
  # node, its level and whether the threshold was met.
  resolve_role <- function(role, partner_mat, modes = rep(rule$mode, n_nodes)) {
    M <- new_mat()
    info <- vector("list", n_nodes)
    for (k in k_seq) {
      x <- resolve_direct(node_rows[[k]], rx, role, rule, partner_mat, mode = modes[k],
                          cache = cache)
      M[k, ] <- x$rows
      info[[k]] <- list(level = x$level, met = x$met)
    }
    list(M = M, info = info)
  }
  resolve_all <- function(rel) {
    con <- resolve_role("consumer", rel$consumer)
    # As resources, predators do not climb above their starting level
    # (modes "climb" and "ladder"), so that top predators do not inherit the
    # predators of their whole order or class.
    no_climb <- if (rule$mode != "fixed" && !rule$climb_predators) {
      predator_nodes(nodes, con$M, rel$consumer, rule$predator_majority)
    } else {
      rep(FALSE, n_nodes)
    }
    res <- resolve_role("resource", rel$resource, ifelse(no_climb, "fixed", rule$mode))
    list(res_d = res$info, con_d = con$info, predator = no_climb, Rd = res$M, Cd = con$M)
  }
  downward <- function(r) {
    Rw <- new_mat()
    Cw <- new_mat()
    for (k in k_seq) {
      Rw[k, ] <- resolve_downward(node_rows[[k]], rx, "resource", r$Rd[k, ], coarse$resource,
                                  rule, cache)
      Cw[k, ] <- resolve_downward(node_rows[[k]], rx, "consumer", r$Cd[k, ], coarse$consumer,
                                  rule, cache)
    }
    list(Rw = Rw, Cw = Cw)
  }
  levels_of <- function(r) {
    c(vapply(r$res_d, function(x) x$level, character(1)),
      vapply(r$con_d, function(x) x$level, character(1)))
  }
  # First pass: partners counted with the whole ladder of the partner
  # (potential). Next passes: with the level at which the partner was
  # resolved plus downward inheritance, so that the threshold counts the
  # partners that will really produce links. Iterate until levels are stable.
  # partners: for a consumer, the RESOURCE nodes its records bring, and the
  # other way round
  r <- resolve_all(list(resource = relevant_matrix(nodes, rx, "consumer", cache),
                        consumer = relevant_matrix(nodes, rx, "resource", cache)))
  converged <- FALSE
  for (it in seq_len(rule$max_iter)) {
    lev_old <- levels_of(r)
    w <- downward(r)
    rel <- list(resource = r$Cd | w$Cw, consumer = r$Rd | w$Rw)
    rm(w, r)
    r <- resolve_all(rel)
    rm(rel)
    if (identical(lev_old, levels_of(r))) {
      converged <- TRUE
      break
    }
  }
  w <- downward(r)
  list(
    res_level = vapply(r$res_d, function(x) x$level, character(1)),
    con_level = vapply(r$con_d, function(x) x$level, character(1)),
    res_met = vapply(r$res_d, function(x) x$met, logical(1)),
    con_met = vapply(r$con_d, function(x) x$met, logical(1)),
    predator = r$predator, converged = converged,
    Rd = r$Rd, Cd = r$Cd, Rw = w$Rw, Cw = w$Cw
  )
}
