# Taxonomic decay in the package (version 0.3.0): matching_rule(decay = TRUE)
# with the model of R/decay.R (consumer resolved at its closest level with
# records, from the family; one effect for the resource levels within the
# family of a recorded resource; additive effect of generic records).
#
#   1. Cross-validation (the 12 webs other than Tuesday Lake, left out by
#      system) of the default rule with and without the decay, with the Gray
#      registry, with Gray + GloBI and with Gray + GATEWAy (every registry at
#      every level), as in dev/validation/compare_sources.R.
#   2. Tuesday Lake: calibration with the 12 webs, metaweb of its 82 taxa
#      with the registry without Tuesday Lake, scores of each year and
#      indicators of 1000 sampled webs (among taxa and among trophic species)
#      without the decay, with it, and with it and the body masses of each
#      year (tuesday_lake_community$mass_kg; the size term of size_model).
#
# The classes of evidence are free (monotone = FALSE), as when the decay and
# the size term were chosen; the ordered classes, the default from the
# seventh phase, are compared in dev/validation/evidence_amount.R.
#
# Run from the root of the package (about an hour):
#   Rscript dev/validation/decay_package.R <globi_registry.rds>
# or each registry in its own process, then merge:
#   Rscript dev/validation/decay_package.R <globi_registry.rds> gray
#   Rscript dev/validation/decay_package.R <globi_registry.rds> both_all
#   Rscript dev/validation/decay_package.R <globi_registry.rds> gray_gateway
#   Rscript dev/validation/decay_package.R <globi_registry.rds> merge
pkgload::load_all(".", quiet = TRUE)
args <- commandArgs(trailingOnly = TRUE)
g <- readRDS(args[1])
out_dir <- file.path("dev", "validation")
gray <- gray2015[gray2015$study_id != "12", ]
both <- combine_registries(gray, g$globi_test)
src <- setdiff(gray2015_sources$study_id[gray2015_sources$validation], "12")
groups <- list(c("13", "50"))
rule <- matching_rule()
gw_ex <- utils::read.csv(file.path("dev", "gateway", "gateway_exclusions.csv"), colClasses = "character")
setups <- list(gray = list(registry = gray, exclusions = NULL),
               both_all = list(registry = both, exclusions = g$exclusions),
               gray_gateway = list(registry = combine_registries(gray, gateway), exclusions = gw_ex))
years <- c("TL84", "TL86")
present <- split(tuesday_lake_community$node_id, tuesday_lake_community$sample_id)
masses <- lapply(split(tuesday_lake_community, tuesday_lake_community$sample_id),
                 function(d) tapply(d$mass_kg, d$node_id, mean))

# Indicators among taxa (without isolated taxa and without aggregation into
# trophic species), next to those of web_indicators() (trophic species).
ind_cols <- c("S", "L", "C", "prop_basal", "prop_top", "TL_mean", "prop_omnivores")
taxa_ind <- function(A) {
  diag(A) <- 0
  keep <- rowSums(A) > 0 | colSums(A) > 0
  if (!any(keep)) {
    return(data.frame(S_taxa = 0, C_taxa = NA_real_, basal_taxa = NA_real_, top_taxa = NA_real_,
                      TL_taxa = NA_real_))
  }
  A <- A[keep, keep, drop = FALSE]
  n_prey <- colSums(A)
  n_pred <- rowSums(A)
  data.frame(S_taxa = nrow(A), C_taxa = sum(A) / nrow(A)^2, basal_taxa = mean(n_prey == 0),
             top_taxa = mean(n_prey > 0 & n_pred == 0), TL_taxa = mean(trophic_level_prey_avg(A)))
}
all_ind <- function(A) cbind(web_indicators(A)[ind_cols], taxa_ind(A))
# Adjacency of draw d among the nodes `ids`
draw_adjacency <- function(dr, d, ids) {
  A <- matrix(0, length(ids), length(ids), dimnames = list(ids, ids))
  l <- dr$links[dr$z[, d], , drop = FALSE]
  if (!is.null(dr$background)) {
    l <- rbind(l[c("resource", "consumer")],
               dr$background$pairs[dr$background$present[[d]], c("resource", "consumer")])
  }
  l <- l[l$resource %in% ids & l$consumer %in% ids, , drop = FALSE]
  A[cbind(l$resource, l$consumer)] <- 1
  A
}
q <- function(x) sprintf("%.3f [%.3f, %.3f]", mean(x, na.rm = TRUE),
                         stats::quantile(x, 0.025, na.rm = TRUE),
                         stats::quantile(x, 0.975, na.rm = TRUE))

run_setups <- if (length(args) > 1L && args[2] != "merge") args[2] else names(setups)
results <- list()
cv_table <- tl_table <- ind_table <- list()
if (length(args) > 1L && args[2] == "merge") {
  parts <- lapply(names(setups), function(nm) {
    readRDS(file.path(out_dir, sprintf("decay_package_%s.rds", nm)))
  })
  run_setups <- character(0)
  for (pt in parts) {
    results <- c(results, pt$results)
    cv_table <- c(cv_table, list(pt$cv))
    tl_table <- c(tl_table, list(pt$tuesday))
    ind_table <- c(ind_table, list(pt$indicators))
  }
}
for (nm in run_setups) {
  st <- setups[[nm]]
  t0 <- Sys.time()
  v <- validate_rule(st$registry, src, groups = groups, rule = rule,
                     exclusions = st$exclusions, quiet = TRUE)
  v0 <- v
  v0$decay <- NULL
  message(nm, ": validation, ", format(Sys.time() - t0))
  s1 <- score_rule(v, monotone = FALSE, seed = 2026)
  s0 <- score_rule(v0, monotone = FALSE, seed = 2026)
  message(nm, ": scores, ", format(Sys.time() - t0))
  d <- s1$webs$brier - s0$webs$brier
  cv_table[[nm]] <- data.frame(
    source = nm, brier_rule = s0$summary$brier, brier_decay = s1$summary$brier,
    diff = mean(d), diff_se = stats::sd(d) / sqrt(length(d)), webs_better = sum(d < 0),
    skill_rule = s0$summary$skill, skill_decay = s1$summary$skill,
    pairs_not_inferred = sum(v$decay$pairs), missed = s1$summary$n_missed,
    expected_background = s1$summary$expected_background, stringsAsFactors = FALSE
  )
  cal <- calibrate_links(v, monotone = FALSE, seed = 2026, quiet = TRUE)
  nodes <- prepare_nodes(tuesday_lake_taxa, rule = rule)
  mw <- build_metaweb(nodes, st$registry)
  mw0 <- mw
  mw0$decay <- NULL
  dr <- sample_metawebs(mw, cal, n_draws = 1000, seed = 2026)
  dr0 <- sample_metawebs(mw0, cal, n_draws = 1000, seed = 2026)
  dr2 <- lapply(masses, function(m) sample_metawebs(mw, cal, n_draws = 1000, masses = m,
                                                     seed = 2026))
  message(nm, ": calibration and draws, ", format(Sys.time() - t0))
  tl <- ind <- list()
  for (y in years) {
    obs <- tuesday_lake_links[tuesday_lake_links$sample_id == y, ]
    a1 <- score_metaweb(mw, cal, obs, nodes = present[[y]])$summary
    a0 <- score_metaweb(mw0, cal, obs, nodes = present[[y]])$summary
    a2 <- score_metaweb(mw, cal, obs, nodes = present[[y]], masses = masses[[y]])$summary
    tl[[y]] <- data.frame(source = nm, year = y, observed = a1$n_observed,
                          missed = a1$n_missed, expected_background = a1$expected_background,
                          expected_background_size = a2$expected_background,
                          brier_rule = a0$brier, brier_decay = a1$brier, brier_size = a2$brier,
                          skill_rule = a0$skill, skill_decay = a1$skill, skill_size = a2$skill,
                          stringsAsFactors = FALSE)
    ids <- present[[y]]
    o <- obs[obs$resource %in% ids & obs$consumer %in% ids, , drop = FALSE]
    A_obs <- matrix(0, length(ids), length(ids), dimnames = list(ids, ids))
    A_obs[cbind(o$resource, o$consumer)] <- 1
    m1 <- do.call(rbind, lapply(seq_len(1000L), function(k) all_ind(draw_adjacency(dr, k, ids))))
    m0 <- do.call(rbind, lapply(seq_len(1000L), function(k) all_ind(draw_adjacency(dr0, k, ids))))
    m2 <- do.call(rbind, lapply(seq_len(1000L), function(k) {
      all_ind(draw_adjacency(dr2[[y]], k, ids))
    }))
    ind[[y]] <- list(observed = all_ind(A_obs), rule = m0, decay = m1, size = m2)
    ind_table[[paste(nm, y)]] <- data.frame(
      source = nm, year = y, indicator = names(m1),
      observed = sprintf("%.3f", unlist(all_ind(A_obs))),
      rule = vapply(m0, q, character(1)), decay = vapply(m1, q, character(1)),
      size = vapply(m2, q, character(1)), stringsAsFactors = FALSE
    )
  }
  tl_table[[nm]] <- do.call(rbind, tl)
  # light results: the scores by web, the fits and the Tuesday Lake tables
  results[[nm]] <- list(scores = s1$webs, scores_rule = s0$webs, decay_counts = v$decay,
                        calibration = list(coefficients = cal$coefficients, sigma = cal$sigma,
                                           decay = cal$decay[c("effects", "probabilities",
                                                               "sigma", "webs", "diagnostics")]),
                        tuesday = tl, indicators = ind)
  message(nm, ": done, ", format(Sys.time() - t0))
}
cv_table <- do.call(rbind, cv_table)
tl_table <- do.call(rbind, tl_table)
ind_table <- do.call(rbind, ind_table)
rownames(cv_table) <- rownames(tl_table) <- rownames(ind_table) <- NULL
print(cv_table, digits = 4)
print(tl_table, digits = 4)
print(ind_table, row.names = FALSE)
for (nm in names(results)) {
  cat("\n", nm, "\n")
  print(results[[nm]]$calibration$decay$effects, digits = 3)
  cat("sigma", results[[nm]]$calibration$decay$sigma, "\n")
}
out_file <- if (length(args) > 1L && args[2] != "merge") {
  sprintf("decay_package_%s.rds", args[2])
} else {
  "decay_package.rds"
}
saveRDS(list(cv = cv_table, tuesday = tl_table, indicators = ind_table, results = results),
        file.path(out_dir, out_file))
