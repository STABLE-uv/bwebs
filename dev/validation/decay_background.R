# Taxonomic distance decay (experimental, R/decay.R): cross-validation and
# Tuesday Lake.
#
# Every pair gets its evidence profile: for each level of relatedness on the
# consumer side (the consumer itself, genus, subfamily, family, order, class,
# phylum or any consumer), the closest relation on the resource side among the
# records of those consumers (decay_profiles()). Two uses are compared with
# the default rule (climb/family/class/5), with the Gray registry alone and
# with Gray + GloBI (every registry at every level):
#   1. decay instead of the rule: every pair takes the probability of its
#      profile (variants: the best evidence of any consumer level, with or
#      without a discount for how well the diet of the consumer is known, and
#      the consumer resolved at its closest level with records);
#   2. decay as the background of the rule: the links inferred by the rule
#      keep the probability of their class, and the pairs not inferred take
#      the probability of their profile instead of 0 (consumer resolved at its
#      closest level with records, from the family).
# On Tuesday Lake (Gray + GloBI), the scores of each year and the indicators
# of 1000 webs sampled with the rule alone and with the background.
#
# This was the first exploration of the decay, before it became part of the
# package (matching_rule(decay = TRUE), with the model of R/decay.R); its
# functions are in dev/validation/decay_experimental.R.
#
# Run from the root of the package (about 25 minutes):
#   Rscript dev/validation/decay_background.R <globi_registry.rds>
pkgload::load_all(".", quiet = TRUE)
source(file.path("dev", "validation", "decay_experimental.R"))
args <- commandArgs(trailingOnly = TRUE)
g <- readRDS(args[1])
out_dir <- file.path("dev", "validation")
gray <- gray2015[gray2015$study_id != "12", ]
both <- combine_registries(gray, g$globi_test)
src <- setdiff(gray2015_sources$study_id[gray2015_sources$validation], "12")
groups <- list(c("13", "50"))
# the rule as it was then, without the decay of the package
rule <- matching_rule(decay = FALSE)
setups <- list(
  gray = list(registry = gray, exclusions = NULL,
              reference = readRDS(file.path(out_dir, "rule_comparison_gray.rds"))$scores[[rule_label(rule)]]),
  both_all = list(registry = both, exclusions = g$exclusions,
                  reference = readRDS(file.path(out_dir, "source_comparison_both_all.rds"))$scores[[rule_label(rule)]])
)

# Pairs inferred by the default rule in each validation web (as validate_web())
inferred_pairs <- function(registry, exclusions) {
  r <- rule
  r$manual_diet <- character(0)
  left_out <- left_out_sets(groups, src, registry)
  out <- list()
  for (k in seq_along(src)) {
    w <- src[k]
    extra <- if (is.null(exclusions)) character(0) else exclusions$study_id[exclusions$web == w]
    vn <- registry_nodes(registry, w, r)
    keep <- !registry$study_id %in% c(left_out[[k]], extra)
    mw <- metaweb_core(vn$nodes, registry[keep, , drop = FALSE],
                       check_manual_links(NULL, vn$nodes$node_id), prefilter = "both")
    e <- mw$edges[!mw$edges$cannibal & mw$edges$consumer %in% vn$nodes$node_id[vn$studied], ,
                  drop = FALSE]
    out[[w]] <- e[c("resource", "consumer", "evidence_class")]
  }
  out
}

resolved <- function(w) {
  x <- decay_resolve(w$P, w$G, w$Q, min_known = 1, from = "family")
  w$P <- x$P
  w$G <- x$G
  w
}
summary_row <- function(setup, use, variant, method, s, ref) {
  d <- s$webs$brier - ref$webs$brier[match(s$webs$study_id, ref$webs$study_id)]
  data.frame(source = setup, use = use, variant = variant, method = method,
             brier = mean(s$webs$brier), brier_rule = mean(ref$webs$brier), diff = mean(d),
             diff_se = stats::sd(d) / sqrt(length(d)), skill = mean(s$webs$skill),
             stringsAsFactors = FALSE)
}

results <- list()
table <- list()
for (nm in names(setups)) {
  st <- setups[[nm]]
  t0 <- Sys.time()
  inf <- inferred_pairs(st$registry, st$exclusions)
  v <- decay_validate(st$registry, src, groups = groups, exclusions = st$exclusions, inferred = inf)
  message(nm, ": inferred pairs and profiles, ", format(Sys.time() - t0))
  # 1. decay instead of the rule (maximum likelihood)
  vr <- v
  vr$webs <- lapply(v$webs, resolved)
  instead <- list(
    best_evidence = decay_score(v, model = "steps", method = "ml"),
    best_evidence_known_diet = decay_score(v, model = "steps_known", method = "ml"),
    resolved_family = decay_score(vr, model = "steps", method = "ml")
  )
  for (k in names(instead)) {
    table[[length(table) + 1L]] <- summary_row(nm, "instead", k, "ml", instead[[k]], st$reference)
  }
  # 2. decay as the background of the rule (Bayesian)
  background <- decay_score_background(vr, st$reference, model = "steps", method = "bayes",
                                       seed = 2026)
  table[[length(table) + 1L]] <- summary_row(nm, "background", "resolved_family", "bayes",
                                             background, st$reference)
  results[[nm]] <- list(validation = v, inferred = inf, instead = instead, background = background)
  message(nm, ": scores, ", format(Sys.time() - t0))
}
table <- do.call(rbind, table)
print(table, digits = 3, row.names = FALSE)

# Tuesday Lake (Gray + GloBI): the rule with and without the background
v <- results$both_all$validation
fit <- decay_fit(lapply(lapply(v$webs, decay_not_inferred), resolved), model = "steps",
                 method = "bayes", n_chains = 4L, n_iter = 1000L, n_warmup = 1000L, seed = 2026)
fit$draws <- fit$draws[seq(1L, nrow(fit$draws), by = 10L), , drop = FALSE]
cal <- readRDS(file.path(out_dir, "tuesday_lake_sources.rds"))$calibrations[[
  paste("both_all", rule_label(rule), sep = " | ")]]
nodes <- prepare_nodes(tuesday_lake_taxa, rule = rule)
mw <- build_metaweb(nodes, both)
consumers <- nodes$node_id[nodes$role == "consumer"]
prof <- decay_profiles(nodes, both, consumers = consumers)
present <- split(tuesday_lake_community$node_id, tuesday_lake_community$sample_id)
key <- function(r, c) paste(r, c, sep = "\u001f")
co <- cal$coefficients
edges <- mw$edges[!mw$edges$cannibal, , drop = FALSE]
cls <- edges$evidence_class
cls[cls == "beyond"] <- "class"
edges$p <- co$p_new[match(cls, co$evidence_class)]
edges$p[is.na(edges$p)] <- 0
ind_cols <- c("S", "L", "C", "prop_basal", "prop_top", "TL_mean", "prop_omnivores")
# indicators of the taxa (without isolated taxa and without aggregation into
# trophic species) next to those of web_indicators() (trophic species)
taxa_ind <- function(A) {
  diag(A) <- 0
  keep <- rowSums(A) > 0 | colSums(A) > 0
  A <- A[keep, keep, drop = FALSE]
  n_prey <- colSums(A)
  n_pred <- rowSums(A)
  data.frame(S_taxa = nrow(A), C_taxa = sum(A) / nrow(A)^2, basal_taxa = mean(n_prey == 0),
             top_taxa = mean(n_prey > 0 & n_pred == 0), TL_taxa = mean(trophic_level_prey_avg(A)))
}
all_ind <- function(A) cbind(web_indicators(A)[ind_cols], taxa_ind(A))
tl_scores <- list()
tl_ind <- list()
set.seed(2026)
for (y in c("TL84", "TL86")) {
  obs <- tuesday_lake_links[tuesday_lake_links$sample_id == y, ]
  ids <- present[[y]]
  obs <- obs[obs$resource %in% ids & obs$consumer %in% ids, , drop = FALSE]
  pairs <- expand.grid(resource = ids, consumer = intersect(consumers, ids), stringsAsFactors = FALSE)
  pairs <- pairs[pairs$resource != pairs$consumer, , drop = FALSE]
  pk <- key(pairs$resource, pairs$consumer)
  ek <- key(edges$resource, edges$consumer)
  pairs$y <- pk %in% key(obs$resource, obs$consumer)
  pairs$inferred <- pk %in% ek
  pairs$p_rule <- 0
  pairs$p_rule[pairs$inferred] <- edges$p[match(pk[pairs$inferred], ek)]
  pp <- decay_pair_profiles(prof, pairs$resource, pairs$consumer)
  x <- decay_resolve(pp$P, pp$G, pp$Q, min_known = 1, from = "family")
  pairs$p_background <- decay_p_new(fit, x$P, x$G, pp$Q)
  pairs$p <- ifelse(pairs$inferred, pairs$p_rule, pairs$p_background)
  # scores on the pairs of score_metaweb(): consumers with observed prey
  sp <- pairs[pairs$consumer %in% unique(obs$consumer[obs$resource != obs$consumer]), , drop = FALSE]
  sc <- score_metaweb(mw, cal, tuesday_lake_links[tuesday_lake_links$sample_id == y, ], nodes = ids)$summary
  b_rule <- mean((sp$p_rule - sp$y)^2)
  b_bg <- mean((sp$p - sp$y)^2)
  tl_scores[[y]] <- data.frame(
    year = y, n_pairs = nrow(sp), observed = sum(sp$y), brier_score_metaweb = sc$brier,
    brier_rule = b_rule, brier_background = b_bg, skill_rule = 1 - b_rule / sc$brier_reference,
    skill_background = 1 - b_bg / sc$brier_reference, missed = sum(sp$y & !sp$inferred),
    expected_missed = sum(sp$p_background[!sp$inferred]), links_rule = sum(sp$p_rule),
    links_background = sum(sp$p), stringsAsFactors = FALSE
  )
  sample_ind <- function(p) {
    do.call(rbind, lapply(seq_len(1000L), function(d) {
      z <- stats::runif(length(p)) < p
      A <- matrix(0, length(ids), length(ids), dimnames = list(ids, ids))
      A[cbind(pairs$resource[z], pairs$consumer[z])] <- 1
      all_ind(A)
    }))
  }
  A_obs <- matrix(0, length(ids), length(ids), dimnames = list(ids, ids))
  A_obs[cbind(obs$resource, obs$consumer)] <- 1
  tl_ind[[y]] <- list(observed = all_ind(A_obs), rule = sample_ind(pairs$p_rule),
                      background = sample_ind(pairs$p))
}
tl_scores <- do.call(rbind, tl_scores)
print(tl_scores, digits = 4, row.names = FALSE)
for (y in names(tl_ind)) {
  q <- function(m) apply(m, 2, function(x) sprintf("%.3f [%.3f, %.3f]", mean(x),
                                                   stats::quantile(x, 0.025), stats::quantile(x, 0.975)))
  cat("\n", y, "\n")
  print(cbind(observed = sprintf("%.3f", unlist(tl_ind[[y]]$observed)), rule = q(tl_ind[[y]]$rule),
              background = q(tl_ind[[y]]$background)), quote = FALSE)
}
saveRDS(list(table = table, results = results, tuesday = list(fit = fit, scores = tl_scores,
                                                              indicators = tl_ind)),
        file.path(out_dir, "decay_background.rds"))
