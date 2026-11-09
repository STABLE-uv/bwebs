# Scores of matching rules and of probabilistic metawebs --------------------------

#' Cross-validated scores of a matching rule
#'
#' `score_rule()` measures how well a matching rule, with its calibrated link
#' probabilities, predicts the links of webs that took no part in the
#' calibration. Each validation system of [validate_rule()] is left out of
#' the calibration in turn; its inferred links take the probability of their
#' class in a new system under the calibration fitted to the other systems,
#' the pairs without an inferred link take probability 0 or, with
#' `decay = TRUE` in the rule, the probability of their evidence under the
#' decay model fitted to the other systems, and the predictions are compared
#' with the links recorded in the web.
#'
#' @param validation A [validate_rule()] result.
#' @param method,monotone,prior_sd,prior_sigma,n_chains,n_iter,n_warmup,quadrature
#'   Passed to [calibrate_links()] for each refit (by default with fewer
#'   draws than [calibrate_links()], and the chains of the decay model make
#'   2 random-walk steps per draw instead of 20: the scores only need the
#'   posterior mean of the probabilities; those of the ordered classes,
#'   `monotone = TRUE`, make 10, as in [calibrate_links()]).
#' @param seed `NULL` (default) to use the current random number stream, or a
#'   number. When given, the result is reproducible and the global random
#'   number stream is restored afterwards.
#' @param quiet If `FALSE`, report progress with messages.
#'
#' @return An object of class `bwebs_score`, a list with
#' * `webs`: one row per validation web with the pairs evaluated, recorded,
#'   predicted and predicted-and-recorded links, `n_uncalibrated` (inferred
#'   links of a class without data in the other systems, given probability
#'   0), `n_missed` (recorded links not inferred), `expected_background`
#'   (links expected among the pairs not inferred; 0 without the decay),
#'   `brier`, `p_reference`, `brier_reference`, `skill`, `tss`, `precision`
#'   and `recall`;
#' * `summary`: one row with the number of webs and the mean of each score
#'   over webs (`brier_se`: its standard error);
#' * `probabilities`: the probability of each class used for each web;
#' * `settings`: the rule and the options used.
#'
#' @details
#' For web \eqn{w} with \eqn{N_w} pairs evaluated (see [validate_rule()]),
#' the Brier score is \eqn{B_w = N_w^{-1} \sum_{ij} (p_{ij} - y_{ij})^2}{B =
#' mean((p - y)^2)}, where \eqn{y_{ij} = 1} if the link was recorded: a pair
#' of probability \eqn{p} costs \eqn{(1 - p)^2} if recorded and \eqn{p^2}
#' otherwise, so that without the decay a recorded link that was not
#' inferred costs 1. Lower is better. `brier_reference` is the score of a null
#' model that gives every pair the same probability, `p_reference`: the mean
#' connectance (recorded links over pairs evaluated) of the webs of the other
#' systems. `skill`
#' is \eqn{1 - B_w / B_{ref}}{1 - B / B_ref} (higher is better; positive when
#' the rule beats the null model). `tss` is the true skill statistic of the deterministic
#' metaweb (all inferred links present), recall minus false-positive rate, as
#' used by Gray et al. (2015, section 2.3, p. 15). The other studies of the
#' system of each web are left out of both the registry (in
#' [validate_rule()]) and the refit.
#'
#' @references
#' Gray C, Figueroa DH, Hudson LN, Ma A, Perkins D, Woodward G (2015).
#' Joining the dots: an automated method for constructing food webs from
#' compendia of published interactions. *Food Webs*, 5, 11--20.
#' \doi{10.1016/j.fooweb.2015.09.001}
#'
#' @seealso [compare_rules()] to score several rules; [score_metaweb()] for
#'   an independent web.
#' @export
#' @examplesIf interactive()
#' v <- validate_rule(gray2015, sources = c("46", "49", "58", "59", "63"))
#' score_rule(v, seed = 1)
score_rule <- function(validation, method = "bayes", monotone = TRUE, prior_sd = 1.5,
                       prior_sigma = 2.5, n_chains = 2L, n_iter = 500L, n_warmup = 500L,
                       quadrature = 20L, seed = NULL, quiet = TRUE) {
  if (!inherits(validation, "bwebs_validation")) {
    abort_input("`validation` must be created with validate_rule().")
  }
  check_choice(method, c("bayes", "ml"))
  check_flag(monotone)
  check_number(seed, allow_null = TRUE)
  check_flag(quiet)
  webs <- validation$webs
  counts <- validation$counts
  decay <- validation$decay
  run <- function() {
    rows <- lapply(seq_len(nrow(webs)), function(k) {
      w <- webs[k, ]
      out_ids <- strsplit(w$left_out, ", ", fixed = TRUE)[[1]]
      train <- counts[!counts$study_id %in% out_ids, , drop = FALSE]
      test <- counts[counts$study_id == w$study_id, , drop = FALSE]
      if (!quiet) message(sprintf("Scoring web %s (%d/%d)", w$study_id, k, nrow(webs)))
      p <- if (nrow(train) && sum(train$predicted) > 0) {
        fit <- if (method == "ml") {
          fit_calibration(train, quadrature, monotone = monotone)
        } else {
          fit_calibration_bayes(train, prior_sd, prior_sigma, as.integer(n_chains),
                                as.integer(n_iter), as.integer(n_warmup), quadrature, seed = NULL,
                                monotone = monotone, thin = 10L)
        }
        stats::setNames(fit$coefficients$p_new, fit$coefficients$evidence_class)
      } else {
        numeric(0)
      }
      # pairs not inferred: probability of the decay model fitted to the other systems
      background <- NULL
      if (!is.null(decay)) {
        background <- decay[decay$study_id == w$study_id, , drop = FALSE]
        train_d <- decay[!decay$study_id %in% out_ids, , drop = FALSE]
        background$p <- if (nrow(background) && nrow(train_d) && sum(train_d$pairs) > 0) {
          fd <- fit_decay(train_d, method = method, prior_sd = prior_sd, prior_sigma = prior_sigma,
                          n_chains = as.integer(n_chains), n_iter = as.integer(n_iter),
                          n_warmup = as.integer(n_warmup), quadrature = quadrature,
                          thin = 2L)
          decay_p_new_of(fd, background, n_nodes = w$n_nodes)
        } else {
          rep(0, nrow(background))
        }
      }
      others <- webs[!webs$study_id %in% out_ids, , drop = FALSE]
      p_ref <- if (nrow(others)) mean(others$n_recorded / others$n_pairs) else NA_real_
      list(score = web_brier(test, p, w$n_pairs, w$n_recorded, w$n_predicted,
                             w$n_predicted_observed, p_ref, background),
           p = data.frame(study_id = rep(w$study_id, length(p)), evidence_class = names(p),
                          p = unname(p), stringsAsFactors = FALSE))
    })
    list(score = do.call(rbind, lapply(rows, `[[`, "score")),
         p = do.call(rbind, lapply(rows, `[[`, "p")))
  }
  res <- if (is.null(seed)) run() else withr::with_seed(seed, run())
  out <- cbind(webs[c("study_id", "left_out", "n_pairs", "n_recorded", "n_predicted",
                      "n_predicted_observed")], res$score)
  structure(
    list(webs = out, summary = score_summary(out), probabilities = res$p,
         settings = list(rule = validation$settings$rule, decay = !is.null(decay),
                         method = method, monotone = monotone, prior_sd = prior_sd,
                         prior_sigma = prior_sigma, n_chains = n_chains, n_iter = n_iter,
                         n_warmup = n_warmup, seed = seed)),
    class = "bwebs_score"
  )
}

# Brier score and companions of one web from its counts by class (`test`),
# the probability of each class (`p`, named), the web totals and the
# reference probability given to every pair by the null model (`p_ref`).
# `background`: NULL (the pairs not inferred have probability 0) or the
# decay counts of the web with the probability `p` of each profile.
web_brier <- function(test, p, n_pairs, n_recorded, n_predicted, n_predicted_observed, p_ref,
                      background = NULL) {
  pc <- unname(p[test$evidence_class])
  uncal <- is.na(pc)
  pc[uncal] <- 0
  inferred <- sum(test$observed * (1 - pc)^2 + (test$predicted - test$observed) * pc^2)
  missed <- n_recorded - n_predicted_observed
  rest <- if (is.null(background)) {
    missed
  } else {
    sum(background$observed * (1 - background$p)^2 +
          (background$pairs - background$observed) * background$p^2)
  }
  brier <- (inferred + rest) / n_pairs
  ref <- (n_recorded * (1 - p_ref)^2 + (n_pairs - n_recorded) * p_ref^2) / n_pairs
  tpr <- if (n_recorded > 0) n_predicted_observed / n_recorded else NA_real_
  fpr <- (n_predicted - n_predicted_observed) / (n_pairs - n_recorded)
  data.frame(n_uncalibrated = as.integer(sum(test$predicted[uncal])),
             n_missed = as.integer(missed),
             expected_background = if (is.null(background)) 0 else
               sum(background$pairs * background$p),
             brier = brier, p_reference = p_ref, brier_reference = ref,
             skill = if (ref > 0) 1 - brier / ref else NA_real_, tss = tpr - fpr,
             precision = if (n_predicted > 0) n_predicted_observed / n_predicted else NA_real_,
             recall = tpr)
}

score_summary <- function(x) {
  data.frame(n_webs = nrow(x), brier = mean(x$brier),
             brier_se = stats::sd(x$brier) / sqrt(nrow(x)), skill = mean(x$skill),
             tss = mean(x$tss), precision = mean(x$precision, na.rm = TRUE),
             recall = mean(x$recall, na.rm = TRUE), n_predicted = sum(x$n_predicted),
             n_uncalibrated = sum(x$n_uncalibrated), n_missed = sum(x$n_missed),
             expected_background = sum(x$expected_background))
}

#' @export
print.bwebs_score <- function(x, ...) {
  s <- x$summary
  cat(sprintf("<bwebs_score> %s, %d webs left out by system\n",
              rule_label(x$settings$rule), s$n_webs))
  cat(sprintf("  Brier %.4f (se %.4f); skill %.3f; TSS %.3f; precision %.2f; recall %.2f\n",
              s$brier, s$brier_se, s$skill, s$tss, s$precision, s$recall))
  if (isTRUE(x$settings$decay)) {
    cat(sprintf(paste("  pairs not inferred (taxonomic decay): %.0f links expected,",
                      "%d recorded links missed by the rule\n"),
                s$expected_background, as.integer(s$n_missed)))
  }
  if (s$n_uncalibrated > 0) {
    cat(sprintf("  %d inferred link(s) of classes without data in the other systems (p = 0)\n",
                s$n_uncalibrated))
  }
  invisible(x)
}

#' Candidate matching rules
#'
#' `candidate_rules()` returns the grid of matching rules compared by
#' [compare_rules()]: fixed matching at the node's own taxon (`"exact"`),
#' genus, family or order; climbing from genus or family up to order or
#' class with thresholds of 1 or 5 partner nodes; and the WebBuilder ladder
#' up to family, order or class with the same thresholds; plus three rules
#' without downward inheritance. Every other setting is that of
#' [matching_rule()].
#'
#' @param ... Settings passed to every [matching_rule()] (for example
#'   `primary_registry` or `manual_diet`).
#'
#' @return A named list of `bwebs_rule` objects; the names describe each rule
#'   (`mode/base_level[/climb_max/min_partners]`).
#' @seealso [compare_rules()].
#' @export
#' @examples
#' names(candidate_rules())
candidate_rules <- function(...) {
  grid <- list()
  add <- function(...) {
    r <- matching_rule(...)
    grid[[rule_label(r)]] <<- r
  }
  for (b in c("exact", "genus", "family", "order")) {
    add(mode = "fixed", base_level = b, climb_max = b, ...)
  }
  for (b in c("genus", "family")) {
    for (mx in c("order", "class")) {
      for (k in c(1L, 5L)) add(mode = "climb", base_level = b, climb_max = mx, min_partners = k, ...)
    }
  }
  for (mx in c("family", "order", "class")) {
    for (k in c(1L, 5L)) add(mode = "ladder", base_level = "exact", climb_max = mx, min_partners = k, ...)
  }
  add(mode = "fixed", base_level = "family", climb_max = "family", downward_ranks = character(0), ...)
  add(mode = "climb", base_level = "family", climb_max = "class", min_partners = 1L,
      downward_ranks = character(0), ...)
  add(mode = "ladder", base_level = "exact", climb_max = "class", min_partners = 1L,
      downward_ranks = character(0), ...)
  grid
}

#' Compare matching rules by cross-validation
#'
#' `compare_rules()` validates and scores several matching rules on the same
#' validation webs ([validate_rule()] and [score_rule()]) and ranks them by
#' their mean Brier score.
#'
#' @inheritParams validate_rule
#' @param rules A named list of [matching_rule()] objects (default
#'   [candidate_rules()]).
#' @param ... Passed to [score_rule()].
#'
#' @return An object of class `bwebs_comparison`, a list with
#' * `table`: one row per rule, from the lowest mean Brier score, with the
#'   settings of the rule, the mean scores over webs, `brier_diff` (difference
#'   from the best rule) and `brier_diff_se` (its standard error, from the
#'   differences web by web);
#' * `scores`, `validations`: the [score_rule()] and [validate_rule()]
#'   results of every rule (pass a validation to [calibrate_links()] to
#'   calibrate that rule without validating it again);
#' * `settings`: sources and groups.
#'
#' @seealso [candidate_rules()], [score_rule()].
#' @export
#' @examplesIf interactive()
#' rules <- candidate_rules()[c("fixed/genus", "climb/family/class/1")]
#' cmp <- compare_rules(gray2015, sources = c("46", "49", "58", "59", "63"),
#'                      rules = rules, seed = 1)
#' cmp
compare_rules <- function(registry, sources, groups = NULL, rules = candidate_rules(),
                          train_registries = NULL, exclusions = NULL, ..., quiet = FALSE) {
  if (!is.list(rules) || !length(rules) || is.null(names(rules)) || anyNA(names(rules)) ||
      any(names(rules) == "") || anyDuplicated(names(rules)) ||
      !all(vapply(rules, inherits, logical(1), "bwebs_rule"))) {
    abort_input("`rules` must be a named list of matching_rule() objects with unique names.")
  }
  check_flag(quiet)
  validations <- scores <- list()
  for (k in seq_along(rules)) {
    nm <- names(rules)[k]
    if (!quiet) message(sprintf("Rule %s (%d/%d)", nm, k, length(rules)))
    validations[[nm]] <- validate_rule(registry, sources = sources, groups = groups,
                                       rule = rules[[k]], train_registries = train_registries,
                                       exclusions = exclusions, quiet = TRUE)
    scores[[nm]] <- score_rule(validations[[nm]], ...)
  }
  tab <- do.call(rbind, lapply(names(rules), function(nm) {
    r <- rules[[nm]]
    cbind(data.frame(rule = nm, mode = r$mode, base_level = r$base_level, climb_max = r$climb_max,
                     min_partners = r$min_partners, downward = length(r$downward_ranks) > 0L,
                     stringsAsFactors = FALSE),
          scores[[nm]]$summary)
  }))
  best <- tab$rule[which.min(tab$brier)]
  b_best <- scores[[best]]$webs$brier
  tab$brier_diff <- tab$brier - min(tab$brier)
  tab$brier_diff_se <- vapply(tab$rule, function(nm) {
    d <- scores[[nm]]$webs$brier - b_best
    stats::sd(d) / sqrt(length(d))
  }, numeric(1))
  tab <- tab[order(tab$brier, method = "radix"), ]
  rownames(tab) <- NULL
  structure(list(table = tab, scores = scores, validations = validations,
                 settings = list(sources = sources, groups = groups,
                                 train_registries = train_registries)),
            class = "bwebs_comparison")
}

#' @export
print.bwebs_comparison <- function(x, n = 10, ...) {
  t <- x$table
  cat(sprintf("<bwebs_comparison> %d rules, %d validation webs; best: %s\n",
              nrow(t), t$n_webs[1], t$rule[1]))
  out <- data.frame(rule = t$rule, brier = round(t$brier, 4),
                    diff = round(t$brier_diff, 4), diff_se = round(t$brier_diff_se, 4),
                    skill = round(t$skill, 3), tss = round(t$tss, 3),
                    precision = round(t$precision, 2), recall = round(t$recall, 2))
  print(utils::head(out, n), row.names = FALSE)
  if (nrow(out) > n) cat(sprintf("  ... and %d more rules\n", nrow(out) - n))
  invisible(x)
}

#' Score a probabilistic metaweb against an observed food web
#'
#' `score_metaweb()` compares the links of a metaweb, with the probabilities
#' of a calibration, with an independently observed food web of the same
#' taxa: the external test of a rule and its calibration. The metaweb must
#' be built from a registry without the records of the observed web.
#'
#' @param metaweb A metaweb from [build_metaweb()].
#' @param calibration A [calibrate_links()] result for the same matching
#'   rule and registries (a warning is given when the registries differ).
#' @param observed Data frame of observed links with columns `resource` and
#'   `consumer` (node identifiers of the metaweb); cannibalistic links are
#'   ignored.
#' @param nodes Node identifiers evaluated (default: all the nodes of the
#'   metaweb), for example the taxa present in one sample.
#' @param uncalibrated Probability of links without a calibrated probability
#'   (see [sample_metawebs()]): `"exclude"` (default, 0) or `"keep"` (1).
#' @param masses `NULL` (default, no size term) or the body masses of the
#'   nodes, a vector named by `node_id` (see [sample_metawebs()]).
#' @param size_model The effect of the predator/prey mass ratio used with
#'   `masses`; `NULL` (default) uses [size_model].
#'
#' @return An object of class `bwebs_web_score`, a list with `summary` (one
#'   row: pairs evaluated, observed, predicted and predicted-and-observed
#'   links, `n_missed` (observed links not inferred), `expected_background`
#'   (links expected among the pairs without a link; 0 without the decay),
#'   `brier`, `p_reference`, `brier_reference`, `skill`, `tss`, `precision`,
#'   `recall`) and `pairs` (every pair that was inferred or observed and,
#'   with the decay, every pair evaluated, with its class, `source` of its
#'   probability (`"rule"`, `"decay"` or `"none"`), probability `p` and
#'   whether it was `observed`).
#'
#' @details
#' The pairs evaluated are the resource \eqn{\to} consumer pairs among
#' `nodes` whose consumer has at least one observed prey among them, without
#' cannibalism. Each inferred link takes the probability of its class in a
#' new system (`p_new` of the calibration), manual links probability 1, and
#' pairs without an inferred link probability 0 or, with `decay = TRUE` in
#' the rule, the probability of their evidence for the taxonomic decay in a
#' new system (`p_new` of the decay model of the calibration, for a web of
#' as many nodes as are evaluated). With
#' `masses`, the logit of those probabilities is shifted by the size term of
#' each predation pair, with the weights of `size_model` (see
#' [sample_metawebs()]). The scores are those of [score_rule()]; the null
#' model gives every pair the mean connectance of the validation webs of the
#' calibration.
#'
#' @seealso [score_rule()], [tuesday_lake_taxa] for an observed web.
#' @export
#' @examples
#' nodes <- prepare_nodes(tuesday_lake_taxa)
#' # study 12 of gray2015 is Tuesday Lake (gateway has no copy of it)
#' registry <- combine_registries(gray2015[gray2015$study_id != "12", ], gateway)
#' mw <- build_metaweb(nodes, registry)
#' tl84 <- tuesday_lake_links[tuesday_lake_links$sample_id == "TL84", ]
#' present <- tuesday_lake_community$node_id[tuesday_lake_community$sample_id == "TL84"]
#' score_metaweb(mw, gray2015_gateway_calibration, tl84, nodes = present)
score_metaweb <- function(metaweb, calibration, observed, nodes = NULL, uncalibrated = "exclude",
                          masses = NULL, size_model = NULL) {
  if (!inherits(metaweb, "bwebs_metaweb")) abort_input("`metaweb` must be created with build_metaweb().")
  if (!inherits(calibration, "bwebs_calibration")) {
    abort_input("`calibration` must be created with calibrate_links().")
  }
  if (!rules_compatible(metaweb$settings$rule, calibration$settings$rule)) {
    abort_input("The calibration was fitted with a different matching rule from the one that built the metaweb.")
  }
  warn_registries(metaweb, calibration)
  check_data_frame(observed)
  check_columns(observed, c("resource", "consumer"), arg = "observed")
  check_choice(uncalibrated, c("exclude", "keep"))
  ids <- metaweb$nodes$node_id
  if (is.null(nodes)) nodes <- ids
  if (!is.character(nodes) || anyNA(nodes)) abort_input("`nodes` must be a character vector of node identifiers.")
  unknown <- setdiff(nodes, ids)
  if (length(unknown)) {
    abort_input(sprintf("`nodes` has identifiers that are not in the metaweb: %s.", format_items(unknown)))
  }
  obs <- as_plain_data_frame(observed[c("resource", "consumer")])
  bad <- setdiff(c(obs$resource, obs$consumer), ids)
  if (length(bad)) {
    abort_input(sprintf("`observed` has nodes that are not in the metaweb: %s.", format_items(bad)))
  }
  obs <- obs[obs$resource != obs$consumer & obs$resource %in% nodes & obs$consumer %in% nodes, ,
             drop = FALSE]
  obs <- obs[!duplicated(obs), , drop = FALSE]
  consumers <- unique(obs$consumer)
  n_pairs <- length(nodes) * length(consumers) - length(consumers)
  e <- metaweb$edges[!metaweb$edges$cannibal & metaweb$edges$resource %in% nodes &
                       metaweb$edges$consumer %in% consumers, , drop = FALSE]
  st <- link_status(e, calibration, uncalibrated)
  p <- link_class_probability(calibration, st$class_used, "random")
  p[st$status == "manual"] <- 1
  p[st$status == "uncalibrated_kept"] <- 1
  p[st$status == "uncalibrated_excluded"] <- 0
  use_size <- !is.null(masses)
  if (use_size) {
    size_model <- check_size_model(size_model)
    log_masses <- check_masses(masses, metaweb$nodes$node_id)
    animal <- metaweb$nodes$node_id[metaweb$nodes$role == "consumer"]
    shifted <- function(q, resource, consumer, weight) {
      size_shift_probability(q, resource, consumer, weight, size_model, log_masses, animal)
    }
    p <- shifted(p, e$resource, e$consumer, size_model$weights[["rule"]])
  }
  key <- function(r, c) paste(r, c, sep = "\u001f")
  obs_key <- key(obs$resource, obs$consumer)
  e_obs <- key(e$resource, e$consumer) %in% obs_key
  # pairs without a link: probability of their evidence for the decay, or 0
  cand <- NULL
  if (!is.null(metaweb$decay)) {
    if (is.null(calibration$decay)) {
      abort_input(paste(
        "The metaweb gives the pairs without a link a probability from the taxonomic",
        "decay (matching_rule(decay = TRUE)), but the calibration has no decay model;",
        "calibrate the links with that rule, or build the metaweb with",
        "matching_rule(decay = FALSE)."
      ))
    }
    cand <- metaweb$decay[metaweb$decay$resource %in% nodes &
                            metaweb$decay$consumer %in% consumers, , drop = FALSE]
    cand$p <- decay_p_new_of(calibration$decay, cand, n_nodes = length(nodes))
    if (use_size) cand$p <- shifted(cand$p, cand$resource, cand$consumer,
                                    size_model$weights[["decay"]])
  }
  cand_key <- if (is.null(cand)) character(0) else key(cand$resource, cand$consumer)
  missed <- obs[!obs_key %in% c(key(e$resource, e$consumer), cand_key), , drop = FALSE]
  pairs <- rbind(
    data.frame(resource = e$resource, consumer = e$consumer, evidence_class = e$evidence_class,
               source = rep("rule", nrow(e)), p = p, observed = e_obs, stringsAsFactors = FALSE),
    if (!is.null(cand)) {
      data.frame(resource = cand$resource, consumer = cand$consumer,
                 evidence_class = rep(NA_character_, nrow(cand)), source = rep("decay", nrow(cand)),
                 p = cand$p, observed = cand_key %in% obs_key, stringsAsFactors = FALSE)
    },
    data.frame(resource = missed$resource, consumer = missed$consumer,
               evidence_class = rep(NA_character_, nrow(missed)), source = rep("none", nrow(missed)),
               p = rep(0, nrow(missed)), observed = rep(TRUE, nrow(missed)),
               stringsAsFactors = FALSE)
  )
  rownames(pairs) <- NULL
  n_obs <- nrow(obs)
  brier <- sum((pairs$p - pairs$observed)^2) / n_pairs
  cw <- calibration$webs
  p_ref <- if (!is.null(cw$n_pairs) && nrow(cw)) mean(cw$n_recorded / cw$n_pairs) else NA_real_
  ref <- (n_obs * (1 - p_ref)^2 + (n_pairs - n_obs) * p_ref^2) / n_pairs
  n_pred <- nrow(e)
  n_po <- sum(e_obs)
  summary <- data.frame(
    n_nodes = length(nodes), n_consumers = length(consumers), n_pairs = n_pairs,
    n_observed = n_obs, n_predicted = n_pred, n_predicted_observed = n_po,
    n_missed = n_obs - n_po,
    expected_background = if (is.null(cand)) 0 else sum(cand$p),
    brier = brier, p_reference = p_ref, brier_reference = ref,
    skill = if (isTRUE(ref > 0)) 1 - brier / ref else NA_real_,
    tss = (if (n_obs) n_po / n_obs else NA_real_) - (n_pred - n_po) / (n_pairs - n_obs),
    precision = if (n_pred) n_po / n_pred else NA_real_,
    recall = if (n_obs) n_po / n_obs else NA_real_
  )
  structure(list(summary = summary, pairs = pairs,
                 settings = list(rule = metaweb$settings$rule, decay = !is.null(cand),
                                 uncalibrated = uncalibrated,
                                 size_model = if (use_size) size_model)),
            class = "bwebs_web_score")
}

#' @export
print.bwebs_web_score <- function(x, ...) {
  s <- x$summary
  cat(sprintf("<bwebs_web_score> %s; %d nodes, %d consumers, %d pairs, %d observed links\n",
              rule_label(x$settings$rule), s$n_nodes, s$n_consumers, s$n_pairs, s$n_observed))
  cat(sprintf("  Brier %.4f (reference %.4f, skill %.3f); TSS %.3f; precision %.2f; recall %.2f\n",
              s$brier, s$brier_reference, s$skill, s$tss, s$precision, s$recall))
  if (isTRUE(x$settings$decay)) {
    cat(sprintf(paste("  pairs without a link (taxonomic decay): %.1f links expected,",
                      "%d observed links missed by the rule\n"),
                s$expected_background, as.integer(s$n_missed)))
  }
  if (!is.null(x$settings$size_model)) {
    w <- x$settings$size_model$weights
    cat(sprintf("  size term: weight %.2g on the inferred links, %.2g on the pairs of the decay\n",
                w[["rule"]], w[["decay"]]))
  }
  invisible(x)
}
