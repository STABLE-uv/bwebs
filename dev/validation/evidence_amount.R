# Amount of evidence and order of the classes of evidence (version 0.3.0).
#
# Does the probability of a link inferred by the rule depend, beyond its
# class of evidence (the taxonomic level of its closest record), on the
# amount of evidence that supports it: how many studies and records, how
# many relatives of the consumer eat the resource or its relatives, and how
# many relatives of the resource are eaten? And does a closer class transfer
# better when it is never less probable than a more distant one?
#
#   1. validate_rule() of the default rule (the 12 webs other than Tuesday
#      Lake, 13 + 50 together) with Gray, Gray + GloBI and Gray + GATEWAy:
#      the links inferred in each web with their evidence (`links`).
#   2. Folds (each system left out) of logistic models with a normal web
#      effect fitted by maximum marginal likelihood (adaptive Gauss-Hermite,
#      as calibrate_links(method = "ml")): the class alone, or with the log
#      of the number of studies, of records, of consumer taxa and of
#      resource taxa that support the link, the proportion of the relatives
#      of the consumer with a recorded diet that eat the resource or its
#      relatives (consumer taxa / n_consumer_known) and the proportion of
#      the relatives of the resource eaten by those of the consumer
#      (resource taxa / n_resource_known). Each link of the web left out
#      takes its probability in a new system; the score is the Brier part
#      of the inferred links over the pairs of the web (the part of the
#      pairs not inferred is the same for every model).
#   3. score_rule() (Bayesian) with the classes free (monotone = FALSE, the
#      model of 0.2.0) and ordered (monotone = TRUE: exact >= genus >=
#      subfamily >= family >= order >= class on the logit scale; generic
#      free), without the taxonomic decay and with it (the same fit of the
#      decay for both class models).
#
#   4. Tuesday Lake with the ordered classes (scores and indicators of the
#      sampled webs, next to those of the free classes in decay_package.rds).
#
# Run from the root of the package, each registry in its own process (about
# 10 minutes each, and 20 for Tuesday Lake), then merge:
#   Rscript dev/validation/evidence_amount.R <globi_registry.rds> gray
#   Rscript dev/validation/evidence_amount.R <globi_registry.rds> both_all
#   Rscript dev/validation/evidence_amount.R <globi_registry.rds> gray_gateway
#   Rscript dev/validation/evidence_amount.R <globi_registry.rds> tuesday gray
#   (and both_all, gray_gateway)
#   Rscript dev/validation/evidence_amount.R <globi_registry.rds> merge
pkgload::load_all(".", quiet = TRUE)
args <- commandArgs(trailingOnly = TRUE)
out_dir <- file.path("dev", "validation")
setup_names <- c("gray", "both_all", "gray_gateway")

# Logistic model with a normal web effect, by maximum marginal likelihood
# (adaptive Gauss-Hermite quadrature at the conditional mode of each web
# effect), with the gradient of the approximation. Links with the same row
# of X in a web are pooled.
gh <- gauss_hermite(20)
fit_re <- function(X, y, web) {
  key <- paste(web, do.call(paste, c(as.data.frame(X), sep = ",")))
  first <- !duplicated(key)
  ys <- rowsum(y, key, reorder = FALSE)[, 1]
  ms <- rowsum(rep(1, length(y)), key, reorder = FALSE)[, 1]
  Xa <- X[first, , drop = FALSE]
  wa <- web[first]
  groups <- split(seq_along(ys), wa)
  K <- ncol(X)
  lw <- log(gh$w) + gh$x^2
  fg <- function(th) {
    beta <- th[seq_len(K)]
    s <- exp(th[[K + 1L]])
    tot <- 0
    grad <- numeric(K + 1L)
    for (idx in groups) {
      Xw <- Xa[idx, , drop = FALSE]
      eta <- drop(Xw %*% beta)
      md <- web_mode(eta, s, ys[idx], ms[idx])
      sd_ <- 1 / sqrt(-md$curvature)
      u <- md$u + sqrt(2) * sd_ * gh$x
      E <- outer(eta, u, "+")
      g <- colSums(ys[idx] * E - ms[idx] * softplus(E)) - u^2 / (2 * s^2)
      a <- lw + g
      amax <- max(a)
      tot <- tot + log(sqrt(2) * sd_) - log(sqrt(2 * pi) * s) + amax + log(sum(exp(a - amax)))
      om <- exp(a - amax)
      om <- om / sum(om)
      d_eta <- drop((ys[idx] - ms[idx] * stats::plogis(E)) %*% om)
      grad[seq_len(K)] <- grad[seq_len(K)] + drop(crossprod(Xw, d_eta))
      grad[K + 1L] <- grad[K + 1L] + sum(om * (u^2 / s^2 - 1))
    }
    list(value = -tot, gradient = -grad)
  }
  cache <- new.env()
  f <- function(th) {
    if (!identical(th, cache$th)) {
      cache$th <- th
      cache$v <- fg(th)
    }
    cache$v
  }
  opt <- stats::optim(rep(0, K + 1L), function(t) f(t)$value, function(t) f(t)$gradient,
                      method = "L-BFGS-B", lower = c(rep(-30, K), -5), upper = c(rep(30, K), 3),
                      control = list(maxit = 2000L, factr = 1e6))
  list(beta = stats::setNames(opt$par[seq_len(K)], colnames(X)), sigma = exp(opt$par[[K + 1L]]),
       loglik = -opt$value, converged = opt$convergence == 0L)
}
gh30 <- gauss_hermite(30)
p_new_eta <- function(eta, s) {
  P <- 0
  for (k in seq_along(gh30$x)) P <- P + gh30$w[k] * stats::plogis(eta + sqrt(2) * s * gh30$x[k])
  P / sqrt(pi)
}

classes <- evidence_classes[1:7]   # "beyond" links are treated as "class", as in sampling
features <- function(L) {
  data.frame(
    y = as.numeric(L$observed), web = L$study_id,
    class = factor(ifelse(L$evidence_class == "beyond", "class", L$evidence_class), levels = classes),
    studies = log(L$n_sources), records = log(L$n_records),
    consumer_taxa = log(L$n_consumer_taxa), resource_taxa = log(L$n_resource_taxa),
    consumer_share = L$n_consumer_taxa / L$n_consumer_known,
    resource_share = L$n_resource_taxa / L$n_resource_known,
    stringsAsFactors = FALSE
  )
}
models <- list(
  class = character(0),
  studies = "studies",
  records = "records",
  consumer_taxa = "consumer_taxa",
  resource_taxa = "resource_taxa",
  taxa = c("consumer_taxa", "resource_taxa"),
  consumer_share = "consumer_share",
  resource_share = "resource_share",
  shares = c("consumer_share", "resource_share"),
  studies_taxa = c("studies", "consumer_taxa", "resource_taxa"),
  all = c("studies", "records", "consumer_taxa", "resource_taxa", "consumer_share",
          "resource_share")
)
design <- function(d, vars) {
  X <- stats::model.matrix(~ 0 + class, d)
  colnames(X) <- sub("^class", "", colnames(X))
  if (length(vars)) X <- cbind(X, as.matrix(d[vars]))
  X
}

if (length(args) > 1L && args[2] == "tuesday") {
  # 4. Tuesday Lake with the ordered classes: calibration of the 12 webs
  # (calibrate_links(monotone = TRUE), the settings of decay_package.R),
  # metaweb of the 82 taxa with the registry without Tuesday Lake, scores of
  # each year without the decay, with it and with it and the body masses of
  # each year, and the indicators of 1000 sampled webs (among taxa and among
  # trophic species). The free classes are those of decay_package.rds (the
  # same validation, seed and decay fit).
  g <- readRDS(args[1])
  nm <- args[3]
  gray <- gray2015[gray2015$study_id != "12", ]
  src <- setdiff(gray2015_sources$study_id[gray2015_sources$validation], "12")
  gw_ex <- utils::read.csv(file.path("dev", "gateway", "gateway_exclusions.csv"),
                           colClasses = "character")
  setups <- list(gray = list(registry = gray, exclusions = NULL),
                 both_all = list(registry = combine_registries(gray, g$globi_test),
                                 exclusions = g$exclusions),
                 gray_gateway = list(registry = combine_registries(gray, gateway),
                                     exclusions = gw_ex))
  st <- setups[[nm]]
  years <- c("TL84", "TL86")
  present <- split(tuesday_lake_community$node_id, tuesday_lake_community$sample_id)
  masses <- lapply(split(tuesday_lake_community, tuesday_lake_community$sample_id),
                   function(d) tapply(d$mass_kg, d$node_id, mean))
  ind_cols <- c("S", "L", "C", "prop_basal", "prop_top", "TL_mean", "prop_omnivores")
  taxa_ind <- function(A) {
    diag(A) <- 0
    keep <- rowSums(A) > 0 | colSums(A) > 0
    if (!any(keep)) {
      return(data.frame(S_taxa = 0, C_taxa = NA_real_, basal_taxa = NA_real_,
                        top_taxa = NA_real_, TL_taxa = NA_real_))
    }
    A <- A[keep, keep, drop = FALSE]
    n_prey <- colSums(A)
    n_pred <- rowSums(A)
    data.frame(S_taxa = nrow(A), C_taxa = sum(A) / nrow(A)^2, basal_taxa = mean(n_prey == 0),
               top_taxa = mean(n_prey > 0 & n_pred == 0), TL_taxa = mean(trophic_level_prey_avg(A)))
  }
  all_ind <- function(A) cbind(web_indicators(A)[ind_cols], taxa_ind(A))
  draw_adjacency <- function(dr, d, ids) {
    A <- matrix(0, length(ids), length(ids), dimnames = list(ids, ids))
    l <- rbind(dr$links[dr$z[, d], c("resource", "consumer")],
               dr$background$pairs[dr$background$present[[d]], c("resource", "consumer")])
    l <- l[l$resource %in% ids & l$consumer %in% ids, , drop = FALSE]
    A[cbind(l$resource, l$consumer)] <- 1
    A
  }
  q <- function(x) sprintf("%.3f [%.3f, %.3f]", mean(x, na.rm = TRUE),
                           stats::quantile(x, 0.025, na.rm = TRUE),
                           stats::quantile(x, 0.975, na.rm = TRUE))
  t0 <- Sys.time()
  v <- validate_rule(st$registry, src, groups = list(c("13", "50")), rule = matching_rule(),
                     exclusions = st$exclusions, quiet = TRUE)
  cal <- calibrate_links(v, monotone = TRUE, seed = 2026, quiet = TRUE)
  message(nm, ": calibration, ", format(Sys.time() - t0))
  mw <- build_metaweb(prepare_nodes(tuesday_lake_taxa), st$registry)
  mw0 <- mw
  mw0$decay <- NULL
  dp <- readRDS(file.path(out_dir, "decay_package.rds"))
  scores <- ind <- list()
  for (y in years) {
    obs <- tuesday_lake_links[tuesday_lake_links$sample_id == y, ]
    ids <- present[[y]]
    a0 <- score_metaweb(mw0, cal, obs, nodes = ids)$summary
    a1 <- score_metaweb(mw, cal, obs, nodes = ids)$summary
    a2 <- score_metaweb(mw, cal, obs, nodes = ids, masses = masses[[y]])$summary
    ref <- dp$tuesday[dp$tuesday$source == nm & dp$tuesday$year == y, ]
    scores[[y]] <- data.frame(
      source = nm, year = y,
      brier_rule_free = ref$brier_rule, brier_rule_ordered = a0$brier,
      brier_decay_free = ref$brier_decay, brier_decay_ordered = a1$brier,
      brier_size_free = ref$brier_size, brier_size_ordered = a2$brier,
      skill_decay_free = ref$skill_decay, skill_decay_ordered = a1$skill,
      skill_size_free = ref$skill_size, skill_size_ordered = a2$skill,
      stringsAsFactors = FALSE)
    o <- obs[obs$resource %in% ids & obs$consumer %in% ids, , drop = FALSE]
    A_obs <- matrix(0, length(ids), length(ids), dimnames = list(ids, ids))
    A_obs[cbind(o$resource, o$consumer)] <- 1
    dr <- sample_metawebs(mw, cal, n_draws = 1000, seed = 2026)
    dr2 <- sample_metawebs(mw, cal, n_draws = 1000, masses = masses[[y]], seed = 2026)
    m1 <- do.call(rbind, lapply(seq_len(1000L), function(k) all_ind(draw_adjacency(dr, k, ids))))
    m2 <- do.call(rbind, lapply(seq_len(1000L), function(k) all_ind(draw_adjacency(dr2, k, ids))))
    free <- dp$indicators[dp$indicators$source == nm & dp$indicators$year == y, ]
    ind[[y]] <- data.frame(source = nm, year = y, indicator = names(m1),
                           observed = sprintf("%.3f", unlist(all_ind(A_obs))),
                           decay_free = free$decay[match(names(m1), free$indicator)],
                           decay_ordered = vapply(m1, q, character(1)),
                           size_free = free$size[match(names(m1), free$indicator)],
                           size_ordered = vapply(m2, q, character(1)), stringsAsFactors = FALSE)
  }
  out <- list(scores = do.call(rbind, scores), indicators = do.call(rbind, ind),
              calibration = list(coefficients = cal$coefficients, sigma = cal$sigma,
                                 diagnostics = cal$diagnostics[c("rhat", "ess")]))
  rownames(out$scores) <- rownames(out$indicators) <- NULL
  print(out$scores, digits = 4, row.names = FALSE)
  print(out$indicators, row.names = FALSE)
  saveRDS(out, file.path(out_dir, sprintf("evidence_amount_tuesday_%s.rds", nm)))
  message(nm, ": done, ", format(Sys.time() - t0))
  quit(save = "no")
}

if (length(args) > 1L && args[2] == "merge") {
  parts <- lapply(setup_names, function(nm) {
    readRDS(file.path(out_dir, sprintf("evidence_amount_%s.rds", nm)))
  })
  res <- list(amount = do.call(rbind, lapply(parts, `[[`, "amount")),
              order = do.call(rbind, lapply(parts, `[[`, "order")),
              effects = do.call(rbind, lapply(parts, `[[`, "effects")),
              classes = do.call(rbind, lapply(parts, `[[`, "classes")),
              results = do.call(c, lapply(parts, `[[`, "results")))
  tl_files <- file.path(out_dir, sprintf("evidence_amount_tuesday_%s.rds", setup_names))
  if (all(file.exists(tl_files))) {
    tl <- lapply(tl_files, readRDS)
    res$tuesday <- do.call(rbind, lapply(tl, `[[`, "scores"))
    res$tuesday_indicators <- do.call(rbind, lapply(tl, `[[`, "indicators"))
    res$tuesday_calibrations <- stats::setNames(lapply(tl, `[[`, "calibration"), setup_names)
    print(res$tuesday, digits = 4, row.names = FALSE)
  }
  print(res$amount, digits = 3, row.names = FALSE)
  print(res$order, digits = 4, row.names = FALSE)
  print(res$effects, digits = 3, row.names = FALSE)
  print(res$classes, digits = 3, row.names = FALSE)
  saveRDS(res, file.path(out_dir, "evidence_amount.rds"))
  quit(save = "no")
}

g <- readRDS(args[1])
gray <- gray2015[gray2015$study_id != "12", ]
src <- setdiff(gray2015_sources$study_id[gray2015_sources$validation], "12")
groups <- list(c("13", "50"))
gw_ex <- utils::read.csv(file.path("dev", "gateway", "gateway_exclusions.csv"),
                         colClasses = "character")
setups <- list(gray = list(registry = gray, exclusions = NULL),
               both_all = list(registry = combine_registries(gray, g$globi_test),
                               exclusions = g$exclusions),
               gray_gateway = list(registry = combine_registries(gray, gateway),
                                   exclusions = gw_ex))
run_setups <- if (length(args) > 1L) args[2] else setup_names
out <- list()
for (nm in run_setups) {
  st <- setups[[nm]]
  t0 <- Sys.time()
  v <- validate_rule(st$registry, src, groups = groups, rule = matching_rule(),
                     exclusions = st$exclusions, quiet = TRUE)
  message(nm, ": validation, ", format(Sys.time() - t0))
  d <- features(v$links)
  webs <- v$webs
  # 2. folds of the models of the amount of evidence
  B <- do.call(rbind, lapply(seq_len(nrow(webs)), function(k) {
    ids <- strsplit(webs$left_out[k], ", ", fixed = TRUE)[[1]]
    tr <- d[!d$web %in% ids, ]
    te <- d[d$web == webs$study_id[k], ]
    vapply(models, function(vars) {
      Xtr <- design(tr, vars)
      seen <- colSums(Xtr != 0) > 0
      Xtr <- Xtr[, seen, drop = FALSE]
      f <- fit_re(Xtr, tr$y, tr$web)
      Xte <- design(te, vars)
      p <- p_new_eta(drop(Xte[, seen, drop = FALSE] %*% f$beta), f$sigma)
      # links of a class without data in the other systems: probability 0
      p[rowSums(Xte[, !seen, drop = FALSE] != 0) > 0] <- 0
      sum((p - te$y)^2) / webs$n_pairs[k]
    }, numeric(1))
  }))
  rownames(B) <- webs$study_id
  message(nm, ": folds of the amount of evidence, ", format(Sys.time() - t0))
  dB <- B - B[, "class"]
  amount <- data.frame(source = nm, model = colnames(B),
                       covariates = vapply(models, paste, character(1), collapse = " + "),
                       brier_inferred = colMeans(B), diff = colMeans(dB),
                       diff_se = apply(dB, 2, stats::sd) / sqrt(nrow(B)),
                       webs_better = colSums(dB < 0), stringsAsFactors = FALSE)
  # effects with every web (the model with all the covariates and each alone)
  effects <- do.call(rbind, lapply(names(models)[-1], function(m) {
    f <- fit_re(design(d, models[[m]]), d$y, d$web)
    data.frame(source = nm, model = m, covariate = models[[m]],
               estimate = unname(f$beta[models[[m]]]), stringsAsFactors = FALSE)
  }))
  # 3. Bayesian cross-validation with the classes free and ordered. The
  # decay model does not depend on the classes: it is fitted once (with the
  # free classes) and its part of the score is added to the part of the
  # inferred links of each class model, scored without the decay, so that
  # the Monte Carlo error of the decay does not enter the comparison.
  s_free <- score_rule(v, monotone = FALSE, seed = 2026)
  v0 <- v
  v0$decay <- NULL
  s_free0 <- score_rule(v0, monotone = FALSE, seed = 2026)
  s_ord0 <- score_rule(v0, monotone = TRUE, seed = 2026)
  message(nm, ": scores, ", format(Sys.time() - t0))
  inferred_part <- function(s) {
    vapply(webs$study_id, function(w) {
      ct <- v$counts[v$counts$study_id == w, , drop = FALSE]
      pr <- s$probabilities[s$probabilities$study_id == w, , drop = FALSE]
      p <- pr$p[match(ct$evidence_class, pr$evidence_class)]
      p[is.na(p)] <- 0
      sum(ct$observed * (1 - p)^2 + (ct$predicted - ct$observed) * p^2)
    }, numeric(1))
  }
  background <- s_free$webs$brier * webs$n_pairs - inferred_part(s_free)
  brier_ord <- (inferred_part(s_ord0) + background) / webs$n_pairs
  dd <- brier_ord - s_free$webs$brier
  d0 <- s_ord0$webs$brier - s_free0$webs$brier
  order_tab <- data.frame(source = nm,
                          brier_free_rule = s_free0$summary$brier,
                          brier_ordered_rule = s_ord0$summary$brier, diff_rule = mean(d0),
                          diff_rule_se = stats::sd(d0) / sqrt(length(d0)),
                          webs_better_rule = sum(d0 < 0),
                          brier_free = s_free$summary$brier, brier_ordered = mean(brier_ord),
                          diff = mean(dd), diff_se = stats::sd(dd) / sqrt(length(dd)),
                          webs_better = sum(dd < 0), stringsAsFactors = FALSE)
  # class probabilities in a new system with every web (maximum likelihood)
  f_free <- fit_calibration(v$counts, 20)
  f_ord <- fit_calibration(v$counts, 20, monotone = TRUE)
  cls <- data.frame(source = nm, evidence_class = f_free$coefficients$evidence_class,
                    links = f_free$coefficients$predicted, observed = f_free$coefficients$observed,
                    p_new_free = f_free$coefficients$p_new, p_new_ordered = f_ord$coefficients$p_new,
                    stringsAsFactors = FALSE)
  out[[nm]] <- list(amount = amount, order = order_tab, effects = effects, classes = cls,
                    results = stats::setNames(list(list(
                      folds = B, scores_free = s_free$webs, brier_ordered = brier_ord,
                      scores_free_rule = s_free0$webs, scores_ordered_rule = s_ord0$webs,
                      probabilities_free = s_free0$probabilities,
                      probabilities_ordered = s_ord0$probabilities)), nm))
  print(amount, digits = 3, row.names = FALSE)
  print(order_tab, digits = 4, row.names = FALSE)
  print(cls, digits = 3, row.names = FALSE)
  message(nm, ": done, ", format(Sys.time() - t0))
}
if (length(args) > 1L) {
  saveRDS(out[[args[2]]], file.path(out_dir, sprintf("evidence_amount_%s.rds", args[2])))
} else {
  res <- list(amount = do.call(rbind, lapply(out, `[[`, "amount")),
              order = do.call(rbind, lapply(out, `[[`, "order")),
              effects = do.call(rbind, lapply(out, `[[`, "effects")),
              classes = do.call(rbind, lapply(out, `[[`, "classes")),
              results = do.call(c, lapply(out, `[[`, "results")))
  saveRDS(res, file.path(out_dir, "evidence_amount.rds"))
}
