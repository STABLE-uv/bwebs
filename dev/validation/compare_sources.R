# Comparison of sources of interactions: Gray, GloBI and both together, with
# the main rules of the comparison with Gray (dev/validation/compare_rules.R),
# and final test on Tuesday Lake.
#
# Cross-validation by system: the 12 validation webs of Gray without Tuesday
# Lake (study 12), 13 and 50 together; when each web is rebuilt, the GloBI
# studies holding at least 10% of its links are left out too (copies of the
# web in other compilations; dev/globi/process_globi.R). Neither study 12 nor
# the GloBI studies with Tuesday Lake are in the registry.
#
# Sources:
#   gray         Gray registry (result of compare_rules.R);
#   globi        GloBI alone (GloBI is the primary registry);
#   both_family  Gray and GloBI; GloBI up to family (Gray primary,
#                supplementary_max_level = "family");
#   both_class   Gray and GloBI; GloBI at every level (Gray primary);
#   both_all     Gray and GloBI at every level and for downward inheritance
#                (default of matching_rule(), primary_registry = NULL).
# With a single registry, primary_registry = NULL is the same as naming it.
#
# Run from the root of the package, with the registry of process_globi.R:
#   Rscript dev/validation/compare_sources.R <globi_registry.rds> source globi
#   Rscript dev/validation/compare_sources.R <globi_registry.rds> source both_family
#   Rscript dev/validation/compare_sources.R <globi_registry.rds> source both_class
#   Rscript dev/validation/compare_sources.R <globi_registry.rds> source both_all 1 2
#   Rscript dev/validation/compare_sources.R <globi_registry.rds> source both_all 2 2
#   Rscript dev/validation/compare_sources.R <globi_registry.rds> join both_all 2
#   Rscript dev/validation/compare_sources.R <globi_registry.rds> merge
#   Rscript dev/validation/compare_sources.R <globi_registry.rds> tuesday
pkgload::load_all(".", quiet = TRUE)
args <- commandArgs(trailingOnly = TRUE)
g <- readRDS(args[1])
mode <- args[2]
out_dir <- file.path("dev", "validation")
gray <- gray2015[gray2015$study_id != "12", ]
both <- combine_registries(gray, g$globi_test)
src <- setdiff(gray2015_sources$study_id[gray2015_sources$validation], "12")
groups <- list(c("13", "50"))
selected <- c("climb/family/class/5", "climb/family/class/1", "climb/genus/class/1",
              "climb/family/order/1", "fixed/family", "fixed/genus", "fixed/exact",
              "ladder/class/1")
# rules that use levels above family: the only ones that change when GloBI is
# used at every level
above_family <- c("climb/family/class/5", "climb/family/class/1", "climb/genus/class/1",
                  "climb/family/order/1", "ladder/class/1")
sources <- list(
  gray = list(rules = candidate_rules(primary_registry = "gray2015")[selected], train = "gray2015"),
  globi = list(rules = candidate_rules(primary_registry = "globi")[selected], train = "globi"),
  both_family = list(rules = candidate_rules(primary_registry = "gray2015")[selected], train = NULL),
  both_class = list(rules = candidate_rules(primary_registry = "gray2015",
                                            supplementary_max_level = "class")[above_family],
                    train = NULL),
  # both at every level and also for downward inheritance (the default of
  # matching_rule(), primary_registry = NULL)
  both_all = list(rules = candidate_rules(primary_registry = NULL)[selected], train = NULL)
)
file_of <- function(s) file.path(out_dir, sprintf("source_comparison_%s.rds", s))

if (mode == "source") {
  # optionally, a part of the rules: source <s> <part> <n_parts>
  s <- args[3]
  cfg <- sources[[s]]
  rules <- cfg$rules
  out_file <- file_of(s)
  if (length(args) >= 5L) {
    part <- as.integer(args[4])
    n_parts <- as.integer(args[5])
    rules <- rules[seq_along(rules) %% n_parts == part %% n_parts]
    out_file <- sub("\\.rds$", sprintf("_part%d.rds", part), out_file)
  }
  t0 <- Sys.time()
  cmp <- compare_rules(both, sources = src, groups = groups, rules = rules,
                       train_registries = cfg$train, exclusions = g$exclusions, seed = 2026)
  cat("time:", format(Sys.time() - t0), "\n")
  print(cmp, n = 30)
  saveRDS(cmp, out_file)
  quit(save = "no")
}

if (mode == "join") {
  # joins the parts of a source into a single result, as compare_rules()
  s <- args[3]
  parts <- lapply(seq_len(as.integer(args[4])), function(k) {
    readRDS(sub("\\.rds$", sprintf("_part%d.rds", k), file_of(s)))
  })
  tab <- do.call(rbind, lapply(parts, `[[`, "table"))
  scores <- do.call(c, lapply(parts, `[[`, "scores"))
  validations <- do.call(c, lapply(parts, `[[`, "validations"))
  best <- tab$rule[which.min(tab$brier)]
  tab$brier_diff <- tab$brier - min(tab$brier)
  tab$brier_diff_se <- vapply(tab$rule, function(nm) {
    d <- scores[[nm]]$webs$brier - scores[[best]]$webs$brier
    stats::sd(d) / sqrt(length(d))
  }, numeric(1))
  tab <- tab[order(tab$brier, method = "radix"), ]
  rownames(tab) <- NULL
  cmp <- structure(list(table = tab, scores = scores[tab$rule], validations = validations[tab$rule],
                        settings = parts[[1]]$settings), class = "bwebs_comparison")
  print(cmp, n = 30)
  saveRDS(cmp, file_of(s))
  quit(save = "no")
}

if (mode == "merge") {
  # Gray: results of compare_rules.R (same webs, groups and seeds; without
  # exclusions, because the Gray registry has no copies of its webs)
  gray_cmp <- readRDS(file.path(out_dir, "rule_comparison_gray.rds"))
  cmps <- c(list(gray = list(table = gray_cmp$table[gray_cmp$table$rule %in% selected, ],
                             scores = gray_cmp$scores[selected])),
            lapply(stats::setNames(nm = setdiff(names(sources), "gray")), function(s) readRDS(file_of(s))))
  tab <- do.call(rbind, lapply(names(cmps), function(s) {
    t <- cmps[[s]]$table
    cbind(source = s, t[c("rule", "brier", "brier_se", "skill", "tss", "precision", "recall",
                          "n_predicted", "n_uncalibrated")], stringsAsFactors = FALSE)
  }))
  web_brier <- function(s, r) {
    w <- cmps[[s]]$scores[[r]]$webs
    stats::setNames(w$brier, w$study_id)[src]
  }
  best <- which.min(tab$brier)
  b_best <- web_brier(tab$source[best], tab$rule[best])
  tab$brier_diff <- tab$brier - tab$brier[best]
  tab$brier_diff_se <- vapply(seq_len(nrow(tab)), function(i) {
    d <- web_brier(tab$source[i], tab$rule[i]) - b_best
    stats::sd(d) / sqrt(length(d))
  }, numeric(1))
  # difference from the same rule with Gray (effect of the source)
  tab$diff_vs_gray <- tab$diff_vs_gray_se <- NA_real_
  for (i in seq_len(nrow(tab))) {
    d <- web_brier(tab$source[i], tab$rule[i]) - web_brier("gray", tab$rule[i])
    tab$diff_vs_gray[i] <- mean(d)
    tab$diff_vs_gray_se[i] <- stats::sd(d) / sqrt(length(d))
  }
  tab <- tab[order(tab$brier, method = "radix"), ]
  rownames(tab) <- NULL
  print(tab, digits = 3, row.names = FALSE)
  saveRDS(list(table = tab, scores = lapply(cmps, `[[`, "scores"),
               settings = list(sources = src, groups = groups,
                               exclusions = g$exclusions, download = g$download)),
          file.path(out_dir, "source_comparison.rds"))
  quit(save = "no")
}

if (mode == "tuesday") {
  # Final test: calibration of each combination with its 12 webs, metaweb of
  # the 82 taxa of Tuesday Lake and score against the web of each year.
  # With sources after "tuesday", only those (they replace their rows in the
  # previous result): tuesday both_all
  res <- readRDS(file.path(out_dir, "source_comparison.rds"))
  gray_cmp <- readRDS(file.path(out_dir, "rule_comparison_gray.rds"))
  tl_file <- file.path(out_dir, "tuesday_lake_sources.rds")
  run <- if (length(args) > 2L) args[-(1:2)] else names(sources)
  years <- c("TL84", "TL86")
  present <- split(tuesday_lake_community$node_id, tuesday_lake_community$sample_id)
  build_combo <- function(s, nm) {
    cmp <- if (s == "gray") gray_cmp else readRDS(file_of(s))
    reg_s <- if (is.null(sources[[s]]$train)) both else
      both[both$registry %in% sources[[s]]$train, , drop = FALSE]
    v <- cmp$validations[[nm]]
    list(calibration = calibrate_links(v, seed = 2026, quiet = TRUE),
         metaweb = build_metaweb(prepare_nodes(tuesday_lake_taxa, rule = v$settings$rule), reg_s))
  }
  rows <- list()
  keep <- list()
  for (s in run) {
    for (nm in names(sources[[s]]$rules)) {
      t0 <- Sys.time()
      k <- build_combo(s, nm)
      for (y in years) {
        obs <- tuesday_lake_links[tuesday_lake_links$sample_id == y, ]
        sc <- score_metaweb(k$metaweb, k$calibration, obs, nodes = present[[y]])$summary
        rows[[length(rows) + 1L]] <- cbind(source = s, rule = nm, year = y, sc,
                                           converged = k$metaweb$diagnostics$converged,
                                           stringsAsFactors = FALSE)
      }
      keep[[paste(s, nm, sep = " | ")]] <- k
      message(s, " | ", nm, ": ", format(Sys.time() - t0))
    }
  }
  tl <- do.call(rbind, rows)
  if (length(args) > 2L && file.exists(tl_file)) {
    tl <- rbind(readRDS(tl_file)$scores, tl)
    tl <- tl[!duplicated(tl[c("source", "rule", "year")], fromLast = TRUE), , drop = FALSE]
  }
  tl_mean <- stats::aggregate(cbind(brier, skill, tss, precision, recall) ~ source + rule, tl, mean,
                              na.rm = TRUE, na.action = stats::na.pass)
  out <- merge(res$table[c("source", "rule", "brier", "brier_diff", "brier_diff_se", "skill", "tss")],
               tl_mean, by = c("source", "rule"), suffixes = c("_cv", "_tl"))
  out <- out[order(out$brier_cv, method = "radix"), ]
  out$rank_cv <- seq_len(nrow(out))
  out$rank_tl <- rank(out$brier_tl)
  print(out, digits = 3, row.names = FALSE)
  cat("Spearman (Brier CV vs Tuesday Lake):",
      round(stats::cor(out$brier_cv, out$brier_tl, method = "spearman"), 2), "\n")

  # Indicators of 1000 sampled webs: the default rule with each source and the
  # best combination of the cross-validation
  observed_ind <- do.call(rbind, lapply(years, function(y) {
    ids <- present[[y]]
    A <- matrix(0, length(ids), length(ids), dimnames = list(ids, ids))
    l <- tuesday_lake_links[tuesday_lake_links$sample_id == y, ]
    A[cbind(l$resource, l$consumer)] <- 1
    cbind(sample_id = y, web_indicators(A))
  }))
  ind_cols <- c("S", "L", "LS", "C", "prop_basal", "prop_top", "TL_mean", "prop_omnivores",
                "MxSim_trophic_species")
  default <- rule_label(matching_rule())
  combos <- unique(c(paste(out$source[1], out$rule[1], sep = " | "),
                     paste(intersect(c("gray", "globi", "both_all"), out$source), default,
                           sep = " | ")))
  compare_ind <- do.call(rbind, lapply(combos, function(k) {
    if (is.null(keep[[k]])) {
      p <- strsplit(k, " | ", fixed = TRUE)[[1]]
      keep[[k]] <<- build_combo(p[1], p[2])
    }
    draws <- sample_metawebs(keep[[k]]$metaweb, keep[[k]]$calibration, n_draws = 1000, seed = 2026)
    ind <- calc_indicators(draws, tuesday_lake_community)
    sm <- summarize_indicators(ind, indicators = ind_cols, na_rm = TRUE)
    sm$observed <- vapply(seq_len(nrow(sm)), function(i) {
      observed_ind[observed_ind$sample_id == sm$sample_id[i], sm$indicator[i]]
    }, numeric(1))
    sm$covered <- sm$observed >= sm$lower & sm$observed <= sm$upper
    cbind(combination = k, sm[c("sample_id", "indicator", "estimate", "lower", "upper",
                                "observed", "covered", "n_omitted")])
  }))
  print(compare_ind, digits = 3, row.names = FALSE)
  saveRDS(list(scores = tl, summary = out, indicators = compare_ind, observed = observed_ind,
               calibrations = lapply(keep[combos], `[[`, "calibration")),
          tl_file)
}
