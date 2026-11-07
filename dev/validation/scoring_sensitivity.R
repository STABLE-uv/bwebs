# Sensitivity of the comparison of rules (dev/validation/compare_rules.R):
#  - skill against the fair null model (connectance of the other systems);
#  - Brier score with a background probability p0 for the pairs without an
#    inferred link, estimated with the other systems (mean share of recorded
#    links among their pairs without an inferred link), with maximum
#    likelihood fits.
# Run from the root of the package: Rscript dev/validation/scoring_sensitivity.R
pkgload::load_all(".", quiet = TRUE)
cmp <- readRDS(file.path("dev", "validation", "rule_comparison_gray.rds"))
rows <- lapply(cmp$table$rule, function(nm) {
  v <- cmp$validations[[nm]]
  webs <- v$webs
  counts <- v$counts
  bayes <- cmp$scores[[nm]]$webs
  per_web <- lapply(seq_len(nrow(webs)), function(k) {
    w <- webs[k, ]
    out_ids <- strsplit(w$left_out, ", ", fixed = TRUE)[[1]]
    others <- webs[!webs$study_id %in% out_ids, ]
    train <- counts[!counts$study_id %in% out_ids, ]
    test <- counts[counts$study_id == w$study_id, ]
    fit <- fit_calibration(train, 20)
    p <- stats::setNames(fit$coefficients$p_new, fit$coefficients$evidence_class)
    p_ref <- mean(others$n_recorded / others$n_pairs)
    p0 <- mean((others$n_recorded - others$n_predicted_observed) / (others$n_pairs - others$n_predicted))
    b <- web_brier(test, p, w$n_pairs, w$n_recorded, w$n_predicted, w$n_predicted_observed, p_ref)
    missed <- w$n_recorded - w$n_predicted_observed
    non_inf <- w$n_pairs - w$n_predicted
    b_bg <- b$brier - missed / w$n_pairs + (missed * (1 - p0)^2 + (non_inf - missed) * p0^2) / w$n_pairs
    ref_b <- bayes[bayes$study_id == w$study_id, ]
    data.frame(brier_ml = b$brier, brier_bg = b_bg, p0 = p0,
               skill_bayes = 1 - ref_b$brier / b$brier_reference)
  })
  pw <- do.call(rbind, per_web)
  data.frame(rule = nm, brier_bayes = mean(bayes$brier), brier_ml = mean(pw$brier_ml),
             brier_background = mean(pw$brier_bg), p0_mean = mean(pw$p0),
             skill_bayes = mean(pw$skill_bayes))
})
tab <- do.call(rbind, rows)
tab$rank_bayes <- rank(tab$brier_bayes)
tab$rank_ml <- rank(tab$brier_ml)
tab$rank_background <- rank(tab$brier_background)
print(tab, digits = 4, row.names = FALSE)
saveRDS(tab, file.path("dev", "validation", "scoring_sensitivity.rds"))
