# Validation of the predictive model of the calibration: for each validation
# web, the model is fitted with the other 12 and the observed precision of that
# web is checked against the 95% predictive interval (unknown system: web
# effect drawn and parameter uncertainty).
pkgload::load_all(".", quiet = TRUE)
cal <- gray2015_calibration
webs <- unique(cal$counts$study_id)
fit_calibration <- utils::getFromNamespace("fit_calibration", "bwebs")
n_sim <- 4000
res <- do.call(rbind, lapply(webs, function(w) {
  train <- cal$counts[cal$counts$study_id != w, ]
  test <- cal$counts[cal$counts$study_id == w, ]
  fit <- fit_calibration(train, quadrature = 20)
  cls <- intersect(test$evidence_class, fit$coefficients$evidence_class)
  test <- test[test$evidence_class %in% cls, ]
  R <- chol(fit$vcov)
  withr::with_seed(1, {
    th <- fit$theta + t(R) %*% matrix(stats::rnorm(length(fit$theta) * n_sim), length(fit$theta))
    u <- stats::rnorm(n_sim, 0, exp(th[nrow(th), ]))
    prec <- vapply(seq_len(n_sim), function(s) {
      p <- stats::plogis(th[match(test$evidence_class, rownames(th)), s] + u[s])
      sum(stats::rbinom(length(p), test$predicted, p)) / sum(test$predicted)
    }, numeric(1))
  })
  obs <- sum(test$observed) / sum(test$predicted)
  q <- stats::quantile(prec, c(0.025, 0.5, 0.975), names = FALSE)
  data.frame(study_id = w, observed = obs, lower = q[1], median = q[2], upper = q[3],
             covered = obs >= q[1] & obs <= q[3])
}))
print(res, digits = 3)
p <- mean(res$covered)
cat(sprintf("Coverage of the 95%% predictive interval: %d/%d = %.2f (Monte Carlo error %.2f)\n",
            sum(res$covered), nrow(res), p, sqrt(p * (1 - p) / nrow(res))))
saveRDS(res, file.path("dev", "validation", "web_coverage.rds"))
