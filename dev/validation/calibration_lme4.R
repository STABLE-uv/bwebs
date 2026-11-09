# Fit of the calibration of Gray and comparison with lme4::glmer (nAGQ = 20).
pkgload::load_all(".", quiet = TRUE)
src <- gray2015_sources$study_id[gray2015_sources$validation]
t0 <- Sys.time()
# Maximum likelihood fit (method = "ml"): the one compared with lme4.
# (the class model only: without the taxonomic decay)
cal <- calibrate_links(gray2015, sources = src, groups = gray2015_calibration$settings$groups,
                       rule = matching_rule(decay = FALSE), method = "ml", monotone = FALSE,
                       quiet = TRUE)
cat("time:", format(Sys.time() - t0), "\n")
print(cal)
print(cal$webs)
if (requireNamespace("lme4", quietly = TRUE)) {
  d <- cal$counts
  d$evidence_class <- factor(d$evidence_class, levels = cal$coefficients$evidence_class)
  f <- lme4::glmer(cbind(observed, predicted - observed) ~ 0 + evidence_class + (1 | study_id),
                   data = d, family = stats::binomial, nAGQ = 20)
  cmp <- data.frame(class = cal$coefficients$evidence_class, bwebs = cal$coefficients$estimate,
                    lme4 = unname(lme4::fixef(f)), se_bwebs = cal$coefficients$std_error,
                    se_lme4 = unname(sqrt(diag(as.matrix(stats::vcov(f))))))
  print(cmp)
  cat("sigma bwebs:", cal$sigma, " lme4:", attr(lme4::VarCorr(f)$study_id, "stddev"), "\n")
  cat("logLik lme4:", as.numeric(stats::logLik(f)), "\n")
}
# Posterior (bundled calibration) against the maximum likelihood fit
bay <- gray2015_calibration$coefficients
print(data.frame(class = bay$evidence_class, ml = cal$coefficients$estimate,
                 posterior_mean = bay$estimate, posterior_sd = bay$std_error,
                 se_ml = cal$coefficients$std_error))
cat("sigma: ml", cal$sigma, " posterior mean", gray2015_calibration$sigma, "\n")
