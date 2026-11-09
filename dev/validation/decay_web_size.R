# Taxonomic decay with and without the size of the web as a covariate
# (version 0.3.0): cross-validation of the default rule (the 12 webs other
# than Tuesday Lake, left out by system) with the Gray registry and with
# Gray + GloBI (every registry at every level). The covariate is
# log(number of nodes of the web / 50); without it, the decay of the
# fifth phase.
#
# Run from the root of the package (about 30 minutes):
#   Rscript dev/validation/decay_web_size.R <globi_registry.rds>
pkgload::load_all(".", quiet = TRUE)
g <- readRDS(commandArgs(trailingOnly = TRUE)[1])
out_dir <- file.path("dev", "validation")
gray <- gray2015[gray2015$study_id != "12", ]
setups <- list(gray = list(registry = gray, exclusions = NULL),
               both_all = list(registry = combine_registries(gray, g$globi_test),
                               exclusions = g$exclusions))
src <- setdiff(gray2015_sources$study_id[gray2015_sources$validation], "12")
out <- list()
webs <- list()
for (nm in names(setups)) {
  st <- setups[[nm]]
  v <- validate_rule(st$registry, src, groups = list(c("13", "50")), rule = matching_rule(),
                     exclusions = st$exclusions, quiet = TRUE)
  v0 <- v
  v0$decay$web_nodes <- NULL
  s_size <- score_rule(v, monotone = FALSE, seed = 2026)
  s_none <- score_rule(v0, monotone = FALSE, seed = 2026)
  d <- s_size$webs$brier - s_none$webs$brier
  out[[nm]] <- data.frame(source = nm, brier_decay = s_none$summary$brier,
                          brier_decay_web_size = s_size$summary$brier, diff = mean(d),
                          diff_se = stats::sd(d) / sqrt(length(d)), webs_better = sum(d < 0),
                          expected_decay = s_none$summary$expected_background,
                          expected_decay_web_size = s_size$summary$expected_background,
                          missed = s_size$summary$n_missed, stringsAsFactors = FALSE)
  webs[[nm]] <- data.frame(source = nm, web = s_size$webs$study_id, nodes = v$webs$n_nodes,
                           missed = s_size$webs$n_missed,
                           expected_decay = s_none$webs$expected_background,
                           expected_decay_web_size = s_size$webs$expected_background,
                           brier_decay = s_none$webs$brier, brier_decay_web_size = s_size$webs$brier)
  message(nm, " done")
}
out <- do.call(rbind, out)
webs <- do.call(rbind, webs)
rownames(out) <- rownames(webs) <- NULL
print(out, digits = 4)
print(webs, digits = 4)
saveRDS(list(summary = out, webs = webs), file.path(out_dir, "decay_web_size.rds"))
