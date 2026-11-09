# GATEWAy as a third registry (version 0.3.0): the default rule (with the
# taxonomic decay) with the registries combined, as in
# dev/validation/decay_package.R:
#
#   gray_gateway   Gray + GATEWAy
#   all            Gray + GloBI + GATEWAy
#   both_nowb      Gray + GloBI without the five chalk streams of Gray et al.
#                  (2015) and Thompson et al. (2017), partly extrapolated
#                  with WebBuilder
#   all_nowb       Gray + GloBI without them + GATEWAy
#
# Gray and Gray + GloBI are those of dev/validation/decay_package.rds. For
# each: cross-validation (the 12 webs other than Tuesday Lake, left out by
# system, with the studies of GloBI and GATEWAy that copy each web left out
# with it) and Tuesday Lake (calibration with the 12 webs; scores of each
# year without and with the body masses of that year).
#
# Run from the root of the package after dev/gateway/process_gateway.R
# (about 30 minutes per registry):
#   Rscript dev/validation/compare_gateway.R <globi_registry.rds> <gateway_registry.rds> [registries]
#   Rscript dev/validation/compare_gateway.R <globi_registry.rds> <gateway_registry.rds> merge
pkgload::load_all(".", quiet = TRUE)
args <- commandArgs(trailingOnly = TRUE)
g <- readRDS(args[1])
gw <- readRDS(args[2])
out_dir <- file.path("dev", "validation")
gray <- gray2015[gray2015$study_id != "12", ]
wb <- "Gray et al. (2015), Thompson et al. (2017)"
globi_nowb <- g$globi_test[!g$globi_test$study_id %in% wb, ]
src <- setdiff(gray2015_sources$study_id[gray2015_sources$validation], "12")
groups <- list(c("13", "50"))
rule <- matching_rule()
setups <- list(
  gray_gateway = list(registry = combine_registries(gray, gw$gateway_test), exclusions = gw$exclusions),
  all = list(registry = combine_registries(gray, g$globi_test, gw$gateway_test),
             exclusions = rbind(g$exclusions, gw$exclusions)),
  both_nowb = list(registry = combine_registries(gray, globi_nowb),
                   exclusions = g$exclusions[!g$exclusions$study_id %in% wb, ]),
  all_nowb = list(registry = combine_registries(gray, globi_nowb, gw$gateway_test),
                  exclusions = rbind(g$exclusions[!g$exclusions$study_id %in% wb, ], gw$exclusions))
)
years <- c("TL84", "TL86")
present <- split(tuesday_lake_community$node_id, tuesday_lake_community$sample_id)
masses <- lapply(split(tuesday_lake_community, tuesday_lake_community$sample_id),
                 function(d) tapply(d$mass_kg, d$node_id, mean))

run <- if (length(args) > 2L) strsplit(args[3], ",", fixed = TRUE)[[1]] else names(setups)
if (identical(run, "merge")) {
  parts <- lapply(names(setups), function(nm) readRDS(file.path(out_dir, sprintf("compare_gateway_%s.rds", nm))))
  dp <- readRDS(file.path(out_dir, "decay_package.rds"))
  cv <- rbind(
    data.frame(source = dp$cv$source, brier = dp$cv$brier_decay, skill = dp$cv$skill_decay,
               missed = dp$cv$missed, expected_background = dp$cv$expected_background),
    do.call(rbind, lapply(parts, function(p) p$cv[names(p$cv) %in% c("source", "brier", "skill", "missed",
                                                                      "expected_background")]))
  )
  tl <- rbind(
    data.frame(source = dp$tuesday$source, year = dp$tuesday$year, brier = dp$tuesday$brier_decay,
               brier_size = dp$tuesday$brier_size, skill = dp$tuesday$skill_decay,
               skill_size = dp$tuesday$skill_size, missed = dp$tuesday$missed,
               expected_background = dp$tuesday$expected_background),
    do.call(rbind, lapply(parts, function(p) p$tuesday))
  )
  # differences by web from Gray + GloBI
  ref <- dp$results$both_all$scores
  webs <- do.call(rbind, lapply(parts, function(p) {
    s <- p$scores
    d <- s$brier - ref$brier[match(s$study_id, ref$study_id)]
    data.frame(source = p$cv$source, diff = mean(d), diff_se = stats::sd(d) / sqrt(length(d)),
               webs_better = sum(d < 0))
  }))
  print(cv, digits = 4)
  print(webs, digits = 3)
  print(tl, digits = 4)
  parts <- lapply(parts, function(p) p[c("cv", "scores", "tuesday")])
  saveRDS(list(cv = cv, versus_both_all = webs, tuesday = tl, parts = parts),
          file.path(out_dir, "compare_gateway.rds"))
  quit(save = "no")
}
for (nm in run) {
  st <- setups[[nm]]
  t0 <- Sys.time()
  v <- validate_rule(st$registry, src, groups = groups, rule = rule, exclusions = st$exclusions,
                     quiet = TRUE)
  s <- score_rule(v, monotone = FALSE, seed = 2026)
  message(nm, ": validation and scores, ", format(Sys.time() - t0))
  cv <- data.frame(source = nm, brier = s$summary$brier, skill = s$summary$skill,
                   recall = s$summary$recall, precision = s$summary$precision,
                   missed = s$summary$n_missed, expected_background = s$summary$expected_background,
                   records = nrow(st$registry), stringsAsFactors = FALSE)
  cal <- calibrate_links(v, monotone = FALSE, seed = 2026, quiet = TRUE)
  mw <- build_metaweb(prepare_nodes(tuesday_lake_taxa, rule = rule), st$registry)
  tl <- do.call(rbind, lapply(years, function(y) {
    obs <- tuesday_lake_links[tuesday_lake_links$sample_id == y, ]
    a1 <- score_metaweb(mw, cal, obs, nodes = present[[y]])$summary
    a2 <- score_metaweb(mw, cal, obs, nodes = present[[y]], masses = masses[[y]])$summary
    data.frame(source = nm, year = y, brier = a1$brier, brier_size = a2$brier, skill = a1$skill,
               skill_size = a2$skill, missed = a1$n_missed, expected_background = a1$expected_background)
  }))
  message(nm, ": Tuesday Lake, ", format(Sys.time() - t0))
  print(cv, digits = 4)
  print(tl, digits = 4)
  saveRDS(list(cv = cv, scores = s$webs, tuesday = tl, calibration = cal),
          file.path(out_dir, sprintf("compare_gateway_%s.rds", nm)))
}
