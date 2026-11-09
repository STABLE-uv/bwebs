# Sensitivity of the comparison of sources (dev/validation/compare_sources.R)
# to the secondary compilations of GloBI: metawebs and databases that gather
# interactions from other sources (among them Gray et al. 2015 and GloBI) or
# deduce them from traits (functional feeding groups). They may contain the
# validation webs generalised to family or order, which the exclusions by web
# (links with the same names) do not detect.
#
# The cross-validation of the two best rules is repeated with GloBI without
# those compilations and, for comparison, with an exclusion threshold by web
# of 0.05 instead of 0.1.
#
# Run from the root of the package:
#   Rscript dev/validation/source_sensitivity.R <globi_registry.rds>
# and, for the default of matching_rule() (both sources at every level,
# primary_registry = NULL), which is added to the result:
#   Rscript dev/validation/source_sensitivity.R <globi_registry.rds> default
pkgload::load_all(".", quiet = TRUE)
# Without the taxonomic decay (decay = FALSE), as when these results were
# computed (before version 0.3.0).
args <- commandArgs(trailingOnly = TRUE)
g <- readRDS(args[1])
out_dir <- file.path("dev", "validation")
gray <- gray2015[gray2015$study_id != "12", ]
src <- setdiff(gray2015_sources$study_id[gray2015_sources$validation], "12")
groups <- list(c("13", "50"))
compilations <- c(
  "Reji Chacko",                              # trophiCH (metaweb of Switzerland)
  "^Ho, H\\.-C\\. \\(2022\\)",                # metawebs of Ho et al. (2022)
  "Poelen, J\\. ?H\\.? Global Biotic Interactions",  # products derived from GloBI
  "Poelen.{0,60}\\(2014\\) Global biotic interactions",  # GloBI cited as the source
  "^Globi$", "^Gray$", "Gray et al\\.? \\(2015\\)",
  "Maiorano",                                 # Tetra-EU (metaweb of vertebrates)
  "freshwaterecology\\.info|Schmidt-Kloiber", # traits (feeding groups)
  "AQEM", "Schmedtje"                         # traits (feeding groups)
)
is_comp <- grepl(paste(compilations, collapse = "|"), g$globi_test$study_id, ignore.case = TRUE)
cat("Compilations:", length(unique(g$globi_test$study_id[is_comp])), "studies,",
    sum(is_comp), "records of", nrow(g$globi_test), "\n")
print(sort(table(substr(g$globi_test$study_id[is_comp], 1, 70)), decreasing = TRUE))
globi_primary <- g$globi_test[!is_comp, , drop = FALSE]

if (length(args) > 1L && args[2] == "default") {
  both <- combine_registries(gray, globi_primary)
  sh <- shared_studies(both, src, threshold = 0.1)
  ex <- sh[sh$study_id %in% globi_primary$study_id, c("web", "study_id")]
  t0 <- Sys.time()
  cmp <- compare_rules(both, sources = src, groups = groups,
                       rules = stats::setNames(list(matching_rule(decay = FALSE)),
                                               rule_label(matching_rule())),
                       exclusions = ex, monotone = FALSE, seed = 2026, quiet = TRUE)
  message("both_all without compilations: ", format(Sys.time() - t0))
  old <- readRDS(file.path(out_dir, "source_sensitivity.rds"))
  row <- cbind(analysis = "both_all_no_compilations",
               cmp$table[c("rule", "brier", "brier_se", "skill", "tss", "precision", "recall",
                           "n_predicted")])
  old$table <- rbind(old$table[old$table$analysis != "both_all_no_compilations", ], row)
  old$comparisons$both_all_no_compilations <- cmp
  print(old$table, digits = 3, row.names = FALSE)
  saveRDS(old, file.path(out_dir, "source_sensitivity.rds"))
  quit(save = "no")
}

rules <- c("climb/family/class/1", "climb/family/class/5")
configs <- list(
  globi = list(rules = candidate_rules(primary_registry = "globi", decay = FALSE)[rules],
               train = "globi"),
  both_family = list(rules = candidate_rules(primary_registry = "gray2015", decay = FALSE)[rules],
                     train = NULL),
  both_class = list(rules = candidate_rules(primary_registry = "gray2015",
                                            supplementary_max_level = "class", decay = FALSE)[rules],
                    train = NULL)
)
res <- list()
# 1. GloBI without compilations (exclusions by web computed again)
both <- combine_registries(gray, globi_primary)
sh <- shared_studies(both, src, threshold = 0.1)
ex <- sh[sh$study_id %in% globi_primary$study_id, c("web", "study_id")]
for (s in names(configs)) {
  t0 <- Sys.time()
  res[[paste0(s, "_no_compilations")]] <- compare_rules(
    both, sources = src, groups = groups, rules = configs[[s]]$rules,
    train_registries = configs[[s]]$train, exclusions = ex, monotone = FALSE, seed = 2026,
    quiet = TRUE)
  message(s, " without compilations: ", format(Sys.time() - t0))
}
# 2. Exclusion threshold of 0.05 (whole GloBI registry)
both <- combine_registries(gray, g$globi_test)
sh <- shared_studies(both, src, threshold = 0.05)
ex05 <- sh[sh$study_id %in% g$globi_test$study_id, c("web", "study_id")]
for (s in c("globi", "both_family")) {
  t0 <- Sys.time()
  res[[paste0(s, "_threshold_0.05")]] <- compare_rules(
    both, sources = src, groups = groups, rules = configs[[s]]$rules,
    train_registries = configs[[s]]$train, exclusions = ex05, monotone = FALSE,
    seed = 2026, quiet = TRUE)
  message(s, " threshold 0.05: ", format(Sys.time() - t0))
}
tab <- do.call(rbind, lapply(names(res), function(k) {
  cbind(analysis = k, res[[k]]$table[c("rule", "brier", "brier_se", "skill", "tss", "precision",
                                         "recall", "n_predicted")])
}))
print(tab, digits = 3, row.names = FALSE)
saveRDS(list(table = tab, comparisons = res, compilations = compilations,
             n_compilation_records = sum(is_comp)),
        file.path(out_dir, "source_sensitivity.rds"))
