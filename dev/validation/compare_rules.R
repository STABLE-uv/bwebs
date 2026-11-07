# Comparison of matching rules with the Gray registry (cross-validation by
# system, Brier score). Tuesday Lake (study 12) is left out of the registry
# and of the validation: it is kept for the final test
# (dev/validation/tuesday_lake.R). Groups: 13 and 50 are the same system.
#
# Run from the root of the package. A part of the grid can be run, to share it
# among processes:
#   Rscript dev/validation/compare_rules.R 1 2   # odd rules
#   Rscript dev/validation/compare_rules.R 2 2   # even rules
#   Rscript dev/validation/compare_rules.R merge 2
pkgload::load_all(".", quiet = TRUE)
args <- commandArgs(trailingOnly = TRUE)
out_dir <- file.path("dev", "validation")
rules <- candidate_rules()
if (length(args) == 2L && args[1] == "merge") {
  parts <- lapply(seq_len(as.integer(args[2])), function(k) {
    readRDS(file.path(out_dir, sprintf("rule_comparison_gray_part%d.rds", k)))
  })
  tab <- do.call(rbind, lapply(parts, `[[`, "table"))
  scores <- do.call(c, lapply(parts, `[[`, "scores"))
  validations <- do.call(c, lapply(parts, `[[`, "validations"))
  # Null model of score_rule() (mean connectance of the other systems): the
  # parts were computed with the previous reference (the connectance of the
  # web itself); the Brier score does not change, only the reference and the
  # skill.
  for (nm in names(scores)) {
    w <- scores[[nm]]$webs
    for (k in seq_len(nrow(w))) {
      out <- strsplit(w$left_out[k], ", ", fixed = TRUE)[[1]]
      o <- w[!w$study_id %in% out, ]
      p <- mean(o$n_recorded / o$n_pairs)
      w$p_reference[k] <- p
      w$brier_reference[k] <- (w$n_recorded[k] * (1 - p)^2 + (w$n_pairs[k] - w$n_recorded[k]) * p^2) /
        w$n_pairs[k]
      w$skill[k] <- 1 - w$brier[k] / w$brier_reference[k]
    }
    scores[[nm]]$webs <- w
    scores[[nm]]$summary <- score_summary(w)
    tab$skill[tab$rule == nm] <- scores[[nm]]$summary$skill
  }
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
  saveRDS(cmp, file.path(out_dir, "rule_comparison_gray.rds"))
  quit(save = "no")
}
part <- if (length(args)) as.integer(args[1]) else 1L
n_parts <- if (length(args) > 1L) as.integer(args[2]) else 1L
rules <- rules[seq_along(rules) %% n_parts == part %% n_parts]
reg <- gray2015[gray2015$study_id != "12", ]
src <- setdiff(gray2015_sources$study_id[gray2015_sources$validation], "12")
groups <- list(c("13", "50"))
t0 <- Sys.time()
cmp <- compare_rules(reg, sources = src, groups = groups, rules = rules, seed = 2026)
cat("time:", format(Sys.time() - t0), "\n")
print(cmp, n = 30)
saveRDS(cmp, file.path(out_dir, if (n_parts == 1L) "rule_comparison_gray.rds" else
  sprintf("rule_comparison_gray_part%d.rds", part)))
