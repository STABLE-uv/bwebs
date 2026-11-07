# How much does the calibration change when the validation webs of the same
# system are left out together? And how much does the example of the README
# change with each option (and with system = "average")?
#
# Sign of shared systems: share of the links of each validation web that are
# also recorded in another study of the registry.
#
# Run from the root of the package: Rscript dev/validation/validation_groups.R
# (about 6 minutes).
pkgload::load_all(".", quiet = TRUE)

r <- gray2015
val <- gray2015_sources$study_id[gray2015_sources$validation]
links <- tapply(seq_len(nrow(r)), r$study_id, function(i) {
  unique(paste(tolower(r$resource[i]), tolower(r$consumer[i])))
})
shared <- do.call(rbind, lapply(val, function(w) {
  o <- setdiff(names(links), w)
  s <- vapply(o, function(x) mean(links[[w]] %in% links[[x]]), numeric(1))
  k <- order(-s)[1:2]
  data.frame(web = w, n_links = length(links[[w]]), other_1 = o[k[1]], share_1 = round(s[k[1]], 2),
             other_2 = o[k[2]], share_2 = round(s[k[2]], 2), row.names = NULL)
}))
print(shared, row.names = FALSE)

options_groups <- list(
  alone = NULL,
  shared_links = list(c("13", "50")),
  shared_links_and_authors = list(c("13", "50"), c("46", "49"), c("15", "59"))
)
cals <- lapply(options_groups, function(g) calibrate_links(gray2015, sources = val, groups = g,
                                                            quiet = TRUE))
coef_table <- do.call(rbind, lapply(names(cals), function(nm) {
  co <- cals[[nm]]$coefficients
  data.frame(groups = nm, class = co$evidence_class, p_average = round(co$p_average, 3),
             sigma = round(cals[[nm]]$sigma, 2))
}))
print(reshape(coef_table, idvar = c("class"), timevar = "groups", direction = "wide"),
      row.names = FALSE)

nodes <- prepare_nodes(example_taxa)
mw <- build_metaweb(nodes, gray2015)
example_for <- function(cal, system) {
  dr <- sample_metawebs(mw, cal, n_draws = 1000, system = system, seed = 2026)
  ind <- calc_indicators(dr, example_community, permanent = c("cpom", "fpom"))
  s <- summarize_indicators(ind, indicators = c("S", "C", "TL_mean"), samples = example_samples,
                            by = "treatment", na_rm = TRUE)
  k <- contrast_indicators(ind, example_samples, by = "treatment", reference = "control",
                           indicators = c("S", "C", "TL_mean"), na_rm = TRUE)
  list(summary = s, contrast = k, links = stats::median(colSums(dr$z)))
}
runs <- list(
  alone_random = example_for(cals$alone, "random"),
  alone_average = example_for(cals$alone, "average"),
  shared_links_random = example_for(cals$shared_links, "random"),
  shared_links_and_authors_random = example_for(cals$shared_links_and_authors, "random")
)
fmt <- function(e, l, u) sprintf("%.3f [%.3f; %.3f]", e, l, u)
example_table <- do.call(rbind, lapply(names(runs), function(nm) {
  s <- runs[[nm]]$summary
  k <- runs[[nm]]$contrast
  data.frame(
    option = nm, median_links = runs[[nm]]$links,
    C_control = fmt(s$estimate[s$treatment == "control" & s$indicator == "C"],
                    s$lower[s$treatment == "control" & s$indicator == "C"],
                    s$upper[s$treatment == "control" & s$indicator == "C"]),
    TL_control = fmt(s$estimate[s$treatment == "control" & s$indicator == "TL_mean"],
                     s$lower[s$treatment == "control" & s$indicator == "TL_mean"],
                     s$upper[s$treatment == "control" & s$indicator == "TL_mean"]),
    dC = fmt(k$estimate[k$indicator == "C"], k$lower[k$indicator == "C"], k$upper[k$indicator == "C"]),
    P_dC_pos = k$prob_positive[k$indicator == "C"],
    dTL = fmt(k$estimate[k$indicator == "TL_mean"], k$lower[k$indicator == "TL_mean"],
              k$upper[k$indicator == "TL_mean"]),
    P_dTL_pos = k$prob_positive[k$indicator == "TL_mean"]
  )
}))
print(example_table, row.names = FALSE)
saveRDS(list(shared = shared, coefficients = coef_table, example = example_table,
             webs = lapply(cals, function(x) x$webs)),
        file.path("dev", "validation", "validation_groups.rds"))
