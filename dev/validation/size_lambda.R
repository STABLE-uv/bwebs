# Body-mass term: a weaker effect on the links inferred by the rule
# (version 0.3.0; for the weights of size_model). From the pairs scored by
# dev/validation/size_folds.R: the size term g(x) is added in full to the
# pairs of the decay and multiplied by lambda for the inferred links
# (lambda = 0: decay only; lambda = 1: every predation pair).
#
# Run from the root of the package after size_folds.R:
#   Rscript dev/validation/size_lambda.R
r <- readRDS(file.path("dev", "validation", "size_folds.rds"))
size_effect <- function(beta, x) {
  b1 <- beta[["x"]]; b2 <- beta[["I(x^2)"]]; x_opt <- -b1 / (2 * b2)
  b1 * (x - x_opt) + b2 * (x^2 - x_opt^2)
}
web_systems <- list(`13` = "frome", `50` = "frome", `46` = "broadstone", `49` = "broadstone",
                    `27` = c("layer", "broadstone", "frome"), `12` = "tuesday")
fit_for <- function(w) {
  key <- paste(sort(c("tuesday", if (is.null(web_systems[[w]])) character(0) else web_systems[[w]])), collapse = "+")
  r$size_fits[[key]]$quadratic
}
out <- list()
for (rn in names(r$pairs)) for (w in names(r$pairs[[rn]])) {
  p <- r$pairs[[rn]][[w]]
  g <- ifelse(p$animal & !is.na(p$x), size_effect(fit_for(w), p$x), 0)
  lg <- stats::qlogis(pmin(pmax(p$p, 1e-12), 1 - 1e-12))
  b <- function(q) sum((q - p$y)^2) / nrow(p)
  row <- c(none = b(p$p))
  for (lam in c(0, 0.25, 0.5, 0.75, 1)) {
    q <- stats::plogis(lg + ifelse(p$inferred, lam, 1) * g)
    row[paste0("rule_lambda_", lam)] <- b(q)
  }
  out[[paste(rn, w)]] <- data.frame(registry = rn, web = w, t(row))
}
tab <- do.call(rbind, out); rownames(tab) <- NULL
print(tab, digits = 4)
for (rn in unique(tab$registry)) {
  d <- tab[tab$registry == rn & tab$web != "12", ]
  cat("\n", rn, ": mean difference from no size (se), webs better\n")
  for (v in grep("lambda", names(d), value = TRUE)) {
    dd <- d[[v]] - d$none
    cat(sprintf("  %-18s %8.4f (%.4f) %d/%d\n", v, mean(dd), sd(dd) / sqrt(length(dd)), sum(dd < 0), length(dd)))
  }
}

# Grid of weights on the inferred links (lambda) and on the pairs of the
# decay, and the choice of the weights by the validation by folds: the
# weights with the best mean difference over the webs, chosen with every
# web and without each web in turn (nested), over the whole grid and over
# the grid with lambda <= 0.5. Both registries together: the mean of each
# web over the two registries.
grid <- expand.grid(rule = c(0, 0.25, 0.5, 0.75, 1), decay = c(1, 3, 10, 30))
rows <- list()
for (rn in names(r$pairs)) for (w in setdiff(names(r$pairs[[rn]]), "12")) {
  p <- r$pairs[[rn]][[w]]
  g <- ifelse(p$animal & !is.na(p$x), size_effect(fit_for(w), p$x), 0)
  g[is.na(g)] <- 0
  lg <- stats::qlogis(pmin(pmax(p$p, 1e-12), 1 - 1e-12))
  none <- sum((p$p - p$y)^2) / nrow(p)
  for (k in seq_len(nrow(grid))) {
    q <- stats::plogis(lg + ifelse(p$inferred, grid$rule[k], grid$decay[k]) * g)
    rows[[length(rows) + 1L]] <- data.frame(registry = rn, web = w, id = k,
                                            diff = sum((q - p$y)^2) / nrow(p) - none)
  }
}
gt <- do.call(rbind, rows)
choose <- function(M, ids) {
  best <- function(rows_used) {
    ids[which.min(colMeans(M[rows_used, as.character(ids), drop = FALSE]))]
  }
  nested <- vapply(rownames(M), function(w) {
    M[w, as.character(best(rownames(M) != w))]
  }, numeric(1))
  k <- best(rep(TRUE, nrow(M)))
  data.frame(rule = grid$rule[k], decay = grid$decay[k], cv = mean(M[, as.character(k)]),
             nested = mean(nested), nested_se = stats::sd(nested) / sqrt(length(nested)),
             nested_better = sum(nested < 0))
}
sel <- list()
for (rn in c(names(r$pairs), "both registries")) {
  d <- if (rn == "both registries") gt else gt[gt$registry == rn, ]
  M <- tapply(d$diff, list(d$web, d$id), mean)
  for (sub in c("lambda <= 1", "lambda <= 0.5")) {
    ids <- if (sub == "lambda <= 1") seq_len(nrow(grid)) else which(grid$rule <= 0.5)
    sel[[paste(rn, sub)]] <- cbind(data.frame(registry = rn, grid = sub), choose(M, ids))
  }
}
sel <- do.call(rbind, sel)
rownames(sel) <- NULL
cat("\nChoice of the weights (curve fitted to every web with masses of the other systems)\n")
print(sel, digits = 3)
