# Predator/prey body-mass ratio (first exploration of the size term of
# version 0.3.0, sample_metawebs(masses = )): size effect fitted to the food
# webs of Brose et al. (2005), tested on Tuesday Lake with the probabilities
# of the package (rule + taxonomic decay).
#
# Data: Brose U, et al. (2005) Body sizes of consumers and their resources.
# Ecology 86: 2545 (Ecological Archives E086-135-D1, bodysizes_2008.txt; free
# for scientific use, not redistributed here).
#
#   1. Webs with the mass of consumers and resources, without Tuesday Lake
#      and without the compilations (Ruess): every predator x animal pair of a
#      web is a link or not. Logistic models of the link with a consumer
#      effect (within-consumer size effect) and the log10 mass ratio x =
#      log10(M_consumer / M_resource): quadratic (hump-shaped), and with the
#      optimum depending on the mass of the consumer (Gravel et al. 2013).
#   2. Tuesday Lake, both years: the probabilities of the package for the
#      Gray + GloBI registry without Tuesday Lake (calibration of the 12
#      webs, from dev/validation/decay_package.R) and the size effect
#      g(x) <= 0 (0 at the optimum ratio of each consumer) added on the logit
#      scale to the predation pairs (animal resources), in several ways.
#
# Run from the root of the package after decay_package.R:
#   Rscript dev/validation/size_brose.R <bodysizes_2008.txt> <globi_registry.rds>
pkgload::load_all(".", quiet = TRUE)
args <- commandArgs(trailingOnly = TRUE)
out_dir <- file.path("dev", "validation")
b <- utils::read.delim(args[1], check.names = FALSE, stringsAsFactors = FALSE,
                       na.strings = c("-999", "-999.0"), quote = "\"", comment.char = "",
                       fileEncoding = "latin1")
names(b) <- trimws(names(b))
g <- readRDS(args[2])

# 1. Webs of Brose et al. -------------------------------------------------------
web_of <- b[["Geographic location"]]
keep_webs <- c("Lake Neuch", "Celtic Sea", "Benguela", "Broadstone", "River Frome")
web_id <- rep(NA_character_, nrow(b))
for (w in keep_webs) web_id[grepl(w, web_of, fixed = TRUE)] <- w
animal_cat <- c("invertebrate", "Invertebrate", "ectotherm vertebrate", "endotherm vertebrate")
name_of <- function(tax, common) ifelse(is.na(tax) | tax == "", common, tax)
rec <- data.frame(
  web = web_id,
  consumer = name_of(b[["Taxonomy consumer"]], b[["Common name(s) consumer"]]),
  resource = name_of(b[["Taxonomy resource"]], b[["Common name(s) resource"]]),
  mc = b[["Mean mass (g) consumer"]], mr = b[["Mean mass (g) resource"]],
  cat_c = b[["Metabolic category consumer"]], cat_r = b[["Metabolic category resource"]],
  stringsAsFactors = FALSE
)
rec <- rec[!is.na(rec$web) & !is.na(rec$consumer) & !is.na(rec$resource), ]
pairs <- do.call(rbind, lapply(split(rec, rec$web), function(d) {
  # one mass per taxon: geometric mean of its masses in both roles
  m <- c(stats::setNames(d$mc, d$consumer), stats::setNames(d$mr, d$resource))
  m <- m[!is.na(m) & m > 0]
  mass <- tapply(log10(m), names(m), mean)
  cats <- c(stats::setNames(d$cat_c, d$consumer), stats::setNames(d$cat_r, d$resource))
  animal <- unique(names(cats)[cats %in% animal_cat])
  pred <- unique(d$consumer[d$cat_r %in% animal_cat])
  animal <- intersect(animal, names(mass))
  pred <- intersect(pred, names(mass))
  p <- expand.grid(consumer = pred, resource = animal, stringsAsFactors = FALSE)
  p <- p[p$consumer != p$resource, ]
  p$y <- paste(p$resource, p$consumer) %in% paste(d$resource, d$consumer)
  p$log_mc <- unname(mass[p$consumer])
  p$x <- p$log_mc - unname(mass[p$resource])
  p$web <- d$web[1]
  p
}))
pairs$cons <- paste(pairs$web, pairs$consumer)
# consumers with at least one animal prey and one non-prey
ok <- tapply(pairs$y, pairs$cons, function(v) any(v) && !all(v))
pairs <- pairs[pairs$cons %in% names(ok)[ok], ]
print(aggregate(cbind(pairs = 1, links = y) ~ web, pairs, sum))

m1 <- stats::glm(y ~ 0 + factor(cons) + x + I(x^2), binomial, data = pairs)
m2 <- stats::glm(y ~ 0 + factor(cons) + x + I(x^2) + x:log_mc, binomial, data = pairs)
co <- function(m) stats::coef(m)[!grepl("^factor", names(stats::coef(m)))]
print(rbind(quadratic = c(co(m1), `x:log_mc` = NA), with_consumer_mass = co(m2)))
cat("AIC:", stats::AIC(m1), stats::AIC(m2), "\n")
# the same by web (quadratic)
by_web <- do.call(rbind, lapply(split(pairs, pairs$web), function(d) {
  m <- stats::glm(y ~ 0 + factor(cons) + x + I(x^2), binomial, data = d)
  c(web = d$web[1], round(co(m), 3), pairs = nrow(d), links = sum(d$y))
}))
print(by_web, quote = FALSE)

# size effect relative to the optimum of each consumer (<= 0)
size_effect <- function(beta, x, log_mc) {
  b1 <- beta[["x"]] + if ("x:log_mc" %in% names(beta)) beta[["x:log_mc"]] * log_mc else 0
  b2 <- beta[["I(x^2)"]]
  x_opt <- -b1 / (2 * b2)
  b1 * (x - x_opt) + b2 * (x^2 - x_opt^2)
}

# 2. Tuesday Lake -------------------------------------------------------------
pk <- readRDS(file.path(out_dir, "decay_package.rds"))$results$both_all$calibration
gray <- gray2015[gray2015$study_id != "12", ]
both <- combine_registries(gray, g$globi_test)
nodes <- prepare_nodes(tuesday_lake_taxa)
mw <- build_metaweb(nodes, both)
cal_co <- pk$coefficients
dec <- pk$decay
key <- function(r, c) paste(r, c, sep = "\u001f")
edges <- mw$edges[!mw$edges$cannibal, ]
cls <- edges$evidence_class
cls[cls == "beyond"] <- "class"
edges$p <- cal_co$p_new[match(cls, cal_co$evidence_class)]
edges$p[is.na(edges$p)] <- 0
cand <- mw$decay
cand$p <- dec$probabilities$p_new[decay_profile_id(cand$consumer_level, cand$resource_level,
                                                   cand$generic)]
animal_ids <- tuesday_lake_taxa$node_id[tuesday_lake_taxa$role == "consumer"]
variants <- list()
for (y in c("TL84", "TL86")) {
  com <- tuesday_lake_community[tuesday_lake_community$sample_id == y, ]
  ids <- com$node_id
  lm <- stats::setNames(log10(com$mass_kg * 1000), com$node_id)   # log10 g
  obs <- tuesday_lake_links[tuesday_lake_links$sample_id == y, ]
  obs <- obs[obs$resource %in% ids & obs$consumer %in% ids & obs$resource != obs$consumer, ]
  consumers <- unique(obs$consumer)
  pr <- expand.grid(resource = ids, consumer = consumers, stringsAsFactors = FALSE)
  pr <- pr[pr$resource != pr$consumer, ]
  k <- key(pr$resource, pr$consumer)
  pr$y <- k %in% key(obs$resource, obs$consumer)
  pr$source <- ifelse(k %in% key(edges$resource, edges$consumer), "rule",
                      ifelse(k %in% key(cand$resource, cand$consumer), "decay", "none"))
  pr$p <- 0
  pr$p[pr$source == "rule"] <- edges$p[match(k[pr$source == "rule"], key(edges$resource, edges$consumer))]
  pr$p[pr$source == "decay"] <- cand$p[match(k[pr$source == "decay"], key(cand$resource, cand$consumer))]
  pr$predation <- pr$resource %in% animal_ids
  pr$log_mc <- unname(lm[pr$consumer])
  pr$x <- pr$log_mc - unname(lm[pr$resource])
  for (mname in c("quadratic", "with_consumer_mass")) {
    beta <- co(if (mname == "quadratic") m1 else m2)
    gx <- ifelse(pr$predation, size_effect(beta, pr$x, pr$log_mc), 0)
    shift <- function(p, s) stats::plogis(stats::qlogis(pmin(pmax(p, 1e-12), 1 - 1e-12)) + s)
    p_bg <- ifelse(pr$source == "decay", shift(pr$p, gx), pr$p)
    p_all <- ifelse(pr$source != "none", shift(pr$p, gx), pr$p)
    # background renormalised to its expected number of links (one constant)
    bg <- pr$source == "decay"
    target <- sum(pr$p[bg])
    kk <- stats::uniroot(function(c) sum(shift(pr$p[bg], gx[bg] + c)) - target, c(-20, 20))$root
    p_bg_norm <- pr$p
    p_bg_norm[bg] <- shift(pr$p[bg], gx[bg] + kk)
    brier <- function(p) sum((p - pr$y)^2) / nrow(pr)
    variants[[length(variants) + 1L]] <- data.frame(
      year = y, model = mname, pairs = nrow(pr), observed = sum(pr$y),
      brier_no_size = brier(pr$p), brier_background = brier(p_bg),
      brier_background_normalised = brier(p_bg_norm), brier_all = brier(p_all),
      links_no_size = sum(pr$p), links_background = sum(p_bg), links_all = sum(p_all),
      stringsAsFactors = FALSE
    )
  }
}
variants <- do.call(rbind, variants)
print(variants, digits = 4, row.names = FALSE)
saveRDS(list(pairs_summary = aggregate(cbind(pairs = 1, links = y) ~ web, pairs, sum),
             models = list(quadratic = co(m1), with_consumer_mass = co(m2)), by_web = by_web,
             tuesday = variants), file.path(out_dir, "size_brose.rds"))
