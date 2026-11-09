# Body-mass term, validated by folds (version 0.3.0; the size term of
# sample_metawebs(masses = ) and score_metaweb(masses = )).
#
# The validation webs of the registry that have body masses are left out by
# system, as in score_rule(): their pairs take the probabilities of the class
# and decay models fitted to the other systems (the rule infers the link, or
# the taxonomic decay gives the pair its probability), and the size term is
# fitted to the webs with masses of the other systems. Tuesday Lake (study
# 12) is never used to fit anything and is scored with the fits of the 12
# webs.
#
# Size term: within-consumer effect of x = log10(M_consumer / M_resource) on
# the links of each predator with the animals of its web (consumer fixed
# effects; quadratic in x, and with the optimum depending on log10
# M_consumer), relative to the optimum of each consumer (g(x) <= 0), added on
# the logit scale to the predation pairs (animal resources) with the masses
# of both nodes. Variants: to the pairs of the decay only, keeping or not
# their expected number of links in the web; to every predation pair; and
# the same with a one-sided effect (only prey larger than the optimum of the
# consumer are penalised). Every validation web whose consumers have masses
# (by name, or by genus) for at least 60% of them is scored.
#
# Masses (and links for the size term): Brose U et al. (2005), Body sizes of
# consumers and their resources, Ecology 86: 2545 (Ecological Archives
# E086-135, bodysizes_2008.txt); the cheddar package (Hudson et al. 2013;
# data files BroadstoneStream, Millstream, pHWebs, TL84 and TL86 of its
# source). Neither is redistributed here.
#
# Run from the root of the package (about 30 minutes; then
# dev/validation/size_lambda.R and dev/validation/size_sources.R):
#   Rscript dev/validation/size_folds.R <bodysizes_2008.txt> <cheddar data folder> \
#     <globi_registry.rds>
pkgload::load_all(".", quiet = TRUE)
args <- commandArgs(trailingOnly = TRUE)
out_dir <- file.path("dev", "validation")
set.seed(2026)

# 1. Webs with body masses ---------------------------------------------------------
clean <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x <- gsub("\\s*\\(.*\\)$", "", x)
  x <- gsub("\\s+(sp|spp|agg|sp\\.|spp\\.|agg\\.)\\.?$", "", x)
  x[!is.na(x) & x == ""] <- NA_character_
  x
}
genus_of <- function(x) sub(" .*$", "", x)
animal_cat <- c("invertebrate", "Invertebrate", "ectotherm vertebrate", "endotherm vertebrate",
                "vert.ecto", "vert.endo")

# Brose et al. (2005): one web per location with masses
b <- utils::read.delim(args[1], check.names = FALSE, stringsAsFactors = FALSE,
                       na.strings = c("-999", "-999.0"), quote = "\"", comment.char = "",
                       fileEncoding = "latin1")
names(b) <- trimws(names(b))
brose_webs <- c(caricaie = "Lake Neuch", celtic = "Celtic Sea", benguela = "Benguela",
                broadstone = "Broadstone", frome = "River Frome", tuesday = "Tuesday Lake")
name_of <- function(tax, common) ifelse(is.na(tax) | tax == "", common, tax)
webs <- list()
for (sys in names(brose_webs)) {
  d <- b[grepl(brose_webs[[sys]], b[["Geographic location"]], fixed = TRUE), ]
  con <- name_of(d[["Taxonomy consumer"]], d[["Common name(s) consumer"]])
  res <- name_of(d[["Taxonomy resource"]], d[["Common name(s) resource"]])
  m <- c(stats::setNames(d[["Mean mass (g) consumer"]], con),
         stats::setNames(d[["Mean mass (g) resource"]], res))
  m <- m[!is.na(m) & m > 0 & !is.na(names(m))]
  cats <- c(stats::setNames(d[["Metabolic category consumer"]], con),
            stats::setNames(d[["Metabolic category resource"]], res))
  nodes <- unique(c(con, res))
  nodes <- nodes[!is.na(nodes)]
  webs[[paste0("brose_", sys)]] <- list(
    systems = sys,
    nodes = data.frame(name = nodes, log_mass = unname(tapply(log10(m), names(m), mean)[nodes]),
                       animal = nodes %in% names(cats)[cats %in% animal_cat],
                       genus = NA_character_, stringsAsFactors = FALSE),
    links = data.frame(resource = res, consumer = con, stringsAsFactors = FALSE)
  )
}
# cheddar: masses in mg (Tuesday Lake in kg)
ld <- function(f) {
  e <- new.env()
  load(file.path(args[2], paste0(f, ".RData")), envir = e)
  get(ls(e)[1], e)
}
from_cheddar <- function(x, systems) {
  n <- x$nodes
  to_g <- if (identical(x$properties$M.units, "kg")) 1000 else 1e-3
  list(systems = systems,
       nodes = data.frame(name = n$node, log_mass = log10(n$M * to_g),
                          animal = n$category %in% animal_cat,
                          genus = if ("genus" %in% names(n)) n$genus else NA_character_,
                          stringsAsFactors = FALSE),
       links = x$trophic.links[c("resource", "consumer")])
}
webs$cheddar_broadstone <- from_cheddar(ld("BroadstoneStream"), "broadstone")
ms <- ld("Millstream")
for (k in names(ms)) webs[[paste0("cheddar_millstream_", k)]] <- from_cheddar(ms[[k]], "frome")
ph <- ld("pHWebs")
for (k in names(ph)) {
  sys <- c("layer", if (k == "Broadstone") "broadstone", if (k %in% c("Mill Stream", "Bere Stream")) "frome")
  webs[[paste0("cheddar_ph_", gsub("[^A-Za-z]", "", k))]] <- from_cheddar(ph[[k]], sys)
}
webs$cheddar_tl84 <- from_cheddar(ld("TL84"), "tuesday")
webs$cheddar_tl86 <- from_cheddar(ld("TL86"), "tuesday")

# Predation pairs of each web: every predator (a node with an animal prey)
# with every animal with a mass, except itself
size_pairs <- function(w) {
  n <- w$nodes[!is.na(w$nodes$log_mass), ]
  animals <- n$name[n$animal]
  l <- w$links[w$links$resource %in% animals & w$links$resource != w$links$consumer, ]
  pred <- intersect(unique(l$consumer), n$name)
  if (!length(pred)) return(NULL)
  p <- expand.grid(consumer = pred, resource = animals, stringsAsFactors = FALSE)
  p <- p[p$consumer != p$resource, ]
  p$y <- paste(p$resource, p$consumer) %in% paste(l$resource, l$consumer)
  lm <- stats::setNames(n$log_mass, n$name)
  p$log_mc <- unname(lm[p$consumer])
  p$x <- p$log_mc - unname(lm[p$resource])
  p
}
train_pairs <- lapply(webs, size_pairs)

# Mass of every name: the mean over the webs where it is a node (log10 g)
mass_rows <- do.call(rbind, lapply(webs, function(w) {
  data.frame(name = clean(w$nodes$name), genus = clean(ifelse(is.na(w$nodes$genus) | w$nodes$genus == "",
                                                            genus_of(w$nodes$name), w$nodes$genus)),
             log_mass = w$nodes$log_mass, stringsAsFactors = FALSE)
}))
mass_rows <- mass_rows[!is.na(mass_rows$log_mass) & is.finite(mass_rows$log_mass), ]
by_name <- tapply(mass_rows$log_mass, mass_rows$name, mean)
by_genus <- tapply(mass_rows$log_mass, mass_rows$genus, mean)
mass_lookup <- function(taxon, genus) {
  k <- clean(taxon)
  out <- unname(by_name[k])
  g <- clean(ifelse(is.na(genus), genus_of(k), genus))
  miss <- is.na(out)
  out[miss] <- unname(by_genus[g[miss]])
  out
}

# Size term fitted to the webs that share no system with `exclude`
size_fits <- new.env()
size_fit <- function(exclude) {
  key <- paste(sort(exclude), collapse = "+")
  if (!is.null(size_fits[[key]])) return(size_fits[[key]])
  use <- !vapply(webs, function(w) any(w$systems %in% exclude), logical(1))
  d <- do.call(rbind, Map(function(p, nm) if (!is.null(p)) cbind(p, web = nm), train_pairs[use],
                          names(webs)[use]))
  d$cons <- paste(d$web, d$consumer)
  ok <- tapply(d$y, d$cons, function(v) any(v) && !all(v))
  d <- d[d$cons %in% names(ok)[ok], ]
  co <- function(m) stats::coef(m)[!grepl("^factor", names(stats::coef(m)))]
  fit <- list(
    quadratic = co(stats::glm(y ~ 0 + factor(cons) + x + I(x^2), stats::binomial, data = d)),
    consumer_mass = co(stats::glm(y ~ 0 + factor(cons) + x + I(x^2) + x:log_mc, stats::binomial,
                                  data = d)),
    pairs = nrow(d), links = sum(d$y), webs = sum(use)
  )
  assign(key, fit, envir = size_fits)
  fit
}
size_effect <- function(beta, x, log_mc, one_sided = FALSE) {
  b1 <- beta[["x"]] + if ("x:log_mc" %in% names(beta)) beta[["x:log_mc"]] * log_mc else 0
  b2 <- beta[["I(x^2)"]]
  x_opt <- -b1 / (2 * b2)
  g <- b1 * (x - x_opt) + b2 * (x^2 - x_opt^2)
  if (one_sided) g[!is.na(x) & x > x_opt] <- 0
  g
}

# 2. Validation webs, left out by system ----------------------------------------
g <- readRDS(args[3])
gray <- gray2015[gray2015$study_id != "12", ]
registries <- list(gray = list(registry = gray, exclusions = NULL),
                   both_all = list(registry = combine_registries(gray, g$globi_test),
                                   exclusions = g$exclusions))
rule <- matching_rule()
src <- setdiff(gray2015_sources$study_id[gray2015_sources$validation], "12")
groups <- list(c("13", "50"))
# systems of the validation webs with webs of their own among those with
# body masses (left out of the fit of the size term); the others share no
# system with them
web_systems <- list(`13` = "frome", `50` = "frome", `46` = "broadstone", `49` = "broadstone",
                    `27` = c("layer", "broadstone", "frome"), `12` = "tuesday")
systems_of <- function(w) if (is.null(web_systems[[w]])) character(0) else web_systems[[w]]
key <- function(r, c) paste(r, c, sep = "\u001f")

pairs_of <- function(reg, w, train) {
  vn <- registry_nodes(reg, w, rule)
  nodes <- vn$nodes
  mw <- metaweb_core(nodes, train, check_manual_links(NULL, nodes$node_id), prefilter = "both")
  studied <- nodes$node_id[vn$studied]
  e <- mw$edges[!mw$edges$cannibal & mw$edges$consumer %in% studied, ]
  p <- expand.grid(resource = nodes$node_id, consumer = studied, stringsAsFactors = FALSE)
  p <- p[p$resource != p$consumer, ]
  k <- key(p$resource, p$consumer)
  p$y <- k %in% key(vn$records$resource, vn$records$consumer)
  p$inferred <- k %in% key(e$resource, e$consumer)
  p$evidence_class <- e$evidence_class[match(k, key(e$resource, e$consumer))]
  ev <- decay_evidence(decay_profiles(nodes, train, rule, consumers = studied),
                       p$resource[!p$inferred], p$consumer[!p$inferred])
  p$consumer_level <- p$resource_level <- NA_character_
  p$generic <- NA
  p$consumer_level[!p$inferred] <- ev$consumer_level
  p$resource_level[!p$inferred] <- ev$resource_level
  p$generic[!p$inferred] <- ev$generic
  role <- stats::setNames(nodes$role, nodes$node_id)
  p$animal <- unname(role[p$resource]) == "consumer"
  lm <- stats::setNames(mass_lookup(nodes$taxon, nodes$genus), nodes$node_id)
  p$log_mc <- unname(lm[p$consumer])
  p$x <- p$log_mc - unname(lm[p$resource])
  attr(p, "mass_coverage") <- mean(!is.na(lm[nodes$node_id[nodes$role == "consumer"]]))
  p
}

score_web <- function(p, fit_size) {
  pred <- p$animal & !is.na(p$x)
  logit <- stats::qlogis(pmin(pmax(p$p, 1e-12), 1 - 1e-12))
  shift <- function(s, idx, c = 0) stats::plogis(logit[idx] + s[idx] + c)
  brier <- function(q) sum((q - p$y)^2) / nrow(p)
  out <- list(none = brier(p$p))
  for (m in c("quadratic", "consumer_mass", "one_sided")) {
    beta <- fit_size[[if (m == "one_sided") "quadratic" else m]]
    s <- ifelse(pred, size_effect(beta, p$x, p$log_mc, one_sided = m == "one_sided"), 0)
    bg <- !p$inferred
    q <- p$p
    q[bg] <- shift(s, bg)
    out[[paste0("decay_", m)]] <- brier(q)
    target <- sum(p$p[bg])
    k <- if (any(bg)) {
      stats::uniroot(function(c) sum(shift(s, bg, c)) - target, c(-30, 30))$root
    } else {
      0
    }
    q[bg] <- shift(s, bg, k)
    out[[paste0("decay_kept_", m)]] <- brier(q)
    q <- p$p
    q[pred] <- shift(s, pred)
    out[[paste0("all_", m)]] <- brier(q)
  }
  out
}

results <- list()
for (rn in names(registries)) {
  st <- registries[[rn]]
  reg <- st$registry
  v <- validate_rule(reg, src, groups = groups, rule = rule, exclusions = st$exclusions, quiet = TRUE)
  left_out <- stats::setNames(strsplit(v$webs$left_out, ", ", fixed = TRUE), v$webs$study_id)
  # the probabilities of a web: fits of the class and decay models to the
  # other systems (the 12 webs for Tuesday Lake)
  fold_p <- function(out_ids) {
    tr <- v$counts[!v$counts$study_id %in% out_ids, ]
    td <- v$decay[!v$decay$study_id %in% out_ids, ]
    fc <- fit_calibration_bayes(tr, 1.5, 2.5, 2L, 500L, 500L, 20L, seed = NULL)
    fd <- fit_decay(td, n_chains = 2L, n_iter = 500L, n_warmup = 500L, thin = 2L)
    list(class = stats::setNames(fc$coefficients$p_new, fc$coefficients$evidence_class), decay = fd)
  }
  tests <- c(src, "12")
  rows <- list()
  pair_data <- list()
  for (w in tests) {
    out_ids <- if (w == "12") "12" else left_out[[w]]
    extra <- if (is.null(st$exclusions)) character(0) else st$exclusions$study_id[st$exclusions$web == w]
    # Tuesday Lake: its records are needed for its nodes and links (study 12),
    # and are left out of the training registry
    reg_w <- if (w != "12") reg else if (rn == "gray") gray2015 else
      combine_registries(gray2015, g$globi_test)
    train <- reg_w[!reg_w$study_id %in% c(out_ids, extra), ]
    p <- pairs_of(reg_w, w, train)
    if (attr(p, "mass_coverage") < 0.6) {
      message(rn, " web ", w, ": masses for ", round(100 * attr(p, "mass_coverage")),
              "% of the consumers, not scored")
      next
    }
    fp <- fold_p(if (w == "12") character(0) else out_ids)
    cls <- p$evidence_class
    cls[!is.na(cls) & cls == "beyond"] <- "class"
    p$p <- ifelse(p$inferred, unname(fp$class[cls]), 0)
    p$p[is.na(p$p)] <- 0
    p$p[!p$inferred] <- decay_p_new_of(fp$decay, p[!p$inferred, ],
                                       n_nodes = length(unique(p$resource)))
    fs <- size_fit(c("tuesday", systems_of(w)))
    sc <- score_web(p, fs)
    rows[[w]] <- cbind(data.frame(registry = rn, web = w, pairs = nrow(p), recorded = sum(p$y),
                                  mass_coverage = round(attr(p, "mass_coverage"), 2),
                                  predation_pairs = sum(p$animal & !is.na(p$x)),
                                  predation_links = sum(p$y & p$animal & !is.na(p$x)),
                                  size_training_pairs = fs$pairs, stringsAsFactors = FALSE),
                       as.data.frame(sc))
    pair_data[[w]] <- p[c("resource", "consumer", "y", "inferred", "p", "animal", "x", "log_mc")]
    message(rn, " web ", w, " done")
  }
  results[[rn]] <- do.call(rbind, rows)
  attr(results[[rn]], "pairs") <- pair_data
}
pairs_all <- lapply(results, attr, "pairs")
tab <- do.call(rbind, results)
rownames(tab) <- NULL
options(width = 200)
print(tab, digits = 4)
variants <- setdiff(names(tab), c("registry", "web", "pairs", "recorded", "mass_coverage",
                                  "predation_pairs", "predation_links", "size_training_pairs",
                                  "none"))
summ <- do.call(rbind, lapply(split(tab[tab$web != "12", ], tab$registry[tab$web != "12"]), function(d) {
  do.call(rbind, lapply(variants, function(vn) {
    dd <- d[[vn]] - d$none
    data.frame(registry = d$registry[1], variant = vn, webs = nrow(d), diff = mean(dd),
               diff_se = stats::sd(dd) / sqrt(nrow(d)), webs_better = sum(dd < 0))
  }))
}))
print(summ, digits = 3, row.names = FALSE)
size_coefs <- lapply(as.list(size_fits), function(f) f[c("quadratic", "consumer_mass", "pairs", "links", "webs")])
saveRDS(list(webs = tab, summary = summ, size_fits = size_coefs, pairs = pairs_all),
        file.path(out_dir, "size_folds.rds"))
