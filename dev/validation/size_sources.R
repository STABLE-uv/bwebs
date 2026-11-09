# Body-mass term: which webs to fit the size curve to, and the weights of
# size_model (version 0.3.0). The curves of the folds of
# dev/validation/size_folds.R are refitted with every web with masses (as
# there) and without Grand Cariçaie (Brose et al. 2005; a wetland web of
# arthropods that holds most of the pairs), and the pairs of the validation
# webs scored by size_folds.R are scored again with a grid of weights of the
# size term on the inferred links and on the pairs of the decay. The
# configuration with the best mean score over the webs is chosen for each
# registry and for both together; the same choice made without each web
# (nested) estimates how well the chosen configuration does on a new web.
#
# Run from the root of the package after size_folds.R:
#   Rscript dev/validation/size_sources.R <bodysizes_2008.txt> <cheddar data folder>
pkgload::load_all(".", quiet = TRUE)
args <- commandArgs(trailingOnly = TRUE)
out_dir <- file.path("dev", "validation")
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


r <- readRDS(file.path(out_dir, "size_folds.rds"))
web_systems <- list(`13` = "frome", `50` = "frome", `46` = "broadstone", `49` = "broadstone",
                    `27` = c("layer", "broadstone", "frome"), `12` = "tuesday")
systems_of <- function(w) if (is.null(web_systems[[w]])) character(0) else web_systems[[w]]
weights <- expand.grid(rule = c(0, 0.25, 0.5), decay = c(1, 3, 10, 30))
rows <- list()
for (curve in c("all webs", "without Grand Cariçaie")) {
  extra <- if (curve == "all webs") character(0) else "caricaie"
  for (rn in names(r$pairs)) for (w in names(r$pairs[[rn]])) {
    p <- r$pairs[[rn]][[w]]
    beta <- size_fit(c("tuesday", extra, systems_of(w)))$quadratic
    g <- ifelse(p$animal & !is.na(p$x), size_effect(beta, p$x, p$log_mc), 0)
    g[is.na(g)] <- 0
    lg <- stats::qlogis(pmin(pmax(p$p, 1e-12), 1 - 1e-12))
    for (k in seq_len(nrow(weights))) {
      q <- stats::plogis(lg + ifelse(p$inferred, weights$rule[k], weights$decay[k]) * g)
      rows[[length(rows) + 1L]] <- data.frame(curve = curve, registry = rn, web = w,
                                              rule = weights$rule[k], decay = weights$decay[k],
                                              none = sum((p$p - p$y)^2) / nrow(p),
                                              brier = sum((q - p$y)^2) / nrow(p))
    }
  }
}
tab <- do.call(rbind, rows)
tab$diff <- tab$brier - tab$none
cv <- tab[tab$web != "12", ]
summ <- do.call(rbind, lapply(split(cv, list(cv$curve, cv$registry, cv$rule, cv$decay), drop = TRUE),
                              function(d) {
  data.frame(curve = d$curve[1], registry = d$registry[1], rule = d$rule[1], decay = d$decay[1],
             diff = mean(d$diff), diff_se = stats::sd(d$diff) / sqrt(nrow(d)),
             webs_better = sum(d$diff < 0),
             tuesday = tab$brier[tab$web == "12" & tab$curve == d$curve[1] &
                                   tab$registry == d$registry[1] & tab$rule == d$rule[1] &
                                   tab$decay == d$decay[1]])
}))
summ <- summ[order(summ$registry, summ$diff), ]
rownames(summ) <- NULL
print(summ, digits = 3)
curves <- lapply(as.list(size_fits), function(f) f$quadratic)
print(do.call(rbind, curves), digits = 3)

# Nested selection: for each web left out, the curve and weights with the
# best mean score over the other webs, and its score on that web
cfg <- unique(tab[c("curve", "rule", "decay")])
cfg$id <- seq_len(nrow(cfg))
tab <- merge(tab, cfg)
nested <- do.call(rbind, lapply(unique(cv$registry), function(rn) {
  d <- tab[tab$registry == rn & tab$web != "12", ]
  M <- tapply(d$diff, list(d$web, d$id), mean)
  out <- do.call(rbind, lapply(rownames(M), function(w) {
    k <- as.integer(names(which.min(colMeans(M[rownames(M) != w, , drop = FALSE]))))
    cbind(data.frame(registry = rn, web = w, diff = M[w, as.character(k)]),
          cfg[cfg$id == k, c("curve", "rule", "decay")])
  }))
  k_all <- as.integer(names(which.min(colMeans(M))))
  message(sprintf("%s: nested mean %.5f (se %.5f), %d of %d webs better; chosen with every web: %s, rule %.2f, decay %g",
                  rn, mean(out$diff), stats::sd(out$diff) / sqrt(nrow(out)), sum(out$diff < 0),
                  nrow(out), cfg$curve[k_all], cfg$rule[k_all], cfg$decay[k_all]))
  out
}))
print(nested, digits = 3, row.names = FALSE)

# Both registries together (the size model does not depend on the
# calibration): the mean difference of each web over the two registries,
# the configuration with the best mean over the webs, and the same choice
# made without each web left out (nested)
d <- tab[tab$web != "12", ]
M <- tapply(d$diff, list(d$web, d$id), mean)
both <- cfg[match(as.integer(colnames(M)), cfg$id), c("curve", "rule", "decay")]
both$diff <- colMeans(M)
both$diff_se <- apply(M, 2L, stats::sd) / sqrt(nrow(M))
both$webs_better <- colSums(M < 0)
for (rn in unique(d$registry)) {
  dr <- d[d$registry == rn, ]
  both[[rn]] <- tapply(dr$diff, dr$id, mean)[colnames(M)]
}
both <- both[order(both$diff), ]
rownames(both) <- NULL
print(both, digits = 3)
nested_both <- do.call(rbind, lapply(rownames(M), function(w) {
  k <- as.integer(names(which.min(colMeans(M[rownames(M) != w, , drop = FALSE]))))
  out <- data.frame(web = w, diff = M[w, as.character(k)])
  for (rn in unique(d$registry)) {
    out[[rn]] <- d$diff[d$registry == rn & d$web == w & d$id == k]
  }
  cbind(out, cfg[cfg$id == k, c("curve", "rule", "decay")])
}))
k_all <- as.integer(names(which.min(colMeans(M))))
message(sprintf("both registries: nested mean %.5f (se %.5f), %d of %d webs better; chosen with every web: %s, rule %.2f, decay %g",
                mean(nested_both$diff), stats::sd(nested_both$diff) / sqrt(nrow(nested_both)),
                sum(nested_both$diff < 0), nrow(nested_both), cfg$curve[k_all], cfg$rule[k_all],
                cfg$decay[k_all]))
print(nested_both, digits = 3, row.names = FALSE)
# The same nested choice with the curve fixed (every web, or without Grand
# Cariçaie) and only the weights chosen, against choosing the curve too
by_curve <- do.call(rbind, lapply(list("both registries", "gray", "both_all"), function(rn) {
  dd <- if (rn == "both registries") d else d[d$registry == rn, ]
  do.call(rbind, lapply(c("all webs", "without Grand Cariçaie", "either"), function(cv) {
    ids <- if (cv == "either") cfg$id else cfg$id[cfg$curve == cv]
    M <- tapply(dd$diff, list(dd$web, dd$id), mean)[, as.character(ids), drop = FALSE]
    nest <- vapply(rownames(M), function(w) {
      M[w, which.min(colMeans(M[rownames(M) != w, , drop = FALSE]))]
    }, numeric(1))
    k <- as.integer(colnames(M)[which.min(colMeans(M))])
    data.frame(registry = rn, curve = cv, cv = min(colMeans(M)), nested = mean(nest),
               nested_se = stats::sd(nest) / sqrt(length(nest)), nested_better = sum(nest < 0),
               chosen_curve = cfg$curve[cfg$id == k], chosen_rule = cfg$rule[cfg$id == k],
               chosen_decay = cfg$decay[cfg$id == k], stringsAsFactors = FALSE)
  }))
}))
print(by_curve, digits = 3, row.names = FALSE)
saveRDS(list(summary = summ, webs = tab, curves = curves, nested = nested, both = both,
             nested_both = nested_both, by_curve = by_curve),
        file.path(out_dir, "size_sources.rds"))
