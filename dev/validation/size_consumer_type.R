# Body-mass term by kind of consumer (version 0.3.0): one curve for every
# predator (size_model), the size term only for invertebrate predators, or
# one curve for invertebrate and another for vertebrate predators (fish), each
# fitted to the webs with masses of the other systems, as in
# dev/validation/size_sources.R. The pairs of the seven validation webs with
# masses scored by dev/validation/size_folds.R are scored again with the
# weights of size_model (0.5 on the inferred links, 10 on the pairs of the
# decay) and with the grid of weights of size_sources.R; Tuesday Lake (study
# 12) is scored with the fits without it. Vertebrate consumers: those of
# category vert.ecto in the registry (fish).
#
# Run from the root of the package after size_folds.R:
#   Rscript dev/validation/size_consumer_type.R <bodysizes_2008.txt> <cheddar data folder>
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
vertebrate_cat <- c("ectotherm vertebrate", "endotherm vertebrate", "vert.ecto", "vert.endo")

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
                       vertebrate = nodes %in% names(cats)[cats %in% vertebrate_cat],
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
                          vertebrate = n$category %in% vertebrate_cat,
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


# Predation pairs of each web, with the kind of consumer
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
  p$x <- unname(lm[p$consumer]) - unname(lm[p$resource])
  p$vertebrate <- p$consumer %in% n$name[n$vertebrate]
  p
}
train_pairs <- lapply(webs, size_pairs)

# Curves fitted to the webs that share no system with `exclude`: every
# predator, invertebrate predators and vertebrate predators
fits <- new.env()
fit_curves <- function(exclude) {
  key <- paste(sort(exclude), collapse = "+")
  if (!is.null(fits[[key]])) return(fits[[key]])
  use <- !vapply(webs, function(w) any(w$systems %in% exclude), logical(1))
  d <- do.call(rbind, Map(function(p, nm) if (!is.null(p)) cbind(p, web = nm), train_pairs[use],
                          names(webs)[use]))
  d$cons <- paste(d$web, d$consumer)
  ok <- tapply(d$y, d$cons, function(v) any(v) && !all(v))
  d <- d[d$cons %in% names(ok)[ok], ]
  co <- function(dd) {
    m <- stats::glm(y ~ 0 + factor(cons) + x + I(x^2), stats::binomial, data = dd)
    stats::coef(m)[c("x", "I(x^2)")]
  }
  out <- list(all = co(d), invertebrate = co(d[!d$vertebrate, ]), vertebrate = co(d[d$vertebrate, ]),
              pairs = c(invertebrate = sum(!d$vertebrate), vertebrate = sum(d$vertebrate)),
              webs_vertebrate = length(unique(d$web[d$vertebrate])))
  assign(key, out, envir = fits)
  out
}
g_of <- function(beta, x) {
  b1 <- beta[["x"]]
  b2 <- beta[["I(x^2)"]]
  x_opt <- -b1 / (2 * b2)
  g <- b1 * (x - x_opt) + b2 * (x^2 - x_opt^2)
  g[is.na(g)] <- 0
  g
}

r <- readRDS(file.path(out_dir, "size_folds.rds"))
web_systems <- list(`13` = "frome", `50` = "frome", `46` = "broadstone", `49` = "broadstone",
                    `27` = c("layer", "broadstone", "frome"), `12` = "tuesday")
systems_of <- function(w) if (is.null(web_systems[[w]])) character(0) else web_systems[[w]]
vert_names <- unique(gray2015$consumer[gray2015$con_category == "vert.ecto"])
weights <- expand.grid(rule = c(0, 0.25, 0.5), decay = c(1, 3, 10, 30))
variants <- c("one curve", "invertebrate predators only", "one curve per kind of predator")
rows <- list()
for (rn in names(r$pairs)) for (w in names(r$pairs[[rn]])) {
  p <- r$pairs[[rn]][[w]]
  cf <- fit_curves(c("tuesday", systems_of(w)))
  pred <- p$animal & !is.na(p$x)
  vert <- p$consumer %in% vert_names
  g_all <- ifelse(pred, g_of(cf$all, p$x), 0)
  g_inv <- ifelse(pred & !vert, g_of(cf$invertebrate, p$x), 0)
  g_type <- ifelse(pred & vert, g_of(cf$vertebrate, p$x), g_inv)
  gs <- list(g_all, ifelse(vert, 0, g_all), g_type)
  lg <- stats::qlogis(pmin(pmax(p$p, 1e-12), 1 - 1e-12))
  none <- sum((p$p - p$y)^2) / nrow(p)
  for (v in seq_along(variants)) for (k in seq_len(nrow(weights))) {
    q <- stats::plogis(lg + ifelse(p$inferred, weights$rule[k], weights$decay[k]) * gs[[v]])
    rows[[length(rows) + 1L]] <- data.frame(
      variant = variants[v], registry = rn, web = w, rule = weights$rule[k],
      decay = weights$decay[k], vertebrate_links = sum(p$y & pred & vert), none = none,
      brier = sum((q - p$y)^2) / nrow(p), stringsAsFactors = FALSE)
  }
}
tab <- do.call(rbind, rows)
tab$diff <- tab$brier - tab$none
curves <- lapply(as.list(fits), function(f) rbind(all = f$all, invertebrate = f$invertebrate,
                                                  vertebrate = f$vertebrate))
print(curves, digits = 3)
print(lapply(as.list(fits), function(f) c(f$pairs, webs_vertebrate = f$webs_vertebrate)))
cat("\nWith the weights of size_model (0.5 and 10):\n")
d0 <- tab[tab$rule == 0.5 & tab$decay == 10, ]
print(stats::reshape(d0[c("variant", "registry", "web", "brier")], idvar = c("registry", "web"),
                     timevar = "variant", direction = "wide"), digits = 4, row.names = FALSE)
cfg <- unique(tab[c("variant", "rule", "decay")])
cfg$id <- seq_len(nrow(cfg))
tab <- merge(tab, cfg)
d <- tab[tab$web != "12", ]
summ <- do.call(rbind, lapply(split(d, list(d$variant, d$rule, d$decay), drop = TRUE), function(x) {
  M <- tapply(x$diff, x$web, mean)
  data.frame(variant = x$variant[1], rule = x$rule[1], decay = x$decay[1], diff = mean(M),
             diff_se = stats::sd(M) / sqrt(length(M)), webs_better = sum(M < 0),
             gray = mean(x$diff[x$registry == "gray"]),
             both_all = mean(x$diff[x$registry == "both_all"]),
             tuesday_gray = tab$brier[tab$web == "12" & tab$registry == "gray" & tab$id == x$id[1]],
             tuesday_both_all = tab$brier[tab$web == "12" & tab$registry == "both_all" &
                                            tab$id == x$id[1]])
}))
summ <- summ[order(summ$diff), ]
rownames(summ) <- NULL
print(summ, digits = 3)
# nested choice (both registries together), by variant and over every variant
M <- tapply(d$diff, list(d$web, d$id), mean)
nested_of <- function(ids) {
  M <- M[, as.character(ids), drop = FALSE]
  nest <- vapply(rownames(M), function(w) {
    M[w, which.min(colMeans(M[rownames(M) != w, , drop = FALSE]))]
  }, numeric(1))
  k <- as.integer(colnames(M)[which.min(colMeans(M))])
  data.frame(cv = min(colMeans(M)), nested = mean(nest),
             nested_se = stats::sd(nest) / sqrt(length(nest)), nested_better = sum(nest < 0),
             chosen = paste0(cfg$variant[k], "; ", cfg$rule[k], ", ", cfg$decay[k]))
}
nested <- do.call(rbind, lapply(c(variants, "any"), function(v) {
  ids <- if (v == "any") cfg$id else cfg$id[cfg$variant == v]
  cbind(data.frame(variant = v), nested_of(ids))
}))
print(nested, digits = 3, row.names = FALSE)
saveRDS(list(webs = tab, summary = summ, nested = nested, curves = curves),
        file.path(out_dir, "size_consumer_type.rds"))
