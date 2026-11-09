# Body-mass term with GATEWAy (version 0.3.0): the validation by folds of
# dev/validation/size_folds.R with the masses and the webs of GATEWAy
# version 1.0 (Brose et al. 2019; 290 webs with the body masses of consumers
# and resources; CC BY 4.0).
#
# - Masses of the nodes of the validation webs: by name, or by genus, from
#   every web with masses (GATEWAy, Brose et al. 2005, cheddar). Every
#   validation web whose consumers have masses for at least 60% of them is
#   scored (left out by system, with the class and decay models fitted to the
#   other systems, as in score_rule()).
# - Curve of the size term (within-consumer quadratic effect of
#   x = log10(M_consumer / M_resource), size_curve_fit()), fitted to the webs
#   of the other systems: those of Brose et al. (2005) and cheddar (as in
#   size_model), or those of GATEWAy (freshwater, aquatic or every
#   ecosystem); one curve for every predator, or one for invertebrate and
#   another for vertebrate predators. Never fitted to Tuesday Lake (the
#   independent test, scored with the fits without it) nor to the five chalk
#   streams of Gray et al. (2015) and Thompson et al. (2017), whose links were
#   partly extrapolated with WebBuilder. When a validation web is left out,
#   the webs of GATEWAy that copy it (dev/gateway/gateway_exclusions.csv) and
#   the webs of its system are left out of the curve.
# - Weights on the logit of the inferred links (0 to 1) and of the pairs of
#   the decay (1 to 30), chosen as in dev/validation/size_sources.R: the
#   configuration with the best mean difference in the Brier score over the
#   webs, and the same choice made without each web (nested).
#
# Run from the root of the package after dev/gateway/process_gateway.R:
#   Rscript dev/validation/size_gateway.R <GATEWAy csv> <gateway_registry.rds> \
#     <bodysizes_2008.txt> <cheddar data folder> <globi_registry.rds> [registries]
# registries: comma-separated among gray, both_all, gray_gateway, all
# (default: all of them); one process per registry, then
#   ... merge
pkgload::load_all(".", quiet = TRUE)
args <- commandArgs(trailingOnly = TRUE)
out_dir <- file.path("dev", "validation")
set.seed(2026)

clean <- function(x) {
  x <- tolower(trimws(gsub("[\t Ê]+", " ", as.character(x))))
  x <- gsub("\\s*\\(.*\\)$", "", x)
  x <- gsub("\\s+(sp|spp|agg|sp\\.|spp\\.|agg\\.)\\.?$", "", x)
  x[!is.na(x) & x == ""] <- NA_character_
  x
}
genus_of <- function(x) sub(" .*$", "", x)
animal_cat <- c("invertebrate", "Invertebrate", "ectotherm vertebrate", "endotherm vertebrate",
                "vert.ecto", "vert.endo")
vertebrate_cat <- c("ectotherm vertebrate", "endotherm vertebrate", "vert.ecto", "vert.endo")
vertebrate_classes <- c("Actinopterygii", "Cephalaspidomorphi", "Chondrichthyes", "Amphibia",
                        "Reptilia", "Aves", "Mammalia", "Crocodylia")

# 1. Webs with body masses ---------------------------------------------------------
webs <- list()
# GATEWAy: every food web; nodes are taxa at a life stage
gw_raw <- utils::read.csv(args[1], stringsAsFactors = FALSE, na.strings = c("NA", "-999", ""))
gw_raw$con_key <- paste(trimws(gsub("[\t Ê]+", " ", gw_raw$con.taxonomy)),
                        ifelse(is.na(gw_raw$con.lifestage), "", gw_raw$con.lifestage), sep = "|")
gw_raw$res_key <- paste(trimws(gsub("[\t Ê]+", " ", gw_raw$res.taxonomy)),
                        ifelse(is.na(gw_raw$res.lifestage), "", gw_raw$res.lifestage), sep = "|")
for (fw in unique(gw_raw$foodweb.name)) {
  d <- gw_raw[gw_raw$foodweb.name == fw, ]
  m <- c(stats::setNames(d$con.mass.mean.g., d$con_key), stats::setNames(d$res.mass.mean.g., d$res_key))
  m <- m[!is.na(m) & m > 0]
  type <- c(stats::setNames(d$con.metabolic.type, d$con_key), stats::setNames(d$res.metabolic.type, d$res_key))
  type <- type[!duplicated(names(type))]
  nodes <- unique(c(d$con_key, d$res_key))
  webs[[paste0("gw: ", fw)]] <- list(
    source = "gateway", ecosystem = d$ecosystem.type[1], study = paste("gateway:", fw),
    tuesday = grepl("^Tuesday Lake", fw),
    webbuilder = d$link.citation[1] %in% "Gray et al. (2015), Thompson et al. (2017)",
    systems = character(0),
    nodes = data.frame(name = sub("\\|.*$", "", nodes),
                       log_mass = unname(tapply(log10(m), names(m), mean)[nodes]),
                       animal = unname(type[nodes]) %in% animal_cat,
                       vertebrate = unname(type[nodes]) %in% vertebrate_cat,
                       genus = NA_character_, key = nodes, stringsAsFactors = FALSE),
    links = data.frame(resource = d$res_key, consumer = d$con_key, stringsAsFactors = FALSE)
  )
}
# Brose et al. (2005): one web per location with masses
b <- utils::read.delim(args[3], check.names = FALSE, stringsAsFactors = FALSE,
                       na.strings = c("-999", "-999.0"), quote = "\"", comment.char = "",
                       fileEncoding = "latin1")
names(b) <- trimws(names(b))
brose_webs <- c(caricaie = "Lake Neuch", celtic = "Celtic Sea", benguela = "Benguela",
                broadstone = "Broadstone", frome = "River Frome", tuesday = "Tuesday Lake")
name_of <- function(tax, common) ifelse(is.na(tax) | tax == "", common, tax)
for (sys in names(brose_webs)) {
  d <- b[grepl(brose_webs[[sys]], b[["Geographic location"]], fixed = TRUE), ]
  con <- name_of(d[["Taxonomy consumer"]], d[["Common name(s) consumer"]])
  res <- name_of(d[["Taxonomy resource"]], d[["Common name(s) resource"]])
  m <- c(stats::setNames(d[["Mean mass (g) consumer"]], con), stats::setNames(d[["Mean mass (g) resource"]], res))
  m <- m[!is.na(m) & m > 0 & !is.na(names(m))]
  cats <- c(stats::setNames(d[["Metabolic category consumer"]], con),
            stats::setNames(d[["Metabolic category resource"]], res))
  nodes <- unique(c(con, res))
  nodes <- nodes[!is.na(nodes)]
  webs[[paste0("brose: ", sys)]] <- list(
    source = "brose_cheddar", ecosystem = NA, study = NA, tuesday = sys == "tuesday", webbuilder = FALSE,
    systems = sys,
    nodes = data.frame(name = nodes, log_mass = unname(tapply(log10(m), names(m), mean)[nodes]),
                       animal = nodes %in% names(cats)[cats %in% animal_cat],
                       vertebrate = nodes %in% names(cats)[cats %in% vertebrate_cat],
                       genus = NA_character_, key = nodes, stringsAsFactors = FALSE),
    links = data.frame(resource = res, consumer = con, stringsAsFactors = FALSE)
  )
}
# cheddar: masses in mg (Tuesday Lake in kg)
ld <- function(f) {
  e <- new.env()
  load(file.path(args[4], paste0(f, ".RData")), envir = e)
  get(ls(e)[1], e)
}
from_cheddar <- function(x, systems) {
  n <- x$nodes
  to_g <- if (identical(x$properties$M.units, "kg")) 1000 else 1e-3
  list(source = "brose_cheddar", ecosystem = NA, study = NA, tuesday = "tuesday" %in% systems,
       webbuilder = FALSE, systems = systems,
       nodes = data.frame(name = n$node, log_mass = log10(n$M * to_g), animal = n$category %in% animal_cat,
                          vertebrate = n$category %in% vertebrate_cat,
                          genus = if ("genus" %in% names(n)) n$genus else NA_character_,
                          key = n$node, stringsAsFactors = FALSE),
       links = x$trophic.links[c("resource", "consumer")])
}
webs[["cheddar: Broadstone"]] <- from_cheddar(ld("BroadstoneStream"), "broadstone")
ms <- ld("Millstream")
for (k in names(ms)) webs[[paste("cheddar: Mill Stream", k)]] <- from_cheddar(ms[[k]], "frome")
ph <- ld("pHWebs")
for (k in names(ph)) {
  sys <- c("layer", if (k == "Broadstone") "broadstone", if (k %in% c("Mill Stream", "Bere Stream")) "frome")
  webs[[paste("cheddar:", k)]] <- from_cheddar(ph[[k]], sys)
}
webs[["cheddar: TL84"]] <- from_cheddar(ld("TL84"), "tuesday")
webs[["cheddar: TL86"]] <- from_cheddar(ld("TL86"), "tuesday")
message(length(webs), " webs with masses")

# Predation pairs of each web: every predator (a node with an animal prey)
# with every animal with a mass, except itself
size_pairs <- function(w) {
  n <- w$nodes[!is.na(w$nodes$log_mass), ]
  animals <- n$key[n$animal]
  l <- w$links[w$links$resource %in% animals & w$links$resource != w$links$consumer, ]
  pred <- intersect(unique(l$consumer), n$key)
  if (!length(pred)) return(NULL)
  p <- expand.grid(consumer = pred, resource = animals, stringsAsFactors = FALSE)
  p <- p[p$consumer != p$resource, ]
  p$y <- paste(p$resource, p$consumer) %in% paste(l$resource, l$consumer)
  lm <- stats::setNames(n$log_mass, n$key)
  p$x <- unname(lm[p$consumer]) - unname(lm[p$resource])
  p$vertebrate <- p$consumer %in% n$key[n$vertebrate]
  p
}
train_pairs <- lapply(webs, size_pairs)

# Mass of every name: the mean over the webs where it is a node (log10 g)
mass_rows <- do.call(rbind, lapply(webs, function(w) {
  data.frame(name = clean(w$nodes$name),
             genus = clean(ifelse(is.na(w$nodes$genus) | w$nodes$genus == "", genus_of(clean(w$nodes$name)),
                                  w$nodes$genus)),
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

# Curves fitted to the webs of a source, without the webs of `exclude`
# (systems of Brose et al. 2005 and cheddar, studies of GATEWAy)
gwx <- readRDS(args[2])$exclusions
web_systems <- list(`13` = "frome", `50` = "frome", `46` = "broadstone", `49` = "broadstone",
                    `27` = c("layer", "broadstone", "frome"), `12` = "tuesday")
systems_of <- function(w) if (is.null(web_systems[[w]])) character(0) else web_systems[[w]]
sources <- list(
  brose_cheddar = function(w) w$source == "brose_cheddar",
  gateway_freshwater = function(w) w$source == "gateway" && w$ecosystem %in% c("lakes", "streams"),
  gateway_aquatic = function(w) w$source == "gateway" && w$ecosystem %in% c("lakes", "streams", "marine"),
  gateway_all = function(w) w$source == "gateway"
)
curves <- new.env()
curve_for <- function(src_name, left_out) {
  key <- paste(src_name, left_out, sep = "@")
  if (!is.null(curves[[key]])) return(curves[[key]])
  ex_sys <- c("tuesday", systems_of(left_out))
  ex_gw <- gwx$study_id[gwx$web == left_out]
  use <- vapply(webs, function(w) {
    sources[[src_name]](w) && !w$tuesday && !w$webbuilder && !any(w$systems %in% ex_sys) &&
      !(w$source == "gateway" && w$study %in% ex_gw)
  }, logical(1))
  d <- do.call(rbind, Map(function(p, nm) if (!is.null(p)) cbind(p, web = nm), train_pairs[use], names(webs)[use]))
  d$cons <- paste(d$web, d$consumer)
  fit <- function(dd) if (sum(dd$y) >= 30) size_curve_fit(dd$y, dd$x, dd$cons) else NULL
  out <- list(all = fit(d), invertebrate = fit(d[!d$vertebrate, ]), vertebrate = fit(d[d$vertebrate, ]),
              webs = sum(use))
  # without enough vertebrate (or invertebrate) predators, the curve of every predator
  if (is.null(out$vertebrate)) out$vertebrate <- out$all
  if (is.null(out$invertebrate)) out$invertebrate <- out$all
  assign(key, out, envir = curves)
  out
}
g_of <- function(f, x) {
  b <- f$coefficients
  g <- b[["x2"]] * (x - f$optimum)^2
  g[is.na(g)] <- 0
  g
}

# 2. Validation webs, left out by system ----------------------------------------
g <- readRDS(args[5])
gw <- readRDS(args[2])
gray <- gray2015[gray2015$study_id != "12", ]
wb <- "Gray et al. (2015), Thompson et al. (2017)"
registries <- list(
  gray = list(registry = gray, exclusions = NULL, full = gray2015),
  both_all = list(registry = combine_registries(gray, g$globi_test), exclusions = g$exclusions,
                  full = combine_registries(gray2015, g$globi_test)),
  gray_gateway = list(registry = combine_registries(gray, gw$gateway_test), exclusions = gw$exclusions,
                      full = combine_registries(gray2015, gw$gateway_test)),
  all = list(registry = combine_registries(gray, g$globi_test, gw$gateway_test),
             exclusions = rbind(g$exclusions, gw$exclusions),
             full = combine_registries(gray2015, g$globi_test, gw$gateway_test))
)
rule <- matching_rule()
src <- setdiff(gray2015_sources$study_id[gray2015_sources$validation], "12")
groups <- list(c("13", "50"))
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
  vert <- stats::setNames(nodes$class %in% vertebrate_classes | nodes$phylum %in% "Chordata", nodes$node_id)
  p$vertebrate <- unname(vert[p$consumer])
  lm <- stats::setNames(mass_lookup(nodes$taxon, nodes$genus), nodes$node_id)
  p$x <- unname(lm[p$consumer]) - unname(lm[p$resource])
  attr(p, "mass_coverage") <- mean(!is.na(lm[nodes$node_id[nodes$role == "consumer"]]))
  p
}

weights <- expand.grid(rule = c(0, 0.25, 0.5, 0.75, 1), decay = c(1, 3, 10, 30))
kinds <- c("one curve", "by kind of predator")
score_web <- function(p, w) {
  pred <- p$animal & !is.na(p$x)
  lg <- stats::qlogis(pmin(pmax(p$p, 1e-12), 1 - 1e-12))
  none <- sum((p$p - p$y)^2) / nrow(p)
  out <- list()
  for (s in names(sources)) {
    cf <- curve_for(s, w)
    g_one <- ifelse(pred, g_of(cf$all, p$x), 0)
    g_kind <- ifelse(pred, ifelse(p$vertebrate, g_of(cf$vertebrate, p$x), g_of(cf$invertebrate, p$x)), 0)
    for (kd in kinds) {
      gg <- if (kd == "one curve") g_one else g_kind
      for (k in seq_len(nrow(weights))) {
        q <- stats::plogis(lg + ifelse(p$inferred, weights$rule[k], weights$decay[k]) * gg)
        out[[length(out) + 1L]] <- data.frame(curve = s, kind = kd, rule = weights$rule[k],
                                              decay = weights$decay[k], none = none,
                                              brier = sum((q - p$y)^2) / nrow(p), stringsAsFactors = FALSE)
      }
    }
  }
  do.call(rbind, out)
}

run <- if (length(args) > 5L) strsplit(args[6], ",", fixed = TRUE)[[1]] else names(registries)
if (identical(run, "merge")) {
  parts <- lapply(names(registries), function(rn) {
    f <- file.path(out_dir, sprintf("size_gateway_%s.rds", rn))
    if (file.exists(f)) readRDS(f) else NULL
  })
  parts <- parts[!vapply(parts, is.null, logical(1))]
  tab <- do.call(rbind, lapply(parts, `[[`, "webs"))
  tab$diff <- tab$brier - tab$none
  cfg <- unique(tab[c("curve", "kind", "rule", "decay")])
  cfg$id <- seq_len(nrow(cfg))
  tab <- merge(tab, cfg)
  d <- tab[tab$web != "12", ]
  summarise <- function(dd) {
    M <- tapply(dd$diff, list(dd$web, dd$id), mean)
    data.frame(id = as.integer(colnames(M)), diff = colMeans(M),
               diff_se = apply(M, 2L, stats::sd) / sqrt(nrow(M)), webs_better = colSums(M < 0), webs = nrow(M))
  }
  nested_of <- function(dd, ids) {
    M <- tapply(dd$diff, list(dd$web, dd$id), mean)[, as.character(ids), drop = FALSE]
    nest <- vapply(rownames(M), function(w) {
      M[w, which.min(colMeans(M[rownames(M) != w, , drop = FALSE]))]
    }, numeric(1))
    k <- as.integer(colnames(M)[which.min(colMeans(M))])
    cbind(cfg[cfg$id == k, c("curve", "kind", "rule", "decay")],
          data.frame(cv = min(colMeans(M)), nested = mean(nest),
                     nested_se = stats::sd(nest) / sqrt(length(nest)), nested_better = sum(nest < 0)))
  }
  out <- list()
  for (rn in c(unique(d$registry), "registries together")) {
    dd <- if (rn == "registries together") d else d[d$registry == rn, ]
    s <- merge(cfg, summarise(dd))
    s <- s[order(s$diff), ]
    tl <- tab[tab$web == "12" & (rn == "registries together" | tab$registry == rn), ]
    s$tuesday <- tapply(tl$brier, tl$id, mean)[as.character(s$id)]
    s$tuesday_none <- mean(tl$none[!duplicated(tl$registry)])
    cat("\n", rn, "\n")
    print(utils::head(s, 15), digits = 3, row.names = FALSE)
    nest <- do.call(rbind, lapply(c(names(sources), "any"), function(cs) {
      do.call(rbind, lapply(c(kinds, "any"), function(kd) {
        ids <- cfg$id[(cs == "any" | cfg$curve == cs) & (kd == "any" | cfg$kind == kd)]
        cbind(data.frame(curves = cs, kinds = kd), nested_of(dd, ids))
      }))
    }))
    print(nest, digits = 3, row.names = FALSE)
    out[[rn]] <- list(summary = s, nested = nest)
  }
  saveRDS(list(webs = tab, configurations = cfg, results = out,
               curves = do.call(c, lapply(parts, `[[`, "curves")), coverage = do.call(rbind, lapply(parts, `[[`, "coverage"))),
          file.path(out_dir, "size_gateway.rds"))
  quit(save = "no")
}

for (rn in run) {
  st <- registries[[rn]]
  reg <- st$registry
  v <- validate_rule(reg, src, groups = groups, rule = rule, exclusions = st$exclusions, quiet = TRUE)
  left_out <- stats::setNames(strsplit(v$webs$left_out, ", ", fixed = TRUE), v$webs$study_id)
  fold_p <- function(out_ids) {
    tr <- v$counts[!v$counts$study_id %in% out_ids, ]
    td <- v$decay[!v$decay$study_id %in% out_ids, ]
    fc <- fit_calibration_bayes(tr, 1.5, 2.5, 2L, 500L, 500L, 20L, seed = NULL)
    fd <- fit_decay(td, n_chains = 2L, n_iter = 500L, n_warmup = 500L, thin = 2L)
    list(class = stats::setNames(fc$coefficients$p_new, fc$coefficients$evidence_class), decay = fd)
  }
  rows <- list()
  cover <- list()
  for (w in c(src, "12")) {
    out_ids <- if (w == "12") "12" else left_out[[w]]
    extra <- if (is.null(st$exclusions)) character(0) else st$exclusions$study_id[st$exclusions$web == w]
    reg_w <- if (w != "12") reg else st$full
    train <- reg_w[!reg_w$study_id %in% c(out_ids, extra), ]
    p <- pairs_of(reg_w, w, train)
    cover[[w]] <- data.frame(registry = rn, web = w, mass_coverage = attr(p, "mass_coverage"))
    if (attr(p, "mass_coverage") < 0.6) {
      message(rn, " web ", w, ": masses for ", round(100 * attr(p, "mass_coverage")), "% of the consumers, not scored")
      next
    }
    fp <- fold_p(if (w == "12") character(0) else out_ids)
    cls <- p$evidence_class
    cls[!is.na(cls) & cls == "beyond"] <- "class"
    p$p <- ifelse(p$inferred, unname(fp$class[cls]), 0)
    p$p[is.na(p$p)] <- 0
    p$p[!p$inferred] <- decay_p_new_of(fp$decay, p[!p$inferred, ], n_nodes = length(unique(p$resource)))
    sc <- score_web(p, w)
    rows[[w]] <- cbind(data.frame(registry = rn, web = w, pairs = nrow(p), recorded = sum(p$y),
                                  predation_links = sum(p$y & p$animal & !is.na(p$x)),
                                  vertebrate_links = sum(p$y & p$animal & !is.na(p$x) & p$vertebrate),
                                  stringsAsFactors = FALSE), sc)
    message(rn, " web ", w, " done")
  }
  cf <- lapply(as.list(curves), function(x) {
    lapply(x[c("all", "invertebrate", "vertebrate")], `[`, c("coefficients", "optimum", "pairs", "links", "consumers"))
  })
  saveRDS(list(webs = do.call(rbind, rows), curves = cf, coverage = do.call(rbind, cover)),
          file.path(out_dir, sprintf("size_gateway_%s.rds", rn)))
}
