# size_model: the effect of the predator/prey body-mass ratio on links.
#
# Within-consumer effect of x = log10(M_consumer / M_resource) on the links of
# each predator with the animals of its web: logistic regression with one
# intercept per consumer and a quadratic term in x (size_curve_fit()), fitted
# to the aquatic webs (lakes, streams and marine) of GATEWAy version 1.0
# (Brose et al. 2019; CC BY 4.0) except Tuesday Lake (kept apart as the
# independent test of the package) and the five chalk streams of Gray et al.
# (2015) and Thompson et al. (2017), whose links were partly extrapolated
# with WebBuilder. Nodes are taxa at a life stage, with the mean of their
# masses in the web. Only the coefficients are bundled.
#
# The webs of the curve and the weights of the effect on the links inferred
# by the rule (0.25) and on the pairs of the taxonomic decay (1) are those
# chosen by the validation by folds of dev/validation/size_gateway.R: the 12
# validation webs of the registry with body masses for at least 60% of their
# consumers, left out by system, against the curves of Brose et al. (2005)
# and cheddar, of the freshwater webs of GATEWAy and of all its webs, one
# curve or one per kind of predator, and weights from 0 to 1 and from 1 to
# 30 (see dev/validation/results.md, "Size term").
#
# Source: Brose U, et al. (2019) Predator traits determine food-web
# architecture across ecosystems. Nature Ecology & Evolution 3: 919-927.
# GATEWAy version 1.0, iDiv Data Repository, dataset 283
# (283_2_FoodWebDataBase_2018_12_10.csv).
#
# Use, from the root of the package:
#   Rscript data-raw/size_model.R <283_2_FoodWebDataBase_2018_12_10.csv>
# Packages used only here: pkgload, usethis.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Use: Rscript data-raw/size_model.R <GATEWAy csv>")
pkgload::load_all(".", quiet = TRUE)
animal <- c("invertebrate", "ectotherm vertebrate", "endotherm vertebrate")

g <- utils::read.csv(args[1], stringsAsFactors = FALSE, na.strings = c("NA", "-999", ""))
# the ecosystem of a web is that of its first record (six webs mix records of
# two ecosystems)
eco <- tapply(g$ecosystem.type, g$foodweb.name, `[`, 1L)
g$ecosystem.type <- unname(eco[g$foodweb.name])
g <- g[g$ecosystem.type %in% c("lakes", "streams", "marine"), ]
g <- g[!grepl("^Tuesday Lake", g$foodweb.name), ]
g <- g[!g$link.citation %in% "Gray et al. (2015), Thompson et al. (2017)", ]
node <- function(name, stage) {
  paste(trimws(gsub("[\t Ê]+", " ", name)), ifelse(is.na(stage), "", stage), sep = "|")
}
g$con <- node(g$con.taxonomy, g$con.lifestage)
g$res <- node(g$res.taxonomy, g$res.lifestage)

# Predation pairs of each web: every predator (a node with an animal prey)
# with every animal with a mass, except itself
pairs <- do.call(rbind, lapply(split(g, g$foodweb.name), function(d) {
  m <- c(stats::setNames(d$con.mass.mean.g., d$con), stats::setNames(d$res.mass.mean.g., d$res))
  m <- m[!is.na(m) & m > 0]
  lm <- tapply(log10(m), names(m), mean)
  type <- c(stats::setNames(d$con.metabolic.type, d$con), stats::setNames(d$res.metabolic.type, d$res))
  type <- type[!duplicated(names(type))]
  animals <- intersect(names(type)[type %in% animal], names(lm))
  l <- d[d$res %in% animals & d$res != d$con, c("res", "con")]
  pred <- intersect(unique(l$con), names(lm))
  if (!length(pred)) return(NULL)
  p <- expand.grid(consumer = pred, resource = animals, stringsAsFactors = FALSE)
  p <- p[p$consumer != p$resource, ]
  p$y <- paste(p$resource, p$consumer) %in% paste(l$res, l$con)
  p$x <- unname(lm[p$consumer]) - unname(lm[p$resource])
  p$web <- d$foodweb.name[1]
  p$ecosystem <- d$ecosystem.type[1]
  p
}))
fit <- size_curve_fit(pairs$y, pairs$x, paste(pairs$web, pairs$consumer))
stopifnot(fit$converged)
used <- pairs[paste(pairs$web, pairs$consumer) %in%
                names(which(tapply(pairs$y, paste(pairs$web, pairs$consumer), function(v) any(v) && !all(v)))), ]

size_model <- structure(
  list(
    coefficients = fit$coefficients,
    std_error = fit$std_error,
    optimum = fit$optimum,
    weights = c(rule = 0.25, decay = 1),
    data = data.frame(
      ecosystem = c("lakes", "streams", "marine"),
      webs = vapply(c("lakes", "streams", "marine"), function(e) length(unique(used$web[used$ecosystem == e])),
                    integer(1)),
      pairs = vapply(c("lakes", "streams", "marine"), function(e) sum(used$ecosystem == e), integer(1)),
      links = vapply(c("lakes", "streams", "marine"), function(e) sum(used$y[used$ecosystem == e]), integer(1)),
      row.names = NULL, stringsAsFactors = FALSE
    ),
    sources = c(
      "Brose U, et al. (2019) Predator traits determine food-web architecture across ecosystems. Nature Ecology & Evolution 3: 919-927. https://doi.org/10.1038/s41559-019-0899-x",
      "Brose U, et al. GlobAL daTabasE of traits and food Web Architecture (GATEWAy) version 1.0. iDiv Data Repository, dataset 283 (CC BY 4.0)."
    )
  ),
  class = "bwebs_size_model"
)
print(size_model)
print(size_model$data)
usethis::use_data(size_model, overwrite = TRUE)
