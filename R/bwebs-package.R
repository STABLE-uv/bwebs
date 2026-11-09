#' *bwebs*: Bayesian Food Webs from Compiled Trophic Interactions
#'
#' *bwebs* (Amador Crespo and Rico 2026) infers the potential trophic links
#' among the taxa of a community from a registry of published interactions
#' with the WebBuilder logic of Gray et al. (2015) ([build_metaweb()]),
#' calibrates how often links inferred with each class of evidence are
#' observed in independent webs ([calibrate_links()]), gives every pair its
#' probability ([link_probabilities()]), samples plausible food webs from
#' those probabilities ([sample_metawebs()]), summarises their topological
#' indicators for every sample of a community ([calc_indicators()],
#' [summarize_indicators()]) and computes their energy fluxes with the
#' approach of *fluxweb* (Gauzens et al. 2019; [calc_fluxes()]).
#'
#' @references
#' Amador Crespo P, Rico A (2026) *bwebs*: Bayesian Food Webs from Compiled
#' Trophic Interactions. R package version 0.4.0.
#' <https://github.com/STABLE-uv/bwebs>
#'
#' Gauzens B, Barnes A, Giling DP, Hines J, Jochum M, Lefcheck JS, Rosenbaum
#' B, Wang S, Brose U (2019) *fluxweb*: an R package to easily estimate energy
#' fluxes in food webs. *Methods in Ecology and Evolution* 10: 270-279.
#' \doi{10.1111/2041-210X.13109}
#'
#' Gray C, Figueroa DH, Hudson LN, Ma A, Perkins D, Woodward G (2015) Joining
#' the dots: an automated method for constructing food webs from compendia of
#' published interactions. *Food Webs* 5: 11-20.
#' \doi{10.1016/j.fooweb.2015.09.001}
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL
