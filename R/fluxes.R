# Energy fluxes -------------------------------------------------------------------
#
# The food-web energetics approach of fluxweb (Gauzens et al. 2019): at steady
# state, the energy a consumer assimilates from its resources balances its
# metabolic losses and the energy its consumers take from it. calc_fluxes()
# solves it in every realization of sample_metawebs() (its links) with the
# populations of the same draw (biomasses and body masses, for example the
# draws of bmassa), with metabolic losses from body mass and temperature
# (Brown et al. 2004) and assimilation efficiencies by type of resource.

# Parameters -----------------------------------------------------------------------

flux_metabolic_types <- c("unicellular", "invertebrate", "fish", "ectotherm_vertebrate",
                          "endotherm")
flux_resource_types <- c("animal", "producer", "detritus")
boltzmann_ev <- 8.617333262e-5   # Boltzmann constant, eV K-1

flux_parameter_set <- function(set) {
  metabolism <- data.frame(
    metabolic_type = flux_metabolic_types,
    intercept = c(15.85, 17.17, 18.47, 18.18, 19.50),
    exponent = rep(-0.29, 5L),
    source = c(
      "Polazzo et al. (2023), from Brown et al. (2004): unicellular organisms",
      "Polazzo et al. (2023) and Gauzens et al. (2019, Table 3), from Brown et al. (2004)",
      "Polazzo et al. (2023), from Brown et al. (2004): fish",
      "Gauzens et al. (2019, Table 3), from Brown et al. (2004): ectotherm vertebrates",
      "Gauzens et al. (2019, Table 3), from Brown et al. (2004): endotherm vertebrates"
    ),
    stringsAsFactors = FALSE
  )
  efficiencies <- data.frame(
    resource_type = flux_resource_types,
    efficiency = c(0.906, 0.77, 0.158),
    source = c(
      "Lang et al. (2017), as in Gauzens et al. (2019) and Polazzo et al. (2023)",
      "Polazzo et al. (2023), from Landry et al. (1984)",
      "Lang et al. (2017), as in Gauzens et al. (2019)"
    ),
    stringsAsFactors = FALSE
  )
  if (set == "polazzo_2023") {
    e <- 0.69
    src <- "Polazzo et al. (2023): the slope across groups of Brown et al. (2004, Fig. 1A)"
  } else {
    e <- 0.63
    src <- paste("Brown et al. (2004): the value with which the rates were corrected for",
                 "temperature before the intercepts were fitted (Fig. 1B)")
  }
  list(metabolism = metabolism, activation_energy = e, activation_source = src,
       efficiencies = efficiencies)
}

#' Parameters of the energy fluxes
#'
#' `flux_parameters()` gives the parameters with which [calc_fluxes()]
#' computes the metabolic losses of the consumers, from their body mass and
#' the temperature, and the efficiencies with which they assimilate each type
#' of resource: a published set, with any of its values replaced.
#'
#' @param set Name of the set: `"brown_2004"` (default), the metabolic
#'   parameters of Brown et al. (2004), intercepts and exponent with the
#'   activation energy with which they were fitted (0.63 eV), and the
#'   efficiencies of Polazzo et al. (2023); or `"polazzo_2023"`, the
#'   parameters of Polazzo et al. (2023) for the food webs of freshwater
#'   mesocosms, the same with 0.69 eV (see Details). Both add the intercepts
#'   of Gauzens et al. (2019) for the other vertebrates and an efficiency for
#'   detritus.
#' @param metabolism `NULL` (default) or a data frame with columns
#'   `metabolic_type`, `intercept` and `exponent` (and optionally `source`)
#'   whose rows replace those of the same metabolic types of the set or add
#'   new types.
#' @param activation_energy `NULL` (default) or the activation energy
#'   \eqn{E} of metabolism, in eV.
#' @param efficiencies `NULL` (default) or a vector of assimilation
#'   efficiencies (in (0, 1]) named by type of resource (`"animal"`,
#'   `"producer"`, `"detritus"`) that replace those of the set.
#'
#' @return An object of class `bwebs_flux_parameters`: a list with `set`,
#'   `metabolism` (one row per metabolic type: `metabolic_type`,
#'   `intercept`, `exponent`, `source`), `activation_energy`,
#'   `activation_source`, `efficiencies` (`resource_type`, `efficiency`,
#'   `source`) and `boltzmann`, the Boltzmann constant in eV
#'   K\eqn{^{-1}}{^-1}.
#'
#' @details
#' The metabolic rate per gram of wet mass of a population of mean individual
#' wet mass \eqn{M} (g) at temperature \eqn{T} (K) is
#' \deqn{x = \exp(c + b \ln M - E / (k T))}{x = exp(c + b log(M) - E / (k T))}
#' in W g\eqn{^{-1}}{^-1}: \eqn{c} is the intercept of its metabolic type on
#' the natural-log scale (\eqn{\ln i_0}{log(i0)} of Brown et al. 2004, for
#' the whole-organism rate in W and the mass in g), \eqn{b} the exponent of
#' the rate per unit of mass (the whole-organism exponent 0.71 minus 1) and
#' \eqn{k} the Boltzmann constant. The intercepts are logarithms, the
#' normalization constants \eqn{C = \ln i_0}{C = log(i0)} of Brown et al.
#' (2004): Polazzo et al. (2023) write the equation with them as multipliers
#' of \eqn{M^b}{M^b}, as the code of the vignette of *fluxweb* does, which
#' would give losses about a million times lower and other ratios between
#' groups. The losses of a population are \eqn{x} times its biomass.
#'
#' Metabolic types and intercepts: `"unicellular"` (heterotrophic protists,
#' 15.85), `"invertebrate"` (17.17), `"fish"` (18.47), used by Polazzo et
#' al. (2023), and `"ectotherm_vertebrate"` (amphibians and reptiles, 18.18)
#' and `"endotherm"` (birds and mammals, 19.50), from Table 3 of Gauzens et
#' al. (2019); all from Brown et al. (2004), with the exponent -0.29.
#' Polazzo et al. (2023) gave 15.85 to the unicellular producers, whose
#' losses do not change the fluxes (they only add to the outflow of the
#' producers), so here it applies to heterotrophic protists. Brown et al.
#' (2004) give an intercept for each group of vertebrates: use `metabolism`
#' to give amphibians or reptiles their own.
#'
#' Assimilation efficiencies by type of resource: animals 0.906 (Lang et al.
#' 2017, as in Gauzens et al. 2019 and Polazzo et al. 2023), primary
#' producers 0.77 (Polazzo et al. 2023, from Landry et al. 1984; Lang et al.
#' 2017 give 0.545 for plants) and detritus 0.158 (Lang et al. 2017;
#' Polazzo et al. 2023 give none).
#'
#' Activation energy: 0.63 eV by default. Brown et al. (2004) estimated
#' 0.69 eV as the common slope of all groups (Fig. 1A), the value of Polazzo
#' et al. (2023), but corrected the rates for temperature with 0.63 eV, the
#' mean without the endotherms, before fitting the exponent (0.71, so -0.29
#' per gram) and the intercepts of each group (Fig. 1B), so the intercepts
#' go with 0.63 eV. With 0.69 eV (`set = "polazzo_2023"`) the rates are 10
#' to 13 times lower between 0 and 30 °C (for an invertebrate of 1 mg at
#' 20 °C, 0.3 instead of 3 mW per g). The fluxes are proportional to the
#' losses, so this scales all of them, and changes their response to
#' temperature slightly, but not their relative values within a web at one
#' temperature.
#'
#' @references
#' Brown JH, Gillooly JF, Allen AP, Savage VM, West GB (2004) Toward a
#' metabolic theory of ecology. *Ecology* 85: 1771-1789.
#' \doi{10.1890/03-9000}
#'
#' Gauzens B, Barnes A, Giling DP, Hines J, Jochum M, Lefcheck JS, Rosenbaum
#' B, Wang S, Brose U (2019) *fluxweb*: an R package to easily estimate energy
#' fluxes in food webs. *Methods in Ecology and Evolution* 10: 270-279.
#' \doi{10.1111/2041-210X.13109}
#'
#' Landry MR, Hassett RP, Fagerness V, Downs J, Lorenzen CJ (1984) Effect of
#' food acclimation on assimilation efficiency of *Calanus pacificus*.
#' *Marine Ecology Progress Series* 16. \doi{10.3354/meps016127}
#'
#' Lang B, Ehnes RB, Brose U, Rall BC (2017) Temperature and consumer type
#' dependencies of energy flows in natural communities. *Oikos* 126:
#' 1717-1725. \doi{10.1111/oik.04419}
#'
#' Polazzo F, Hermann M, Crettaz-Minaglia M, Rico A (2023) Impacts of extreme
#' climatic events on trophic network complexity and multidimensional
#' stability. *Ecology* 104: e3951. \doi{10.1002/ecy.3951}
#'
#' @seealso [calc_fluxes()].
#' @export
#' @examples
#' flux_parameters()
#' # losses of the invertebrates doubled, for a sensitivity analysis, and
#' # another efficiency for detritus
#' flux_parameters(metabolism = data.frame(metabolic_type = "invertebrate",
#'                                         intercept = 17.17 + log(2), exponent = -0.29),
#'                 efficiencies = c(detritus = 0.2))
flux_parameters <- function(set = "brown_2004", metabolism = NULL, activation_energy = NULL,
                            efficiencies = NULL) {
  check_choice(set, c("brown_2004", "polazzo_2023"))
  p <- flux_parameter_set(set)
  if (!is.null(metabolism)) {
    check_data_frame(metabolism, allow_empty = FALSE)
    check_columns(metabolism, c("metabolic_type", "intercept", "exponent"), arg = "metabolism")
    check_numeric_column(metabolism, "intercept", arg = "metabolism")
    check_numeric_column(metabolism, "exponent", arg = "metabolism")
    type <- as.character(metabolism$metabolic_type)
    if (anyNA(type) || any(!nzchar(type)) || anyDuplicated(type)) {
      abort_input("`metabolism` must have one row per metabolic type, without missing types.")
    }
    m <- data.frame(metabolic_type = type, intercept = as.numeric(metabolism$intercept),
                    exponent = as.numeric(metabolism$exponent),
                    source = if ("source" %in% names(metabolism)) {
                      as.character(metabolism$source)
                    } else {
                      rep("user", length(type))
                    },
                    stringsAsFactors = FALSE)
    i <- match(type, p$metabolism$metabolic_type)
    p$metabolism[i[!is.na(i)], ] <- m[!is.na(i), , drop = FALSE]
    p$metabolism <- rbind(p$metabolism, m[is.na(i), , drop = FALSE])
    rownames(p$metabolism) <- NULL
  }
  if (!is.null(activation_energy)) {
    check_number(activation_energy, lower = 0)
    p$activation_energy <- activation_energy
    p$activation_source <- "user"
  }
  if (!is.null(efficiencies)) {
    if (!is.numeric(efficiencies) || is.null(names(efficiencies)) ||
        !all(names(efficiencies) %in% flux_resource_types) ||
        anyDuplicated(names(efficiencies)) || anyNA(efficiencies) ||
        any(efficiencies <= 0 | efficiencies > 1)) {
      abort_input(sprintf(
        "`efficiencies` must be numbers in (0, 1] named by type of resource: %s.",
        format_items(quote_items(flux_resource_types), max = Inf)
      ))
    }
    i <- match(names(efficiencies), p$efficiencies$resource_type)
    p$efficiencies$efficiency[i] <- unname(efficiencies)
    p$efficiencies$source[i] <- "user"
  }
  structure(c(list(set = set), p, list(boltzmann = boltzmann_ev)),
            class = "bwebs_flux_parameters")
}

#' @export
print.bwebs_flux_parameters <- function(x, ...) {
  cat(sprintf("<bwebs_flux_parameters> %s\n", x$set))
  cat("  metabolic rate per g of wet mass: exp(c + b log(M) - E / (k T)) W/g; M in g, T in K\n")
  cat(sprintf("  activation energy E = %.2f eV (%s)\n", x$activation_energy,
              x$activation_source))
  m <- x$metabolism
  for (i in seq_len(nrow(m))) {
    cat(sprintf("  %-21s c = %.2f, b = %.2f\n", m$metabolic_type[i], m$intercept[i],
                m$exponent[i]))
  }
  e <- x$efficiencies
  cat(sprintf("  assimilation efficiencies: %s\n",
              paste(sprintf("%s %.3g", e$resource_type, e$efficiency), collapse = ", ")))
  invisible(x)
}

# Metabolic rate per gram of wet mass (W g-1) of populations of mean wet mass
# `mass` (g) at `temperature` (degrees C), of metabolic types `type`.
metabolic_rate <- function(mass, temperature, type, parameters) {
  m <- parameters$metabolism
  i <- match(type, m$metabolic_type)
  exp(m$intercept[i] + m$exponent[i] * log(mass) -
        parameters$activation_energy / (parameters$boltzmann * (temperature + 273.15)))
}

# Types of the nodes ---------------------------------------------------------------

protist_phyla <- c("Amoebozoa", "Apicomplexa", "Bigyra", "Cercozoa", "Choanozoa", "Ciliophora",
                   "Cryptophyta", "Dinoflagellata", "Dinophyta", "Discosea", "Euglenophyta",
                   "Euglenozoa", "Evosea", "Foraminifera", "Haptophyta", "Heliozoa", "Loukozoa",
                   "Metamonada", "Myzozoa", "Ochrophyta", "Percolozoa", "Protozoa", "Retaria",
                   "Sarcomastigophora", "Tubulinea")
protist_classes <- c("Armophorea", "Choanoflagellatea", "Colpodea", "Dinophyceae",
                     "Heterotrichea", "Kinetoplastea", "Litostomatea", "Nassophorea",
                     "Oligohymenophorea", "Phyllopharyngea", "Prostomatea", "Spirotrichea",
                     "Tubulinea")
fish_classes <- c("Actinopteri", "Actinopterygii", "Cephalaspidomorphi", "Chondrichthyes",
                  "Coelacanthi", "Dipnoi", "Elasmobranchii", "Holocephali", "Myxini",
                  "Petromyzonti", "Pisces", "Sarcopterygii", "Teleostei")
ectotherm_vertebrate_classes <- c("Amphibia", "Crocodylia", "Reptilia", "Squamata", "Testudines")
endotherm_classes <- c("Aves", "Mammalia")

# Default type of resource of every node (as prey: animal for consumers,
# producer or detritus by role) and metabolic type of the consumers (NA for
# the other nodes), from the class and phylum.
default_node_types <- function(nodes) {
  role <- nodes$role
  column <- function(n) {
    if (n %in% names(nodes)) as.character(nodes[[n]]) else rep(NA_character_, nrow(nodes))
  }
  cls <- column("class")
  phy <- column("phylum")
  resource <- ifelse(role == "producer", "producer",
                     ifelse(role == "detritus", "detritus", "animal"))
  met <- rep("invertebrate", nrow(nodes))
  met[phy %in% protist_phyla | cls %in% protist_classes] <- "unicellular"
  met[cls %in% fish_classes] <- "fish"
  met[cls %in% ectotherm_vertebrate_classes] <- "ectotherm_vertebrate"
  met[cls %in% endotherm_classes] <- "endotherm"
  met[role != "consumer"] <- NA_character_
  list(resource = stats::setNames(resource, nodes$node_id),
       metabolic = stats::setNames(met, nodes$node_id))
}

check_named_types <- function(x, ids, allowed, arg, call = sys.call(-1)) {
  if (!is.character(x) || is.null(names(x)) || anyNA(x) || anyDuplicated(names(x)) ||
      !all(names(x) %in% ids) || !all(x %in% allowed)) {
    abort_input(sprintf(
      "`%s` must be a character vector named by node_id of the metaweb, with values in %s.",
      arg, format_items(quote_items(allowed), max = Inf)
    ), call = call)
  }
  invisible(x)
}

# Populations ----------------------------------------------------------------------

mass_unit_g <- c(pg = 1e-12, ng = 1e-9, ug = 1e-6, mg = 1e-3, g = 1)

# The populations of calc_fluxes() and node_masses(): one source or a list of
# sources (data frames and simulations of bmassa). Each source gives its
# samples, its number of draws (NA when it is the same in every draw) and
# `get(d)`: the node (index in `ids`), sample, biomass and mean individual
# mass (g of wet mass) of its rows in draw d.
population_sources <- function(populations, ids, sample_col, node_col, draw_col, biomass_col,
                               mass_col, call = sys.call(-1)) {
  single <- is.data.frame(populations) || inherits(populations, "bmassa_simulation")
  sources <- if (single) list(populations) else populations
  if (!is.list(sources) || length(sources) == 0L) {
    abort_input(paste("`populations` must be a data frame, a simulation of bmassa or a list",
                      "of them."), call = call)
  }
  lapply(seq_along(sources), function(i) {
    arg <- if (single) "populations" else sprintf("populations[[%d]]", i)
    p <- sources[[i]]
    if (inherits(p, "bmassa_simulation")) {
      bmassa_source(p, ids, sample_col, node_col, arg, call)
    } else if (is.data.frame(p)) {
      table_source(p, ids, sample_col, node_col, draw_col, biomass_col, mass_col, arg, call)
    } else {
      abort_input(sprintf("`%s` must be a data frame or a simulation of bmassa, not %s.", arg,
                          describe_value(p)), call = call)
    }
  })
}

table_source <- function(p, ids, sample_col, node_col, draw_col, biomass_col, mass_col, arg,
                         call) {
  check_data_frame(p, arg = arg, allow_empty = FALSE, call = call)
  check_columns(p, c(sample_col = sample_col, node_col = node_col, biomass_col = biomass_col),
                arg = arg, call = call)
  check_numeric_column(p, biomass_col, lower = 0, allow_na = TRUE, arg = arg, call = call)
  has_mass <- mass_col %in% names(p)
  if (has_mass) {
    check_numeric_column(p, mass_col, lower = 0, include_lower = FALSE, allow_na = TRUE,
                         arg = arg, call = call)
  }
  samp <- as.character(p[[sample_col]])
  node <- as.character(p[[node_col]])
  if (anyNA(samp) || anyNA(node)) {
    abort_input(sprintf("Columns '%s' and '%s' of `%s` must not have missing values.",
                        sample_col, node_col, arg), call = call)
  }
  unknown <- setdiff(node, ids)
  if (length(unknown)) {
    abort_input(sprintf("`%s` has nodes that are not in the metaweb: %s.", arg,
                        format_items(unknown)), call = call)
  }
  node <- match(node, ids)
  biomass <- as.numeric(p[[biomass_col]])
  mass <- if (has_mass) as.numeric(p[[mass_col]]) else rep(NA_real_, nrow(p))
  rows <- NULL
  n_draws <- NA_integer_
  if (draw_col %in% names(p)) {
    check_numeric_column(p, draw_col, lower = 1, whole_numbers = TRUE, arg = arg, call = call)
    draw <- as.integer(p[[draw_col]])
    n_draws <- max(draw)
    rows <- split(seq_along(draw), factor(draw, levels = seq_len(n_draws)))
    empty <- which(lengths(rows) == 0L)
    if (length(empty)) {
      abort_input(sprintf("`%s` has no rows for draws %s; give every draw from 1 to %d.", arg,
                          format_items(empty), n_draws), call = call)
    }
  }
  list(arg = arg, samples = unique(samp), n_draws = n_draws, has_mass = has_mass,
       get = function(d) {
         r <- if (is.null(rows)) seq_along(node) else rows[[d]]
         list(node = node[r], sample = samp[r], biomass = biomass[r], mass = mass[r])
       })
}

bmassa_source <- function(p, ids, sample_col, node_col, arg, call) {
  if (!is.list(p) || !is.data.frame(p$rows) || !is.matrix(p$ind_mass) ||
      nrow(p$ind_mass) != nrow(p$rows)) {
    abort_input(sprintf("`%s` is not a valid simulation of bmassa (rows and ind_mass).", arg),
                call = call)
  }
  r <- p$rows
  rarg <- paste0(arg, "$rows")
  check_columns(r, c(sample_col = sample_col, node_col = node_col), arg = rarg, call = call)
  need <- c("abundance", "currency", "unit")
  if (!all(need %in% names(r))) {
    abort_input(sprintf("`%s` lacks the columns %s of a simulation of bmassa.", rarg,
                        format_items(quote_items(setdiff(need, names(r))), max = Inf)),
                call = call)
  }
  samp <- as.character(r[[sample_col]])
  node <- as.character(r[[node_col]])
  if (anyNA(samp) || anyNA(node)) {
    abort_input(sprintf(paste(
      "Columns '%s' and '%s' of `%s` must not have missing values: give every row of the",
      "abundance table its sample and its node, and leave out the rows that are not nodes",
      "of the food web."), sample_col, node_col, rarg), call = call)
  }
  unknown <- setdiff(node, ids)
  if (length(unknown)) {
    abort_input(sprintf("`%s` has nodes that are not in the metaweb: %s.", rarg,
                        format_items(unknown)), call = call)
  }
  drawn <- rowSums(!is.na(p$ind_mass)) > 0
  currency <- as.character(r$currency)
  if (any(drawn & (is.na(currency) | currency != "wet_mass"))) {
    abort_input(sprintf(paste(
      "The masses of `%s` must be wet mass, the currency of the metabolic parameters:",
      "convert them with convert_currency(to = \"wet_mass\") of bmassa before simulate()."),
      arg), call = call)
  }
  f <- unname(mass_unit_g[as.character(r$unit)])
  if (any(drawn & is.na(f))) {
    abort_input(sprintf("Units of mass not recognised in `%s`: %s.", rarg,
                        format_items(unique(r$unit[drawn & is.na(f)]))), call = call)
  }
  if ("biomass_unit" %in% names(r)) {
    per <- sub("^[^/]*/", "", as.character(r$biomass_unit))
    ok <- drawn & !is.na(r$biomass_unit)
    n_basis <- tapply(per[ok], samp[ok], function(z) length(unique(z)))
    mixed <- names(n_basis)[n_basis > 1L]
    if (length(mixed)) {
      abort_input(sprintf(paste(
        "Samples %s of `%s` have abundances on different bases (%s), whose biomasses cannot",
        "be added in a web: express all the abundances on one basis (for instance per",
        "mesocosm, with `per` in bmassa)."), format_items(mixed), arg,
        format_items(unique(per[ok & samp %in% mixed]))), call = call)
    }
  }
  node <- match(node, ids)
  abundance <- as.numeric(r$abundance)
  absent <- !is.na(abundance) & abundance == 0
  f[!drawn] <- 1
  list(arg = arg, samples = unique(samp), n_draws = ncol(p$ind_mass), has_mass = TRUE,
       get = function(d) {
         m <- p$ind_mass[, d] * f
         b <- abundance * m
         b[absent] <- 0
         list(node = node, sample = samp, biomass = b, mass = m)
       })
}

# Rows of all the sources in draw d (the sources without draws are the same in
# every draw).
population_draw <- function(sources, d) {
  parts <- lapply(sources, function(s) s$get(if (is.na(s$n_draws)) 1L else d))
  list(node = unlist(lapply(parts, `[[`, "node"), use.names = FALSE),
       sample = unlist(lapply(parts, `[[`, "sample"), use.names = FALSE),
       biomass = unlist(lapply(parts, `[[`, "biomass"), use.names = FALSE),
       mass = unlist(lapply(parts, `[[`, "mass"), use.names = FALSE))
}

# Number of draws of the populations that can be used: the smallest number
# among the sources with draws (NA if none has draws).
population_n_draws <- function(sources) {
  nd <- vapply(sources, function(s) as.integer(s$n_draws), integer(1))
  if (all(is.na(nd))) NA_integer_ else min(nd, na.rm = TRUE)
}

#' Body masses of the nodes from their populations
#'
#' `node_masses()` gives the mean individual mass of every node in every
#' draw of its populations, pooling all the samples: the total biomass
#' divided by the total number of individuals. It is the `masses` argument
#' of [sample_metawebs()] (and of [link_probabilities()] and
#' [score_metaweb()], with one draw), so that realization `d` of the web gets
#' the size term of draw `d` of the populations, the same draw that
#' [calc_fluxes()] then uses for its biomasses.
#'
#' @param x A [build_metaweb()] or [sample_metawebs()] result: the nodes.
#' @inheritParams calc_fluxes
#' @param n_draws `NULL` (default: all the draws of the populations) or the
#'   number of draws to return (the first ones), such as the `n_draws` of
#'   [sample_metawebs()].
#'
#' @return A matrix nodes (rows, named by `node_id`, all the nodes of `x`) x
#'   draws of masses in g of wet mass, or a vector named by `node_id` if the
#'   populations have no draws. Nodes without populations with a positive
#'   biomass and a mass, such as producers and detritus, get `NA`.
#'
#' @details
#' The number of individuals of a population is its biomass divided by its
#' mean mass (its abundance, for a simulation of *bmassa*, Amador Crespo and
#' Rico 2026); populations with missing values are left out.
#'
#' @references
#' Amador Crespo P, Rico A (2026) *bmassa*: Biomass of Aquatic Taxa from
#' Abundances and Body Sizes. R package version 0.5.0.
#' <https://github.com/STABLE-uv/bmassa>
#'
#' @seealso [calc_fluxes()], [sample_metawebs()].
#' @export
#' @examples
#' nodes <- prepare_nodes(example_taxa)
#' mw <- build_metaweb(nodes, combine_registries(gray2015, gateway))
#' m <- node_masses(mw, example_populations)
#' head(m[!is.na(m)])
node_masses <- function(x, populations, sample_col = "sample_id", node_col = "node_id",
                        draw_col = "draw", biomass_col = "biomass", mass_col = "mass",
                        n_draws = NULL) {
  if (!inherits(x, c("bwebs_draws", "bwebs_metaweb"))) {
    abort_input("`x` must be created with build_metaweb() or sample_metawebs().")
  }
  for (a in c("sample_col", "node_col", "draw_col", "biomass_col", "mass_col")) {
    check_string(get(a), arg = a)
  }
  ids <- x$nodes$node_id
  sources <- population_sources(populations, ids, sample_col, node_col, draw_col, biomass_col,
                                mass_col)
  available <- population_n_draws(sources)
  if (!is.null(n_draws)) {
    check_number(n_draws, lower = 1)
    if (is.na(available) || n_draws > available) {
      abort_input(sprintf("`n_draws` is %d but the populations have %s draws.",
                          as.integer(n_draws), if (is.na(available)) "no" else available))
    }
  }
  nd <- if (is.na(available)) 1L else if (is.null(n_draws)) available else as.integer(n_draws)
  out <- matrix(NA_real_, length(ids), nd, dimnames = list(ids, NULL))
  for (d in seq_len(nd)) {
    pd <- population_draw(sources, d)
    ok <- !is.na(pd$biomass) & pd$biomass > 0 & !is.na(pd$mass)
    if (!any(ok)) next
    b <- rowsum(pd$biomass[ok], pd$node[ok])
    n <- rowsum(pd$biomass[ok] / pd$mass[ok], pd$node[ok])
    out[as.integer(rownames(b)), d] <- b[, 1L] / n[, 1L]
  }
  if (is.na(available)) stats::setNames(out[, 1L], ids) else out
}

# Steady-state fluxes ----------------------------------------------------------------

# Fluxes of the food-web energetics approach with efficiencies defined by the
# resource (Gauzens et al. 2019, eq. 6): F = (diag(W'e + b) - W)^-1 X, where W
# holds the preferences of each consumer for its resources (the columns of
# `A`, scaled by the biomass of the resources if `scale`, summing to 1), e the
# efficiencies of the resources, b is 1 for the nodes without resources and X
# the losses; the flux from i to j is W[i, j] F[j]. NULL when the system is
# singular or has no non-negative solution.
flux_matrix_internal <- function(A, biomass, losses, efficiency, scale = TRUE) {
  W <- A
  if (scale) W <- W * biomass
  consumer <- colSums(A) > 0
  W[, consumer] <- W[, consumer, drop = FALSE] /
    rep(colSums(W)[consumer], each = nrow(W))
  M <- -W
  diag(M) <- diag(M) + as.vector(crossprod(W, efficiency)) + as.numeric(!consumer)
  Fv <- tryCatch(solve(M, losses), error = function(e) NULL)
  if (is.null(Fv) || any(!is.finite(Fv)) || any(Fv < -1e-10 * max(abs(Fv)))) {
    return(NULL)
  }
  W * rep(pmax(Fv, 0), each = nrow(W))
}

flux_matrix_fluxweb <- function(A, biomass, losses, efficiency, scale = TRUE) {
  out <- tryCatch(
    fluxweb::fluxing(A, biomasses = as.numeric(biomass), losses = as.numeric(losses),
                     efficiencies = as.numeric(efficiency), bioms.prefs = scale,
                     bioms.losses = FALSE, ef.level = "prey"),
    error = function(e) NULL
  )
  if (is.null(out)) return(NULL)
  out <- unname(as.matrix(out))
  dimnames(out) <- dimnames(A)
  out
}

# Adjacency resource x consumer of realization d of `x` (all the inferred
# links of a metaweb), without cannibalism.
realization_adjacency <- function(x, d) {
  ids <- x$nodes$node_id
  if (!inherits(x, "bwebs_draws")) {
    a <- x$adjacency[ids, ids, drop = FALSE]
    storage.mode(a) <- "double"
    diag(a) <- 0
    return(a)
  }
  a <- matrix(0, length(ids), length(ids), dimnames = list(ids, ids))
  on <- x$z[, d]
  if (any(on)) a[cbind(x$links$resource[on], x$links$consumer[on])] <- 1
  bg <- x$background
  if (!is.null(bg) && length(bg$present[[d]])) {
    b <- bg$present[[d]]
    a[cbind(bg$pairs$resource[b], bg$pairs$consumer[b])] <- 1
  }
  diag(a) <- 0
  a
}

#' Energy fluxes of the sampled food webs
#'
#' `calc_fluxes()` computes the energy fluxes among the nodes of every sample
#' in every realization of the food web, with the food-web energetics
#' approach of the package *fluxweb* (Gauzens et al. 2019): at steady state, the
#' energy a consumer assimilates from its resources balances its metabolic
#' losses and the energy its consumers take from it. Realization `d` has the
#' links of draw `d` of [sample_metawebs()] and the biomasses and body
#' masses of draw `d` of the populations, such as the draws of
#' `simulate()` of *bmassa* (Amador Crespo and Rico 2026), so that the
#' uncertainty of the links and
#' that of the biomasses are both carried into the fluxes.
#'
#' @param x A [sample_metawebs()] result, or a [build_metaweb()] result to
#'   use all its inferred links in every realization.
#' @param populations The populations of every sample (for example a taxon
#'   at a life stage in a habitat): a data frame, a simulation of
#'   *bmassa* (`simulate()` on an estimate of biomass in wet mass), or a
#'   list of them, such as the simulation of the animals and a table of the
#'   biomasses of producers and detritus. A data frame has one row per
#'   population and sample, and per draw if the biomasses vary between
#'   draws, with the sample, the node, the biomass (**g of wet mass** per
#'   unit of the sampling basis, such as per mesocosm or per m\eqn{^2}{^2};
#'   0 is an observed absence), the mean individual mass of consumers (**g
#'   of wet mass**; optional for producers and detritus) and optionally the
#'   draw (1, 2, ...; without it, the same populations in every
#'   realization). In a simulation of *bmassa*, `sample_col` and
#'   `node_col` name columns of its abundance table, and the biomass of a
#'   row is its abundance times its mean mass in each draw.
#' @param temperature Temperature in °C: one number for all the samples, or
#'   a data frame with the sample and the temperature (Polazzo et al. 2023
#'   used the mean temperature of each treatment between the previous
#'   sampling and that of the sample).
#' @param parameters A [flux_parameters()] result (default: the metabolic
#'   parameters of Brown et al. 2004 and the efficiencies of Polazzo et al.
#'   2023).
#' @param metabolic_types `NULL` (default) or a vector named by `node_id`
#'   with the metabolic type of consumers (one of
#'   `parameters$metabolism$metabolic_type`), replacing the default:
#'   unicellular for protists, fish, ectotherm vertebrates (amphibians and
#'   reptiles) and endotherms by class, invertebrate otherwise.
#' @param resource_types `NULL` (default) or a vector named by `node_id` with
#'   the type of resource of nodes (`"animal"`, `"producer"` or
#'   `"detritus"`), replacing the default given by their role.
#' @param sample_col,node_col,draw_col,biomass_col,mass_col Columns of the
#'   populations (the draw and mass columns may be absent).
#' @param temperature_col Column of the temperature in `temperature`.
#' @param permanent Nodes present in every sample (for example detritus);
#'   with `preferences = "biomass"` they still need a biomass in every
#'   sample.
#' @param preferences `"biomass"` (default): a consumer takes from each of
#'   its resources in proportion to the biomass of the resource
#'   (`bioms.prefs = TRUE` of *fluxweb*); `"equal"`: in equal parts.
#' @param na_rm If `FALSE` (default), a sample with a missing biomass or a
#'   consumer without a mass gets no fluxes (`status = "missing_values"`).
#'   If `TRUE`, the nodes with missing values are left out of that sample
#'   (counted in `n_missing`).
#' @param engine `"fluxweb"` (default): `fluxing()` of *fluxweb*;
#'   `"internal"`: the same equations solved within this package, which gives
#'   the same fluxes without *fluxweb*.
#' @param node_fluxes,link_fluxes If `TRUE`, also return the fluxes of every
#'   node or of every link in every sample and realization (large for many
#'   realizations).
#' @param quiet If `FALSE`, report progress.
#'
#' @return An object of class `bwebs_fluxes` (and `bwebs_indicators`, so that
#'   [summarize_indicators()] and [contrast_indicators()] summarise it), a
#'   list with:
#'   * `values`: one row per realization and sample (ordered by realization
#'     and sample): `draw` (0 for a metaweb and populations without draws),
#'     `sample_id`, `n_present` (nodes with a biomass), `n_missing` (nodes
#'     with missing values), `n_links` (links among the nodes of the
#'     fluxes), `total_flux`, `herbivory`, `detritivory` and `carnivory`
#'     (the sum of all the fluxes and of those from producers, detritus and
#'     animals), `metabolism` (the losses of the consumers), `n_unsupported`
#'     and `biomass_unsupported` (consumers left out because no path of
#'     links reaches them from a producer or detritus in that realization,
#'     and their biomass) and `status` (`"ok"`; `"no_links"`, all fluxes 0;
#'     `"unsolved"`, fluxes `NA`, when the web of that realization has no
#'     steady state, see Details; `"missing_values"`);
#'   * `nodes` (with `node_fluxes = TRUE`, otherwise `NULL`): `draw`,
#'     `sample_id`, `node_id`, `biomass`, `losses`, `intake` (fluxes from its
#'     resources), `consumed` (fluxes to its consumers) and `supported`;
#'   * `links` (with `link_fluxes = TRUE`, otherwise `NULL`): `draw`,
#'     `sample_id`, `resource`, `consumer` and `flux`;
#'   * `node_types`: the type of resource, metabolic type and efficiency (as
#'     a resource) of every node;
#'   * `settings`: the parameters, the engine, the preferences, the units
#'     and the settings of `x`.
#'
#'   The fluxes and the losses are in W (J s\eqn{^{-1}}{^-1}) per unit of
#'   the sampling basis of the biomasses.
#'
#' @details
#' Losses: the metabolic rate per gram of each population
#' ([flux_parameters()]), from its mean mass and the temperature of its
#' sample, times its biomass, summed by node (so that the life stages of a
#' taxon keep their own masses). Producers and detritus have no losses: the
#' fluxes from them are set by their consumers. Efficiencies depend on the
#' resource (`ef.level = "prey"` of *fluxweb*), so an omnivore
#' assimilates animals and algae with different efficiencies.
#'
#' The web of a sample in a realization has the nodes with a positive
#' biomass (and the `permanent` ones) and the links of the realization among
#' them. A consumer that no path of links reaches from a producer or
#' detritus (one left without resources in that realization, or eating only
#' such consumers) cannot be sustained by the web: it is left out of the
#' fluxes and counted in `n_unsupported`, rather than being treated as a
#' source of energy, as *fluxweb* would treat a node without resources.
#'
#' Webs without a steady state. Energy is lost at every step of
#' assimilation, so consumers that take much of their diet from one another,
#' in loops, cannot sustain their losses: the system then has no
#' non-negative solution and the realization gets `status = "unsolved"`.
#' This happens in realizations with few inferred links, where the links
#' that the taxonomic decay adds among consumers make up most of their
#' diets (grazers that only eat predators). The share of unsolved
#' realizations measures how many sampled webs are energetically
#' implausible; `summarize_indicators(na_rm = TRUE)` summarises the others.
#'
#' @section Assumptions:
#' The web is at steady state and its biomasses are those sampled. The
#' masses and biomasses are wet mass, the currency of the metabolic
#' intercepts (Brown et al. 2004), and every population has the metabolic
#' rate of an individual of its mean mass (with sizes spread around the
#' mean, the losses of a population differ from that). Fluxes are energy;
#' other currencies (carbon, nutrients) would need the content of each
#' resource and their own efficiencies.
#'
#' @references
#' Amador Crespo P, Rico A (2026) *bmassa*: Biomass of Aquatic Taxa from
#' Abundances and Body Sizes. R package version 0.5.0.
#' <https://github.com/STABLE-uv/bmassa>
#'
#' Gauzens B, Barnes A, Giling DP, Hines J, Jochum M, Lefcheck JS, Rosenbaum
#' B, Wang S, Brose U (2019) *fluxweb*: an R package to easily estimate energy
#' fluxes in food webs. *Methods in Ecology and Evolution* 10: 270-279.
#' \doi{10.1111/2041-210X.13109}
#'
#' Polazzo F, Hermann M, Crettaz-Minaglia M, Rico A (2023) Impacts of extreme
#' climatic events on trophic network complexity and multidimensional
#' stability. *Ecology* 104: e3951. \doi{10.1002/ecy.3951}
#'
#' @seealso [flux_parameters()], [node_masses()], [summarize_indicators()].
#' @export
#' @examples
#' nodes <- prepare_nodes(example_taxa)
#' mw <- build_metaweb(nodes, combine_registries(gray2015, gateway))
#' # the size term with the masses of the populations
#' draws <- sample_metawebs(mw, gray2015_gateway_calibration, n_draws = 20,
#'                          masses = node_masses(mw, example_populations), seed = 1)
#' # engine = "fluxweb", the default, gives the same fluxes with the package fluxweb
#' fl <- calc_fluxes(draws, example_populations, temperature = 15, engine = "internal")
#' fl
#' # na_rm = TRUE: the realizations without a steady state are left out
#' summarize_indicators(fl, indicators = c("total_flux", "herbivory", "detritivory"),
#'                      samples = example_samples, by = "treatment", na_rm = TRUE)
calc_fluxes <- function(x, populations, temperature, parameters = flux_parameters(),
                        metabolic_types = NULL, resource_types = NULL,
                        sample_col = "sample_id", node_col = "node_id", draw_col = "draw",
                        biomass_col = "biomass", mass_col = "mass",
                        temperature_col = "temperature", permanent = character(0),
                        preferences = "biomass", na_rm = FALSE, engine = "fluxweb",
                        node_fluxes = FALSE, link_fluxes = FALSE, quiet = TRUE) {
  if (!inherits(x, c("bwebs_draws", "bwebs_metaweb"))) {
    abort_input("`x` must be created with sample_metawebs() or build_metaweb().")
  }
  if (!inherits(parameters, "bwebs_flux_parameters")) {
    abort_input("`parameters` must be created with flux_parameters().")
  }
  for (a in c("sample_col", "node_col", "draw_col", "biomass_col", "mass_col",
              "temperature_col")) {
    check_string(get(a), arg = a)
  }
  check_choice(preferences, c("biomass", "equal"))
  check_flag(na_rm)
  check_choice(engine, c("fluxweb", "internal"))
  check_flag(node_fluxes)
  check_flag(link_fluxes)
  check_flag(quiet)
  if (engine == "fluxweb") {
    check_installed("fluxweb", "compute the fluxes with engine = \"fluxweb\"")
  }
  ids <- x$nodes$node_id
  role <- x$nodes$role
  if (!is.character(permanent) || anyNA(permanent) || length(setdiff(permanent, ids))) {
    abort_input(sprintf("`permanent` must contain node identifiers of the metaweb; unknown: %s.",
                        format_items(setdiff(permanent, ids))))
  }

  # types of the nodes
  types <- default_node_types(x$nodes)
  if (!is.null(resource_types)) {
    check_named_types(resource_types, ids, flux_resource_types, "resource_types")
    types$resource[names(resource_types)] <- resource_types
  }
  if (!is.null(metabolic_types)) {
    check_named_types(metabolic_types, ids, parameters$metabolism$metabolic_type,
                      "metabolic_types")
    not_consumer <- names(metabolic_types)[role[match(names(metabolic_types), ids)] != "consumer"]
    if (length(not_consumer)) {
      abort_input(sprintf(
        "Only consumers have metabolic losses; `metabolic_types` names other nodes: %s.",
        format_items(not_consumer)
      ))
    }
    types$metabolic[names(metabolic_types)] <- metabolic_types
  }
  unknown_types <- setdiff(stats::na.omit(types$metabolic), parameters$metabolism$metabolic_type)
  if (length(unknown_types)) {
    abort_input(sprintf("Metabolic types without parameters in `parameters`: %s.",
                        format_items(quote_items(unknown_types))))
  }
  efficiency <- parameters$efficiencies$efficiency[
    match(types$resource, parameters$efficiencies$resource_type)
  ]

  # populations and realizations
  sources <- population_sources(populations, ids, sample_col, node_col, draw_col, biomass_col,
                                mass_col)
  available <- population_n_draws(sources)
  pop_draws <- !is.na(available)
  if (inherits(x, "bwebs_draws")) {
    n_real <- ncol(x$z)
    if (pop_draws && available < n_real) {
      abort_input(sprintf(paste(
        "The populations have %d draws and `x` has %d realizations: give at least as many",
        "draws of the populations (or none, for the same populations in every realization)."),
        available, n_real))
    }
    labels <- seq_len(n_real)
  } else if (pop_draws) {
    n_real <- available
    labels <- seq_len(n_real)
  } else {
    n_real <- 1L
    labels <- 0L
  }
  samples <- sort(unique(unlist(lapply(sources, `[[`, "samples"), use.names = FALSE)),
                  method = "radix")
  n_s <- length(samples)
  n_ids <- length(ids)

  # temperature of each sample (degrees C)
  temp <- if (is.numeric(temperature) && length(temperature) == 1L) {
    check_number(temperature, lower = -273.15, include_lower = FALSE)
    rep(temperature, n_s)
  } else {
    check_data_frame(temperature, allow_empty = FALSE)
    check_columns(temperature, c(sample_col = sample_col, temperature_col = temperature_col),
                  arg = "temperature")
    check_numeric_column(temperature, temperature_col, lower = -273.15, include_lower = FALSE,
                         arg = "temperature")
    check_unique_key(temperature, sample_col, arg = "temperature")
    i <- match(samples, as.character(temperature[[sample_col]]))
    if (anyNA(i)) {
      abort_input(sprintf("`temperature` lacks samples of the populations: %s.",
                          format_items(samples[is.na(i)])))
    }
    temperature[[temperature_col]][i]
  }

  basal <- role != "consumer"
  perm <- ids %in% permanent
  rtype <- unname(types$resource)
  mtype <- unname(types$metabolic)
  scale <- preferences == "biomass"
  solve_flux <- if (engine == "fluxweb") flux_matrix_fluxweb else flux_matrix_internal
  n_out <- n_real * n_s
  v_present <- v_missing <- v_links <- v_uns <- integer(n_out)
  v_total <- v_herb <- v_detr <- v_carn <- v_met <- v_buns <- numeric(n_out)
  v_status <- character(n_out)
  nd_rows <- if (node_fluxes) vector("list", n_out)
  lk_rows <- if (link_fluxes) vector("list", n_out)
  adj <- if (!inherits(x, "bwebs_draws")) realization_adjacency(x, 1L)

  k <- 0L
  for (r in seq_len(n_real)) {
    if (!quiet) message(sprintf("Realization %d of %d", r, n_real))
    if (inherits(x, "bwebs_draws")) adj <- realization_adjacency(x, r)
    pd <- population_draw(sources, r)
    si <- match(pd$sample, samples)
    mt <- mtype[pd$node]
    living <- !is.na(mt) & !is.na(pd$biomass) & pd$biomass > 0
    rate <- numeric(length(pd$node))
    rate[living] <- metabolic_rate(pd$mass[living], temp[si[living]], mt[living], parameters)
    loss <- rate * pd$biomass
    cell <- pd$node + (si - 1L) * n_ids
    Bm <- Lm <- matrix(0, n_ids, n_s)
    sb <- rowsum(pd$biomass, cell)
    Bm[as.integer(rownames(sb))] <- sb[, 1L]
    sl <- rowsum(loss, cell)
    Lm[as.integer(rownames(sl))] <- sl[, 1L]
    for (s in seq_len(n_s)) {
      k <- k + 1L
      B <- Bm[, s]
      L <- Lm[, s]
      miss <- is.na(B) | is.na(L)
      detected <- !miss & B > 0
      v_present[k] <- sum(detected)
      v_missing[k] <- sum(miss)
      if (any(miss) && !na_rm) {
        v_status[k] <- "missing_values"
        v_links[k] <- NA_integer_
        v_uns[k] <- NA_integer_
        v_total[k] <- v_herb[k] <- v_detr[k] <- v_carn[k] <- v_met[k] <- v_buns[k] <- NA_real_
        next
      }
      pres <- which(detected | (perm & !miss))
      if (scale) {
        no_b <- pres[B[pres] <= 0]
        if (length(no_b)) {
          abort_input(sprintf(paste(
            "With preferences = \"biomass\" every node of a web needs a biomass: %s in",
            "sample '%s' have none. Give them a biomass in `populations` or use",
            "preferences = \"equal\"."), format_items(ids[no_b]), samples[s]))
        }
      }
      A <- adj[pres, pres, drop = FALSE]
      # consumers reached by a path of links from a producer or detritus
      sup <- basal[pres]
      repeat {
        new <- !sup & colSums(A[sup, , drop = FALSE]) > 0
        if (!any(new)) break
        sup <- sup | new
      }
      v_uns[k] <- sum(!sup)
      v_buns[k] <- sum(B[pres][!sup])
      A[!sup, ] <- 0
      A[, !sup] <- 0
      keep <- rowSums(A) > 0 | colSums(A) > 0
      kp <- pres[keep]
      Ak <- A[keep, keep, drop = FALSE]
      v_links[k] <- as.integer(sum(Ak))
      v_met[k] <- sum(L[pres][sup & !basal[pres]])
      Fm <- NULL
      if (v_links[k] == 0L) {
        v_status[k] <- "no_links"
      } else {
        Lk <- L[kp]
        Lk[basal[kp]] <- 0
        Fm <- solve_flux(Ak, B[kp], Lk, efficiency[kp], scale)
        if (is.null(Fm)) {
          v_status[k] <- "unsolved"
          v_total[k] <- v_herb[k] <- v_detr[k] <- v_carn[k] <- NA_real_
        } else {
          v_status[k] <- "ok"
          out <- rowSums(Fm)
          v_total[k] <- sum(out)
          v_herb[k] <- sum(out[rtype[kp] == "producer"])
          v_detr[k] <- sum(out[rtype[kp] == "detritus"])
          v_carn[k] <- sum(out[rtype[kp] == "animal"])
        }
      }
      if (node_fluxes) {
        intake <- consumed <- numeric(length(pres))
        if (!is.null(Fm)) {
          intake[keep] <- colSums(Fm)
          consumed[keep] <- rowSums(Fm)
        } else if (v_status[k] == "unsolved") {
          intake[keep] <- consumed[keep] <- NA_real_
        }
        nd_rows[[k]] <- list(node = pres, biomass = B[pres], losses = L[pres], intake = intake,
                             consumed = consumed, supported = sup)
      }
      if (link_fluxes && !is.null(Fm)) {
        w <- which(Ak > 0, arr.ind = TRUE)
        lk_rows[[k]] <- list(resource = kp[w[, 1L]], consumer = kp[w[, 2L]], flux = Fm[w])
      }
    }
  }
  values <- data.frame(
    draw = rep(labels, each = n_s), sample_id = rep(samples, times = n_real),
    n_present = v_present, n_missing = v_missing, n_links = v_links, total_flux = v_total,
    herbivory = v_herb, detritivory = v_detr, carnivory = v_carn, metabolism = v_met,
    n_unsupported = v_uns, biomass_unsupported = v_buns, status = v_status,
    stringsAsFactors = FALSE
  )
  draw_of <- rep(labels, each = n_s)
  sample_of <- rep(samples, times = n_real)
  stack_rows <- function(rows, cols, node_cols) {
    n <- vapply(rows, function(z) if (is.null(z)) 0L else length(z[[cols[1L]]]), integer(1))
    out <- data.frame(draw = rep(draw_of, n), sample_id = rep(sample_of, n),
                      stringsAsFactors = FALSE)
    for (cl in cols) {
      val <- unlist(lapply(rows, `[[`, cl), use.names = FALSE)
      if (is.null(val)) val <- if (cl %in% node_cols) integer(0) else numeric(0)
      out[[cl]] <- if (cl %in% node_cols) ids[val] else val
    }
    out
  }
  nodes_out <- if (node_fluxes) {
    z <- stack_rows(nd_rows, c("node", "biomass", "losses", "intake", "consumed", "supported"),
                    "node")
    names(z)[names(z) == "node"] <- "node_id"
    if (!nrow(z)) z$supported <- logical(0)
    z
  }
  links_out <- if (link_fluxes) {
    stack_rows(lk_rows, c("resource", "consumer", "flux"), c("resource", "consumer"))
  }
  xs <- x$settings
  draws_x <- inherits(x, "bwebs_draws")
  structure(
    list(values = values, nodes = nodes_out, links = links_out,
         node_types = data.frame(node_id = ids, resource_type = rtype, metabolic_type = mtype,
                                 efficiency = efficiency, stringsAsFactors = FALSE),
         settings = list(source = if (draws_x) "draws" else "metaweb", n_draws = n_real,
                         system = if (draws_x) xs$system else NA,
                         parameter_uncertainty = if (draws_x) xs$parameter_uncertainty else NA,
                         population_draws = pop_draws, parameters = parameters,
                         engine = engine, preferences = preferences, permanent = permanent,
                         na_rm = na_rm,
                         units = "W (J/s) per unit of the sampling basis of the biomasses")),
    class = c("bwebs_fluxes", "bwebs_indicators")
  )
}

#' @export
print.bwebs_fluxes <- function(x, ...) {
  v <- x$values
  s <- x$settings
  cat(sprintf("<bwebs_fluxes> %d sample(s) x %d realization(s) of %s\n",
              length(unique(v$sample_id)), s$n_draws,
              if (s$source == "draws") "sampled metawebs" else "the deterministic metaweb"))
  cat(sprintf("  biomasses: %s; parameters: %s; engine: %s\n",
              if (isTRUE(s$population_draws)) "drawn with each realization" else
                "the same in every realization",
              s$parameters$set, s$engine))
  st <- table(factor(v$status, levels = c("ok", "no_links", "unsolved", "missing_values")))
  st <- st[st > 0 | names(st) == "ok"]
  cat("  status: ", paste(names(st), st, sep = " = ", collapse = ", "), "\n", sep = "")
  ok <- v$status == "ok"
  if (any(ok)) {
    q <- stats::quantile(v$total_flux[ok], c(0.5, 0.025, 0.975), names = FALSE)
    cat(sprintf("  total flux: median %.3g W (%.3g to %.3g) per unit of the sampling basis\n",
                q[1], q[2], q[3]))
    cat(sprintf("  consumers left out without a path from a basal node: %.2g per web on average\n",
                mean(v$n_unsupported[ok])))
  }
  cat("  use summarize_indicators() for posterior summaries; values in $values\n")
  invisible(x)
}
