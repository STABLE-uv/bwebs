# A web whose fluxes can be computed by hand: an alga and detritus eaten by a
# grazer (an insect), which is eaten by a fish.
flux_web <- function() {
  ids <- c("alga", "det", "graz", "fish")
  nodes <- data.frame(node_id = ids, taxon = ids,
                      role = c("producer", "detritus", "consumer", "consumer"),
                      class = c("Bacillariophyceae", NA, "Insecta", "Actinopterygii"),
                      phylum = c("Ochrophyta", NA, "Arthropoda", "Chordata"),
                      stringsAsFactors = FALSE)
  A <- matrix(0L, 4, 4, dimnames = list(ids, ids))
  A["alga", "graz"] <- A["det", "graz"] <- A["graz", "fish"] <- 1L
  structure(list(nodes = nodes, adjacency = A, settings = list()), class = "bwebs_metaweb")
}

flux_populations <- function() {
  data.frame(sample_id = "s1", node_id = c("alga", "det", "graz", "fish"),
             biomass = c(2, 6, 0.5, 1), mass = c(NA, NA, 1e-3, 1), stringsAsFactors = FALSE)
}

# Draws of the web of flux_web(): draw 1 with every link, draw 2 without the
# link grazer -> fish.
flux_draws <- function() {
  mw <- flux_web()
  links <- data.frame(resource = c("alga", "det", "graz"), consumer = c("graz", "graz", "fish"),
                      stringsAsFactors = FALSE)
  z <- cbind(c(TRUE, TRUE, TRUE), c(TRUE, TRUE, FALSE))
  structure(list(links = links, z = z, background = NULL, nodes = mw$nodes,
                 settings = list(n_draws = 2L, system = "random", parameter_uncertainty = TRUE)),
            class = "bwebs_draws")
}

# Fluxes of flux_web() by hand, at temperature `t` (degrees C), with the
# activation energy `e` (that of the default parameters).
hand_fluxes <- function(t = 20, equal = FALSE, e = 0.63) {
  kT <- 8.617333262e-5 * (t + 273.15)
  l_graz <- exp(17.17 - 0.29 * log(1e-3) - e / kT) * 0.5
  l_fish <- exp(18.47 - 0.29 * log(1) - e / kT) * 1
  w_alga <- if (equal) 0.5 else 2 / 8
  f_fish <- l_fish / 0.906
  f_graz <- (l_graz + f_fish) / (w_alga * 0.77 + (1 - w_alga) * 0.158)
  list(total = f_graz + f_fish, herbivory = w_alga * f_graz, detritivory = (1 - w_alga) * f_graz,
       carnivory = f_fish, metabolism = l_graz + l_fish, l_graz = l_graz, f_graz = f_graz)
}

test_that("flux_parameters() gives the published sets and replaces values", {
  p <- flux_parameters()
  expect_s3_class(p, "bwebs_flux_parameters")
  # by default the activation energy with which Brown et al. (2004) fitted
  # the intercepts
  expect_identical(p$set, "brown_2004")
  expect_identical(p$activation_energy, 0.63)
  expect_equal(p$metabolism$intercept[match(c("unicellular", "invertebrate", "fish"),
                                            p$metabolism$metabolic_type)],
               c(15.85, 17.17, 18.47))
  expect_true(all(p$metabolism$exponent == -0.29))
  expect_equal(p$efficiencies$efficiency, c(0.906, 0.77, 0.158))
  q <- flux_parameters(metabolism = data.frame(metabolic_type = c("invertebrate", "amphibian"),
                                               intercept = c(17, 18), exponent = -0.3),
                       activation_energy = 0.65, efficiencies = c(detritus = 0.2))
  expect_identical(q$metabolism$intercept[q$metabolism$metabolic_type == "invertebrate"], 17)
  expect_true("amphibian" %in% q$metabolism$metabolic_type)
  expect_identical(q$activation_energy, 0.65)
  expect_identical(q$efficiencies$efficiency[3], 0.2)
  expect_identical(q$efficiencies$source[3], "user")
  pz <- flux_parameters("polazzo_2023")
  expect_identical(pz$activation_energy, 0.69)
  expect_identical(pz$metabolism, p$metabolism)
  expect_identical(pz$efficiencies, p$efficiencies)
  # the intercepts of Brown et al. (2004) with their activation energy: about
  # 3 mW per g for an invertebrate of 1 mg at 20 degrees C; 0.3 with 0.69 eV
  expect_equal(metabolic_rate(1e-3, 20, "invertebrate", p), 3.13e-3, tolerance = 0.01)
  expect_equal(metabolic_rate(1e-3, 20, "invertebrate", pz) /
                 metabolic_rate(1e-3, 20, "invertebrate", p),
               exp(-0.06 / (8.617333262e-5 * 293.15)))
  expect_error(flux_parameters("other"), class = "bwebs_error_input")
  expect_error(flux_parameters(efficiencies = c(plant = 0.5)), class = "bwebs_error_input")
  expect_error(flux_parameters(efficiencies = c(animal = 1.5)), class = "bwebs_error_input")
  expect_error(flux_parameters(activation_energy = -1), class = "bwebs_error_input")
  expect_error(flux_parameters(metabolism = data.frame(metabolic_type = c("a", "a"),
                                                       intercept = 1, exponent = 1)),
               class = "bwebs_error_input")
})

test_that("the metabolic rate follows Brown et al. (2004) on the log scale", {
  p <- flux_parameters()
  kT <- 8.617333262e-5 * (15 + 273.15)
  expect_equal(metabolic_rate(c(1e-3, 2), 15, c("invertebrate", "fish"), p),
               exp(c(17.17, 18.47) - 0.29 * log(c(1e-3, 2)) - 0.63 / kT))
  # warmer is faster, larger is slower per gram
  expect_gt(metabolic_rate(1e-3, 25, "invertebrate", p), metabolic_rate(1e-3, 15, "invertebrate", p))
  expect_gt(metabolic_rate(1e-3, 15, "invertebrate", p), metabolic_rate(1, 15, "invertebrate", p))
})

test_that("the internal solver gives the fluxes of a chain and an omnivore by hand", {
  ids <- c("p", "c1", "c2")
  A <- matrix(0, 3, 3, dimnames = list(ids, ids))
  A["p", "c1"] <- A["p", "c2"] <- A["c1", "c2"] <- 1
  B <- c(10, 2, 1)
  L <- c(0, 0.3, 0.2)
  e <- c(0.77, 0.906, 0.906)
  Fm <- flux_matrix_internal(A, B, L, e, scale = TRUE)
  wp <- 10 / 12
  f2 <- 0.2 / (wp * 0.77 + (1 - wp) * 0.906)
  f1 <- (0.3 + (1 - wp) * f2) / 0.77
  expect_equal(Fm["p", "c2"], wp * f2)
  expect_equal(Fm["c1", "c2"], (1 - wp) * f2)
  expect_equal(Fm["p", "c1"], f1)
  # balance: what each consumer assimilates is its losses plus what is eaten of it
  assim <- colSums(Fm * e)
  expect_equal(unname(assim[2:3]), c(0.3 + Fm["c1", "c2"], 0.2))
  eq <- flux_matrix_internal(A, B, L, e, scale = FALSE)
  expect_equal(eq["p", "c2"], eq["c1", "c2"])
})

test_that("calc_fluxes() reproduces the fluxes of a web by hand", {
  fl <- calc_fluxes(flux_web(), flux_populations(), temperature = 20, engine = "internal",
                    node_fluxes = TRUE, link_fluxes = TRUE)
  expect_s3_class(fl, c("bwebs_fluxes", "bwebs_indicators"))
  v <- fl$values
  h <- hand_fluxes(20)
  expect_identical(v$draw, 0L)
  expect_identical(v$status, "ok")
  expect_identical(v$n_present, 4L)
  expect_identical(v$n_links, 3L)
  expect_equal(v$total_flux, h$total)
  expect_equal(v$herbivory, h$herbivory)
  expect_equal(v$detritivory, h$detritivory)
  expect_equal(v$carnivory, h$carnivory)
  expect_equal(v$metabolism, h$metabolism)
  expect_identical(fl$node_types$metabolic_type, c(NA, NA, "invertebrate", "fish"))
  expect_equal(sum(fl$links$flux), v$total_flux)
  nd <- fl$nodes
  expect_equal(nd$intake[nd$node_id == "graz"], h$f_graz)
  expect_equal(nd$losses[nd$node_id == "graz"], h$l_graz)
  expect_equal(nd$consumed[nd$node_id == "alga"], h$herbivory)
  equal <- calc_fluxes(flux_web(), flux_populations(), temperature = 20, engine = "internal",
                       preferences = "equal")
  expect_equal(equal$values$herbivory, hand_fluxes(20, equal = TRUE)$herbivory)
  pz <- calc_fluxes(flux_web(), flux_populations(), temperature = 20, engine = "internal",
                    parameters = flux_parameters("polazzo_2023"))
  expect_equal(pz$values$total_flux, hand_fluxes(20, e = 0.69)$total)
  expect_equal(pz$values$metabolism, hand_fluxes(20, e = 0.69)$metabolism)
})

test_that("losses are summed over the populations of a node", {
  pops <- flux_populations()
  stages <- rbind(pops[pops$node_id != "graz", ],
                  data.frame(sample_id = "s1", node_id = "graz", biomass = c(0.2, 0.3),
                             mass = c(1e-4, 1e-3)))
  fl <- calc_fluxes(flux_web(), stages, temperature = 20, engine = "internal",
                    node_fluxes = TRUE)
  kT <- 8.617333262e-5 * 293.15
  l_graz <- sum(exp(17.17 - 0.29 * log(c(1e-4, 1e-3)) - 0.63 / kT) * c(0.2, 0.3))
  expect_equal(fl$nodes$losses[fl$nodes$node_id == "graz"], l_graz)
  expect_equal(fl$nodes$biomass[fl$nodes$node_id == "graz"], 0.5)
})

test_that("the temperature of each sample sets its losses", {
  pops <- rbind(flux_populations(), transform(flux_populations(), sample_id = "s2"))
  temp <- data.frame(sample_id = c("s2", "s1"), temperature = c(25, 20))
  fl <- calc_fluxes(flux_web(), pops, temperature = temp, engine = "internal")
  expect_equal(fl$values$total_flux, c(hand_fluxes(20)$total, hand_fluxes(25)$total))
  expect_error(calc_fluxes(flux_web(), pops, temperature = temp[1, ], engine = "internal"),
               "lacks samples", class = "bwebs_error_input")
  expect_error(calc_fluxes(flux_web(), pops, temperature = -300, engine = "internal"),
               class = "bwebs_error_input")
})

test_that("consumers without a path from a basal node are left out", {
  mw <- flux_web()
  ids <- c(mw$nodes$node_id, "orph")
  mw$nodes <- rbind(mw$nodes, data.frame(node_id = "orph", taxon = "orph", role = "consumer",
                                         class = "Insecta", phylum = "Arthropoda"))
  A <- matrix(0L, 5, 5, dimnames = list(ids, ids))
  A[1:4, 1:4] <- mw$adjacency
  A["orph", "fish"] <- 1L
  mw$adjacency <- A
  pops <- rbind(flux_populations(),
                data.frame(sample_id = "s1", node_id = "orph", biomass = 0.4, mass = 1e-3))
  fl <- calc_fluxes(mw, pops, temperature = 20, engine = "internal")
  v <- fl$values
  expect_identical(v$n_unsupported, 1L)
  expect_equal(v$biomass_unsupported, 0.4)
  expect_identical(v$n_links, 3L)
  expect_equal(v$total_flux, hand_fluxes(20)$total)
  expect_equal(v$metabolism, hand_fluxes(20)$metabolism)
})

test_that("missing values block a sample unless na_rm = TRUE", {
  pops <- flux_populations()
  pops$biomass[3] <- NA
  blocked <- calc_fluxes(flux_web(), pops, temperature = 20, engine = "internal")
  expect_identical(blocked$values$status, "missing_values")
  expect_true(is.na(blocked$values$total_flux))
  expect_identical(blocked$values$n_missing, 1L)
  kept <- calc_fluxes(flux_web(), pops, temperature = 20, engine = "internal", na_rm = TRUE)
  # without the grazer the fish has no path from a basal node
  expect_identical(kept$values$status, "no_links")
  expect_identical(kept$values$n_unsupported, 1L)
  expect_identical(kept$values$total_flux, 0)
  # a consumer without a mass is a missing value too
  pops <- flux_populations()
  pops$mass[4] <- NA
  expect_identical(calc_fluxes(flux_web(), pops, temperature = 20,
                               engine = "internal")$values$status, "missing_values")
  # a zero biomass is an absence
  pops <- flux_populations()
  pops$biomass[4] <- 0
  pops$mass[4] <- NA
  z <- calc_fluxes(flux_web(), pops, temperature = 20, engine = "internal")
  expect_identical(z$values$status, "ok")
  expect_identical(z$values$n_present, 3L)
  expect_identical(z$values$carnivory, 0)
})

test_that("realization d has the links of draw d and the populations of draw d", {
  dr <- flux_draws()
  fl <- calc_fluxes(dr, flux_populations(), temperature = 20, engine = "internal")
  v <- fl$values
  expect_identical(v$draw, 1:2)
  h <- hand_fluxes(20)
  expect_equal(v$total_flux[1], h$total)
  # draw 2: the fish has no resources and is left out
  expect_identical(v$n_unsupported, c(0L, 1L))
  expect_equal(v$total_flux[2], h$l_graz / (0.25 * 0.77 + 0.75 * 0.158))
  # populations with draws: the grazer doubles its biomass in draw 2
  p1 <- flux_populations()
  p2 <- flux_populations()
  p2$biomass[3] <- 1
  pops <- rbind(cbind(draw = 1L, p1), cbind(draw = 2L, p2), cbind(draw = 3L, p1))
  fd <- calc_fluxes(dr, pops, temperature = 20, engine = "internal")
  expect_equal(fd$values$total_flux[1], h$total)
  expect_equal(fd$values$total_flux[2], 2 * h$l_graz / (0.25 * 0.77 + 0.75 * 0.158))
  expect_true(fd$settings$population_draws)
  expect_error(calc_fluxes(dr, pops[pops$draw == 1, ], temperature = 20, engine = "internal"),
               "1 draws", class = "bwebs_error_input")
  expect_error(calc_fluxes(dr, pops[pops$draw != 2, ], temperature = 20, engine = "internal"),
               "no rows for draws", class = "bwebs_error_input")
  # a metaweb with populations drawn: one realization per draw of the populations
  fm <- calc_fluxes(flux_web(), pops, temperature = 20, engine = "internal")
  expect_identical(fm$values$draw, 1:3)
  s <- summarize_indicators(fm, indicators = "total_flux")
  expect_identical(s$uncertainty, "biomass")
  expect_false(is.na(s$sd))
  s2 <- summarize_indicators(fd, indicators = "total_flux")
  expect_identical(s2$uncertainty, "link_inference+between_system+parameter+biomass")
})

test_that("simulations of bmassa give the populations of each draw", {
  sim <- structure(list(
    rows = data.frame(sample_id = "s1", node_id = c("graz", "fish"), abundance = c(500, 1),
                      currency = "wet_mass", unit = c("mg", "g"),
                      biomass_unit = c("mg/tank", "g/tank"), stringsAsFactors = FALSE),
    ind_mass = matrix(c(1, 1, 2, 1, 1, 1), 2, 3)
  ), class = "bmassa_simulation")
  basal <- flux_populations()[1:2, c("sample_id", "node_id", "biomass")]
  fl <- calc_fluxes(flux_draws(), list(sim, basal), temperature = 20, engine = "internal")
  h <- hand_fluxes(20)
  expect_equal(fl$values$total_flux[1], h$total)
  # draw 2: grazers of 2 mg (biomass 1 g), without the link to the fish
  l_graz2 <- exp(17.17 - 0.29 * log(2e-3) - 0.63 / (8.617333262e-5 * 293.15)) * 1
  expect_equal(fl$values$total_flux[2], l_graz2 / (0.25 * 0.77 + 0.75 * 0.158))
  m <- node_masses(flux_web(), list(sim, basal))
  expect_identical(dim(m), c(4L, 3L))
  expect_equal(unname(m[c("graz", "fish"), 2]), c(2e-3, 1))
  expect_true(all(is.na(m[c("alga", "det"), ])))
  dry <- sim
  dry$rows$currency <- "dry_mass"
  expect_error(calc_fluxes(flux_web(), list(dry, basal), temperature = 20, engine = "internal"),
               "wet mass", class = "bwebs_error_input")
  mixed <- sim
  mixed$rows$biomass_unit <- c("mg/tank", "g/m2")
  expect_error(calc_fluxes(flux_web(), list(mixed, basal), temperature = 20,
                           engine = "internal"),
               "different bases", class = "bwebs_error_input")
  bad <- sim
  bad$rows$node_id[1] <- NA
  expect_error(calc_fluxes(flux_web(), list(bad, basal), temperature = 20, engine = "internal"),
               class = "bwebs_error_input")
})

test_that("node_masses() pools the individuals of every population", {
  pops <- rbind(flux_populations(),
                data.frame(sample_id = "s2", node_id = "graz", biomass = c(0.2, 0),
                           mass = c(2e-3, 5e-3)))
  m <- node_masses(flux_web(), pops)
  expect_named(m, c("alga", "det", "graz", "fish"))
  expect_equal(m[["graz"]], 0.7 / (0.5 / 1e-3 + 0.2 / 2e-3))
  expect_true(is.na(m[["alga"]]))
  pd <- rbind(cbind(draw = 1L, flux_populations()), cbind(draw = 2L, flux_populations()))
  pd$mass[pd$draw == 2 & pd$node_id == "fish"] <- 2
  md <- node_masses(flux_web(), pd)
  expect_identical(dim(md), c(4L, 2L))
  expect_equal(md["fish", ], c(1, 2))
  expect_identical(ncol(node_masses(flux_web(), pd, n_draws = 1)), 1L)
  expect_error(node_masses(flux_web(), pd, n_draws = 3), class = "bwebs_error_input")
})

test_that("metabolic and resource types have defaults and can be replaced", {
  nodes <- data.frame(node_id = c("cil", "fish", "frog", "duck", "snail", "alg", "det"),
                      role = c(rep("consumer", 5), "producer", "detritus"),
                      class = c("Spirotrichea", "Actinopterygii", "Amphibia", "Aves", "Gastropoda",
                                "Chlorophyceae", NA),
                      phylum = c("Ciliophora", "Chordata", "Chordata", "Chordata", "Mollusca",
                                 "Chlorophyta", NA), stringsAsFactors = FALSE)
  ty <- default_node_types(nodes)
  expect_identical(unname(ty$metabolic),
                   c("unicellular", "fish", "ectotherm_vertebrate", "endotherm", "invertebrate",
                     NA, NA))
  expect_identical(unname(ty$resource),
                   c(rep("animal", 5), "producer", "detritus"))
  fl <- calc_fluxes(flux_web(), flux_populations(), temperature = 20, engine = "internal",
                    metabolic_types = c(fish = "invertebrate"),
                    resource_types = c(det = "producer"))
  expect_identical(fl$node_types$metabolic_type[4], "invertebrate")
  expect_equal(fl$node_types$efficiency[2], 0.77)
  expect_error(calc_fluxes(flux_web(), flux_populations(), temperature = 20, engine = "internal",
                           metabolic_types = c(alga = "invertebrate")),
               "Only consumers", class = "bwebs_error_input")
  expect_error(calc_fluxes(flux_web(), flux_populations(), temperature = 20, engine = "internal",
                           metabolic_types = c(fish = "shark")), class = "bwebs_error_input")
})

test_that("permanent nodes need a biomass for preferences scaled by biomass", {
  pops <- flux_populations()[-2, ]
  expect_error(calc_fluxes(flux_web(), pops, temperature = 20, engine = "internal",
                           permanent = "det"),
               "needs a biomass", class = "bwebs_error_input")
  fl <- calc_fluxes(flux_web(), pops, temperature = 20, engine = "internal", permanent = "det",
                    preferences = "equal")
  expect_equal(fl$values$total_flux, hand_fluxes(20, equal = TRUE)$total)
  expect_identical(fl$values$n_present, 3L)
})

test_that("a steady state without solution is reported, not an error", {
  local_mocked_bindings(flux_matrix_internal = function(...) NULL)
  fl <- calc_fluxes(flux_web(), flux_populations(), temperature = 20, engine = "internal",
                    node_fluxes = TRUE)
  expect_identical(fl$values$status, "unsolved")
  expect_true(is.na(fl$values$total_flux))
  expect_false(is.na(fl$values$metabolism))
  expect_true(all(is.na(fl$nodes$intake[fl$nodes$node_id %in% c("graz", "fish")])))
})

test_that("calc_fluxes() rejects invalid input with classed errors", {
  mw <- flux_web()
  pops <- flux_populations()
  expect_error(calc_fluxes(list(), pops, 20, engine = "internal"), class = "bwebs_error_input")
  expect_error(calc_fluxes(mw, pops, 20, parameters = list(), engine = "internal"),
               class = "bwebs_error_input")
  expect_error(calc_fluxes(mw, pops[0, ], 20, engine = "internal"), class = "bwebs_error_input")
  bad <- pops
  bad$node_id[1] <- "nope"
  expect_error(calc_fluxes(mw, bad, 20, engine = "internal"), "not in the metaweb",
               class = "bwebs_error_input")
  bad <- pops
  bad$biomass[1] <- -1
  expect_error(calc_fluxes(mw, bad, 20, engine = "internal"), class = "bwebs_error_input")
  expect_error(calc_fluxes(mw, pops, 20, engine = "other"), class = "bwebs_error_input")
  expect_error(calc_fluxes(mw, pops, 20, engine = "internal", permanent = "nope"),
               class = "bwebs_error_input")
  expect_error(calc_fluxes(mw, list(pops, "a"), 20, engine = "internal"),
               class = "bwebs_error_input")
  local_mocked_bindings(is_installed = function(pkg) FALSE)
  expect_error(calc_fluxes(mw, pops, 20), class = "bwebs_error_missing_package")
})

test_that("the draws keep the classification of the nodes", {
  mw <- build_metaweb(prepare_nodes(tiny_taxa(), rule = tiny_rule()), tiny_registry())
  cal <- fake_calibration(c(genus = 0, family = 0, order = 0))
  d <- sample_metawebs(mw, cal, n_draws = 3, seed = 1)
  expect_true(all(c("node_id", "role", "class", "family") %in% names(d$nodes)))
  pops <- data.frame(sample_id = "s1", node_id = c("cocc", "baet", "perl"),
                     biomass = c(1, 0.2, 0.1), mass = c(NA, 1e-3, 1e-2))
  fl <- calc_fluxes(d, pops, temperature = 15, engine = "internal")
  expect_identical(fl$values$draw, 1:3)
  expect_true(all(fl$values$status %in% c("ok", "no_links")))
})

test_that("the internal engine and fluxweb give the same fluxes", {
  skip_if_not_installed("fluxweb")
  nodes <- prepare_nodes(example_taxa)
  mw <- build_metaweb(nodes, combine_registries(gray2015, gateway))
  d <- sample_metawebs(mw, gray2015_gateway_calibration, n_draws = 5, seed = 1)
  a <- calc_fluxes(d, example_populations, temperature = 15, engine = "internal",
                   link_fluxes = TRUE)
  b <- calc_fluxes(d, example_populations, temperature = 15, engine = "fluxweb",
                   link_fluxes = TRUE)
  expect_equal(a$values[a$values$status == "ok", c("total_flux", "herbivory", "carnivory")],
               b$values[a$values$status == "ok", c("total_flux", "herbivory", "carnivory")])
  expect_equal(a$links$flux, b$links$flux)
})

test_that("print methods are stable", {
  expect_snapshot(print(flux_parameters()))
  expect_snapshot(print(calc_fluxes(flux_draws(), flux_populations(), temperature = 20,
                                    engine = "internal")))
})
