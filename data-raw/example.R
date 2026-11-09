# Example community -> data/example_taxa.rda, data/example_community.rda,
# data/example_samples.rda
#
# Source: invented for the package examples (no field data). The taxa are
# freshwater taxa present in the registry of Gray et al. (2015), classified as
# in that registry; Dugastella valentina (an Iberian atyid shrimp without
# records in the registry) takes the records of Asellus through a surrogate.
# The abundances are invented: six samples, three controls and three exposed
# to a stressor that removes some taxa. Licence: the same as the package.
# Packages used only here: usethis.

example_taxa <- data.frame(
  node_id = c("cocconeis_placentula", "gomphonema", "achnanthes_minutissima", "fragilaria",
              "cpom", "fpom",
              "baetis_rhodani", "ecdyonurus", "leuctra_nigra", "nemurella_pictetii",
              "simulium", "orthocladiinae", "tanypodinae", "asellus_aquaticus",
              "gammarus_pulex", "dugastella_valentina", "ancylus_fluviatilis", "pisidium",
              "hydropsyche", "rhyacophila_dorsalis", "plectrocnemia_conspersa",
              "sialis_fuliginosa", "cordulegaster_boltonii"),
  taxon = c("Cocconeis placentula", "Gomphonema", "Achnanthes minutissima", "Fragilaria",
            "CPOM", "FPOM",
            "Baetis rhodani", "Ecdyonurus", "Leuctra nigra", "Nemurella pictetii",
            "Simulium", "Orthocladiinae", "Tanypodinae", "Asellus aquaticus",
            "Gammarus pulex", "Dugastella valentina", "Ancylus fluviatilis", "Pisidium",
            "Hydropsyche", "Rhyacophila dorsalis", "Plectrocnemia conspersa",
            "Sialis fuliginosa", "Cordulegaster boltonii"),
  rank = c("species", "genus", "species", "genus", "non_taxonomic", "non_taxonomic",
           "species", "genus", "species", "species", "genus", "subfamily", "subfamily",
           "species", "species", "species", "species", "genus", "genus", "species",
           "species", "species", "species"),
  role = c(rep("producer", 4), rep("detritus", 2), rep("consumer", 17)),
  node_type = c(rep("taxon", 4), rep("non_taxonomic", 2), rep("taxon", 17)),
  genus = c("Cocconeis", "Gomphonema", "Achnanthes", "Fragilaria", NA, NA,
            "Baetis", "Ecdyonurus", "Leuctra", "Nemurella", "Simulium", NA, NA, "Asellus",
            "Gammarus", "Dugastella", "Ancylus", "Pisidium", "Hydropsyche", "Rhyacophila",
            "Plectrocnemia", "Sialis", "Cordulegaster"),
  subfamily = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, "Orthocladiinae", "Tanypodinae",
                NA, NA, NA, NA, NA, NA, NA, NA, NA, NA),
  family = c("Cocconeidaceae", "Gomphonemataceae", "Achnanthaceae", "Fragilariaceae", NA, NA,
             "Baetidae", "Heptageniidae", "Leuctridae", "Nemouridae", "Simuliidae",
             "Chironomidae", "Chironomidae", "Asellidae", "Gammaridae", "Atyidae",
             "Ancylidae", "Pisidiidae", "Hydropsychidae", "Rhyacophilidae",
             "Polycentropodidae", "Sialidae", "Cordulegastridae"),
  order = c("Achnanthales", "Cymbellales", "Achnanthales", "Fragilariales", NA, NA,
            "Ephemeroptera", "Ephemeroptera", "Plecoptera", "Plecoptera", "Diptera",
            "Diptera", "Diptera", "Isopoda", "Amphipoda", "Decapoda", "Basommatophora",
            "Veneroida", "Trichoptera", "Trichoptera", "Trichoptera", "Megaloptera",
            "Odonata"),
  class = c("Bacillariophyceae", "Bacillariophyceae", "Bacillariophyceae",
            "Coscinodiscophyceae", NA, NA, rep("Insecta", 7), "Malacostraca", "Malacostraca",
            "Malacostraca", "Gastropoda", "Bivalvia", rep("Insecta", 5)),
  phylum = c(rep("Ochrophyta", 4), NA, NA, rep("Arthropoda", 10), "Mollusca", "Mollusca",
             rep("Arthropoda", 5)),
  alias = c(NA, NA, NA, NA, "detritus", "detritus", rep(NA, 17)),
  consumer_surrogate = c(rep(NA, 15), "Asellus", rep(NA, 7)),
  resource_surrogate = c(rep(NA, 15), "Asellus", rep(NA, 7)),
  stringsAsFactors = FALSE
)

# Invented abundances (individuals per sample, or cells for diatoms);
# exposed samples lose mayflies, stoneflies, amphipods and the shrimp.
abund <- rbind(
  control_1 = c(120, 80, 60, 40, NA, NA, 12, 6, 9, 4, 15, 30, 8, 10, 14, 6, 3, 5, 4, 3, 2, 1, 1),
  control_2 = c(100, 90, 50, 45, NA, NA, 10, 5, 7, 3, 12, 28, 9, 12, 11, 4, 2, 6, 5, 2, 2, 1, 0),
  control_3 = c(130, 70, 65, 35, NA, NA, 14, 4, 8, 5, 18, 33, 7, 9, 16, 5, 4, 4, 3, 4, 1, 2, 1),
  exposed_1 = c(150, 60, 80, 50, NA, NA, 0, 0, 2, 0, 10, 45, 10, 8, 0, 0, 3, 6, 4, 2, 2, 1, 0),
  exposed_2 = c(160, 55, 85, 40, NA, NA, 1, 0, 0, 0, 9, 50, 12, 6, 0, 0, 2, 5, 3, 1, 1, 0, 0),
  exposed_3 = c(140, 65, 75, 55, NA, NA, 0, 1, 1, 0, 11, 41, 9, 7, 2, 0, 3, 7, 5, 2, 3, 1, 1)
)
colnames(abund) <- example_taxa$node_id
# Detritus is not sampled: it is a permanent resource (calc_indicators(permanent = ))
abund <- abund[, !colnames(abund) %in% c("cpom", "fpom")]
example_community <- data.frame(
  sample_id = rep(rownames(abund), times = ncol(abund)),
  node_id = rep(colnames(abund), each = nrow(abund)),
  value = as.vector(abund),
  stringsAsFactors = FALSE
)
example_community <- example_community[order(example_community$sample_id,
                                             example_community$node_id, method = "radix"), ]
rownames(example_community) <- NULL

example_samples <- data.frame(
  sample_id = rownames(abund),
  treatment = rep(c("control", "exposed"), each = 3),
  replicate = rep(1:3, times = 2),
  stringsAsFactors = FALSE
)

usethis::use_data(example_taxa, example_community, example_samples, overwrite = TRUE)

# Populations for the energy fluxes -> data/example_populations.rda
#
# Invented, like the abundances: biomass (g of wet mass per sample) and mean
# individual wet mass (g) of every node of each sample. Consumers: the
# abundances times a typical wet mass of each taxon; diatoms: 2 mg of
# periphyton per counted unit; detritus: 5 g of CPOM and 2 g of FPOM in every
# sample. Producers and detritus have no individual mass.
mean_mass <- c(baetis_rhodani = 3e-3, ecdyonurus = 1.5e-2, leuctra_nigra = 1e-3,
               nemurella_pictetii = 2e-3, simulium = 1e-3, orthocladiinae = 2e-4,
               tanypodinae = 5e-4, asellus_aquaticus = 5e-3, gammarus_pulex = 2e-2,
               dugastella_valentina = 0.1, ancylus_fluviatilis = 2e-2, pisidium = 1e-2,
               hydropsyche = 2e-2, rhyacophila_dorsalis = 3e-2, plectrocnemia_conspersa = 1.5e-2,
               sialis_fuliginosa = 5e-2, cordulegaster_boltonii = 0.5)
producers <- example_taxa$node_id[example_taxa$role == "producer"]
pop <- example_community
pop$mass <- unname(mean_mass[pop$node_id])
pop$biomass <- ifelse(pop$node_id %in% producers, pop$value * 2e-3, pop$value * pop$mass)
detritus <- data.frame(sample_id = rep(example_samples$sample_id, each = 2),
                       node_id = rep(c("cpom", "fpom"), times = nrow(example_samples)),
                       biomass = rep(c(5, 2), times = nrow(example_samples)),
                       mass = NA_real_, stringsAsFactors = FALSE)
example_populations <- rbind(pop[c("sample_id", "node_id", "biomass", "mass")], detritus)
example_populations <- example_populations[order(example_populations$sample_id,
                                                 example_populations$node_id,
                                                 method = "radix"), ]
rownames(example_populations) <- NULL

usethis::use_data(example_populations, overwrite = TRUE)
