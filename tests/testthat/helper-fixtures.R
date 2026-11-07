# Small hand-made registry and nodes whose metaweb can be derived by hand
# (see test-metaweb.R for the expected links).

tiny_registry <- function() {
  rec <- data.frame(
    study_id = c("s1", "s1", "s2", "s2", "s3", "s3"),
    resource = c("Cocconeis placentula", "Baetis rhodani", "Navicula gregaria", "Diptera",
                 "CPOM", "CPOM"),
    consumer = c("Baetis rhodani", "Perla bipunctata", "Ecdyonurus venosus",
                 "Gerris lacustris", "Asellus aquaticus", "Simulium ornatum"),
    res_genus = c("Cocconeis", "Baetis", "Navicula", NA, NA, NA),
    res_subfamily = NA_character_,
    res_family = c("Cocconeidaceae", "Baetidae", "Naviculaceae", NA, NA, NA),
    res_order = c("Achnanthales", "Ephemeroptera", "Naviculales", "Diptera", NA, NA),
    res_class = c("Bacillariophyceae", "Insecta", "Bacillariophyceae", "Insecta", NA, NA),
    con_genus = c("Baetis", "Perla", "Ecdyonurus", "Gerris", "Asellus", "Simulium"),
    con_subfamily = NA_character_,
    con_family = c("Baetidae", "Perlidae", "Heptageniidae", "Gerridae", "Asellidae", "Simuliidae"),
    con_order = c("Ephemeroptera", "Plecoptera", "Ephemeroptera", "Hemiptera", "Isopoda", "Diptera"),
    con_class = c("Insecta", "Insecta", "Insecta", "Insecta", "Malacostraca", "Insecta"),
    res_category = c("producer", "invertebrate", "producer", "invertebrate", NA, NA),
    stringsAsFactors = FALSE
  )
  as_registry(rec, registry = "gray2015")
}

tiny_taxa <- function() {
  data.frame(
    node_id = c("cocc", "navi", "cpom", "baet", "ecdy", "perl", "simu", "gerr", "asel"),
    taxon = c("Cocconeis pediculus", "Navicula", "CPOM", "Baetis muticus", "Rhithrogena",
              "Perla marginata", "Simulium", "Gerris", "Asellus"),
    rank = c("species", "genus", "non_taxonomic", "species", "genus", "species", "genus",
             "genus", "genus"),
    role = c("producer", "producer", "detritus", rep("consumer", 6)),
    node_type = c("taxon", "taxon", "non_taxonomic", rep("taxon", 6)),
    genus = c("Cocconeis", "Navicula", NA, "Baetis", "Rhithrogena", "Perla", "Simulium",
              "Gerris", "Asellus"),
    subfamily = NA_character_,
    family = c("Cocconeidaceae", "Naviculaceae", NA, "Baetidae", "Heptageniidae", "Perlidae",
               "Simuliidae", "Gerridae", "Asellidae"),
    order = c("Achnanthales", "Naviculales", NA, "Ephemeroptera", "Ephemeroptera", "Plecoptera",
              "Diptera", "Hemiptera", "Isopoda"),
    class = c("Bacillariophyceae", "Bacillariophyceae", NA, "Insecta", "Insecta", "Insecta",
              "Insecta", "Insecta", "Malacostraca"),
    stringsAsFactors = FALSE
  )
}

# The rule of the reanalysis of Amador et al. (2024), the one the hand-derived
# expectations of the tests were written for (independent of the defaults).
tiny_rule <- function(...) {
  matching_rule(mode = "climb", base_level = "family", climb_max = "class", min_partners = 1L, ...)
}

edge_keys <- function(edges) {
  sort(paste(edges$resource, edges$consumer, sep = ">"), method = "radix")
}

# A calibration object with known coefficients (no fitting), for sampling tests.
fake_calibration <- function(alpha, sigma = 1, se = 0.01, rule = tiny_rule(),
                             registries = "gray2015") {
  theta <- c(alpha, log_sigma = log(sigma))
  V <- diag(se^2, length(theta))
  dimnames(V) <- list(names(theta), names(theta))
  structure(
    list(
      coefficients = data.frame(evidence_class = names(alpha), estimate = unname(alpha),
                                std_error = se, lower = unname(alpha) - 2 * se,
                                upper = unname(alpha) + 2 * se,
                                p_average = stats::plogis(unname(alpha)),
                                p_new = p_new_system(unname(alpha), sigma),
                                webs = 3L, predicted = 100L, observed = 50L,
                                stringsAsFactors = FALSE),
      sigma = sigma, theta = theta,
      vcov = V, draws = NULL,
      counts = data.frame(), webs = data.frame(),
      settings = list(rule = rule, sources = c("a", "b", "c"), registries = registries,
                      n_records = 1L, quadrature = 20L),
      diagnostics = list(converged = TRUE, information_pd = TRUE, message = NA_character_,
                         classes_without_data = character(0))
    ),
    class = "bwebs_calibration"
  )
}

# Records of two webs: A (Baetis rhodani eats Cocconeis placentula; Perla
# bipunctata eats Baetis rhodani) and B (the same genera, plus Perla eating
# Cocconeis), for the validation tests.
two_web_records <- function() {
  rec <- data.frame(
    study_id = c("A", "A", "B", "B", "B"),
    resource = c("Cocconeis placentula", "Baetis rhodani", "Cocconeis pediculus",
                 "Baetis muticus", "Cocconeis pediculus"),
    consumer = c("Baetis rhodani", "Perla bipunctata", "Baetis muticus", "Perla marginata",
                 "Perla marginata"),
    res_genus = c("Cocconeis", "Baetis", "Cocconeis", "Baetis", "Cocconeis"), res_subfamily = NA,
    res_family = c("Cocconeidaceae", "Baetidae", "Cocconeidaceae", "Baetidae", "Cocconeidaceae"),
    res_order = c("Achnanthales", "Ephemeroptera", "Achnanthales", "Ephemeroptera", "Achnanthales"),
    res_class = c("Bacillariophyceae", "Insecta", "Bacillariophyceae", "Insecta", "Bacillariophyceae"),
    con_genus = c("Baetis", "Perla", "Baetis", "Perla", "Perla"), con_subfamily = NA,
    con_family = c("Baetidae", "Perlidae", "Baetidae", "Perlidae", "Perlidae"),
    con_order = c("Ephemeroptera", "Plecoptera", "Ephemeroptera", "Plecoptera", "Plecoptera"),
    con_class = "Insecta",
    res_category = c("producer", "invertebrate", "producer", "invertebrate", "producer")
  )
  rec
}
