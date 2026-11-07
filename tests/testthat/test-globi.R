test_that("paths keep empty elements, also at the end", {
  expect_identical(split_path("a | b |"), c("a", "b", ""))
  expect_identical(rank_from_path("Animalia | Arthropoda | ", "kingdom | phylum | ", "phylum"),
                   "Arthropoda")
  expect_identical(rank_from_path("A | B", "kingdom", "kingdom"), NA_character_)
})

test_that("placeholders of unresolved ranks become NA", {
  ranks <- "kingdom | class | order | family"
  path <- c("Animalia | Not assigned | Ploima | Brachionidae",
            "Animalia | Branchiopoda | Diplostraca | unclassified Cladoceromorpha",
            "Animalia | Teleostei | Perciformes | Ovalentaria incertae sedis")
  expect_identical(rank_from_path(path, rep(ranks, 3), "class"),
                   c(NA, "Branchiopoda", "Teleostei"))
  expect_identical(rank_from_path(path, rep(ranks, 3), "family"),
                   c("Brachionidae", NA, NA))
  expect_identical(is_placeholder(c("unknown_eukaryote", "Unionidae", NA, "unidentified_bacteria")),
                   c(TRUE, FALSE, FALSE, TRUE))
})

test_that("ranks_from_paths() gives the ranks of each path once per distinct path", {
  ranks <- "kingdom | class | order | family | genus"
  path <- c("Animalia | Insecta | Ephemeroptera | Baetidae | Baetis",
            "Animalia | Malacostraca | Isopoda | Asellidae | Asellus")[c(1, 2, 1, 1, 2)]
  out <- ranks_from_paths(path, rep(ranks, 5))
  expect_named(out, c("genus", "subfamily", "family", "order", "class"))
  for (r in names(out)) {
    expect_identical(out[[r]], rank_from_path(path, rep(ranks, 5), r))
  }
  expect_identical(out$family, c("Baetidae", "Asellidae", "Baetidae", "Baetidae", "Asellidae"))
  expect_true(all(is.na(out$subfamily)))
})

test_that("globi_to_registry() keeps trophic records of the taxon queried", {
  path_a <- "Animalia | Arthropoda | Malacostraca | Isopoda | Asellidae | Asellus | Asellus aquaticus"
  ranks <- "kingdom | phylum | class | order | family | genus | species"
  raw <- data.frame(
    source_taxon_name = c("Asellus aquaticus", "Hadra sp.", "Asellus aquaticus", "Asellus aquaticus"),
    source_taxon_path = c(path_a, "Animalia | Mollusca | Hadra", path_a, path_a),
    source_taxon_path_ranks = c(ranks, "kingdom | phylum | genus", ranks, ranks),
    interaction_type = c("eats", "eats", "parasiteOf", "eats"),
    target_taxon_name = "Cocconeis placentula",
    target_taxon_path = "Chromista | Ochrophyta | Bacillariophyceae | Achnanthales | Cocconeidaceae | Cocconeis | Cocconeis placentula",
    target_taxon_path_ranks = ranks,
    study_citation = c("Study 1", "Study 2", "Study 3", "Gray copy"),
    study_source_id = c(NA, NA, NA, "globalbioticinteractions/gray2015"),
    query = "Asellus", query_type = "eats"
  )
  reg <- suppressMessages(globi_to_registry(raw))
  expect_identical(nrow(reg), 1L)
  expect_identical(c(reg$consumer, reg$resource), c("Asellus aquaticus", "Cocconeis placentula"))
  expect_identical(c(reg$con_family, reg$res_class), c("Asellidae", "Bacillariophyceae"))
  expect_identical(c(reg$registry, reg$study_id), c("globi", "Study 1"))
  expect_message(expect_message(globi_to_registry(raw), "compilations"), "other taxa")
})

test_that("fetch_globi() validates its input before going online", {
  expect_error(fetch_globi(character(0)), class = "bwebs_error_input")
  expect_error(fetch_globi("Asellus", cache_dir = file.path(tempdir(), "no_such_dir")), "does not exist")
})
