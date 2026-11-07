web_from_links <- function(nodes, links) {
  a <- matrix(0L, length(nodes), length(nodes), dimnames = list(nodes, nodes))
  for (l in links) a[l[1], l[2]] <- 1L
  a
}

test_that("a three-level chain gives the indicators computed by hand", {
  a <- web_from_links(c("alga", "grazer", "predator"),
                      list(c("alga", "grazer"), c("grazer", "predator")))
  w <- web_indicators(a)
  expect_identical(c(w$S_taxa, w$n_isolated, w$S, w$L), c(3L, 0L, 3L, 2L))
  expect_equal(c(w$LS, w$C), c(2 / 3, 2 / 9))
  expect_equal(c(w$prop_basal, w$prop_intermediate, w$prop_top), c(1, 1, 1) / 3)
  # generality: prey per consumer (1, 1); sd of (0, 1, 1) / (2/3) = sqrt(0.75)
  expect_equal(c(w$gen_mean, w$gen_sd, w$vul_mean, w$vul_sd), c(1, sqrt(0.75), 1, sqrt(0.75)))
  expect_equal(c(w$TL_mean, w$TL_max, w$TL_short_mean), c(2, 3, 2))
  expect_identical(w$TL_status, "solve")
  expect_equal(c(w$prop_omnivores, w$omnivory_mean), c(0, 0))
  expect_equal(c(w$MxSim_taxa, w$MxSim_trophic_species), c(0, 0))
  expect_identical(w$status, "ok")
})

test_that("an omnivore feeding on two trophic levels is detected", {
  a <- web_from_links(c("alga", "grazer", "predator"),
                      list(c("alga", "grazer"), c("alga", "predator"), c("grazer", "predator")))
  w <- web_indicators(a)
  # TL: alga 1, grazer 2, predator 1 + (1 + 2) / 2 = 2.5; shortest chain 1, 2, 2
  expect_equal(c(w$TL_mean, w$TL_max, w$TL_short_mean), c(5.5 / 3, 2.5, 5 / 3))
  # omnivory of the predator: sd(c(1, 2)); the grazer has one prey
  expect_equal(c(w$prop_omnivores, w$omnivory_mean), c(0.5, sd(c(1, 2)) / 2))
  # each species shares one of three partners with another: similarity 1/3
  expect_equal(w$MxSim_taxa, 1 / 3)
})

test_that("identical taxa are aggregated into trophic species and isolated nodes counted", {
  a <- web_from_links(c("alga", "g1", "g2", "pred", "lonely"),
                      list(c("alga", "g1"), c("alga", "g2"), c("g1", "pred"), c("g2", "pred"),
                           c("lonely", "lonely")))
  w <- web_indicators(a)
  expect_identical(c(w$S_taxa, w$n_isolated, w$S, w$L), c(5L, 1L, 3L, 2L))
  expect_equal(c(w$MxSim_taxa, w$MxSim_trophic_species), c(0.5, 0))
})

test_that("R50 of a chain with a permanent resource is 0.5", {
  a <- web_from_links(c("alga", "grazer", "predator"),
                      list(c("alga", "grazer"), c("grazer", "predator")))
  w <- web_indicators(a, permanent = "alga", robustness = TRUE, robustness_reps = 20, seed = 1)
  expect_equal(c(w$R50_random, w$R50_connected, w$R50_random_living, w$R50_connected_living),
               rep(0.5, 4))
})

test_that("a loop without basal input has unresolved trophic levels", {
  a <- web_from_links(c("a", "b"), list(c("a", "b"), c("b", "a")))
  w <- web_indicators(a)
  expect_identical(w$TL_status, "unresolved")
  expect_true(is.na(w$TL_mean))
  expect_true(is.na(w$TL_short_mean))
})

test_that("a web without links is empty, not zero", {
  a <- web_from_links(c("a", "b"), list())
  w <- web_indicators(a)
  expect_identical(w$status, "empty_web")
  expect_identical(c(w$S_taxa, w$n_isolated), c(2L, 2L))
  expect_true(is.na(w$C))
})

test_that("indicators do not depend on the order of the nodes", {
  withr::local_seed(5)
  n <- 12
  a <- matrix(stats::rbinom(n * n, 1, 0.2), n, n, dimnames = list(letters[1:n], letters[1:n]))
  a[, 1:3] <- 0L   # three basal species
  o <- sample.int(n)
  x <- web_indicators(a)
  y <- web_indicators(a[o, o])
  num <- vapply(x, is.numeric, logical(1))
  expect_equal(unlist(x[num]), unlist(y[num]))
})

test_that("shortest-chain trophic levels agree with igraph distances", {
  skip_if_not_installed("igraph")
  withr::local_seed(7)
  n <- 15
  B <- matrix(stats::rbinom(n * n, 1, 0.15), n, n, dimnames = list(letters[1:n], letters[1:n]))
  diag(B) <- 0L
  B[, 1:3] <- 0L
  tl <- trophic_level_short(B)
  g <- igraph::graph_from_adjacency_matrix(B, mode = "directed")
  d <- igraph::distances(g, v = which(colSums(B) == 0), mode = "out")
  ref <- 1 + apply(d, 2, min)
  ref[is.infinite(ref)] <- NA
  expect_equal(unname(tl), unname(ref))
})

test_that("modularity needs igraph and returns a value", {
  skip_if_not_installed("igraph")
  a <- web_from_links(c("a1", "a2", "c1", "c2", "b1", "b2", "d1", "d2"),
                      list(c("a1", "c1"), c("a2", "c1"), c("a1", "c2"), c("a2", "c2"),
                           c("b1", "d1"), c("b2", "d1"), c("b1", "d2"), c("b2", "d2")))
  w <- web_indicators(a, modularity = TRUE, seed = 1)
  expect_identical(w$n_modules, 2L)
  expect_gt(w$modularity, 0.4)
})

test_that("modularity explains how to install igraph when it is missing", {
  local_mocked_bindings(is_installed = function(pkg) FALSE)
  a <- web_from_links(c("a", "b"), list(c("a", "b")))
  expect_error(web_indicators(a, modularity = TRUE), class = "bwebs_error_missing_package")
})

test_that("web_indicators() rejects invalid matrices", {
  expect_error(web_indicators(matrix(1, 2, 3)), class = "bwebs_error_input")
  expect_error(web_indicators(matrix(c(0, NA, 0, 0), 2, 2, dimnames = list(c("a", "b"), c("a", "b")))),
               "NA")
  expect_error(web_indicators(matrix(0, 2, 2)), "names")
})
