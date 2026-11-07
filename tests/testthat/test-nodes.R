test_that("own_match_level() maps ranks to matching levels (rule R1)", {
  expect_identical(own_match_level(c("species", "genus", "tribe", "subfamily", "family",
                                     "suborder", "order", "class", "phylum", "informal")),
                   c("exact", "genus", "subfamily", "subfamily", "family", "order", "order",
                     "class", "exact", "exact"))
})

test_that("prepare_nodes() applies rules R2, R3, R5, NT and basal", {
  tx <- tiny_taxa()
  tx <- rbind(tx, data.frame(node_id = "ostr", taxon = "Ostracoda", rank = "class",
                             role = "consumer", node_type = "taxon", genus = NA, subfamily = NA,
                             family = NA, order = NA, class = "Ostracoda"))
  sur <- data.frame(node_id = "asel", role = "consumer", surrogate = "Gammarus")
  nodes <- prepare_nodes(tx, surrogates = sur)
  get <- function(id, col) nodes[[col]][nodes$node_id == id]
  # R2: start at family, climb up to class
  expect_identical(c(get("baet", "consumer_method"), get("baet", "consumer_base"),
                     get("baet", "consumer_max"), get("baet", "consumer_rule")),
                   c("exact", "family", "class", "R2"))
  # R3: identified coarser than family
  expect_identical(c(get("ostr", "consumer_base"), get("ostr", "consumer_rule")), c("class", "R3"))
  # R5: surrogate replaces the ladder in that role only
  expect_identical(get("asel", "consumer_rule"), "R5")
  expect_identical(get("asel", "resource_rule"), "R2")
  # NT and basal
  expect_identical(get("cpom", "resource_rule"), "NT")
  expect_identical(get("cpom", "consumer_method"), "none")
  expect_identical(get("cocc", "consumer_rule"), "basal")
  expect_s3_class(nodes, "bwebs_nodes")
  expect_identical(nodes$node_id, tx$node_id)
})

test_that("exceptions and manual-diet nodes change the ladder", {
  rule <- matching_rule(manual_diet = "simu",
                        exceptions = data.frame(role = "consumer", rank = "family",
                                                taxon = "Baetidae", base = "order", rule = "R3",
                                                reason = "test"))
  nodes <- prepare_nodes(tiny_taxa(), rule = rule)
  expect_identical(nodes$consumer_base[nodes$node_id == "baet"], "order")
  expect_identical(nodes$resource_base[nodes$node_id == "baet"], "family")
  expect_identical(nodes$consumer_rule[nodes$node_id == "simu"], "manual")
})

test_that("prepare_nodes() rejects invalid tables with classed errors", {
  tx <- tiny_taxa()
  expect_error(prepare_nodes(tx[0, ]), class = "bwebs_error_input")
  bad <- tx; bad$node_id[2] <- bad$node_id[1]
  expect_error(prepare_nodes(bad), "one row per combination")
  bad <- tx; bad$role[1] <- "plant"
  expect_error(prepare_nodes(bad), "role")
  bad <- tx; bad$family <- NULL
  expect_error(prepare_nodes(bad), "not found")
  bad <- tx; bad$taxon[3] <- NA
  expect_error(prepare_nodes(bad), "missing values")
  sur <- data.frame(node_id = "nope", role = "consumer", surrogate = "Gammarus")
  expect_error(prepare_nodes(tx, surrogates = sur), "not in `taxa`")
  expect_error(prepare_nodes(tx, rule = list()), "matching_rule")
})

test_that("prepare_nodes() harmonises the classification like the registry", {
  tx <- rbind(tiny_taxa(), data.frame(
    node_id = "salm", taxon = "Salmo trutta", rank = "species", role = "consumer",
    node_type = "taxon", genus = "Salmo", subfamily = NA, family = "Salmonidae",
    order = "Salmoniformes", class = "Teleostei"
  ))
  tx$order[tx$node_id == "gerr"] <- "Heteroptera"
  nodes <- prepare_nodes(tx)
  expect_identical(nodes$class[nodes$node_id == "salm"], "Actinopterygii")
  expect_identical(nodes$order[nodes$node_id == "gerr"], "Hemiptera")
  expect_identical(nodes$taxon, tx$taxon)                        # names are not changed
  kept <- prepare_nodes(tx, taxonomy_fixes = character(0))
  expect_identical(kept$class[kept$node_id == "salm"], "Teleostei")
  expect_error(prepare_nodes(tx, taxonomy_fixes = "Hemiptera"), "named character")
  # the fixes are idempotent and the bundled Gray et al. (2015) registry is
  # already harmonised with them
  fx <- default_taxonomy_fixes()
  expect_false(any(fx %in% names(fx)))
  in_gray <- unique(unlist(gray2015[c(paste0("res_", rank_cols), paste0("con_", rank_cols),
                                      "resource_std", "consumer_std")]))
  expect_false(any(names(fx) %in% in_gray))
})

test_that("factors and tibbles are accepted as text columns", {
  tx <- tiny_taxa()
  tx$role <- factor(tx$role)
  expect_identical(prepare_nodes(tx)$role, as.character(tiny_taxa()$role))
})

test_that("print method is stable", {
  expect_snapshot(print(prepare_nodes(tiny_taxa(), rule = tiny_rule())))
})
