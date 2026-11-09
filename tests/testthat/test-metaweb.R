test_that("build_metaweb() infers the links derived by hand for the tiny registry", {
  # Expected links with the default rule (start at family, climb to class,
  # one partner, predators do not climb as resources):
  #   cocc -> baet  Baetis rhodani eats Cocconeis placentula (same genera)
  #   baet -> perl  Perla bipunctata eats Baetis rhodani (same genera)
  #   navi -> ecdy  Ecdyonurus (Heptageniidae) eats Navicula (same family / genus)
  #   ecdy -> perl  ecdy has no predators in its family and is not a predator,
  #                 so as a resource it climbs to Ephemeroptera (Perla eats Baetis)
  #   cpom -> asel, cpom -> simu  by name (CPOM) and genus
  #   simu -> gerr  simu climbs as resource to Diptera: "Gerris eats Diptera"
  # perl and gerr are predators: they do not climb as resources (no predators).
  nodes <- prepare_nodes(tiny_taxa(), rule = tiny_rule())
  mw <- build_metaweb(nodes, tiny_registry())
  expect_s3_class(mw, "bwebs_metaweb")
  expect_identical(edge_keys(mw$edges), sort(c(
    "cocc>baet", "baet>perl", "navi>ecdy", "ecdy>perl", "cpom>asel", "cpom>simu", "simu>gerr"
  ), method = "radix"))
  cls <- stats::setNames(mw$edges$evidence_class, paste(mw$edges$resource, mw$edges$consumer, sep = ">"))
  expect_identical(cls[["cocc>baet"]], "genus")
  expect_identical(cls[["baet>perl"]], "genus")
  expect_identical(cls[["navi>ecdy"]], "family")
  expect_identical(cls[["ecdy>perl"]], "order")
  expect_identical(cls[["cpom>asel"]], "genus")
  expect_identical(cls[["simu>gerr"]], "order")
  expect_true(all(mw$edges$support == "direct"))
  lev <- stats::setNames(mw$nodes$resource_level_used, mw$nodes$node_id)
  expect_identical(lev[["ecdy"]], "order")
  expect_identical(lev[["perl"]], "no_match")
  expect_true(mw$nodes$predator[mw$nodes$node_id == "perl"])
  expect_identical(sum(mw$adjacency), 7L)
  expect_identical(mw$adjacency["simu", "gerr"], 1L)
})

test_that("each link keeps the amount of evidence that supports it", {
  nodes <- prepare_nodes(example_taxa)
  mw <- build_metaweb(nodes, gray2015)
  e <- mw$edges[mw$edges$support != "manual", ]
  expect_true(all(e$n_sources <= e$n_records))
  expect_true(all(e$n_consumer_taxa >= 1L & e$n_consumer_taxa <= e$n_consumer_known))
  expect_true(all(e$n_resource_taxa >= 1L & e$n_resource_taxa <= e$n_resource_known))
  expect_true(all(e$n_consumer_taxa <= e$n_records & e$n_resource_taxa <= e$n_records))
  # the relatives with records are a property of each node
  known <- tapply(e$n_consumer_known, e$consumer, function(v) length(unique(v)))
  expect_true(all(known == 1L))
  tiny <- build_metaweb(prepare_nodes(tiny_taxa(), rule = tiny_rule()), tiny_registry(),
                        manual_links = data.frame(resource = "cocc", consumer = "asel",
                                                  action = "add"))
  m <- tiny$edges[tiny$edges$support == "manual", ]
  expect_true(all(is.na(unlist(m[c("n_consumer_taxa", "n_resource_taxa", "n_consumer_known",
                                   "n_resource_known")]))))
})

test_that("downward inheritance extends coarse records to the nodes of that taxon", {
  # Without climbing, simu has no records as a resource; the record
  # "Gerris eats Diptera" (resolved only to order) reaches it downwards.
  nodes <- prepare_nodes(tiny_taxa(), rule = matching_rule(mode = "fixed", base_level = "family"))
  mw <- build_metaweb(nodes, tiny_registry())
  e <- mw$edges[mw$edges$resource == "simu" & mw$edges$consumer == "gerr", ]
  expect_identical(nrow(e), 1L)
  expect_identical(e$support, "downward")
  expect_identical(e$evidence_class, "order")
  # and without downward inheritance the link disappears
  nodes2 <- prepare_nodes(tiny_taxa(), rule = matching_rule(mode = "fixed", base_level = "family", downward_ranks = character(0)))
  mw2 <- build_metaweb(nodes2, tiny_registry())
  expect_false("simu>gerr" %in% edge_keys(mw2$edges))
})

test_that("manual links are added and removed, and basal nodes never consume", {
  nodes <- prepare_nodes(tiny_taxa(), rule = tiny_rule())
  ml <- data.frame(resource = c("cocc", "cpom", "baet"), consumer = c("asel", "asel", "cocc"),
                   action = c("add", "remove", "add"), source = "test")
  mw <- build_metaweb(nodes, tiny_registry(), manual_links = ml)
  keys <- edge_keys(mw$edges)
  expect_true("cocc>asel" %in% keys)
  expect_false("cpom>asel" %in% keys)
  expect_false("baet>cocc" %in% keys)   # a producer does not consume
  expect_identical(mw$edges$evidence_class[mw$edges$resource == "cocc" & mw$edges$consumer == "asel"],
                   "manual")
})

test_that("a surrogate replaces the ladder and links sharing no rank are extrapolated", {
  tx <- tiny_taxa()
  tx <- rbind(tx, data.frame(
    node_id = "duga", taxon = "Dugastella valentina", rank = "species", role = "consumer",
    node_type = "taxon", genus = "Dugastella", subfamily = NA, family = "Atyidae",
    order = "Decapoda", class = "Malacostraca"
  ))
  sur <- data.frame(node_id = "duga", role = "consumer", surrogate = "Asellus")
  nodes <- prepare_nodes(tx, surrogates = sur, rule = tiny_rule())
  expect_identical(nodes$consumer_rule[nodes$node_id == "duga"], "R5")
  mw <- build_metaweb(nodes, tiny_registry())
  e <- mw$edges[mw$edges$consumer == "duga", ]
  expect_identical(e$resource, "cpom")
  # CPOM by name; Dugastella and Asellus share only the class Malacostraca
  expect_identical(e$evidence_class, "class")
  expect_false(e$extrapolated)
})

test_that("build_metaweb() rejects invalid input with classed errors", {
  nodes <- prepare_nodes(tiny_taxa(), rule = tiny_rule())
  expect_error(build_metaweb(tiny_taxa(), tiny_registry()), class = "bwebs_error_input")
  expect_error(build_metaweb(nodes, data.frame(a = 1)), class = "bwebs_error_input")
  ml <- data.frame(resource = "zzz", consumer = "asel", action = "add")
  expect_error(build_metaweb(nodes, tiny_registry(), manual_links = ml), "not in `nodes`")
})

test_that("discarding records that match no node does not change the metaweb", {
  nodes <- prepare_nodes(example_taxa)
  keep <- relevant_rows(gray2015, nodes, attr(nodes, "rule"))
  both <- relevant_rows(gray2015, nodes, attr(nodes, "rule"), sides = "both")
  expect_lt(sum(keep), nrow(gray2015))
  expect_true(all(keep[both]))
  filtered <- build_metaweb(nodes, gray2015)
  local_mocked_bindings(relevant_rows = function(registry, nodes, rule, sides) rep(TRUE, nrow(registry)))
  full <- build_metaweb(nodes, gray2015)
  expect_identical(filtered$edges, full$edges)
  expect_identical(filtered$nodes, full$nodes)
  expect_identical(filtered$adjacency, full$adjacency)
  expect_identical(filtered$settings, full$settings)
})

test_that("the metaweb does not depend on the order of the registry records", {
  nodes <- prepare_nodes(tiny_taxa(), rule = tiny_rule())
  reg <- tiny_registry()
  a <- build_metaweb(nodes, reg)
  b <- build_metaweb(nodes, reg[rev(seq_len(nrow(reg))), ])
  expect_identical(a$edges[c("resource", "consumer", "evidence_class", "n_records", "n_sources")],
                   b$edges[c("resource", "consumer", "evidence_class", "n_records", "n_sources")])
})

test_that("print method is stable", {
  mw <- build_metaweb(prepare_nodes(tiny_taxa(), rule = tiny_rule()), tiny_registry())
  expect_snapshot(print(mw))
})
