test_that("matching_rule() validates its settings", {
  expect_s3_class(matching_rule(), "bwebs_rule")
  expect_error(matching_rule(mode = "climbing"), class = "bwebs_error_input")
  expect_error(matching_rule(base_level = "tribe"), "must be one of")
  expect_error(matching_rule(min_partners = 0), class = "bwebs_error_input")
  expect_error(matching_rule(min_partners = 1.5), "whole number")
  expect_error(matching_rule(downward_ranks = "class"), "downward")
  expect_error(matching_rule(exceptions = data.frame(role = "x", rank = "order", taxon = "A",
                                                     base = "order", rule = "R3", reason = "r")),
               "role")
  expect_identical(nrow(matching_rule(generic_groups = NULL)$generic_groups), 0L)
})

test_that("rules that differ only in manual-diet nodes are compatible", {
  a <- matching_rule()
  b <- matching_rule(manual_diet = "x")
  expect_false(same_rule(a, b))
  expect_true(rules_compatible(a, b))
  expect_false(rules_compatible(a, matching_rule(min_partners = 2)))
})

test_that("every registry is used at every level unless a primary one is named", {
  # Web A (registry gray2015) is rebuilt from study C of the registry "other",
  # in which a heptageniid mayfly eats Cocconeis: Baetis (Baetidae) only
  # reaches that record by climbing to its order, Ephemeroptera
  rec <- two_web_records()
  rec_c <- data.frame(
    study_id = "C", resource = "Cocconeis pediculus", consumer = "Ecdyonurus venosus",
    res_genus = "Cocconeis", res_subfamily = NA, res_family = "Cocconeidaceae",
    res_order = "Achnanthales", res_class = "Bacillariophyceae",
    con_genus = "Ecdyonurus", con_subfamily = NA, con_family = "Heptageniidae",
    con_order = "Ephemeroptera", con_class = "Insecta", res_category = "producer"
  )
  a <- as_registry(rec[rec$study_id == "A", ], registry = "gray2015")
  mixed <- rbind(a, as_registry(rec_c, registry = "other"))
  expect_null(matching_rule()$primary_registry)
  expect_error(matching_rule(primary_registry = 1), class = "bwebs_error_input")
  all_levels <- validate_web(mixed, "A", tiny_rule())
  expect_true("order" %in% all_levels$counts$evidence_class)
  # with gray2015 as primary, "other" is only used up to family
  primary <- validate_web(mixed, "A", tiny_rule(primary_registry = "gray2015"))
  expect_identical(nrow(primary$counts), 0L)
  expect_identical(validate_web(mixed, "A", tiny_rule(primary_registry = "gray2015",
                                                      supplementary_max_level = "class")),
                   all_levels)
  # with a single registry, naming it primary changes nothing
  reg <- as_registry(rec, registry = "gray2015")
  expect_identical(validate_web(reg, "A", tiny_rule(primary_registry = "gray2015")),
                   validate_web(reg, "A", tiny_rule()))
})

test_that("print method is stable", {
  expect_snapshot(print(matching_rule()))
  expect_snapshot(print(matching_rule(primary_registry = "gray2015")))
})
