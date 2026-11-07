test_that("as_registry() harmonises names, stages and keeps the record order", {
  rec <- data.frame(
    study_id = c("a", "a", "b"), resource = c(" Corixa", "Diatom", NA),
    consumer = c("Notonecta", "Baetis", "X"),
    res_genus = c("Corixa", NA, NA), res_subfamily = NA, res_family = c("Corixidae", NA, NA),
    res_order = c("Heteroptera", NA, NA), res_class = c("Insecta", NA, NA),
    con_genus = c("Notonecta", "Baetis", NA), con_subfamily = NA,
    con_family = c("Notonectidae", "Baetidae", NA), con_order = c("Hemiptera", "Ephemeroptera", NA),
    con_class = c("Insecta", "Insecta", NA),
    con_stage = c("larvae", "adults; juveniles", "adult"), stringsAsFactors = FALSE
  )
  expect_message(reg <- as_registry(rec, registry = "test"), "dropped")
  expect_identical(nrow(reg), 2L)
  expect_identical(reg$resource, c("Corixa", "Diatom"))          # trimmed
  expect_identical(reg$res_order[1], "Hemiptera")                 # synonym harmonised
  expect_identical(reg$con_stage, c("larva", NA))                 # mixed stage -> NA
  expect_identical(reg$registry, c("test", "test"))
  expect_identical(names(reg), registry_columns)
})

test_that("normalise_stage() follows the original rules", {
  expect_identical(normalise_stage(c("larvae", "nymph", "adult", "adults", "larvae and adults",
                                     "adults; juveniles", "0+", NA, "pupae")),
                   c("larva", "larva", "adult", "adult", NA, NA, NA, NA, NA))
})

test_that("as_registry() reports missing columns and studies", {
  expect_error(as_registry(data.frame(resource = "a"), registry = "x"), "not found")
  rec <- tiny_registry()
  rec$study_id[1] <- NA
  expect_error(as_registry(rec, registry = "x"), "study")
  expect_error(as_registry(tiny_registry(), registry = ""), class = "bwebs_error_input")
})

test_that("combine_registries() keeps studies of different registries apart", {
  a <- tiny_registry()
  b <- a
  b$registry <- "other"
  ab <- combine_registries(a, b)
  expect_identical(nrow(ab), 2L * nrow(a))
  expect_identical(length(unique(ab$study_id)), 2L * length(unique(a$study_id)))
  expect_error(combine_registries(a, a), "different names")
})

test_that("the bundled registry has the expected size and checksum-derived content", {
  expect_identical(nrow(gray2015), 20822L)
  expect_identical(length(unique(gray2015$study_id)), 50L)
  expect_identical(names(gray2015), registry_columns)
  expect_identical(sum(gray2015_sources$validation), 13L)
  expect_identical(sum(gray2015_sources$n_records), nrow(gray2015))
})
