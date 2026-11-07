test_that("fetch_globi() downloads records from GloBI", {
  skip_if_not_installed("curl")
  skip_if_offline()
  raw <- fetch_globi("Asellus aquaticus", page_size = 20, max_pages = 1, quiet = TRUE)
  expect_true(all(c("source_taxon_name", "target_taxon_path", "query") %in% names(raw)))
  reg <- suppressMessages(globi_to_registry(raw))
  expect_true(all(reg$registry == "globi"))
})
