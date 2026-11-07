test_that("shared_studies() finds studies of the same system", {
  x <- shared_studies(gray2015, c("13", "50", "46"), threshold = 0.2)
  expect_identical(x$web, c("13", "50"))
  expect_identical(x$study_id, c("50", "13"))
  expect_equal(x$share[1], x$n_shared[1] / x$n_links[1])
  expect_gt(x$share[1], 0.6)
  expect_identical(attr(x, "groups"), list(c("13", "50")))
  expect_error(shared_studies(gray2015, "zz"), "not in `registry`")
})

test_that("Tuesday Lake is study 12 of the Gray registry", {
  tl <- tuesday_lake_links
  name <- stats::setNames(tuesday_lake_taxa$taxon, tuesday_lake_taxa$node_id)
  links <- data.frame(web = tl$sample_id, resource = unname(name[tl$resource]),
                      consumer = unname(name[tl$consumer]))
  x <- shared_studies(gray2015, links, threshold = 0.1)
  expect_setequal(x$web, c("TL84", "TL86"))
  expect_true(all(x$study_id == "12"))
  expect_true(all(x$share > 0.6))
})

test_that("connected groups join chains of pairs", {
  expect_identical(connected_groups(c("a", "b", "d"), c("b", "c", "e")),
                   list(c("a", "b", "c"), c("d", "e")))
  expect_identical(connected_groups(character(), character()), list())
})
