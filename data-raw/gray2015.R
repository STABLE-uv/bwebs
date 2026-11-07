# Registry of trophic interactions of Gray et al. (2015) -> data/gray2015.rda
#
# Source: Gray C, Ma A, Perkins D, Hudson L, Figueroa D, Woodward G (2015).
#   Database of trophic interactions (version 1). Zenodo,
#   doi:10.5281/zenodo.13751, file trophic.links.2014-11-10.csv
#   (record checked on 2026-11-03).
# Licence: Creative Commons Zero v1.0 Universal (stated on the Zenodo record).
# Original file checksum (tools::md5sum): cb2d75df4291e2214bda9fc3ab743d64.
# The file is windows-1252 encoded, not UTF-8.
# Packages used only here: readr, pkgload, usethis.
#
# Run from the package root: Rscript data-raw/gray2015.R
# The CSV is not kept in the repository; give its path in BWEBS_GRAY_CSV or
# put it in data-raw/.

path <- Sys.getenv("BWEBS_GRAY_CSV", file.path("data-raw", "trophic.links.2014-11-10.csv"))
stopifnot(file.exists(path))
stopifnot(unname(tools::md5sum(path)) == "cb2d75df4291e2214bda9fc3ab743d64")
pkgload::load_all(".", quiet = TRUE)

raw <- readr::read_csv(
  path, na = c("", "NA"), show_col_types = FALSE,
  locale = readr::locale(encoding = "windows-1252"),
  col_types = readr::cols(.default = readr::col_character())
)
raw <- as.data.frame(raw, stringsAsFactors = FALSE)
for (col in names(raw)) raw[[col]] <- trimws(raw[[col]])

# Test rows of the original file (names starting with "zz")
test_row <- grepl("^zz", raw$resource) | grepl("^zz", raw$consumer)
raw <- raw[!test_row, , drop = FALSE]

reg <- data.frame(
  record_id = raw$link.id, study_id = raw$source.id,
  resource = raw$resource, consumer = raw$consumer,
  res_stage = raw$resource.lifestage, con_stage = raw$consumer.lifestage,
  res_genus = raw$res.genus, res_subfamily = raw$res.subfamily,
  res_family = raw$res.family, res_order = raw$res.order, res_class = raw$res.class,
  con_genus = raw$con.genus, con_subfamily = raw$con.subfamily,
  con_family = raw$con.family, con_order = raw$con.order, con_class = raw$con.class,
  res_category = raw$res.category, con_category = raw$con.category,
  evidence = raw$link.evidence, citation = raw$full.source,
  stringsAsFactors = FALSE
)

# Row-level corrections of the original analysis: names whose hierarchy was
# resolved wrongly by the Global Names Resolver ("Cladocera" as a lepidopteran
# genus; "Diatoma" and "Cyclops" as homonyms in Rhizophoraceae and
# Nassariidae), a missing class ("Oligochaeta") and a missing family
# ("Caenis"). The whole hierarchy of every row whose resource or consumer has
# that name is replaced.
name_fixes <- data.frame(
  name      = c("Cladocera", "Oligochaeta", "Diatoma", "Cyclops", "Caenis"),
  genus     = c(NA, NA, "Diatoma", "Cyclops", "Caenis"),
  subfamily = c(NA, NA, NA, NA, NA),
  family    = c(NA, NA, "Fragilariaceae", "Cyclopidae", "Caenidae"),
  order     = c("Diplostraca", NA, "Fragilariales", "Cyclopoida", "Ephemeroptera"),
  class     = c("Branchiopoda", "Clitellata", "Coscinodiscophyceae", "Maxillopoda", "Insecta"),
  stringsAsFactors = FALSE
)
for (side in c("resource", "consumer")) {
  pre <- if (side == "resource") "res_" else "con_"
  for (k in seq_len(nrow(name_fixes))) {
    hit <- !is.na(reg[[side]]) & tolower(reg[[side]]) == tolower(name_fixes$name[k])
    for (r in c("genus", "subfamily", "family", "order", "class")) {
      reg[[paste0(pre, r)]][hit] <- name_fixes[[r]][k]
    }
  }
}

# Misspellings and synonyms of higher ranks (default_taxonomy_fixes()) and
# life stages are harmonised by as_registry().
gray2015 <- as_registry(reg, registry = "gray2015")

# One row per study, to choose the validation webs of calibrate_links().
by_study <- split(gray2015, gray2015$study_id)
gray2015_sources <- do.call(rbind, lapply(names(by_study), function(id) {
  d <- by_study[[id]]
  pairs <- unique(d[c("resource", "consumer")])
  data.frame(
    study_id = id,
    citation = d$citation[1],
    n_records = nrow(d),
    n_links = nrow(pairs),
    n_consumers = length(unique(d$consumer)),
    n_resources = length(unique(d$resource)),
    n_taxa = length(unique(c(d$consumer, d$resource))),
    n_repeated_records = nrow(d) - nrow(pairs),
    stringsAsFactors = FALSE
  )
}))
gray2015_sources <- gray2015_sources[order(as.integer(gray2015_sources$study_id)), ]
# Provisional validation set (whole communities with at least 8 consumers
# whose diets were recorded; single-taxon diet studies and compilations of
# one group, such as water mites, are excluded). To be confirmed by the user.
validation_ids <- c("27", "13", "50", "47", "12", "59", "63", "11", "15", "58", "46",
                    "49", "23")
gray2015_sources$validation <- gray2015_sources$study_id %in% validation_ids
rownames(gray2015_sources) <- NULL

usethis::use_data(gray2015, gray2015_sources, overwrite = TRUE, compress = "xz")
