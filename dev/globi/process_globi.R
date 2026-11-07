# Processing of the GloBI download (dev/globi/download_globi.R) for the
# comparison of sources: GloBI registry, GloBI studies copied from Tuesday Lake
# and studies that share links with each validation web of Gray (exclusions by
# web).
#
# Run from the root of the package:
#   Rscript dev/globi/process_globi.R <path>/globi_raw.rds [<output folder>]
# Output (in the output folder; by default, that of globi_raw.rds):
#   globi_registry.rds: list with the GloBI registry (`globi`), the registry
#     without the studies of Tuesday Lake (`globi_test`), the exclusions by web
#     (`exclusions`) and the data of the download.
# And in dev/globi/ (small, for version control):
#   globi_tuesday_lake.csv: GloBI studies removed for the final test;
#   globi_exclusions.csv: GloBI studies left out when each web is validated.
# The GloBI registry is not included in the package: the licences of its data
# depend on each source (see the citations).
pkgload::load_all(".", quiet = TRUE)
args <- commandArgs(trailingOnly = TRUE)
raw_path <- args[1]
out_dir <- if (length(args) > 1L) args[2] else dirname(raw_path)
raw <- readRDS(raw_path)
# attribute written by dev/globi/download_globi.R ("descarga", with the fields
# fecha, campos and tipos, in downloads made before 2026-11-07)
info <- attr(raw, "download")
if (is.null(info)) {
  info <- attr(raw, "descarga")
  old_names <- c(fecha = "date", campos = "fields", tipos = "types")
  hit <- names(info) %in% names(old_names)
  names(info)[hit] <- old_names[names(info)[hit]]
}
cat("Records downloaded:", nrow(raw), "\n")

# 1. Registry (without the copies of Gray and Brose in GloBI and without the
#    records of other taxa than the one queried)
globi <- globi_to_registry(raw)
rm(raw)
cat("GloBI registry:", nrow(globi), "records of", length(unique(globi$study_id)), "studies\n")
clash <- intersect(globi$study_id, gray2015$study_id)
if (length(clash)) stop("Studies with the same identifier in Gray and GloBI: ", toString(clash))

# 2. GloBI studies that contain Tuesday Lake: those holding at least 10% of
#    the links of TL84 or TL86 (with the names of the taxa) and those citing
#    its authors (Cohen, Jonsson and Carpenter; Brose et al. 2005 includes
#    them)
name <- stats::setNames(tuesday_lake_taxa$taxon, tuesday_lake_taxa$node_id)
tl <- data.frame(web = tuesday_lake_links$sample_id,
                 resource = unname(name[tuesday_lake_links$resource]),
                 consumer = unname(name[tuesday_lake_links$consumer]),
                 stringsAsFactors = FALSE)
sh_tl <- shared_studies(globi, tl, threshold = 0.1)
tl_pattern <- paste0("Tuesday|Jonsson et al\\.? \\(200[45]\\)|Cohen.{0,40}Jonsson|",
                     "Brose.{0,60}Jonsson|Cohen_2003")
cit_tl <- unique(globi$study_id[grepl(tl_pattern, globi$study_id, ignore.case = TRUE)])
drop_tl <- sort(union(sh_tl$study_id, cit_tl), method = "radix")
tl_table <- data.frame(study_id = drop_tl,
                       share_tl84 = vapply(drop_tl, function(s) {
                         v <- sh_tl$share[sh_tl$web == "TL84" & sh_tl$study_id == s]
                         if (length(v)) v else 0
                       }, numeric(1)),
                       share_tl86 = vapply(drop_tl, function(s) {
                         v <- sh_tl$share[sh_tl$web == "TL86" & sh_tl$study_id == s]
                         if (length(v)) v else 0
                       }, numeric(1)),
                       by_citation = drop_tl %in% cit_tl,
                       n_records = as.integer(table(globi$study_id)[drop_tl]),
                       stringsAsFactors = FALSE, row.names = NULL)
cat("GloBI studies with Tuesday Lake:", length(drop_tl), "(", sum(tl_table$n_records),
    "records)\n")
globi_test <- globi[!globi$study_id %in% drop_tl, , drop = FALSE]

# 3. Exclusions by web: GloBI studies holding at least 10% of the links of
#    each validation web (without Tuesday Lake, study 12)
val <- setdiff(gray2015_sources$study_id[gray2015_sources$validation], "12")
both <- combine_registries(gray2015[gray2015$study_id != "12", ], globi_test)
sh_val <- shared_studies(both, val, threshold = 0.05)
sh_val <- sh_val[sh_val$study_id %in% globi_test$study_id, , drop = FALSE]
exclusions <- sh_val[sh_val$share >= 0.1, , drop = FALSE]
rownames(exclusions) <- NULL
cat("Exclusions by web (GloBI studies):\n")
print(table(exclusions$web))
cat("Web-study pairs with threshold 0.05 / 0.1 / 0.2 / 0.5:",
    vapply(c(0.05, 0.1, 0.2, 0.5), function(th) sum(sh_val$share >= th), integer(1)), "\n")

utils::write.csv(tl_table, file.path("dev", "globi", "globi_tuesday_lake.csv"), row.names = FALSE)
utils::write.csv(exclusions, file.path("dev", "globi", "globi_exclusions.csv"), row.names = FALSE)
saveRDS(list(globi = globi, globi_test = globi_test, exclusions = exclusions[c("web", "study_id")],
             shared_validation = sh_val, tuesday_lake = tl_table, download = info),
        file.path(out_dir, "globi_registry.rds"))
cat("Saved in", file.path(out_dir, "globi_registry.rds"), "\n")
