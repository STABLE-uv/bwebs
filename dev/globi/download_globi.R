# =============================================================================
# GloBI download for the comparison of methods of bwebs
# =============================================================================
# Run it from the root of the package on a computer with internet access (it
# only uses base R):
#   Rscript dev/globi/download_globi.R [<output folder>]
# or open it in RStudio, with the working directory at the root of the
# package, and press Source. The output folder (by default
# ../data/GloBI/bwebs, outside the repository) receives the cache of the
# queries, globi_raw.rds and query_summary.csv.
#
# - For each query of globi_queries.csv (consumers of the validation webs of
#   Gray, of Tuesday Lake and of the azoxystrobin experiment: name, genus and
#   family), it downloads the trophic interactions in which the taxon is the
#   consumer, record by record and with their citation
#   (includeObservations=true), so that the studies of Tuesday Lake and of
#   the validation webs can be removed.
# - Resumable cache: each query is saved in its own CSV; if the run stops, run
#   it again and it only downloads what is missing. Failed queries are not
#   saved and are tried again in the next run.
# - At the end it writes globi_raw.rds (all the records) and
#   query_summary.csv, to be processed with dev/globi/process_globi.R.
# - It takes one to three hours (about 500 queries).
# - It does not use the cache of the original analysis (records without
#   citations).
# =============================================================================

# ---- Paths and parameters (adjust them if needed) ----------------------------
args <- commandArgs(trailingOnly = TRUE)
queries_csv <- file.path("dev", "globi", "globi_queries.csv")
output      <- if (length(args)) args[1] else file.path("..", "data", "GloBI", "bwebs")
if (!file.exists(queries_csv)) {
  stop("Run the script from the root of the package (", queries_csv, " not found).",
       call. = FALSE)
}
api         <- "https://api.globalbioticinteractions.org/interaction"
page_size   <- 1000L   # records per page
max_pages   <- 100L    # limit per query (a warning is given when it is reached)
pause       <- 0.5     # seconds between requests
fields <- c("source_taxon_name", "source_taxon_path", "source_taxon_path_ranks",
            "interaction_type", "target_taxon_name", "target_taxon_path",
            "target_taxon_path_ranks", "study_citation", "study_source_citation",
            "study_source_id", "study_doi")
minimal_fields <- fields[1:8]   # those that worked in the original analysis

# ---- Functions -----------------------------------------------------------------
options(timeout = 300)
cache_dir <- file.path(output, "cache")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(output, "download_log.txt")
note <- function(...) {
  line <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "  ", paste0(...))
  cat(line, "\n")
  cat(line, "\n", file = log_file, append = TRUE)
}

page_url <- function(taxon, type, offset, fields, api) {
  paste0(api, "?sourceTaxon=", utils::URLencode(taxon, reserved = TRUE),
         "&interactionType=", type, "&includeObservations=true&type=csv",
         "&limit=", page_size, "&offset=", offset,
         "&fields=", paste(fields, collapse = ","))
}

# One page; NULL if it fails three times or the answer is not the expected table.
read_page <- function(u) {
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp), add = TRUE)
  for (attempt in 1:3) {
    ok <- tryCatch(withCallingHandlers({
      utils::download.file(u, tmp, mode = "wb", quiet = TRUE)
      TRUE
    }, warning = function(w) invokeRestart("muffleWarning")), error = function(e) FALSE)
    if (isTRUE(ok) && file.exists(tmp)) {
      d <- tryCatch(utils::read.csv(tmp, colClasses = "character", na.strings = c("", "NA"),
                                    encoding = "UTF-8", check.names = FALSE),
                    error = function(e) {
                      if (grepl("no lines available", conditionMessage(e))) data.frame() else NULL
                    })
      if (!is.null(d) && (nrow(d) == 0L || "source_taxon_name" %in% names(d))) return(d)
    }
    Sys.sleep(5 * attempt)
  }
  NULL
}

cache_file <- function(taxon, type) {
  file.path(cache_dir, paste0(gsub("[^A-Za-z0-9]+", "_", taxon), "_", type, ".csv"))
}

# Marker of a truncated query (".truncada" in caches made before 2026-11-07).
is_truncated <- function(f) file.exists(paste0(f, ".truncated")) || file.exists(paste0(f, ".truncada"))

# Downloads a whole query; returns list(data, truncated) or NULL if it fails.
download_query <- function(taxon, type, fields) {
  pages <- list()
  truncated <- FALSE
  for (p in seq_len(max_pages)) {
    d <- read_page(page_url(taxon, type, (p - 1L) * page_size, fields, api))
    if (is.null(d)) return(NULL)
    pages[[p]] <- d
    if (pause > 0) Sys.sleep(pause)
    if (nrow(d) < page_size) break
    if (p == max_pages) truncated <- TRUE
  }
  pages <- pages[vapply(pages, nrow, integer(1)) > 0]
  data <- if (length(pages)) do.call(rbind, pages) else {
    as.data.frame(stats::setNames(rep(list(character()), length(fields)), fields))
  }
  list(data = data, truncated = truncated)
}

# ---- 1. Test with a small query ------------------------------------------------
note("Start. Queries: ", queries_csv, "; output: ", output)
test <- download_query("Chaoborus", "eats", fields)
if (is.null(test)) {
  note("The test with all the fields failed; trying with the minimal fields")
  fields <- minimal_fields
  test <- download_query("Chaoborus", "eats", fields)
}
if (is.null(test) || nrow(test$data) == 0) {
  stop("GloBI does not answer or returns no records for 'Chaoborus'. Is there a connection? ",
       "See ", log_file, call. = FALSE)
}
missing_fields <- setdiff(minimal_fields, names(test$data))
if (length(missing_fields)) {
  stop("The API did not return the fields: ", paste(missing_fields, collapse = ", "),
       ". Fields received: ", paste(names(test$data), collapse = ", "), call. = FALSE)
}
n_cited <- sum(!is.na(test$data$study_citation))
note("Test 'Chaoborus': ", nrow(test$data), " records; with citation: ", n_cited,
     "; types: ", paste(names(table(test$data$interaction_type)),
                        table(test$data$interaction_type), collapse = ", "))
if (n_cited == 0) {
  stop("The records come without citation (empty study_citation): stopping before the download.",
       call. = FALSE)
}
# Does the server limit the size of the page? If the first page brings fewer
# than page_size records but there are more after it, the limit is the server's.
p1 <- read_page(page_url("Daphnia", "eats", 0L, fields, api))
if (!is.null(p1) && nrow(p1) > 0L && nrow(p1) < page_size) {
  p2 <- read_page(page_url("Daphnia", "eats", nrow(p1), fields, api))
  if (!is.null(p2) && nrow(p2) > 0L) {
    page_size <- nrow(p1)
    note("The server returns at most ", page_size, " records per page: that size is used")
  }
}
# 'eats' should include 'preysOn' (a subtype); if not, both types are queried.
e_eats <- read_page(page_url("Esox lucius", "eats", 0L, fields, api))
e_prey <- read_page(page_url("Esox lucius", "preysOn", 0L, fields, api))
included <- !is.null(e_eats) && !is.null(e_prey) && nrow(e_prey) > 0L &&
  "preysOn" %in% e_eats$interaction_type
types <- if (included) "eats" else c("eats", "preysOn")
note("Interaction types queried: ", paste(types, collapse = ", "),
     "; fields: ", paste(fields, collapse = ", "), "; records per page: ", page_size)

# ---- 2. Download -------------------------------------------------------------
queries <- utils::read.csv(queries_csv, stringsAsFactors = FALSE)$query
tasks <- expand.grid(query = queries, type = types, stringsAsFactors = FALSE)
summary_table <- data.frame(query = tasks$query, type = tasks$type, records = NA_integer_,
                            truncated = FALSE, failed = FALSE, stringsAsFactors = FALSE)
t0 <- Sys.time()
n_new <- 0L
for (k in seq_len(nrow(tasks))) {
  q <- tasks$query[k]
  type <- tasks$type[k]
  f <- cache_file(q, type)
  if (file.exists(f)) {
    d <- utils::read.csv(f, colClasses = "character", na.strings = c("", "NA"),
                         encoding = "UTF-8", check.names = FALSE)
    summary_table$records[k] <- nrow(d)
    summary_table$truncated[k] <- is_truncated(f)
    next
  }
  res <- download_query(q, type, fields)
  n_new <- n_new + 1L
  if (is.null(res)) {
    summary_table$failed[k] <- TRUE
    note(sprintf("[%d/%d] %s (%s): FAILED, it will be tried again in the next run",
                 k, nrow(tasks), q, type))
    next
  }
  utils::write.csv(res$data, f, row.names = FALSE, fileEncoding = "UTF-8")
  if (res$truncated) {
    writeLines("truncated", paste0(f, ".truncated"))
    note(sprintf("[%d/%d] %s (%s): WARNING, limit of %d records reached",
                 k, nrow(tasks), q, type, page_size * max_pages))
  }
  summary_table$records[k] <- nrow(res$data)
  summary_table$truncated[k] <- res$truncated
  min_per_query <- as.numeric(difftime(Sys.time(), t0, units = "mins")) / n_new
  left <- sum(!file.exists(vapply(seq_len(nrow(tasks)),
                                  function(i) cache_file(tasks$query[i], tasks$type[i]),
                                  character(1))))
  note(sprintf("[%d/%d] %s (%s): %d records; about %.0f min left",
               k, nrow(tasks), q, type, nrow(res$data), min_per_query * left))
}

# ---- 3. Join and summary -----------------------------------------------------
utils::write.csv(summary_table, file.path(output, "query_summary.csv"), row.names = FALSE)
if (any(summary_table$failed)) {
  note(sum(summary_table$failed), " queries failed: run the script again to retry them ",
       "(the downloaded ones are not repeated).")
}
parts <- lapply(seq_len(nrow(tasks)), function(k) {
  f <- cache_file(tasks$query[k], tasks$type[k])
  if (!file.exists(f)) return(NULL)
  d <- utils::read.csv(f, colClasses = "character", na.strings = c("", "NA"),
                       encoding = "UTF-8", check.names = FALSE)
  if (!nrow(d)) return(NULL)
  d$query <- tasks$query[k]
  d$query_type <- tasks$type[k]
  d
})
parts <- parts[!vapply(parts, is.null, logical(1))]
cols <- unique(unlist(lapply(parts, names)))
parts <- lapply(parts, function(d) {
  for (col in setdiff(cols, names(d))) d[[col]] <- NA_character_
  d[cols]
})
raw <- do.call(rbind, parts)
attr(raw, "download") <- list(date = Sys.time(), api = api, fields = fields, types = types,
                              R = R.version.string)
saveRDS(raw, file.path(output, "globi_raw.rds"))
note("Finished: ", nrow(raw), " records of ", length(parts), " queries in ",
     file.path(output, "globi_raw.rds"), if (any(summary_table$failed)) " (incomplete: some failed)" else "")
