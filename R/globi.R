# Global Biotic Interactions (GloBI) -------------------------------------------

#' Download trophic interactions from GloBI
#'
#' `fetch_globi()` queries the web API of Global Biotic Interactions for the
#' trophic interactions in which each taxon is the consumer, page by page,
#' and returns the raw records. Convert them with [globi_to_registry()].
#' Requires an internet connection.
#'
#' @param taxa Character vector of taxon names to query as consumers (for
#'   example the names, genera and families of the consumer nodes).
#' @param interaction_types GloBI interaction types to query (default
#'   `"eats"`). In GloBI, `"eats"` also returns the records typed as
#'   `"preysOn"` and as parasitism; [globi_to_registry()] keeps the trophic
#'   types.
#' @param cache_dir Optional directory where every query is saved as a CSV
#'   file and read from on later calls (delete a file to download it again).
#'   `NULL` (default) does not write anything.
#' @param page_size,max_pages Records per page and maximum pages per query.
#' @param pause Seconds to wait between requests.
#' @param api Base URL of the interaction endpoint.
#' @param quiet If `FALSE`, report progress with messages.
#'
#' @return A data frame with one row per record and the character columns
#'   returned by the API (`source_taxon_name`, `source_taxon_path`,
#'   `source_taxon_path_ranks`, `interaction_type`, `target_taxon_name`,
#'   `target_taxon_path`, `target_taxon_path_ranks`, `study_citation`,
#'   `study_source_citation`, `study_source_id`, `study_doi`), plus `query`
#'   (the taxon queried) and `query_type`. Records are requested one by one
#'   (`includeObservations=true`), which is needed for the citations. A
#'   query that fails is reported with a warning and is not cached, so it can
#'   be retried.
#'
#' @section Assumptions:
#' A query matches the name at any rank of any classification linked by
#' GloBI, so it can return records of homonyms in other groups;
#' [globi_to_registry()] keeps, by default, only the records whose consumer
#' path contains the name queried. GloBI is rebuilt about weekly: record the
#' date of the download. GloBI indexes the Gray et al. (2015) and Brose et
#' al. (2005) compilations; [globi_to_registry()] drops them by default to
#' avoid counting the same records twice when combining with [gray2015].
#'
#' @references
#' Poelen JH, Simons JD, Mungall CJ (2014). Global Biotic Interactions: an
#' open infrastructure to share and analyze species-interaction datasets.
#' *Ecological Informatics*.
#' \doi{10.1016/j.ecoinf.2014.08.005}
#'
#' @export
#' @examplesIf interactive()
#' raw <- fetch_globi(c("Asellus", "Gammarus"))
#' globi <- globi_to_registry(raw)
fetch_globi <- function(taxa, interaction_types = "eats", cache_dir = NULL,
                        page_size = 1000L, max_pages = 30L, pause = 0.5,
                        api = "https://api.globalbioticinteractions.org/interaction",
                        quiet = FALSE) {
  if (!is.character(taxa) || length(taxa) == 0L || anyNA(taxa)) {
    abort_input("`taxa` must be a character vector of taxon names.")
  }
  if (!is.character(interaction_types) || length(interaction_types) == 0L) {
    abort_input("`interaction_types` must be a character vector, such as c(\"eats\", \"preysOn\").")
  }
  check_string(cache_dir, allow_null = TRUE)
  check_number(page_size, lower = 1)
  check_number(max_pages, lower = 1)
  check_number(pause, lower = 0)
  check_string(api)
  check_flag(quiet)
  if (!is.null(cache_dir) && !dir.exists(cache_dir)) {
    abort_input(sprintf("`cache_dir` does not exist: %s.", cache_dir))
  }
  out <- list()
  for (q in unique(taxa)) {
    for (type in interaction_types) {
      if (!quiet) message(sprintf("GloBI: %s (%s)", q, type))
      res <- fetch_globi_one(q, type, cache_dir, as.integer(page_size), as.integer(max_pages),
                             pause, api)
      if (!is.null(res) && nrow(res)) {
        res$query <- rep(q, nrow(res))
        res$query_type <- rep(type, nrow(res))
        out[[length(out) + 1L]] <- res
      }
    }
  }
  if (length(out) == 0L) {
    empty <- as.data.frame(stats::setNames(rep(list(character()), length(globi_fields) + 2L),
                                           c(globi_fields, "query", "query_type")),
                           stringsAsFactors = FALSE)
    return(empty)
  }
  cols <- unique(unlist(lapply(out, names)))
  out <- lapply(out, function(d) {
    for (col in setdiff(cols, names(d))) d[[col]] <- rep(NA_character_, nrow(d))
    d[cols]
  })
  do.call(rbind, out)
}

globi_fields <- c("source_taxon_name", "source_taxon_path", "source_taxon_path_ranks",
                  "interaction_type", "target_taxon_name", "target_taxon_path",
                  "target_taxon_path_ranks", "study_citation", "study_source_citation",
                  "study_source_id", "study_doi")

fetch_globi_one <- function(taxon, type, cache_dir, page_size, max_pages, pause, api) {
  cache <- if (is.null(cache_dir)) NULL else {
    file.path(cache_dir, paste0(gsub("[^A-Za-z0-9]+", "_", taxon), "_", type, ".csv"))
  }
  if (!is.null(cache) && file.exists(cache)) {
    return(utils::read.csv(cache, colClasses = "character", na.strings = c("", "NA")))
  }
  pages <- list()
  for (p in seq_len(max_pages)) {
    url <- paste0(api, "?sourceTaxon=", utils::URLencode(taxon, reserved = TRUE),
                  "&interactionType=", type, "&includeObservations=true",
                  "&type=csv&limit=", page_size,
                  "&offset=", (p - 1L) * page_size,
                  "&fields=", paste(globi_fields, collapse = ","))
    res <- NULL
    for (attempt in 1:3) {
      res <- tryCatch(utils::read.csv(url, colClasses = "character", na.strings = c("", "NA")),
                      error = function(e) NULL)
      if (!is.null(res)) break
      Sys.sleep(2 * attempt)
    }
    if (is.null(res)) {
      warning(sprintf("GloBI: query %s (%s) failed at page %d; it is not cached.",
                      taxon, type, p), call. = FALSE)
      return(NULL)
    }
    pages[[p]] <- res
    if (pause > 0) Sys.sleep(pause)
    if (nrow(res) < page_size) break
  }
  out <- do.call(rbind, pages)
  if (!is.null(cache)) utils::write.csv(out, cache, row.names = FALSE)
  out
}

#' Convert GloBI records into a registry
#'
#' `globi_to_registry()` turns the records of [fetch_globi()] into the
#' registry format of [as_registry()]: the target of each interaction is the
#' resource and the source is the consumer; the classification of both is
#' read from the GloBI taxon paths.
#'
#' @param raw Records from [fetch_globi()].
#' @param registry Name of the registry (default `"globi"`).
#' @param types Interaction types kept (default `"eats"` and `"preysOn"`;
#'   parasitism and other types are dropped).
#' @param exclude_sources Values of `study_source_id` dropped: by default the
#'   GloBI copies of the Gray et al. (2015) and Brose et al. (2005)
#'   compilations, which [gray2015] already contains.
#' @param match_query If `TRUE` (default), keep only records whose consumer
#'   path contains the name queried (GloBI matches names approximately).
#' @param exclude Regular expressions; records whose citation matches any of
#'   them are dropped (for example inferred webs or generic diet matrices).
#' @param taxonomy_fixes Passed to [as_registry()].
#'
#' @return A registry (see [as_registry()]). The study of each record is its
#'   citation (`study_citation`, or the dataset citation when missing);
#'   records without any citation share the study
#'   `"globi_without_citation"`. Duplicated records (same resource, consumer
#'   and study) are kept once.
#'
#' @export
#' @examples
#' ranks <- "kingdom | phylum | class | order | family | genus | species"
#' raw <- data.frame(
#'   source_taxon_name = "Asellus aquaticus",
#'   source_taxon_path = paste("Animalia | Arthropoda | Malacostraca | Isopoda |",
#'                             "Asellidae | Asellus | Asellus aquaticus"),
#'   source_taxon_path_ranks = ranks,
#'   interaction_type = "eats", target_taxon_name = "Cocconeis placentula",
#'   target_taxon_path = paste("Chromista | Ochrophyta | Bacillariophyceae | Achnanthales |",
#'                             "Cocconeidaceae | Cocconeis | Cocconeis placentula"),
#'   target_taxon_path_ranks = ranks,
#'   study_citation = "An invented study", query = "Asellus", query_type = "eats"
#' )
#' globi_to_registry(raw)
globi_to_registry <- function(raw, registry = "globi", types = c("eats", "preysOn"),
                              exclude_sources = c("globalbioticinteractions/gray2015",
                                                  "globalbioticinteractions/brose"),
                              match_query = TRUE, exclude = character(0),
                              taxonomy_fixes = default_taxonomy_fixes()) {
  check_data_frame(raw)
  base_fields <- c("source_taxon_name", "source_taxon_path", "source_taxon_path_ranks",
                   "interaction_type", "target_taxon_name", "target_taxon_path",
                   "target_taxon_path_ranks", "study_citation", "query")
  check_columns(raw, base_fields, arg = "raw")
  if (!is.character(types) || !length(types)) abort_input("`types` must be a character vector.")
  if (!is.character(exclude_sources)) abort_input("`exclude_sources` must be a character vector.")
  check_string(registry)
  check_flag(match_query)
  if (!is.character(exclude)) abort_input("`exclude` must be a character vector of regular expressions.")
  d <- as_plain_data_frame(raw)
  for (col in c("study_source_citation", "study_source_id")) {
    if (!col %in% names(d)) d[[col]] <- rep(NA_character_, nrow(d))
  }
  d <- d[!is.na(d$source_taxon_name) & !is.na(d$target_taxon_name) &
           !is.na(d$source_taxon_path) & d$interaction_type %in% types, , drop = FALSE]
  dup_src <- !is.na(d$study_source_id) & d$study_source_id %in% exclude_sources
  if (any(dup_src)) {
    message(sprintf("globi_to_registry(): %d record(s) of excluded compilations dropped.",
                    sum(dup_src)))
    d <- d[!dup_src, , drop = FALSE]
  }
  if (match_query && nrow(d)) {
    key <- paste(tolower(d$query), d$source_taxon_path, sep = "\u001f")
    first <- which(!duplicated(key))
    ok <- vapply(first, function(i) {
      tolower(d$query[i]) %in% tolower(split_path(d$source_taxon_path[i]))
    }, logical(1))
    in_path <- ok[match(key, key[first])]
    if (any(!in_path)) {
      message(sprintf("globi_to_registry(): %d record(s) of other taxa than the one queried dropped.",
                      sum(!in_path)))
    }
    d <- d[in_path, , drop = FALSE]
  }
  cit <- ifelse(is.na(d$study_citation), d$study_source_citation, d$study_citation)
  if (length(exclude) && nrow(d)) {
    drop <- !is.na(cit) & grepl(paste(exclude, collapse = "|"), cit)
    if (any(drop)) message(sprintf("globi_to_registry(): %d record(s) excluded by `exclude`.", sum(drop)))
    d <- d[!drop, , drop = FALSE]
    cit <- cit[!drop]
  }
  res <- ranks_from_paths(d$target_taxon_path, d$target_taxon_path_ranks)
  con <- ranks_from_paths(d$source_taxon_path, d$source_taxon_path_ranks)
  reg <- data.frame(
    study_id = ifelse(is.na(cit), "globi_without_citation", cit),
    resource = d$target_taxon_name, consumer = d$source_taxon_name,
    res_genus = res$genus, res_subfamily = res$subfamily, res_family = res$family,
    res_order = res$order, res_class = res$class,
    con_genus = con$genus, con_subfamily = con$subfamily, con_family = con$family,
    con_order = con$order, con_class = con$class,
    evidence = sprintf("GloBI: %s", d$interaction_type), citation = cit,
    stringsAsFactors = FALSE
  )
  reg <- reg[!duplicated(reg[c("resource", "consumer", "study_id")]), , drop = FALSE]
  as_registry(reg, registry = registry, taxonomy_fixes = taxonomy_fixes)
}

# Genus, subfamily, family, order and class of each path, computed once per
# distinct path (GloBI downloads repeat the same paths many times).
ranks_from_paths <- function(path, ranks) {
  key <- paste(path, ranks, sep = "\u001f")
  first <- which(!duplicated(key))
  m <- match(key, key[first])
  out <- lapply(c("genus", "subfamily", "family", "order", "class"), function(r) {
    rank_from_path(path[first], ranks[first], r)[m]
  })
  stats::setNames(out, c("genus", "subfamily", "family", "order", "class"))
}

# Value of a rank from GloBI paths ("Animalia | Arthropoda | ...") and their
# ranks ("kingdom | phylum | ..."); NA when the rank is absent, holds a
# placeholder ("Not assigned", "unclassified ...", "... incertae sedis") or
# the two strings do not have the same number of elements.
rank_from_path <- function(path, ranks, rank) {
  out <- vapply(seq_along(path), function(i) {
    p <- path[i]
    r <- ranks[i]
    if (is.na(p) || is.na(r)) {
      return(NA_character_)
    }
    pv <- split_path(p)
    rv <- tolower(split_path(r))
    if (length(pv) != length(rv)) {
      return(NA_character_)
    }
    v <- pv[rv == rank]
    if (length(v) == 0L || v[1] == "") NA_character_ else v[1]
  }, character(1))
  out[is_placeholder(out)] <- NA_character_
  out
}

# Names that stand for an unresolved rank in GloBI classifications.
is_placeholder <- function(x) {
  !is.na(x) & grepl("^not assigned$|^unclassified\\b|incertae sedis|^unidentified|^unknown",
                    x, ignore.case = TRUE)
}

# Split a GloBI path on "|" keeping empty elements, also the trailing one
# (strsplit() drops it), so that paths and ranks stay aligned.
split_path <- function(x) {
  parts <- strsplit(paste0(x, "|\u0001"), "|", fixed = TRUE)[[1]]
  trimws(parts[-length(parts)])
}
