# Registries of trophic interactions -------------------------------------------

#' Prepare a registry of published trophic interactions
#'
#' `as_registry()` validates a table of trophic records (one row per record
#' of a consumer eating a resource in a study) and harmonises it into the
#' format used by [build_metaweb()] and [calibrate_links()]. The bundled
#' [gray2015] registry is already in this format.
#'
#' @param data A data frame with columns `study_id` (identifier of the study
#'   or data set the record comes from), `resource` and `consumer` (taxon
#'   names as recorded), and the classification of both partners:
#'   `res_genus`, `res_subfamily`, `res_family`, `res_order`, `res_class`,
#'   `con_genus`, `con_subfamily`, `con_family`, `con_order`, `con_class`
#'   (`NA` where not resolved). Optional columns: `record_id`, `res_stage`
#'   and `con_stage` (life stage as recorded), `res_category` and
#'   `con_category` (for example `"producer"`), `evidence` and `citation`.
#' @param registry Name of the registry (a single string, such as
#'   `"globi"`), stored in the `registry` column. Must not be empty.
#' @param taxonomy_fixes Named character vector of names to replace in the
#'   rank columns (names: as found; values: replacement), used to align
#'   misspellings and synonyms of higher ranks with the other registries and
#'   the node classification. `default_taxonomy_fixes()` gives the
#'   harmonisation used for the Gray et al. (2015) registry; `character(0)`
#'   changes nothing.
#'
#' @return A data frame with one row per record in the input order and
#'   columns `record_id`, `registry`, `study_id`, `resource`, `consumer`,
#'   `resource_std` and `consumer_std` (names after `taxonomy_fixes`, used
#'   to recognise records resolved only to a coarse rank), `res_stage` and
#'   `con_stage` (`"larva"`, `"adult"` or `NA` when absent or mixed), the
#'   ten rank columns, `res_category`, `con_category`, `evidence` and
#'   `citation` (all character). Records without a resource or consumer name
#'   are dropped with a message.
#'
#' @section Missing values and zeros:
#' A missing rank (`NA`) means the rank was not resolved in the source; it is
#' never filled in.
#'
#' @seealso [combine_registries()], [fetch_globi()].
#' @export
#' @examples
#' rec <- data.frame(
#'   study_id = "s1", resource = "Baetis rhodani", consumer = "Perla bipunctata",
#'   res_genus = "Baetis", res_subfamily = NA, res_family = "Baetidae",
#'   res_order = "Ephemeroptera", res_class = "Insecta",
#'   con_genus = "Perla", con_subfamily = NA, con_family = "Perlidae",
#'   con_order = "Plecoptera", con_class = "Insecta"
#' )
#' as_registry(rec, registry = "my_registry")
as_registry <- function(data, registry, taxonomy_fixes = default_taxonomy_fixes()) {
  check_data_frame(data)
  check_string(registry)
  if (!is.character(taxonomy_fixes) || (length(taxonomy_fixes) && is.null(names(taxonomy_fixes)))) {
    abort_input("`taxonomy_fixes` must be a named character vector (names: as found; values: replacement).")
  }
  required <- c("study_id", "resource", "consumer",
                paste0("res_", rank_cols), paste0("con_", rank_cols))
  check_columns(data, stats::setNames(required, rep("", length(required))))
  optional <- c("record_id", "res_stage", "con_stage", "res_category", "con_category",
                "evidence", "citation")
  reg <- as_plain_data_frame(data)
  for (col in optional) {
    if (!col %in% names(reg)) reg[[col]] <- rep(NA_character_, nrow(reg))
  }
  for (col in c(required, optional)) {
    v <- reg[[col]]
    if (is.factor(v) || is.numeric(v) || is.logical(v)) v <- as.character(v)
    if (!is.character(v)) {
      abort_input(sprintf("Column '%s' of `data` must be text, not %s.", col, class(v)[1L]))
    }
    v <- trimws(v)
    v[!is.na(v) & v == ""] <- NA_character_
    reg[[col]] <- v
  }
  drop <- is.na(reg$resource) | is.na(reg$consumer)
  if (any(drop)) {
    message(sprintf("as_registry(): %d record(s) without resource or consumer name dropped.",
                    sum(drop)))
    reg <- reg[!drop, , drop = FALSE]
  }
  if (anyNA(reg$study_id)) {
    abort_input(sprintf(
      "Column 'study_id' has missing values in %s; every record needs the study it comes from.",
      format_rows(which(is.na(reg$study_id)))
    ))
  }
  for (col in c(paste0("res_", rank_cols), paste0("con_", rank_cols))) {
    reg[[col]] <- apply_fixes(reg[[col]], taxonomy_fixes)
  }
  reg$resource_std <- apply_fixes(reg$resource, taxonomy_fixes)
  reg$consumer_std <- apply_fixes(reg$consumer, taxonomy_fixes)
  reg$res_stage <- normalise_stage(reg$res_stage)
  reg$con_stage <- normalise_stage(reg$con_stage)
  if (all(is.na(reg$record_id))) reg$record_id <- paste0(registry, "_", seq_len(nrow(reg)))
  reg$registry <- rep(registry, nrow(reg))
  out <- reg[registry_columns]
  rownames(out) <- NULL
  out
}

registry_columns <- c(
  "record_id", "registry", "study_id", "resource", "consumer", "resource_std",
  "consumer_std", "res_stage", "con_stage", paste0("res_", c("genus", "subfamily",
  "family", "order", "class")), paste0("con_", c("genus", "subfamily", "family",
  "order", "class")), "res_category", "con_category", "evidence", "citation"
)

#' Combine registries
#'
#' Binds registries prepared with [as_registry()] (or the bundled
#' [gray2015]). Records keep the name of their registry, which
#' [matching_rule()] uses to decide which registry is primary.
#'
#' @param ... Registries (data frames in the format of [as_registry()]).
#'
#' @return A data frame in the same format, with the records of each
#'   registry in the order given. Study identifiers are prefixed with the
#'   registry name when two registries share an identifier, so that studies
#'   of different registries are never merged.
#'
#' @section Assumptions:
#' A study that appears in two registries (for example a data set indexed
#' both in the Gray et al. (2015) compilation and in GloBI) is counted as two
#' studies; remove the duplicate records before combining if this matters.
#'
#' @export
#' @examples
#' extra <- gray2015[gray2015$study_id == "1", ]
#' extra$registry <- "copy"
#' nrow(combine_registries(gray2015[1:10, ], extra))
combine_registries <- function(...) {
  regs <- list(...)
  if (length(regs) == 0L) {
    abort_input("Give at least one registry.")
  }
  for (k in seq_along(regs)) {
    check_registry(regs[[k]], arg = sprintf("registry %d", k))
  }
  ids <- unlist(lapply(regs, function(r) unique(r$study_id)), use.names = FALSE)
  names_reg <- unlist(lapply(regs, function(r) unique(r$registry)), use.names = FALSE)
  if (anyDuplicated(names_reg)) {
    abort_input(sprintf(
      "Registries must have different names in their `registry` column; repeated: %s.",
      format_items(unique(names_reg[duplicated(names_reg)]))
    ))
  }
  if (anyDuplicated(ids)) {
    regs <- lapply(regs, function(r) {
      r$study_id <- paste(r$registry, r$study_id, sep = ":")
      r
    })
  }
  out <- do.call(rbind, lapply(regs, function(r) r[registry_columns]))
  rownames(out) <- NULL
  out
}

check_registry <- function(reg, arg = "registry", call = sys.call(-1)) {
  check_data_frame(reg, arg = arg, allow_empty = FALSE, call = call)
  check_columns(reg, registry_columns, arg = arg, call = call)
  if (anyNA(reg$registry) || anyNA(reg$study_id)) {
    abort_input(sprintf("`%s` has records without registry or study_id; prepare it with as_registry().",
                        arg), call = call)
  }
  invisible(reg)
}

#' Harmonisation of rank names used for the Gray et al. (2015) registry
#'
#' Misspellings in the registry and synonyms of higher ranks used by other
#' taxonomies (such as the classifications that GloBI attaches to its
#' records), aligned with the names of the Gray et al. (2015) registry.
#' [as_registry()] applies them to the classification and names of the
#' records, and [prepare_nodes()] to the classification of the nodes.
#'
#' @return A named character vector: names as found, values as used.
#'
#' @details
#' Only names that are not used in the Gray et al. (2015) registry are
#' replaced, so the registry itself is not changed. Synonyms are given when
#' the two names cover practically the same taxa in freshwater webs (for
#' example `"Teleostei"` and `"Actinopteri"` for `"Actinopterygii"`); broader
#' or ambiguous names (`"Osteichthyes"`, `"Hexapoda"`) are left as they are.
#' @export
#' @examples
#' default_taxonomy_fixes()
default_taxonomy_fixes <- function() {
  c(
    "Acnanthales" = "Achnanthales",           # misspelling in the registry
    "Gomophonemataceae" = "Gomphonemataceae", # misspelling in the registry
    "Cladocera" = "Diplostraca", "Anomopoda" = "Diplostraca",
    "Ctenopoda" = "Diplostraca", "Onychopoda" = "Diplostraca",
    "Hexanauplia" = "Maxillopoda", "Copepoda" = "Maxillopoda",
    "Thecostraca" = "Maxillopoda",
    "Monogononta" = "Eurotatoria",
    "Monogonta" = "Eurotatoria",              # misspelling found in GloBI
    "Euglenoidea" = "Euglenida", "Euglenophyceae" = "Euglenida",
    "Heteroptera" = "Hemiptera",
    "Rhabditophora" = "Turbellaria",
    # classes of other classifications (for example, in GloBI)
    "Teleostei" = "Actinopterygii", "Actinopteri" = "Actinopterygii",
    "Chondrostei" = "Actinopterygii",
    "Petromyzonti" = "Cephalaspidomorphi", "Hyperoartia" = "Cephalaspidomorphi",
    "Lepidosauria" = "Reptilia",
    "Conjugatophyceae" = "Zygnematophyceae", "Zygnemophyceae" = "Zygnematophyceae",
    "Synurophyceae" = "Chrysophyceae", "Mediophyceae" = "Coscinodiscophyceae",
    "Cyanobacteriia" = "Cyanophyceae"
  )
}

apply_fixes <- function(x, fixes) {
  if (length(fixes) == 0L) {
    return(x)
  }
  hit <- !is.na(x) & x %in% names(fixes)
  x[hit] <- unname(fixes[x[hit]])
  x
}

# Life stage as larva / adult / NA (mixed or absent), as in the original
# analysis.
normalise_stage <- function(x) {
  x <- tolower(x)
  out <- rep(NA_character_, length(x))
  adult <- !is.na(x) & grepl("adult", x)
  larva <- !is.na(x) & grepl("larv|nymph", x)
  mixed <- adult & (!is.na(x) & grepl("larv|juven|nymph", x))
  out[larva] <- "larva"
  out[adult] <- "adult"
  out[mixed] <- NA_character_
  out
}
