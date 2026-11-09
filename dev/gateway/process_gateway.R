# GATEWAy as a registry: the freshwater links (lakes and streams) of the
# GlobAL daTabasE of traits and food Web Architecture (GATEWAy) version 1.0
# (Brose et al. 2019, Nature Ecology & Evolution 3: 919-927; iDiv Data
# Repository 283, CC BY 4.0) in the format of as_registry(). GATEWAy gives
# only the names of the taxa; their genus, subfamily, family, order and class
# are taken from the Gray et al. (2015) and GloBI registries (by name, then
# by genus or higher rank) or, for the rest, from the GBIF backbone (cached
# in dev/gateway/gateway_taxonomy.csv, so that a new run needs no network).
#
# Left out:
#   - Tuesday Lake (Jonsson et al. 2005), the independent test of the
#     package, and any study that shares at least 10% of its links (as for
#     GloBI, dev/globi/process_globi.R);
#   - the five chalk streams of Gray et al. (2015) and Thompson et al.
#     (2017), whose links were partly extrapolated from similar taxa with
#     WebBuilder ("gut content, extrapolated from similar taxa");
#   - parasitic and parasitoid links.
# Every food web is a study (study_id "gateway: <food web>"). The studies that
# share at least 10% of the links of a validation web of Gray are left out
# when that web is validated (exclusions, as for GloBI).
#
# Outputs:
#   <folder of the csv>/gateway_registry.rds: list(gateway = the registry,
#     gateway_test = without the studies of Tuesday Lake, exclusions = data
#     frame web / study_id, shared_validation, tuesday_lake, source);
#   dev/gateway/gateway_taxonomy.csv: the taxonomy of every name and where it
#     comes from (registry, genus, family, order, class, gbif or none);
#   dev/gateway/gateway_tuesday_lake.csv, dev/gateway/gateway_exclusions.csv.
#
# Use, from the root of the package:
#   Rscript dev/gateway/process_gateway.R <283_2_FoodWebDataBase_2018_12_10.csv> <globi_registry.rds>
# Packages used only here: pkgload, jsonlite (GBIF, only for names not yet in
# the cache).
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Use: Rscript dev/gateway/process_gateway.R <GATEWAy csv> <globi_registry.rds>")
}
pkgload::load_all(".", quiet = TRUE)
out_dir <- file.path("dev", "gateway")
tax_file <- file.path(out_dir, "gateway_taxonomy.csv")

g <- utils::read.csv(args[1], stringsAsFactors = FALSE, na.strings = c("NA", "-999", ""))
fw <- g[g$ecosystem.type %in% c("lakes", "streams"), ]
fw <- fw[!grepl("^Tuesday Lake", fw$foodweb.name), ]
fw <- fw[!fw$link.citation %in% "Gray et al. (2015), Thompson et al. (2017)", ]
fw <- fw[!fw$interaction.type %in% c("parasitic", "parasitoid"), ]
message(sprintf("%d freshwater links of %d food webs", nrow(fw), length(unique(fw$foodweb.name))))

# Names as recorded, without tabs, encoding artefacts or repeated spaces
clean <- function(x) {
  x <- gsub("[\t Ê]+", " ", x)
  trimws(gsub(" +", " ", x))
}
fw$con <- clean(fw$con.taxonomy)
fw$res <- clean(fw$res.taxonomy)

# Name to look up: without qualifiers, stages and morphotypes, with the
# English forms of family names ("Hydrobiosid" -> "Hydrobiosidae")
query_name <- function(x) {
  # groups resolved to their class
  groups <- c("fish fry" = "Actinopterygii", "fish eggs" = "Actinopterygii",
              "copepod nauplii" = "Maxillopoda")
  x <- ifelse(tolower(x) %in% names(groups), groups[tolower(x)], x)
  q <- gsub("\\(.*?\\)", " ", x)
  q <- gsub("\\b(indet|spp|sp|cf|new|type|larvae|larva|adults|adult|nymphs|nymph|pupae|juveniles|nauplii)\\b", " ",
            q, ignore.case = TRUE)
  q <- gsub("[.?]", " ", q)
  q <- trimws(gsub(" +", " ", q))
  q <- gsub("\\s(I|II|III|IV|[A-Z])$", "", q)
  q <- sub("/.*$", "", q)
  fixes <- c(Oligochaete = "Oligochaeta", Pelicypod = "Bivalvia", Pelecypod = "Bivalvia",
             Scrtid = "Scirtidae", Potamopurgus = "Potamopyrgus", Psilachorema = "Psilochorema",
             Archicauliodes = "Archichauliodes", Anguila = "Anguilla", Acnanthidium = "Achnanthidium",
             Tanytarsini = "Chironomidae", Hexatomini = "Limoniidae", Eriopterini = "Limoniidae",
             Phaoniinae = "Muscidae", Cocconeia = "Cocconeis", Acnantidium = "Achnanthidium",
             Sephanodiscus = "Stephanodiscus", Aucalosira = "Aulacoseira", Illybius = "Ilybius",
             Scapaloberis = "Scapholeberis", Polypedellum = "Polypedilum",
             Eukiefidrella = "Eukiefferiella", Zealolessica = "Zelolessica", Keriphyron = "Kephyrion")
  w <- strsplit(q, " ", fixed = TRUE)
  vapply(w, function(v) {
    v <- v[nzchar(v)]
    if (!length(v)) return(NA_character_)
    if (v[1] %in% names(fixes)) v[1] <- fixes[[v[1]]]
    if (length(v) == 1L && grepl("id$", v[1])) v[1] <- paste0(v[1], "ae")
    v <- v[!grepl("^(broad|skinny|Broad|Skinny|confusus)$", v) | seq_along(v) == 1L]
    paste(v[seq_len(min(2L, length(v)))], collapse = " ")
  }, character(1))
}

# Taxonomy of the Gray et al. (2015) and GloBI registries
reg <- combine_registries(gray2015, readRDS(args[2])$globi)
ranks <- c("genus", "subfamily", "family", "order", "class")
side <- function(prefix, name_col) {
  out <- reg[c(name_col, paste0(prefix, "_", ranks))]
  names(out) <- c("name", ranks)
  out
}
lk <- rbind(side("con", "consumer"), side("res", "resource"), side("con", "consumer_std"),
            side("res", "resource_std"))
lk$name <- tolower(lk$name)
lk <- lk[!is.na(lk$class) | !is.na(lk$order) | !is.na(lk$family), ]
# the classification of each name (or rank) with the most ranks resolved and,
# among those, the most frequent in the records
most_common <- function(d, by, cols) {
  d <- d[!is.na(d[[by]]), unique(c(by, cols))]
  key <- do.call(paste, c(d, sep = "\u001f"))
  n <- table(key)
  d$n <- as.integer(n[key])
  d$resolved <- rowSums(!is.na(d[cols]))
  d <- d[order(d[[by]], -d$resolved, -d$n), ]
  d[!duplicated(d[[by]]), unique(c(by, cols))]
}
by_name <- most_common(lk, "name", ranks)
by_genus <- most_common(lk, "genus", ranks)
by_subfamily <- most_common(lk, "subfamily", c("family", "order", "class"))
by_family <- most_common(lk, "family", c("order", "class"))
by_order <- most_common(lk, "order", "class")
known_class <- unique(stats::na.omit(lk$class))

# GBIF backbone for the rest (cached)
cache <- if (file.exists(tax_file)) utils::read.csv(tax_file, stringsAsFactors = FALSE) else NULL
# Accepted: exact matches, and fuzzy matches (spelling) within two letters of
# the name looked up, of the kingdom and phylum expected from the metabolic
# type of GATEWAy and resolved at least to class. The name is looked up in
# the kingdoms expected (Animalia for animals; Plantae, Chromista, Bacteria
# and Protozoa for producers), then by its first word (the genus).
gbif_one <- function(q, type, kingdom) {
  animal <- type %in% c("invertebrate", "ectotherm vertebrate", "endotherm vertebrate")
  url <- paste0("https://api.gbif.org/v1/species/match?verbose=false&name=",
                utils::URLencode(q, reserved = TRUE), "&kingdom=", kingdom)
  x <- tryCatch(jsonlite::fromJSON(url), error = function(e) NULL)
  get <- function(f) if (is.null(x[[f]])) NA_character_ else as.character(x[[f]])
  if (is.null(x) || !get("matchType") %in% c("EXACT", "FUZZY")) return(NULL)
  if (get("matchType") == "FUZZY" &&
      (is.na(get("canonicalName")) || utils::adist(tolower(q), tolower(get("canonicalName"))) > 2 ||
       as.numeric(get("confidence")) < 80)) {
    return(NULL)
  }
  vertebrate <- type %in% c("ectotherm vertebrate", "endotherm vertebrate")
  if (is.na(get("class")) || (animal && !get("kingdom") %in% "Animalia") ||
      (!animal && get("kingdom") %in% "Animalia") ||
      (vertebrate && !get("phylum") %in% "Chordata") ||
      (type %in% "invertebrate" && get("phylum") %in% "Chordata")) {
    return(NULL)
  }
  rk <- toupper(get("rank"))
  c(genus = if (rk %in% c("GENUS", "SPECIES", "SUBSPECIES", "VARIETY", "FORM")) get("genus") else NA,
    subfamily = NA_character_, family = get("family"), order = get("order"), class = get("class"),
    source = "gbif")
}
gbif_match <- function(q, type) {
  animal <- type %in% c("invertebrate", "ectotherm vertebrate", "endotherm vertebrate")
  kingdoms <- if (animal) "Animalia" else c("Plantae", "Chromista", "Bacteria", "Protozoa")
  for (qq in unique(c(q, sub(" .*", "", q)))) {
    for (k in kingdoms) {
      out <- gbif_one(qq, type, k)
      if (!is.null(out)) return(out)
    }
  }
  c(stats::setNames(rep(NA_character_, 5L), ranks), source = "none")
}

names_tab <- unique(data.frame(
  name = c(fw$con, fw$res),
  type = c(fw$con.metabolic.type, fw$res.metabolic.type),
  level = c(fw$con.taxonomy.level, fw$res.taxonomy.level), stringsAsFactors = FALSE
))
names_tab <- names_tab[!is.na(names_tab$name) & !duplicated(names_tab$name), ]
names_tab$query <- query_name(names_tab$name)
resolve <- function(name, query, type, level) {
  for (n in unique(stats::na.omit(tolower(c(name, query))))) {
    i <- match(n, by_name$name)
    if (!is.na(i)) return(c(unlist(by_name[i, ranks]), source = "registry"))
  }
  if (is.na(query)) return(c(stats::setNames(rep(NA_character_, 5L), ranks), source = "none"))
  first <- sub(" .*", "", query)
  if (grepl(" ", query)) {
    i <- match(first, by_genus$genus)
    if (!is.na(i)) return(c(unlist(by_genus[i, ranks]), source = "genus"))
  }
  i <- match(first, by_subfamily$subfamily)
  if (!is.na(i)) return(c(genus = NA, subfamily = first,
                          unlist(by_subfamily[i, c("family", "order", "class")]), source = "subfamily"))
  i <- match(first, by_family$family)
  if (!is.na(i)) return(c(genus = NA, subfamily = NA, unlist(by_family[i, c("family", "order", "class")]),
                          source = "family"))
  i <- match(first, by_order$order)
  if (!is.na(i)) return(c(genus = NA, subfamily = NA, family = NA, unlist(by_order[i, c("order", "class")]),
                          source = "order"))
  if (first %in% known_class) return(c(genus = NA, subfamily = NA, family = NA, order = NA, class = first,
                                       source = "class"))
  # names without a taxonomic level in GATEWAy (common names, such as the
  # local names of the fish of Lake Malawi, and groups such as zooplankton)
  # are left unresolved
  if (is.na(level)) return(c(stats::setNames(rep(NA_character_, 5L), ranks), source = "none"))
  i <- match(name, cache$name)
  if (!is.na(i) && cache$source[i] %in% c("gbif", "none")) {
    return(stats::setNames(as.character(unlist(cache[i, c(ranks, "source")])), c(ranks, "source")))
  }
  gbif_match(query, type)
}
res <- t(vapply(seq_len(nrow(names_tab)), function(k) {
  resolve(names_tab$name[k], names_tab$query[k], names_tab$type[k], names_tab$level[k])
}, character(6)))
tax <- cbind(names_tab[c("name", "query", "type", "level")], as.data.frame(res, stringsAsFactors = FALSE))
# higher ranks missing in the source, from the lower ones
fill <- function(tax, from, lookup, cols) {
  i <- match(tax[[from]], lookup[[from]])
  for (cl in cols) {
    miss <- is.na(tax[[cl]]) & !is.na(i)
    tax[[cl]][miss] <- lookup[[cl]][i[miss]]
  }
  tax
}
tax <- fill(tax, "genus", by_genus, c("family", "order", "class"))
tax <- fill(tax, "family", by_family, c("order", "class"))
tax <- fill(tax, "order", by_order, "class")
utils::write.csv(tax, tax_file, row.names = FALSE)
print(table(tax$source))

# Records in the format of as_registry()
stage <- function(x) {
  out <- rep(NA_character_, length(x))
  out[x %in% c("larvae", "nymphs", "nauplii")] <- "larva"
  out[x %in% "adults"] <- "adult"
  out
}
category <- function(x) {
  map <- c("primary producer" = "producer", invertebrate = "invertebrate",
           "ectotherm vertebrate" = "vert.ecto", "endotherm vertebrate" = "vert.endo",
           detritus = "detritus", "dead organic material" = "detritus")
  unname(map[x])
}
evidence <- function(x) {
  out <- rep(NA_character_, length(x))
  obs <- grepl("gut|field|feeding trial|stomach", x, ignore.case = TRUE)
  out[obs] <- paste("observed -", x[obs])
  out[!obs & grepl("published account", x)] <- "published account"
  out
}
ti <- match(fw$con, tax$name)
tr <- match(fw$res, tax$name)
rec <- data.frame(
  study_id = paste("gateway:", fw$foodweb.name),
  resource = fw$res, consumer = fw$con,
  res_stage = stage(fw$res.lifestage), con_stage = stage(fw$con.lifestage),
  res_genus = tax$genus[tr], res_subfamily = tax$subfamily[tr], res_family = tax$family[tr],
  res_order = tax$order[tr], res_class = tax$class[tr],
  con_genus = tax$genus[ti], con_subfamily = tax$subfamily[ti], con_family = tax$family[ti],
  con_order = tax$order[ti], con_class = tax$class[ti],
  res_category = category(fw$res.metabolic.type), con_category = category(fw$con.metabolic.type),
  evidence = evidence(fw$link.methodology),
  citation = paste0(fw$link.citation, ", in GATEWAy version 1.0 (Brose et al. 2019)"),
  stringsAsFactors = FALSE
)
rec <- rec[!is.na(rec$resource) & !is.na(rec$consumer), ]
gateway <- as_registry(rec, registry = "gateway")
message(sprintf("registry: %d records, %d studies", nrow(gateway), length(unique(gateway$study_id))))

# Studies that copy Tuesday Lake
name <- stats::setNames(tuesday_lake_taxa$taxon, tuesday_lake_taxa$node_id)
tl <- data.frame(web = tuesday_lake_links$sample_id,
                 resource = unname(name[tuesday_lake_links$resource]),
                 consumer = unname(name[tuesday_lake_links$consumer]), stringsAsFactors = FALSE)
sh_tl <- shared_studies(gateway, tl, threshold = 0.1)
drop_tl <- unique(sh_tl$study_id)
gateway_test <- gateway[!gateway$study_id %in% drop_tl, , drop = FALSE]
message(sprintf("studies sharing >= 10%% of the links of Tuesday Lake: %d", length(drop_tl)))

# Studies to leave out when each validation web of Gray is validated
val <- gray2015_sources$study_id[gray2015_sources$validation]
val <- setdiff(val, "12")
both <- combine_registries(gray2015[gray2015$study_id != "12", ], gateway_test)
sh_val <- shared_studies(both, val, threshold = 0.05)
sh_val <- sh_val[sh_val$study_id %in% gateway_test$study_id, , drop = FALSE]
exclusions <- sh_val[sh_val$share >= 0.1, c("web", "study_id"), drop = FALSE]
# The webs of a publication pooled: a validation web can be spread over
# several webs of GATEWAy (the 20 streams of Layer et al. 2010 are web 27 of
# Gray), so every web of a publication that holds at least 10% of the links
# of a validation web, pooled, is left out with it
cit <- sub(", in GATEWAy.*$", "", gateway_test$citation)
pooled <- gateway_test
pooled$study_id <- paste("gateway publication:", cit)
both_p <- combine_registries(gray2015[gray2015$study_id != "12", ], pooled)
sh_pub <- shared_studies(both_p, val, threshold = 0.05)
sh_pub <- sh_pub[sh_pub$study_id %in% pooled$study_id & sh_pub$share >= 0.1, , drop = FALSE]
print(sh_pub[c("web", "study_id", "n_links", "n_shared", "share")], digits = 2)
by_pub <- do.call(rbind, lapply(seq_len(nrow(sh_pub)), function(k) {
  st <- unique(gateway_test$study_id[pooled$study_id == sh_pub$study_id[k]])
  data.frame(web = sh_pub$web[k], study_id = st, stringsAsFactors = FALSE)
}))
exclusions <- unique(rbind(exclusions, by_pub))
exclusions <- exclusions[order(exclusions$web, exclusions$study_id), ]
rownames(exclusions) <- NULL
print(table(exclusions$web))

utils::write.csv(sh_tl, file.path(out_dir, "gateway_tuesday_lake.csv"), row.names = FALSE)
utils::write.csv(exclusions, file.path(out_dir, "gateway_exclusions.csv"), row.names = FALSE)
saveRDS(list(gateway = gateway, gateway_test = gateway_test, exclusions = exclusions,
             shared_validation = sh_val, shared_publications = sh_pub, tuesday_lake = sh_tl,
             source = list(file = basename(args[1]), date = Sys.time(), R = R.version.string)),
        file.path(dirname(args[1]), "gateway_registry.rds"))
