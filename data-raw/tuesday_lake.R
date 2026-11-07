# Tuesday Lake 1984 and 1986 -> data/tuesday_lake_taxa.rda,
# data/tuesday_lake_community.rda, data/tuesday_lake_links.rda
#
# Source: communities TL84 and TL86 of the cheddar package (version 0.1-639,
# BSD 2-Clause licence, see inst/COPYRIGHTS); data of Jonsson, Cohen &
# Carpenter (2005) and Cohen, Jonsson & Carpenter (2003). Packages used only
# here: cheddar, pkgload, usethis.
#
# The node table is aligned with the registry of Gray et al. (2015): the
# rotifer class (spelt "Monogonata" in cheddar) is Monogononta, and the rank
# names go through default_taxonomy_fixes(), as the registry does. Study 12
# of the registry ("from Brose et al 2005 (Jonsson et al. (2004))") is
# Tuesday Lake: leave it out of the registry when Tuesday Lake is used to
# evaluate a method.
#
# Run from the package root: Rscript data-raw/tuesday_lake.R

pkgload::load_all(".", quiet = TRUE)
stopifnot(requireNamespace("cheddar", quietly = TRUE))
e <- new.env()
utils::data("TL84", "TL86", package = "cheddar", envir = e)
webs <- list(TL84 = e$TL84, TL86 = e$TL86)

node_id <- function(x) gsub("^_|_$", "", gsub("[^a-z0-9]+", "_", tolower(x)))
fixes <- c(Monogonata = "Monogononta", default_taxonomy_fixes())

nodes <- do.call(rbind, lapply(names(webs), function(w) {
  n <- cheddar::NPS(webs[[w]])
  data.frame(sample_id = w, node = n$node, category = n$category, M = n$M, N = n$N,
             phylum = n$phylum, class = n$class, order = n$order, family = n$family,
             genus = n$genus, species = n$species, stringsAsFactors = FALSE)
}))
blank <- function(x) ifelse(is.na(x) | trimws(x) == "", NA_character_, trimws(x))
for (col in c("phylum", "class", "order", "family", "genus", "species")) nodes[[col]] <- blank(nodes[[col]])
for (col in c("class", "order", "family")) nodes[[col]] <- apply_fixes(apply_fixes(nodes[[col]], fixes), fixes)
nodes$node_id <- node_id(nodes$node)

taxa <- nodes[!duplicated(nodes$node_id), ]
tuesday_lake_taxa <- data.frame(
  node_id = taxa$node_id,
  taxon = ifelse(!is.na(taxa$species), paste(taxa$genus, taxa$species),
                 ifelse(!is.na(taxa$genus), taxa$genus, taxa$node)),
  rank = ifelse(!is.na(taxa$species), "species", ifelse(!is.na(taxa$genus), "genus", "non_taxonomic")),
  role = ifelse(taxa$category == "producer", "producer", "consumer"),
  node_type = ifelse(is.na(taxa$genus), "non_taxonomic", "taxon"),
  genus = taxa$genus, subfamily = NA_character_, family = taxa$family, order = taxa$order,
  class = taxa$class, phylum = taxa$phylum, category = taxa$category, name = taxa$node,
  stringsAsFactors = FALSE
)
tuesday_lake_taxa <- tuesday_lake_taxa[order(tuesday_lake_taxa$role, tuesday_lake_taxa$node_id,
                                             method = "radix"), ]
rownames(tuesday_lake_taxa) <- NULL

tuesday_lake_community <- data.frame(sample_id = nodes$sample_id, node_id = nodes$node_id,
                                     value = nodes$N, mass_kg = nodes$M,
                                     stringsAsFactors = FALSE)

tuesday_lake_links <- do.call(rbind, lapply(names(webs), function(w) {
  l <- cheddar::TLPS(webs[[w]])
  data.frame(sample_id = w, resource = node_id(l$resource), consumer = node_id(l$consumer),
             stringsAsFactors = FALSE)
}))
stopifnot(all(c(tuesday_lake_links$resource, tuesday_lake_links$consumer) %in%
                tuesday_lake_taxa$node_id))

print(table(tuesday_lake_community$sample_id))
print(table(tuesday_lake_links$sample_id))
print(table(tuesday_lake_taxa$role, tuesday_lake_taxa$rank, useNA = "ifany"))
usethis::use_data(tuesday_lake_taxa, tuesday_lake_community, tuesday_lake_links,
                  overwrite = TRUE, compress = "xz")
