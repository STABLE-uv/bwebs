# Phylum of each class (internal table `class_phylum`, R/sysdata.rda)
#
# The taxonomic distances of the decay model (decay_profiles()) go up to the
# phylum, but the registries and the node tables classify only up to the
# class. The phylum of each class is taken from the taxonomic paths of the
# GloBI download (dev/globi/download_globi.R; the records are not
# redistributed, only this table of class names): the most frequent phylum
# of each class, after giving the phyla the names of the GBIF backbone, which
# Tuesday Lake and default_generic_groups() use (Ochrophyta for diatoms and
# chrysophytes, Myzozoa for dinoflagellates...). Classes of the registries or
# of Tuesday Lake that are missing in GloBI are added by hand. The class
# names are those of default_taxonomy_fixes().
#
# Run from the package root:
#   Rscript data-raw/class_phylum.R <globi_raw.rds>
pkgload::load_all(".", quiet = TRUE)
args <- commandArgs(trailingOnly = TRUE)
raw <- readRDS(args[1])
paths <- unique(rbind(
  data.frame(path = raw$source_taxon_path, ranks = raw$source_taxon_path_ranks),
  data.frame(path = raw$target_taxon_path, ranks = raw$target_taxon_path_ranks)
))
paths <- paths[!is.na(paths$path) & !is.na(paths$ranks), , drop = FALSE]
d <- data.frame(class = apply_fixes(rank_from_path(paths$path, paths$ranks, "class"),
                                    default_taxonomy_fixes()),
                phylum = trimws(sub(";.*", "", rank_from_path(paths$path, paths$ranks, "phylum"))),
                n = 1L)
d <- d[!is.na(d$class) & !is.na(d$phylum), , drop = FALSE]
# Phylum names of the GBIF backbone
backbone <- c(Bacillariophyta = "Ochrophyta", Heterokontophyta = "Ochrophyta",
              Chrysophyta = "Ochrophyta", Cryptista = "Cryptophyta",
              Cyanobacteriota = "Cyanobacteria", Miozoa = "Myzozoa", Dinoflagellata = "Myzozoa",
              Magnoliophyta = "Tracheophyta", Pseudomonadota = "Proteobacteria",
              Actinomycetota = "Actinobacteria", Bryophytes = "Bryophyta", Radiolaria = "Radiozoa")
d$phylum <- ifelse(d$phylum %in% names(backbone), backbone[d$phylum], d$phylum)
tab <- stats::aggregate(n ~ class + phylum, d, sum)
# Streptophyta is a clade above the phyla of the backbone
tab <- tab[!(tab$phylum == "Streptophyta" & tab$class %in% tab$class[tab$phylum != "Streptophyta"]), ]
tab <- tab[order(tab$class, -tab$n, tab$phylum, method = "radix"), ]
tab <- tab[!duplicated(tab$class), c("class", "phylum")]
# Classes of the registries and of Tuesday Lake that are not in the GloBI paths
by_hand <- c(
  Acantharia = "Radiozoa", Adenophorea = "Nematoda", Anopla = "Nemertea",
  Dicotyledones = "Tracheophyta", Eustigmatophyceae = "Ochrophyta",
  Glomeromycetes = "Glomeromycota", Isokontae = "Chlorophyta", Nuda = "Ctenophora",
  Pennatae = "Ochrophyta", Polythalamea = "Foraminifera", Priapulomorpha = "Priapulida",
  Pyramimonadophyceae = "Chlorophyta", Raphidophyceae = "Ochrophyta",
  Synurophyceae = "Ochrophyta"
)
by_hand <- by_hand[!names(by_hand) %in% tab$class]
class_phylum <- rbind(tab, data.frame(class = names(by_hand), phylum = unname(by_hand)))
class_phylum <- class_phylum[order(class_phylum$class, method = "radix"), ]
rownames(class_phylum) <- NULL
print(table(class_phylum$phylum))
usethis::use_data(class_phylum, internal = TRUE, overwrite = TRUE)
