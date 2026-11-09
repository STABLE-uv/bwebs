# gateway: the freshwater links of GATEWAy version 1.0 as a registry
# -> data/gateway.rda
#
# Source: Brose U, et al. (2019) Predator traits determine food-web
# architecture across ecosystems. Nature Ecology & Evolution 3: 919-927.
# GlobAL daTabasE of traits and food Web Architecture (GATEWAy) version 1.0,
# iDiv Data Repository, dataset 283, version 3 (CC BY 4.0;
# 283_2_FoodWebDataBase_2018_12_10.csv).
#
# The lakes and streams of GATEWAy, processed by dev/gateway/process_gateway.R
# (taxonomy of the names from the Gray et al. 2015 and GloBI registries or
# the GBIF backbone, dev/gateway/gateway_taxonomy.csv), without Tuesday Lake
# (the independent test of the package), the five chalk streams of Gray et
# al. (2015) and Thompson et al. (2017) (partly extrapolated with
# WebBuilder) and parasitic and parasitoid links. The studies of GATEWAy that
# copy each validation web of Gray are in dev/gateway/gateway_exclusions.csv.
# Packages used only here: usethis.
#
# Run from the package root after dev/gateway/process_gateway.R:
#   Rscript data-raw/gateway.R <path>/gateway_registry.rds
x <- readRDS(commandArgs(trailingOnly = TRUE)[1])
gateway <- x$gateway_test
rownames(gateway) <- NULL
cat(nrow(gateway), "records,", length(unique(gateway$study_id)), "studies\n")
usethis::use_data(gateway, overwrite = TRUE, compress = "xz")
