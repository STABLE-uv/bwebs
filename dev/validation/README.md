# Analysis of the method (dev/validation)

Scripts that produce the results of the method of bwebs: choice of the
matching rule, comparison of sources, test on Tuesday Lake and the tables
included in the package. They are run from the root of the package, in this
order. The times are those of a laptop with two cores (Ubuntu, R 4.3).

| Step | Script | Input | Output | Time |
|---|---|---|---|---|
| 1 | `dev/globi/globi_queries.R` (needs the local copies of `dev/reference`) | validation webs, Tuesday Lake and nodes of the azoxy experiment | `dev/globi/globi_queries.csv` | seconds |
| 2 | `dev/globi/download_globi.R` (on a computer with access to GloBI) | queries | `globi_raw.rds` (832,003 records, 500 queries) | ~1 h |
| 3 | `dev/globi/process_globi.R <globi_raw.rds>` | download | `globi_registry.rds` (next to the download), `dev/globi/globi_tuesday_lake.csv`, `dev/globi/globi_exclusions.csv` | 1 min |
| 4 | `compare_rules.R 1 2`, `2 2`, `merge 2` | `gray2015` | `rule_comparison_gray.rds` (21 rules, Gray only) | 30 min in two processes |
| 5 | `compare_sources.R <globi_registry.rds> source globi`, `source both_family`, `source both_class`, `source both_all 1 2` and `2 2`, `join both_all 2`, `merge` | steps 3 and 4 | `source_comparison_*.rds`, `source_comparison.rds` | ~1 h 40 min in two processes |
| 6 | `compare_sources.R <globi_registry.rds> tuesday` (or `tuesday <sources>` to add only those) | step 5 | `tuesday_lake_sources.rds` | ~30 min |
| 6b | `source_sensitivity.R <globi_registry.rds>` and `... default` | step 3 | `source_sensitivity.rds` | ~35 min |
| 7 | `tuesday_lake.R` | step 4 | `tuesday_lake.rds` (21 rules, Gray only) | 15 min |
| 8 | `data-raw/gray2015_rule_comparison.R`, `data-raw/source_comparison.R` | steps 4 to 7 | `data/gray2015_rule_comparison.rda`, `data/source_comparison.rda` | seconds |
| 9 | `data-raw/gray2015_calibration.R`, `data-raw/gray2015_globi_calibration.R <globi_registry.rds>` | default rule | `data/gray2015_calibration.rda`, `data/gray2015_globi_calibration.rda` | ~25 and ~30 min |
| 10 | `decay_background.R <globi_registry.rds>` (taxonomic distance decay, experimental) | steps 4 and 6 | `decay_background.rds` | ~25 min |
| 11 | `decay_package.R <globi_registry.rds> gray`, `... both_all`, `... gray_gateway` (after step 17) and `... merge` (the taxonomic decay of the package; Tuesday Lake also with the body masses of `tuesday_lake_community`) | step 3 | `decay_package.rds` | ~15 min per registry |
| 11b | `decay_role.R <globi_registry.rds>` (needs lme4) and `decay_web_size.R <globi_registry.rds>` (variants of the decay: role of the resource, size of the web) | step 3 | `decay_role.rds`, `decay_web_size.rds` | ~10 and ~30 min |
| 12 | `size_brose.R <bodysizes_2008.txt> <globi_registry.rds>` (first exploration of the size term) | Brose et al. (2005), step 3 | `size_brose.rds` | minutes |
| 13 | `size_folds.R <bodysizes_2008.txt> <cheddar data folder> <globi_registry.rds>`, then `size_lambda.R`, `size_sources.R <bodysizes_2008.txt> <cheddar data folder>` and `size_consumer_type.R <bodysizes_2008.txt> <cheddar data folder>` (validation by folds of the size term, its weights, the webs of its curve and curves by kind of predator) | Brose et al. (2005), the data files of the source of cheddar, step 3 | `size_folds.rds`, `size_sources.rds`, `size_consumer_type.rds` | ~35 min, then ~15 min each |
| 14 | `dev/gateway/process_gateway.R <GATEWAy csv> <globi_registry.rds>` (GATEWAy version 1.0 as a registry; the csv of the iDiv Data Repository, dataset 283) | step 3 | `gateway_registry.rds` (next to the csv), `dev/gateway/gateway_taxonomy.csv`, `dev/gateway/gateway_exclusions.csv`, `dev/gateway/gateway_tuesday_lake.csv` | 2 min |
| 15 | `compare_gateway.R <globi_registry.rds> <gateway_registry.rds> [registries]`, then `... merge` (registries with GATEWAy) | steps 3, 11 and 14 | `compare_gateway.rds` | ~10 min per registry |
| 16 | `size_gateway.R <GATEWAy csv> <gateway_registry.rds> <bodysizes_2008.txt> <cheddar data folder> <globi_registry.rds> [registries]`, then `... merge` (the size term with the masses and webs of GATEWAy) | steps 3 and 14 | `size_gateway.rds` | ~5 min per registry |
| 17 | `data-raw/gateway.R <gateway_registry.rds>`, `data-raw/gray2015_gateway_calibration.R`, `data-raw/size_model.R <GATEWAy csv>` | steps 14 to 16 | `data/gateway.rda`, `data/gray2015_gateway_calibration.rda`, `data/size_model.rda` | ~15 min |
| 18 | `evidence_amount.R <globi_registry.rds> gray`, `... both_all`, `... gray_gateway`, `... tuesday gray` (and the other two), then `... merge` (amount of evidence of the links and ordered classes; Tuesday Lake with the ordered classes) | steps 3, 11 and 17 | `evidence_amount.rds` | ~10 min per registry and ~20 min per registry for Tuesday Lake |
| 19 | `data-raw/gray2015_calibration.R`, `data-raw/gray2015_globi_calibration.R <globi_registry.rds>`, `data-raw/gray2015_gateway_calibration.R` again (the bundled calibrations with the ordered classes) | step 18 | `data/gray2015_calibration.rda`, `data/gray2015_globi_calibration.rda`, `data/gray2015_gateway_calibration.rda` | ~15, ~25 and ~15 min |

The validation scripts of steps 4 to 16 fit the free classes of evidence
(`monotone = FALSE`), the model of those phases; the scripts of `data-raw/`
use the defaults of the package, whose classes are ordered from step 18
(step 19 refits the bundled calibrations with them).

Other scripts: `scoring_sensitivity.R` (scoring with a background probability
and with maximum likelihood fits), `validation_groups.R` (webs of the same
system in Gray), `web_coverage.R` (predictive coverage of the calibration),
`calibration_lme4.R` (maximum likelihood fit against lme4), `example_system.R`
(example) and `reanalysis_metaweb.R` (webs sampled with the metaweb of the
reanalysis; needs the local copies of `dev/reference`).

The results and their interpretation are in `results.md`; the decisions, in
`../decisions.md`.

The GloBI registry is not included in the package (the licences of its data
depend on the source); the package includes only the scores. Neither are
the body masses of Brose et al. (2005) and of cheddar. GATEWAy (CC BY 4.0)
is included as a registry (`gateway`); of its body masses, only the
coefficients of the size term (`size_model`).
