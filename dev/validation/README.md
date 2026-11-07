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
| 9 | `data-raw/gray2015_calibration.R`, `data-raw/gray2015_globi_calibration.R <globi_registry.rds>` | default rule | `data/gray2015_calibration.rda`, `data/gray2015_globi_calibration.rda` | 5 and 6 min |

Other scripts: `scoring_sensitivity.R` (scoring with a background probability
and with maximum likelihood fits), `validation_groups.R` (webs of the same
system in Gray), `web_coverage.R` (predictive coverage of the calibration),
`calibration_lme4.R` (maximum likelihood fit against lme4), `example_system.R`
(example) and `reanalysis_metaweb.R` (webs sampled with the metaweb of the
reanalysis; needs the local copies of `dev/reference`).

The results and their interpretation are in `results.md`; the decisions, in
`../decisions.md`.

The GloBI registry is not included in the package (the licences of its data
depend on the source); the package includes only the scores.
