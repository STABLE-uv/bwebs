# bwebs 0.2.0

## Changes that alter results

* `calibrate_links()` fits the calibration model as a Bayesian model by
  default (`method = "bayes"`): normal priors on the class coefficients and a
  half-normal prior on the standard deviation of the web effect, posterior
  sampled by Metropolis-Hastings. `sample_metawebs()` takes the parameters of
  each draw from the posterior draws. The previous fit is `method = "ml"`.
* The calibration reports `p_new`, the probability of a class in a new
  system averaged over web effects, used to score rules and metawebs.
* `matching_rule()` uses every registry at every level by default
  (`primary_registry = NULL`; downward inheritance from every registry): with
  the Gray et al. (2015) registry combined with GloBI this was the best of
  the sources compared (`source_comparison`). Results with a single registry
  do not change; `primary_registry = "gray2015"` gives the previous
  behaviour with several registries.
* Data: `gray2015_globi_calibration`, the calibration of the Gray et al.
  (2015) registry combined with GloBI, and `source_comparison`, the scores of
  the sources compared.

## New features

* `validate_rule()`: the leave-one-system-out validation of a matching rule
  (previously inside `calibrate_links()`, which now accepts its result).
* `score_rule()`: cross-validated Brier score, Brier skill score and TSS of a
  rule, refitting the calibration without the system of each web.
* `candidate_rules()` and `compare_rules()`: grid of matching rules (fixed,
  climbing and WebBuilder ladder, several starting levels and thresholds)
  ranked by their cross-validated Brier score.
* `score_metaweb()`: scores of a probabilistic metaweb against an observed
  food web.
* Data: `tuesday_lake_taxa`, `tuesday_lake_community` and
  `tuesday_lake_links`, the Tuesday Lake food webs of 1984 and 1986 (from the
  cheddar package), as an independent test. Study 12 of `gray2015` is Tuesday
  Lake and must be left out of the registry for that test.
* `validate_rule()` and `compare_rules()` gain `train_registries` (validate
  one registry, such as GloBI, on the webs of another) and `exclusions`
  (studies left out when a given web is rebuilt, such as copies of that web
  in other compilations).
* `shared_studies()`: proportion of the links of each web recorded in each
  study of a registry, to find studies of the same system and copies.

## Other changes

* `prepare_nodes()` applies `default_taxonomy_fixes()` to the classification
  of the nodes (new argument `taxonomy_fixes`), the same harmonisation the
  registry receives. `default_taxonomy_fixes()` adds synonyms of classes used
  by other classifications, such as those of GloBI (`"Teleostei"` and
  `"Actinopteri"` for `"Actinopterygii"`, `"Petromyzonti"` for
  `"Cephalaspidomorphi"`); the bundled registry does not change.
* `globi_to_registry()` reads placeholders of unresolved ranks
  (`"Not assigned"`, `"unclassified ..."`, `"... incertae sedis"`) as
  missing, and reads the ranks of each distinct path once.
* `globi_to_registry()` and `as_registry()` return an empty registry when
  no record is left (for example when no record of GloBI belongs to the
  taxon queried), instead of failing.
* The matching is faster and needs less memory with large registries:
  records that match no node are dropped first, the records of each name and
  rank are matched once, the node x record matrices are filled in place and
  the links are assembled without dense matrix products. The metawebs are
  identical (the comparison with the original code still gives the same
  1319 links, levels and records).
* `sample_metawebs()` and `score_metaweb()` warn when the calibration was
  fitted with other registries than those of the metaweb.
* Andreu Rico is a coauthor of the package; `CITATION.cff` describes how to
  cite it.

# bwebs 0.1.0

First usable version, written from the code of the reanalysis of the
azoxystrobin mesocosm experiment of Amador et al. (2024).

## Changes that alter results

* Compared with the original analysis, the food webs are no longer only the
  projection of one deterministic metaweb: `sample_metawebs()` keeps each
  inferred link with a probability calibrated by its class of evidence, so
  the webs are sparser than the metaweb. The deterministic metaweb is still
  available with `build_metaweb()` and `calc_indicators()`.

## New features

* `matching_rule()`, `prepare_nodes()` and `build_metaweb()`: the matching
  engine of the original analysis (WebBuilder ladder with climbing, evidence
  threshold, surrogates, aliases, life stages, downward inheritance, generic
  resources and manual links), with the class of evidence of each link.
* `calibrate_links()`: leave-one-web-out calibration of link probabilities
  by class of evidence, with a logistic model with a random effect of the web;
  `groups` leaves out together the studies of the same system. The bundled
  calibration leaves out studies 13 and 50 together (62% of the links of
  study 13 are also recorded in study 50).
* `sample_metawebs()`, `calc_indicators()`, `web_indicators()`,
  `summarize_indicators()` and `contrast_indicators()`.
* `as_registry()`, `combine_registries()`, `fetch_globi()` and
  `globi_to_registry()`.
* Data: `gray2015` (registry of Gray et al. 2015), `gray2015_sources`,
  `gray2015_calibration` and the invented example community.
