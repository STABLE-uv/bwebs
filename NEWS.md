# bwebs 0.3.0

## Changes that alter results

* Taxonomic distance decay of the pairs without an inferred link
  (`matching_rule(decay = TRUE)`, the default). Instead of probability 0,
  each such pair takes the probability of its evidence: the closest level
  (family, order, class, phylum or any consumer) at which relatives of the
  consumer have recorded diets, and the closest level at which the resource
  is related to those diets (within the family of a recorded resource,
  order, class, phylum or not related), plus whether a record of a generic
  group applies to it. Closer evidence never gives a lower probability.
  `validate_rule()` returns the evidence of the pairs not inferred of each
  validation web (`decay`); `calibrate_links()` fits a second hierarchical
  model to them, with its own web effect (`decay`); `build_metaweb()` lists
  the pairs without a link with their evidence (`decay`);
  `sample_metawebs()` draws them in every draw (`background`), and
  `calc_indicators()`, `score_rule()` and `score_metaweb()` include them.
  `decay = FALSE` gives the results of 0.2.0. The decay model includes the
  size of the web (log of its number of nodes): with Gray + GloBI, fewer of
  the pairs not inferred are links in larger webs. In the cross-validation
  of the default rule with Gray + GloBI, the Brier score went from 0.1988
  without the decay to 0.1983 with it and to 0.1971 with the size of the
  web; with the Gray et al. (2015) registry alone, from 0.223 to 0.208 and
  0.209.
* Body masses (`sample_metawebs(masses = )`, `score_metaweb(masses = )`): the
  probability of every predation pair (an animal resource) whose two nodes
  have a mass is lowered by how far the predator/prey mass ratio is from the
  optimum of `size_model`, a within-consumer effect fitted to the 206
  aquatic webs of GATEWAy (Brose et al. 2019; Tuesday Lake left out), with
  one weight for the links inferred by the rule and another for the pairs
  of the taxonomic decay. The masses can be one per node or one set per draw
  (the draws of a biomass model): each draw uses its own masses and keeps
  them (`masses` of the result), so that a draw holds the links of a web and
  the masses of its nodes. The curve and the weights (0.25 for the inferred
  links and 1 for the pairs of the decay) are those that gave the best
  validation by folds of the 12 validation webs with body masses, left out
  by system, among curves fitted to the webs of Brose et al. (2005) and
  cheddar or to those of GATEWAy, one curve or one per kind of predator:
  their mean Brier score went down by 0.006 (9 of 12 webs better; 0.008 with
  the Gray registry, 0.005 with Gray + GloBI, 0.006 with Gray + GATEWAy), and
  that of Tuesday Lake from 0.140 to 0.137 with Gray + GloBI and from 0.149
  to 0.139 with Gray + GATEWAy. Without `masses` nothing changes.
* GATEWAy (Brose et al. 2019, CC BY 4.0) as a registry: `gateway`, the
  19,185 links of its 81 freshwater webs (lakes and streams, without
  Tuesday Lake), classified with the taxonomy of the other registries or of
  the GBIF backbone, and `gray2015_gateway_calibration`, the calibration of
  `combine_registries(gray2015, gateway)`. In the cross-validation it
  improved the Gray registry (Brier 0.2093 to 0.2052, 7 of 12 webs better)
  but not Gray + GloBI (0.1971 to 0.2013), which stays the best combination;
  being distributable, Gray + GATEWAy is the registry of the examples.
* Ordered classes of evidence (`calibrate_links(monotone = TRUE)`, the
  default; also in `score_rule()`): a link supported by a closer relative
  (the taxon itself, its genus, subfamily, family, order, class) is never
  less probable than one supported by a more distant relative; `generic`
  is free. The free coefficients of 0.2.0 (`monotone = FALSE`) broke that
  order: links of the same genus were more probable than those of the same
  taxon, links that needed a class more than those that needed an order,
  and with Gray + GloBI links that needed an order (0.51) more than those
  of the same family (0.30). In the cross-validation of the 12 webs the
  ordered classes scored the same (Brier score with the decay 0.2104
  against 0.2093 with Gray, 0.1978 against 0.1971 with Gray + GloBI and
  0.2044 against 0.2052 with Gray + GATEWAy; standard errors of the
  differences 0.003 to 0.005); on Tuesday Lake they did better with Gray +
  GloBI and Gray + GATEWAy and worse with the Gray registry alone. With
  Gray + GATEWAy (`gray2015_gateway_calibration`) the probability in a new
  system is 0.51 for the taxon and its genus, 0.38 for the subfamily and
  the family, and 0.31 for the order and the class.
* The bundled calibrations `gray2015_calibration`,
  `gray2015_globi_calibration` and `gray2015_gateway_calibration` have the
  ordered classes and the decay model (with the size of the web).

## New features

* `link_probabilities()`: the matrix of the probability of every resource →
  consumer pair of a metaweb, from which `sample_metawebs()` draws its webs
  (the inferred links by their class of evidence, the other pairs by the
  taxonomic decay and, with `masses`, the size term), averaged over the
  posterior and over the effect of the system.
* The amount of evidence of each link: `build_metaweb()` keeps, besides its
  records and studies, the consumer and resource taxa among them
  (`n_consumer_taxa`, `n_resource_taxa`: how many relatives of the consumer
  eat the resource or its relatives, and how many relatives of the
  resource are eaten) and the relatives of each node with records in the
  network (`n_consumer_known`, `n_resource_known`); `validate_rule()`
  returns them for every link inferred in each validation web (`links`).
  None of these amounts, nor the shares of the relatives that support a
  link, improved the prediction of the links of the webs left out in the
  cross-validation of the default rule with Gray, Gray + GloBI or Gray +
  GATEWAy (their effects changed sign from one registry to another), so the
  probability of a link still depends on its class only.

## Other changes

* Logo, shown in the README (with a version for dark backgrounds).
* `score_rule()` and `score_metaweb()` report the recorded links that the
  rule missed (`n_missed`) and the links expected among the pairs without a
  link (`expected_background`).
* Data: `size_model`, the effect of the predator/prey body-mass ratio
  (coefficients only); `gateway` and `gray2015_gateway_calibration`.

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
