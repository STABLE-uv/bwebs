# Validation — bwebs 0.2.0

Results obtained on 2026-11-03 with the development version 0.1.0.9000 (the
same code as version 0.2.0), on Ubuntu 24.04 with R 4.3.3.

## Refactoring: equality with the original code

`dev/reference/run_original.R` runs `01_derive_foodwebs.R` and
`02_network_properties.R` of the original project (Azoxy_fw) without editing them
(only `setwd()` neutralised). The regenerated metaweb matches the metaweb of
the reanalysis byte for byte (`metaweb_adjacency.csv`, `metaweb_edges.csv`).

| Check | Result |
|---|---|
| `build_metaweb()` against `metaweb_edges.csv` (108 nodes) | 1319 links, same order; `support`, `n_records`, `n_sources`, `res_level`, `con_level`, `evidence`, `source_ids`, `examples`, `citations`, `cannibal` and adjacency identical; levels used per node identical |
| `web_indicators()` against the regenerated `web_metrics.csv` (120 webs) | 27 indicators equal with `all.equal(tolerance = 1e-10)`, including modularity and R50 with the seed 20240101 |

## Calibration model

| Check | Result |
|---|---|
| Maximum likelihood fit (`method = "ml"`) against `lme4::glmer(nAGQ = 20)` on the counts of the bundled calibration (default rule, 13 webs; `dev/validation/calibration_lme4.R`) | coefficients equal to within 3e-6; standard errors to within 7e-6; σ 1.399165 against 1.399169 |
| Posterior (bundled calibration) against maximum likelihood | posterior means 0.13 above the estimates in every class; posterior standard deviations 0.36 against standard errors 0.40; posterior mean of σ 1.58 against 1.40 |
| Marginal likelihood against numerical integration (`integrate()`) | equal with tolerance 1e-8 (test) |
| Recovery of simulated parameters (12 webs, 3 classes) | within 4 standard errors (test) |
| Predictive coverage: each validation web against the model fitted to the other 12 (unknown system, 4000 simulations; `dev/validation/web_coverage.R`, bundled calibration) | 12/13 webs within the 95% interval (0.92; Monte Carlo error 0.07). Outside: Jones (1951), precision 0.03 against [0.05; 0.87]. The intervals are very wide: the precision of a particular system can hardly be anticipated |

## Webs of the same system (`dev/validation/validation_groups.R`)

Proportion of the links of each validation web that are also recorded in
another study of the registry (the two studies with the most):

| Web | Links | Study 1 | Proportion | Study 2 | Proportion |
|---|---|---|---|---|---|
| 11 | 191 | 18 | 0.05 | 47 | 0.02 |
| 12 | 382 | 1 | 0.00 | 10 | 0.00 |
| 13 | 367 | **50** | **0.62** | 27 | 0.02 |
| 15 | 145 | 59 | 0.05 | 32 | 0.01 |
| 23 | 180 | 29 | 0.01 | 1 | 0.00 |
| 27 | 5498 | 50 | 0.00 | 13 | 0.00 |
| 46 | 156 | 41 | 0.09 | 47 | 0.03 |
| 47 | 481 | 38 | 0.04 | 29 | 0.01 |
| 49 | 31 | 46 | 0.06 | 47 | 0.06 |
| 50 | 768 | **13** | **0.30** | 27 | 0.02 |
| 58 | 70 | 50 | 0.01 | 1 | 0.00 |
| 59 | 316 | 15 | 0.02 | 1 | 0.00 |
| 63 | 154 | 50 | 0.03 | 13 | 0.01 |

13 (Ledger, Edwards and Woodward, unpublished, in Brose et al. 2005) and 50
(Ledger et al. 2012) describe the same system. When each was left out alone,
the other predicted it: 640 of the 1189 predictions of class `exact` came
from there, 442 of them observed. Other pairs by the same authors share few
links under the names of the registry: 46 and 49 (Woodward, Hildrew et al.),
15 and 59 (P. H. Warren); 41 (Townsend and Hildrew 1979) holds 9% of those of
46.

Probability in a web with effect 0 (`p_average`) according to the groups:

| Class | Each web alone | 13 + 50 (bundled) | 13 + 50, 46 + 49, 15 + 59 |
|---|---|---|---|
| exact | 0.687 | 0.407 | 0.418 |
| genus | 0.496 | 0.514 | 0.530 |
| subfamily | 0.269 | 0.286 | 0.302 |
| family | 0.308 | 0.313 | 0.324 |
| order | 0.221 | 0.249 | 0.268 |
| class | 0.390 | 0.375 | 0.391 |
| generic | 0.696 | 0.715 | 0.730 |
| σ (logit) | 1.50 | 1.47 | 1.49 |

Example in the README (1000 iterations, seed 2026, unknown system):

| Groups | Links per iteration (median) | C control | Δ C (exposed − control) | P(Δ C > 0) | Δ mean TL | P(Δ TL > 0) |
|---|---|---|---|---|---|---|
| each web alone | 75.5 | 0.155 [0.051; 0.290] | 0.014 [−0.002; 0.035] | 0.95 | −0.074 [−0.200; 0.115] | 0.10 |
| 13 + 50 (bundled) | 64 | 0.143 [0.051; 0.288] | 0.014 [−0.002; 0.040] | 0.94 | −0.069 [−0.178; 0.111] | 0.11 |
| 13 + 50, 46 + 49, 15 + 59 | 65 | 0.145 [0.052; 0.289] | 0.015 [−0.002; 0.037] | 0.94 | −0.069 [−0.181; 0.112] | 0.11 |

## Unknown or typical system (`dev/validation/example_system.R`)

Example in the README with the bundled calibration, 1000 iterations:

| `system` | Links (median) | S control | C control | Δ S | Δ C, P(> 0) | Δ mean TL, P(> 0) |
|---|---|---|---|---|---|---|
| `"random"` (default) | 64 | 20.2 [6.0; 22.7] | 0.143 [0.051; 0.288] | −4.6 [−6.3; −0.7] | 0.014 [−0.002; 0.041], 0.94 | −0.069 [−0.178; 0.111], 0.11 |
| `"average"` | 66 | 22.3 [20.3; 22.7] | 0.128 [0.075; 0.187] | −5.0 [−5.7; −4.0] | 0.011 [−0.003; 0.026], 0.94 | −0.088 [−0.185; 0.020], 0.05 |

## Bundled calibration (13 webs; 13 and 50 left out together)

Precision of the links inferred with the rule of the reanalysis in the web
left out: 36% overall, from 2% (Jones 1951) to 84% (Harper-Smith, in Brose et
al. 2005). In a web with effect 0: exact 0.41; genus 0.51; subfamily 0.29;
family 0.31; order 0.25; class 0.38; generic groups 0.72. Genus beats exact
(difference 0.43 on the logit scale, standard error of the difference 0.12),
mainly because of webs 27 (exact 23/85, genus 375/561) and 50 (31/114 against
43/91); the model assumes the same class effect in every web. In the
exploratory prototype, without groups, the number of studies that support a
link did not improve the prediction (coefficient −0.06, standard error 0.05,
with lme4).

## Second phase: Bayesian calibration

| Check | Result |
|---|---|
| Posterior of a small problem (1 class, 3 webs) against the grid posterior | means of \(\alpha\) and \(\sigma\) within 4 Monte Carlo errors (test) |
| Split R-hat and ESS | ≈ 1 and ≈ 4000 in independent chains; R-hat > 1.1 with shifted chains; ESS < 600 with AR(0.9) (test) |
| `p_new` against `integrate()` | equal with tolerance 1e-8 (test) |
| Bundled calibration (13 webs, default rule, 4 × 4000) | maximum R-hat 1.001; minimum ESS 2982; σ = 1.58 |

Bundled calibration, probability in a new system (`p_new`): exact 0.46;
genus 0.53; subfamily 0.37; family 0.39; order 0.31; class 0.51; generic
0.67.

## Comparison of rules (`dev/validation/compare_rules.R`)

12 validation webs of Gray without Tuesday Lake (study 12), left out by
system (13 and 50 together); registry without study 12; 21 rules of
`candidate_rules()`; Bayesian refits (2 chains × 500) with seed 2026. About
30 min in two processes.

| Rule | Brier | Difference (se) | Skill | TSS | Precision | Recall |
|---|---|---|---|---|---|---|
| climb/family/class/5 | 0.2225 | — | −0.025 | 0.17 | 0.39 | 0.67 |
| climb/family/class/1 (reanalysis) | 0.2229 | 0.0004 (0.0024) | −0.022 | 0.18 | 0.39 | 0.65 |
| climb/family/class/1 without downward inheritance | 0.2236 | 0.0011 (0.0024) | −0.023 | 0.19 | 0.40 | 0.64 |
| climb/family/order/1 | 0.2276 | 0.0051 (0.0027) | −0.049 | 0.18 | 0.39 | 0.63 |
| fixed/order | 0.2284 | 0.0059 (0.0042) | −0.061 | 0.04 | 0.34 | 0.75 |
| fixed/family | 0.2338 | 0.0113 (0.0067) | −0.074 | 0.18 | 0.40 | 0.57 |
| climb/genus/class/1 | 0.2360 | 0.0134 (0.0073) | −0.063 | 0.21 | 0.44 | 0.52 |
| fixed/genus | 0.2562 | 0.0337 (0.0196) | −0.121 | 0.18 | 0.47 | 0.34 |
| ladder/class/1 (WebBuilder) | 0.2625 | 0.0400 (0.0215) | −0.141 | 0.16 | 0.43 | 0.42 |
| fixed/exact | 0.2964 | 0.0739 (0.0326) | −0.247 | 0.08 | 0.43 | 0.17 |

(Full table: `gray2015_rule_comparison`.) Conclusions: climbing from family
to class is best; the threshold (1 or 5) and the downward inheritance make no
distinguishable difference; starting at genus, the WebBuilder ladder and the
fixed level score worse, and not generalising scores worst. The skill against
the null model (mean connectance of the other systems) is slightly negative
for every rule: a third of the recorded links are not inferred and each one
costs 1.

Sensitivity (`dev/validation/scoring_sensitivity.R`):

- With maximum likelihood refits the ranking hardly changes (the best is
  still climb/family/class/5).
- If the pairs not inferred receive the background probability of the other
  systems (mean proportion of recorded links among their pairs not inferred,
  ≈ 0.19, almost the connectance), the Brier score falls to 0.208–0.211 and
  the 11 climb rules and fixed/family-order end up tied; the ladder and
  fixed/exact remain worse. In the webs of the registry, the fact that a pair
  is not inferred says little about the absence of a link.

## Final test on Tuesday Lake (`dev/validation/tuesday_lake.R`)

Calibration of each rule with the 12 webs; metaweb of the 82 taxa with the
registry without study 12; observed web of each year (only the taxa
present).

| Rule | Brier | Skill | TSS | Precision | Recall |
|---|---|---|---|---|---|
| climb/family/class/5 (default) | 0.171 | −0.007 | 0.15 | 0.44 | 0.22 |
| climb/family/class/1 | 0.173 | −0.019 | 0.15 | 0.45 | 0.22 |
| climb/genus/class/5 | 0.174 | −0.025 | 0.13 | 0.42 | 0.20 |
| climb/family/order/1 | 0.190 | −0.118 | 0.01 | 0.23 | 0.06 |
| ladder/class/1 | 0.179 | −0.056 | 0.11 | 0.44 | 0.16 |
| fixed/genus | 0.192 | −0.130 | 0.00 | 0 | 0 |

The default rule is also the best of the 21 on Tuesday Lake (Spearman
between the two rankings, 21 rules: 0.53). Only about a fifth of the
observed links are inferred: Gray (mostly British rivers) has few diets of
plankton; the rules that do not climb to class hardly link any rotifers or
cladocerans.

Indicators of 1000 sampled webs (default rule) against the observed web: S,
C and proportion of omnivores within the 95% interval in both years; L within
in 1984 and below in 1986; proportion of basal species too high (0.62
against 0.36 and 0.29) and mean trophic level too low (1.43 against 2.03 and
2.08): links are missing, and the sampling can only remove links from the
metaweb, not add the missing ones.

## Third phase: sources (Gray, GloBI and both)

### Download and processing of GloBI (`dev/globi/`)

- 500 queries: the consumers of the 13 validation webs, of Tuesday Lake and
  of the azoxy experiment, by name, genus and family; 832,003 `eats` records
  with citations (2026-11-03, on the maintainer's computer), without
  failures; one truncated query (Cyprinidae, 100,100 records).
- `globi_to_registry()`: removes 59,737 records of the copies of Gray and
  Brose in GloBI and 39,090 of taxa other than the one queried; 245,360
  records of 17,579 studies remain.
- Classification: 20% of the GloBI records with an order present in Gray had
  a different class (Teleostei and Actinopteri instead of Actinopterygii,
  Petromyzonti, Lepidosauria…); those synonyms were added to
  `default_taxonomy_fixes()`, which is now also applied to the nodes. The
  differences in family by genus (2%) and in order by family (3%;
  Centrarchiformes, Cottidae…) are matters of classification and are left as
  they are.
- Tuesday Lake in GloBI: 6 studies (908 records) with at least 10% of the
  links of TL84 or TL86 or that cite their authors (copies of Brose et al.
  2005, Jonsson et al. 2005 and Cohen et al. 2003): left out of the registry
  for the selection and the final test (`dev/globi/globi_tuesday_lake.csv`).
- Exclusions by web (`dev/globi/globi_exclusions.csv`): 43 web–study pairs
  of GloBI with at least 10% of the links of the web (51 with 5%, 31 with
  20%, 17 with 50%): copies of Brose et al. 2005 (up to 100% of web 46),
  Layer et al. 2010, Gray et al. 2015 / Thompson et al. 2017, “Globi”,
  “Gray” and the products derived from GloBI.
- After the exclusions, GloBI contains, with the same names, 4–18% of the
  links of each validation web (Gray, 0.4–13%).

### Cross-validation (`dev/validation/compare_sources.R`)

12 validation webs without Tuesday Lake, 13 and 50 together, exclusions by
web; 8 main rules with each source; Bayesian refits (2 × 500) with seed 2026.
About 80 min in two processes.

| Source | Use of the registries | Rule | Brier | Diff. from Gray (se) | Skill | TSS | Precision | Recall |
|---|---|---|---|---|---|---|---|---|
| Gray + GloBI | every registry at every level (**default**) | climb/family/class/5 | 0.1988 | −0.024 (0.013) | 0.058 | 0.28 | 0.39 | 0.91 |
| Gray + GloBI | every registry at every level | climb/family/class/1 | 0.2001 | −0.023 (0.014) | 0.054 | 0.28 | 0.40 | 0.91 |
| GloBI | — | climb/family/class/1 | 0.2003 | −0.023 (0.014) | 0.054 | 0.28 | 0.40 | 0.90 |
| GloBI | — | climb/family/class/5 | 0.2010 | −0.022 (0.013) | 0.046 | 0.27 | 0.39 | 0.91 |
| Gray + GloBI | Gray primary; GloBI at every level | climb/family/class/5 | 0.2034 | −0.019 (0.010) | 0.043 | 0.25 | 0.39 | 0.86 |
| Gray + GloBI | Gray primary; GloBI up to family | climb/family/class/5 | 0.2072 | −0.015 (0.010) | 0.024 | 0.23 | 0.39 | 0.83 |
| Gray | — | climb/family/class/5 | 0.2225 | — | −0.025 | 0.17 | 0.39 | 0.67 |

(Full table: `source_comparison`.) With any source, climbing from family to
class (or to order, or fixed family) is best; starting at genus, the ladder
and matching without generalisation score worse. GloBI improves the Brier
score of every rule (except the ladder when Gray is the primary registry);
the improvement of the default rule is 1.9 standard errors (9 of the 12 webs
improve, mainly 27, 11, 50 and 15; 47, 58 and 49 get somewhat worse). Among
the configurations with GloBI the differences are within the error (GloBI at
every level against only up to family, with Gray primary: −0.0039, se 0.0018,
8 of 12 webs). The precision does not change (0.39): GloBI raises the recall.

Sensitivity (`dev/validation/source_sensitivity.R`):

- Without the secondary compilations of GloBI (23 studies, 135,946 records:
  trophiCH, Ho 2022, products and copies of GloBI and Gray, Tetra-EU and
  trait databases such as freshwaterecology.info), which could contain
  generalised versions of the validation webs: default 0.199 (skill 0.062);
  GloBI alone 0.203; Gray primary 0.206–0.207. The improvement does not
  depend on them.
- Exclusion threshold of 5%: GloBI alone 0.203; Gray primary with GloBI up
  to family 0.210.

### Final test on Tuesday Lake

| Source | Rule | Brier | Skill | TSS | Precision | Recall |
|---|---|---|---|---|---|---|
| Gray + GloBI, every registry at every level | climb/family/class/5 (default) | 0.153 | 0.099 | 0.34 | 0.31 | 0.71 |
| GloBI | climb/family/class/1 | 0.150 | 0.116 | 0.35 | 0.33 | 0.67 |
| GloBI | fixed/family | 0.148 | 0.127 | 0.36 | 0.36 | 0.62 |
| Gray primary, GloBI at every level | climb/family/class/5 | 0.175 | −0.031 | 0.34 | 0.31 | 0.71 |
| Gray primary, GloBI up to family | climb/family/class/5 | 0.164 | 0.037 | 0.36 | 0.36 | 0.62 |
| Gray | climb/family/class/5 | 0.171 | −0.007 | 0.15 | 0.44 | 0.22 |

Spearman between the rankings of the cross-validation and of Tuesday Lake
(37 combinations): 0.73. With GloBI, two thirds of the observed links are
inferred (a fifth with Gray). With Gray as the primary registry and GloBI at
every level, the links of Tuesday Lake are the same as with GloBI alone, but
the Brier score is worse: the class `generic` (“Algae” extended to all algae)
is calibrated at 0.61 with Gray (frequent in the validation rivers) and at
0.31 with GloBI, and on Tuesday Lake only 13% of those 182 links are
observed. The transfer of the probability of a class between types of system
is the weak point.

Indicators of 1000 sampled webs (default, Gray + GloBI) against the observed
web: C, L, proportion of omnivores, of top predators and trophic similarity
within the 95% interval in both years; S (trophic species without isolated
nodes) too high (50 against 22 and 21: the observed web aggregates many
producers with the same consumers), proportion of basal species too high
(0.53–0.57 against 0.36 and 0.29) and mean trophic level too low (1.55–1.62
against 2.03 and 2.08), although closer than with Gray alone (0.62 and
1.43).

### Bundled calibrations (13 webs)

| | `gray2015_calibration` | `gray2015_globi_calibration` |
|---|---|---|
| Inferred links (observed) | 14,497 (36%) | 20,214 (40%) |
| σ (logit) | 1.58 | 1.23 |
| p in a new system: exact, genus, subfamily, family, order, class, generic | 0.46; 0.53; 0.37; 0.39; 0.31; 0.51; 0.67 | 0.45; 0.55; 0.36; 0.33; 0.52; 0.42; 0.33 |
| Maximum R-hat; minimum ESS | 1.001; 2982 | 1.001; 3520 |

### Matching engine with large registries

| Check | Result |
|---|---|
| Metaweb of the reanalysis against the original code (`dev/reference/compare_metaweb.R`) after each change of the engine | identical (1319 links, fields, adjacency and levels) |
| Counts of the validation of the bundled calibration | identical |
| Prefilter of records (`relevant_rows()`): with `build_metaweb()`, the records that match no node in any role are removed | metaweb identical in every field (test) |
| In the validation, the records that do not match in both roles are removed too | same links when the levels converge (test with web 49) |
| Validation of the Gray webs (12–13 webs) | 2.4 min → 37 s |
| Web 27 (221 nodes) with Gray + GloBI | 1 min, 1.5 GB of peak memory (before, more than 3 GB: two processes did not fit) |

## Sampling

| Check | Result |
|---|---|
| Frequency of each link in 4000 iterations against its probability | within 4 Monte Carlo errors (test) |
| Seed | reproducible and without altering the global random number stream (test) |

The sections “Webs of the same system”, “Unknown or typical system” and
“Bundled calibration (13 webs; 13 and 50 left out together)” are from the
first phase (maximum likelihood fit and rule with `min_partners = 1`); the
others, from the second (default rule with `min_partners = 5`, Bayesian fit)
and from the third (sources; by default, every registry at every level).

Scripts: `dev/validation/calibration_lme4.R`, `dev/validation/web_coverage.R`,
`dev/validation/validation_groups.R`, `dev/validation/example_system.R`,
`dev/validation/compare_rules.R`, `dev/validation/scoring_sensitivity.R`,
`dev/validation/tuesday_lake.R`, `dev/globi/process_globi.R`,
`dev/validation/compare_sources.R`, `dev/validation/source_sensitivity.R`
(order and run times in `dev/validation/README.md`). The comparisons with
the original code (`dev/reference/compare_metaweb.R`,
`dev/reference/compare_indicators.R`) need a local copy of the original
project, which is not distributed.
