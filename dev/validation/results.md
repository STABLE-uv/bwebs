# Validation — *bwebs* 0.2.0

<div align="justify">

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
| Maximum likelihood fit (`method = "ml"`) against `glmer(nAGQ = 20)` of *lme4* (Bates et al. 2015) on the counts of the bundled calibration (default rule, 13 webs; `dev/validation/calibration_lme4.R`) | coefficients equal to within 3e-6; standard errors to within 7e-6; σ 1.399165 against 1.399169 |
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
with *lme4*).

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
1.43). Among taxa (without the aggregation into trophic species) the
differences are much smaller: see “Fourth phase”.

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

## Fourth phase: taxonomic distance decay (experimental; `dev/validation/decay_background.R`)

Every pair gets an evidence profile (`decay_profiles()`, `R/decay.R`): for
each level of relatedness on the consumer side (the consumer itself, genus,
subfamily, family, order, class, phylum or any consumer), the closest
relation on the resource side among the records of those consumers, and
whether they eat a generic group that applies to the resource. The phylum of
each class comes from the taxonomic paths of GloBI, with the names of the
GBIF backbone (internal table `class_phylum`, `data-raw/class_phylum.R`). The
probability of a profile: logit p = the best of its evidence, a[consumer
level] + b[resource level] (d for a generic group instead of b), plus the
web effect, with effects that never increase with the distance. Profiles of
the 12 webs: 3 s with Gray, 29 s with Gray + GloBI.

### Cross-validation (12 webs; difference from the default rule web by web)

| Registries | Use | Variant | Fit | Brier | Default | Difference (se) |
|---|---|---|---|---|---|---|
| Gray | instead of the rule | best evidence of any consumer level | ML | 0.247 | 0.2225 | +0.025 (0.033) |
| Gray | instead of the rule | best evidence, discounted by the resources known for closer consumers | ML | 0.248 | 0.2225 | +0.025 (0.033) |
| Gray | instead of the rule | consumer resolved at its closest level with records, from the family | ML | 0.220 | 0.2225 | −0.003 (0.010) |
| Gray | background of the rule | consumer resolved, from the family | Bayesian | 0.206 | 0.2225 | −0.017 (0.006) |
| Gray + GloBI | instead of the rule | best evidence of any consumer level | ML | 0.242 | 0.1988 | +0.043 (0.028) |
| Gray + GloBI | instead of the rule | best evidence, discounted | ML | 0.238 | 0.1988 | +0.039 (0.027) |
| Gray + GloBI | instead of the rule | consumer resolved, from the family | ML | 0.226 | 0.1988 | +0.027 (0.018) |
| Gray + GloBI | background of the rule | consumer resolved, from the family | Bayesian | 0.197 | 0.1988 | −0.0021 (0.0010) |

- Instead of the rule, the decay is worse: the evidence of distant relatives
  reaches every consumer (some insect of the registry eats “Diatom”), and the
  webs where few consumers eat diatoms are over-predicted (web 23 with Gray:
  1100 links expected against 179 recorded). The discount for the resources
  known for closer consumers improves the fit (log likelihood +376 with Gray)
  but not the scores.
- As the background of the rule, it improves the scores: the inferred links
  keep the probability of their class and the pairs not inferred, which had
  probability 0, take that of their profile. With Gray, the score improves in
  10 of the 12 webs; with Gray + GloBI the gain is small, because the rule
  already infers 91% of the recorded links (477 recorded links among the
  8,814 pairs not inferred).

### Tuesday Lake (Gray + GloBI; background fitted to the 12 webs)

| Year | Brier: rule; rule + background | Skill | Missed links (expected by the background) | Links: observed; rule; rule + background |
|---|---|---|---|---|
| 1984 | 0.147; 0.142 | 0.137; 0.168 | 78 (63) | 264; 228; 291 |
| 1986 | 0.160; 0.154 | 0.061; 0.093 | 69 (52) | 236; 242; 294 |

Background (posterior means, logit scale): consumer level, from the consumer
itself to any consumer, −0.10, −0.38, −0.85, −1.08, −1.13, −1.35, −1.85,
−2.16; resource level, from the resource itself to unrelated, 0, −0.36,
−0.68, −0.98, −1.22, −1.26, −1.27, −1.29; generic group +1.03; σ 1.15. The
resource effects stop decreasing beyond the order: pairs that are not
related up to the phylum keep a probability of about 0.08 (in the validation
webs, 5% of the pairs not inferred are recorded links).

Indicators of 1000 sampled webs (independent links with their probability
in a new system), observed and mean of the sampled webs:

| Indicator | 1984: observed; rule; rule + background | 1986: observed; rule; rule + background |
|---|---|---|
| C (taxa) | 0.106; 0.076; 0.093 | 0.091; 0.078; 0.092 |
| Mean trophic level (taxa) | 1.68; 1.61; 1.70 | 1.56; 1.51; 1.55 |
| Proportion of basal taxa | 0.50; 0.55; 0.55 | 0.57; 0.61; 0.61 |
| Proportion of top taxa | 0.020; 0.032; 0.005 | 0.039; 0.029; 0.006 |
| S (trophic species) | 22; 55; 56 | 21; 55; 57 |
| Proportion of basal trophic species | 0.36; 0.55; 0.55 | 0.29; 0.60; 0.61 |
| Mean trophic level (trophic species) | 2.03; 1.61; 1.70 | 2.08; 1.51; 1.55 |
| Proportion of omnivores (trophic species) | 0.50; 0.62; 0.94 | 0.47; 0.63; 0.93 |

- Among taxa the sampled webs are close to the observed ones, and with the
  background their connectance and mean trophic level reach the observed
  values. The excess of basal species and the low trophic level found with
  `web_indicators()` come mostly from the aggregation into trophic species:
  the observed webs collapse into 21 and 22 trophic species (many producers
  share their consumers), the sampled webs do not (their links are drawn
  independently), so their indicators of trophic species are not comparable.
- The background spreads its probability over many unlikely pairs: almost
  every consumer becomes an omnivore and almost no taxon remains a top
  predator. The taxonomy of the resource does not separate plausible from
  implausible pairs well; the predator/prey size ratio is meant for that.

## Fifth phase: the taxonomic decay in the package (version 0.3.0; `dev/validation/decay_package.R`)

`matching_rule(decay = TRUE)`, the default since 0.3.0, gives the pairs
without an inferred link the probability of their evidence (`R/decay.R`):
the consumer resolved at its closest level with recorded diets, from the
family (family, order, class, phylum, any consumer), and the closest level
at which the resource is related to that diet (within the family of a
recorded resource, order, class, phylum, not related), plus an additive
effect when a record of a generic group applies and an effect of the size
of the web; logit p = a[consumer] + b[resource] + γ·generic +
β·log(nodes / 50) + v[web], with steps ≥ 0 so that closer evidence never
gives a lower probability. Changes from the exploration (fourth phase): the
resource levels within the family share one effect (with Gray only 129 of
the 14,102 pairs not inferred are related within the family: the rule
infers almost all of them); the steps are on the natural scale with
half-normal priors; the generic effect is added instead of competing with
the resource level (the maximum left a flat region); the MCMC reflects the
steps at 0 and makes 20 random-walk steps per draw kept (with 10, the R-hat
of the σ of the web effect was 1.014 in the bundled Gray calibration). With the model of
the exploration the chains did not mix (R-hat up to 1.24, effective sample
size down to 12); in the bundled calibrations the maximum R-hat of the decay
is 1.004 (Gray) and 1.003 (Gray + GloBI), with effective sample sizes above
4800. The class model is unchanged (the bundled class draws are identical).

### Cross-validation (12 webs; the same validation scored without and with the decay)

| Registries | Brier: rule; rule + decay | Difference (se) | Webs better | Skill: rule; rule + decay | Pairs not inferred | Recorded among them | Expected by the decay |
|---|---|---|---|---|---|---|---|
| Gray | 0.2225; 0.2093 | −0.0132 (0.0085) | 7 of 12 | −0.025; 0.015 | 14,102 | 3,248 | 4,746 |
| Gray + GloBI | 0.1988; 0.1971 | −0.0017 (0.0018) | 6 of 12 | 0.058; 0.066 | 8,814 | 477 | 846 |
| Gray + GATEWAy (sixth phase) | 0.2142; 0.2052 | −0.0091 (0.0046) | 9 of 12 | 0.004; 0.034 | 11,570 | 1,421 | 1,601 |

- With Gray the decay improves the score as in the exploration (0.206
  there). Without the size of the web its expected number of links among
  the pairs not inferred matched the recorded one (3,369 and 3,248); with it
  the decay expects more in the largest webs left out (27: 2,486 expected,
  2,114 recorded; 47: 888 and 109), whose size the other webs do not cover.
- With Gray + GloBI the gain is small: the rule already infers 91% of the
  recorded links. Without the size of the web the decay over-predicted the
  largest webs, where few pairs not inferred are recorded links (web 27,
  3,623 pairs not inferred: 731 expected, 55 recorded; web 13: 315 and 64),
  while it helped the small predator webs (15, 59) and 50; the web effect of
  the decay absorbed this in the fit (−1.83 for web 27, +1.38 for web 59) but
  not for a web left out. With the size of the web, 289 and 168 links are
  expected in webs 27 and 13.
- The size of the web explains part of it, as expected from the decrease of
  connectance with the number of species: with log(number of nodes) as a
  covariate (`dev/validation/decay_role.R`: `glmer()` of *lme4*, leave one system
  out, probabilities in a new system) the coefficient is −0.62 with Gray +
  GloBI, the standard deviation of the web effect falls from 0.98 to 0.73
  and the squared error over the pairs not inferred from 473 to 445 (477
  with probability 0); with Gray it has no effect (+0.10).
- Added to the decay model of the package (effect of log(nodes / 50); in
  `sample_metawebs()` the number of nodes of the metaweb, in
  `score_metaweb()` the nodes evaluated), because it gave the best
  validation by folds with the default sources
  (`dev/validation/decay_web_size.R`, the same validation scored with the
  decay without and with the size of the web): Gray + GloBI 0.1983 → 0.1971
  (difference −0.0012, se 0.0019; 6 of 12 webs better; links expected among
  the pairs not inferred 1,536 → 846, against 477 recorded), Gray 0.2077 →
  0.2093 (+0.0016, se 0.0027). The tables of this section are with the size
  of the web.

### Tuesday Lake (calibration with the 12 webs; registries without Tuesday Lake)

The scores of each year with the rule alone, with the decay, and with the
decay and the body masses of that year (`tuesday_lake_community$mass_kg`;
the size term of `size_model`, fitted to the aquatic webs of GATEWAy, sixth
phase):

| Registries | Year | Brier: rule; + decay; + masses | Skill: rule; + decay; + masses | Observed links missed by the rule (expected by the decay; with masses) |
|---|---|---|---|---|
| Gray | 1984 | 0.1685; 0.1433; 0.1413 | 0.009; 0.157; 0.169 | 198 (159; 131) |
| Gray | 1986 | 0.1737; 0.1434; 0.1426 | −0.022; 0.156; 0.161 | 190 (136; 117) |
| Gray + GloBI | 1984 | 0.1467; 0.1413; 0.1384 | 0.137; 0.169; 0.186 | 78 (68; 43) |
| Gray + GloBI | 1986 | 0.1595; 0.1539; 0.1506 | 0.061; 0.094; 0.113 | 69 (55; 37) |
| Gray + GATEWAy | 1984 | 0.1581; 0.1490; 0.1393 | 0.070; 0.124; 0.181 | 111 (153; 113) |
| Gray + GATEWAy | 1986 | 0.1539; 0.1454; 0.1368 | 0.094; 0.144; 0.194 | 101 (143; 114) |

Indicators of 1000 webs sampled with `sample_metawebs()` (default: system
effect drawn in each draw), observed and mean of the sampled webs, rule
alone, rule + decay and rule + decay with the masses of the year (among
taxa: without isolated taxa and without aggregation into trophic species):

| Indicator | Gray, 1984 | Gray + GloBI, 1984 | Gray + GloBI, 1986 | Gray + GATEWAy, 1984 |
|---|---|---|---|---|
| C (taxa) | 0.106; 0.074; 0.080; 0.069 | 0.106; 0.080; 0.093; 0.086 | 0.091; 0.083; 0.092; 0.086 | 0.106; 0.079; 0.113; 0.099 |
| Mean trophic level (taxa) | 1.68; 1.32; 1.66; 1.52 | 1.68; 1.60; 1.70; 1.60 | 1.56; 1.52; 1.57; 1.51 | 1.68; 1.81; 1.82; 1.68 |
| Proportion of basal taxa | 0.50; 0.71; 0.58; 0.59 | 0.50; 0.55; 0.55; 0.55 | 0.57; 0.60; 0.60; 0.60 | 0.50; 0.53; 0.55; 0.55 |
| Proportion of top taxa | 0.020; 0.178; 0.035; 0.085 | 0.020; 0.077; 0.034; 0.070 | 0.039; 0.063; 0.028; 0.062 | 0.020; 0.101; 0.016; 0.062 |
| Proportion of omnivores (trophic species) | 0.50; 0.49; 0.76; 0.58 | 0.50; 0.58; 0.84; 0.69 | 0.47; 0.58; 0.83; 0.68 | 0.50; 0.69; 0.95; 0.86 |

- On Tuesday Lake the decay improves every score, with every registry.
  Among taxa, the sampled webs come close to the observed ones in
  connectance, trophic level and top taxa; with Gray alone, which has few
  diets of plankton, the rule leaves 47 taxa without links and the decay
  fills them. With Gray + GloBI it expects about as many links among the
  pairs not inferred as were missed (68 and 55, against 78 and 69); with
  Gray + GATEWAy more (153 and 143, against 111 and 101).
- The decay spreads probability over pairs that size would exclude
  (rotifers eating fish), so most consumers become omnivores. The draws are
  also much more variable in their number of links (95% of draws between 54
  and 568 links among trophic species with Gray + GloBI in 1984), because
  the system effect of the decay is drawn in every draw.
- With the body masses every score improves, most with Gray + GATEWAy
  (0.1490 → 0.1393, 0.1454 → 0.1368), and the omnivores fall (0.84 → 0.69
  with Gray + GloBI, 0.95 → 0.86 with Gray + GATEWAy, 0.76 → 0.58 with
  Gray; 0.47–0.50 observed), but the top taxa rise above the observed ones
  (0.06–0.09 against 0.02–0.04) and the mean trophic level falls slightly.
  With the first size term (curve of Brose et al. 2005 and *cheddar*, optimum
  120 times, weights 0.5 and 10) the omnivores fell further (0.61 with Gray
  + GloBI) but the links of fish to zooplankton (median log10 mass ratio
  4.7 in 1984, against 2.1 for the invertebrate predators) were removed,
  the top taxa rose to 0.09 and the Brier score of Gray got worse (0.1433 →
  0.1468).

### Role of the resource (`dev/validation/decay_role.R`)

Recorded links among the pairs not inferred of the 12 webs, by the role of
the resource and whether the rule infers links of the consumer to resources
of that kind (basal or animal):

| Resource | The consumer has inferred links to that kind | Gray | Gray + GloBI |
|---|---|---|---|
| Animal | yes | 4.4% of 6,664 pairs | 1.8% of 5,887 |
| Animal | no, only to basal resources | 0.4% of 1,977 | 0.1% of 1,370 |
| Producer | yes | 65% of 3,727 | 18% of 824 |
| Producer | no, only to animals | 39% of 919 | 30% of 208 |
| Detritus | yes | 36% of 240 | 42% of 300 |
| Detritus | no, only to animals | 17% of 266 | 14% of 225 |
| Any | the consumer has no inferred link | 17% of 309 | none |

Logistic models with a random effect of the web, each system left out in
turn and predicted with its probability in a new system; squared error over
the pairs not inferred (3,248 with Gray and 477 with Gray + GloBI when they
all have probability 0):

| Model | Gray | Gray + GloBI |
|---|---|---|
| Taxonomy (the evidence of the decay of the package) | 2,651 | 473 |
| Taxonomy + role of the resource | 1,932 | 547 |
| Taxonomy + role + inferred diet of the consumer | 1,849 | 576 |
| Role + inferred diet only | 1,874 | 540 |

- With Gray alone the role of the resource predicts the missed links much
  better than taxonomy: the rule misses most of the algae eaten by grazers.
  With Gray + GloBI it makes the webs left out worse than probability 0: the
  missed links are few and concentrated in a few webs (in the predator webs
  15, 46 and 59 a sixth of the animal pairs are recorded links, in the
  whole-community webs 13, 27 and 47 almost none), so the effects do not
  transfer from web to web. Not used in the package.

## Size term (version 0.3.0; `dev/validation/size_brose.R`, `size_folds.R`, `size_lambda.R`, `size_sources.R`)

Webs with the body masses of consumers and resources: those of Brose et al.
(2005; `bodysizes_2008.txt`: Grand Cariçaie, Celtic Sea, Benguela,
Broadstone, the Frome, Tuesday Lake) and of the *cheddar* package (Hudson et
al. 2013; Broadstone
Stream, the two Mill Stream channels of Ledger et al. 2012, the ten streams
of Layer et al. 2010 with masses, Tuesday Lake 1984 and 1986). Size term:
within-consumer effect of x = log10(M_consumer / M_resource) on the links
of each predator with the animals of its web (logistic, consumer fixed
effects, quadratic in x), relative to the optimum of each consumer, so that
g(x) ≤ 0; added on the logit scale to the predation pairs (animal
resources) whose two nodes have a mass, with a weight for the links
inferred by the rule and another for the pairs of the decay. Fitted without
the webs of the system left out, the curve hardly changes: optimum x
between 1.7 and 2.1 (a consumer 50 to 120 times heavier than its prey),
quadratic coefficient between −0.20 and −0.16. Without Grand Cariçaie,
which holds 84,456 of the 90,376 pairs, the curve is flatter and its
optimum much larger (0.89 x − 0.10 x²: optimum 4.3, a consumer 21,000 times
heavier than its prey), set by the marine webs of fish.

Folds (`size_folds.R`): every validation web whose consumers have masses
(by name or genus, from any of those webs) for at least 60% of them is left
out by system, as in `score_rule()`; its pairs take the probabilities of
the class and decay models fitted to the other systems (the decay with the
size of the web), and the size term is fitted to the webs with masses of
the other systems (never to Tuesday Lake, scored with the fits of the 12
webs). Seven webs qualify: 13 and 50 (the Frome), 46 and 49 (Broadstone),
27 (the streams of Layer et al. 2010), 47 and 23; 11, 15, 58, 59 and 63
have masses for 18–53% of their consumers. Difference in the Brier score
from no size term (mean of the 7 webs, its standard error and the webs
that improve), and Tuesday Lake, with the weight on the pairs of the decay
fixed at 1 (`size_folds.R`, `size_lambda.R`):

| Size term applied to | Gray | Gray + GloBI | Tuesday Lake (Gray + GloBI; 0.1402 without) |
|---|---|---|---|
| the pairs of the decay | −0.0068 (0.0022), 7 of 7 | −0.0014 (0.0009), 7 of 7 | 0.1386 |
| the pairs of the decay, keeping their expected number of links | −0.0029 (0.0026), 5 of 7 | +0.0005 (0.0005), 4 of 7 | 0.1414 |
| the pairs of the decay, and ¼ of it to the inferred links | −0.0096 (0.0041), 6 of 7 | −0.0051 (0.0031), 6 of 7 | 0.1360 |
| the pairs of the decay, and ½ of it to the inferred links | −0.0107 (0.0065), 6 of 7 | −0.0067 (0.0057), 6 of 7 | 0.1363 |
| every predation pair | −0.0102 (0.0112), 6 of 7 | −0.0071 (0.0106), 6 of 7 | 0.1374 |

- On the pairs of the decay the size term helps in every web, with both
  registries and on Tuesday Lake: it lowers the probability of the pairs
  that size excludes, which the decay over-predicts. Keeping the expected
  number of links of the decay (moving that probability to plausible pairs)
  loses the gain.
- On the inferred links it helps a lot where the rule infers predation
  links between taxa of implausible sizes (Gray + GloBI, full effect: web 13
  0.121 → 0.103, 47 0.155 → 0.131, 49 0.200 → 0.163, 50 0.185 → 0.170), but
  in Broadstone 46, a predator web where 64% of the predator–prey pairs are
  links, it lowers the probability of many recorded links (0.214 → 0.265).
  A one-sided effect (only prey heavier than the optimum penalised) does not
  avoid it (0.263).
- Letting the optimum depend on the mass of the consumer (Gravel et al.
  2013) changes little: the mean differences are within 0.001 of those of
  the quadratic term.

Weights and curve (`size_lambda.R`, `size_sources.R`): the pairs of the
seven webs scored again with a grid of weights, 0 to 1 on the inferred
links and 1, 3, 10 or 30 on the pairs of the decay, and with the curve
fitted to every web with masses of the other systems or also without Grand
Cariçaie. The configuration with the best mean difference over the seven
webs is the choice of the validation by folds; the same choice made
without each web in turn and scored on that web (nested) estimates how the
chosen configuration does on a new web. Both registries together (mean of
each web over Gray and Gray + GloBI):

| Choice among | Chosen (curve; weight on inferred links, on the decay) | Mean difference over the 7 webs | Nested: mean (se), webs better |
|---|---|---|---|
| weights up to 1, curve of every web | every web; 0.75, 10 | −0.0100 | −0.0050 (0.0100), 6 of 7 |
| weights up to ½, curve of every web | every web; 0.5, 10 | −0.0097 | −0.0089 (0.0062), 6 of 7 |
| weights up to ½, curve without Grand Cariçaie | without; 0.25, 1 | −0.0098 | −0.0031 (0.0141), 6 of 7 |
| weights up to ½, either curve | without; 0.25, 1 | −0.0098 | −0.0009 (0.0137), 6 of 7 |

- The best configurations are tied within 0.0003 of each other, but only
  with the curve of every web and weights on the inferred links up to ½ is
  the choice stable: chosen without each web it does almost as well on that
  web (−0.0089). Larger weights, or the curve without Grand Cariçaie, win
  in the six other webs and lose heavily on web 46 when it is left out. The
  same holds with each registry (nested, with the curve of every web and
  weights up to ½: Gray −0.0097, Gray + GloBI −0.0072; their choices 0.5
  and 10, and 0.5 and 30, within 0.0001 of 0.5 and 10).
- First default (replaced in the sixth phase, with the masses and webs of
  GATEWAy): the curve of the 18 webs with masses other than Tuesday Lake (0.668 x − 0.161 x²; optimum 2.08, a
  consumer 119 times heavier than its prey), weight 0.5 on the inferred
  links and 10 on the pairs of the decay. With it, Gray + GloBI: web 13
  0.121 → 0.108, 23 0.187 → 0.185, 27 0.186 → 0.183, 46 0.214 → 0.237, 47
  0.155 → 0.139, 49 0.200 → 0.173, 50 0.185 → 0.174 (mean −0.0071, se
  0.0059); Gray: mean −0.0122 (se 0.0075), 6 of 7. Tuesday Lake (fits of
  the 12 webs): Gray + GloBI 0.1402 → 0.1361, Gray 0.1395 → 0.1388.

Kind of predator (`size_consumer_type.R`): on Tuesday Lake the curve fits
the invertebrate predators (median log10 mass ratio of their recorded prey
2.1) but not the fish (4.7 in 1984; see the fifth phase), so a size term
only for invertebrate predators, and one curve for invertebrates and another
for vertebrates (each fitted to the predators of that kind of the webs with
masses of the other systems), were scored in the same folds (both
registries together; fish are predators in webs 23, 27, 47 and Tuesday
Lake):

| Size term | Best weights | Mean difference over the 7 webs | Nested: mean (se) | Tuesday Lake (Gray + GloBI; 0.1402 without) |
|---|---|---|---|---|
| one curve (`size_model`) | 0.5, 10 | −0.0097 | −0.0089 (0.0062) | 0.1361 |
| invertebrate predators only | 0.5, 10 | −0.0095 | −0.0088 (0.0062) | 0.1375 |
| one curve per kind of predator | 0.5, 30 | −0.0092 | −0.0089 (0.0060) | 0.1372 |

- The curve of the vertebrate predators of the other webs (optimum 2.75, a
  fish 560 times heavier than its prey; invertebrates 1.49, 31 times) does
  not reach the mass ratios of the fish of Tuesday Lake, and neither variant
  improves the folds or Tuesday Lake: the single curve was kept (see the
  sixth phase for the curve of GATEWAy).

## Sixth phase: GATEWAy (version 0.3.0; `dev/gateway/process_gateway.R`, `dev/validation/compare_gateway.R`, `size_gateway.R`)

GATEWAy version 1.0 (Brose et al. 2019; iDiv Data Repository 283, CC BY
4.0): 290 food webs with the body masses of consumers and resources
(222,151 links; masses for 99.6% of the consumers and 91% of the
resources), 88 of them freshwater (60 lakes, 28 streams). Several are the
validation webs of Gray (Tuesday Lake, Broadstone, Mill Stream, the streams
of Layer et al. 2010, Skipwith Pond, the Icelandic streams), and 14,889 of
its records were already in the GloBI download (those of the queried
taxa). New: the 50 Adirondack lakes of Havens (1992; plankton and fish,
146 fish predators), the ten streams of the Taieri (Townsend et al. 1998)
and Lake Malawi, about 5,000 freshwater name pairs absent from Gray and
GloBI. The five chalk streams of Gray et al. (2015) and Thompson et al.
(2017) ("gut content, extrapolated from similar taxa") were built in part
with WebBuilder and the Gray registry.

### Registry (`dev/gateway/process_gateway.R`)

The lakes and streams without Tuesday Lake, the five chalk streams and
parasitic and parasitoid links: 19,185 records of 81 webs (one study per
web). Taxonomy of the names: from the Gray and GloBI registries by name
(582 names), by genus or higher rank (78), from the GBIF backbone (30:
exact matches, or misspellings within two letters, in the kingdom and
phylum expected from the metabolic type of GATEWAy), and unresolved (56:
detritus, groups such as "Phytoplankton", the local names of the fish of
Lake Malawi); 98.6% of the records have the class of the consumer and 92%
that of the resource. No study shares 10% of the links of Tuesday Lake.
Webs left out with each validation web: those holding at least 10% of its
links, alone or pooled by publication (the 20 streams of Layer et al. 2010
are web 27 of Gray, spread over ten webs of GATEWAy): 13, 27, 46, 47, 49
and 50 lose the ten webs of Layer et al. (2010), 59 Skipwith Pond and 63
the four Icelandic streams.

### Cross-validation and Tuesday Lake (default rule with the decay)

| Registries | Brier | Skill | Missed | Expected by the decay | Tuesday Lake 1984; 1986 | With the masses of each year |
|---|---|---|---|---|---|---|
| Gray | 0.2093 | 0.015 | 3,248 | 4,746 | 0.1433; 0.1434 | 0.1468; 0.1455 |
| Gray + GATEWAy | 0.2052 | 0.034 | 1,421 | 1,601 | 0.1490; 0.1454 | 0.1412; 0.1367 |
| Gray + GloBI | 0.1971 | 0.066 | 477 | 846 | 0.1413; 0.1539 | 0.1399; 0.1493 |
| Gray + GloBI + GATEWAy | 0.2013 | 0.038 | 443 | 1,115 | 0.1433; 0.1518 | 0.1410; 0.1460 |
| Gray + GloBI without the chalk streams | 0.2029 | 0.043 | 496 | 734 | 0.1456; 0.1582 | 0.1435; 0.1539 |
| Gray + GloBI + GATEWAy, without them | 0.2015 | 0.048 | 461 | 725 | 0.1481; 0.1563 | 0.1452; 0.1508 |

(Tuesday Lake with the masses: the size term of the first selection, curve
of Brose et al. 2005 and *cheddar* with weights 0.5 and 10.)

- With the Gray registry, GATEWAy improves the cross-validation (−0.0041,
  se 0.0030, 7 of 12 webs better), mostly where its webs have relatives of
  the validation taxa (web 47: 0.166 → 0.144; 27: 0.245 → 0.231; 11: 0.269
  → 0.245). Gray + GATEWAy is bundled (`gateway`,
  `gray2015_gateway_calibration`) and used in the examples [USER].
- Added to Gray + GloBI it makes the score worse (+0.0042, se 0.0025, 5 of
  12 webs better): it adds coarse diets (the Adirondack lakes have links
  between 59% of their predator–prey pairs) and the decay expects more
  links among the pairs not inferred (1,115 against 443 recorded). Gray +
  GloBI stays the best combination.
- Without the chalk streams the score of Gray + GloBI is worse (0.2029; 2 of
  12 webs better) and so is Tuesday Lake. Part of their gain may come from
  the validation webs themselves (they were built with the Gray registry);
  they are already left out with the webs that share at least 10% of their
  links (13, 23, 47, 50 and 63). Kept [USER].

### Size term with the masses of GATEWAy (`size_gateway.R`)

With the masses of GATEWAy, Brose et al. (2005) and *cheddar* (by name, then
genus), all 12 validation webs other than Tuesday Lake have masses for at
least 60% of their consumers (the fewest, web 11, 71%; seven webs in the
first selection). Curves fitted with
`size_curve_fit()` (one intercept per consumer profiled out; the estimates
of `glm()`) to the webs of the other systems: Brose et al. (2005) and
*cheddar* (as in the first selection), the freshwater webs of GATEWAy, its
aquatic webs (with the marine ones) or all its webs; one curve or one per
kind of predator (invertebrate, vertebrate); weights 0 to 1 on the inferred
links and 1 to 30 on the pairs of the decay. Curves fitted to all the webs
of each source but Tuesday Lake:

| Curve | Pairs (links) | x | x² | Optimum (consumer / prey) |
|---|---|---|---|---|
| Brose et al. (2005) and *cheddar* | 90,097 (7,965) | 0.66 | −0.162 | 110 |
| GATEWAy, freshwater | 17,112 (6,125) | 1.05 | −0.094 | 400,000 |
| GATEWAy, aquatic | 491,829 (47,467) | 0.50 | −0.059 | 17,000 |
| GATEWAy, every ecosystem | 1,093,988 (154,750) | 0.40 | −0.075 | 440 |

Choice of the configuration (mean difference over the 12 webs of the three
registries, Gray, Gray + GloBI and Gray + GATEWAy; nested: the same choice
made without each web and scored on it):

| Choice among | Chosen | Mean difference | Nested: mean (se), webs better |
|---|---|---|---|
| every configuration | GATEWAy aquatic, one curve; 0.5, 1 | −0.0068 | +0.0047 (0.0122), 9 of 12 |
| GATEWAy aquatic, one curve, weights on the inferred links up to 0.25 | 0.25, 1 | −0.0061 | −0.0061 (0.0045), 9 of 12 |
| GATEWAy aquatic, one per kind of predator, up to 0.5 | 0.5, 1 | −0.0058 | −0.0048 (0.0053), 9 of 12 |
| GATEWAy, every ecosystem, one curve, up to 0.5 | 0.5, 3 | −0.0054 | −0.0047 (0.0048), 9 of 12 |
| decay only (no weight on the inferred links), GATEWAy aquatic | 0, 1 | −0.0029 | −0.0029 (0.0014), 9 of 12 |
| Brose et al. (2005) and *cheddar*, up to 0.25 | 0.25, 1 | −0.0035 | −0.0028 (0.0041), 8 of 12 |
| the first default (Brose et al. and *cheddar*; 0.5, 10), fixed | | −0.0021 (se 0.0076), 8 of 12 | |

- Three webs whose recorded links are mostly predation (11, a pond with
  fish; 46, Broadstone; 59, Skipwith Pond: 96, 148 and 275 of their 186,
  149 and 316 links) get worse with any weight on the inferred links, the
  others better; with larger weights the choice made without one of them
  does badly on it.
- Default of the package (`size_model`, `data-raw/size_model.R`): the curve
  of the 206 aquatic webs of GATEWAy (49 lakes, 23 streams, 134 marine;
  0.503 x − 0.0593 x², optimum 4.24), weight 0.25 on the inferred links and
  1 on the pairs of the decay: −0.0061 over the 12 webs (Gray −0.0078, Gray
  + GloBI −0.0049, Gray + GATEWAy −0.0056; 9 of 12 webs better with each);
  Tuesday Lake 0.139 → 0.137 (Gray), 0.140 → 0.137 (Gray + GloBI), 0.149 →
  0.139 (Gray + GATEWAy). The curve is flat (−0.24 on the logit for a ratio
  100 times off the optimum, −1.1 for a consumer as heavy as its prey) and
  its optimum close to that of fish (Tuesday Lake: 4.7 in 1984), so it no
  longer removes the links of fish to zooplankton.

## Seventh phase: amount of evidence and ordered classes (version 0.3.0; `dev/validation/evidence_amount.R`)

Question of the maintainer: should the probability of an inferred link
depend, besides the taxonomic category of its evidence (the class: the
level of the closest record), on how much evidence supports it, so that a
link eaten by the taxon and its relatives, and recorded for the prey and
its relatives, is more probable than one recorded once for a distant
relative? And should a closer relative never give a lower probability?

Each link of the metaweb keeps its records and studies (`n_records`,
`n_sources`), the consumer and resource taxa among them (`n_consumer_taxa`,
`n_resource_taxa`) and the relatives of each node with records in the
network (`n_consumer_known`, `n_resource_known`); `validate_rule()` returns
them for every link inferred in the validation webs (`links`).

### Amount of evidence (folds by system; logistic model with the class and a web effect, maximum likelihood)

Difference in the Brier score of the inferred links of each web left out
(over its pairs; the part of the pairs not inferred is the same for every
model) from the class alone, and webs better, of 12:

| Covariates (with the class) | Gray | Gray + GloBI | Gray + GATEWAy |
|---|---|---|---|
| log studies | +0.0000 (6) | +0.0008 (4) | +0.0072 (6) |
| log records | +0.0003 (5) | +0.0015 (5) | +0.0023 (6) |
| log consumer taxa | +0.0003 (5) | +0.0010 (2) | +0.0026 (6) |
| log resource taxa | +0.0013 (4) | −0.0001 (5) | +0.0030 (6) |
| both | +0.0027 (3) | −0.0015 (7) | +0.0030 (6) |
| share of the relatives of the consumer | +0.0002 (7) | −0.0003 (4) | +0.0024 (8) |
| share of the relatives of the resource | +0.0030 (4) | +0.0004 (6) | +0.0021 (3) |
| both shares | +0.0032 (4) | +0.0004 (6) | +0.0042 (6) |
| studies and both taxa | +0.0041 (3) | −0.0013 (6) | +0.0074 (6) |
| all | +0.0063 (2) | −0.0007 (6) | +0.0106 (3) |

(standard errors of the differences 0.0003 to 0.005)

- No amount of evidence improves the prediction of the webs left out with
  every registry; the best with Gray + GloBI (consumer and resource taxa,
  −0.0015, se 0.0026) is worse with the other two.
- Fitted to all the webs, the effects are clear but change sign from one
  registry to another: log studies −0.32 with Gray, +0.03 with Gray +
  GloBI, +0.53 with Gray + GATEWAy; share of the relatives of the resource
  +0.69, −1.25 and +0.41. Within a registry the number of records of a link
  grows with how much its taxa were studied and with how coarse the levels
  at which its nodes were matched, not with how probable the link is in a
  new system.
- The probability of a link stays a function of its class only.

### Ordered classes (`calibrate_links(monotone = TRUE)`)

Coefficients of the classes exact, genus, subfamily, family, order, class
on the logit scale, never increasing in that order (non-negative steps,
half-normal priors; the reflecting random walk of the decay model, 10 steps
per draw kept); generic free. Probability in a new system with the 12 webs
(maximum likelihood):

| Class | Gray: free | ordered | Gray + GloBI: free | ordered | Gray + GATEWAy: free | ordered |
|---|---|---|---|---|---|---|
| exact | 0.432 | 0.494 | 0.416 | 0.482 | 0.438 | 0.507 |
| genus | 0.508 | 0.494 | 0.525 | 0.482 | 0.525 | 0.507 |
| subfamily | 0.338 | 0.363 | 0.333 | 0.340 | 0.356 | 0.386 |
| family | 0.362 | 0.363 | 0.295 | 0.340 | 0.388 | 0.386 |
| order | 0.282 | 0.309 | 0.508 | 0.340 | 0.309 | 0.327 |
| class | 0.471 | 0.309 | 0.402 | 0.340 | 0.462 | 0.327 |
| generic | 0.661 | 0.667 | 0.290 | 0.275 | 0.552 | 0.554 |

Cross-validation (Bayesian, `score_rule()`; the same fit of the decay for
both class models, so that its Monte Carlo error does not enter the
difference):

| Registries | Without the decay: free | ordered | difference (se), webs better | With the decay: free | ordered | difference (se), webs better |
|---|---|---|---|---|---|---|
| Gray | 0.2225 | 0.2245 | +0.0020 (0.0045), 8 | 0.2093 | 0.2104 | +0.0012 (0.0047), 8 |
| Gray + GloBI | 0.1988 | 0.1995 | +0.0007 (0.0033), 6 | 0.1971 | 0.1978 | +0.0007 (0.0032), 6 |
| Gray + GATEWAy | 0.2142 | 0.2138 | −0.0004 (0.0034), 7 | 0.2052 | 0.2044 | −0.0008 (0.0036), 7 |

- The two class models score the same; the ordered one makes 6 to 8 of the
  12 webs better and a few much worse (with Gray, Broadstone, web 46,
  +0.052: when it is left out the ordered fit lowers every class, and its
  links are recorded more often than in the other webs).
- With 2 steps per draw kept (the setting of the decay in `score_rule()`)
  the chains of the ordered classes had not mixed in some folds (with Gray,
  web 46: 0.42 instead of 0.46 for the exact class) and the scores were
  worse; `score_rule()` uses 10, as `calibrate_links()`.
- Default: the ordered classes, since the folds do not separate the two
  models and the maintainer asked that a more distant relative never give a
  higher probability [USER, confirmed on 2026-11-09].

Bundled calibrations (13 webs; `data-raw/`), probability in a new system
with the ordered classes (free classes in brackets):

| Class | `gray2015_calibration` | `gray2015_globi_calibration` | `gray2015_gateway_calibration` |
|---|---|---|---|
| exact | 0.493 | 0.497 | 0.510 (0.478) |
| genus | 0.488 | 0.495 | 0.506 (0.552) |
| subfamily | 0.367 | 0.360 | 0.383 (0.379) |
| family | 0.360 | 0.351 | 0.377 (0.411) |
| order | 0.313 | 0.351 | 0.309 (0.324) |
| class | 0.312 | 0.345 | 0.307 (0.434) |
| generic | 0.656 | 0.302 | 0.535 (0.561) |

### Tuesday Lake (calibration with the 12 webs; registries without Tuesday Lake)

Brier score, free → ordered classes (skill with the decay and the masses in
the last column):

| Registries | Year | Without the decay | With the decay | With the decay and the masses | Skill (masses) |
|---|---|---|---|---|---|
| Gray | 1984 | 0.1685 → 0.1725 | 0.1433 → 0.1472 | 0.1413 → 0.1457 | 0.169 → 0.143 |
| Gray | 1986 | 0.1737 → 0.1774 | 0.1434 → 0.1471 | 0.1426 → 0.1468 | 0.161 → 0.136 |
| Gray + GloBI | 1984 | 0.1467 → 0.1403 | 0.1413 → 0.1350 | 0.1384 → 0.1324 | 0.186 → 0.221 |
| Gray + GloBI | 1986 | 0.1595 → 0.1468 | 0.1539 → 0.1413 | 0.1506 → 0.1386 | 0.113 → 0.184 |
| Gray + GATEWAy | 1984 | 0.1581 → 0.1533 | 0.1490 → 0.1443 | 0.1393 → 0.1355 | 0.181 → 0.203 |
| Gray + GATEWAy | 1986 | 0.1539 → 0.1496 | 0.1454 → 0.1411 | 0.1368 → 0.1333 | 0.194 → 0.215 |

Indicators among taxa of the webs sampled with the decay and the masses
(1984; mean of 1000 draws, free → ordered):

| Registries | Connectance (observed 0.106) | Omnivores (0.50) | Top taxa (0.02) | Mean trophic level (1.68) |
|---|---|---|---|---|
| Gray | 0.069 → 0.064 | 0.58 → 0.56 | 0.085 → 0.105 | 1.52 → 1.51 |
| Gray + GloBI | 0.086 → 0.079 | 0.69 → 0.65 | 0.070 → 0.076 | 1.60 → 1.60 |
| Gray + GATEWAy | 0.099 → 0.092 | 0.86 → 0.85 | 0.062 → 0.069 | 1.68 → 1.67 |

- With Gray + GloBI and Gray + GATEWAy every score improves (−0.004 to
  −0.013); with the Gray registry alone they get worse (+0.004). With the
  combined registries the free fits gave coarse classes high probabilities
  (Gray + GloBI: order 0.51, family 0.30), and most of the links of Tuesday
  Lake inferred from an order or a class are not observed.
- The ordered classes give a little fewer links (connectance lower,
  further from the observed one), fewer omnivores and a few more top taxa.
- Tuesday Lake was not used to choose the default (the tie in the folds and
  the requirement of the maintainer were); it is the final test.

## Sampling

| Check | Result |
|---|---|
| Frequency of each link in 4000 iterations against its probability | within 4 Monte Carlo errors (test) |
| Seed | reproducible and without altering the global random number stream (test) |

The sections “Webs of the same system”, “Unknown or typical system” and
“Bundled calibration (13 webs; 13 and 50 left out together)” are from the
first phase (maximum likelihood fit and rule with `min_partners = 1`); the
others, from the second (default rule with `min_partners = 5`, Bayesian fit),
the third (sources; by default, every registry at every level) and the
fourth (taxonomic distance decay, experimental).

Scripts: `dev/validation/calibration_lme4.R`, `dev/validation/web_coverage.R`,
`dev/validation/validation_groups.R`, `dev/validation/example_system.R`,
`dev/validation/compare_rules.R`, `dev/validation/scoring_sensitivity.R`,
`dev/validation/tuesday_lake.R`, `dev/globi/process_globi.R`,
`dev/validation/compare_sources.R`, `dev/validation/source_sensitivity.R`,
`dev/validation/decay_background.R`
(order and run times in `dev/validation/README.md`). The comparisons with
the original code (`dev/reference/compare_metaweb.R`,
`dev/reference/compare_indicators.R`) need a local copy of the original
project, which is not distributed.

## R packages cited

Amador Crespo P, Rico A (2026). *bwebs*: Bayesian Food Webs from Compiled
Trophic Interactions. R package version 0.4.0.
<https://github.com/STABLE-uv/bwebs>

Bates D, Mächler M, Bolker B, Walker S (2015). Fitting linear mixed-effects
models using *lme4*. *Journal of Statistical Software*, 67(1), 1–48.
<https://doi.org/10.18637/jss.v067.i01>

Hudson LN, Emerson R, Jenkins GB, Layer K, Ledger ME, Pichler DE, Thompson
MSA, O'Gorman EJ, Woodward G, Reuman DC (2013). *Cheddar*: analysis and
visualisation of ecological communities in R. *Methods in Ecology and
Evolution*, 4(1), 99–104. <https://doi.org/10.1111/2041-210X.12005>

</div>
