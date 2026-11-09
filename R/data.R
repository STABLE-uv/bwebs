#' Registry of trophic interactions of Gray et al. (2015)
#'
#' The database of published trophic interactions among freshwater taxa,
#' mostly from the United Kingdom, compiled by Gray et al. (2015), in the
#' format of [as_registry()]. The original file has 20,823 rows and 51 source
#' identifiers; one of them is a test record ("Test"), which is removed.
#' Five names whose hierarchy was resolved wrongly (homonyms or missing
#' ranks) are corrected, and misspellings and synonyms of higher ranks are
#' harmonised with [default_taxonomy_fixes()] (see `data-raw/gray2015.R`).
#'
#' @format A data frame with 20,822 rows (one per record) and 23 character
#'   columns: `record_id` (link identifier of the original file),
#'   `registry` (`"gray2015"`), `study_id` (source identifier of the
#'   original file), `resource` and `consumer` (names as recorded),
#'   `resource_std` and `consumer_std` (names after harmonisation),
#'   `res_stage` and `con_stage` (`"larva"`, `"adult"` or `NA`), the
#'   classification of the resource (`res_genus`, `res_subfamily`,
#'   `res_family`, `res_order`, `res_class`) and of the consumer (`con_*`),
#'   `NA` where not resolved, `res_category` and `con_category` (`"producer"`,
#'   `"invertebrate"`, `"vert.ecto"`), `evidence` (type of evidence) and
#'   `citation` (full reference of the source).
#' @source Gray C, Ma A, Perkins D, Hudson L, Figueroa D, Woodward G (2015).
#'   Database of trophic interactions. Zenodo, \doi{10.5281/zenodo.13751},
#'   file `trophic.links.2014-11-10.csv`, licence CC0 1.0. Described in Gray
#'   et al. (2015), *Food Webs*, 5, 11--20, \doi{10.1016/j.fooweb.2015.09.001}.
"gray2015"

#' Sources of the Gray et al. (2015) registry
#'
#' One row per source (study) of [gray2015], with the counts used to choose
#' the validation webs of [validate_rule()] and [calibrate_links()]. Study 12
#' ("from Brose et al 2005 (Jonsson et al. (2004))") is Tuesday Lake: its 382
#' links are those of [tuesday_lake_links]. Leave it out of the registry when
#' Tuesday Lake is used to test a method. Studies 13 and 50 describe the same
#' system (see [gray2015_calibration]).
#'
#' @format A data frame with 50 rows and 9 columns: `study_id`, `citation`,
#'   `n_records`, `n_links` (distinct resource-consumer pairs), `n_consumers`,
#'   `n_resources`, `n_taxa`, `n_repeated_records` (records that repeat a
#'   pair, for example from several sites of the same study) and `validation`
#'   (logical: provisional set of whole communities with at least eight
#'   consumers whose diets were recorded, used for [gray2015_calibration]).
#' @source Derived from [gray2015]; see `data-raw/gray2015.R`.
"gray2015_sources"

#' Calibration of link probabilities for the Gray et al. (2015) registry
#'
#' The Bayesian calibration ([calibrate_links()], 4 chains of 4000 draws,
#' ordered classes of evidence) of [gray2015] with the default
#' [matching_rule()], including the model of the taxonomic decay of the
#' pairs not inferred (`decay`), and the 13 validation
#' webs marked in [gray2015_sources]; Tuesday Lake (study 12) was kept apart
#' only while the rule was chosen and tested ([gray2015_rule_comparison]).
#' Studies 13 and 50 are left out together (`groups`): 62% of the links of
#' study 13 are also recorded in study 50 (see [shared_studies()]), so they
#' describe the same system, and when each was left out alone the other
#' predicted it. Use it with [sample_metawebs()] for metawebs built with the
#' same rule from the same registry; to test a method on Tuesday Lake,
#' calibrate it without study 12.
#'
#' @format An object of class `bwebs_calibration` (see [calibrate_links()]).
#' @source `data-raw/gray2015_calibration.R`.
"gray2015_calibration"

#' Calibration of link probabilities for Gray et al. (2015) and GloBI
#'
#' The Bayesian calibration ([calibrate_links()], 4 chains of 4000 draws,
#' ordered classes of evidence) of [gray2015] combined with records of
#' Global Biotic Interactions (GloBI),
#' with the default [matching_rule()] (every registry at every level, with
#' the taxonomic decay of the pairs not inferred) and the 13 validation webs
#' marked in [gray2015_sources]. Studies 13 and 50 are left
#' out together, and when each web is rebuilt the GloBI studies holding at
#' least 10% of its links are left out too (copies of the web in other
#' compilations; for Tuesday Lake, study 12, also the GloBI studies that copy
#' it). Use it with [sample_metawebs()] for metawebs built with the same rule
#' from [gray2015] combined with GloBI records downloaded with
#' [fetch_globi()] and converted with [globi_to_registry()] (consumers
#' queried by name, genus and family, as here).
#'
#' The GloBI records were downloaded on 3 November 2026 for the consumers of
#' the validation webs, of Tuesday Lake and of the azoxystrobin experiment of
#' Amador et al. (2024); they are not included in the package (their licences
#' depend on each source). Only the counts of the validation (by class of
#' evidence, and by evidence profile for the decay) and the posterior draws
#' are.
#'
#' @format An object of class `bwebs_calibration` (see [calibrate_links()]).
#' @source `data-raw/gray2015_globi_calibration.R`.
"gray2015_globi_calibration"

#' Example community
#'
#' Invented data for the examples: 23 freshwater nodes (four diatoms,
#' detritus and 17 invertebrates, classified as in [gray2015]) and their
#' abundances in six samples, three controls and three exposed to a stressor
#' that removes mayflies, stoneflies, amphipods and the shrimp *Dugastella
#' valentina*, and the biomasses and masses of their populations for the
#' energy fluxes. The shrimp has no records in the registry and takes those
#' of *Asellus* through a surrogate.
#'
#' @format
#' `example_taxa`: a data frame with 23 rows (one per node) in the format of
#' [prepare_nodes()]: `node_id`, `taxon`, `rank`, `role`, `node_type`,
#' `genus`, `subfamily`, `family`, `order`, `class`, `phylum`, `alias`,
#' `consumer_surrogate` and `resource_surrogate`.
#'
#' `example_community`: a data frame with 126 rows (one per sample and
#' sampled node) and columns `sample_id`, `node_id` and `value` (individuals,
#' or cells for diatoms, per sample; 0 is an observed absence). Detritus is
#' not sampled: use it as a permanent resource.
#'
#' `example_samples`: a data frame with 6 rows and columns `sample_id`,
#' `treatment` (`"control"` or `"exposed"`) and `replicate`.
#'
#' `example_populations`: a data frame with 138 rows (one per sample and
#' node, detritus included) for [calc_fluxes()] and [node_masses()], with
#' columns `sample_id`, `node_id`, `biomass` (g of wet mass per sample) and
#' `mass` (mean individual wet mass of the consumers, g; `NA` for producers
#' and detritus): the abundances of `example_community` times a typical
#' mass of each taxon, 2 mg of periphyton per counted unit of a diatom and
#' 5 g of CPOM and 2 g of FPOM in every sample.
#' @source Invented for the package; see `data-raw/example.R`.
#' @name example_community
NULL

#' @rdname example_community
"example_taxa"

#' @rdname example_community
"example_community"

#' @rdname example_community
"example_samples"

#' @rdname example_community
"example_populations"

#' Tuesday Lake food webs, 1984 and 1986
#'
#' The pelagic communities of Tuesday Lake (Michigan, USA) sampled in 1984 and
#' 1986: 82 taxa (50 producers, 29 invertebrates and 3 fish), their numerical
#' abundance and body mass, and the observed trophic links of each year (269
#' in 1984, 241 in 1986), as compiled in the *cheddar* package (Hudson et al.
#' 2013). In *bwebs* they
#' serve as an independent test of the construction of food webs: the
#' registry used to build them must not contain study 12 of [gray2015], which
#' is Tuesday Lake.
#'
#' @format
#' `tuesday_lake_taxa`: a data frame with 82 rows (one per node) in the format
#' of [prepare_nodes()]: `node_id`, `taxon` (name matched with the registry),
#' `rank`, `role` (`"producer"` or `"consumer"`), `node_type`, `genus`,
#' `subfamily`, `family`, `order`, `class` and `phylum`, aligned with the
#' taxonomy of [gray2015] (the rotifer class is Eurotatoria, as in the
#' registry), plus `category` (`"producer"`, `"invertebrate"`, `"vert.ecto"`)
#' and `name` (the node name in *cheddar*). Two producers, unclassified
#' flagellates and microflagellates, are non-taxonomic nodes.
#'
#' `tuesday_lake_community`: a data frame with 113 rows (one per year and
#' node present) and columns `sample_id` (`"TL84"`, `"TL86"`), `node_id`,
#' `value` (numerical abundance, individuals per m\out{<sup>3</sup>}) and
#' `mass_kg` (mean body mass, kg).
#'
#' `tuesday_lake_links`: a data frame with 510 rows (one per observed link)
#' and columns `sample_id`, `resource` and `consumer` (node identifiers);
#' it includes cannibalistic links.
#' @source Communities TL84 and TL86 of the *cheddar* R package (Hudson et al.
#'   2013; version 0.1-639, BSD 2-Clause licence, see the file `COPYRIGHTS`
#'   of *bwebs*), from Cohen et al. (2003) and Jonsson et al. (2005); see
#'   `data-raw/tuesday_lake.R`.
#' @references
#' Cohen JE, Jonsson T, Carpenter SR (2003). Ecological community description
#' using the food web, species abundance, and body size. *Proceedings of the
#' National Academy of Sciences of the USA*, 100(4), 1781--1786.
#' \doi{10.1073/pnas.232715699}
#'
#' Hudson LN, Emerson R, Jenkins GB, Layer K, Ledger ME, Pichler DE, Thompson
#' MSA, O'Gorman EJ, Woodward G, Reuman DC (2013). *Cheddar*: analysis and
#' visualisation of ecological communities in R. *Methods in Ecology and
#' Evolution*, 4(1), 99--104. \doi{10.1111/2041-210X.12005}
#'
#' Jonsson T, Cohen JE, Carpenter SR (2005). Food webs, body size, and
#' species abundance in ecological community description. *Advances in
#' Ecological Research*, 36, 1--84. \doi{10.1016/S0065-2504(05)36001-6}
#' @name tuesday_lake
NULL

#' @rdname tuesday_lake
"tuesday_lake_taxa"

#' @rdname tuesday_lake
"tuesday_lake_community"

#' @rdname tuesday_lake
"tuesday_lake_links"

#' Comparison of matching rules on the Gray et al. (2015) registry
#'
#' The cross-validated scores of the [candidate_rules()] on the 12 validation
#' webs of [gray2015] other than Tuesday Lake (left out by system, studies 13
#' and 50 together), from [compare_rules()], and their scores on the
#' Tuesday Lake webs of 1984 and 1986 ([score_metaweb()]) with a registry and
#' a calibration without Tuesday Lake. With this registry alone, climbing
#' from family to class with a threshold of five or one partners scored best,
#' indistinguishably; the defaults of [matching_rule()] were chosen with the
#' comparison of sources ([source_comparison]).
#'
#' @format A data frame with 21 rows (one per rule, from the lowest mean
#'   Brier score) and columns `rule`, `mode`, `base_level`, `climb_max`,
#'   `min_partners`, `downward` (downward inheritance), the cross-validated
#'   mean scores over webs (`brier`, its standard error `brier_se`,
#'   `brier_diff` and `brier_diff_se`, the difference from the best rule and
#'   its standard error from the differences web by web, `skill`, `tss`,
#'   `precision`, `recall`) and `n_predicted` (links inferred in the 12
#'   webs), and the mean scores on Tuesday Lake (`tl_brier`, `tl_skill`,
#'   `tl_tss`, `tl_precision`, `tl_recall`; `NA` precision when no link was
#'   inferred). See [score_rule()] for the definitions.
#' @source `data-raw/gray2015_rule_comparison.R`, from
#'   `dev/validation/compare_rules.R` and `dev/validation/tuesday_lake.R`
#'   (Bayesian refits with seed 2026).
"gray2015_rule_comparison"

#' Comparison of sources of trophic interactions
#'
#' The cross-validated scores of the main matching rules with the registry of
#' Gray et al. (2015) ([gray2015]), with records of Global Biotic Interactions
#' (GloBI) and with both, on the 12 validation webs of [gray2015] other than
#' Tuesday Lake (left out by system, studies 13 and 50 together; when each web
#' is rebuilt, the GloBI studies holding at least 10% of its links are left
#' out too), from [compare_rules()], and their scores on the Tuesday Lake
#' webs of 1984 and 1986 ([score_metaweb()]) with registries and calibrations
#' without Tuesday Lake. The combination with the lowest mean Brier score sets
#' the defaults of [matching_rule()].
#'
#' The GloBI records were downloaded on 3 November 2026 for the consumers of
#' the validation webs, of Tuesday Lake and of the azoxystrobin experiment of
#' Amador et al. (2024), by name, genus and family (500 queries, 832,003
#' records; 244,452 trophic records of 17,573 studies after removing the
#' copies of the Gray et al. 2015 and Brose et al. 2005 compilations, the
#' records of other taxa than the one queried and the studies that copy
#' Tuesday Lake). They are not included in the package.
#'
#' @format A data frame with one row per source and rule (from the lowest
#'   mean Brier score) and columns `source` (`"Gray"`, `"GloBI"` or
#'   `"Gray + GloBI"`), `primary_registry` and `supplementary_max_level`
#'   (how the two registries were combined, see [matching_rule()]; `NA`:
#'   every registry at every level), `rule`, the cross-validated mean scores
#'   over webs (`brier`, `brier_se`, `brier_diff` and `brier_diff_se`, the
#'   difference from the best combination and its standard error from the
#'   differences web by web, `diff_vs_gray` and `diff_vs_gray_se`, the same
#'   for the same rule with the Gray registry alone, `skill`, `tss`,
#'   `precision`, `recall`), `n_predicted` (links inferred in the 12 webs)
#'   and the mean scores on Tuesday Lake (`tl_brier`, `tl_skill`, `tl_tss`,
#'   `tl_precision`, `tl_recall`). See [score_rule()] for the definitions.
#' @source `data-raw/source_comparison.R`, from
#'   `dev/validation/compare_sources.R` (Bayesian refits with seed 2026).
"source_comparison"

#' Effect of the predator/prey body-mass ratio on links
#'
#' The size term of [sample_metawebs()] and [score_metaweb()] with
#' `masses`: the within-consumer effect of \eqn{x = \log_{10}(M_{consumer} /
#' M_{resource})}{x = log10(M_consumer / M_resource)} on the links of each
#' predator with the animals of its web, a logistic regression with one
#' intercept per consumer and a quadratic term in \eqn{x}, fitted to the 206
#' aquatic webs of GATEWAy version 1.0 (Brose et al. 2019; 49 lakes, 23
#' streams and 134 marine webs; 491,829 predator x animal pairs with 47,467
#' links). Tuesday Lake was left out, so that it stays an independent test,
#' and so were the five chalk streams of Gray et al. (2015) and Thompson et
#' al. (2017), whose links were partly extrapolated from similar taxa. Only
#' the coefficients are bundled.
#'
#' The shift of the logit of a predation pair is \eqn{w\, b\, (x -
#' x^*)^2}{w b (x - x*)^2}, 0 at the optimum ratio \eqn{x^*} (a consumer
#' about 17,000 times heavier than its prey; the curve is flat: \eqn{b\, (x -
#' x^*)^2}{b (x - x*)^2} is -0.24 for a ratio 100 times smaller or larger than
#' the optimum, and -1.1 for a consumer as heavy as its prey), with weight
#' \eqn{w} = 0.25 for the links inferred by the rule and 1 for the pairs of
#' the taxonomic decay. The curve and the weights are those that maximised the
#' validation by folds: the 12 validation webs of the registry with body
#' masses for at least 60% of their consumers, left out by system, with the
#' probabilities and the curve fitted to the other systems, against curves
#' fitted to the webs of Brose et al. (2005) and *cheddar* (Hudson et al. 2013),
#' to the freshwater
#' webs of GATEWAy or to all of its webs, one curve or one per kind of
#' predator, and weights from 0 to 1 and from 1 to 30. The mean Brier score
#' of those webs went down by 0.006 (9 of 12 webs better; 0.008 with the
#' Gray registry, 0.005 with Gray + GloBI and 0.006 with Gray + GATEWAy), and
#' that of Tuesday Lake from 0.140 to 0.137 with Gray + GloBI and from 0.149
#' to 0.139 with Gray + GATEWAy. The weight of the inferred links was
#' limited to 0.25: larger weights gave about the same mean score but did
#' worse on a web when they were chosen without it, because in webs whose
#' recorded links are mostly predation (Broadstone Stream, Skipwith Pond)
#' the size term lowers many recorded links (`dev/validation/results.md` in
#' the GitHub repository of the package). Change `weights` to use other
#' weights.
#'
#' @format An object of class `bwebs_size_model`, a list with
#'   `coefficients` (`x`, `x2`: the linear and quadratic terms in \eqn{x}),
#'   `std_error`, `optimum` (\eqn{x^*}), `weights` (`rule`, `decay`),
#'   `data` (webs, pairs and links by ecosystem) and `sources`.
#' @source Brose U, et al. (2019) Predator traits determine food-web
#'   architecture across ecosystems. *Nature Ecology & Evolution* 3:
#'   919--927, \doi{10.1038/s41559-019-0899-x}. GlobAL daTabasE of traits and
#'   food Web Architecture (GATEWAy) version 1.0, iDiv Data Repository,
#'   dataset 283 (CC BY 4.0). Fitted by `data-raw/size_model.R`.
#' @references
#' Hudson LN, Emerson R, Jenkins GB, Layer K, Ledger ME, Pichler DE, Thompson
#' MSA, O'Gorman EJ, Woodward G, Reuman DC (2013). *Cheddar*: analysis and
#' visualisation of ecological communities in R. *Methods in Ecology and
#' Evolution*, 4(1), 99--104. \doi{10.1111/2041-210X.12005}
"size_model"

#' Registry of the freshwater trophic links of GATEWAy
#'
#' The links of the 81 freshwater food webs (lakes and streams) of the
#' GlobAL daTabasE of traits and food Web Architecture (GATEWAy) version 1.0
#' (Brose et al. 2019) in the format of [as_registry()]: 19,185 records,
#' one study per food web (`study_id` `"gateway: <food web>"`). Left out:
#' Tuesday Lake (the independent test of the package), the five chalk
#' streams of Gray et al. (2015) and Thompson et al. (2017), whose links were
#' partly extrapolated from similar taxa with WebBuilder, and parasitic and
#' parasitoid links. GATEWAy gives only the names of the taxa: their genus,
#' subfamily, family, order and class are those of the same names in
#' [gray2015] and in the GloBI records of [gray2015_globi_calibration] (by
#' name, then by genus or higher rank), or of the GBIF backbone for the rest
#' (98.6% of the consumers and 92% of the resources resolved to a class; the
#' rest are mostly detritus, groups such as "Phytoplankton" and the local
#' names of the fish of Lake Malawi). Combine it with [gray2015] with
#' [combine_registries()], as [gray2015_gateway_calibration] does.
#'
#' Several webs of GATEWAy copy validation webs of [gray2015_sources]
#' (Broadstone Stream, Mill Stream, Skipwith Pond, the Icelandic streams and
#' the streams of Layer et al. 2010): when such a web is validated, leave
#' them out with it (`exclusions` of [validate_rule()];
#' `dev/gateway/gateway_exclusions.csv` in the GitHub repository of the
#' package, made by `dev/gateway/process_gateway.R`).
#'
#' @format A data frame with 19,185 rows and the 23 columns of [gray2015]
#'   (`registry` `"gateway"`); `evidence` gives how the link was established
#'   in GATEWAy and `citation` its original source.
#' @source Brose U, et al. (2019) Predator traits determine food-web
#'   architecture across ecosystems. *Nature Ecology & Evolution* 3:
#'   919--927, \doi{10.1038/s41559-019-0899-x}. GlobAL daTabasE of traits and
#'   food Web Architecture (GATEWAy) version 1.0, iDiv Data Repository,
#'   dataset 283, file `283_2_FoodWebDataBase_2018_12_10.csv`, licence CC BY
#'   4.0. See `data-raw/gateway.R`.
"gateway"

#' Calibration of link probabilities for Gray et al. (2015) and GATEWAy
#'
#' The Bayesian calibration ([calibrate_links()], 4 chains of 4000 draws,
#' ordered classes of evidence) of [gray2015] combined with [gateway], with
#' the default [matching_rule()]
#' (every registry at every level, with the taxonomic decay of the pairs not
#' inferred) and the 13 validation webs marked in [gray2015_sources].
#' Studies 13 and 50 are left out together, and when each web is rebuilt the
#' webs of GATEWAy that copy it are left out too (those holding at least 10%
#' of its links, alone or pooled by publication). In the cross-validation of
#' the 12 webs other than Tuesday Lake, Gray + GATEWAy had a Brier score of
#' 0.2052, against 0.2093 with the Gray registry alone and 0.1971 with Gray
#' + GloBI (free classes of evidence; GloBI cannot be bundled; adding
#' GATEWAy to Gray + GloBI gave 0.2013), and with the ordered classes of
#' this calibration 0.2044, against 0.2104 and 0.1978. Use it with
#' [sample_metawebs()] or [link_probabilities()] for metawebs built with the
#' same rule from `combine_registries(gray2015, gateway)`.
#'
#' @format An object of class `bwebs_calibration` (see [calibrate_links()]).
#' @source `data-raw/gray2015_gateway_calibration.R`.
"gray2015_gateway_calibration"
