__MEDseq: Mixtures of Exponential-Distance Models with Covariates__   
===================================================================

### Improvements, Bug Fixes & Miscellaneous Edits
* Function `seqdef` added as an exact copy of `TraMineR::seqdef`, to enable experienced  
  users of `MEDseq` & `TraMineR` to use the former without needing to explicitly load the latter.
* `wKModes` now also returns `x$tot.withindiff` (i.e. `sum(x$withindiff)`).
* Fixed rare bug in tie-breaking for modal sequence estimate.
* New function `MEDseq_AvePP` added.

## MEDseq v1.3.3 - (_10<sup>th</sup> release [patch update]: 2022-03-28_)
### Improvements, Bug Fixes & Miscellaneous Edits
* Major speed-ups to E-steps for all model types when `G>1`.
* Minor speed-ups to distance calculations for all model types when `G>1`.
* `MEDseq_meantime` gains the `map.size` arg. and a related `print` method.
* Added `summary` (and related `print`) methods for `MEDCriterion` objects.
* New function `MEDseq_entropy` added.
* Fixed mismatched plotting symbols for models with noise in model-selection criteria plot legends.
* Minor fix to handle (rare) empty components.
* Minor edits for compatibility w/ latest `TraMineR` release, w.r.t. `"mt"` and `"ms"` plots.

## MEDseq v1.3.2 - (_9<sup>th</sup> release [patch update]: 2021-12-19_)
### Bug Fixes & Miscellaneous Edits
* Modifications to `WKModes` (& thus related `MEDseq_control` `init.z`  
  options `"kmodes"`/`"kmodes2"`), by further altering `klaR::kmodes`:  
  * Ties for modes now broken randomly, using new `wKModes` arg. `random` (defaults to `TRUE`).
  * All tie-breaks for cluster assignments now biased towards previous iteration's assignments.
  * Fixed rare bug when `modes` is supplied as a number with aggregated data, e.g. `"kmodes2"`.
* `MEDseq_fit` & other functions now work for sequence alphabets of any size;  
  previously, only sequences with fewer than 10 states/categories were accommodated.
* Minor fix to `dbs` function when supplying `clusters` with a noise component.
* `sapply` replaced with `vapply`, with other negligible speed-ups.
* Updated citation info after final publication in _JRSSA_.

## MEDseq v1.3.1 - (_8<sup>th</sup> release [patch update]: 2021-10-14_)
### Bug Fixes & Miscellaneous Edits
* Fixes for `init.z` options `"kmodes"` & `"kmodes2"` in `MEDseq_control`, with new function `wKModes`  
  provided for running the k-modes algorithm on _weighted_ data: previously, k-modes initialisation  
  was only available for _unweighted_ sequences via the now-replaced `klaR::kmodes` function  
  (consequently, the `klaR` package has been removed from the `DESCRIPTION` `Suggests:` field).
* `plot.MEDseq` gains the arg. `subset`, for use with the `TraMineR` `type` plots:  
  allows plotting some but not all components, e.g. only the noise component (see documentation).
* Fixed minor bug causing `MEDseq_fit` to crash when `weights` are supplied and `unique=FALSE`.
* Fixed ASW calculation when _unweighted_ sequences are aggregated (i.e. `unique=TRUE`, the default).
* Fixed small bug for `type="ms"` plots for models with a noise component when `SPS=TRUE`.
* Fixed printing of `noise.gate` in `MEDseq_compare` for `G=2` models w/ noise & gating covariates.
* Improved checks on `G` in `MEDseq_fit`.

## MEDseq v1.3.0 - (_7<sup>th</sup> release [minor update]: 2021-07-15_)
### New Features & Improvements
* `plot.MEDseq` gains a number of new arguments:  
    * `soft` allows soft cluster membership probabilities to be used for the `"d"`, `"f"`, `"Ht"`, `"ms"`,  
    & `"mt"` `type` plots (default: `soft=TRUE`) + the `"i"` & `"I"` plots (default: `soft=FALSE`), in a  
    manner akin to  `WeightedCluster::fuzzyseqplot()`:  previously, all but the `"ms"` plot used the  
    hard MAP partition and discarded the soft assignment information (i.e. `soft=FALSE`, implicitly).
    * `sortv` allows overriding the `smeth` arg. to instead order observations in certain plots  
    (where `seriated` is one of `"observations"` or `"both"`) by the `"dbs"` or `"asw"` values;  
    additionally, and for consistency with `WeightedCluster::fuzzyseqplot()`,  
    `sortv="membership"` is provided for `soft=TRUE` `type="I"` plots.
    * `weighted` (`TRUE`, by default) allows control over whether the weights (if any) are used;  
    relevant only for `"d"`, `"f"`, `"Ht"`, `"i"`, `"I"`, `"ms"`, & `"mt"` `type` plots.
* Exported `MEDseq_clustnames` & `MEDseq_nameclusts` functions and added `SPS` arg. to `plot.MEDseq`,  
  `MEDseq_meantime`, `MEDseq_stderr`, & various/more `print`/`summary` methods: now certain plots &  
  outputs can be (or are by default) labelled with the central sequences in SPS format, as per the paper.
* `seriated` options `"observations"` & `"both"` can now be used for `"i"` type plots,  
  with related minor fixes for `"i"` & `"I"` type plots for weighted data with seriated observations.
* Added `predict`, `fitted`, & `residuals` methods for `"MEDgating"` objects, i.e. `x$gating`.
* `MEDseq_meantime` gains the arg. `wt.size` (defaults to `FALSE`).
* Minor speed-ups to model-fitting for `modtype="CU"`.

### Bug Fixes & Miscellaneous Edits
* A warning message is now printed if the gating network's MLR ever fails to converge, prompting users to  
  modify the `itmax` arg. to `MEDseq_control`: the 2<sup>nd</sup> element of this arg. governs the maximum number of  
  MLR iterations --- consequently, its default has been modified from `100` to `1000`, which is liable to slow  
  down internal calls to `nnet::multinom`, but generally reduces the required number of EM iterations. 
* Changes to default colour palettes & plotting symbols for certain plot types:  
  `Suggests:` package `viridisLite` now only invoked if available.
* Minor fixes to returned `x$gating` object, especially for `equalPro` models  
  with a noise component and weighted models _without_ any gating covariates at all.
* Stronger checks to ensure `weights` arg. is explicitly supplied to `MEDseq_fit`  
  in cases where the `"stslist"` object passed via `seqs` has the `"weights"` attribute.
* Added error message to `MEDseq_fit` when the number of states exceeds 9,   
  to better inform of this bug which will be rectified in future updates.
* Fixed bug preventing inclusion of higher-order terms in `gating` formulas when there are duplicates.
* Minor fixes to `get_MEDseq_results` and how its optional args. are internally handled by `plot.MEDseq`.
* Stronger checks for variables in `gating` formula which are not found in `covars`.
* `type="mean"` option renamed to `type="central"` in `plot.MEDseq`.
* `type="ms"` plots now work properly when `seriated="clusters"` or `seriated="both"`.
* Removed some superfluous warnings for all but the `"mt"` `TraMineR` `type` plots.
* Fixed small bug in `MEDseq_meantime` when `MAP=FALSE`.
* Further robustifications to handle empty components.
* Minor fixes to `print.MEDseq` for models where DBS &/or ASW statistics weren't computed.
* Minor vignette edits and documentation clarifications.
* Updated citation info after online publication in _JRSSA_.

## MEDseq v1.2.1 - (_6<sup>th</sup> release [patch update]: 2020-12-29_)
### Bug Fixes & Miscellaneous Edits
* The `"d"`, `"f"`, `"Ht"`, `"i"`, & `"I"` plot types now properly account for sampling weights.
* Layout and legend-placement has been improved for these same types of plots.
* Mimicking `TraMineR` further, `plot.MEDseq` also gains the `type` options `"ms"` & `"mt"`.
* Minor speed-ups associated with the `opti="medoid"` setting.
* Added ORCID iDs to DESCRIPTION.
* Minor CRAN compliance edits to the vignette.

## MEDseq v1.2.0 - (_5<sup>th</sup> release [minor update]: 2020-11-20_)
### Significant User-Visible Changes
* Corrected the parameter count penalty for the BIC, ICL, and AIC model selection criteria,  
  specifically, the count is now greater for the central sequence position estimates.
* Hence, `criterion="bic"` is now the default for `MEDseq_control`, `MEDseq_compare`, and  
  `get_MEDseq_results` (previously `"dbs"`), with modifications to `print` & `summary` functions.
* Non-noise components' central sequence positions associated with precision parameters of zero  
  are now printed (`print.MEDseqtheta`) & plotted (`plot.MEDseq(..., type="mean")`) always:  
  the `preczero` argument has thus been removed from both functions.
  
### New Features & Improvements
* `MEDseq_meantime` gains two new arguments (see documentation for more details):  
    * `weighted` (default: `TRUE`, old: `FALSE`) allows the sampling weights to be used,  
    with or without the cluster assignment probabilities, in the computation of the weighted averages.
    * `prop` (default: `FALSE`) divides the output when `norm=TRUE` by the sequence length.
* `MEDseq_control` gains the arg. `random=TRUE`, governing tie-breaking of estimated central sequence  
  positions: old behaviour (always choosing the first candidate state) recoverable via `random=FALSE`.
* `plot.MEDseq` arg. `quant.scale=FALSE` replaces old arg. `log.scale`: quantiles now used  
  to determine non-linear colour breakpoints when invoked with `type="precision"`.
* Sped-up `init.z="kmedoids"` initialisation via `pam` for _unweighted_ sequences, by using the  
  _highest available_ value for the `pamonce` option,  based on the `cluster` package's version number.
* `init.z` gains the options `"kmodes"` & `"kmodes2"`, though only for _unweighted_ sequences:  
  both require the newly _suggested_ `klaR (>= 0.6-13)` package.
* `plot.MEDseq` gains the arg. `smeth`, governing the seriation method to be used (`"TSP"`, by default).
* For weighted sequences, `init.z="kmedoids"` is now itself initialised by Ward's hierarchical clustering.
* Significant speed-ups to computation of central sequences for all `opti` settings (esp. `"mode"`).
* Added `SPS` arg. (default=`FALSE`) to `print.MEDtheta` & `summary.MEDseq`.
* `dbs` gains the optional/experimental arg. `clusters` - no change to default.
  
### Bug Fixes & Miscellaneous Edits
* Various fixes to the `seriated` arg. to `plot.MEDseq`:  
    * Arg. name changed from `seriate` to avoid conflict with function `seriation::seriate`.  
    * Fixed `seriated` options `"clusters"`/`"both"` for models with no noise component.  
    * `seriated="observations"` (the default) now also works for `type="I"` plots.  
    * `seriated="clusters"` now also works for `type="dbsvals"` & `type="aswvals"` plots. 
* `MEDseq_fit` now always internally normalises the `weights` to sum to the sample size.
* Minor fixes to properly account for weighted sequences &/or duplicates when `noise.gate=FALSE`.
* Minor fix to gathering of results to account for `noise.gate=FALSE` when `G=2`.
* `MEDseq_stderr` now respects the `algo`, `opti`, & `noise.gate` settings of the original model.
* `MEDseq_compare` now returns & prints `opti` info where relevant.
* Fixes to `print` & `summary` methods for `MEDgating` objects if `equalPro=TRUE`.
* `MEDseq_fit` now coerces `"character"` covariates to `"factor"`.
* Minor fixes to `print` method for `MEDlambda` objects also.
* Additional minor edits to `plot.MEDseq(..., type="gating")`.
* `print.MEDseqCompare` gains the args. `maxi` & `rerank=FALSE`.
* Minor speed-ups for `G=1` models.
* Added `viridisLite (>= 0.2.0)` to `Suggests:` (for `plot.MEDseq(..., type="precision")`).
* Ensured `matrixStats (>= 0.53.1)` and `TraMineR (>= 1.6)` in `Imports:`.
* Package startup message now checks if newer version of package is available from CRAN.
* Significant vignette edits.
* Updated maintainer e-mail address.
* Minor documentation, examples, and CRAN compliance edits.

## MEDseq v1.1.1 - (_4<sup>th</sup> release [patch update]: 2020-05-12_)
### Bug Fixes & Miscellaneous Edits
* Maintenance release for compatibility with R 4.0.0 - minor edits.
* `summary.MEDseq` gains the printing-related arguments  
  `classification=TRUE`, `parameters=FALSE`, and `gating=FALSE`.
* `x$params$lambda` now inherits the `MEDlambda` class,  
  with its own `print` method as per `x$params$theta`.
* `x$params$tau` now has informative `dimnames`.
* Minor changes when supplying `x.axis` to `plot.MEDseq(..., type="gating")`.
* Documentation, vignette, examples, and references improvements.
* Added `rmarkdown` to `Suggests:`.
* Reformatted package startup message.

## MEDseq v1.1.0 - (_3<sup>rd</sup> release [minor update]: 2020-03-30_)
### New Features, Improvements, & Bug Fixes
* Significant efficiency gains when ignoring duplicates in the presence of weights:  
    * before, unique cases were defined as unique sequence/covariates/weight combinations,  
    * now, cases with different weights that are otherwise duplicates are treated as duplicates.
* `MEDseq_stderr` is provided for computing the standard errors of the  
  coefficients for the covariates in the gating network via either the  
  weighted likelihood bootstrap or jackknife methods.
* Small robustifications in the presence of empty components.
* Fixed `get_MEDseq_results` when `what="MAP"` and non-noise models are chosen.
* Fixed bug related to the colours used in the vignette plots.
* Odds ratios now returned (and printed) when calling `summary` on `x$gating`.
* Cosmetic fixes to `plot.MEDseq` when `type="clusters"` for small sample sizes.
* Other small cosmetic plotting  & reference-formatting changes.
* Spell-checking of documentation and fixes to `donttest` examples.

## MEDseq v1.0.1 - (_2<sup>nd</sup> release [patch update]: 2019-12-10_)
### Bug Fixes & Miscellaneous Edits
* Speed-ups to E-step, especially for models with a noise component.
* Clarifications and improvements to documentation and examples.

## MEDseq v1.0.0 - (_1<sup>st</sup> release: 2019-08-24_)
