__MEDseq: Mixtures of Exponential-Distance Models with Covariates__   
===================================================================

### Bug Fixes & Miscellaneous Edits
* Added `predict`, `fitted`, & `residuals` methods for `"MEDgating"` objects, i.e. `x$gating`.
* Minor fixes to returned `x$gating` object, especially for `equalPro` models  
with a noise component and weighted models _without_ any gating covariates at all.
* Stronger checks to ensure `weights` arg. is explicitly supplied to `MEDseq_fit`  
in cases where the `"stslist"` object passed via `seqs` has the `"weights"` attribute.
* Stronger checks for variables in `gating` formula which are not found in `covars`.
* Minor fixes to `print.MEDseq` for models where DBS &/or ASW statistics weren't computed.
* Removed some superfluous warnings for all but the `"mt"` `TraMineR` type plots.
* Minor vignette edits.

## MEDseq v1.2.1 - (_6<sup>th</sup> release [patch update]: 2020-12-29_)
### Bug Fixes & Miscellaneous Edits
* The `"d"`, `"f"`, `"Ht"`, `"i"`, & `"I"` plot types now properly account for sampling weights.
* Layout and legend-placement has been improved for these same types of plots.
* Mimicking `TraMineR` further, `plot.MEDseq` also gains the `type` options `"ms"` & `"mt"`.
* Minor speed-ups associated with the `opti="medoid"` setting.
* Added ORCID iDs to DESCRIPTION.
* Minor CRAN compliance edits to the vignette.

## MEDseq v1.2.0 - (_5<sup>th</sup> release [minor update]: 2020-11-20_)
### Significant user-visible changes
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
### New Features, Improvements, and Bug Fixes
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