__MEDseq: Mixtures of Exponential-Distance Models with Covariates__   
===================================================================

### Bug Fixes & Miscellaneous Edits
* Sped-up `"kmedoids"` initialisation via `pam` for _unweighted_ sequences  
  by using the _highest available_ value for the `pamonce` option,  
  based on the version number of the loaded `cluster` package.
* Ensured `matrixStats (>= 0.53.1)` and `TraMineR (>= 1.6)` in `Imports:`.

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