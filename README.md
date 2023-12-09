[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/MEDseq)](https://cran.r-project.org/package=MEDseq)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/MEDseq?)](https://github.com/r-hub/cranlogs.app)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/grand-total/MEDseq?color=82b4e8)](https://github.com/r-hub/cranlogs.app)

# MEDseq R Package
## Mixtures of Exponential-Distance Models
## for Clustering Longitudinal Life-Course Sequences
## with Gating Covariates and Sampling Weights
### Written by Keefe Murphy
## Description

Fits _MEDseq_ models introduced by Murphy et al. (2021) <[doi:10.1111/rssa.12712](https://rss.onlinelibrary.wiley.com/doi/abs/10.1111/rssa.12712)>, i.e. fits mixtures of exponential-distance models for clustering longitudinal/categorical life-course sequence data via the EM/CEM algorithm. A family of parsimonious precision parameter constraints are accommodated. So too are sampling weights. Gating covariates can be supplied via formula interfaces. Visualisation of the results of such models is also facilitated.

The most important function in the __MEDseq__ package is: `MEDseq_fit`, for fitting the models via EM/CEM. This function requires the data to be in `"stslist"` format; the function `seqdef` is conveniently reexported from the __TraMineR__ package for this purpose.

`MEDseq_control` allows supplying additional arguments which govern, among other things, controls on the initialisation of the allocations for the EM/CEM algorithm and the various model selection options. `MEDseq_compare` is provided for conducting model selection between different results from using different covariate combinations &/or initialisation strategies, etc. `MEDseq_stderr` is provided for computing the standard errors of the coefficients for the covariates in the gating network.

A dedicated plotting function exists for visualising various aspects of the results, using new methods as well as some existing methods adapted from the __TraMineR__ package. Finally, the package also contains two data sets: `biofam` and `mvad`.

## Installation

You can install the latest stable official release of the `MEDseq` package from CRAN:

```
install.packages("MEDseq")
```

or the development version from GitHub:

```
# If required install devtools:  
# install.packages('devtools')  
devtools::install_github('Keefe-Murphy/MEDseq')
```

In either case, you can then explore the package with:

```
library(MEDseq)  
help(MEDseq_fit) # Help on the main modelling function
```

For a more thorough intro, the vignette document is available as follows:

```
vignette("MEDseq", package="MEDseq")
```

However, if the package is installed from GitHub the vignette is not automatically created. It can be accessed when installing from GitHub with the code:

```
devtools::install_github('Keefe-Murphy/MEDseq', build_vignettes = TRUE)
```

Alternatively, the vignette is available on the package's CRAN page.

### References
Murphy, K., Murphy, T. B., Piccarreta, R., and Gormley, I. C. (2021). Clustering longitudinal life-course sequences using mixtures of exponential-distance models. _Journal of the Royal Statistical Society: Series A (Statistics in Society)_, 184(4): 1414--1451. <[doi:10.1111/rssa.12712](https://rss.onlinelibrary.wiley.com/doi/abs/10.1111/rssa.12712)>.
