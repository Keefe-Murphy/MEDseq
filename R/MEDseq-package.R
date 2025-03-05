#' MEDseq: Mixtures of Exponential-Distance Models with Covariates
#'
#' Fits MEDseq models: mixtures of Exponential-Distance models with gating covariates and sampling weights. Typically used for clustering categorical/longitudinal life-course sequences.
#' @section Usage:
#' Fits _MEDseq_ models introduced by Murphy et al. (2021) <\doi{10.1111/rssa.12712}>, i.e. fits mixtures of exponential-distance models for clustering longitudinal life-course sequence data via the EM/CEM algorithm. 
#' 
#' A family of parsimonious precision parameter constraints are accommodated. So too are sampling weights. Gating covariates can be supplied via formula interfaces.
#' 
#' The most important function in the \pkg{MEDseq} package is: \code{\link{MEDseq_fit}}, for fitting the models via EM/CEM. This function requires the data to be in \code{"stslist"} format; the function \code{\link[TraMineR]{seqdef}} is conveniently reexported from the \pkg{TraMineR} package for this purpose.
#' 
#' \code{\link{MEDseq_control}} allows supplying additional arguments which govern, among other things, controls on the initialisation of the allocations for the EM/CEM algorithm and the various model selection options. 
#' 
#' \code{\link{MEDseq_compare}} is provided for conducting model selection between different results from using different covariate combinations &/or initialisation strategies, etc. 
#' 
#' \code{\link{MEDseq_stderr}} is provided for computing the standard errors of the coefficients for the covariates in the gating network.
#' 
#' A dedicated plotting function \code{\link{plot.MEDseq}} exists for visualising various aspects of the results, using new methods as well as some existing methods adapted from the \pkg{TraMineR} package.
#' 
#' Finally, the package also contains two data sets: \code{\link{biofam}} and \code{\link{mvad}}.
#' 
#' @details
#' \describe{
#' \item{Type: }{Package}
#' \item{Package: }{\pkg{MEDseq}}
#' \item{Version: }{1.4.2}
#' \item{Date: }{2025-03-10 (this version), 2019-08-24 (original release)}
#' \item{Licence: }{GPL (>= 3)}
#' }
#'
#' @seealso Useful links:
#' \itemize{
#' \item \url{https://cran.r-project.org/package=MEDseq}
#' \item Report bugs at \url{https://github.com/Keefe-Murphy/MEDseq}
#' }
#' 
#' Further details and examples are given in the associated vignette document:
#' \preformatted{vignette("MEDseq", package = "MEDseq")}
#' @author
#' Keefe Murphy [aut, cre], Thomas Brendan Murphy [ctb], Raffaella Piccarreta [ctb], Isobel Claire Gormley [ctb]
#'
#' \strong{Maintainer}: Keefe Murphy - <\email{keefe.murphy@@mu.ie}>
#' @references Murphy, K., Murphy, T. B., Piccarreta, R., and Gormley, I. C. (2021). Clustering longitudinal life-course sequences using mixtures of exponential-distance models. \emph{Journal of the Royal Statistical Society: Series A (Statistics in Society)}, 184(4): 1414-1451. <\doi{10.1111/rssa.12712}>.
#' @examples
#' # Load the MVAD data
#' data(mvad)
#' mvad$Location <- factor(apply(mvad[,5:9], 1L, function(x) 
#'                  which(x == "yes")), labels = colnames(mvad[,5:9]))
#' mvad          <- list(covariates = mvad[c(3:4,10:14,87)],
#'                       sequences = mvad[,15:86], 
#'                       weights = mvad[,2])
#' mvad.cov      <- mvad$covariates
#' 
#' # Create a state sequence object with the first two (summer) time points removed
#' states        <- c("EM", "FE", "HE", "JL", "SC", "TR")
#' labels        <- c("Employment", "Further Education", "Higher Education", 
#'                    "Joblessness", "School", "Training")
#' mvad.seq      <- seqdef(mvad$sequences[-c(1,2)], states=states, labels=labels)
#' \donttest{                         
#' # Fit a range of unweighted models without covariates
#' # Only consider models with a noise component
#' # Supply some MEDseq_control() arguments
#' mod1          <- MEDseq_fit(mvad.seq, G=9:10, modtype=c("CCN", "CUN", "UCN", "UUN"),
#'                             algo="CEM", init.z="kmodes", criterion="icl")
#' 
#' # Fit a model with weights and gating covariates
#' # Have the probability of noise-component membership be constant
#' mod2          <- MEDseq_fit(mvad.seq, G=11, modtype="UUN", weights=mvad$weights, 
#'                             gating=~ gcse5eq, covars=mvad.cov, noise.gate=FALSE)
#'                             
#' # Examine this model and its gating network
#' summary(mod2, network=TRUE)
#' plot(mod2, "clusters")}
#' @docType package
#' @keywords package
"_PACKAGE"

.onAttach <- function(lib, pkg) {
  path    <- file.path(lib, pkg, "DESCRIPTION")
  version <- read.dcf(path, "Version")
  name    <- read.dcf(path, "Package")
  if(interactive()) {
    packageStartupMessage(paste("\nMixtures of Exponential-Distance Models with Covariates\n __  __ _____ _____\t\t\t  version", version, "\n|  \\/  |  ___|  __ \\\n| \\  / | |__ | |  \\ |___  ___  __ _\n| |\\/| |  __|| |  | / __|/ _ \\/ _` |\n| |  | | |___| |__/ \\__ \\  __/ (_| |\n|_|  |_|_____|_____/|___/\\___|\\__, |\n                                 | |\n                                 |/\n"))                 
  } else   {
    packageStartupMessage("\nPackage ", sQuote(name), " version ", version, ".\n")
  }
    packageStartupMessage(paste("See", sQuote("?MEDseq"), "to see a brief guide to how to use this R package.\nSee", sQuote(paste0("citation(", dQuote(name),")")) ,"for citing the package in publications.\nSee", sQuote("MEDseq_news()"), "to see new features, changes, and bug fixes.\n"))
  if(interactive() &&
     name %in% utils::old.packages()[,1L]) {
    packageStartupMessage("\n !!! A newer version of this package is available from CRAN !!!")
  }
}
