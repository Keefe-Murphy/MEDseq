#' MEDseq: Mixtures of Exponential-Distance Models with Covariates
#'
#' Fits MEDseq models: mixtures of Exponential-Distance models with gating covariates and sampling weights. Typically used for clustering categorical/longitudinal life-course sequences
#' @section Usage:
#' Fits _MEDseq_ models introduced by Murphy et al. (2019) <\href{https://arxiv.org/abs/1908.07963}{arXiv:1908.07963}>, i.e. fits mixtures of exponential-distance models for clustering longitudinal life-course sequence data via the EM/CEM algorithm. 
#' 
#' A family of parsimonious precision parameter constraints are accommodated. So too are sampling weights. Gating covariates can be supplied via formula interfaces.
#' 
#' The most important function in the \pkg{MEDseq} package is: \code{\link{MEDseq_fit}}, for fitting the models via EM/CEM. 
#' 
#' \code{\link{MEDseq_control}} allows supplying additional arguments which govern, among other things, controls on the initialisation of the allocations for the EM/CEM algorithm and the various model selection options. 
#' 
#' \code{\link{MEDseq_compare}} is provided for conducting model selection between different results from using different covariate combinations &/or initialisation strategies, etc. 
#' 
#' A dedicated plotting function exists for visualising various aspects of results, using new methods as well as some existing methods from the \pkg{TraMineR} package.
#' 
#' Finally, the package also contains two data sets: \code{\link{biofam}} and \code{\link{mvad}}.
#' 
#' @section Details:
#' \itemize{
#' \item{Type: }{Package}
#' \item{Package: }{MEDseq}
#' \item{Version: }{1.0.1}
#' \item{Date: }{2019-12-10 (this version), 2019-08-24 (original release)}
#' \item{Licence: }{GPL (>=2)}
#' }
#'
#' @section See Also:
#' Further details and examples are given in the associated vignette document:\cr
#' \code{vignette("MEDseq", package = "MEDseq")}
#'
#' @author
#' Keefe Murphy [aut, cre], Thomas Brendan Murphy [ctb], Raffaella Piccarreta [ctb], Isobel Claire Gormley [ctb]
#'
#' \strong{Maintainer}: Keefe Murphy - <\email{keefe.murphy@@ucd.ie}>
#' @references Keefe Murphy, T. Brendan Murphy, Raffaella Piccarreta, and I. Claire Gormley (2019). Clustering longitudinal life-course sequences using mixtures of exponential-distance models. \emph{To appear}. <\href{https://arxiv.org/abs/1908.07963}{arXiv:1908.07963}>.
#' @examples
#' \dontshow{library(TraMineR)}
#' # Load the MVAD data
#' data(mvad)
#' mvad$Location <- factor(apply(mvad[,5:9], 1L, function(x) 
#'                  which(x == "yes")), labels = colnames(mvad[,5:9]))
#' mvad          <- list(covariates = mvad[c(3:4,10:14,87)],
#'                       sequences = mvad[,15L:86L], 
#'                       weights = mvad[,2])
#' mvad.cov      <- mvad$covariates
#' states        <- c("EM", "FE", "HE", "JL", "SC", "TR")
#' labels        <- c("Employment", "FE", "HE", "Joblessness", "School", "Training")
#' mvad.seq      <- seqdef(mvad$sequences, states=states, labels=labels)
#' \donttest{                         
#' # Fit a range of unweighted models without covariates
#' # Only consider models with a noise component
#' # Supply some MEDseq_control() arguments
#' mod1          <- MEDseq_fit(mvad.seq, G=9:10, modtype=c("CCN", "CUN", "UCN", "UUN"),
#'                             algo="CEM", init.z="hc", criterion="asw")
#' 
#' # Fit a model with weights and gating covariates
#' # Drop the 1st time point which was used to define the weights
#' mvad.seq2     <- seqdef(mvad$sequences[,-1], states=states, labels=labels)
#' mod2          <- MEDseq_fit(mvad.seq2, G=10, modtype="UCN", weights=mvad$weights, 
#'                             gating=~ fmpr + gcse5eq + livboth, covars=mvad.cov)
#'                             
#' # Examine this model in greater detail
#' summary(mod2)
#' summary(mod2$gating)
#' plot(mod2, "clusters")}
#' @docType package
#' @keywords package
"_PACKAGE"

.onAttach <- function(lib, pkg) {
  version <- read.dcf(file.path(lib, pkg, "DESCRIPTION"), "Version")
  if(interactive()) {
    packageStartupMessage(paste("\nMixtures of Exponential-Distance Models with Covariates\n___  ___ ___________\n|  \\/  ||  ___|  _  \\\n| .  . || |__ | | | |___  ___  __ _\n| |\\/| ||  __|| | | / __|/ _ \\/ _` |\n| |  | || |___| |/ /\\__ \\  __/ (_| |\n\\_|  |_/\\____/|___/ |___/\\___|\\__, |\n                                 | |\n                                 |/       version", version, "\n"))                 
  } else   {
    packageStartupMessage("\nPackage ", sQuote("MEDseq"), " version ", version, ".\n")
  }
    packageStartupMessage(paste("See", sQuote("?MEDseq"), "for a brief guide to how to use this R package.\nSee", sQuote(paste0("citation(", dQuote("MEDseq"),")")) ,"for citing the package in publications.\nSee", sQuote("MEDseq_news()"), "for new features, recent changes, and bug fixes.\n"))
}
