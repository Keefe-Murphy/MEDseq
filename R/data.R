#' MVAD: Transition from school to work
#'
#' The data comes from a study by McVicar and Anyadike-Danes on transition from school to work. The data consist of static background characteristics and a time series sequence of 72 monthly labour market activities for each of a cohort of 712 individuals in the Status Zero Survey. The individuals were followed up from July 1993 to June 1999. The monthly states are recorded in columns 15 (\code{Jul.93}) to 86 (\code{Jun.99}).
#' @format A data frame containing 712 rows, 72 state variables, 1 id variable and 13 covariates.
#' @details States are:\cr
#' 
#'\tabular{ll}{
#'  \code{employment}  \tab (EM) \cr
#'  \code{FE}          \tab further education (FE) \cr
#'  \code{HE}          \tab higher education (HE) \cr
#'  \code{joblessness} \tab (JL) \cr
#'  \code{school}      \tab (SC) \cr
#'  \code{training}    \tab (TR) \cr
#'}

#'The data set contains also ids (\code{id}) and sample weights (\code{weights}) as well as the following binary covariates:\cr
#'\cr
#'\code{male}\cr
#'\code{catholic}\cr
#'\code{Belfast}, \code{N.Eastern}, \code{Southern}, \code{S.Eastern}, \code{Western} (location of school, one of five Education and Library Board areas in Northern Ireland)\cr
#'\code{Grammar} (type of secondary education, 1=grammar school)\cr
#'\code{funemp} (father's employment status at time of survey, 1=father unemployed)\cr
#'\code{gcse5eq} (qualifications gained by the end of compulsory education, 1=5+ GCSEs at grades A-C, or equivalent)\cr
#'\code{fmpr} (SOC code of father's current or most recent job at time of survey, 1=SOC1 (professional, managerial or related))\cr
#'\code{livboth} (living arrangements at time of first sweep of survey (June 1995), 1=living with both parents)
#' @references McVicar, D. (2000). Status 0 four years on: young people and social exclusion in Northern Ireland. \emph{Labour Market Bulletin}, 14, 114-119.
#' 
#' McVicar, D. and Anyadike-Danes, M. (2002). Predicting successful and unsuccessful transitions from school to work by using sequence methods. \emph{Journal of the Royal Statistical Society: Series A (Statistics in Society)}, 165(2): 317-334.
#' @note The first two months of the observation period coincide with summer holidays from school. Hence, documented examples throughout this package extract only the states in columns 17 to 86; i.e. sequences of length 70 from \code{Sep.93} to \code{Jun.99}.
#' @examples
#' data(mvad, package="MEDseq")
#' 
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
#' mvad.seq      <- seqdef(mvad$sequences[-c(1,2)], states=states, 
#'                         labels=labels, weights=mvad$weights)
#' @docType data
#' @importFrom TraMineR "seqdef"
#' @source McVicar and Anyadike-Danes (2002)
#' @keywords datasets
#' @usage data(mvad)
"mvad"

#' Family life states from the Swiss Household Panel biographical survey
#'
#' 2000 16 year-long family life sequences built from the retrospective biographical survey carried out by the Swiss Household Panel (SHP) in 2002.
#' @format A data frame with 2000 rows, 16 state variables, 1 id variable and 7 covariates and 2 weights variables.
#' @details The \emph{biofam} data set was constructed by Müller et al. (2007) from the data of the retrospective biographical survey carried out by the Swiss Household Panel (SHP) in 2002. 
#' 
#' The data set contains (in columns 10 to 25) sequences of family life states from age 15 to 30 (sequence length is 16) and a series of covariates. The sequences are a sample of 2000 sequences of those created from the SHP biographical survey. It includes only individuals who were at least 30 years old at the time of the survey. The \emph{biofam} data set describes family life courses of 2000 individuals born between 1909 and 1972. 
#' 
#' The states numbered from 0 to 7 are defined from the combination of five basic states, namely Living with parents (Parent), Left home (Left), Married (Marr), Having Children (Child), Divorced:\cr
#' \cr
#' 0 = "Parent" \cr
#' 1 = "Left" \cr
#' 2 = "Married" \cr
#' 3 = "Left+Marr" \cr
#' 4 = "Child" \cr
#' 5 = "Left+Child" \cr
#' 6 = "Left+Marr+Child" \cr
#' 7 = "Divorced" \cr
#' 
#' The covariates are: \cr
#' \tabular{ll}{
#' \code{sex}      \tab \cr
#' \code{birthyr}  \tab (birth year) \cr
#' \code{nat_1_02} \tab (first nationality) \cr
#' \code{plingu02} \tab (language of questionnaire) \cr
#' \code{p02r01}   \tab (religion) \cr
#' \code{p02r04}   \tab (religious participation) \cr
#' \code{cspfaj}   \tab (father's social status) \cr
#' \code{cspmoj}   \tab (mother's social status) \cr
#' }
#' Two additional weights variables are inserted for illustrative purpose ONLY (since \code{biofam} is a subsample of the original data, these weights are not adapted to the actual data): \cr
#' \tabular{ll}{
#' \code{wp00tbgp} \tab (weights inflating to the Swiss population) \cr
#' \code{wp00tbgs} \tab (weights respecting sample size) \cr
#'}
#' @references Müller, N. S., Studer, M. and Ritschard, G. (2007). Classification de parcours de vie à l'aide de l'optimal matching. In \emph{XIVe Rencontre de la Société francophone de classification (SFC 2007)}, Paris, 5-7 septembre 2007, pp. 157-160.
#' @examples
#' data(biofam, package="MEDseq")
#' 
#' biofam         <- list(covariates = biofam[2L:9L], sequences = biofam[10L:25L] + 1L)
#' biofam.cov     <- biofam$covariates[,colSums(is.na(biofam$covariates)) == 0]
#' biofam.seq     <- seqdef(biofam$sequences,
#'                         states = c("P", "L", "M", "L+M", "C", "L+C", "L+M+C", "D"),
#'                         labels = c("Parent", "Left", "Married", "Left+Marr", "Child", 
#'                                    "Left+Child", "Left+Marr+Child", "Divorced"))
#' biofam.cov$age <- 2002 - biofam.cov$birthyr
#' @docType data
#' @importFrom TraMineR "seqdef"
#' @source Swiss Household Panel \url{https://forscenter.ch/projects/swiss-household-panel/}
#' @keywords datasets
#' @usage data(biofam)
"biofam"