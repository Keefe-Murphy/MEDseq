#' Compute the Density-based Silhouette
#' 
#' Computes the Density-based Silhouette for a 'soft' clustering assignment matrix.
#' @param z A numeric matrix such that rows correspond to observations, columns correspond to clusters, and rows sum to \code{1}.
#' @param ztol A small (single, numeric, non-negative) tolerance parameter governing whether small assignment probabilities are treated instead as crisp assignments. Defaults to \code{1E-100}.
#' @param weights An optional numeric vector giving observation-specific weights for computing the (weighted) mean/median DBS (see \code{summ}).
#' @param summ A single character string indicating whether the (possibly weighted) \code{"mean"} (the default) or \code{"median"} DBS should be computed.
#' @param clusters Optional/experimental argument for giving the indicator labels of the cluster assignments. Defaults to the MAP assignment derived from \code{z} when not supplied. Note that actually supplying the MAP assignment here is slightly less efficient than the \code{NULL} default and \strong{not} advised.
#' @param ... Catches unused arguments.
#'
#' @return A list with the following elements:
#' \describe{
#' \item{\code{silvals}}{A matrix where each row contains the cluster to which each observation belongs in the first column and the observation-specific DBS width in the second column.}
#' \item{\code{msw}}{Depending on the value of \code{summ}, either the mean or median DBS width.}
#' \item{\code{wmsw}}{Depending on the value of \code{summ}, either the weighted mean or weighted median DBS width.}}
#' @note When calling \code{\link{MEDseq_fit}}, the \code{summ} argument can be passed via the \code{...} construct, in which case it governs both the \code{dbs} and \code{asw} criteria.
#' @importFrom matrixStats "rowSums2" "weightedMedian" "weightedMean"
#' @importFrom TraMineR "seqdef"
#' @references Menardi, G. (2011). Density-based Silhouette diagnostics for clustering methods. \emph{Statistics and Computing}, 21(3): 295-308.
#' @export 
#' @seealso \code{\link{MEDseq_fit}}
#' @author Keefe Murphy - <\email{keefe.murphy@@mu.ie}>
#' @keywords utility
#' @usage
#' dbs(z,
#'     ztol = 1E-100,
#'     weights = NULL,
#'     summ = c("mean", "median"),
#'     clusters = NULL,
#'     ...)
#' @examples
#' # Generate a toy z matrix
#' z <- abs(matrix(rnorm(50), ncol=2))
#' z <- z/rowSums(z)
#' 
#' # Return the median DBS width
#' dbs(z, summ="median")$msw
#' 
#' # For real sequence data
#' data(mvad)
#' \donttest{
#' mod <- MEDseq_fit(seqdef(mvad[,17:86]), G=11, modtype="UUN", weights=mvad$weight)
#' 
#' dbs(mod$z, weights=mvad$weight)}
dbs               <- function(z, ztol = 1E-100, weights = NULL, summ = c("mean", "median"), clusters = NULL, ...) {
  if(any(!is.matrix(z), !is.numeric(z)) ||
     ncol(z)      <= 1     ||
     nrow(z)      <= 1)          stop("'z' must be a numeric matrix with 2 or more columns & 2 or more rows", call.=FALSE)
  z               <- .renorm_z(z)
  if(length(ztol)  > 1     ||
     !is.numeric(ztol)     ||
     ztol          < 0)          stop("Invalid 'ztol'", call.=FALSE)
  if(!is.null(weights))     {
    N             <- nrow(z)
   if(!is.numeric(weights) ||
      length(weights)      != N) stop(paste0("'weights' must be a numeric vector of length N=", N), call.=FALSE)
   if(any(weights < 0)     || 
      any(!is.finite(weights)))  stop("'weights' must be non-negative and finite", call.=FALSE)
  }
  if(!missing(summ)        && 
    (length(summ)  > 1     ||
     !is.character(summ)))       stop("'summ' must be a single character string",  call.=FALSE)
  summ            <- match.arg(summ)
  if(is.null(clusters))     {
    MAP           <- if(any(names(list(...)) == "MAP")) list(...)$MAP else max.col(z)
    zz            <- matrix(z[order(row(z), -z)], nrow(z), byrow=TRUE)
    l2            <- log(zz[,2L])
    zz            <- log(zz[,1L])    - l2
  } else           {
    if(length(clusters)    != nrow(z)    ||
       any(clusters        !=
       floor(clusters)))         stop("Invalid 'clusters'", call.=FALSE)
    cmax          <- length(unique(clusters))
    if(ncol(z)    != cmax)  {    warning("Number of groups in 'clusters' differs from the number of columns in 'z'\n", call.=FALSE)
      if(ncol(z)   < cmax)  {
        z         <- cbind(z, matrix(0L, nrow=nrow(z), ncol=cmax - ncol(z)))
      }
    }
    MAP           <- clusters
    ordered       <- apply(z, 1L, order, decreasing = TRUE)
    indX          <- .misclass(clusters, ordered[1L,])$misclassified
    ind           <- setdiff(seq_len(nrow(z)), indX)
    l2            <-
    zz            <- integer(nrow(z))
    l2[indX]      <- log(diag(z[indX,ordered[1L,indX]]))
    zz[indX]      <- log(diag(z[indX,clusters]))       - l2[indX]
    l2[ind]       <- log(diag(z[ind,ordered[2L,ind]]))
    zz[ind]       <- log(diag(z[ind,ordered[1L,ind]])) - l2[ind]
  }
  zz.inf          <- is.infinite(zz) | l2 < log(ztol)
  ds              <- zz/max(1L, abs(zz[!zz.inf]))
  ds[zz.inf & zz  >= 0]    <- 1L
  ds[is.nan(ds)]  <- 0L
  DS              <- cbind(cluster=MAP, dbs_width=ds)
  class(DS)       <- "MEDsil"
  msw             <- switch(EXPR=summ, median=stats::median(ds), mean=mean(ds))
  dbs_res         <- list(silvals = DS, msw = msw, wmsw = ifelse(is.null(weights), msw, 
                          switch(EXPR=summ, median=weightedMedian(ds, weights), mean=weightedMean(ds, weights))))
  attr(dbs_res, "summ")    <- summ
    return(dbs_res)
}

#' Extract results from a MEDseq model
#'
#' Utility function for extracting results of submodels from \code{"MEDseq"} objects when a range of models were run via \code{\link{MEDseq_fit}}.
#' @param x An object of class \code{"MEDseq"} generated by \code{\link{MEDseq_fit}} or an object of class \code{"MEDseqCompare"} generated by \code{\link{MEDseq_compare}}.
#' @param what A character string indicating the desired results to extract.
#' @param rank A number indicating what \code{rank} model results should be extracted from, where the \code{rank} is determined by \code{criterion}. Defaults to \code{1}, i.e. the best model.
#' @param criterion The \code{criterion} used to determine the ranking. Defaults to \code{"bic"}.
#' @param G Optional argument giving the number of components in the model for which results are desired. Can be supplied with or without also specifying \code{modtype}.
#' @param modtype Optional argument giving the desired model type for which results are desired. Can be supplied with or without also specifying \code{G}.
#' @param noise A logical indicating whether models with a noise component should be considered. Defaults to \code{TRUE}.
#' @param ... Catches unused arguments.
#'
#' @return The desired results extracted from the \code{MEDseq} model.
#' @details The arguments \code{rank} and \code{criterion} are invoked when one or more of the arguments \code{G} and \code{modtype} are missing. Thus, supplying \code{G} and \code{modtype} allows \code{rank} and \code{criterion} to be bypassed entirely.
#' @note Arguments to this function can be supplied to \code{\link{plot.MEDseq}} via the \code{...} construct.
#' @export
#' @importFrom TraMineR "seqdef"
#' @seealso \code{\link{MEDseq_fit}}, \code{\link{plot.MEDseq}}
#' @author Keefe Murphy - <\email{keefe.murphy@@mu.ie}>
#' @keywords utility
#' @usage
#' get_MEDseq_results(x,
#'                    what = c("z", "MAP", "DBS", "ASW"), 
#'                    rank = 1L, 
#'                    criterion = c("bic", "icl", "aic", "dbs", 
#'                                  "asw", "cv", "nec", "loglik"), 
#'                    G = NULL, 
#'                    modtype = NULL, 
#'                    noise = TRUE, 
#'                    ...)
#' @examples
#' \donttest{data(biofam)
#' # mod <- MEDseq_fit(seqdef(biofam[10:25] + 1L), G=9:10)
#' 
#' # Extract the MAP clustering of the best 9-cluster model according to the asw criterion
#' # get_MEDseq_results(mod, what="MAP", G=9, criterion="asw")
#' 
#' # Extract the DBS values of the best UUN model according to the dbs criterion
#' # get_MEDseq_results(mod, what="DBS", modtype="UUN", criterion="dbs")
#' 
#' # Plot the DBS values of this same model, by passing get_MEDseq_results arguments through plot
#' # plot(mod, type="dbsvals", modtype="UUN", criterion="dbs")}
get_MEDseq_results            <- function(x, what = c("z", "MAP", "DBS", "ASW"), rank = 1L, criterion = c("bic", "icl", "aic", "dbs", "asw", "cv", "nec", "loglik"), G = NULL, modtype = NULL, noise = TRUE, ...) {
    UseMethod("get_MEDseq_results")
}

#' @method get_MEDseq_results MEDseq
#' @export
get_MEDseq_results.MEDseq     <- function(x, what = c("z", "MAP", "DBS", "ASW"), rank = 1L, criterion = c("bic", "icl", "aic", "dbs", "asw", "cv", "nec", "loglik"), G = NULL, modtype = NULL, noise = TRUE, ...) {
  x               <- if(inherits(x, "MEDseqCompare")) x$optimal else x
  if(!missing(what)           && 
     (length(what) > 1        ||
     !is.character(what)))       stop("'what' must be a single character string",      call.=FALSE)
  what            <- match.arg(what)
  minG            <- 1L  + is.element(what, c("DBS", "ASW"))
  if(!(missing(G) -> m.G)     &&
    (length(G)    != 1        ||
     !is.numeric(G)           ||
     (G < minG    || floor(G) != G))) {
    if(is.element(what, 
       c("dbs", "asw"))) {       stop(paste0("'G' must be a single integer > 1 when 'what'=", what), call.=FALSE)
    } else                       stop("'G' must be a single integer >= 1",             call.=FALSE)
  }  
  if(!(missing(modtype) ->
       m.M) &&
    (length(modtype) > 1      ||
     !is.character(modtype)))    stop("'modtype' must be a single character string",   call.=FALSE)
  if(length(noise) > 1        ||
     !is.logical(noise))         stop("'noise' must be a single logical indicator",    call.=FALSE)
  if(!missing(criterion)      ||
     !missing(rank)           ||
     any(m.G, m.M))            {
    if(!missing(criterion)    && 
      (length(criterion)  > 1 ||
      !is.character(criterion))) stop("'criterion' must be a single character string", call.=FALSE)
    criterion     <- match.arg(criterion)
    if((criterion == "nec"    ||
       criterion  == "dbs"    ||
       criterion  == "asw")   &&
      !m.G  &&  G == 1)          stop(paste0("Can't select based on the ", toupper(criterion), " criterion when G=1"), call.=FALSE)
    if(criterion  == "dbs"    &&
       attr(x, "Algo") == "CEM") stop(paste0("Can't select based on the ", toupper(criterion), " criterion as the CEM algorithm was used to fit the model"), call.=FALSE)
    if(criterion  == "cv"     &&
       !attr(x, "CV"))           stop("Can't select based on the CV criterion as cross-validated likelihood wasn't performed", call.=FALSE)
    tmp           <- switch(EXPR=criterion, bic=x$BIC, icl=x$ICL, aic=x$AIC, cv=x$CV, nec=x$NEC, dbs=x$DBS, asw=x$ASW, loglik=x$LOGLIK)
    if(is.null(tmp))             stop(paste0("Can't select based on the ", toupper(criterion), " criterion when all models contain only 1 component"),       call.=FALSE)
    if(!noise)     {
      n.names     <- colnames(tmp) %in% c("CC", "UC", "CU", "UU")
      if(!any(n.names))          stop("No non-noise models to select", call.=FALSE)
      tmp         <- tmp[,n.names, drop=FALSE]
    }
    if(!m.G)       {
      Gallow      <- as.numeric(rownames(tmp))
      if(!(G %in% Gallow))       stop("Invalid 'G'",       call.=FALSE)
      tmp         <- tmp[as.character(G),, drop=FALSE]
    }
    if(!m.M)       {
      Mallow      <- colnames(tmp)
      if(!(modtype %in% Mallow)) stop("Invalid 'modtype'", call.=FALSE)
      tmp         <- tmp[,modtype, drop=FALSE]
    }
    class(tmp)    <- "MEDcriterion"
    attr(tmp, "Criterion")  <- toupper(criterion)
    if(length(rank) > 1     ||
       !is.numeric(rank)    ||
       rank != floor(rank)  ||
       rank <= 0  ||
       rank  > sum(!is.na(tmp))) stop("Invalid 'rank'",    call.=FALSE)
    best          <- strsplit(names(.pick_MEDCrit(tmp, pick=rank)$crits[rank]), ",")[[1L]]
    modtype       <- best[1L]
    G             <- as.numeric(best[2L])
  }
  switch(EXPR=what, 
         DBS=, ASW=         {
    SILS          <- switch(EXPR=what, DBS=x$DBSvals, ASW=x$ASWvals)
    summ          <- attr(SILS, "Summ")
    if(!(G  %in%
       as.numeric(names(SILS)))) stop("Invalid 'G' value", call.=FALSE)
    S             <- SILS[[as.character(G)]]
    s.ind         <- if(noise) names(S) else names(S)[names(S) %in% c("CC", "UC", "CU", "UU")]
    if(!(modtype %in% s.ind))    stop("Invalid 'modtype'", call.=FALSE)
    res           <- S[[modtype]]
    if(anyNA(res))               message("Selected model didn't converge: no silhouettes available\n")
    attr(res, "Summ")      <- summ
  }, {
    ZS            <- x$ZS
    if(!(G  %in%
       as.numeric(names(ZS))))   stop("Invalid 'G' value", call.=FALSE)
    Z             <- ZS[[as.character(G)]]
    z.ind         <- if(noise) names(Z) else names(Z)[names(Z) %in% c("CC", "UC", "CU", "UU")]
    if(!(modtype %in% z.ind))    stop("Invalid 'modtype'", call.=FALSE)
    res           <- Z[[modtype]]
    if(anyNA(res))               message("Selected model didn't converge: no partition available\n")
    noise         <- is.element(modtype, c("CCN", "UCN", "CUN", "UUN"))
    if(what == "MAP")       {
      MAP         <- max.col(res)
      res         <- if(noise) replace(MAP, MAP == G, 0L) else MAP
    }
  })
  attr(res, "G")           <- G
  attr(res, "ModelType")   <- modtype
  attr(res, "Noise")       <- noise
    return(res)
}

#' Choose the best MEDseq model
#'
#' Takes one or more sets of \code{"MEDseq"} models fitted by \code{\link{MEDseq_fit}} and ranks them according to a specified model selection criterion. It's possible to respect the internal ranking within each set of models, or to discard models within each set which were already deemed sub-optimal. This function can help with model selection via exhaustive or stepwise searches.
#' @param ... One or more objects of class \code{"MEDseq"} outputted by \code{\link{MEDseq_fit}}. All models must have been fit to the same data set. A single \emph{named} list of such objects can also be supplied. Additionally, objects of class \code{"MEDseqCompare"} outputted by this very function can also be supplied here.
#' 
#' This argument is only relevant for the \code{\link{MEDseq_compare}} function and will be ignored for the associated \code{print} function.
#' @param criterion The criterion used to determine the ranking. Defaults to \code{"bic"}.
#' @param pick The (integer) number of models to be ranked and compared. Defaults to \code{10L}. Will be constrained by the number of models within the \code{"MEDseq"} objects supplied via \code{...} if \code{optimal.only} is \code{FALSE}, otherwise constrained simply by the number of \code{"MEDseq"} objects supplied. Setting \code{pick=Inf} is a valid way to select all models.
#' @param optimal.only Logical indicating whether to only rank models already deemed optimal within each \code{"MEDeq"} object (\code{TRUE}), or to allow models which were deemed suboptimal enter the final ranking (\code{FALSE}, the default). See \code{Details}.
#' @param x,index,rerank,digits,maxi Arguments required for the associated \code{print} function:
#' \describe{
#' \item{\code{x}}{An object of class \code{"MEDseqCompare"} resulting from a call to \code{\link{MEDseq_compare}}.}
#' \item{\code{index}}{A logical or numeric vector giving the indices of the rows of the table of ranked models to print. This defaults to the full set of ranked models. It can be useful when the table of ranked models is large to examine a subset via this \code{index} argument, for display purposes. See \code{rerank}.}
#' \item{\code{rerank}}{A logical indicating whether the ranks should be recomputed when subsetting using \code{index}. Defaults to \code{FALSE}.}
#' \item{\code{digits}}{The number of decimal places to round model selection criteria to (defaults to \code{3}).}
#' \item{\code{maxi}}{A number specifying the maximum number of rows/models to print. Defaults to \code{length(index)}.}}
#' @note The \code{criterion} argument here need not comply with the criterion used for model selection within each \code{"MEDseq"} object, but be aware that a mismatch in terms of \code{criterion} \emph{may} require the optimal model to be re-fit in order to be extracted, thereby slowing down \code{\link{MEDseq_compare}}.
#' 
#' If random starts had been used via \code{init.z="random"} the \code{optimal} model may not necessarily correspond to the highest-ranking model in the presence of a criterion mismatch, due to the randomness of the initialisation. 
#'
#' A dedicated \code{print} function exists for objects of class \code{"MEDseqCompare"} and \code{\link{plot.MEDseq}} can also be called on objects of class \code{"MEDseqCompare"}.
#' @return A list of class \code{"MEDseqCompare"}, for which a dedicated print function exists, containing the following elements (each of length \code{pick}, and ranked according to \code{criterion}, where appropriate):
#' \item{\code{data}}{The name of the data set to which the models were fitted.}
#' \item{\code{optimal}}{The single optimal model (an object of class \code{"MEDseq"}) among those supplied, according to the chosen \code{criterion}.}
#' \item{\code{pick}}{The final number of ranked models. May be different (i.e. less than) the supplied \code{pick} value.}
#' \item{\code{MEDNames}}{The names of the supplied \code{"MEDseq"} objects.}
#' \item{\code{modelNames}}{The MEDseq model names (denoting the constraints or lack thereof on the precision parameters).}
#' \item{\code{G}}{The optimal numbers of components.}
#' \item{\code{df}}{The numbers of estimated parameters.}
#' \item{\code{iters}}{The numbers of EM/CEM iterations.}
#' \item{\code{bic}}{BIC values, ranked according to \code{criterion}.}
#' \item{\code{icl}}{ICL values, ranked according to \code{criterion}.}
#' \item{\code{aic}}{AIC values, ranked according to \code{criterion}.}
#' \item{\code{dbs}}{(Weighted) mean/median DBS values, ranked according to \code{criterion}.}
#' \item{\code{asw}}{(Weighted) mean/median ASW values, ranked according to \code{criterion}.}
#' \item{\code{cv}}{Cross-validated log-likelihood values, ranked according to \code{criterion}.}
#' \item{\code{nec}}{NEC values, ranked according to \code{criterion}.}
#' \item{\code{loglik}}{Maximal log-likelihood values, ranked according to \code{criterion}.}
#' \item{\code{gating}}{The gating formulas.}
#' \item{\code{algo}}{The algorithm used for fitting the model - either \code{"EM"}, \code{"CEM"}, \code{"cemEM"}.}
#' \item{\code{equalPro}}{Logical indicating whether mixing proportions were constrained to be equal across components.}
#' \item{\code{opti}}{The method used for estimating the central sequence(s).}
#' \item{\code{weights}}{Logical indicating whether the given model was fitted with sampling weights.}
#' \item{\code{noise}}{Logical indicating the presence/absence of a noise component. Only displayed if at least one of the compared models has a noise component.}
#' \item{\code{noise.gate}}{Logical indicating whether gating covariates were allowed to influence the noise component's mixing proportion. Only printed for models with a noise component, when at least one of the compared models has gating covariates.}
#' \item{\code{equalNoise}}{Logical indicating whether the mixing proportion of the noise component for \code{equalPro} models is also equal (\code{TRUE}) or estimated (\code{FALSE}).}
#' @details The purpose of this function is to conduct model selection on \code{"MEDseq"} objects, fit to the same data set, with different combinations of gating network covariates or different initialisation settings.
#'
#' Model selection will have already been performed in terms of choosing the optimal number of components and MEDseq model type within each supplied set of results, but \code{\link{MEDseq_compare}} will respect the internal ranking of models when producing the final ranking if \code{optimal.only} is \code{FALSE}: otherwise only those models already deemed optimal within each \code{"MEDseq"} object will be ranked.
#'
#' As such if two sets of results are supplied when \code{optimal.only} is \code{FALSE}, the 1st, 2nd, and 3rd best models could all belong to the first set of results, meaning a model deemed suboptimal according to one set of covariates could be superior to one deemed optimal under another set of covariates.
#' @export
#' @keywords clustering main
#' @importFrom TraMineR "seqdef"
#' @author Keefe Murphy - <\email{keefe.murphy@@mu.ie}>
#' @references Murphy, K., Murphy, T. B., Piccarreta, R., and Gormley, I. C. (2021). Clustering longitudinal life-course sequences using mixtures of exponential-distance models. \emph{Journal of the Royal Statistical Society: Series A (Statistics in Society)}, 184(4): 1414-1451. <\href{https://rss.onlinelibrary.wiley.com/doi/abs/10.1111/rssa.12712}{doi:10.1111/rssa.12712}>.
#' @seealso \code{\link{MEDseq_fit}}, \code{\link{plot.MEDseq}}
#' @usage
#' MEDseq_compare(...,
#'                criterion = c("bic", "icl", "aic", 
#'                              "dbs", "asw", "cv", "nec"),
#'                pick = 10L,
#'                optimal.only = FALSE)
#' @examples
#' data(biofam)
#' seqs <- seqdef(biofam[10:25] + 1L,
#'                states = c("P", "L", "M", "L+M", "C", 
#'                           "L+C", "L+M+C", "D"))
#' covs <- cbind(biofam[2:3], age=2002 - biofam$birthyr)
#' \donttest{ 
#' # Fit a range of models
#' # m1   <- MEDseq_fit(seqs, G=9:10)
#' # m2   <- MEDseq_fit(seqs, G=9:10, gating=~sex,       covars=covs, noise.gate=FALSE)
#' # m3   <- MEDseq_fit(seqs, G=9:10, gating=~age,       covars=covs, noise.gate=FALSE)
#' # m4   <- MEDseq_fit(seqs, G=9:10, gating=~sex + age, covars=covs, noise.gate=FALSE)
#' 
#' # Rank only the optimal models (according to the dbs criterion)
#' # Examine the best model in more detail
#' # (comp <- MEDseq_compare(m1, m2, m3, m4, criterion="dbs", optimal.only=TRUE))
#' # (best <- comp$optimal)
#' # (summ <- summary(best, parameters=TRUE))
#' 
#' # Examine all models visited, including those already deemed suboptimal
#' # Only print models with gating covariates & 10 components
#' # comp2 <- MEDseq_compare(comp, m1, m2, m3, m4, criterion="dbs", pick=Inf)
#' # print(comp2, index=comp2$gating != "None" & comp2$G == 10)}
MEDseq_compare    <- function(..., criterion = c("bic", "icl", "aic", "dbs", "asw", "cv", "nec"), pick = 10L, optimal.only = FALSE) {
  crit.miss       <- missing(criterion)
  if(!missing(criterion)   && (length(criterion) > 1 ||
     !is.character(criterion)))  stop("'criterion' must be a single character string", call.=FALSE)
  criterion       <- match.arg(criterion)
  num.miss        <- missing(pick)
  opt.miss        <- missing(optimal.only)
  if(length(pick) != 1     ||
     !is.numeric(pick))          stop("'pick' must be a single number", call.=FALSE)
  if(floor(pick)  != pick  ||
     pick          < 1)          stop("'pick' must be a strictly positive integer", call.=FALSE)
  if(length(optimal.only)   > 1 ||
     !is.logical(optimal.only))  stop("'optimal.only' must be a single logical indicator", call.=FALSE)
  call            <- match.call(expand.dots=TRUE)[-1L]
  call            <- if(crit.miss) call else call[-which(names(call) == "criterion")]
  call            <- if(num.miss)  call else call[-which(names(call) == "pick")]
  call            <- if(opt.miss)  call else call[-which(names(call) == "optimal.only")]
  len.call        <- length(as.list(call))
  if(len.call     == 1     && inherits(..., "list") && !inherits(..., "MEDseq")) {
    dots          <- as.list(...)
    mod.names     <- unique(names(dots))
    comparison    <- vapply(dots, inherits,  logical(1L), "MEDseqCompare")
    dat.name      <- if(any(comparison)) dots[[1L]]$data
    dots[comparison]       <- sapply(dots[comparison], "[", "optimal")
    MEDs          <- dots[mod.names]
    if(is.null(mod.names))       stop("When supplying models as a list, every element of the list must be named", call.=FALSE)
  } else           {
    dots          <- list(...)
    mod.names     <- vapply(call, deparse, character(1L))
    comparison    <- vapply(dots, inherits,  logical(1L), "MEDseqCompare")
    dat.name      <- if(any(comparison)) dots[[1L]]$data
    dots[comparison]       <- sapply(dots[comparison], "[", "optimal")
    MEDs          <- stats::setNames(dots, mod.names)
    mod.names     <- unique(mod.names)
    MEDs          <- MEDs[mod.names]
  }
  Mclass          <- vapply(MEDs, class,          character(1L))
  if(any(Mclass   != "MEDseq"))  stop("All models must be of class 'MEDseq'!", call.=FALSE)
  data            <- lapply(MEDs, "[[", "data")
  data            <- lapply(data, unname)
  if(length(data)  > 1     && 
     !.unique_list(data))        stop("All models being compared must have been fit to the same data set!", call.=FALSE)
  dat.name        <- if(is.null(dat.name)) deparse(MEDs[[1L]]$call$seqs) else dat.name
  gate.x          <- lapply(MEDs,   "[[", "gating")
  algo            <- vapply(MEDs,   attr, character(1L), "Algo")
  opti            <- vapply(MEDs,   attr, character(1L), "Opti")
  equalNoise      <- vapply(MEDs,   attr, logical(1L),   "EqualNoise")
  equalPro        <- vapply(MEDs,   attr, logical(1L),   "EqualPro")
  noise.gate      <- vapply(MEDs,   attr, logical(1L),   "NoiseGate")
  weights         <- vapply(MEDs,   attr, logical(1L),   "Weighted")
  gating          <- lapply(gate.x, attr, "Formula")
  BICs            <- lapply(MEDs, "[[", "BIC")
  ICLs            <- lapply(MEDs, "[[", "ICL")
  AICs            <- lapply(MEDs, "[[", "AIC")
  LLxs            <- lapply(MEDs, "[[", "LOGLIK")
  DFxs            <- lapply(MEDs, "[[", "DF")
  ITxs            <- lapply(MEDs, "[[", "ITERS")
  CVs             <- lapply(MEDs, "[[", "CV")
  NECs            <- lapply(MEDs, "[[", "NEC")
  DBSs            <- lapply(MEDs, "[[", "DBS")
  ASWs            <- lapply(MEDs, "[[", "ASW")
  dbsnull         <- vapply(DBSs, is.null, logical(1L))
  aswnull         <- vapply(ASWs, is.null, logical(1L))
  cvnull          <- vapply(CVs,  is.null, logical(1L))
  necnull         <- vapply(NECs, is.null, logical(1L))
  if(all(dbsnull) &&             
     criterion    == "dbs")      stop(paste0("'criterion' cannot be 'dbs' when all models being compared ", ifelse(all(algo == "CEM"), "were fitted via CEM", "contain only 1 component")), call.=FALSE)
  if(all(aswnull) &&             
     criterion    == "asw")      stop("'criterion' cannot be 'asw' when all models being compared contain only 1 component", call.=FALSE)
  if(all(cvnull)  && 
     criterion    == "cv")       stop("'criterion' cannot be 'cv' when cross-validation was not performed for any of the supplied models", call.=FALSE)
  if(all(necnull) &&             
     criterion    == "nec")      stop("'criterion' cannot be 'nec' when all models being compared contain only 1 component", call.=FALSE)
  choice          <- max(lengths(BICs))
  bics            <- lapply(BICs, function(x) .pick_MEDCrit(x, choice)$crits)
  icls            <- lapply(ICLs, function(x) .pick_MEDCrit(x, choice)$crits)
  aics            <- lapply(AICs, function(x) .pick_MEDCrit(x, choice)$crits)
  llxs            <- lapply(LLxs, function(x) .pick_MEDCrit(x, choice)$crits)
  dfxs            <- lapply(DFxs, function(x) .pick_MEDCrit(x, choice)$crits)
  itxs            <- lapply(ITxs, function(x) .pick_MEDCrit(x, choice)$crits)
  dbss            <- lapply(DBSs, function(x) if(!is.null(x)) .pick_MEDCrit(x, choice)$crits)[!dbsnull]
  asws            <- lapply(ASWs, function(x) if(!is.null(x)) .pick_MEDCrit(x, choice)$crits)[!aswnull]
  cvs             <- lapply(CVs,  function(x) if(!is.null(x)) .pick_MEDCrit(x, choice)$crits)[!cvnull]
  necs            <- lapply(NECs, function(x) if(!is.null(x)) .pick_MEDCrit(x, choice)$crits)[!necnull]
  if(optimal.only) {
    opt.names     <- names(.crits_names(lapply(switch(EXPR=criterion, bic=bics, icl=icls, aic=aics, cv=cvs, nec=necs, dbs=dbss, asw=asws), "[", 1L)))
  }
  bics            <- .crits_names(bics)
  icls            <- .crits_names(icls)
  aics            <- .crits_names(aics)
  llxs            <- .crits_names(llxs)
  dfxs            <- .crits_names(dfxs)
  itxs            <- .crits_names(itxs)
  dbss            <- .crits_names(dbss)
  asws            <- .crits_names(asws)
  cvs             <- .crits_names(cvs)
  necs            <- .crits_names(necs)
  if(criterion    == "cv"  &&
    (length(cvs)  != 
     length(bics)))              warning("Discarding models for which the CV criterion was not computed\n",  call.=FALSE, immediate.=TRUE)
  if(criterion    == "nec" &&
    (length(necs) != 
     length(bics)))              warning("Discarding models for which the NEC criterion was not computed\n", call.=FALSE, immediate.=TRUE)
  if(criterion    == "dbs" &&
    (length(dbss) != 
     length(bics)))              warning("Discarding models for which the DBS criterion was not computed\n", call.=FALSE, immediate.=TRUE)
  if(criterion    == "asw" &&
    (length(asws) != 
     length(bics)))              warning("Discarding models for which the ASW criterion was not computed\n", call.=FALSE, immediate.=TRUE)
  if(optimal.only) {
    bics          <- bics[names(bics) %in% opt.names]
    icls          <- icls[names(icls) %in% opt.names]
    aics          <- aics[names(aics) %in% opt.names]
    llxs          <- llxs[names(llxs) %in% opt.names]
    dfxs          <- dfxs[names(dfxs) %in% opt.names]
    itxs          <- itxs[names(itxs) %in% opt.names]
    dbss          <- dbss[names(dbss) %in% opt.names]
    asws          <- asws[names(asws) %in% opt.names]
    cvs           <- cvs[names(cvs)   %in% opt.names]
    necs          <- necs[names(necs) %in% opt.names]
  }
  crits           <- switch(EXPR=criterion, bic=bics, icl=icls, aic=aics, dbs=dbss, asw=asws, cv=cvs, nec=necs)
  pick            <- min(pick, length(crits))
  max.crits       <- sort(crits, decreasing=criterion != "nec")[seq_len(pick)]
  if(length(unique(max.crits))  < pick) {
    ties          <- max.crits == max.crits[1L]
    if(any(ties[-1L]))      {     warning(paste0("Ties for the optimal model exist according to the '", criterion, "' criterion: choosing the most parsimonious model\n"), call.=FALSE, immediate.=TRUE)
      df.ties     <- dfxs[names(max.crits)][which(ties)]
      max.crits[ties]      <- max.crits[order(df.ties)]
      if(any((df.ties      == df.ties[1L])[-1L])) {
        max.crits[ties]    <- max.crits[order(as.numeric(gsub(".*,", "", names(max.crits[ties]))))]
      }
    } else                       warning(paste0("Ties exist according to the '", criterion, "' criterion\n"), call.=FALSE, immediate.=TRUE)
  }
  max.names       <- names(max.crits)
  crit.names      <- gsub("\\|.*", "",          max.names)
  G               <- as.numeric(gsub(".*,", "", max.names))
  gating          <- unname(unlist(gating[crit.names]))
  modelNames      <- gsub(",.*", "", gsub(".*\\|", "", max.names))
  best.model      <- MEDs[[crit.names[1L]]]
  if(best.model$modtype != modelNames[1L] || best.model$G != G[1L]) {
    message("Re-fitting optimal model due to mismatched 'criterion'...\n\n")
    old.call    <- best.model$call
    old.call    <- c(as.list(old.call)[1L], list(criterion=criterion), as.list(old.call)[-1L])
    old.call    <- as.call(old.call[!duplicated(names(old.call))])
    if(!is.null(old.call$init.z)    &&
       old.call$init.z  == 
       "random")                 warning("Optimal model may differ slightly due to criterion mismatch and random starts used in the initialisation:\nPrinted output intended only as a guide", call.=FALSE, immediate.=TRUE)
    best.call   <- c(list(data=best.model$data, modtype=modelNames[1L], G=G[1L], criterion="bic", verbose=FALSE, do.cv=FALSE, do.nec=FALSE), as.list(old.call[-1L]))
    best.mod    <- try(do.call(MEDseq_fit, best.call[!duplicated(names(best.call))]), silent=TRUE)
    if(!inherits(best.model, "try-error")) {
      best.model$call               <- old.call
      best.model$modtype            <- best.mod$modtype
      best.model$G                  <- best.mod$G
      best.model$bic                <- best.mod$bic
      best.model$icl                <- best.mod$icl
      best.model$aic                <- best.mod$aic
      best.model$dbs                <- if(attr(best.model, "DBS")) best.model$DBS[which(best.mod$G == as.numeric(rownames(best.model$DBS))),best.mod$modtype]
      best.model$asw                <- if(attr(best.model, "ASW")) best.model$ASW[which(best.mod$G == as.numeric(rownames(best.model$ASW))),best.mod$modtype]
      best.model$cv                 <- if(attr(best.model, "CV"))  best.model$CV[best.mod$G,best.mod$modtype]
      best.model$nec                <- if(attr(best.model, "NEC")) best.model$NEC[which(best.mod$G == as.numeric(rownames(best.model$NEC))),best.mod$modtype]
      best.model$gating             <- best.mod$gating
      best.model$loglik             <- best.mod$loglik
      best.model$df                 <- best.mod$df
      best.model$iters              <- best.mod$iters
      best.model$params             <- best.mod$params
      best.model$z                  <- best.mod$z
      best.model$MAP                <- best.mod$MAP
      best.model$uncert             <- best.mod$uncert
      attributes(best.model)        <- attributes(best.mod)
      attr(best.model, "Criterion") <- criterion
    } else best.model               <- paste0("Failed to re-fit the optimal model: ", gsub("\"", "'", deparse(old.call, width.cutoff=500L), fixed=TRUE))
  }
  gating[gating == "~1" | G   == 1] <- "None"
  noise         <- modelNames %in% c("CCN", "UCN", "CUN", "UUN")
  noise.gate    <- ifelse(!noise, NA, noise.gate[crit.names])
  equalPro      <- replace(unname(equalPro[crit.names]), gating != "None" | G == 1, NA)
  equalNoise    <- ifelse(!noise | G == 1, NA, equalNoise[crit.names] & vapply(equalPro, isTRUE, logical(1L)))
  comp          <- list(data = dat.name, optimal = best.model, pick = pick, MEDNames = crit.names, modelNames = modelNames, G = as.integer(G), 
                        df = as.integer(unname(dfxs[max.names])), iters = as.integer(unname(itxs[max.names])), 
                        bic = unname(bics[max.names]), icl = unname(icls[max.names]), aic = unname(aics[max.names]), 
                        dbs = replace(unname(dbss[max.names]), G == 1, NA), asw = replace(unname(asws[max.names]), G == 1, NA), 
                        cv = unname(cvs[max.names]), nec = replace(unname(necs[max.names]), G == 1, NA), loglik = unname(llxs[max.names]), 
                        gating = gating, algo = unname(algo[crit.names]), opti = unname(opti[crit.names]), weights = unname(weights[crit.names]),
                        equalPro = equalPro, noise = unname(noise), noise.gate = unname(replace(noise.gate, gating == "None" | G <= 2L - noise, NA)), 
                        equalNoise = unname(replace(equalNoise, !equalPro | is.na(equalPro), NA)))
  class(comp)   <- c("MEDseqCompare", "MEDseq")
  bic.tmp       <- sapply(BICs, as.vector)
  attr(comp, "Crit")   <- criterion
  attr(comp, "Opt")    <- optimal.only
  attr(comp, "NMods")  <- c(tried = sum(vapply(bic.tmp, function(x) length(x[!is.na(x)]),    numeric(1L))),
                            ran   = sum(vapply(bic.tmp, function(x) length(x[is.finite(x)]), numeric(1L))))
    comp
}

#' Set control values for use with MEDseq_fit
#'
#' Supplies a list of arguments (with defaults) for use with \code{\link{MEDseq_fit}}.
#' @param algo Switch controlling whether models are fit using the \code{"EM"} (the default) or \code{"CEM"} algorithm. The option \code{"cemEM"} allows running the EM algorithm starting from convergence of the CEM algorithm.
#' @param init.z The method used to initialise the cluster labels. All options respect the presence of sampling \code{weights}, if any. Defaults to \code{"kmedoids"}. Other options include \code{"kmodes"}, \code{"kmodes2"}, Ward's hierarchical clustering (\code{"hc"}, via \code{\link[stats]{hclust}}), \code{"random"} initialisation, and a user-supplied \code{"list"} (see \code{z.list} below). For weighted sequences, \code{"kmedoids"} is itself initialised using Ward's hierarchical clustering.
#' 
#' The \code{"kmodes"} and \code{"kmodes2"} options both internally call the function \code{\link{wKModes}}, which \emph{typically} uses random initial modes. Under \code{"kmodes"}, the algorithm is instead initialised via the medoids of the clusters obtained from a call to \code{\link[stats]{hclust}}. The option \code{"kmodes2"} is slightly faster, by virtue of using the \emph{random} initial medoids. However, final results are by default still subject to randomness under both options (unless \code{\link{set.seed}} is invoked), as ties for modes and cluster assignments are \emph{typically} broken at random throughout the algorithm (see the \code{random} argument below, and in \code{\link{wKModes}} itself).
#' @param z.list A user supplied list of initial cluster allocation matrices, with number of rows given by the number of observations, and numbers of columns given by the range of component numbers being considered. Only relevant if \code{init.z == "z.list"}. These matrices are allowed correspond to both soft or hard clusterings, and will be internally normalised so that the rows sum to 1.
#' @param dist.mat An optional distance matrix to use for initialisation when \code{init.z} is one of \code{"kmedoids"} or \code{"hc"}. Defaults to a Hamming distance matrix. This is an experimental feature and should only be tampered with by expert users.
#' @param unique A logical indicating whether the model is fit only to the unique observations (defaults to \code{TRUE}). When there are covariates, this means all unique combinations of covariate and sequence patterns, otherwise only the sequence patterns. 
#' 
#' When \code{weights} \emph{are not} supplied to \code{\link{MEDseq_fit}} and \code{isTRUE(unique)}, weights are given by the occurrence frequency of the corresponding sequences, and the model is then fit to the unique observations only.
#' 
#' When \code{weights} \emph{are} supplied and \code{isTRUE(unique)}, the weights are summed for each set of duplicate observations and assigned to one retained copy of each corresponding unique sequence. Hence, observations with different weights that are otherwise duplicates are treated as duplicates and significant computational gains can be made. 
#' 
#' In both cases, the results will be unchanged, but setting \code{unique} to \code{TRUE} can often be much faster. 
#' @param criterion When either \code{G} or \code{modtype} is a vector, \code{criterion} governs how the 'best' model is determined when gathering output. Defaults to \code{"bic"}. Note that all criteria will be returned in any case, if possible.
#' @param tau0 Prior mixing proportion for the noise component. If supplied, a noise component will be added to the model in the estimation, with \code{tau0} giving the prior probability of belonging to the noise component for \emph{all} observations. Typically supplied as a scalar in the interval (0, 1), e.g. \code{0.1}. Can be supplied as a vector when gating covariates are present and \code{noise.gate} is \code{TRUE}.
#' @param noise.gate A logical indicating whether gating network covariates influence the mixing proportion for the noise component, if any. Defaults to \code{TRUE}, but leads to greater parsimony if \code{FALSE}. Only relevant in the presence of a noise component (i.e. the \code{"CCN"}, \code{"UCN"}, \code{"CUN"}, and \code{"UUN"} models); only affects estimation in the presence of gating covariates.
#' @param random A logical governing how ties for estimated central sequence positions are handled. When \code{TRUE} (the default), such ties are broken at random. When \code{FALSE} (the implied default prior to version \code{1.2.0} of this package), the first candidate state is always chosen. This argument affects all \code{opti} options. If \code{verbose} is \code{TRUE} and there are tie-breaking operations performed, a warning message is printed once per model, regardless of the number of such operations. 
#' 
#' Note that this argument is \emph{also} passed to \code{\link{wKModes}} if \code{init.z} is \code{"kmodes"} or \code{"kmodes2"} and that, in certain rare cases when the \code{"CEM"} \code{algo} is invoked when \code{equalPro} is \code{TRUE} and the precision parameter(s) are somehow constrained across clusters, this argument also governs ties for cluster assignments within \code{MEDseq_fit} as well.
#' @param do.cv A logical indicating whether cross-validated log-likelihood scores should also be computed (see \code{nfolds}). Defaults to \code{FALSE} due to significant computational burden incurred.
#' @param do.nec A logical indicating whether the normalised entropy criterion (NEC) should also be computed (for models with more than one component). Defaults to \code{FALSE}. When \code{TRUE}, models with \code{G=1} are fitted always.
#' @param nfolds The number of folds to use when \code{isTRUE{do.cv}}.
#' @param nstarts The number of random initialisations to use when \code{init.z="random"}. Defaults to \code{1}. Results will be based on the random start yielding the highest estimated log-likelihood.
#' @param stopping The criterion used to assess convergence of the EM/CEM algorithm. The default (\code{"aitken"}) uses Aitken's acceleration method, otherwise the \code{"relative"} change in log-likelihood is monitored (which may be less strict).
#' @param equalPro Logical variable indicating whether or not the mixing proportions are to be constrained to be equal in the model. Default: \code{equalPro = FALSE}. Only relevant when \code{gating} covariates are \emph{not} supplied within \code{\link{MEDseq_fit}}, otherwise ignored. In the presence of a noise component, only the mixing proportions for the non-noise components are constrained to be equal (by default, see \code{equalNoise}), after accounting for the noise component.
#' @param equalNoise Logical which is \strong{only} invoked when \code{isTRUE(equalPro)} and gating covariates are not supplied. Under the default setting (\code{FALSE}), the mixing proportion for the noise component is estimated, and remaining mixing proportions are equal; when \code{TRUE} all components, including the noise component, have equal mixing proportions.
#' @param tol A vector of length two giving \emph{relative} convergence tolerances for 1) the log-likelihood of the EM/CEM algorithm, and 2) optimisation in the multinomial logistic regression in the gating network, respectively. The default is \code{c(1e-05, 1e-08)}. If only one number is supplied, it is used as the tolerance in both cases.
#' @param itmax A vector of length two giving integer limits on the number of iterations for 1) the EM/CEM algorithm, and 2) the multinomial logistic regression in the gating network, respectively. The default is \code{c(.Machine$integer.max, 1000)}. This allows termination of the EM/CEM algorithm to be completely governed by \code{tol[1]}. If only one number is supplied, it is used as the iteration limit for the EM/CEM algorithm only and the other element of \code{itmax} retains its usual default.
#' 
#' If, for any model with gating covariates, the multinomial logistic regression in the gating network fails to converge in \code{itmax[2]} iterations at any stage of the EM/CEM algorithm, an appropriate warning will be printed, prompting the user to modify this argument.
#' @param opti Character string indicating how central sequence parameters should be estimated. The default \code{"mode"} is exact and thus this experimental argument should only be tampered with by expert users. The option \code{"medoid"} fixes the central sequence(s) to be one of the observed sequences (like k-medoids). The other options \code{"first"} and \code{"GA"} use the first-improvement and genetic algorithms, respectively, to mutate the medoid. Pre-computation of the Hamming distance matrix for the observed sequences speeds-up computation of all options other than \code{"mode"}.
#' @param ordering Experimental feature that should only be tampered with by experienced users. Allows sequences to be reordered on the basis of the column-wise entropy when \code{opti} is \code{"first"} or \code{"GA"}.
#' @param MaxNWts The maximum allowable number of weights in the call to \code{\link[nnet]{multinom}} for the multinomial logistic regression in the gating network. There is no intrinsic limit in the code, but increasing \code{MaxNWts} will probably allow fits that are very slow and time-consuming. It may be necessary to increase \code{MaxNWts} when categorical concomitant variables with many levels are included or the number of components is high.
#' @param verbose Logical indicating whether to print messages pertaining to progress to the screen during fitting. By default is \code{TRUE} if the session is interactive, and \code{FALSE} otherwise. If \code{FALSE}, warnings and error messages will still be printed to the screen, but everything else will be suppressed.
#' @param ... Catches unused arguments, and also allows the optional arguments \code{ztol} and \code{summ} to be passed to \code{\link{dbs}} (\code{ztol} and \code{summ}) as well as the ASW computation (\code{summ}), and the optional \code{\link{wKModes}} arguments \code{iter.max}, \code{freq.weighted}, and \code{fast} (provided \code{init.z} is one of \code{"kmodes"} or \code{"kmodes2"}). In such cases, the \code{wKModes} argument \code{random} is already controlled by \code{random} above here.
#'
#' @return A named list in which the names are the names of the arguments and the values are the values supplied to the arguments.
#' @details \code{\link{MEDseq_control}} is provided for assigning values and defaults within \code{\link{MEDseq_fit}}. While the \code{criterion} argument controls the choice of the optimal number of components and MEDseq model type (in terms of the constraints or lack thereof on the precision parameters), \code{\link{MEDseq_compare}} is provided for choosing between fits with different combinations of covariates or different initialisation settings.
#' @importFrom nnet "multinom"
#' @importFrom TraMineR "seqdef"
#' @importFrom WeightedCluster "wcKMedoids"
#' @keywords control
#' @author Keefe Murphy - <\email{keefe.murphy@@mu.ie}>
#' @seealso \code{\link{MEDseq_fit}}, \code{\link{dbs}}, \code{\link[WeightedCluster]{wcKMedoids}}, \code{\link[cluster]{pam}}, \code{\link{wKModes}}, \code{\link[stats]{hclust}}, \code{\link[TraMineR]{seqdist}}, \code{\link[nnet]{multinom}}, \code{\link{MEDseq_compare}}
#' @references Murphy, K., Murphy, T. B., Piccarreta, R., and Gormley, I. C. (2021). Clustering longitudinal life-course sequences using mixtures of exponential-distance models. \emph{Journal of the Royal Statistical Society: Series A (Statistics in Society)}, 184(4): 1414-1451. <\href{https://rss.onlinelibrary.wiley.com/doi/abs/10.1111/rssa.12712}{doi:10.1111/rssa.12712}>.
#' 
#' Menardi, G. (2011). Density-based Silhouette diagnostics for clustering methods. \emph{Statistics and Computing}, 21(3): 295-308.
#' @export
#' @usage 
#' MEDseq_control(algo = c("EM", "CEM", "cemEM"), 
#'                init.z = c("kmedoids", "kmodes", "kmodes2", "hc", "random", "list"), 
#'                z.list = NULL, 
#'                dist.mat = NULL, 
#'                unique = TRUE, 
#'                criterion = c("bic", "icl", "aic", "dbs", "asw", "cv", "nec"), 
#'                tau0 = NULL, 
#'                noise.gate = TRUE, 
#'                random = TRUE,
#'                do.cv = FALSE, 
#'                do.nec = FALSE, 
#'                nfolds = 10L, 
#'                nstarts = 1L, 
#'                stopping = c("aitken", "relative"), 
#'                equalPro = FALSE, 
#'                equalNoise = FALSE, 
#'                tol = c(1E-05, 1E-08), 
#'                itmax = c(.Machine$integer.max, 1000L), 
#'                opti = c("mode", "medoid", "first", "GA"), 
#'                ordering = c("none", "decreasing", "increasing"), 
#'                MaxNWts = 1000L, 
#'                verbose = TRUE, 
#'                ...)
#' @examples
#' # The CC MEDseq model is almost equivalent to k-medoids when the
#' # CEM algorithm is employed, mixing proportions are constrained,
#' # and the central sequences are restricted to the observed sequences
#' ctrl   <- MEDseq_control(algo="CEM", equalPro=TRUE, opti="medoid", criterion="asw")
#' \donttest{
#' data(mvad)
#' # Note that ctrl must be explicitly named 'ctrl'
#' mod   <- MEDseq_fit(seqdef(mvad[,17:86]), G=11, modtype="CC", weights=mvad$weight, ctrl=ctrl)
#' 
#' # Alternatively, specify the control arguments directly
#' mod   <- MEDseq_fit(seqdef(mvad[,17:86]), G=11, modtype="CC", weights=mvad$weight,
#'                     algo="CEM", equalPro=TRUE, opti="medoid", criterion="asw")
#' 
#' # Note that supplying control arguments via a mix of the ... construct and the named argument 
#' # 'control' or supplying MEDseq_control output without naming it 'control' can throw an error}
MEDseq_control    <- function(algo = c("EM", "CEM", "cemEM"), init.z = c("kmedoids", "kmodes", "kmodes2", "hc", "random", "list"), z.list = NULL, dist.mat = NULL, unique = TRUE, 
                              criterion = c("bic", "icl", "aic", "dbs", "asw", "cv", "nec"), tau0 = NULL, noise.gate = TRUE, random = TRUE, do.cv = FALSE, do.nec = FALSE, nfolds = 10L, 
                              nstarts = 1L, stopping = c("aitken", "relative"), equalPro = FALSE, equalNoise = FALSE, tol = c(1E-05, 1E-08), itmax = c(.Machine$integer.max, 1000L), 
                              opti = c("mode", "medoid", "first", "GA"), ordering = c("none", "decreasing", "increasing"), MaxNWts = 1000L, verbose = TRUE, ...) {
  miss.args                <- list(tau0=missing(tau0), init.z = missing(init.z), z.list = missing(z.list))
  if(!missing(algo)        &&
    (length(algo)      > 1 ||
     !is.character(algo)))       stop("'algo' must be a character vector of length 1",      call.=FALSE)
  if(!missing(init.z)      &&
    (length(init.z)    > 1 ||
     !is.character(init.z)))     stop("'init.z' must be a character vector of length 1",    call.=FALSE)
  init.z                   <- match.arg(init.z)
  if(!missing(dist.mat))    {
    dist.mat               <- tryCatch(suppressWarnings(stats::as.dist(dist.mat)), error=function(e)     {
                                 stop("'dist.mat' must be coercible to the class 'dist'",   call.=FALSE) })
  }
  if(length(unique)    > 1 ||
     !is.logical(unique))        stop("'unique' must be a single logical indicator",        call.=FALSE)
  if(init.z == "random")    {
   if(length(nstarts) != 1 ||
      !is.numeric(nstarts) ||
      (nstarts         < 1 ||
       floor(nstarts) !=
       nstarts))                 stop(paste0("'nstarts' must be a single integer >= 1 if when 'init.z'=", init.z), call.=FALSE)
  }
  if(!missing(criterion)   &&
    (length(criterion) > 1 ||
     !is.character(criterion)))  stop("'criterion' must be a character vector of length 1", call.=FALSE)
  if(length(random)    > 1 ||
     !is.logical(random))        stop("'random' must be a single logical indicator",        call.=FALSE)
  if(length(do.cv)     > 1 ||
     !is.logical(do.cv))         stop("'do.cv' must be a single logical indicator",         call.=FALSE)
  if(length(do.nec)    > 1 ||
     !is.logical(do.nec))        stop("'do.nec' must be a single logical indicator",        call.=FALSE)
  if(!missing(stopping)    &&
    (length(stopping)  > 1 ||
     !is.character(stopping)))   stop("'stopping' must be a character vector of length 1",  call.=FALSE)
  if(!miss.args$tau0       &&
    (!is.numeric(tau0)     ||
     any(tau0  < 0)        || 
     any(tau0 >= 1)))            stop("'tau0' must lie in the interval [0, 1)",             call.=FALSE)
  if(!missing(opti)        &&
    (length(opti)      > 1 ||
     !is.character(opti)))       stop("'opti' must be a character vector of length 1",      call.=FALSE)
  if(!missing(ordering)    &&
    (length(ordering)  > 1 ||
     !is.character(ordering)))   stop("'ordering' must be a character vector of length 1",  call.=FALSE)
  if(length(noise.gate)     > 1 ||
     !is.logical(noise.gate))    stop("'noise.gate' must be a single logical indicator",    call.=FALSE)
  if(length(MaxNWts)   > 1 ||
     !is.numeric(MaxNWts)  ||
     MaxNWts      <= 0)          stop("'MaxNWts' must a strictly positive scalar",          call.=FALSE)
  if(length(equalPro)  > 1 ||
     !is.logical(equalPro))      stop("'equalPro' must be a single logical indicator",      call.=FALSE)
  if(length(equalNoise)     > 1 ||
     !is.logical(equalNoise))    stop("'equalNoise' must be a single logical indicator",    call.=FALSE)
  if(equalNoise   && !equalPro  &&
     isTRUE(verbose))            message("'equalNoise' forced to FALSE as 'equalPro' is FALSE\n")
  equalNoise      <- equalPro   && equalNoise
  if((len.tol     <- 
      length(tol)) > 2     ||
     !is.numeric(tol))           stop("'tol' must be a numeric vector of length at most 2", call.=FALSE)
  if(any(tol   < 0,
         tol  >= 1))             stop("'tol' must be in the interval [0, 1)",               call.=FALSE)
  if(len.tol  == 1)    tol <- rep(tol, 2L)
  if(length(itmax)         == 1) {
    itmax     <- c(itmax, 1000L)
  } else if(length(itmax)  != 2) stop("'itmax' must be of length 2",                        call.=FALSE)
  if(!is.numeric(itmax)    ||
     any(floor(itmax) != itmax) ||
     any(itmax    <= 0))         stop("'itmax' must contain strictly positive integers",    call.=FALSE)
  inf         <- is.infinite(itmax)
  if(any(inf))   itmax[inf]     <- .Machine$integer.max
  itmax[1L]   <- ifelse(itmax[1L] == .Machine$integer.max, itmax[1L], itmax[1L] + 2L)
  if(length(verbose)   > 1 ||
     !is.logical(verbose))       stop("'verbose' must be a single logical indicator",       call.=FALSE)
  pamonce     <- ifelse(!is.null(list(...)$pamonce), list(...)$pamonce, ifelse(.version_above("cluster", "2.0.8"), 5, ifelse(.version_above("cluster", "1.14.2"), 2, 0)))
  control                  <- list(algo = match.arg(algo), init.z = init.z, dist.mat = dist.mat, nstarts = nstarts, criterion = match.arg(criterion), nfolds = nfolds, random = random, do.cv = do.cv, 
                                   do.nec = do.nec, MaxNWts = MaxNWts, stopping = match.arg(stopping), tau0 = tau0, opti = match.arg(opti), ordering = match.arg(ordering), noise.gate = noise.gate, unique = unique,
                                   equalPro = equalPro, equalNoise = equalNoise, tol = tol[1L], g.tol = tol[2L], itmax = itmax[1L], g.itmax = itmax[2L], verbose = verbose, z.list = z.list, pamonce = pamonce)
  attr(control, "missing") <- miss.args
    return(control)
}

#' MEDseq: Mixtures of Exponential-Distance Models with Covariates
#'
#' Fits MEDseq models: mixtures of Exponential-Distance models with gating covariates and sampling weights. Typically used for clustering categorical/longitudinal life-course sequences. Additional arguments are available via the function \code{\link{MEDseq_control}}.
#' @param seqs A state-sequence object of class \code{"stslist"} as created by the \code{\link[TraMineR]{seqdef}} function in the \pkg{TraMineR} package. Note that the data set must have equal sequence lengths, the intervals are assumed to be evenly spaced, and missingness is not allowed.
#' @param G A positive integer vector specifying the numbers of mixture components (clusters) to fit. Defaults to \code{G=1:9}.
#' @param modtype A vector of character strings indicating the type of MEDseq models to be fitted, in terms of the constraints or lack thereof on the precision parameters. By default, all valid model types are fitted (except some only where \code{G > 1} or \code{G > 2}, see \code{Note}). 
#' The models are named \code{"CC"}, \code{"CU"}, \code{"UC"}, \code{"UU"}, \code{"CCN"}, \code{"CUN"}, \code{"UCN"}, and \code{"UUN"}. The first letter denotes whether the precision parameters are constrained/unconstrained across clusters. The second letter denotes whether the precision parameters are constrained/unconstrained across sequence positions (i.e. time points). The third letter denotes whether one of the components is constrained to have zero-precision/infinite variance. Such a noise component assumes sequences in that cluster follow a uniform distribution.
#' @param gating A \code{\link[stats]{formula}} for determining the model matrix for the multinomial logistic regression in the gating network when fixed covariates enter the mixing proportions. Defaults to \code{~1}, i.e. no covariates. This will be ignored where \code{G=1}. Continuous, categorical, and/or ordinal covariates are allowed. Logical covariates will be coerced to factors. Interactions, transformations, and higher order terms are permitted: the latter \strong{must} be specified explicitly using the \code{AsIs} operator (\code{\link{I}}). The specification of the LHS of the formula is ignored. Intercept terms are included by default.
#' @param covars An optional data frame (or a matrix with named columns) in which to look for the covariates in the \code{gating} network formula, if any. If not found in \code{covars}, any supplied \code{gating} covariates are taken from the environment from which \code{MEDseq_fit} is called. Try to ensure the names of variables in \code{covars} do not match any of those in \code{seqs}.
#' @param weights Optional numeric vector containing observation-specific sampling weights, which are accounted for in the model fitting and other functions where applicable. \code{weights} are always internally normalised to sum to the sample size. See the \code{unique} argument to \code{\link{MEDseq_control}} to see how incorporating weights also yields computational benefits. Note that \code{weights} must \strong{always} be explicitly supplied here; it is not enough to use weights when constructing the state sequence object via \code{\link[TraMineR]{seqdef}}. If you \emph{are} using a weighted \code{"stslist"} state sequence object and do not specify \code{weights}, you will be prompted to explicitly specify \code{weights=attr(seqs, "weights")} for a weighted model or \code{weights=NULL} for an unweighted model.
#' @param ctrl A list of control parameters for the EM/CEM and other aspects of the algorithm. The defaults are set by a call to \code{\link{MEDseq_control}}.
#' @param ... Catches unused arguments (see \code{\link{MEDseq_control}}).
#' @param x,object,digits,classification,parameters,network,SPS Arguments required for the \code{print} and \code{summary} functions: \code{x} and \code{object} are objects of class \code{"MEDseq"} resulting from a call to \code{\link{MEDseq_fit}}, while \code{digits} gives the number of decimal places to round to for printing purposes (defaults to \code{3}). \code{classification}, \code{parameters}, and \code{network} are logicals which govern whether a table of the MAP classification of observations, the mixture component parameters, and the gating network coefficients are printed, respectively. \code{SPS} governs the printing of the relevant quantities in \code{"summaryMEDseq"} objects when any of \code{classification}, \code{parameters}, &/or \code{network} are \code{TRUE} (see \code{\link{MEDseq_clustnames}} and \code{\link[TraMineR]{seqformat}}).
#'
#' @return A list (of class \code{"MEDseq"}) with the following named entries (of which some may be missing, depending on the \code{criterion} employed), mostly corresponding to the chosen optimal model (as determined by the \code{criterion} within \code{\link{MEDseq_control}}):
#' \item{\code{call}}{The matched call.}
#' \item{\code{data}}{The input data, \code{seqs}.}
#' \item{\code{modtype}}{A character string denoting the MEDseq model type at which the optimal \code{criterion} occurs.}
#' \item{\code{G}}{The optimal number of mixture components according to \code{criterion}.}
#' \item{\code{params}}{A list with the following named components:
#' \describe{
#' \item{\code{theta}}{A matrix with \code{G} rows and T columns, where T is the number of sequence positions, giving the central sequences of each cluster. The mean of the noise component is not reported, as it does not contribute in any way to the likelihood. A dedicated \code{print} function is provided.}
#' \item{\code{lambda}}{A matrix of precision parameters. Will contain \code{1} row if the 1st letter of \code{modtype} is "C" and \code{G} columns otherwise. Will contain \code{1} column if the 2nd letter of \code{modtype} is "C" and T columns otherwise, where T is the number of sequence positions. Precision parameter values of zero are reported for the noise component, if any. Note that values of \code{Inf} are also possible, corresponding to zero-variance, which is most likely under the \code{"UU"} or \code{"UUN"} models. A dedicated \code{print} function is provided.}
#' \item{\code{tau}}{The mixing proportions: either a vector of length \code{G} or, if \code{gating} covariates were supplied, a matrix with an entry for each observation (rows) and component (columns).}}
#' }
#' \item{\code{gating}}{An object of class \code{"MEDgating"} (for which dedicated \code{print}, \code{summary}, and \code{\link[=predict.MEDgating]{predict}} methods exist) and either \code{"multinom"} or \code{"glm"} (only for single-component models) giving the \code{\link[nnet]{multinom}} regression coefficients of the \code{gating} network. If \code{gating} covariates were \emph{NOT} supplied (or the best model has just one component), this corresponds to a RHS of \code{~1}, otherwise the supplied \code{gating} formula. As such, a fitted \code{gating} network is always returned even in the absence of supplied covariates or clusters. If there is a noise component (and the option \code{noise.gate=TRUE} is invoked), its coefficients are those for the \emph{last} component. \strong{Users are cautioned against making inferences about statistical significance from summaries of the coefficients in the gating network. Users are instead advised to use the function \code{\link{MEDseq_stderr}}}.}
#' \item{\code{z}}{The final responsibility matrix whose \code{[i,k]}-th entry is the probability that observation \emph{i} belongs to the \emph{k}-th component. If there is a noise component, its values are found in the \emph{last} column.}
#' \item{\code{MAP}}{The vector of cluster labels for the chosen model corresponding to \code{z}, i.e. \code{max.col(z)}. Observations belonging to the noise component, if any, will belong to component \code{0}.}
#' \item{\code{BIC}}{A matrix of \emph{all} BIC values with \code{length{G}} rows and \code{length(modtype)} columns. See \code{Note}.}
#' \item{\code{ICL}}{A matrix of \emph{all} ICL values with \code{length{G}} rows and \code{length(modtype)} columns. See \code{Note}.}
#' \item{\code{AIC}}{A matrix of \emph{all} AIC values with \code{length{G}} rows and \code{length(modtype)} columns. See \code{Note}.}
#' \item{\code{DBS}}{A matrix of \emph{all} (weighted) mean/median DBS values with \code{length{G}} rows and \code{length(modtype)} columns. See \code{Note} and \code{\link{dbs}}.}
#' \item{\code{DBSvals}}{A list of lists giving the observation-specific DBS values for \emph{all} fitted models. The first level of the list corresponds to numbers of components, the second to the MEDseq model types.}
#' \item{\code{dbs}}{The (weighted) mean/median DBS value corresponding to the optimal model. May not necessarily be the optimal DBS.}
#' \item{\code{dbsvals}}{Observation-specific DBS values corresponding to the optimum model, which may not be optimal in terms of DBS.}
#' \item{\code{ASW}}{A matrix of \emph{all} (weighted) mean/median ASW values with \code{length{G}} rows and \code{length(modtype)} columns. See \code{Note}.}
#' \item{\code{ASWvals}}{A list of lists giving the observation-specific ASW values for \emph{all} fitted models. The first level of the list corresponds to numbers of components, the second to the MEDseq model types.}
#' \item{\code{asw}}{The (weighted) mean/median ASW value corresponding to the optimal model. May not necessarily be the optimal ASW.}
#' \item{\code{aswvals}}{Observation-specific ASW values corresponding to the optimum model, which may not be optimal in terms of ASW.}
#' \item{\code{LOGLIK}}{A matrix of \emph{all} maximal log-likelihood values with \code{length{G}} rows and \code{length(modtype)} columns. See \code{Note}.}
#' \item{\code{DF}}{A matrix giving the numbers of estimated parameters (i.e. the number of 'used' degrees of freedom) for \emph{all} visited models, with \code{length{G}} rows and \code{length(modtype)} columns. Subtract these numbers from the sample size to get the degrees of freedom. See \code{Note}.}
#' \item{\code{ITERS}}{A matrix giving the total number of EM/CEM iterations for \emph{all} visited models, with \code{length{G}} rows and \code{length(modtype)} columns. See \code{Note}.}
#' \item{\code{CV}}{A matrix of \emph{all} cross-validated log-likelihood values with \code{length{G}} rows and \code{length(modtype)} columns, if available. See \code{Note} and the arguments \code{do.cv} and \code{nfolds} to \code{\link{MEDseq_control}}.}
#' \item{\code{NEC}}{A matrix of \emph{all} NEC values with \code{length{G}} rows and \code{length(modtype)} columns, if available. See \code{Note} and the argument \code{do.nec} to \code{\link{MEDseq_control}}.}
#' \item{\code{bic}}{The BIC value corresponding to the optimal model. May not necessarily be the optimal BIC.}
#' \item{\code{icl}}{The ICL value corresponding to the optimal model. May not necessarily be the optimal ICL.}
#' \item{\code{aic}}{The AIC value corresponding to the optimal model. May not necessarily be the optimal AIC.}
#' \item{\code{loglik}}{The vector of increasing log-likelihood values for every EM/CEM iteration under the optimal model. The last element of this vector is the maximum log-likelihood achieved by the parameters returned at convergence.}
#' \item{\code{df}}{The number of estimated parameters in the optimal model (i.e. the number of 'used' degrees of freedom). Subtract this number from the sample size to get the degrees of freedom.}
#' \item{\code{iters}}{The total number of EM/CEM iterations for the optimal model.}
#' \item{\code{cv}}{The cross-validated log-likelihood value corresponding to the optimal model, if available. May not necessarily be the optimal one.}
#' \item{\code{nec}}{The NEC value corresponding to the optimal model, if available. May not necessarily be the optimal NEC.}
#' \item{\code{ZS}}{A list of lists giving the \code{z} matrices for \emph{all} fitted models. The first level of the list corresponds to numbers of components, the second to the MEDseq model types.}
#' \item{\code{uncert}}{The uncertainty associated with the \code{classification}.}
#' \item{\code{covars}}{A data frame gathering the set of covariates used in the \code{gating} network, if any. Will contain zero columns in the absence of gating covariates. Supplied gating covariates will be excluded if the optimal model has only one component. May have fewer columns than covariates supplied via the \code{covars} argument also, as only the included covariates are gathered here.}
#' Dedicated \code{\link[=plot.MEDseq]{plot}}, \code{print}, and \code{summary} functions exist for objects of class \code{"MEDseq"}. 
#' @details The function effectively allows 8 different MEDseq precision parameter settings for models with or without gating network covariates. By constraining the mixing proportions to be equal (see \code{equalPro} in \code{\link{MEDseq_control}}) an extra special case is facilitated in the latter case. 
#' 
#' While model selection in terms of choosing the optimal number of components and the MEDseq model type is performed within \code{\link{MEDseq_fit}}, using one of the \code{criterion} options within \code{\link{MEDseq_control}}, choosing between multiple fits with different combinations of covariates or different initialisation settings can be done by supplying objects of class \code{"MEDseq"} to \code{\link{MEDseq_compare}}.
#' @note Where \code{BIC}, \code{ICL}, \code{AIC}, \code{DBS}, \code{ASW}, \code{LOGLIK}, \code{DF}, \code{ITERS}, \code{CV}, and \code{NEC} contain \code{NA} entries, this corresponds to a model which was not run; for instance a UU model is never run for single-component models as it is equivalent to CU, while a UCN model is never run for two-component models as it is equivalent to CCN. As such, one can consider the value as not really missing, but equivalent to the corresponding value. On the other hand, \code{-Inf} represents models which were terminated due to error, for which a log-likelihood could not be estimated. These objects all inherit the class \code{"MEDCriterion"} for which dedicated \code{print} and \code{summary} methods exist. For plotting, please see \code{\link[=plot.MEDseq]{plot}}.

#' @importFrom cluster "agnes" "pam"
#' @importFrom matrixStats "colSums2" "logSumExp" "rowLogSumExps" "rowMaxs" "rowMeans2" "rowSums2" "weightedMedian" "weightedMean"
#' @importFrom nnet "multinom"
#' @importFrom stringdist "stringdistmatrix"
#' @importFrom TraMineR "disscenter" "seqdef" "seqformat"
#' @importFrom WeightedCluster "wcKMedoids" "wcSilhouetteObs"
#' @export
#' @author Keefe Murphy - <\email{keefe.murphy@@mu.ie}>
#' @references Murphy, K., Murphy, T. B., Piccarreta, R., and Gormley, I. C. (2021). Clustering longitudinal life-course sequences using mixtures of exponential-distance models. \emph{Journal of the Royal Statistical Society: Series A (Statistics in Society)}, 184(4): 1414-1451. <\href{https://rss.onlinelibrary.wiley.com/doi/abs/10.1111/rssa.12712}{doi:10.1111/rssa.12712}>.
#' @keywords clustering main
#' @seealso \code{\link[TraMineR]{seqdef}}, \code{\link{MEDseq_control}}, \code{\link{MEDseq_compare}}, \code{\link{plot.MEDseq}}, \code{\link{predict.MEDgating}}, \code{\link{MEDseq_stderr}}, \code{\link{I}}, \code{\link{MEDseq_clustnames}}, \code{\link[TraMineR]{seqformat}}
#' @usage
#' MEDseq_fit(seqs, 
#'            G = 1L:9L, 
#'            modtype = c("CC", "UC", "CU", "UU", 
#'                        "CCN", "UCN", "CUN", "UUN"), 
#'            gating = NULL, 
#'            weights = NULL, 
#'            ctrl = MEDseq_control(...), 
#'            covars = NULL, 
#'            ...)
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
#' 
#' # Fit a range of exponential-distance models without clustering
#' mod0          <- MEDseq_fit(mvad.seq, G=1)
#' 
#' # Fit a range of unweighted mixture models without covariates
#' # Only consider models with a noise component
#' # Supply some MEDseq_control() arguments
#' \donttest{
#' # mod1        <- MEDseq_fit(mvad.seq, G=9:10, modtype=c("CCN", "CUN", "UCN", "UUN"),
#' #                           algo="CEM", init.z="kmodes", criterion="icl")
#' 
#' # Fit a model with weights and a gating covariate
#' # Have the probability of noise-component membership be constant
#' mod2          <- MEDseq_fit(mvad.seq, G=11, modtype="UUN", weights=mvad$weights, 
#'                             gating=~ gcse5eq, covars=mvad.cov, noise.gate=FALSE)
#'                             
#' # Examine this model in greater detail
#' summary(mod2, classification=TRUE, parameters=TRUE)
#' summary(mod2$gating, SPS=TRUE)
#' print(mod2$params$theta, SPS=TRUE)
#' plot(mod2, "clusters")}
MEDseq_fit        <- function(seqs, G = 1L:9L, modtype = c("CC", "UC", "CU", "UU", "CCN", "UCN", "CUN", "UUN"), 
                              gating = NULL, weights = NULL, ctrl = MEDseq_control(...), covars = NULL, ...) {
  call            <- cX  <- match.call()
  if(!inherits(seqs, "stslist")) stop("'seqs' must be of class 'stslist'",        call.=FALSE)
  if(any(seqs     ==
         attr(seqs, "nr")))      stop("Missing values in 'seqs' are not allowed", call.=FALSE)
  SEQ             <- apply(.fac_to_num(seqs), 1L, .num_to_char)
  seqX            <- seqs
  levs            <- attr(seqs, "alphabet")
  attr(SEQ, "N")  <- N2  <- attr(SEQ, "W") <- N <- nrow(seqs)
  attr(SEQ, "T")  <- P   <- ncol(seqs)
  attr(SEQ, "V")  <- V   <- length(levs)
  attr(SEQ, "V1") <- V1  <- V - 1L
  attr(SEQ, "V1V")       <- V1/V
  attr(SEQ, "VTS")       <- sum(lengths(apply(seqs, 2L, unique)) - 1L)
  attr(SEQ, "logV1")     <- log(V1)
  attr(SEQ, "lPV")       <- P * log(V)
  if(any(c(N, P, V)      <= 1))  stop("The number of sequences, the sequence length, and the size of the sequence alphabet must all be > 1", call.=FALSE)
  if(V            >= 10)         warning(paste0("Sequence alphabet size (V=", V, ") represents a large number of states/categories: are you sure?\n"), call.=FALSE, immediate.=TRUE)
  if(!is.null(modtype)   &&
     !is.character(modtype))     stop("'modtype' must be a character vector", call.=FALSE)
  modtype         <- if(is.null(modtype)) c("CC", "UC", "CU", "UU", "CCN", "UCN", "CUN", "UUN") else toupper(modtype)
  l.meth          <- match.arg(modtype, several.ok=TRUE)
  dots            <- list(...)
  if(!missing(ctrl)      &&
     length(dots[names(dots) %in%
     names(ctrl)]) > 0)          stop("Arguments cannot be supplied via the '...' construct when the named argument 'ctrl' is supplied", call.=FALSE) 
  algo            <- ctrl$algo
  opti            <- ctrl$opti
  criterion       <- ctrl$criterion
  init.z          <- ctrl$init.z
  z.list          <- ctrl$z.list
  dist.mat        <- ctrl$dist.mat
  nstarts         <- switch(EXPR=init.z, random=ctrl$nstarts, 1L)
  startseq        <- seq_len(nstarts)
  equalPro        <- ctrl$equalPro
  equalNoise      <- ctrl$equalNoise
  noise.gate      <- ctrl$noise.gate
  verbose         <- ctrl$verbose
  pamonce         <- ctrl$pamonce
  ctrl$warn       <- TRUE
  x.ctrl          <- list(equalPro=equalPro, noise.gate=noise.gate, equalNoise=equalNoise)
  ctrl$ordering   <- ifelse(opti == "first", ctrl$ordering, "none")
  miss.args       <- attr(ctrl, "missing")
  zin.miss        <- miss.args$init.z
  zli.miss        <- miss.args$z.list
  if(!zli.miss)    {
    if(!inherits(z.list,
       "list")    ||
       !all(vapply(z.list, inherits,
       logical(1L), "matrix")))  stop("'z.list' must be a list of matrices if supplied", call.=FALSE)
    if(zin.miss   &&
       init.z     != "list")   { 
      init.z      <- "list"
      if(isTRUE(verbose))        message("'init.z' set to 'list' as 'z.list' was supplied\n")
    }
  }
  covmiss         <- missing(covars)
  mt1             <- unique(vapply(l.meth, function(lx) switch(EXPR=lx, CC=, UC="CC", CU=, UU="CU", "CCN"), character(1L)))
  mt2             <- unique(vapply(l.meth, function(lx) switch(EXPR=lx, UCN="CCN", UUN="CUN", lx),          character(1L)))
  mtg             <- unique(l.meth)
  n.meths         <- c("CCN", "UCN", "CUN", "UUN")
  n.meth          <- l.meth %in% n.meths
  if(!(tmiss      <- miss.args$tau0) && 
    all(ctrl$tau0 == 0))       {
    if(all(n.meth))              stop("'tau0' is zero: models with noise component cannot be fitted",        call.=FALSE)
    if(any(n.meth))              warning("'tau0' is zero: models with noise component will not be fitted\n", call.=FALSE, immediate.=TRUE)
    l.meths       <- c("CC", "UC", "CU", "UU")
  } else l.meths  <- c("CC", "UC", "CU", "UU", "CCN", "UCN", "CUN", "UUN")
  HAM.mat         <- suppressMessages(seqdist(seqs, "HAM", full.matrix=FALSE))
  if(is.null(dist.mat))        { 
    dist.mat      <- HAM.mat
  } else if((sqrt(length(dist.mat) * 2L + N) !=
             N))                 stop("Invalid 'dist.mat' dimensions", call.=FALSE)
  HAM.mat         <- as.matrix(HAM.mat)
  
  if((gate.x      <- !missing(gating))) {
    if(!inherits(gating, 
                 "formula"))     stop("'gating' must be a formula", call.=FALSE)
    if(!covmiss)   {   
     if((!is.data.frame(covars) &&
         !is.matrix(covars))    ||
        (ncol(covars)  > 0      &&
     is.null(colnames(covars)))) stop("'covars' must be a data.frame or matrix with named columns if supplied", call.=FALSE)
      call$covars <- covars   <- as.data.frame(covars)
      gvars       <- setdiff(all.vars(gating), ".")
      if(length(gvars) > 0    && 
         !all(gvars %in% 
              names(covars)))    stop("One or more variables in gating formula not found in covars", call.=FALSE)
    }
    if(inherits(try(stats::terms(gating), silent=TRUE), "try-error")) {
      if(covmiss)                stop("Can't use '.' in 'gating' formula without supplying 'covars' argument", call.=FALSE)
      gating      <- attr(stats::terms(gating, data=covars), "term.labels")
      gating      <- stats::reformulate(if(length(gating) == 0) "1" else gating, response="z")
    }
    gating        <- tryCatch(stats::update.formula(stats::as.formula(gating), z ~ .), error=function(e) {
                                 stop("Invalid 'gating' network formula supplied", call.=FALSE) })
    if(gating[[3L]]      == 1) { 
      if(verbose)                message("Not including gating network covariates with only intercept on gating formula RHS\n")
      gate.x      <- FALSE
    }
    if(gating[[3L]]      == 
       "1 - 1")                  stop("'gating' formula must include an intercept when it doesn't include covariates", call.=FALSE)
    gate.names    <- stats::terms(gating)
    gate.names    <- labels(gate.names)[attr(gate.names, "order") <= 1]
  } 
  if(gate.x)       {
    covars        <- eval(bquote(stats::model.frame(.(stats::update.formula(gating, NULL ~ .)), data=.(call$covars), drop.unused.levels=TRUE)), envir=parent.frame(), enclos=environment())
    gate.names    <- colnames(covars)
    covars        <- cbind(covars, eval(bquote(stats::model.frame(.(stats::as.formula(paste("~", paste(eval(bquote(all.vars(.(gating))), envir=parent.frame())[-1L], collapse="+")))), data=.(call$covars), drop.unused.levels=TRUE)), envir=parent.frame(), enclos=environment()))
    covars        <- covars[,unique(colnames(covars)), drop=FALSE]
    covch         <- vapply(covars, is.character, logical(1L))
    covars[covch] <- lapply(covars[covch], factor)
    if(any(covch) && verbose)    message("Character covariates coerced to factors\n")
  }
  gate.names      <- if(gate.x)  gate.names[!is.na(gate.names)]
  if(!covmiss)     {
    if(!all(gate.names  %in% 
            colnames(covars)))   stop("Supplied gating covariates not found in supplied 'covars'", call.=FALSE)
    covars        <- if(gate.x)  covars                                                               else as.data.frame(matrix(0L, nrow=N, ncol=0L))
  } else {
    if(any(grepl("\\$", 
                 gate.names)))   stop("Don't supply covariates to the gating network using the $ operator: use the 'covars' argument instead", call.=FALSE)
    covars        <- if(gate.x)  data.frame(stats::model.frame(gating[-2L], drop.unused.levels=TRUE)) else as.data.frame(matrix(0L, nrow=N, ncol=0L))
  }
  if(nrow(covars) != N)          stop("'gating' covariates must contain the same number of rows as 'seqs'", call.=FALSE)
  glogi           <- vapply(covars, is.logical, logical(1L))
  covars[,glogi]  <- sapply(covars[,glogi], as.factor)
  if(covmiss)      {
    covars        <- data.frame(covars, stringsAsFactors=TRUE)
  }
  
  if(ctrl$do.wts  <- do.wts   <- 
     !is.null(weights))   {
    if(!is.numeric(weights)   ||
       length(weights)   != N)   stop(paste0("'weights' must be a numeric vector of length N=", N), call.=FALSE)
    if(any(weights < 0)  || 
       any(!is.finite(weights))) stop("'weights' must be non-negative and finite", call.=FALSE)
    if(!is.null(attr(seqs, "weights"))        &&
       !identical(unname(attr(seqs, "weights")), 
        weights))                warning(paste0("Supplied 'weights' differ from the weights attribute of the sequence object: attr(", deparse(substitute(seqs)), ", \"weights\")\n"), call.=FALSE, immediate.=TRUE)
    if(ctrl$do.wts       <- 
       do.wts     <- (length(unique(weights))  > 1))  {
      weights     <- replace(weights, weights == 0, .Machine$double.eps)
      weights     <- WEIGHTS  <- 
      attr(SEQ, "Weights")    <- weights/sum(weights) * N
    }
  } else           {
    if(missing(weights)       &&
       !is.null(attr(seqs, 
                "weights")))     stop(paste0("Weighted sequences supplied, but you must explicitly specify the 'weights' argument:\n\tIf you do desire a weighted model, use 'weights=attr(", deparse(substitute(seqs)), ", \"weights\")'\n\tIf you desire an unweighted model, use 'weights=NULL'"), call.=FALSE)
    WEIGHTS       <- rep(1L, N2)
  } 
  do.uni          <- ctrl$unique
  DF              <- seqs   
  uni.sum         <- sum(!duplicated(DF))
  DF              <- if(gate.x)  cbind(DF, covars) else DF
  uni.ind         <- !duplicated(DF)
  sum.uni         <- sum(uni.ind)
  if(sum.uni < N  && !do.uni  &&
     verbose)                   message(paste0("Number of unique observations (", sum.uni, ") is less than N (", N, ")", ifelse(uni.sum < sum.uni, paste0(" - \nNumber of unique sequences only, ignoring covariates, is ", uni.sum, ":\n"), ": "), "Consider setting 'unique'=TRUE\n"))
  if(do.uni       <- do.uni   &&
     sum.uni < N)  {            
    if(verbose)                 message(paste0("Proceeding with ", sum.uni, " unique observation", ifelse(sum.uni == 1, "", "s"), " out of N=", N, "\n"))
    if(do.wts)     {
      agg.DF      <- merge(data.frame(cbind(id=seq_len(N)), DF, check.names=FALSE, wts=weights), data.frame(stats::aggregate(cbind(DF[0L], count=1L), DF, length), check.names=FALSE), by=colnames(DF), sort=FALSE)
      agg.id      <- order(agg.DF$id)
      c2          <- agg.DF$count[agg.id]
      counts      <- c2[uni.ind]
      grp         <- apply(agg.DF[,seq_len(ncol(DF))], 1L, paste, collapse="")
      weights     <- stats::aggregate(agg.DF$wts, by=list(grp=grp), sum)
      weights     <- weights[match(grp, weights$grp),]$x[agg.id[uni.ind]]
    } else         {
      agg.DF      <- merge(data.frame(cbind(id=seq_len(N)), DF, check.names=FALSE),              data.frame(stats::aggregate(cbind(DF[0L], count=1L), DF, length), check.names=FALSE), by=colnames(DF), sort=FALSE)
      agg.id      <- order(agg.DF$id)
      c2          <- agg.DF$count[agg.id]
      weights     <- 
      counts      <- c2[uni.ind]
    }
    if(ctrl$do.wts       <- (length(unique(weights)) < N)) {
      dist.mat2   <- dist.mat
      HAM.mat2    <- HAM.mat[uni.ind,uni.ind]
      w2          <- integer(N)
      w2[uni.ind] <- weights
      dist.mat    <- stats::as.dist(as.matrix(dist.mat)[uni.ind,uni.ind])
      seqs        <- seqs[uni.ind,,   drop=FALSE]
      atts        <- attributes(SEQ)
      SEQ         <- apply(.fac_to_num(seqs), 1L, .num_to_char)
      attributes(SEQ)         <- atts
      attr(SEQ, "N")     <- N <- nrow(seqs)
      attr(SEQ, "Weights")    <- weights
      covars      <- covars[uni.ind,, drop=FALSE]
      dis.agg     <- rep(seq_len(N), counts)[agg.id]
    }
  } 
  if(!do.uni || !ctrl$do.wts)  {
    uni.ind       <- rep(TRUE, N)
    w2            <- weights
    dist.mat2     <- dist.mat
    HAM.mat2      <- HAM.mat
  }
  
  if(any(floor(G) != G)       ||
     any(G         < 1))         stop("'G' must be strictly positive", call.=FALSE)
  if(any(G        >= sum.uni)) {
    if(length(G)   > 1)        { warning(paste0("Removing G values >= the number of ",    ifelse(do.uni, "unique ", " "), "observations\n"), call.=FALSE, immediate.=TRUE)
      G           <- G[G <= sum.uni]
    } else                       stop(paste0("G values must be less than the number of ", ifelse(do.uni, "unique ", " "), "observations"),   call.=FALSE)
  }
  G               <- rG  <- unique(sort(as.integer(G)))
  do.nec          <- ctrl$do.nec
  if(!do.nec      &&
     algo         != "CEM"    &&
     criterion    == "nec"    &&
     !all(G == 1))             { 
    if(verbose)                  message("Forcing 'do.nec' to TRUE as criterion='nec'\n")
    do.nec        <- TRUE
  }
  if(do.nec &&
     algo         != "CEM"    &&
     !any(G == 1))             { 
    if(verbose)                  message("Forcing G=1 models to be fitted for NEC criterion computation\n")
    G             <- rG  <- unique(c(1L, G))  
  }
  if(all(G  == 1)) { if(verbose) message("Silhouettes not computed as only single component models are being fitted\n")
    do.dbs        <- 
    do.asw        <- FALSE
    if(criterion  == "dbs")    { 
      if(isTRUE(verbose))        message("DBS criterion cannot be used to select among only single-component models: defaulting to 'criterion'=\"bic\"\n")
      criterion   <- "bic"
    }
    if(criterion  == "asw")    { 
      if(isTRUE(verbose))        message("ASW criterion cannot be used to select among only single-component models: defaulting to 'criterion'=\"bic\"\n")
      criterion   <- "bic"
    }
    if(do.nec     && 
       !(do.nec   <- FALSE))   {
      if(verbose)                message("Forcing 'do.nec' to FALSE as only single component models are being fit\n")
    }   
    if(criterion  == "nec")      stop("NEC criterion cannot be used to select among only single-component models", call.=FALSE)
  } else if(algo  == "CEM")    { 
    if(verbose)                  message("Density-based silhouettes not computed as the CEM algorithm is employed\n")
    do.asw        <- TRUE
    do.dbs        <- FALSE
    if(criterion  == "dbs")      stop("DBS criterion cannot be used to select among models fitted via CEM", call.=FALSE)
    if(do.nec     &&
       !(do.nec   <- FALSE))   {
      if(verbose)                message("Forcing 'do.nec' to FALSE as models are being fit via CEM\n")
    }  
    if(criterion  == "nec")      stop("NEC criterion cannot be used to select among models fitted via CEM", call.=FALSE)
  } else do.dbs   <- do.asw   <- TRUE
  len.G           <- length(G)
  if(any(G > 1L))  {
    G1            <- any(G == 1L)
    G2            <- any(G == 2L)
    GG            <- any(G  > 2L)
    all.mod       <- if(all(G1, G2, GG)) unique(c(mtg, mt2, mt1)) else if(all(G1, G2)) unique(c(mt2, mt1)) else if(all(G1, GG)) unique(c(mtg, mt1)) else if(all(G2, GG)) unique(c(mtg, mt2)) else if(G2) mt2 else mtg
    all.mod       <- l.meths[l.meths %in% all.mod]
    if(is.element(init.z, 
       c("hc", "kmedoids", "kmodes")))  {
      hcZ         <- if(do.wts) stats::hclust(dist.mat2, method="ward.D2", members=w2) else agnes(dist.mat2, diss=TRUE, method="ward", keep.diss=FALSE, keep.data=FALSE, trace.lev=0L)
    }
    if(!zli.miss)  {
     if(length(z.list) != len.G) stop(paste0("'z.list' must be a list of length ", len.G), call.=FALSE)  
     if(!all(vapply(z.list, ncol, numeric(1L)) == 
             G))                 stop("Each element of 'z.list' must have 'G' columns",    call.=FALSE)
     if(!all(vapply(z.list, nrow, numeric(1L)) == 
             N2))                stop(paste0("Each element of 'z.list' must have N=", N2, " rows"),     call.=FALSE) 
      z.list      <- lapply(seq_along(G), function(g) z.list[[g]][uni.ind,, drop=FALSE])
    }
    if(all(zli.miss, init.z   == 
           "list"))              stop(paste0("'z.list' must be supplied if 'init.z' is set to 'list'"), call.=FALSE)
  } else all.mod  <- l.meths[l.meths %in% mt1]
  multi           <- length(all.mod) > 1
  cvsel           <- ctrl$do.cv
  ctrl$do.cv      <- FALSE
  if(!cvsel       && 
     criterion    == "cv")     { 
    if(verbose)                  message("Forcing 'do.cv' to TRUE as criterion='cv'\n")
    cvsel         <- TRUE
  }
  if(ctrl$numseq  <- any(c("CU", "UU", "CUN", "UUN") %in% all.mod, opti == "mode", ctrl$ordering != "none")) {
    numseq        <- unname(sapply(SEQ, .char_to_num))
    attr(numseq, "T")    <- P
    attr(numseq, "V")    <- V
  } else numseq   <- NULL
  
  BICs            <-
  ICLs            <-
  AICs            <-
  DBSs            <- 
  ASWs            <- 
  NECs            <- 
  LL.x            <- 
  DF.x            <-
  IT.x            <-
  Nzero.x         <-
  Ninfty.x        <- provideDimnames(matrix(NA, nrow=len.G, ncol=length(all.mod)), base=list(as.character(G), all.mod))
  ZS              <- 
  DBSvals         <- 
  ASWvals         <- replicate(len.G, list())
  if(isTRUE(cvsel))         {
    nfolds        <- pmin(ctrl$nfolds, N)
    if(length(nfolds) != 1 ||
       !is.numeric(nfolds) ||
       (nfolds    <= 1     ||
        floor(nfolds) !=
        nfolds))                 stop("'nfolds' must be a single integer > 1 if cross-validated likelihood is invoked", call.=FALSE)
    CV.x          <- BICs
    cv.ind        <- sample(N)
    cv.SEQ        <- SEQ[cv.ind]
    gCV           <- covars[cv.ind,,        drop=FALSE]
    hmCV          <- HAM.mat[cv.ind,cv.ind, drop=FALSE]
    if(ctrl$numseq)         {
      cv.numseq   <- numseq[,cv.ind, drop=FALSE]
    }
    if(ctrl$do.wts)         {
      cv.wts      <- weights[cv.ind]
    }
    cv.folds      <- cut(seq_len(N), breaks=nfolds, labels=FALSE)
    fcount        <- tabulate(cv.folds, nfolds)
    Nfcount       <- N - fcount
    foldseq       <- seq_len(nfolds)
  }
  crit.tx         <-
  crit.gx         <- -sqrt(.Machine$double.xmax)
  noise           <- all.mod %in% c("CCN", "UCN", "CUN", "UUN")
  nonoise         <- any(!noise)
  noise           <- any(noise)

  gate.G          <- matrix(rG > 1  & gate.x, nrow=2L, ncol=length(rG), byrow=TRUE)
  if(gate.x)       {
    Gn            <- G - !noise.gate
    if(gate.x     &&
      ((any(Gn    <= 1)  && noise) ||
       any(G      <= 1))) {      
      if(verbose)                message(paste0("Can't include gating network covariates ", ifelse(noise.gate, "in a single component mixture", "where G is less than 3 when 'noise.gate' is FALSE\n")))
     gate.G[2L,Gn <= 1]  <- FALSE
    }
  } else           {
    Gn            <- G
    gating        <- stats::as.formula(z ~ 1)
    environment(gating)  <- environment()
  }
  noise.gate      <- !gate.G   | noise.gate
  if(all(equalPro, gate.x))    { 
    if(verbose)                  message("Can't constrain mixing proportions to be equal when gating covariates are supplied\n")
    equalPro      <- FALSE
  }
  equal.tau       <- rbind(G  == 1 | equalPro, Gn < 1 | equalPro) & !gate.G
  equal.n0        <- (rbind(G == 1, Gn == 1) | equalNoise) & equal.tau
  attr(covars, "Gating") <- gate.names
  if(!identical(gating, 
                .drop_constants(covars, 
                gating)))        stop("Constant columns exist in gating formula; remove offending gating covariate(s) and try again", call.=FALSE)
  G.last          <- G[len.G]
  MLRconverge     <- TRUE
  if(isTRUE(verbose))            message("\n################\n")
  
  for(g  in G)     {
    if(isTRUE(verbose))     {    message(paste0("\n", g, " cluster model", ifelse(multi, "s", ""), " -\n"))
      last.G      <- g   == G.last
    }
    if(ctrl$numseq)         {
      attr(numseq,   "G")  <- g
    } 
    attr(SEQ, "G")         <- g
    h             <- which(G   == g)
    g0 <- attr(SEQ, "G0")  <- g - noise
    if(isTRUE(noise))   {
      tau0        <- if(tmiss) 1/g else ctrl$tau0
      if(length(tau0)   > 1)   {
        tau0      <- tau0[uni.ind]
        if(anyNA(tau0))          stop("Invalid 'tau0' supplied", call.=FALSE)
        if(all(gate.G[2L,h]   && 
           noise.gate[2L,h]))  {
          if(N != length(tau0))  stop(paste0("'tau0' must be a scalar or a vector of length N=", N), call.=FALSE)
        } else                   stop("'tau0' must be a scalar in the interval (0, 1)", call.=FALSE)
      }
    }
    
    if(g  > 1)     {
      algog       <- algo
      if(init.z   == "random" &&
         nstarts   > 1)     {
        if(isTRUE(nonoise)) {
          zg      <- replicate(nstarts, list(.unMAP(sample(seq_len(g),  size=N, replace=TRUE), groups=seq_len(g))))
        }
        if(isTRUE(noise))   {
          if(g0    > 1)     {
            zg0   <- replicate(nstarts, list(.unMAP(sample(seq_len(g0), size=N, replace=TRUE), groups=seq_len(g0))))
            zg0   <- lapply(zg0, .tau_noise, tau0)
          } else   {
            zg0   <- matrix(tau0, nrow=N, ncol=2L)
          }
        }
      } else       {
        if(isTRUE(nonoise)) {
          switch(EXPR=init.z, 
                 list=      {
          zg      <- .renorm_z(z.list[[h]])
          },       {
          zg      <- .unMAP(switch(EXPR=init.z, 
                                   random=sample(seq_len(g),  size=N, replace=TRUE),
                                   kmodes=              {
                                     Ms <- seqX[disscenter(dist.mat2, group=stats::cutree(hcZ, k=g),  medoids.index="first", weights=if(do.wts) w2),, drop=FALSE]
                                       suppressWarnings(wKModes(seqs,     weights=if(ctrl$do.wts) weights, random=ctrl$random, modes=Ms, ...)$cluster)
                                     },
                                   kmodes2=suppressWarnings(wKModes(seqs, weights=if(ctrl$do.wts) weights, random=ctrl$random, modes=g,  ...)$cluster),
                                   kmedoids= if(do.wts) {
                                     zz <- wcKMedoids(dist.mat, k=g,  weights=weights, cluster.only=TRUE, initialclust=stats::cutree(hcZ, k=g)[uni.ind])
                                       as.numeric(factor(zz, labels=seq_along(unique(zz))))
                                     } else pam(dist.mat2, k=g,  cluster.only=TRUE, pamonce=pamonce)[uni.ind],
                                   hc=stats::cutree(hcZ, k=g)[uni.ind]), groups=seq_len(g))
          })
        }
        if(isTRUE(noise))   {
          switch(EXPR=init.z, 
                 list=      {
          zg0     <- .renorm_z(z.list[[h]])
          },       {
          if(g0    > 1)     {
            zg0   <- .unMAP(switch(EXPR=init.z, 
                                   random=sample(seq_len(g0), size=N, replace=TRUE),
                                   kmodes=              {
                                     Ms <- seqX[disscenter(dist.mat2, group=stats::cutree(hcZ, k=g0), medoids.index="first", weights=if(do.wts) w2),, drop=FALSE]
                                       suppressWarnings(wKModes(seqs,     weights=if(ctrl$do.wts) weights, random=ctrl$random, modes=Ms, ...)$cluster)
                                     },
                                   kmodes2=suppressWarnings(wKModes(seqs, weights=if(ctrl$do.wts) weights, random=ctrl$random, modes=g0, ...)$cluster),
                                   kmedoids= if(do.wts) {
                                     zz <- wcKMedoids(dist.mat, k=g0, weights=weights, cluster.only=TRUE, initialclust=stats::cutree(hcZ, k=g0)[uni.ind])
                                       as.numeric(factor(zz, labels=seq_along(unique(zz))))
                                     } else pam(dist.mat2, k=g0, cluster.only=TRUE, pamonce=pamonce)[uni.ind],
                                   hc=stats::cutree(hcZ, k=g0)[uni.ind]), groups=seq_len(g0))
            zg0   <- .tau_noise(zg0, tau0)
          } else   {
            zg0   <- matrix(tau0, nrow=N, ncol=2L)
          }
          })
        }
      }
      modtypes    <- if(g  == 2L) mt2 else mtg
    } else         {
      algog       <- "EM"
      zg  <- zg0  <- matrix(1L, nrow=N)
      modtypes    <- mt1
    }
      
    T.last        <- modtypes[length(modtypes)]
    for(modtype   in modtypes) {
      if(isTRUE(verbose))      { message(paste0("\n\tModel: ", modtype, "\n"))
        last.T    <- modtype  == T.last
      }
      ctrl$nmeth  <- is.element(modtype, n.meths)
      ctrl$equalNoise    <- equal.n0[ctrl$nmeth    + 1L,h]
      ctrl$equalPro      <- equal.tau[ctrl$nmeth   + 1L,h]
      ctrl$gate.g        <- gate.G[ctrl$nmeth      + 1L,h]
      ctrl$noise.gate    <- !ctrl$nmeth    || noise.gate[2L, h]
      zm          <- if(attr(SEQ, "Noise") <- gN0 <- ctrl$nmeth) zg0 else zg
      if(gN0      &&   !ctrl$noise.gate    && algog  != "EM") {
        if(init.z == "random" &&
           nstarts > 1)   {
          zm      <- lapply(zm, function(x) { x[,-G] <- replace(x[,-G], x[,-G] > 0, 1L); x })
        } else zm[,-G]   <- replace(zm[,-G], zm[,-G] > 0, 1L)
      }
      m           <- which(modtype == modtypes)
      if(init.z   == "random" &&
         g - gN0   > 1        &&
         nstarts   > 1)        {
        EMX       <- list()
        for(i in startseq)     {
         if(isTRUE(verbose))     message(paste0("\tRandom Start: #", i, "...\n"))
         switch(EXPR=algog,
               cemEM=          {
            ctrl$algo      <- "CEM"
            EMX[[i]]       <- .EM_algorithm(SEQ=SEQ, numseq=numseq, g=g, modtype=modtype, z=zm[[i]],  ctrl=ctrl, gating=gating, covars=covars, HAM.mat = HAM.mat2, MLRconverge = MLRconverge)
            if(!EMX[[i]]$ERR)  {
              ctrl$algo    <- "EM"
              tmpEMX       <- EMX[[i]]
              j.i          <- pmax(tmpEMX$j, 2L)
              EMX[[i]]     <- .EM_algorithm(SEQ=SEQ, numseq=numseq, g=g, modtype=modtype, z=tmpEMX$z, ctrl=ctrl, gating=gating, covars=covars, HAM.mat = HAM.mat2, ll=tmpEMX$ll[c(j.i - 1L, j.i)], MLRconverge = MLRconverge)
              if(EMX[[i]]$ERR) {
                EMX[[i]]   <- tmpEMX
                ctrl$algo  <- "CEM"
              }
            }
         }, EMX[[i]]       <- .EM_algorithm(SEQ=SEQ, numseq=numseq, g=g, modtype=modtype, z=zm[[i]],  ctrl=ctrl, gating=gating, covars=covars, HAM.mat = HAM.mat2, MLRconverge = MLRconverge))
         MLRconverge       <- EMX[[i]]$MLRconverge
        }
        EMX       <- EMX[[which.max(vapply(lapply(EMX, "[[", "ll"), max, numeric(1L)))]]
      } else       {
        switch(EXPR=algog,
              cemEM=        {
          ctrl$algo        <- "CEM"
          EMX              <- .EM_algorithm(SEQ=SEQ, numseq=numseq, g=g, modtype=modtype, z=zm,       ctrl=ctrl, gating=gating, covars=covars, HAM.mat = HAM.mat2, MLRconverge = MLRconverge)
          if(!EMX$ERR)      {
            ctrl$algo      <- "EM"
            tmpEMX         <- EMX
            j.i            <- pmax(EMX$j, 2L)
            EMX            <- .EM_algorithm(SEQ=SEQ, numseq=numseq, g=g, modtype=modtype, z=EMX$z,    ctrl=ctrl, gating=gating, covars=covars, HAM.mat = HAM.mat2, ll=EMX$ll[c(j.i - 1L, j.i)],   MLRconverge = MLRconverge)
            if(EMX$ERR)     {
              EMX          <- tmpEMX
              ctrl$algo    <- "CEM"
            }
          }
        }, EMX             <- .EM_algorithm(SEQ=SEQ, numseq=numseq, g=g, modtype=modtype, z=zm,       ctrl=ctrl, gating=gating, covars=covars, HAM.mat = HAM.mat2, MLRconverge = MLRconverge))
        MLRconverge        <- EMX$MLRconverge
      }
      ERR         <- EMX$ERR
      j           <- EMX$j
      Mstep       <- EMX$Mstep
      ll          <- EMX$ll
      z           <- EMX$z
      cvsel.X     <- cvsel && !ERR
      j2          <- max(1L, j - switch(EXPR=algog, cemEM=1L, 2L))
      if(isTRUE(verbose))        message(paste0("\t\t# Iterations: ", ifelse(ERR, "stopped at ", ""), j2, ifelse(last.G && last.T, "\n\n", "\n")))

      if(all((Mstep$lambda -> lambda) == 0) && cvsel.X) {
        CVll      <- -N * P * log(V)
      } else if(cvsel.X)    {
        ctrl$warn <- FALSE
        lCV       <- numeric(nfolds)
        zCV       <- z[cv.ind,,    drop=FALSE]
        for(i in foldseq)   {
          testX   <- which(cv.folds == i)
          CVS     <- cv.SEQ[testX]
          SCV     <- cv.SEQ[-testX]
          CVz     <- zCV[-testX,, drop=FALSE]
          CVg     <- gCV[-testX,, drop=FALSE]
          attributes(CVS)  <-
          attributes(SCV)  <- attributes(SEQ)
          attr(CVS, "N")   <- fcount[i]
          attr(SCV, "N")   <- Nfcount[i]
          if(ctrl$do.wts)   {
           attr(CVS, "Weights") <- cv.wts[testX]
           attr(SCV, "Weights") <- cv.wts[-testX]
           attr(CVS, "W")  <- sum(attr(CVS, "Weights"))
           attr(SCV, "W")  <- sum(attr(SCV, "Weights"))
          }
          if(any(modtype %in% c("CU", "UU", "CUN", "UUN"), opti == "mode", ctrl$ordering != "none")) {
            nCV            <- cv.numseq[,-testX, drop=FALSE]
            CVn            <- cv.numseq[,testX,  drop=FALSE]
            attr(nCV, "G") <- attr(numseq, "G")
            attr(nCV, "T") <- attr(numseq, "T")
            attr(nCV, "V") <- attr(numseq, "V")
          } else   {
            nCV   <-
            CVn   <- NULL
          }
          EMX     <- .EM_algorithm(SEQ=SCV, numseq=nCV, g=g, modtype=modtype, z=CVz, ctrl=ctrl, gating=gating, covars=CVg, HAM.mat = if(opti != "mode") hmCV[-testX,-testX, drop=FALSE], MLRconverge = MLRconverge)
          MCV     <- EMX$Mstep
          MLRconverge      <- EMX$MLRconverge
          if(is.matrix(MCV$tau)) {
            tau.tmp        <- stats::predict(MCV$fitG, newdata=gCV[testX,, drop=FALSE], type="probs")
            MCV$tau        <- if(ctrl$noise.gate) tau.tmp else cbind(tau.tmp * (1 - MCV$tau[1L,g]), MCV$tau[1L,g])
            rm(tau.tmp)
          }
          MCV$dG  <- NULL
          ctrl$do.cv       <- TRUE
          if(is.infinite(lCV[i] <- 
             .E_step(seqs=CVS, params=MCV, modtype=modtype, ctrl=ctrl, 
                     z.old=if(g  > 1 && ctrl$algo == "CEM" && ctrl$equalPro && !ctrl$random &&
                              modtype %in% c("CC", "CU", "CCN", "CUN")) zCV[testX,, drop=FALSE],
             numseq=CVn)))       break
          ctrl$do.cv       <- FALSE
        }
        CVll      <- sum(lCV)
        ctrl$warn <- TRUE
      }
      
      z           <- if(do.uni) z[dis.agg,, drop=FALSE] else z
      log.lik     <- ll[j]
      nzero       <- sum(lambda == 0)
      ninfty      <- sum(is.infinite(lambda))
      Gfit        <- if(ctrl$gate.g) Mstep$fitG
      gate.pen    <- ifelse(ctrl$gate.g, length(stats::coef(Gfit)) + !ctrl$noise.gate, 
                     ifelse(ctrl$equalPro, as.integer(ctrl$nmeth  && !ctrl$equalNoise), g - 1L))
      choice      <- .choice_crit(ll=log.lik, seqs=SEQ, z=z, modtype=modtype, gp=gate.pen)
      bicx        <- choice$bic
      iclx        <- choice$icl
      aicx        <- choice$aic
      dfx         <- choice$df
      tmp.MAP     <- if(g > 1) max.col(z)
      if(do.dbs   && g > 1) {
        dots      <- list(...)
        if(any(names(dots) == "clusters"))   dots$clusters  <- NULL
        ztol      <- ifelse(any(names(dots)  == "ztol"), dots$ztol, 1E-100)
        summ      <- ifelse(any(names(dots)  == "summ")     && dots$summ == "median", "median", "mean")
        DBS       <- dbs(z, ztol=ztol, weights=if(ctrl$do.wts) w2, MAP=tmp.MAP, summ=summ, clusters=NULL, dots)
        dbsx      <- DBS$wmsw
        DBSvals[[h]][[m]]  <- if(ERR)   NA else DBS$silvals
        attr(DBSvals[[h]][[m]], "G")         <- g
        attr(DBSvals[[h]][[m]], "ModelType") <- modtype
      } else dbsx <- NA
      if(do.asw   && g > 1 && length(unique(tmp.MAP)) > 1)   {
        ASWvals[[h]][[m]]  <- ASW <- if(ERR) NA else cbind(tmp.MAP, wcSilhouetteObs(dist.mat2, tmp.MAP, weights=if(ctrl$do.wts) w2, measure=ifelse(do.wts, "ASWw", "ASW")))
        colnames(ASWvals[[h]][[m]])          <- 
        colnames(ASW)      <- c("cluster", "asw_width")
        summ      <- ifelse(any(names(list(...)) == "summ") && list(...)$summ == "median", "median", "mean")
        aswx      <- ifelse(ERR, NA, ifelse(ctrl$do.wts, 
                                            switch(EXPR=summ, median=weightedMedian(ASW[,2L], w=w2),
                                                                mean=weightedMean(ASW[,2L],   w=w2)), 
                                            switch(EXPR=summ, median=stats::median(ASW[,2L]), mean=mean(ASW[,2L]))))
        attr(ASWvals[[h]][[m]], "G")         <- g
        attr(ASWvals[[h]][[m]], "ModelType") <- modtype
      } else aswx <- NA
      necx        <- ifelse(g > 1 && do.nec, -sum(apply(z, 1L, .entropy))/(log.lik - LL.x[1L,switch(EXPR=modtype, CC=, UC="CC", CU=, UU="CU", "CCN")]), NA)
      crit.t      <- switch(EXPR=criterion, bic=bicx, icl=iclx, aic=aicx, dbs=dbsx, asw=aswx, cv=CVll, nec=necx)
      crit.t      <- ifelse(is.na(crit.t) || ERR, -Inf, crit.t)
      if(crit.t    > crit.tx)     {
        crit.tx   <- crit.t
        theta.x   <- Mstep$theta
        lambda.x  <- lambda
        tau.x     <- Mstep$tau
        z.x       <- z
        ll.x      <- ll
        gp.x      <- gate.pen
        if(gcov.x <- ctrl$gate.g) {
          fit.x   <- Gfit  
        }
      }
      BICs[h,modtype]      <- ifelse(ERR, -Inf, bicx)
      ICLs[h,modtype]      <- ifelse(ERR, -Inf, iclx)
      AICs[h,modtype]      <- ifelse(ERR, -Inf, aicx)
      DBSs[h,modtype]      <- ifelse(ERR, -Inf, dbsx)
      ASWs[h,modtype]      <- ifelse(ERR, -Inf, aswx)
      NECs[h,modtype]      <- ifelse(ERR,  Inf, -necx)
      LL.x[h,modtype]      <- ifelse(ERR, -Inf, log.lik)
      DF.x[h,modtype]      <- ifelse(ERR, -Inf, dfx)
      IT.x[h,modtype]      <- ifelse(ERR,  Inf, j2)
      Nzero.x[h,modtype]   <- ifelse(ERR,  NA,  nzero)
      Ninfty.x[h,modtype]  <- ifelse(ERR,  NA,  ninfty)
      ZS[[h]][[m]]         <- if(ERR)      NA   else z
      if(cvsel)    {
        CV.x[h,modtype]    <- ifelse(ERR, -Inf, CVll)
      }
    } # for (modtype)

    if(crit.tx     > crit.gx)  {
      crit.gx     <- crit.tx
      x.theta     <- theta.x
      x.lambda    <- lambda.x
      x.tau       <- tau.x
      x.z         <- z.x
      x.ll        <- ll.x
      x.gp        <- gp.x
      if(x.gcov   <- gcov.x)   {
        fitG      <- fit.x
      }
    }
    ZS[[h]]       <- stats::setNames(ZS[[h]],               modtypes)
    if(do.dbs     && g > 1)    {
      DBSvals[[h]]         <- stats::setNames(DBSvals[[h]], modtypes)  
    }
    if(do.asw     && g > 1 &&
       length(unique(tmp.MAP)) > 1)           {
      ASWvals[[h]]         <- stats::setNames(ASWvals[[h]], modtypes)  
    }
  } # for (g)

  seqs            <- seqX
  if(any(l.warn   <- (x.ll != cummax(x.ll)))) {
    if(which.max(l.warn)   != 
       length(x.ll))             warning("Log-likelihoods are not strictly increasing\n", call.=FALSE)
  }
  if(any(IT.x[!is.na(IT.x)]
         == ctrl$itmax))         warning(paste0("One or more models failed to converge in the maximum number of allowed iterations (", ctrl$itmax, ")\n"), call.=FALSE)
  class(BICs)     <-
  class(ICLs)     <-
  class(AICs)     <-
  class(DF.x)     <- 
  class(IT.x)     <-
  class(LL.x)     <- "MEDcriterion"
  attr(BICs, "Criterion")  <- "BIC"
  attr(ICLs, "Criterion")  <- "ICL"
  attr(AICs, "Criterion")  <- "AIC"
  attr(DF.x, "Criterion")  <- "DF"
  attr(IT.x, "Criterion")  <- "ITERS"
  attr(LL.x, "Criterion")  <- "loglik"
  attr(LL.x, "Weighted")   <- do.wts
  CRITs           <- switch(EXPR=criterion, bic=BICs, icl=ICLs, aic=AICs, dbs=DBSs, asw=ASWs, cv=CV.x, nec=NECs)
  best.ind        <- which(CRITs == switch(EXPR=criterion, nec=-crit.gx, crit.gx), arr.ind=TRUE)
  if(nrow(best.ind) > 1)    {    warning(paste0("Ties for the optimal model exist according to the '", toupper(criterion), "' criterion: choosing the most parsimonious model\n"), call.=FALSE, immediate.=TRUE)
    best.ind      <- which(DF.x  == min(DF.x[best.ind], na.rm=TRUE), arr.ind=TRUE)
    best.ind      <- best.ind[which.min(best.ind[,1L]),]
  }
  best.G          <- best.ind[1L]
  best.mod        <- colnames(CRITs)[best.ind[2L]]
  G               <- G[best.G]
  x.bic           <- BICs[best.ind]
  x.icl           <- ICLs[best.ind]
  x.aic           <- AICs[best.ind]
  x.dbs           <- DBSs[best.ind]
  x.asw           <- ASWs[best.ind]
  x.nec           <- NECs[best.ind]
  attr(BICs, "G")          <-
  attr(ICLs, "G")          <-
  attr(AICs, "G")          <-
  attr(DF.x, "G")          <-
  attr(IT.x, "G")          <-
  attr(LL.x, "G")          <- rownames(BICs)
  attr(BICs, "modelNames") <-
  attr(ICLs, "modelNames") <-
  attr(AICs, "modelNames") <-
  attr(DF.x, "modelNames") <-
  attr(IT.x, "modelNames") <-
  attr(LL.x, "modelNames") <- colnames(BICs)
  if(cvsel)        {
    CV.x          <- CV.x * 2
    x.cv          <- CV.x[best.ind]
    class(CV.x)   <- "MEDcriterion"
    attr(CV.x, "Criterion")     <- "CV"
    attr(CV.x, "G")             <- rownames(BICs)
    attr(CV.x, "modelNames")    <- colnames(BICs)
    attr(CV.x, "Weighted")      <- 
    attr(x.cv, "Weighted")      <- do.wts
  }
  if(len.G > 1    && verbose)    {
    if(G          == min(rG))    message("Best model occurs at the min of the number of components considered\n")
    if(G          == max(rG))    message("Best model occurs at the max of the number of components considered\n")
  }
  Gseq            <- seq_len(G)
  
  noise           <- best.mod %in% c("CCN", "UCN", "CUN", "UUN")
  attr(x.lambda, "G")           <- G
  attr(x.lambda, "Model")       <- best.mod
  attr(x.lambda, "Names")       <- attr(seqs, "names")
  attr(x.lambda, "Nzero")       <- Nzero.x[best.ind]
  attr(x.lambda, "Ninfty")      <- Ninfty.x[best.ind]
  attr(x.lambda, "P")           <- P
  class(x.lambda)               <- "MEDlambda"
  attr(DF.x,     "Nzero")       <- Nzero.x
  attr(DF.x,     "Ninfty")      <- Ninfty.x
  attr(DF.x,     "Gate.Pen")    <- x.gp
  inds                          <- attr(x.theta, "Ind")
  nonunique                     <- attr(x.theta, "NonUnique")
  if(any(apply(x.lambda == 0, 1L, all))) {
    x.theta                     <- if(G > 1) rbind(do.call(rbind, lapply(x.theta[-G], .char_to_num)), NA)                     else matrix(NA, nrow=1L, ncol=P)
  } else x.theta                <- do.call(rbind, lapply(x.theta, .char_to_num))
  dimnames(x.theta)             <- list(if(G == 1  && noise) "Cluster0" else paste0("Cluster", if(noise) replace(Gseq, G, 0L) else Gseq), colnames(seqs))
  storage.mode(x.theta)         <- "integer"
  attr(x.theta, "alphabet")     <- levs
  attr(x.theta, "labels")       <- attr(seqs, "labels")
  attr(x.theta, "Ind")          <- if(ctrl$opti == "medoid") stats::setNames(inds, paste0("Cluster", seq_len(G - noise)))
  attr(x.theta, "Model")        <- best.mod
  attr(x.theta, "NonUnique")    <- nonunique
  class(x.theta)                <- "MEDtheta"
  colnames(x.z)   <- if(G == 1  && noise) "Cluster0"         else paste0("Cluster", if(noise) replace(Gseq, G, 0L)            else Gseq)
  MAP             <- max.col(x.z)
  MAP             <- if(noise) replace(MAP, MAP == G, 0L)    else MAP
  emptywarn       <- length(unique(MAP))   != G
  if(isTRUE(emptywarn))          warning("Final solution contains at least one empty component\n", call.=FALSE, immediate.=TRUE)
  equalPro        <- equal.tau[1L   + noise,best.G] 
  equalNoise      <- equal.n0[1L    + noise,best.G] 
  noise.gate      <- noise.gate[1L  + noise,best.G]
  noise.gate      <- !noise || noise.gate 
  covars          <- if(do.uni) covars[dis.agg,, drop=FALSE] else covars
  rownames(covars)                 <- seq_along(MAP)
  if(!(gate.G[1L   + noise,best.G] -> bG))  {
    z             <- if(do.wts)       x.z   * WEIGHTS        else x.z
    if(G > (noise.gate   + 1L))  {
      if(ctrl$do.wts)            {
        z[apply(z == 0, 1L, all),] <- .Machine$double.eps
      }
      fitG        <- multinom(gating, trace=FALSE, data=covars, maxit=ctrl$g.itmax, reltol=ctrl$g.tol, MaxNWts=ctrl$MaxNWts)
      MLRconverge <- MLRconverge   && fitG$convergence == 0
      if(equalPro && !equalNoise   && noise) {
        tau0      <- mean(z[,G])
        x.tau     <- c(rep((1 - tau0)/(G - 1L), G  - 1L), tau0)
      } else       {
        x.tau     <- if(equalPro   || equalNoise) rep(1/G, G) else fitG$fitted.values[1L,]
      }
      x.tau       <- stats::setNames(x.tau, paste0("Cluster", if(noise) replace(Gseq, G, 0L)    else Gseq))
      if(equalPro) {
        if(noise  && !equalNoise)   {
          fitG$wts[-(G * 2L)]   <- 0L
        } else fitG$wts[]       <- 0L
        fitG$fitted.values      <- matrix(x.tau, nrow=nrow(z), ncol=G, byrow=TRUE)
        fitG$residuals          <- z - fitG$fitted.values
      }
    }   else       {
      if(any(ctrl$do.wts && G > 1,
             !ctrl$do.wts))         {
       if(G > 1)   {
        z[apply(z == 0, 1L, all),] <- .Machine$double.eps
       }
       fitG       <- suppressWarnings(stats::glm(z ~ 1, family=stats::binomial(), maxit=ctrl$g.itmax))
      } else       {
       z          <- rep(1L, N2)
       fitG       <- suppressWarnings(stats::glm(z ~ 1, family=stats::binomial(), maxit=ctrl$g.itmax, weights=WEIGHTS))
      }
      MLRconverge <- MLRconverge   && isTRUE(fitG$converged)
    }
    attr(fitG, "Formula")       <- "~1"
  }     else       {
    x.tau         <- if(do.uni) x.tau[dis.agg,, drop=FALSE]   else x.tau
  }
  if(isFALSE(MLRconverge))       warning(paste0("\tFor one or more models, in one or more ", algo, " iterations, the multinomial logistic regression\n\t\tin the gating network failed to converge in ", ctrl$g.itmax, " iterations:\n\t\tmodify the 2nd element of the 'itmax' argument to MEDseq_control()\n"), call.=FALSE, immediate.=TRUE)
  if(is.matrix(x.tau))           {
    colnames(x.tau)             <- paste0("Cluster",          if(noise) replace(Gseq, G, 0L)    else Gseq)
  } else x.tau    <- stats::setNames(x.tau, paste0("Cluster", if(noise) replace(Gseq, G, 0L)    else Gseq))
  fitG$lab        <- if(noise   && noise.gate && G > 1) c(paste0("Cluster", Gseq[-G]), "Noise") else if(noise && G > 1) paste0("Cluster", Gseq[-G]) else paste0("Cluster", Gseq)
  attr(fitG, "Data")            <- x.z
  attr(fitG, "Disagg")          <- if(do.uni) dis.agg
  attr(fitG, "DoUni")           <- do.uni
  attr(fitG, "EqualNoise")      <- equalNoise
  attr(fitG, "EqualPro")        <- equalPro
  attr(fitG, "Formula")         <- ifelse(is.null(attr(fitG, "Formula")), Reduce(paste, deparse(gating[-2L])), attr(fitG, "Formula"))
  attr(fitG, "G")               <- G
  attr(fitG, "G0")              <- G == 1     && best.mod   == "CCN"
  attr(fitG, "Maxit")           <- ctrl$g.itmax
  attr(fitG, "MaxNWts")         <- ctrl$MaxNWts
  attr(fitG, "NoGate")          <- if(attr(fitG, "Formula") == "~1") x.tau
  attr(fitG, "Noise")           <- noise
  attr(fitG, "NoiseGate")       <- noise.gate
  attr(fitG, "NoisePro")        <- if(noise   && !noise.gate) ifelse(is.matrix(x.tau), x.tau[1L,G], x.tau[G])
  attr(fitG, "Reltol")          <- ctrl$g.tol
  class(fitG)     <- c("MEDgating", class(fitG))
  if(isTRUE(verbose))            message(paste0("\n\t\tBest Model", ifelse(length(CRITs) > 1, paste0(" (according to ", toupper(criterion), "): "), ": "), best.mod, ", with ",  paste0(G, " component", ifelse(G > 1, "s", "")),
                                         ifelse(bG | x.gcov, paste0(" (incl. ", ifelse(do.wts, "weights and ", ""), "gating network covariates)"), ifelse(do.wts, ifelse(x.ctrl$equalPro && G > 1, " (incl. weights and equal mixing proportions)", " (incl. weights)"), ifelse(x.ctrl$equalPro && G > 1, " (and equal mixing proportions)", ""))), "\n\t\t",
                                         "BIC = ",    round(x.bic, 2L), 
                                         " | ICL = ", round(x.icl, 2L), 
                                         " | AIC = ", round(x.aic, 2L), 
                                         ifelse(any(G > 1 && any(do.dbs, do.asw, do.nec), isTRUE(cvsel)), "\n\t\t", ""),
                                         ifelse(G > 1 && do.dbs, paste0("DBS = ", round(x.dbs, 2L)), ""),
                                         ifelse(G > 1 && do.dbs && any(G > 1 && any(do.asw, do.nec), cvsel), " | ", ""),
                                         ifelse(G > 1 && do.asw, paste0("ASW = ", round(x.asw, 2L)), ""),
                                         ifelse(G > 1 && any(do.dbs, do.asw) && any(G > 1 && do.nec, cvsel), " | ", ""),
                                         ifelse(cvsel,           paste0("CV = ",  round(x.cv,  2L)), ""),
                                         ifelse(G > 1 && any(do.dbs, do.asw, do.nec) && all(cvsel, do.nec),  " | ", ""),
                                         ifelse(G > 1 && do.nec, paste0("NEC = ", round(x.nec, 2L)), ""), "\n\n"))
  params          <- list(theta   = x.theta,
                          lambda  = x.lambda,
                          tau     = x.tau)
  attr(seqs, "Weights")         <- if(do.wts) WEIGHTS                    else rep(1L, N2)
  results         <- list(call    = cX,
                          data    = seqs,
                          modtype = best.mod,
                          G       = G,
                          params  = params,
                          gating  = fitG,
                          z       = x.z,
                          MAP     = MAP,
                          BIC     = BICs,
                          ICL     = ICLs,
                          AIC     = AICs)
  if(do.dbs)       {
    DBSs          <- if(any(rG  == 1)) DBSs[-1L,, drop=FALSE]            else DBSs
    class(DBSs)   <- "MEDcriterion"
    attr(DBSs, "Criterion")     <- "DBS"
    attr(DBSs, "G")             <- if(any(rG  == 1)) rownames(BICs)[-1L] else rownames(BICs)
    attr(DBSs, "modelNames")    <- colnames(BICs)
    DBSvals       <- if(any(rG  == 1)) stats::setNames(DBSvals, rG)[-1L] else stats::setNames(DBSvals, rG)
    attr(DBSs,    "Weighted")   <- 
    attr(DBSvals, "Weighted")   <- do.wts
    attr(DBSs,    "Summ")       <-
    attr(DBSvals, "Summ")       <- ifelse(any(names(list(...)) == "summ"), list(...)$summ, "mean")
    results       <- c(results, list(DBS = DBSs, DBSvals = DBSvals))
    if(G > 1)      {
      x.sils      <- DBSvals[[as.character(G)]][[best.mod]]
      attr(x.dbs,  "Weighted")  <-
      attr(x.sils, "Weighted")  <- do.wts
      attr(x.dbs,  "Summ")      <-
      attr(x.sils, "Summ")      <- ifelse(any(names(list(...)) == "summ"), list(...)$summ, "mean")
      results     <- c(results, list(dbs = x.dbs, dbsvals = x.sils))
    }
  }
  if(do.asw)       {
    ASWs          <- if(any(rG  == 1)) ASWs[-1L,, drop=FALSE]            else ASWs
    class(ASWs)   <- "MEDcriterion"
    attr(ASWs, "Criterion")     <- "ASW"
    attr(ASWs, "G")             <- if(any(rG  == 1)) rownames(BICs)[-1L] else rownames(BICs)
    attr(ASWs, "modelNames")    <- colnames(BICs)
    ASWvals       <- if(any(rG  == 1)) stats::setNames(ASWvals, rG)[-1L] else stats::setNames(ASWvals, rG)
    attr(ASWs,    "Weighted")   <- 
    attr(ASWvals, "Weighted")   <- do.wts
    attr(ASWs,    "Summ")       <-
    attr(ASWvals, "Summ")       <- ifelse(any(names(list(...)) == "summ") && list(...)$summ == "median", "median", "mean")
    results       <- c(results, list(ASW = ASWs, ASWvals = ASWvals))
    if(G > 1)      {
      attr(x.asw,  "Weighted")  <- do.wts
      attr(x.asw,  "Summ")      <- ifelse(any(names(list(...)) == "summ") && list(...)$summ == "median", "median", "mean")
      if(length(unique(MAP))     > 1)          {
       x.sils     <- ASWvals[[as.character(G)]][[best.mod]]
       attr(x.sils, "Weighted") <- attr(x.asw,  "Weighted")
       attr(x.sils, "Summ")     <- attr(x.asw,  "Summ")
       results    <- c(results, list(asw = x.asw, aswvals = x.sils))
      } else       {             warning("'asw' not returned as optimal model has only 1 non-empty component\n", call.=FALSE, immediate.=TRUE)
       results    <- c(results, list(asw = x.asw)) 
      }
    }
  }
  results         <- c(results, list(
                          LOGLIK  = LL.x,
                          DF      = DF.x,
                          ITERS   = IT.x))
  results         <- if(cvsel) c(results, list(CV = CV.x))     else results
  if(do.nec)       {
    NECs          <- NECs[-1L,, drop=FALSE]
    class(NECs)   <- "MEDcriterion"
    attr(NECs, "Criterion")     <- "NEC"
    attr(NECs, "G")             <- rownames(BICs)[-1L]
    attr(NECs, "modelNames")    <- colnames(BICs)
    results       <- c(results, list(NEC = NECs))
  }
  x.ll            <- x.ll[if(G > 1) switch(EXPR=algo, cemEM=-1L, -seq_len(2L)) else 1L]
  attr(x.ll, "Weighted")        <- do.wts
  results         <- c(results, list(
                          bic     = x.bic,
                          icl     = x.icl,
                          aic     = x.aic,
                          loglik  = x.ll,
                          df      = DF.x[best.ind],
                          iters   = IT.x[best.ind]))
  results         <- if(cvsel)           c(results, list(cv  = x.cv))  else results
  results         <- if(do.nec && G > 1) c(results, list(nec = x.nec)) else results
  results         <- c(results, list(
                          ZS      = stats::setNames(ZS, rG),
                          uncert  = if(G > 1) 1 - rowMaxs(x.z) else integer(N2),
                          covars  = covars))
  attr(results, "Algo")         <- algo
  attr(results, "ASW")          <- do.asw
  attr(results, "Counts")       <- if(do.uni) c2               else rep(1L, N2)
  attr(results, "Criterion")    <- criterion
  attr(results, "CV")           <- cvsel
  attr(results, "DBS")          <- do.dbs
  attr(results, "DistMat")      <- HAM.mat
  attr(results, "EmptyWarn")    <- emptywarn
  attr(results, "EqualPro")     <- x.ctrl$equalPro
  attr(results, "EqualNoise")   <- x.ctrl$equalNoise
  attr(results, "Gating")       <- x.gcov | bG
  attr(results, "N")            <- N2
  attr(results, "NEC")          <- do.nec
  attr(results, "Noise")        <- noise
  attr(results, "NoiseGate")    <- x.ctrl$noise.gate
  attr(results, "Opti")         <- opti
  attr(results, "T")            <- P
  attr(results, "Unique")       <- do.uni
  attr(results, "V")            <- V
  attr(results, "Weighted")     <- do.wts
  attr(results, "Weights")      <- replace(attr(seqs, "Weights"), attr(seqs, "Weights") == .Machine$double.eps, 0L)
  class(results)  <- "MEDseq"
  Gnames          <- MEDseq_clustnames(results, cluster=TRUE, size=FALSE)
  attributes(Gnames)            <- NULL
  attr(results, "Gname")        <- 
  attr(results$params$lambda, 
       "Gname")                 <- Gnames
  attr(results$gating, "Gname") <- Gnames[seq_len(G - 1L + (!noise || noise.gate))]
    return(results)
}

#' Plot MEDseq results
#'
#' Produces a range of plots of the results of fitted \code{MEDseq} models.
#' @param x An object of class \code{"MEDseq"} generated by \code{\link{MEDseq_fit}} or an object of class \code{"MEDseqCompare"} generated by \code{\link{MEDseq_compare}}.
#' @param type A character string giving the type of plot requested:
#' \describe{
#' \item{\code{"clusters"}}{Visualise the data set with sequences grouped into their respective clusters. See \code{seriated}. Similar to the \code{type="I"} plot (see below). However, \code{type="clusters"} always plots the hard MAP partition and is unaffected by the \code{soft} argument below.}
#' \item{\code{"central"}}{Visualise the central sequences (typically modal sequences, but this depends on the \code{opti} argument to \code{\link{MEDseq_control}} used during model-fitting). See \code{seriated}. The central sequence for the noise component, if any, is not shown as it doesn't contribute in any way to the likelihood. See the \code{type="ms"} option below for an alternative means of displaying the central sequences.}
#' \item{\code{"precision"}}{Visualise the precision parameters in the form of a heatmap. Values of \code{0} and \code{Inf} are shown in \code{"white"} and \code{"black"} respectively (see \code{quant.scale} and \code{seriated}).}
#' \item{\code{"gating"}}{Visualise the gating network, i.e. the observation index (by default) against the mixing proportions for that observation, coloured by cluster. See \code{seriated}. The optional argument \code{x.axis} can be passed via the \code{...} construct to change the x-axis against which mixing proportions are plotted (only advisable for models with a single gating network covariate, when \code{x.axis} is a quantity related to the gating network of the fitted model).}
#' \item{\code{"bic"}}{Plots all BIC values in a fitted \code{MEDseq} object.}
#' \item{\code{"icl"}}{Plots all ICL values in a fitted \code{MEDseq} object.}
#' \item{\code{"aic"}}{Plots all AIC values in a fitted \code{MEDseq} object.}
#' \item{\code{"dbs"}}{Plots all (weighted) mean/median DBS \emph{criterion} values in a fitted \code{MEDseq} object.}
#' \item{\code{"asw"}}{Plots all (weighted) mean/median ASW \emph{criterion} values in a fitted \code{MEDseq} object.}
#' \item{\code{"cv"}}{Plots all cross-validated log-likelihood values in a fitted \code{MEDseq} object.}
#' \item{\code{"nec"}}{Plots all NEC values in a fitted \code{MEDseq} object.}
#' \item{\code{"LOGLIK"}}{Plots all maximal log-likelihood values in a fitted \code{MEDseq} object.}
#' \item{\code{"dbsvals"}}{Silhouette plot using observations-specific DBS values for the optimal model (coloured by cluster). See \code{seriated}.}
#' \item{\code{"aswvals"}}{Silhouette plot using observations-specific ASW values for the optimal model (coloured by cluster). See \code{seriated}.}
#' \item{\code{"similarity"}}{Produces a heatmap of the similarity matrix constructed from the \code{x$z} matrix at convergence, with observations reordered via \code{seriated} for visual clarity. The (potentially seriated) similarity matrix can also be invisibly returned.}
#' \item{\code{"uncert.bar"}}{Plot the observation-specific clustering uncertainties, if any, in the form of a bar plot.}
#' \item{\code{"uncert.profile"}}{Plot the observation-specific clustering uncertainties, if any, in the form of a profile plot.}
#' \item{\code{"loglik"}}{Plot the log-likelihood at every iteration of the EM/CEM algorithm used to fit the model.}
#' }
#' Also available are the following options which act as wrappers to types of plots produced by the \code{\link[TraMineR]{seqplot}} function in the \pkg{TraMineR} package. All are affected by the value of \code{seriated} and all account for the sampling weights (if any) by default (see the \code{weighted} argument and the related \code{Note} below).
#' 
#' Note also that all of the plot types below can be made to either work with the hard MAP partition, or to use the soft cluster membership probabilities, via the \code{soft} argument below. The soft information is used by default for all but the \code{"i"} and \code{"I"} plot types, which (by default) discard this information to instead use the MAP partition: see the \code{soft} argument below for modifying this default behaviour for all of the following plot types.
#' \describe{
#' \item{\code{"d"}}{State distribution plots (chronograms, by cluster).}
#' \item{\code{"dH"}}{State distribution plots (chronograms, by cluster) with overlaid entropy line as per \code{type="Ht"}. Note that this option is only available if version \code{2.2-4} or later of \pkg{TraMineR} is installed.}
#' \item{\code{"f"}}{Sequence frequency plots (by cluster).}
#' \item{\code{"Ht"}}{Transversal entropy plots (by cluster).}
#' \item{\code{"i"}}{Selected sequence index plots (by cluster). By default, bar widths for each observation will be proportional to their weight (if any). However, this can be overruled by specifying \code{weighted=FALSE}.}
#' \item{\code{"I"}}{Whole set index plots (by cluster). This plot effectively contains almost exactly the same information as \code{type="clusters"} plots, and is similarly affected by the \code{seriated} argument, albeit shown on a by-cluster basis rather than stacked in one plot. However, bar widths for each observation will (by default) be proportional to their weight (if any), which is not the case for \code{type="clusters"} plots. However, this can be overruled by specifying \code{weighted=FALSE}.}
#' \item{\code{"ms"}}{Modal state sequence plots (by cluster). This is an alternative way of displaying the central sequences beyond the \code{type="central"} option above. Notably, this option respects arguments passed to \code{\link{get_MEDseq_results}} via the \code{...} construct (see below), while \code{type="central"} does not, although still nothing is shown for the noise component. \strong{Note}: unlike \code{type="central"}, this option always plots \emph{modal} sequences, even if another \code{opti} setting was invoked during model-fitting via \code{\link{MEDseq_control}}, in which case there will be a mismatch between the visualisation and \code{x$params$theta}. Similarly, there may be a mismatch if \code{soft} and/or \code{weighted} are modified from their default values of \code{TRUE}.}
#' \item{\code{"mt"}}{Mean times plots (by cluster). This is equivalent to plotting the results of \code{\link{MEDseq_meantime}(x, MAP=!soft, weighted=weighted, norm=TRUE, prop=FALSE, map.size=FALSE, wt.size=TRUE)}. Other options for \code{norm=FALSE}, \code{prop=TRUE}, \code{map.size=TRUE}, and \code{wt.size=FALSE} may be added in future versions of this package.}
#' }
#' @param seriated Switch indicating whether seriation should be used to improve the visualisation by re-ordering the \code{"observations"} within clusters (the default), the \code{"clusters"}, \code{"both"}, or \code{"none"}. See \code{\link[seriation]{seriate}} and the \code{smeth} and \code{sortv} arguments below. 
#' 
#' The \code{"clusters"} option (and the cluster-related part of \code{"both"}) is only invoked when \code{type} is one of \code{"clusters"}, \code{"central"}, \code{"precision"}, \code{"gating"}, \code{"dbsvals"}, \code{"aswvals"}, \code{"similarity"}, \code{"d"}, \code{"dH"}, \code{"f"}, \code{"Ht"}, \code{"i"}, \code{"I"}, \code{"ms"}, or \code{"mt"}. 
#' 
#' Additionally, the \code{"observations"} option (and the observation-related part of \code{"both"}) is only invoked when \code{type} is one of \code{"clusters"}, \code{"gating"}, \code{"similarity"}, \code{"i"} or \code{"I"}, which are also the only options for which \code{"both"} is relevant.
#' @param soft This argument is a single logical indicator which is only relevant for the \code{"d"}, \code{"dH"}, \code{"f"}, \code{"Ht"}, \code{"i"}, \code{"I"}, \code{"ms"}, and \code{"mt"} plot types borrowed from \pkg{TraMineR}. When \code{soft=TRUE} (the default for all but the \code{"i"} and \code{"I"} \code{type} plots) the soft cluster membership probabilities are used in a manner akin to \code{\link[WeightedCluster]{fuzzyseqplot}}. Otherwise, when \code{FALSE} (the default for \code{"i"} and \code{"I"} \code{type} plots), the soft information is discarded and the hard MAP partition is used instead. 
#' 
#' Note that soft cluster membership probabilities will not be available if \code{x$G=1} or the model was fitted using the \code{algo="CEM"} option to \code{\link{MEDseq_control}}. Plots may still be weighted when \code{soft} is \code{FALSE}, according to the observation-specific sampling weights, when \code{weighted=TRUE}. Note also that \code{type="Ht"} can be used in conjunction with \code{soft=TRUE}, unlike \code{\link[WeightedCluster]{fuzzyseqplot}} for which \code{type="Ht"} is not permissible. Finally, be advised that plotting may be time-consuming when \code{soft=TRUE} for \code{"i"} and \code{"I"} \code{type} plots.
#' @param weighted This argument is a single logical indicator which is only relevant for the \code{"clusters"}, \code{"central"}, and \code{"precision"} plot types, as well as the \code{"d"}, \code{"dH"}, \code{"f"}, \code{"Ht"}, \code{"i"}, \code{"I"}, \code{"ms"}, and \code{"mt"} plot types borrowed from \pkg{TraMineR}. For plots borrowed from \pkg{TraMineR}, when \code{TRUE} (the default), the weights (if any) are accounted for in such plots. Note that when \code{soft} is \code{TRUE}, plots will still be weighted according to the soft cluster membership probabilities; thus \code{weighted=TRUE} and \code{soft=TRUE} allows both these and the observation-specific weights to be used simultaneously (the default behaviour for both arguments). 
#' 
#' Additionally, for these plots and the \code{"clusters"}, \code{"central"}, and \code{"precision"} types, \code{weighted} is passed through to \code{\link{MEDseq_clustnames}} in the rare case where \code{SPS=TRUE} (see below) and the optional \code{\link{MEDseq_clustnames}} argument \code{size=TRUE} is invoked (again, see below).
#' @param SPS A logical indicating whether clusters should be labelled according to the state-permanence-sequence representation of their central sequence. See \code{\link{MEDseq_clustnames}} and \code{\link[TraMineR]{seqformat}}. Defaults to \code{TRUE} for the plot types adapted from \pkg{TraMineR}, i.e. the \code{"d"}, \code{"dH"}, \code{"f"}, \code{"Ht"}, \code{"i"}, \code{"I"}, \code{"ms"}, and \code{"mt"} \code{type} plots. The \code{SPS} argument is also relevant for the following \code{type} plots: \code{"clusters"}, \code{"central"}, and \code{"precision"}, though \code{SPS} defaults to \code{FALSE} in those instances. Note that if \code{SPS=TRUE} for any relevant plot type, the \code{weighted} argument above is relevant if the optional \code{\link{MEDseq_clustnames}} argument \code{size=TRUE} is invoked (see below).
#' @param smeth A character string with the name of the seriation method to be used. Defaults to \code{"TSP"}. See \code{\link[seriation]{seriate}} and \code{seriation::list_seriation_methods("dist")} for further details and the available methods. Only relevant when \code{seriated != "none"}. When \code{seriated == "obs"} or \code{seriated == "both"}, the ordering of observations can be governed by \code{smeth} or \emph{instead} governed by the \code{sortv} argument below.
#' @param sortv A sorting method governing the ordering of observations for \code{"clusters"}, \code{"gating"}, \code{"similarity"}, \code{"i"}, or \code{"I"} \code{type} plots. Potential options include \code{"dbs"} and \code{"asw"} for sorting observations by their DBS or ASW values (if available). Only relevant if \code{seriated} is one of \code{"observations"} or \code{"both"}. Note that the \code{sortv} argument overrides the setting in \code{smeth} as it pertains to the ordering of observations if \code{sortv} is supplied; otherwise \code{sortv} is \code{NULL} and the \code{smeth} is invoked.
#' 
#' Additionally, when (and only when) \code{soft=TRUE} and \code{type="I"}, the additional option \code{sortv="membership"} is provided in accordance with \code{\link[WeightedCluster]{fuzzyseqplot}}, on which such plots are based.
#' @param subset An optional numeric vector giving the indices of the clusters to be plotted. For models with a noise component, values in \code{0:x$G} are admissible, where \code{0} denotes the noise component, otherwise only values in \code{1:x$G}. Only relevant for the \pkg{TraMineR}-\code{type} plots, i.e. \code{"d"}, \code{"dH"}, \code{"f"}, \code{"Ht"}, \code{"i"}, \code{"I"}, \code{"ms"}, and \code{"mt"} \code{type} plots. Note however, that noise components are never plotted for \code{type="ms"} plots, so \code{subset} values of \code{0} will be ignored in this instance.
#' @param quant.scale Logical indicating whether precision parameter heatmaps should use quantiles to determine non-linear colour break-points when \code{type="precision"}. This ensures each colour represents an equal proportion of the data. The behaviour of \code{0} or \code{Inf} values remains unchanged; only strictly-positive finite entries are affected. Heavily imbalanced values are more likely for the \code{"UU"} and \code{"UUN"} model types, thus \code{quant.scale} defaults to \code{TRUE} in those instances and \code{FALSE} otherwise. Note that \code{quant.scale} is \emph{always} \code{FALSE} for the \code{"CC"} and \code{"CCN"} model types.
#' @param ... Catches unused arguments, and allows arguments to \code{\link{get_MEDseq_results}} to be passed when \code{type} is one of \code{"clusters"}, \code{"dbsvals"}, \code{"aswvals"}, \code{"similarity"}, \code{"uncert.bar"}, \code{"uncert.profile"}, \code{"d"}, \code{"dH"}, \code{"f"}, \code{"Ht"}, \code{"i"}, \code{"I"}, \code{"ms"}, or \code{"mt"}, as well as the \code{x.axis} argument when \code{type="gating"}. Also allows select additional arguments to the \pkg{TraMineR} function \code{\link[TraMineR]{seqplot}} to be used for the relevant plot types (e.g. \code{border} and/or \code{ylab}, \code{serr} where \code{type="mt"}, and \code{info} where \code{type="ms"}) and the \code{size} argument to \code{\link{MEDseq_clustnames}}, where relevant.
#'
#' @return The visualisation according to \code{type} of the results of a fitted \code{MEDseq} model.
#' @details The \code{type} options related to model selection criteria plot values for \emph{all} fitted models in the \code{"MEDseq"} object \code{x}. The remaining \code{type} options plot results for the optimal model, by default. However, arguments to \code{get_MEDseq_results} can be passed via the \code{...} construct to plot corresponding results for suboptimal models in \code{x} when \code{type} is one of \code{"clusters"}, \code{"d"}, \code{"dH"}, \code{"f"}, \code{"Ht"}, \code{"i"}, \code{"I"}, \code{"ms"}, or \code{"mt"}. See the examples below.
#' @note Every \code{type} of plot respects the sampling weights, if any. However, those related to \code{\link[TraMineR]{seqplot}} plots from \pkg{TraMineR} (\code{"d"}, \code{"dH"}, \code{"f"}, \code{"Ht"}, \code{"i"}, \code{"I"}, \code{"ms"}, \code{"mt"}) do so only when \code{weighted=TRUE} (the default). 
#' 
#' For these plot types borrowed from \pkg{TraMineR}, when \code{weighted=TRUE}, the y-axis labels (which can be suppressed using \code{ylab=NA}) always display cluster sizes which correspond to the output of \code{\link{MEDseq_meantime}(x, MAP=!soft, weighted=weighted, map.size=FALSE, wt.size=TRUE)}, where \code{wt.size=TRUE} is \strong{NOT} the default behaviour for \code{\link{MEDseq_meantime}}. 
#' 
#' Please note that \code{type="dH"} will be unavailable if versions of \pkg{TraMineR} prior to \code{2.2-4} are in use. The colour of the entropy line(s) will be \code{"blue"} for \code{type="Ht"} and \code{"black"} for \code{type="dH"}. Finally, the plot types borrowed from \pkg{TraMineR} may be too wide to display in the preview panel. The same may also be true when \code{type} is \code{"dbsvals"} or \code{"aswvals"}. 
#' @references Murphy, K., Murphy, T. B., Piccarreta, R., and Gormley, I. C. (2021). Clustering longitudinal life-course sequences using mixtures of exponential-distance models. \emph{Journal of the Royal Statistical Society: Series A (Statistics in Society)}, 184(4): 1414-1451. <\href{https://rss.onlinelibrary.wiley.com/doi/abs/10.1111/rssa.12712}{doi:10.1111/rssa.12712}>.
#' 
#' Studer, M. (2018). Divisive property-based and fuzzy clustering for sequence analysis. In G. Ritschard and M. Studer (Eds.), \emph{Sequence Analysis and Related Approaches: Innovative Methods and Applications}, pp. 223-239. Cham, Switzerland: Springer International Publishing.
#' 
#' Gabadinho, A., Ritschard, G., Mueller, N. S., and Studer, M. (2011). Analyzing and visualizing state sequences in R with \pkg{TraMineR}. \emph{Journal of Statistical Software}, 40(4): 1-37.
#' @usage 
#' \method{plot}{MEDseq}(x,
#'      type = c("clusters", "central", "precision", "gating", 
#'               "bic", "icl", "aic", "dbs", "asw", "cv", 
#'               "nec", "LOGLIK", "dbsvals", "aswvals", "similarity",
#'               "uncert.bar", "uncert.profile", "loglik", 
#'               "d", "dH", "f", "Ht", "i", "I", "ms", "mt"), 
#'      seriated = c("observations", "both", "clusters", "none"), 
#'      soft = NULL,
#'      weighted = TRUE,
#'      SPS = NULL,
#'      smeth = "TSP",
#'      sortv = NULL,
#'      subset = NULL,
#'      quant.scale = FALSE, 
#'      ...)
#' @author Keefe Murphy - <\email{keefe.murphy@@mu.ie}>
#' @keywords plotting main
#' @export
#' @importFrom matrixStats "rowMaxs" "rowSums2" "weightedMedian" "weightedMean"
#' @importFrom seriation "get_order" "list_seriation_methods" "seriate"
#' @importFrom TraMineR "seqdef" "seqdist" "seqformat" "seqplot"
#' @importFrom WeightedCluster "fuzzyseqplot"
#' @seealso \code{\link{MEDseq_fit}}, \code{\link[TraMineR]{seqplot}}, \code{\link{dbs}}, \code{\link{get_MEDseq_results}}, \code{\link[seriation]{seriate}}, \code{\link[seriation]{list_seriation_methods}}, \code{\link[WeightedCluster]{fuzzyseqplot}}, \code{\link{MEDseq_meantime}}, \code{\link{MEDseq_clustnames}}, \code{\link[TraMineR]{seqformat}}
#'
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
#' 
#' # Fit a range of exponential-distance models without clustering
#' mod0          <- MEDseq_fit(mvad.seq, G=1)
#' 
#' # Show the central sequence and precision parameters of the optimal model
#' plot(mod0, type="central")
#' plot(mod0, type="ms")
#' plot(mod0, type="precision")
#' \donttest{
#' # Fit a range of unweighted mixture models without covariates
#' # Only consider models with a noise component
#' # mod1        <- MEDseq_fit(mvad.seq, G=9:11, modtype=c("CCN", "CUN", "UCN", "UUN"))
#' 
#' # Plot the DBS values for all fitted models
#' # plot(mod1, "dbs")
#' 
#' # Plot the clusters of the optimal model (according to the dbs criterion)
#' # plot(mod1, "clusters", criterion="dbs")
#' 
#' # Use seriation to order the observations and the clusters
#' # plot(mod1, "cluster", criterion="dbs", seriated="both")
#' 
#' # Use a different seriation method
#' # seriation::list_seriation_methods("dist")
#' # plot(mod1, "cluster", criterion="dbs", seriated="both", smeth="Spectral")
#' 
#' # Use the DBS values instead to sort the observations, and label the clusters
#' # plot(mod1, "cluster", criterion="dbs", seriated="both", sortv="dbs", SPS=TRUE, size=TRUE)
#' 
#' # Plot the observation-specific ASW values of the best CCN model (according to the asw criterion)
#' # plot(mod1, "aswvals", modtype="CCN", criterion="asw")
#' 
#' # Plot the similarity matrix (as a heatmap) of the best G=9 model (according to the icl criterion)
#' # plot(mod1, "similarity", G=9, criterion="icl")
#' 
#' # Fit a model with weights and gating covariates
#' # mod2        <- MEDseq_fit(mvad.seq, G=10, modtype="UCN", weights=mvad$weights, 
#' #                           gating=~ fmpr + gcse5eq + livboth, covars=mvad.cov)
#' 
#' # Plot the central sequences & precision parameters of this model
#' # plot(mod2, "central")
#' # plot(mod2, "precision")
#' 
#' # Plot the clustering uncertainties in the form of a barplot
#' # plot(mod2, "uncert.bar")
#' 
#' # Plot the observation-specific DBS values
#' # plot(mod2, "dbsvals")
#' 
#' # Plot the  transversal entropies by cluster & then the state-distributions by cluster
#' # Note that these plots may not display properly in the preview panel
#' # plot(mod2, "Ht", ylab=NA)              # suppress the y-axis labels
#' # plot(mod2, "d", border=TRUE)           # add borders
#' # plot(mod2, "dH", ylab=NA, border=TRUE) # both simultaneously (needs TraMineR >=2.2-4)
#' 
#' # The plots above use the soft cluster membership probabilities
#' # Discard this information and reproduce the per-cluster state-distributions plot
#' # plot(mod2, "d", soft=FALSE)
#' 
#' # The plots above use the observation-specific sampling weights
#' # Discard this information and plot the mean times per state per cluster
#' # plot(mod2, "mt", weighted=FALSE)
#' 
#' # Use type="I" and subset=0 to examine the noise component
#' # plot(mod2, "I", subset=0, border=TRUE, weighted=FALSE, seriated="none")}
plot.MEDseq       <- function(x, type = c("clusters", "central", "precision", "gating", "bic", "icl", "aic", "dbs", "asw", "cv", "nec", "LOGLIK", 
                              "dbsvals", "aswvals", "similarity", "uncert.bar", "uncert.profile", "loglik", "d", "dH", "f", "Ht", "i", "I", "ms", "mt"), 
                              seriated = c("observations", "both", "clusters", "none"), soft = NULL, weighted = TRUE, 
                              SPS = NULL, smeth = "TSP", sortv = NULL, subset = NULL, quant.scale = FALSE, ...) {
  x               <- if(inherits(x, "MEDseqCompare")) x$optimal else x
  if(!missing(type)           &&
     (length(type)       > 1  ||
      !is.character(type)))      stop("'type' must be a character vector of length 1",     call.=FALSE)
  if(!missing(seriated)       &&
     (length(seriated)   > 1  ||
      !is.character(seriated)))  stop("'seriated' must be a character vector of length 1", call.=FALSE)
  type            <- match.arg(type)
  if(type == "dH" && 
     !.version_above("TraMineR", 
                     "2.2-4"))   stop("'type'=\"dH\" is only available if version 2.2-4 or later of the TraMineR package is installed", call.=FALSE)
  seriated        <- match.arg(seriated)
  sericlus        <- is.element(seriated, c("both", "clusters"))     && is.element(type, c("clusters", "central", "precision", "gating", "similarity",
                                                                                           "dbsvals", "aswvals", "d", "dH", "f", "Ht", "i", "I", "ms", "mt"))
  seriobs         <- is.element(seriated, c("both", "observations")) && is.element(type, c("clusters", "gating", "similarity", "i", "I"))
  seriated        <- switch(EXPR=seriated, 
                            clusters=ifelse(sericlus, "clusters", "none"),
                            observations=ifelse(seriobs, "observations", "none"),
                            both=ifelse(sericlus && seriobs, "both", ifelse(sericlus, "clusters", ifelse(seriobs, "observations", "none"))))
  if(seriated     != "none")   {
    if(!missing(smeth)        &&
       (length(smeth)    > 1  ||
        !is.character(smeth)))   stop("'smeth' must be a character vector of length 1",    call.=FALSE)
    smeths        <- list_seriation_methods("dist")
    if(!(smeth  %in% smeths))    stop("Invalid 'smeth': see seriation::list_seriation_methods(\"dist\")", call.=FALSE)
  }
  SPS             <- ifelse(is.null(SPS), !is.element(type, c("clusters", "central", "precision")), SPS)
  if(length(SPS)  != 1        ||
     !is.logical(SPS))           stop("'SPS' must be a single logical indicator",          call.=FALSE)
  palette         <- grDevices::palette()
  savepar         <- graphics::par(no.readonly=TRUE)
  on.exit(graphics::par(savepar))
  on.exit(grDevices::palette(palette), add=TRUE)
  G               <- x$G
  N               <- attr(x, "N")
  P               <- attr(x, "T")
  V       <- mV   <- attr(x, "V")
  noise           <- attr(x, "Noise")
  dmat            <- attr(x, "DistMat")
  Weighted        <- attr(x, "Weighted")
  dat             <- x$data
  modtype         <- x$modtype
  Nseq            <- seq_len(N)
  Pseq            <- seq_len(P)
  symbols         <- c(7, 13, 9, 12, 24, 21, 25, 22)
  use.col         <- c("dimgray", "dodgerblue2", "red3", "green3", "slateblue", "violetred4", "gold", "hotpink")
  alpha.x         <- attr(dat, "alphabet")
  cpal.x          <- attr(dat, "cpal")
  label.x <- lab  <- attr(dat, "labels")
  if(type == "central")        {
    theta         <- x$params$theta
    class(theta)  <- NULL
  }
  if(type == "precision")      {
    lambda        <- x$params$lambda
    lambda        <- switch(EXPR=modtype, CCN=, CUN=rbind(matrix(lambda[1L,], nrow=G - 1L, ncol=P, byrow=modtype == "CUN"), 0L), 
                                          matrix(lambda, nrow=G, ncol=P, byrow=modtype == "CU"))
  }
  if(is.element(type, c("aswvals", "dbsvals", "gating", "precision"))) {
    has_viridis   <- all(isTRUE(suppressMessages(requireNamespace("viridisLite", quietly=TRUE))),
                         isTRUE(.version_above("viridisLite", switch(EXPR=type, precision="0.2.0", "0.4.0"))))
    if(isFALSE(has_viridis))     message(paste0("Nicer colour palettes are invoked for the '", type, "' plot types if the suggested \"viridisLite\" library is installed\n"))
  }
  dots            <- list(...)
  dots            <- dots[unique(names(dots))]
  has.dot         <- length(names(dots)[names(dots) %in% setdiff(names(formals(get_MEDseq_results)), c("x", "what", "..."))]) > 0
  if((has.MAP     <- length(dots) > 0 && any(names(dots)  == "MAP")))             {
    MAP           <- dots$MAP
    if(length(MAP)      != N  ||
       !is.numeric(MAP)       ||
       MAP        != floor(MAP)) stop(paste0("'MAP' must be an integer vector of length N=", N),  call.=FALSE)
  } else MAP      <- x$MAP
  if(is.element(type, c("clusters", "similarity", 
                        "d", "dH", "f", "Ht", "i", "I", "ms", "mt"))   ||
    (isTRUE(sericlus) && 
    (is.element(type, c("central", "precision", "dbsvals", "aswvals")) ||
    (type == "gating" && attr(x, "Gating"))))) {
    if(!is.element(type, c("gating", "central", "precision"))) {
      if(has.MAP)  {
        G         <- length(unique(MAP))
        noise     <- attr(x, "Noise")
      } else if(has.dot)       {
        MAP       <- do.call(get_MEDseq_results, c(list(x=x, what="MAP"), dots[!(names(dots) %in% c("x", "what"))]))
        G         <- attr(MAP, "G")
        noise     <- attr(MAP, "Noise")
      }
    }
    if(type       == "similarity") {
      z           <- if(has.dot) do.call(get_MEDseq_results, c(list(x=x, what="z"), dots[!(names(dots) %in% c("x", "what"))])) else x$z
      G           <- ncol(z)
      sim         <- if(G > 1) tcrossprod(z)      else matrix(1L, nrow=N, ncol=N)
      dmat        <- if(G > 1) 1 - sim            else matrix(0L, nrow=N, ncol=N)
    }
    if(G == 1)     {
      sericlus    <- FALSE  
      seriated    <- switch(EXPR=seriated, clusters="none", both="observations", seriated)
    }
    if(isTRUE(sericlus)       &&
       !all((unip <- 
             unique(MAP))     == 0))           {
      meds        <- NULL
      set.seed(100)
      for(g in sort(unip[unip > 0]))           {
        srows     <- Nseq[MAP == g]
        meds      <- c(meds, srows[which.min(rowSums2(dmat[srows,srows, drop=FALSE]))])
      }
      perm        <- get_order(seriate(stats::as.dist(dmat[meds,meds]), method=smeth))
      perm        <- if(isTRUE(noise)) c(perm, G) else perm
    }
    switch(EXPR=type,
          central= {
       theta      <- theta[perm,,  drop=FALSE]
     }, precision= {
       lambda     <- lambda[perm,, drop=FALSE]
    })
  }
  Gseq            <- seq_len(G)
  perm            <- if(isTRUE(sericlus))  perm                         else Gseq
  perm            <- if(length(perm) != G) c(perm, setdiff(Gseq, perm)) else perm
  perm            <- if(isTRUE(noise))     replace(perm, G, ifelse(is.element(type, c("clusters", "gating", "similarity", "i", "I")), 0L, "Noise"))  else perm
  
  if(is.element(type, c("clusters", "similarity", "i", "I")) ||
    (type == "gating"   && 
     attr(x, "Gating"))) {
    if((vsort     <- isTRUE(seriobs)    && 
       !is.null(sortv)  &&
       !identical(sortv, "membership"))) {
      if(length(sortv)   > 1  ||
         !is.character(sortv))   stop("'sortv' must be a character vector of length 1", call.=FALSE)
      if(type     == "I"      &&
         isTRUE(soft))   {
        if(!is.element(sortv,
           c("dbs", "asw", 
             "membership")))     stop("'sortv' must be one of \"dbs\", \"asw\", or \"membership\"", call.=FALSE)
      } else if(!is.element(sortv, 
           c("dbs", "asw")))     stop("'sortv' must be one of \"dbs\" or \"asw\"", call.=FALSE)
      if(sortv    ==  "dbs"   &&
        (!attr(x, "DBS")      ||
         is.null(x$DBSvals)))    stop(paste0("DBS values cannot be used for sorting as ", ifelse(attr(x, "Algo") == "CEM", "the CEM algorithm was used to fit the models", "only 1-component models were fitted")), call.=FALSE)
      if(sortv    ==  "asw"   &&
        (!attr(x, "ASW")      ||
         is.null(x$ASWvals)))    stop("ASW values cannot be plotted as only 1-component models were fitted", call.=FALSE)
      sortv       <- if(has.dot) do.call(get_MEDseq_results, c(list(x=x, what=toupper(sortv)), dots[!(names(dots) %in% c("x", "what"))])) else switch(EXPR=sortv, dbs=x$dbsvals, x$aswvals)
    } else if(!is.null(sortv)) {
      if(isFALSE(seriobs))       stop("'sortv' can only be supplied when 'seriated' is one of \"observations\" or \"both\"", call.=FALSE)
      if(identical(sortv, "membership")) {
        if(any(type     != "I",
               !isTRUE(soft))) { stop("'sortv'=\"membership\" is only permissible when 'type' is one of \"i\" or \"I\" and 'soft'=TRUE", call.=FALSE)
        } else     {
          vsort   <- seriobs  <- FALSE
        }
      }
    }
    glo.order     <-
    num.cl        <- NULL
    set.seed(200)
    for(g in Gseq) {
      if(!any(MAP == perm[g])) {
        num.cl    <- c(num.cl,  0L)
        next
      }
      srows       <- Nseq[MAP == perm[g]]
      num.cl      <- c(num.cl, length(srows))
      glo.order   <- c(glo.order, if(isTRUE(vsort)) srows[order(sortv[sortv[,1L] == replace(perm, G, G)[g],2L])] else if(isTRUE(seriobs)) srows[get_order(seriate(stats::as.dist(dmat[srows,srows]), method=smeth))] else srows)
    }
    cum.cl        <- cumsum(num.cl)
    gcl           <- c(0L,  cum.cl)
  } else if(!is.null(sortv))     message(paste0("'sortv' is not relevant for the \"", type, "\" plot type\n"))
  
  switch(EXPR=type,
         clusters= {
    graphics::layout(rbind(1, 2), heights=c(0.85, 0.15), widths=1)
    OrderedStates <- data.matrix(.fac_to_num(dat))[glo.order,]
    if(isTRUE(SPS)) {
      tmp_labs    <- MEDseq_clustnames(x, cluster=FALSE, weighted=weighted, ...)[replace(perm, perm == 0, G)]
      graphics::par(mar=c(5.1, 4.1, 4.1, .lab_width(tmp_labs)))
    }
    if(any(num.cl == 0))         warning("Model has one or more empty components\n", call.=FALSE, immediate.=TRUE)
    graphics::image(x=Pseq, z=t(OrderedStates), y=Nseq, zlim=c(1L, length(cpal.x)),
                    axes=FALSE, xlab="Time", ylab="Clusters", col=cpal.x, 
                    main=switch(EXPR=seriated, none="Clusters", observations="Observations Ordered Within Clusters", clusters="Ordered Clusters", "Ordered Clusters and Observations"))
    graphics::box(lwd=2)
    y_labs        <- if(noise) replace(perm, G, ifelse(SPS, 0L, "Noise")) else perm
    y_labs        <- y_labs[num.cl   != 0]
    gcl           <- unique(gcl)
    graphics::axis(side=2,   at=gcl[-length(gcl)] + diff(gcl)/2 + 0.5, lwd=1, las=2, tick=FALSE, cex.axis=0.75, labels=y_labs, line=-0.25)
    graphics::axis(side=1,   at=Pseq, labels=attr(x$data, "names"), cex.axis=0.75)
    if(isTRUE(SPS)) {
      y_labs      <- tmp_labs[num.cl != 0]
      graphics::axis(side=4, at=gcl[-length(gcl)] + diff(gcl)/2 + 0.5, lwd=1, las=2, tick=FALSE, cex.axis=0.67, labels=y_labs, line=-0.5)
    }
    for(g in Gseq) {
      graphics::segments(0, cum.cl[g] + 0.5, P + 0.5, cum.cl[g] + 0.5, lwd=2)
    }
    graphics::par(mar=c(1, 1, 0.5, 1) + 0.1, xpd=FALSE)
    graphics::plot.new()
    graphics::legend("bottom", fill=cpal.x, legend=label.x, ncol=ceiling(V/ifelse(V > 6, 3, 2)), cex=0.75)
    graphics::layout(1)
      invisible()
  }, central=      {
    missind       <- setdiff(seq_len(V), which(tabulate(theta, nbins=V) == 0))
    if(noise      <- is.element(modtype, c("CCN", "UCN", "CUN", "UUN"))) {
      lab         <- c(label.x, expression(paste(lambda, " = 0")))
      mV          <- V   + 1L
      theta[G,]   <- NA
    }
    l.ncol        <- ceiling(mV/ifelse(mV > 6, 3, 2))
    dat           <- suppressMessages(seqdef(theta, states=alpha.x[missind], labels=label.x[missind], cpal=cpal.x[missind]))
    attr(dat, "names")  <- attr(x$data, "names")
    graphics::layout(rbind(1, 2), heights=c(0.85, 0.15), widths=1)
    if(isTRUE(SPS)) {
      tmp_labs    <- MEDseq_clustnames(x, cluster=FALSE, weighted=weighted, ...)[as.numeric(replace(perm, G, G))]
      graphics::par(mar=c(5.1, 4.1, 4.1, .lab_width(tmp_labs)))
    }
    seqplot(dat, type="I", with.legend=FALSE, main="Central Sequences Plot", border=NA, missing.color=graphics::par()$bg, 
            yaxis=FALSE, cex.axis=0.75, ylab=switch(EXPR=seriated, clusters=, both="Ordered Clusters", "Clusters"), xlab="Time")
    if(G > 1)        graphics::axis(2, at=seq_len(G) - 0.5, lwd=1, las=2, tick=FALSE, cex.axis=0.75,
                                    line=-0.25, labels=as.character(replace(perm, perm == "Noise", ifelse(isTRUE(SPS), 0L, "Noise"))))
    if(G > 1 && SPS) graphics::axis(4, at=seq_len(G) - 0.5, lwd=1, las=2, tick=FALSE, cex.axis=0.67,
                                    line=-0.5,  labels=tmp_labs)
    graphics::par(mar=c(1, 1, 0.5, 1) + 0.1, xpd=FALSE)
    graphics::plot.new()
    graphics::legend("bottom", fill=if(isTRUE(noise)) c(cpal.x, graphics::par()$bg) else cpal.x, legend=lab, ncol=l.ncol, cex=0.75)
    graphics::layout(1)
      invisible()
  }, precision=    {
    quant.scale   <- ifelse(missing(quant.scale), is.element(modtype, c("UU", "UUN")), !is.element(modtype, c("CC", "CCN")) && quant.scale)
    if(length(quant.scale)  > 1  ||
       !is.logical(quant.scale)) stop("'quant.scale' must be a single logical indicator", call.=FALSE)
    if(all("CCN"  == modtype     &&
           G      == 1))         message("Nothing to plot!\n")
    graphics::layout(rbind(1, 2), heights=c(0.85, 0.15))
    graphics::par(mar=c(4.1, 4.1, 4.1, ifelse(isTRUE(quant.scale), 5.1, 3.1)))
    if(isTRUE(SPS)) {
      tmp_labs    <- MEDseq_clustnames(x, cluster=TRUE, weighted=weighted, ...)[as.numeric(replace(perm, G, G))]
      graphics::par(mar=replace(graphics::par()$mar, 2L, .lab_width(tmp_labs)))
    }
    i.ind         <- is.infinite(lambda)
    num.ind       <- !i.ind &  0  < lambda
    dat           <- lambda[num.ind]
    ncols         <- 13L
    ncols         <- ifelse(length(unique(dat)) > 0, ifelse(length(unique(dat)) > 1, pmin(length(unique(dat)) + 1L, ncols), ncols), 1L)
    if(isTRUE(quant.scale))       {
      breaks      <- stats::quantile(dat, probs=seq(0L, 1L, length.out=ncols))
      breaks      <- unique(breaks)
      ncols       <- length(breaks)
      facs        <- if(ncols > 1) cut(dat, breaks,     include.lowest=TRUE) else 1L
    } else facs   <- if(ncols > 1) cut(dat, ncols - 1L, include.lowest=TRUE) else 1L
    ncols         <- pmax(1L, ncols - 1L)
    cols          <- if(isTRUE(has_viridis)) viridisLite::viridis(ncols, option="D", direction=-1) else grDevices::heat.colors(13L, rev=TRUE)
    cmat          <- matrix("", nrow=G, ncol=P)
    cmat[i.ind]            <- "black"
    cmat[lambda   == 0]    <- "White"
    cmat[num.ind]          <- cols[as.numeric(facs)]
    levels        <- sort(unique(as.vector(cmat)))
    z             <- matrix(unclass(factor(cmat, levels=levels, labels=seq_along(levels))), nrow=P, ncol=G, byrow=TRUE)
    graphics::image(Pseq, Gseq, z, col=levels, axes=FALSE, xlab="Positions", ylab=ifelse(isTRUE(SPS), "", switch(EXPR=seriated, clusters=, both="Ordered Clusters", "Clusters")), main=paste0("Precision Parameters Plot", ifelse(quant.scale, "\n(Quantile Scale)",  "")))
    graphics::axis(1, at=Pseq, tick=FALSE, las=1, cex.axis=0.75,                            labels=colnames(x$params$theta))
    graphics::axis(2, at=Gseq, tick=FALSE, las=1, cex.axis=ifelse(isTRUE(SPS), 0.67, 0.75), labels=if(isTRUE(SPS)) tmp_labs else as.character(perm), line=ifelse(isTRUE(SPS), -0.67, -0.25))
    graphics::box(lwd=2)
    bx            <- graphics::par("usr")
    xpd           <- graphics::par()$xpd
    box.cx        <- c(bx[2L] + (bx[2L]  - bx[1L])/1000, bx[2L] + (bx[2L] - bx[1L])/1000 + (bx[2L] - bx[1L])/50)
    box.cy        <- c(bx[3L],   bx[3L])
    box.sy        <- (bx[4L]  -  bx[3L]) / length(cols)
    xx            <- rep(box.cx, each=2L)
    graphics::par(xpd = TRUE)
    for(i in seq_along(cols)) {
      yy          <- c(box.cy[1L] + (box.sy * (i - 1L)),
                       box.cy[1L] + (box.sy * (i)),
                       box.cy[1L] + (box.sy * (i)),
                       box.cy[1L] + (box.sy * (i - 1L)))
      graphics::polygon(xx, yy, col = cols[i], border = cols[i])
    }
    graphics::par(new=TRUE)
    base::plot(0, 0, type="n", ylim=if(isTRUE(quant.scale) || length(unique(dat)) == 0) c(0L, 100L) else range(dat, na.rm=TRUE), yaxt="n", ylab="", xaxt="n", xlab="", frame.plot=FALSE)
    if(isTRUE(quant.scale))   {
      graphics::axis(side=4, las=2, tick=FALSE, line=-0.125, cex.axis=0.725,
                     at=seq(0L, 100L, length.out=ncols), labels=.tidy_breaks(breaks))
    } else           graphics::axis(side=4, las=2, tick=FALSE, line=0.125)
    graphics::par(mar=c(0, 0, 0, 0))
    graphics::plot.new()
    graphics::legend("center", c(expression(paste(lambda, " = 0")), expression(paste(lambda %->%~infinity))), fill=c("white", "black"), ncol=2, text.width=0.1, cex=1.25)
    graphics::layout(1)
      invisible()
  }, gating=       {
    suppressWarnings(graphics::par(pty="m"))
    Tau           <- .mat_byrow(x$params$tau, nrow=N, ncol=ncol(x$z))
    sericlus      <- isTRUE(sericlus)   && attr(x, "Gating")
    seriobs       <- isTRUE(seriobs)    && attr(x, "Gating")
    perm          <- replace(perm, G, G)
    Tau           <- if(isTRUE(sericlus))  Tau[,perm,      drop=FALSE] else Tau
    vars          <- all.vars(stats::as.formula(attr(x$gating, "Formula")))
    if(miss.x     <- length(dots) > 0   && any(names(dots) == "x.axis")) {
      ncovs       <- length(vars) > 1
      if(isTRUE(ncovs))          warning("Function may produce undesirable plot when 'x.axis' is supplied for a model with multiple gating network covariates\n", call.=FALSE, immediate.=TRUE)
      x.axis      <- dots$x.axis
      if(length(x.axis) != N)    stop(paste0("'x.axis' not of length N=", N), call.=FALSE)
      o.axis      <- order(x.axis, decreasing=FALSE)
      x.axis      <- x.axis[o.axis]
      Tau         <- Tau[o.axis,,        drop=FALSE]
      type        <- ifelse(ncovs, "p", "b")
    } else         {
      Tau         <- if(isTRUE(seriobs))   Tau[glo.order,, drop=FALSE] else Tau
      x.axis      <- seq_len(N)
      type        <- "b"
    }
    xlab          <- ifelse(miss.x, ifelse(is.null(dots$xlab), deparse(match.call()$x.axis), dots$xlab), ifelse(isTRUE(seriobs), "Seriated Observations", ifelse(isTRUE(sericlus), "Observations", "Observation")))
    if(isTRUE(has_viridis)) {
      col         <- if(noise) c(viridisLite::viridis(G - 1L, option="H", direction=1), "grey65") else viridisLite::viridis(G, option="H", direction=1)  
    } else col    <- if(noise) c(grDevices::rainbow(G - 1L), "grey65") else grDevices::rainbow(G)
    col           <- col[perm]
    if(length(x.axis) != N)      stop("'x.axis' must be of length N", call.=FALSE)
    if(x.fac      <- is.factor(x.axis)) {
      xlev        <- levels(x.axis)
      x.axis      <- as.integer(x.axis)
      xaxt        <- "n"
    } else         {
      type        <- ifelse(any(vars %in% names(x$gating$xlevels)), "p", type)
      xaxt        <- ifelse(any(seriobs, sericlus), "n", "s")
    }
    graphics::matplot(x=x.axis, y=Tau, type=type, main="Gating Network", xaxt=xaxt, xlab=xlab, ylab="", col=col, pch=1, lty=perm)
    graphics::mtext(expression(widehat(tau)[g]), side=2, las=2, line=3)
    if(x.fac)      {
      graphics::axis(1, at=unique(x.axis), labels=xlev)
    }
    if(isTRUE(sericlus)    &&
       isFALSE(miss.x))     {
      graphics::abline(v=cum.cl)
      graphics::mtext(replace(perm, G, "Noise"), at=gcl[-length(gcl)] + diff(gcl)/2, side=1, las=1)
      graphics::mtext("Ordered Clusters", side=1, line=1, las=1)
    }
      invisible()
  }, bic=,
     icl=,
     aic=,
     dbs=,
     asw=,
     nec=,
     cv=,
     LOGLIK=       {
    if(all(type   == "cv",
       !attr(x, "CV")))          stop("Cross-validated log-likelihood values cannot be plotted as cross-validation didn't take place during model fitting\n", call.=FALSE)
    dat           <- switch(EXPR=type, bic=x$BIC, icl=x$ICL, aic=x$AIC, dbs=x$DBS, asw=x$ASW, cv=x$CV, nec=x$NEC, LOGLIK=x$LOGLIK)
    if(type ==  "dbs"      &&
       (!attr(x, "DBS")     ||
        is.null(dat)))            stop(paste0("DBS values cannot be plotted as ", ifelse(attr(x, "Algo") == "CEM", "the CEM algorithm was used to fit the models", "only 1-component models were fitted")), call.=FALSE)
    if(type ==  "asw"      &&
       (!attr(x, "ASW")     ||
        is.null(dat)))            stop("ASW values cannot be plotted as only 1-component models were fitted", call.=FALSE)
    if(type ==  "nec"     &&
      (!attr(x, "NEC")     || 
       is.null(dat)))            stop(paste0("NEC values cannot be plotted as ", ifelse(attr(x, "NEC"), "only 1-component models were fitted", "'do.nec' was set to FALSE")), call.=FALSE)
    ms            <- which(c("CC", "UC", "CU", "UU", "CCN", "UCN", "CUN", "UUN") %in% colnames(dat))
    symbols       <- symbols[ms]
    use.col       <- use.col[ms]
    graphics::matplot(dat, type="b", xlab="Number of Components (G)", ylab=switch(EXPR=type, dbs=, asw=, cv=, LOGLIK="", toupper(type)), col=use.col, bg=use.col, pch=symbols, ylim=range(as.vector(dat[!is.na(dat) & is.finite(dat)])), xaxt="n", lty=1)
    if(type == "dbs")    graphics::mtext(paste0(ifelse(Weighted, "Weighted ", ""), switch(EXPR=attr(dat, "Summ"), median="Median", "Mean"), " DBS"), side=2, line=3, las=3)
    if(type == "asw")    graphics::mtext(paste0(ifelse(Weighted, "Weighted ", ""), switch(EXPR=attr(dat, "Summ"), median="Median", "Mean"), " ASW"), side=2, line=3, las=3)
    if(type == "cv")     graphics::mtext(ifelse(Weighted, expression("\u2113"["cv"]^"w"), expression("\u2113"["cv"])), side=2, line=3, las=1, cex=1.5)
    if(type == "LOGLIK") graphics::mtext(paste0(ifelse(Weighted, "Weighted ", ""), "Log-Likelihood"), side=2, line=3, las=3)
    graphics::axis(1, at=seq_len(nrow(dat)), labels=rownames(dat))
    graphics::legend(switch(EXPR=type, nec=, dbs="topright", "bottomright"), ncol=2, cex=1, inset=0.01, legend=colnames(dat), pch=symbols, col=use.col, pt.bg=use.col)
      invisible()
  }, dbsvals=,
     aswvals=      {
    switch(EXPR=type,
          dbsvals= {
      if(!attr(x, "DBS")     ||
         all(type   == "dbs",
         is.null(x$DBSvals)))    stop(paste0("DBS values cannot be plotted as ", ifelse(attr(x, "Algo") == "CEM", "the CEM algorithm was used to fit the model", "only 1-component models were fitted")), call.=FALSE)
      object      <- if(has.dot) do.call(get_MEDseq_results, c(list(x=x, what="DBS"), dots[!(names(dots) %in% c("x", "what"))])) else x$dbsvals
    },    aswvals= {
      if(!attr(x, "ASW")     ||
         all(type   == "asw",
         is.null(x$ASWvals)))    stop("ASW values cannot be plotted as only 1-component models were fitted",      call.=FALSE)
      object      <- if(has.dot) do.call(get_MEDseq_results, c(list(x=x, what="ASW"), dots[!(names(dots) %in% c("x", "what"))])) else x$aswvals
      if(is.null(object))        stop("ASW values cannot be plotted as the model only has 1 non-empty component", call.=FALSE)
    })
    rownames(object)       <- as.character(Nseq)
    G             <- attr(object, "G")
    modtype       <- attr(object, "ModelType")
    cl            <- object[,"cluster"]
    X             <- object[order(cl, -object[,2L]),, drop=FALSE]
    if(isTRUE(sericlus))    {
      perm        <- replace(perm, G, G)
      X           <- do.call(rbind, lapply(Gseq, function(g) X[X[,1L] == perm[g],]))
      col         <- rev(rep(unique(X[,1L]), table(X[,1L])[perm]))
    } else col    <- rev(X[,1L])
    sil           <- X[,2L]
    space         <- c(0L, rev(diff(cli <- X[,"cluster"])))
    space[space   != 0]    <- 0.5
    ng            <- table(cl)
    if(length(ng) != G)          warning("Model has one or more empty components", call.=FALSE, immediate.=TRUE)
    noise         <- modtype %in% c("CCN", "UCN", "CUN", "UUN")
    if(identical(palette, grDevices::palette("default"))) {
      palette     <- if(isTRUE(has_viridis)) viridisLite::viridis(G, option="H", direction=1) else rep(c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999"), length.out=G)
      palette[G]  <- ifelse(noise, "black", palette[G])
      grDevices::palette(if(grDevices::dev.capabilities()$semiTransparency) grDevices::adjustcolor(palette, alpha.f=0.75) else palette)
    }
    dat           <- rev(graphics::barplot(rev(sil), space=space, xlim=c(min(0L, min(sil)), 1L), 
                                           horiz=TRUE, las=1, mgp=c(2.5, 1, 0), col=col, border=NA, 
                                           cex.names=graphics::par("cex.axis"), axisnames=FALSE))
    summ          <- attr(object, "Summ")
    if(Weighted)   {
      weights     <- attr(x, "Weights")[order(cl, -object[,2L])]
      wSIL        <- switch(EXPR=summ, median=weightedMedian(sil, weights), mean=weightedMean(sil, weights))
      switch(EXPR=summ, median=graphics::title(main=paste0("(Weighted) ", switch(EXPR=type, dbsvals="Density-based ", ""), "Silhouette Plot", ifelse(isTRUE(sericlus), " (Ordered by Cluster)", "")), sub=paste0("(Weighted) Median ", switch(EXPR=type, dbsvals="DBS", "Silhouette"), " Width : ", round(wSIL, digits=3)), adj=0),
                          mean=graphics::title(main=paste0("(Weighted) ", switch(EXPR=type, dbsvals="Density-based ", ""), "Silhouette Plot", ifelse(isTRUE(sericlus), " (Ordered by Cluster)", "")), sub=paste0("(Weighted) Mean ",   switch(EXPR=type, dbsvals="DBS", "Silhouette"), " Width : ", round(wSIL, digits=3)), adj=0))
    } else         {
      wSIL        <- switch(EXPR=summ, median=stats::median(sil), mean=mean(sil))
      switch(EXPR=summ, median=graphics::title(main=paste0(switch(EXPR=type, dbsvals="Density-based ", ""), "Silhouette Plot"), sub=paste0("Median ", switch(EXPR=type, dbsvals="DBS", "Silhouette"), " Width : ", round(wSIL, digits=3)), adj=0),
                          mean=graphics::title(main=paste0(switch(EXPR=type, dbsvals="Density-based ", ""), "Silhouette Plot"), sub=paste0("Mean ",   switch(EXPR=type, dbsvals="DBS", "Silhouette"), " Width : ", round(wSIL, digits=3)), adj=0))
    }
    graphics::mtext(paste0("n = ", N),  adj=0)
    graphics::mtext(substitute(G~modtype~"clusters"~C[g], list(G=G, modtype=modtype)), adj=1)
    if(Weighted)   {
      switch(EXPR=summ, median=graphics::mtext(expression(paste(g, " : ", n[g], " | ", w.med[i %in% Cg] ~ ~w[i]~s[i])), adj=1.05, line=-1.2),
                          mean=graphics::mtext(expression(paste(g, " : ", n[g], " | ", w.avg[i %in% Cg] ~ ~w[i]~s[i])), adj=1.05, line=-1.2))
      silcl       <- split(sil,  cli)
      meds        <- switch(EXPR=summ, median=vapply(seq_along(silcl), function(i) weightedMedian(silcl[[i]], weights[cli == i]), numeric(1L)),
                                         mean=vapply(seq_along(silcl), function(i) weightedMean(silcl[[i]],   weights[cli == i]), numeric(1L)))
    } else         {
      switch(EXPR=summ, median=graphics::mtext(expression(paste(g, " : ", n[g], " | ", med[i %in% Cg] ~ ~s[i])), adj=1.05, line=-1.2),
                          mean=graphics::mtext(expression(paste(g, " : ", n[g], " | ", avg[i %in% Cg] ~ ~s[i])), adj=1.05, line=-1.2))
      meds        <- tapply(sil, cli, switch(EXPR=summ, median=stats::median, mean=mean))
    }
    medy          <- tapply(dat, cli, stats::median)
    for(g in seq_len(G)) {
      graphics::text(1, medy[g], paste(ifelse(g == G && noise, "Noise", g), ":  ", ng[g], " | ", format(meds[g], digits=1, nsmall=2)), xpd=NA, adj=0.8)
    }
    w.r           <- cumsum(space + 1L)
    w.l           <- w.r - 1L
    graphics::segments(x0=wSIL, y0=min(w.l), y1=max(w.r), lty=2)
      invisible()
  }, similarity=   {
    sim           <- sim[glo.order,glo.order]
    col           <- grDevices::heat.colors(30L, rev=TRUE)
    oldpar        <- suppressWarnings(graphics::par(no.readonly=TRUE))
    oldpar$new    <- FALSE
    on.exit(suppressWarnings(graphics::par(oldpar)))
    graphics::par(mar=c(5.1, 4.1, 4.1, 3.1))
    graphics::image(t(sim[seq(from=N, to=1L, by=-1L),]), col=col, main="Similarity Matrix",
                    xlab=paste0("Observation 1:N", if(seriated != "none") " (Reordered)"),
                    ylab=paste0("Observation 1:N", if(seriated != "none") " (Reordered)"))
    graphics::title(switch(EXPR=seriated, clusters="(Ordered by Clusters)", observations="(Ordered by Observations)", both="(Ordered by Clusters and Observations)", ""), cex.main=0.8, line=0.5)
    graphics::box(lwd=2)
    suppressWarnings(.heat_legend(data=sim, col=col, cex.lab=0.8))
      invisible(if(seriated == "none") sim else provideDimnames(sim, base=list(as.character(glo.order), as.character(glo.order))))
  }, uncert.profile=,
     uncert.bar=   {
    if(attr(x, "Algo") == "CEM") message("No uncertainties to plot: model was fitted via CEM\n")
    graphics::par(pty="m", mar=c(5.1, 4.1, 4.1, 3.1))
    if(has.dot)    {
      z           <- do.call(get_MEDseq_results, c(list(x=x, what="z"), dots[!(names(dots) %in% c("x", "what"))]))
      G           <- ncol(z)
      uncX        <- if(G > 1) 1 - rowMaxs(z) else integer(N)
    } else uncX   <- x$uncert
    oneG          <- 1/G
    min1G         <- 1 - oneG
    yx            <- unique(c(0, pretty(c(0, min1G))))
    yx            <- replace(yx, length(yx), min1G)
    cm            <- c("dodgerblue2", "red3", "green3")
    switch(EXPR=type,
           uncert.bar=         {
             cu   <- cm[seq_len(2L)][(uncX >= oneG) + 1L]
             cu[uncX == 0] <- NA
             base::plot(uncX, type="h", ylim=range(yx), col=cu, yaxt="n", ylab="", xlab="Observations", lend=1)
             graphics::lines(x=c(0, N), y=c(oneG, oneG), lty=2, col=cm[3L])
             graphics::axis(2, at=yx, labels=replace(yx, length(yx), "1 - 1/G"), las=2, xpd=TRUE)
             graphics::axis(2, at=oneG, labels="1/G", las=2, xpd=TRUE, side=4, xpd=TRUE)
           }, uncert.profile=  {
             ord  <- order(uncX, decreasing=FALSE)
             ucO  <- uncX[ord]
             base::plot(ucO, type="n", ylim=c(-max(uncX)/32, max(yx)), ylab="", xaxt="n", yaxt="n", xlab=paste0("Observations in order of increasing uncertainty"))
             graphics::lines(x=c(0, N), y=c(0, 0), lty=3)
             graphics::lines(ucO)
             graphics::points(ucO, pch=15, cex=0.5, col=1)
             graphics::lines(x=c(0, N), y=c(oneG, oneG), lty=2, col=cm[3L])
             graphics::axis(2, at=yx,   las=2, xpd=TRUE, labels=replace(yx, length(yx), "1 - 1/G"))
             graphics::axis(2, at=oneG, las=2, xpd=TRUE, labels="1/G", side=4)
           })
    graphics::mtext("Uncertainty", side=2, line=3)
    graphics::title(main=list(paste0("Clustering Uncertainty ", switch(EXPR=type, uncert.bar="Barplot", "Profile Plot"))))
      invisible()
  }, loglik=       {
    x             <- x$loglik
    if(all(x      != cummax(x))) warning("Log-likelihoods are not strictly increasing\n", call.=FALSE)
    base::plot(x, type=ifelse(length(x) == 1, "p", "l"), xlab="Iterations", ylab=paste0(ifelse(Weighted, "Weighted ", ""), "Log-Likelihood"), xaxt="n", ...)
    seqll         <- seq_along(x)
    llseq         <- pretty(seqll)
    llseq         <- if(any(llseq  != floor(llseq))) seqll else llseq
    graphics::axis(1, at=llseq, labels=llseq)
      invisible()
  },               {
    if(is.null(soft))       {
      soft        <- !is.element(type, c("i", "I"))
    } else if(length(soft) != 1    || 
              !is.logical(soft)) stop("'soft' must be a single logical indicator",        call.=FALSE)
    soft          <- isTRUE(soft)  && (G > 1 && (attr(x, "Algo") != "CEM"))
    if(length(unique(MAP)) != G) warning("Model contains one or more empty components\n", call.=FALSE, immediate.=TRUE)
    if(!is.null(subset)    &&
      (!is.numeric(subset) ||
      !all(subset == 
           floor(subset))  ||
      !all(subset %in% 
           (!noise):G)))    {    stop(paste0("'subset', if supplied, must be a numeric vector with entries in {", 0L + !noise, ", G=", G-noise, "} for models ", ifelse(noise, "with", "without"), " a noise component"), call.=FALSE)
    } else if(!is.null(subset))     {
      subset      <- if(any(subset == 0)) c(sort(subset[subset != 0]), 0L) else sort(subset)
      MAP[!(MAP   %in% 
              subset)]     <- NA 
    }
    MAP           <- MAP2  <- factor(replace(MAP, MAP == 0, "Noise"), levels=switch(EXPR=type, i=, I=replace(perm, perm == 0, "Noise"), perm))
    if(isTRUE(SPS))         {
      levels(MAP) <- MEDseq_clustnames(x, weighted=weighted, ...)
    }
    MAP           <- switch(EXPR=type, i=, I=MAP[glo.order],  MAP)
    dat           <- switch(EXPR=type, i=, I=dat[glo.order,], dat)
    if(length(weighted) != 1  ||
       !is.logical(weighted))    stop("'weighted' must be a single logical indicator",    call.=FALSE)
    weighted      <- weighted && Weighted
    attr(dat, "weights")      <- if(weighted) switch(EXPR=type, i=, I=attr(dat, "Weights")[glo.order], attr(dat, "Weights"))        else rep(1L, N)
    if(type       == "ms")  {
      noisemap    <- which(MAP2    != "Noise")
      if(sum(noisemap)        == 0 &&
         !is.null(subset))       stop("Can't plot \"ms\" type plot as there are no non-noise components available after invoking 'subset'", call.=FALSE)
      if(any(0    == subset))    warning("'subset' value of '0' ignored for type=\"ms\" plots", call.=FALSE, immediate.=TRUE)
      MAP         <- droplevels(MAP[noisemap])
      dat         <- dat[noisemap,, drop=FALSE]
    }
    if(isTRUE(soft))        {
      if(type     == "ms")  {
        if(has.dot)         {
          x$z     <- do.call(get_MEDseq_results, c(list(x=x, what="z"), dots[!(names(dots) %in% c("x", "what"))]))
        }
        x$z       <- if(isTRUE(weighted)) x$z * attr(dat, "Weights")       else x$z
        x$z       <- x$z[noisemap,as.numeric(replace(perm, perm == "Noise", G)),   drop=FALSE]
        attr(dat, "weights")  <- vapply(seq_along(noisemap), function(i) x$z[i,MAP[i]], numeric(1L))
        if((G - noise) == 0)     stop("Nothing to plot: no non-noise components", call.=FALSE)
      } else       {
        group     <- if(has.dot) do.call(get_MEDseq_results, 
                                         c(list(x=x, what="z"), dots[!(names(dots) %in% c("x", "what"))])) else x$z
        group     <- switch(EXPR=type, i=, I=group[glo.order,,       drop=FALSE], group)
        if(!is.null(subset))   {
         sub.ind  <- !is.na(MAP)
         dat      <- dat[sub.ind,,            drop=FALSE]
         MAP      <- droplevels(MAP[sub.ind])
         attr(dat, "weights") <- attr(dat, "Weights")[sub.ind]
         group    <- group[sub.ind, replace(subset, subset == 0, G), drop=FALSE]
         N        <- nrow(dat)
         G        <- length(unique(MAP))
        }
        dat       <- dat[rep(seq_len(N), G),, drop=FALSE]
        MAP       <- factor(rep(seq_len(G), each=N), labels=levels(MAP))
        attr(dat, "weights")  <- if(isTRUE(weighted)) attr(dat, "weights") * as.vector(group)              else as.vector(group)
        if(type   == "I"      && 
           !is.null(sortv)    && 
           identical(sortv, "membership")) {
          dots$sortv          <- as.vector(group)
        }
      }
      weighted    <- TRUE
    }
    attr(dat, "Weights")      <- NULL
    dots          <- c(list(seqdata=dat, group=MAP, with.legend=type != "Ht", with.missing=FALSE, type=type, weighted=weighted), 
                       dots[!(names(dots) %in% c("G", "modtype", "noise", "cluster", "size"))])
    dots          <- switch(EXPR=type, Ht=dots[names(dots) != "border"], if(!any(names(dots) == "border")) c(list(border=NA), dots) else dots)
    if(type       == "dH"     && 
       !("cols" %in% names(dots))) {
      dots        <- c(dots, list(col="black", cols=NA))
    } 
    dots          <- switch(EXPR=type,  I=c(list(space=0L),   dots), dots)
    dots          <- switch(EXPR=type, mt=c(list(prop=FALSE), dots), dots)
    dots          <- dots[unique(names(dots))]
    suppressWarnings(do.call(seqplot, dots))
      invisible()
  })
}

#' @method summary MEDseq
#' @rdname MEDseq_fit
#' @usage
#' \method{summary}{MEDseq}(object,
#'         classification = TRUE,
#'         parameters = FALSE,
#'         network = FALSE,
#'         SPS = FALSE,
#'         ...)
#' @export
summary.MEDseq  <- function(object, classification = TRUE, parameters = FALSE, network = FALSE, SPS = FALSE, ...) {
  object        <- if(inherits(object, "MEDseqCompare")) object$optimal else object
  if(length(classification)  > 1 || 
    !is.logical(classification)) stop("'classification' must be a single logical indicator", call.=FALSE)
  if(length(parameters)   > 1    ||
    !is.logical(parameters))     stop("'parameters' must be a single logical indicator",     call.=FALSE)
  if(length(network)      > 1    ||
    !is.logical(network))        stop("'network' must be a single logical indicator",        call.=FALSE)
  if(length(SPS)          > 1    ||
     !is.logical(SPS))           stop("'SPS' must be a single logical indicator",            call.=FALSE)
  G             <- object$G
  attr(G, "range")       <- eval(object$call$G)
  params        <- object$params
  equalPro      <- G <= 1 || attr(object$gating, "EqualPro")
  equalN        <- attr(object$gating, "EqualNoise") && equalPro
  summ          <- list(data = deparse(object$call$seqs), N = attr(object, "N"), P = attr(object, "T"), V = attr(object, "V"), G = G, modelName = object$modtype, algo=attr(object, "Algo"), 
                        loglik = object$loglik[length(object$loglik)], df = object$df, iters = object$iters, gating = object$gating, bic = unname(object$bic), icl = unname(object$icl), 
                        aic = unname(object$aic), dbs = unname(object$dbs), asw = unname(object$asw), emptywarn = attr(object, "EmptyWarn"), cv = unname(object$cv), nec = unname(object$nec), 
                        tau = params$tau, theta = params$theta, lambda = params$lambda, z = object$z, equalPro = equalPro, equalNoise = equalN, classification = object$MAP, SPS = SPS, 
                        noise.gate = attr(object, "NoiseGate"), gating = object$gating, printClass = classification, printParams = parameters, printNetwork = network, criterion = attr(object, "Criterion"))
  attr(summ, "Gname")    <- attr(object, "Gname")
  class(summ)   <- "summaryMEDseq"
    summ
}

#' @method summary MEDgating
#' @export
summary.MEDgating  <- function(object, SPS = FALSE, ...) {
  if(length(SPS)   > 1 ||
     !is.logical(SPS))           stop("'SPS' must be a single logical indicator", call.=FALSE)
  equalnoise       <- attr(object, "EqualNoise")
  equalpro         <- attr(object, "EqualPro")
  formula          <- attr(object, "Formula")
  Gname            <- attr(object, "Gname")
  noise            <- attr(object, "Noise")
  noise.gate       <- attr(object, "NoiseGate")
  class(object)    <- class(object)[2L]
  summ             <- summary(object, ...)
  summ$OddsRatios  <- exp(summ$coefficients)
  class(summ)      <- "summaryMEDgate"
  attr(summ, "Class")      <- class(object)
  attr(summ, "EqualNoise") <- equalnoise
  attr(summ, "EqualPro")   <- equalpro
  attr(summ, "Formula")    <- formula
  attr(summ, "Gname")      <- Gname
  attr(summ, "Noise")      <- noise
  attr(summ, "NoiseGate")  <- noise.gate
  attr(summ, "SPS")        <- SPS
    summ
}

#' @method print MEDcriterion
#' @export
print.MEDcriterion       <- function(x, pick = 3L, ...) {
  if(length(pick)        != 1 ||
     !is.numeric(pick))          stop("'pick' must be a single number", call.=FALSE)
  if(floor(pick)  != pick     ||
     pick          < 1)          stop("'pick' be a strictly positive integer", call.=FALSE)
  algo            <- attr(x, "algo")
  crit            <- attr(x, "Criterion")
  choice          <- .pick_MEDCrit(x, pick)
  pick            <- choice$pick
  dim1            <- attr(x, "dim")
  dim2            <- attr(x, "dimnames")
  summ            <- attr(x, "Summ")
  weighted        <- attr(x, "Weighted")
  attributes(x)   <- NULL
  attr(x, "dim")         <- dim1
  attr(x, "dimnames")    <- dim2
  cat(switch(EXPR=crit,
             BIC="Bayesian Information Criterion (BIC):\n",
             ICL="Integrated Completed Likelihood (ICL):\n",
             AIC="Akaike Information Criterion (AIC):\n",
             DBS=paste0(ifelse(weighted, "(Weighted) ", ""), switch(EXPR=summ, median="Median", mean="Mean"), " Density-based Silhouette (DBS):\n"),
             ASW=paste0(ifelse(weighted, "(Weighted) ", ""), switch(EXPR=summ, median="Median", mean="Mean"), " (Average) Silhouette Width (ASW):\n"),
              CV=paste0("Cross-Validated ", ifelse(weighted, "(Weighted) ", ""), "Log-Likelihood (CV):\n"),
             NEC="Normalised Entropy Criterion (NEC):\n",
              DF="Number of Estimated Parameters (Residual DF):\n",
           ITERS=paste0("Number of ", algo, " Iterations:\n"),
          loglik=paste0("Maximal ", ifelse(weighted, "(Weighted) ", ""), "Log-Likelihood:\n")))
  print(unclass(x))
  cat(paste0("\nTop ", ifelse(pick > 1, paste0(pick, " models"), "model"), " based on the ", toupper(crit), " criterion:\n"))
  print(choice$crits)
    invisible()
}

#' @method summary MEDcriterion
#' @export
summary.MEDcriterion <- function(object, G, modtype, ...) {
  attribs            <- attributes(object) 
  if(!missing(G))     {
    object           <- object[rownames(object)  %in% G,,      drop=FALSE]
    attr(object, "dim")        <- dim(object)
    attr(object, "dimnames")   <- list(as.character(G), attribs$dimnames[[2L]])
    attr(object, "G")          <- G
  }
  if(!missing(modtype))         {
    object           <- object[,colnames(object) %in% modtype, drop=FALSE]
    attr(object, "dim")        <- dim(object)
    attr(object, "dimnames")   <- list(attribs$dimnames[[1L]], modtype)
    attr(object, "modelNames") <- modtype
  }
  attr(object, "class")        <- "MEDcriterion"
  attr(object, "Criterion")    <- attribs$Criterion
  object             <- .pick_MEDCrit(object, ...)
  attr(object, "class")        <- "MEDcriterion"
  attr(object, "Criterion")    <- attribs$Criterion
  class(object)      <- "summary.MEDcriterion"
    return(object)
}

#' @method print summary.MEDcriterion
#' @export
print.summary.MEDcriterion <- function(x, digits = 3L, ...){
  if(length(digits) > 1   || !is.numeric(digits) ||
     digits   <= 0)              stop("Invalid 'digits'", call.=FALSE)
  crit        <- attr(x, "Criterion")
  cat(paste0("Best ", crit, " values:\n"))
  x           <- drop(as.matrix(x$crits))
  x           <- rbind(x, x - switch(EXPR=crit, DF=, ITERS=, NEC=min(x), max(x)))
  rownames(x) <- list(crit, paste0(crit, " diff"))
  print(x, digits = digits)
    invisible()
}

#' @method print MEDgating
#' @export
print.MEDgating    <- function(x, call = FALSE, SPS = FALSE, ...) {
  if(length(SPS)   > 1 ||
     !is.logical(SPS))           stop("'SPS' must be a single logical indicator", call.=FALSE)
  equalpro         <- attr(x, "EqualPro")
  formula          <- attr(x, "Formula")
  noise            <- attr(x, "Noise")
  equalNoise       <- noise && equalpro
  gateNoise        <- noise && !equalpro && formula != "~1"
  if(ifelse(inherits(x, "multinom"),
     x$convergence == 1,
     isTRUE(x$converged)))       warning("Multinomial logistic regression failed to converge", call.=FALSE, immediate.=TRUE)
  class(x)         <- class(x)[class(x)  != "MEDgating"]
  if(isTRUE(call)  && 
     !is.null(cl   <- x$call)) {
    cat("Call:\n")
    dput(cl, control = NULL)
  }
  cat("\nCoefficients:\n")
  if(isTRUE(SPS))   {
    x$lab          <- attr(x, "Gname")
  }
  print(stats::coef(x), ...)
  cat(paste("\nFormula:", formula, "\n"))
  cat(paste("Noise:",     noise,   "\n"))
  if(gateNoise)                  cat(paste("Noise Component Gating:", attr(x, "NoiseGate"), "\n"))
  cat(paste("EqualPro:", equalpro, ifelse(equalNoise, "\n", "")))
  if(equalNoise)                 cat(paste("Noise Proportion Estimated:", !attr(x, "EqualNoise")))
  if(equalpro)                   message("\n\nCoefficients set to zero as this is an equal mixing proportion model")
    invisible()
}

#' @method print MEDlambda
#' @export
print.MEDlambda   <- function(x, SPS = FALSE, ...) {
  if(length(SPS)   > 1 ||
     !is.logical(SPS))           stop("'SPS' must be a single logical indicator", call.=FALSE)
  mod             <- attr(x, "Model")
  G               <- attr(x, "G")
  gnames          <- if(isTRUE(SPS)) attr(x, "Gname") else paste0("Cluster", seq_len(G))
  pnames          <- attr(x, "Names")
  attributes(x)[-1L]   <- NULL
  class(x)        <- NULL
  gnames          <- switch(EXPR=mod, CC=,  CU="C",     UC=,  UU=gnames, 
                            CCN=, CUN=c("C", "Noise"), UCN=, UUN=if(isTRUE(SPS)) gnames else replace(gnames, G, "Noise"))
  pnames          <- switch(EXPR=mod, CC=,  UC=, CCN=, UCN="C",  pnames)
  dimnames(x)     <- list(gnames, pnames)
  print(x, ...)
    invisible()
}

#' @method print MEDnames
#' @export
print.MEDnames    <- function(x, ...) {
  attributes(x)   <- 
    class(x)        <- NULL
  print(x, ...)
    invisible()  
}

#' @method print MEDseq
#' @rdname MEDseq_fit
#' @usage
#' \method{print}{MEDseq}(x,
#'       digits = 3L,
#'       ...)
#' @export
print.MEDseq      <- function(x, digits = 3L, ...) {
  cat("Call:\t");  print(x$call)
  if(length(digits)  > 1    || !is.numeric(digits) ||
     digits       <= 0)          stop("Invalid 'digits'", call.=FALSE)
  if(attr(x, "EmptyWarn"))       warning("Model contains at least one empty component\n", call.=FALSE, immediate.=TRUE)
  name            <- x$modtype
  G               <- x$G
  noise           <- attr(x, "Noise")
  gating          <- attr(x$gating, "Formula")
  gate.x          <- !attr(x, "Gating")
  equalP          <- G == 1 || attr(x$gating, "EqualPro")
  equalN          <- noise  && attr(x$gating, "EqualNoise") && equalP
  crit            <- round(c(BIC = x$bic, ICL=x$icl, AIC=x$aic), digits)
  crit2           <- NULL
  crit2           <- if(attr(x, "DBS")) c(crit2, DBS=round(x$dbs, digits)) else crit2
  crit2           <- if(attr(x, "ASW")) c(crit2, ASW=round(x$asw, digits)) else crit2
  crit2           <- if(attr(x, "CV"))  c(crit2, CV =round(x$cv,  digits)) else crit2
  crit2           <- if(attr(x, "NEC")) c(crit2, NEC=round(x$nec, digits)) else crit2
  cat(paste0("\nBest Model", ifelse(length(x$BIC)  > 1, paste0(" (according to ", toupper(attr(x, "Criterion")), "): "), ": "), name, ", with ",
             G, " component",      ifelse(G > 1, "s ", " "),
             paste0("and ",  ifelse(attr(x, "Weighted"), "", "no "), "weights"),
             ifelse(gate.x,        " (no covariates)\n", " (incl. gating covariates)\n"),
             ifelse(!equalP ||
                     G == 1, "",   paste0("Equal Mixing Proportions", ifelse(equalN | G == 1 | !noise, "\n", " (with estimated noise component mixing proportion)\n"))),
             paste(paste0(names(crit),  " = ", crit),  collapse=" | "), "\n",
             paste(paste0(names(crit2), " = ", crit2), collapse=" | "),
             ifelse(gate.x,  "",   paste0("\nGating: ", gating, "\n"))))
    invisible()
}

#' @method print MEDseqCompare
#' @rdname MEDseq_compare
#' @usage
#' \method{print}{MEDseqCompare}(x,
#'       index = seq_len(x$pick),
#'       rerank = FALSE,
#'       digits = 3L,
#'       maxi = length(index),
#'       ...)
#' @export
print.MEDseqCompare    <- function(x, index=seq_len(x$pick), rerank = FALSE, digits = 3L, maxi = length(index), ...) {
  index           <- if(is.logical(index)) which(index) else index
  if(length(index) < 1 || (!is.numeric(index) &&
     (any(index    < 1  | 
          index    > x$pick))))  stop("Invalid 'index'",  call.=FALSE)
  if(length(digits)     > 1 ||
     !is.numeric(digits)    ||
     digits       <= 0)          stop("Invalid 'digits'", call.=FALSE)
  if(length(rerank)     > 1 ||
     !is.logical(rerank))        stop("'rerank' must be a single logical indicator",       call.=FALSE)
  if(length(maxi) != 1 ||
     !is.numeric(maxi) ||
     maxi         <= 0 ||
     floor(maxi)  != maxi)       stop("'maxi' must be a single strictly positive integer", call.=FALSE)
  maxi            <- min(maxi, length(index))
  crit            <- attr(x, "Crit")
  opt             <- attr(x, "Opt")
  x$bic           <- round(x$bic,    digits)
  x$icl           <- round(x$icl,    digits)
  x$aic           <- round(x$aic,    digits)
  x$loglik        <- round(x$loglik, digits)
  na.dbs          <- is.na(x$dbs)
  x$dbs[!na.dbs]  <- if(!all(na.dbs)) round(x$dbs[!na.dbs], digits)
  dbs             <- replace(x$dbs, na.dbs, "")
  x$dbs           <- NULL
  x$dbs           <- if(all(na.dbs))             NULL else dbs
  na.asw          <- is.na(x$asw)
  x$asw[!na.asw]  <- if(!all(na.asw)) round(x$asw[!na.asw], digits)
  asw             <- replace(x$asw, na.asw, "")
  x$asw           <- NULL
  x$asw           <- if(all(na.asw))             NULL else asw
  na.cvs          <- is.na(x$cv)
  x$cv[!na.cvs]   <- if(!all(na.cvs)) round(x$cv[!na.cvs],  digits)
  cvs             <- replace(x$cv, na.cvs,  "")
  x$cv            <- NULL
  x$cv            <- if(all(na.cvs))             NULL else cvs
  na.nec          <- is.na(x$nec)
  x$nec[!na.nec]  <- if(!all(na.nec)) round(x$nec[!na.nec], digits)
  nec             <- replace(x$nec, na.nec, "")
  x$nec           <- NULL
  x$nec           <- if(all(na.nec))             NULL else nec
  x               <- c(x[seq_len(which(names(x) == "loglik") - 1L)], list(dbs=x$dbs, asw=x$asw, cvs=x$cv, nec=x$nec), x[seq(from=which(names(x) == "loglik"), to=length(x) - (2L + !is.null(x$cv) + !is.null(x$nec)), by=1L)])
  x               <- x[unique(names(x))]
  n.all           <- all(!x$noise)
  x$noise         <- if(n.all)                   NULL else x$noise
  noise.gate      <- if(n.all)                   NULL else replace(x$noise.gate, is.na(x$noise.gate), "")
  x$noise.gate    <- NULL
  x$noise.gate    <- if(all(x$gating == "None")) NULL else noise.gate
  equalPro        <- if(all(is.na(x$equalPro)))  NULL else replace(x$equalPro,   is.na(x$equalPro),   "")
  x$equalPro      <- NULL
  x$equalPro      <- equalPro
  na.equalNoise   <- is.na(x$equalNoise)
  equalNoise      <- replace(x$equalNoise, na.equalNoise, "")
  x$equalNoise    <- NULL
  x$equalNoise    <- if(all(na.equalNoise))      NULL else equalNoise
  x$opti          <- if(all(x$opti == "mode"))   NULL else x$opti
  title           <- "Comparison of Mixtures of Exponential-Distance Models with Covariates"
  cat(paste0("---------------------------------------------------------------------\n", 
             title, "\nData: ", x$data, "\nRanking Criterion: ", toupper(crit), "\nOptimal Only: ", opt,
             "\n---------------------------------------------------------------------\n\n"))
  compX           <- data.frame(do.call(cbind, x[-seq_len(3L)]))[index,,           drop=FALSE]
  compX           <- compX[,!vapply(compX, function(x) all(x == ""), logical(1L)), drop=FALSE]
  compX           <- cbind(rank = if(isTRUE(rerank)) seq_along(index) else index, compX)
  rownames(compX) <- NULL
  print(compX[seq_len(maxi),], row.names = FALSE)
    invisible()
}

#' @method print MEDtheta
#' @importFrom TraMineR "seqformat"
#' @export
print.MEDtheta    <- function(x, SPS = FALSE, ...) {
  if(length(SPS)   > 1 ||
     !is.logical(SPS))           stop("'SPS' must be a single logical indicator", call.=FALSE)
  alpha           <- attr(x, "alphabet")
  if(any(attr(x, "NonUnique")))  message("Solution contains at least one non-unique modal central sequence position\n")
  G               <- nrow(x)
  gnames          <- paste0("Cluster", seq_len(G))
  class(x)        <- NULL
  if(is.element(attr(x, "Model"), 
                c("CCN", "UCN", "CUN", "UUN")))    {
    alpha         <- c(alpha, "*")
    gnames        <- replace(gnames, G, "Noise")
    x[G,]         <- length(alpha)
  }
  x               <- as.data.frame(lapply(as.data.frame(x), function(theta) .replace_levels(.num_to_char(theta), alpha)))
  if(isTRUE(SPS))  {
    x             <- provideDimnames(matrix(suppressMessages(seqformat(x, from="STS", to="SPS", compress=TRUE, right=NA, ...)), 
                                            nrow=G, ncol=1L), base=list(gnames, "Sequence"))
    print(x, quote=FALSE, ...)
  } else           {
    rownames(x)   <- gnames
    print(x, ...)
  }
    invisible()
}

#' @method print summaryMEDseq
#' @export
print.summaryMEDseq      <- function(x, digits = 3L, ...) {
  if(length(digits)  > 1 || !is.numeric(digits) ||
     digits     <= 0)            stop("Invalid 'digits'", call.=FALSE)
  tmp           <- data.frame(log.likelihood = round(x$loglik, digits), N = x$N, P = x$P, V = x$V, df = x$df, iters = x$iters)
  SPS           <- x$SPS
  switch(EXPR=toupper(x$criterion),
    BIC={tmp    <- cbind(tmp, BIC=round(x$bic, digits))},
    ICL={tmp    <- cbind(tmp, ICL=round(x$icl, digits))},
    AIC={tmp    <- cbind(tmp, AIC=round(x$aic, digits))},
     CV={tmp    <- cbind(tmp,  CV=round(x$cv,  digits))},
    NEC={tmp    <- cbind(tmp, NEC=round(x$nec, digits))})
  tmp           <- if(is.null(x$dbs)) tmp else cbind(tmp, DBS = round(x$dbs, digits))
  tmp           <- if(is.null(x$asw)) tmp else cbind(tmp, ASW = round(x$asw, digits))
  tmp           <- cbind(tmp, Algo = x$algo)
  rownames(tmp) <- NULL
  name          <- x$modelName
  G             <- x$G
  if(x$emptywarn)                warning("Model contains at least one empty component\n", call.=FALSE, immediate.=TRUE)
  range.G       <- attr(G, "range")
  if(!is.null(range.G)      && length(range.G) > 1 &&
     G          == min(range.G)) message("Best model occurs at the min of the number of components considered\n")
  if(!is.null(range.G)      && length(range.G) > 1 &&
     G          == max(range.G)) message("Best model occurs at the max of the number of components considered\n")
  noise         <- is.element(name, c("CCN", "CUN", "UCN", "UUN"))
  gating        <- attr(x$gating, "Formula")
  gate.x        <- gating   == "~1"
  equalP        <- gate.x   && x$equalPro
  equalN        <- noise    && x$equalNoise && equalP
  title         <- "Mixture of Exponential-Distance Models with Covariates"
  cat(paste0("------------------------------------------------------\n", title, "\nData: ",
             x$data,"\n", "------------------------------------------------------\n\n",
             "MEDseq ", "(", name, "), with ",    paste0(G, " component", ifelse(G > 1, "s", "")),
             paste0("\nGating Network Covariates:  ", ifelse(gate.x, "None", gating)),
             ifelse(G  > 1  && gate.x,            paste0("\nEqual Mixing Proportions:   ", equalP), ""),
             paste0("\nNoise Component:            ", noise, ""),
             ifelse(G  > 1  && !gate.x && noise,  paste0("\nNoise Component Gating:     ", x$noise.gate), ""),
             ifelse(G  > 1  && noise   && equalP, paste0("\nNoise Proportion Estimated: ", !equalN, "\n\n"), "\n\n")))
  print(tmp, row.names = FALSE)
  if(isTRUE(x$printClass))   {
    cat("\nClustering table :")
    CTAB        <- table(x$classification)
    if(isTRUE(SPS))          {
      cat("\n")
      names(CTAB)           <- attr(x, "Gname")[if(noise) c(G, seq_len(G - 1L)) else seq_len(G)]
    }
    print(CTAB, row.names = FALSE)  
  }
  if(isTRUE(x$printParams))  {
    cat("\nMixing proportions :\n")
    if(isTRUE(SPS))          {
      if(is.matrix(x$tau))   {
        colnames(x$tau)     <- attr(x, "Gname")
      } else x$tau          <- stats::setNames(x$tau, attr(x, "Gname"))
    }
    print(x$tau)
    cat("\nComponent central sequences :\n")
    print(x$theta,  SPS=SPS)
    cat("\nComponent precisions :\n")
    print(x$lambda, SPS=SPS)
    cat("\n")
  } else cat("\n")
  if(isTRUE(x$printNetwork)) {
    if(isFALSE(gate.x))      {
      gating    <- list("Gating Network"     = x$gating)
      class(gating)   <- "listof"
      print(gating, SPS=SPS, call = FALSE)
      cat("\n")
    } else                       message("No gating network to display\n")
  }
    invisible()
}

#' @method print summaryMEDgate
#' @export
print.summaryMEDgate  <- function(x, ...) {
  equalpro         <- attr(x, "EqualPro")
  formula          <- attr(x, "Formula")
  noise            <- attr(x, "Noise")
  SPS              <- attr(x, "SPS")
  class(x)         <- "MEDgating"
  if(isTRUE(SPS))   {
    x$lab          <- attr(x, "Gname")
    rownames(x$coefficients)    <-
    rownames(x$OddsRatios)      <-
    rownames(x$standard.errors) <- x$lab[-1L]
  }
  print(x, ...)
  cat("\n\nOddsRatios:\n")
  print(x$OddsRatios)
  cat("\nStd. Errors:\n")
  print(x$standard.errors)
  message("\n\n\nUsers are cautioned against making inferences about statistical significance from summaries of the coefficients in the gating network\nUsers are advised to use the function 'MEDseq_stderr' instead\n")
  class(x)         <- "summaryMEDgate"
    invisible(x)
}

#' MEDseq gating network standard errors
#' 
#' Computes standard errors of the gating network coefficients in a fitted MEDseq model using either the Weighted Likelihood Bootstrap or Jackknife methods.
#' @param mod An object of class \code{"MEDseq"} generated by \code{\link{MEDseq_fit}} or an object of class \code{"MEDseqCompare"} generated by \code{\link{MEDseq_compare}}.
#' @param method The method used to compute the standard errors (defaults to \code{"WLBS"}, the Weighted Likelihood Bootstrap).
#' @param N The (integer) number of samples to use when the \code{"WLBS"} \code{method} is employed. Defaults to \code{1000L}. Not relevant when \code{method="Jackknife"}, in which case \code{N} is always the number of observations. Must be > 1, though \code{N} being greater than or equal to the sample size is recommended under \code{method="WLBS"}.
#' @param symmetric A logical indicating whether symmetric draws from the uniform Dirichlet distribution are used for the \code{WLBS} method in the presence of existing sampling weights. Defaults to \code{TRUE}; when \code{FALSE}, the concentration parameters of the Dirichlet distribution are given by the sampling weights. Only relevant when \code{method="WLBS"} for models with existing sampling weights.
#' @param SPS A logical indicating whether the output should be labelled according to the state-permanence-sequence representation of the central sequences. Defaults to \code{FALSE}. See \code{\link{MEDseq_clustnames}} and \code{\link[TraMineR]{seqformat}}.
#'
#' @return A list with the following two elements:
#' \describe{
#' \item{\code{Coefficients}}{The original matrix of estimated coefficients (\code{coef(mod$gating)}).}
#' \item{\code{Std. Errors}}{The matrix of corresponding standard error estimates.}}
#' @details A progress bar is displayed as the function iterates over the \code{N} samples. The function may take a long time to run for large \code{N}. The function terminates immediately if \code{mod$G == 1}.
#' @note The \code{symmetric} argument is an experimental feature. More generally, caution is advised in interpreting the standard error estimates under \emph{either} the \code{"WLBS"} \strong{or} the \code{"Jackknife"} \code{method} when there are existing sampling weights which arise from complex/stratified sampling designs.
#' @importFrom TraMineR "seqdef" "seqformat"
#' @export
#' @seealso \code{\link{MEDseq_fit}}, \code{\link{MEDseq_clustnames}}, \code{\link[TraMineR]{seqformat}}
#' @references Murphy, K., Murphy, T. B., Piccarreta, R., and Gormley, I. C. (2021). Clustering longitudinal life-course sequences using mixtures of exponential-distance models. \emph{Journal of the Royal Statistical Society: Series A (Statistics in Society)}, 184(4): 1414-1451. <\href{https://rss.onlinelibrary.wiley.com/doi/abs/10.1111/rssa.12712}{doi:10.1111/rssa.12712}>.
#' 
#' O'Hagan, A., Murphy, T. B., Scrucca, L., and Gormley, I. C. (2019). Investigation of parameter uncertainty in clustering using a Gaussian mixture model via jackknife, bootstrap and weighted likelihood bootstrap. \emph{Computational Statistics}, 34(4): 1779-1813.
#' @author Keefe Murphy - <\email{keefe.murphy@@mu.ie}>
#' @keywords clustering main
#' @usage
#' MEDseq_stderr(mod,
#'               method = c("WLBS", "Jackknife"),
#'               N = 1000L,
#'               symmetric = TRUE,
#'               SPS = FALSE)
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
#' # Fit a model with weights and a gating covariate
#' # Have the probability of noise-component membership be constant
#' # mod         <- MEDseq_fit(mvad.seq, G=11, modtype="UUN", weights=mvad$weights, 
#' #                           gating=~ gcse5eq, covars=mvad.cov, noise.gate=FALSE)
#'                             
#' # Estimate standard errors using 100 WLBS samples
#' # (std        <- MEDseq_stderr(mod, N=100))}
MEDseq_stderr <- function(mod, method = c("WLBS", "Jackknife"), N = 1000L, symmetric = TRUE, SPS = FALSE) {
    UseMethod("MEDseq_stderr")
}

#' @method MEDseq_stderr MEDseq
#' @export
MEDseq_stderr.MEDseq <- function(mod, method = c("WLBS", "Jackknife"), N = 1000L, symmetric = TRUE, SPS = FALSE) {
  mod                <- if(inherits(mod, "MEDseqCompare")) mod$optimal else mod
  if(length(SPS)      > 1  ||
     !is.logical(SPS))           stop("'SPS' must be a single logical indicator",   call.=FALSE)
  if(!missing(method)      && 
     (length(method)  > 1  ||
      !is.character(method)))    stop("'method' must be a single character string", call.=FALSE)
  method      <- match.arg(method)
  if(isFALSE(attr(mod, 
                  "Gating")))    warning("No gating covariates in fitted model\n",  call.=FALSE, immediate.=TRUE)
  n           <- nrow(mod$data)
  N           <- switch(EXPR=method, WLBS=N, n)
  if(length(N)       != 1  ||
     !is.numeric(N)        ||
     N        <= 1         ||
     floor(N) != N)              stop("'N' must be a single integer > 1",           call.=FALSE)
  gating      <- mod$gating
  coeffs      <- stats::coef(gating)
  if((mod$G   -> G)  == 1)  {    message("No clustering in fitted model (G=1)\n")
   res        <- list(Coefficients=coeffs, Std.Errors=stats::setNames(0L, names(coeffs)))
   class(res) <- "listof"
    return(res)
  }
  gate.x      <- is.matrix(coeffs)
  coeff       <- array(NA, c(if(gate.x) dim(coeffs) else rep(1L, 2L), N))
  gate        <- stats::as.formula(attr(gating, "Formula"))
  covs        <- mod$covars
  Z           <- switch(EXPR=method, WLBS=list(mod$z), mod$z)
  modtype     <- mod$modtype
  seqdat      <- mod$data
  weights     <- attr(mod, "Weights")
  algo        <- attr(mod, "Algo")
  opti        <- attr(mod, "Opti")
  noise.gate  <- attr(mod$gating, "NoiseGate")
  attr(seqdat, "weights")     <- NULL
  if(method   == "WLBS"       &&
     n         > N)              warning("It is recommended that N match or exceed the sample size\n", call.=FALSE, immediate.=TRUE)
  pb          <- utils::txtProgressBar(min = 0, max = N, style = 3)
  switch(EXPR=method, WLBS= {
    if(length(symmetric)  > 1 ||
       !is.logical(symmetric))   stop("'symmetric' must be a single logical indicator",                call.=FALSE)
    if(isFALSE(symmetric)     &&
       all(weights   == 1))      message("'symmetric=FALSE' has no effect as the fitted model does not use sampling weights\n")
    wts       <- if(isTRUE(symmetric)) replicate(N, weights * .rDirichlet(n)) else replicate(N, .rDirichlet(n, shape=weights))
    for(i in seq_len(N))    {
      STD     <- suppressMessages(MEDseq_fit(seqdat, G=G, modtype=modtype, weights=wts[,i], 
                                             z.list=Z, gating=gate, covars=covs, 
                                             algo=algo, opti=opti, noise.gate=noise.gate, verbose=FALSE))
      coeff[,,i]     <- stats::coef(STD$gating)
      utils::setTxtProgressBar(pb, i)
    }
  },           {
    for(i in seq_len(N))    {
      STD     <- suppressMessages(MEDseq_fit(seqdat[-i,], G=G, modtype=modtype, 
                                             weights=weights[-i], z.list=list(Z[-i,]), 
                                             gating=gate, covars=covs[-i,, drop=FALSE], 
                                             algo=algo, opti=opti, noise.gate=noise.gate, verbose=FALSE))
      coeff[,,i]     <- stats::coef(STD$gating)
      utils::setTxtProgressBar(pb, i)
    }
  })
  message("\n\n")
  std         <- apply(coeff, c(1L, 2L), stats::sd)
  if(isTRUE(SPS))     {
    rownames(coeffs) <- attr(gating, "Gname")[-1L]
  }
  if(gate.x)   {
    dimnames(std)    <- dimnames(coeffs)
  } else       {
    std       <- stats::setNames(as.vector(std), names(coeffs))
  }
  res         <- list(Coefficients=coeffs, Std.Errors=std)
  class(res)  <- "listof"
    res
}

#' Compute the mean time spent in each sequence category
#'
#' Computes the mean time (per cluster) spent in each sequence category (i.e. state value) for a fitted \code{MEDseq} model.
#' @param x An object of class \code{"MEDseq"} generated by \code{\link{MEDseq_fit}} or an object of class \code{"MEDseqCompare"} generated by \code{\link{MEDseq_compare}}.
#' @param MAP A logical indicating whether to use the MAP classification in the computation of the averages, or the 'soft' clustering assignment probabilities given by \code{x$z}. Defaults to \code{FALSE}, but is always \code{TRUE} for models fitted by the CEM algorithm (see \code{\link{MEDseq_control}}). See \code{weighted} for incorporating the sampling weights (regardless of the value of \code{MAP}). See \code{map.size} below.
#' @param weighted A logical indicating whether the sampling weights (if used during model fitting) are used to compute the weighted averages. These can be used alone (when \code{MAP} is \code{TRUE}) or in conjunction with the 'soft' clustering assignment probabilities (when \code{MAP} is \code{FALSE}). Defaults to \code{TRUE}. Note that, \emph{by default}, the first column of the output is not affected by the value of \code{weighted} (see \code{wt.size}).
#' @param norm A logical indicating whether the mean times (outputted values after the first column) are normalised to sum to the sequence length within each cluster (defaults to \code{TRUE}). Otherwise, when \code{FALSE}, entries beyond the first column give the total (weighted) number of times a given sequence category was observed in a given cluster.
#' @param prop A logical (defaulting to \code{FALSE} and only invoked when \code{norm} is also \code{TRUE}) which further normalises the output to give the \emph{proportions} of time spent in each state on average instead of the absolute values.
#' @param map.size A logical (defaulting to \code{FALSE}, unless the model was fitted by the CEM algorithm (see \code{\link{MEDseq_control}})) which overrides \code{MAP} in the \code{Size} column (or \code{Weighted.Size} column, see \code{wt.size}) of the output, e.g. if \code{MAP=FALSE} and \code{map.size=TRUE}, the MAP classification is used to determine the cluster sizes but the soft cluster-membership probabilities are used to calculate quantities in remaining columns. Only relevant when \code{MAP=FALSE} or \code{wt.size=TRUE}.
#' @param wt.size A logical (defaults to \code{FALSE} and only invoked when when \code{weighted} is also \code{TRUE}) which toggles whether the weights are \emph{also} used in the computation of the cluster sizes in the first column of the output (regardless of the values of \code{MAP} or \code{map.size}).
#' @param SPS A logical indicating whether the output should be labelled according to the state-permanence-sequence representation of the central sequences. Defaults to \code{FALSE}. See \code{\link{MEDseq_clustnames}} and \code{\link[TraMineR]{seqformat}}.
#' 
#' @details Models with weights, covariates, &/or a noise component are also accounted for.
#' @note The function \code{\link{plot.MEDseq}} with the option \code{type="mt"} can be used to visualise the mean times (by cluster). However, the results displayed therein (at present) always assume \code{norm=TRUE}, \code{prop=FALSE}, and \code{wt.size=TRUE}, while the \code{MAP} argument is renamed to \code{soft}, where \code{MAP=!soft}.
#' @return A matrix with sequence category and cluster-specific mean times, giving clusters on the rows, corresponding cluster sizes (or weighted cluster sizes) in the first column, and sequence categories in the remaining columns.
#' @importFrom matrixStats "colSums2"
#' @importFrom TraMineR "seqdef" "seqformat" "seqistatd"
#' @export
#' @references Murphy, K., Murphy, T. B., Piccarreta, R., and Gormley, I. C. (2021). Clustering longitudinal life-course sequences using mixtures of exponential-distance models. \emph{Journal of the Royal Statistical Society: Series A (Statistics in Society)}, 184(4): 1414-1451. <\href{https://rss.onlinelibrary.wiley.com/doi/abs/10.1111/rssa.12712}{doi:10.1111/rssa.12712}>.
#' @seealso \code{\link{MEDseq_fit}}, \code{\link{MEDseq_control}}, \code{\link{plot.MEDseq}}
#' @author Keefe Murphy - <\email{keefe.murphy@@mu.ie}>
#' @keywords utility
#' @usage
#' MEDseq_meantime(x,
#'                 MAP = FALSE,
#'                 weighted = TRUE, 
#'                 norm = TRUE,
#'                 prop = FALSE, 
#'                 map.size = FALSE,
#'                 wt.size = FALSE,
#'                 SPS = FALSE)
#' @examples
#' \donttest{data(biofam)
#' seqs <- seqdef(biofam[10:25] + 1L,
#'                states = c("P", "L", "M", "L+M", "C", 
#'                           "L+C", "L+M+C", "D"))
#' mod <- MEDseq_fit(seqs, G=10, modtype="UUN")
#' 
#' MEDseq_meantime(mod)
#' MEDseq_meantime(mod, prop=TRUE)
#' MEDseq_meantime(mod, map.size=TRUE)
#' MEDseq_meantime(mod, MAP=TRUE, norm=FALSE, SPS=TRUE)}
MEDseq_meantime        <- function(x, MAP = FALSE, weighted = TRUE, norm = TRUE, prop = FALSE, 
                                   map.size = FALSE, wt.size = FALSE, SPS = FALSE) {
    UseMethod("MEDseq_meantime")
}

#' @method MEDseq_meantime MEDseq
#' @export
MEDseq_meantime.MEDseq <- function(x, MAP = FALSE, weighted = TRUE, norm = TRUE, prop = FALSE, 
                                   map.size = FALSE, wt.size = FALSE, SPS = FALSE) {
  x               <- if(inherits(x, "MEDseqCompare")) x$optimal else x
  if(length(MAP)   > 1 ||
     !is.logical(MAP))           stop("'MAP' must be a single logical indicator",         call.=FALSE)
  if(length(weighted)   > 1 ||
     !is.logical(weighted))      stop("'weighted' must be a single logical indicator",    call.=FALSE)
  if(length(norm)  > 1 ||
     !is.logical(norm))          stop("'norm' must be a single logical indicator",        call.=FALSE)
  if(length(prop)  > 1 ||
     !is.logical(prop))          stop("'prop' must be a single logical indicator",        call.=FALSE)
  if(length(map.size)   > 1 ||
     !is.logical(map.size))      stop("'map.size' must be a single logical indicator",    call.=FALSE)
  if(length(wt.size)    > 1 ||
     !is.logical(wt.size))       stop("'wt.size' must be a single logical indicator",     call.=FALSE)
  if(length(SPS)   > 1 ||
     !is.logical(SPS))           stop("'SPS' must be a single logical indicator",         call.=FALSE)
  if(attr(x, "EmptyWarn"))       warning("Model contains one or more empty components\n", call.=FALSE, immediate.=TRUE)
  MAP             <- isTRUE(MAP)      && attr(x, "Algo")    != "CEM"
  map.size        <- isTRUE(map.size) || attr(x, "Algo")    == "CEM"
  map.size        <- isTRUE(map.size) && (isTRUE(wt.size)   || isFALSE(MAP)) 
  weighted        <- isTRUE(weighted) && attr(x, "Weighted")
  wt.size         <- isTRUE(weighted) && isTRUE(wt.size)
  weights         <- attr(x, "Weights")
  alph            <- attr(x$data, "alphabet")
  V               <- attr(x, "V")
  P               <- attr(x, "T")
  G               <- x$G
  z               <- x$z
  noise           <- attr(x, "Noise")
  if(isTRUE(SPS))  {
    gnames        <- attr(x, "Gname")
  } else           {
    gnames        <- paste0("Cluster", seq_len(G))
    gnames        <- if(isTRUE(noise))    replace(gnames, G, "Noise") else gnames
  }
  class           <- if(isTRUE(noise))  replace(x$MAP, x$MAP == 0, G) else x$MAP 
  tabMAP          <- if(isTRUE(MAP))         tabulate(class, nbins=G) else colSums2(z)
  z               <- if(isTRUE(weighted))                 z * weights else z
  tabMAPw         <- if(any(weighted, wt.size))           colSums2(z) else tabMAP
  if(isTRUE(wt.size))     {
    wtabMAP       <- if(isTRUE(map.size)) tapply(weights, class, sum) else tabMAPw
  } else wtabMAP  <- if(map.size && 
                        isFALSE(MAP))        tabulate(class, nbins=G) else tabMAP
  if(isTRUE(MAP))  {
    if(isTRUE(weighted))  {
      temp        <- suppressMessages(seqistatd(x$data)) * weights
      temp        <- do.call(rbind, lapply(seq_len(G), function(g) colSums2(temp[class == g,, drop=FALSE])))
    } else temp   <- do.call(rbind, by(x$data, class,  function(x) tabulate(do.call(base::c, x), V)))
    if(nrow(temp) != G)   {
      temp        <- rbind(temp, matrix(0L, nrow=sum(tabMAP == 0), ncol=V))
    }
  }   else         {
    x$data        <- do.call(base::c, .fac_to_num(x$data))
    temp          <- do.call(rbind, lapply(seq_len(G), function(g) tapply(rep(z[,g], P), droplevels(x$data), sum)))  
  }
  if(isTRUE(weighted))    {
    temp          <- temp / tabMAPw
    temp          <- if(isFALSE(norm)) temp * tabMAP                  else if(isTRUE(prop)) temp / P else temp
  } else if(isTRUE(norm)) {
    temp          <- if(isTRUE(prop))  temp / (tabMAP    * P)         else temp / tabMAP
  }
    res           <- provideDimnames(unname(cbind(wtabMAP, temp)), base=list(gnames, c(ifelse(isTRUE(wt.size), "Weighted.Size", "Size"), alph)))
    class(res)    <- "MEDseqMeanTime"
      return(res)
}

#' @method print MEDseqMeanTime
#' @param digits Minimum number of significant digits to be printed in values.
#' @param ... Catches unused arguments.
#' @rdname MEDseq_meantime
#' @usage
#' \method{print}{MEDseqMeanTime}(x,
#'       digits = 3L,
#'       ...)
#' @export
print.MEDseqMeanTime <- function(x, digits = 3L, ...) {
  if(length(digits)   > 1 || 
     !is.numeric(digits)  ||
     digits          <= 0)       stop("Invalid 'digits'", call.=FALSE)
  print(unclass(x), digits = digits)
    invisible()
}

#' Automatic labelling of clusters using central sequences
#' 
#' These functions extract names for clusters according to the SPS representation of their central sequences.
#' @param x An object of class \code{"MEDseq"} generated by \code{\link{MEDseq_fit}} or an object of class \code{"MEDseqCompare"} generated by \code{\link{MEDseq_compare}}.
#' @param cluster A logical indicating whether names should be prepended with the text "\code{Cluster g: }", where \code{g} is the cluster number. Defaults to \code{TRUE}. 
#' @param size A logical indicating whether the (typically 'soft') size of each cluster is appended to the label of each group, expressed as a percentage of the total number of observations. Defaults to \code{FALSE}. 
#' @param weighted A logical indicating whether the sampling weights (if any) are used when appending the \code{size} of each cluster to the labels. Defaults to \code{FALSE}.
#' @param ... Catches unused arguments. 
#' @param names The output of \code{MEDseq_clustnames} to be passed to the convenience function \code{MEDseq_nameclusts} (see \code{Details}).
#' 
#' @details Unlike the \code{\link[WeightedCluster]{seqclustname}} function from the \pkg{WeightedCluster} package which inspired these functions, \code{MEDseq_clustnames} only returns the names themselves, not the \code{factor} variable indicating cluster membership with labels given by those names. Thus, \code{MEDseq_nameclusts} is provided as a convenience function for precisely this purpose (see \code{Examples}). 
#' @return For \code{MEDseq_clustnames}, a character vector containing the names for each component defined by their central sequence, and optionally the cluster name (see \code{cluster} above) and cluster size (see \code{size} above). The name for the noise component, if any, will always be simply \code{"Noise"} (or \code{"Cluster 0: Noise"}).
#' 
#' For \code{MEDseq_nameclusts}, a factor version of \code{x$MAP} with levels given by the output of \code{MEDseq_clustnames}.
#' @note The main \code{MEDseq_clustnames} function is used internally by \code{\link{plot.MEDseq}}, \code{\link{MEDseq_meantime}}, \code{\link{MEDseq_stderr}}, and also other \code{print} and \code{summary} methods, where its invocation can typically controlled via a \code{SPS} logical argument. However, the optional arguments \code{cluster}, \code{size}, and \code{weighted} can only be passed through \code{\link{plot.MEDseq}}; elsewhere \code{cluster=TRUE}, \code{size=FALSE}, and \code{weighted=FALSE} are always assumed.
#' 
#' @references Murphy, K., Murphy, T. B., Piccarreta, R., and Gormley, I. C. (2021). Clustering longitudinal life-course sequences using mixtures of exponential-distance models. \emph{Journal of the Royal Statistical Society: Series A (Statistics in Society)}, 184(4): 1414-1451. <\href{https://rss.onlinelibrary.wiley.com/doi/abs/10.1111/rssa.12712}{doi:10.1111/rssa.12712}>.
#' @author Keefe Murphy - <\email{keefe.murphy@@mu.ie}>
#' @seealso \code{\link[TraMineR]{seqformat}}, \code{\link[WeightedCluster]{seqclustname}}, \code{\link{plot.MEDseq}}, \code{\link{MEDseq_meantime}}, \code{\link{MEDseq_stderr}}
#' @keywords utility
#' @importFrom matrixStats "colMeans2"
#' @importFrom TraMineR "seqformat"
#' @importFrom WeightedCluster "seqclustname"
#' @export
#' @usage
#' MEDseq_clustnames(x,
#'                   cluster = TRUE,
#'                   size = FALSE,
#'                   weighted = FALSE,
#'                   ...)
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
#' # Fit a model with weights and a gating covariate
#' # Have the probability of noise-component membership depend on the covariate
#' mod    <- MEDseq_fit(mvad.seq, G=5, modtype="UUN", weights=mvad$weights, 
#'                      gating=~ gcse5eq, covars=mvad.cov, noise.gate=TRUE)
#'                      
#' # Extract the names
#' names  <- MEDseq_clustnames(mod, cluster=FALSE, size=TRUE)
#' 
#' # Get the renamed MAP cluster membership indicator vector
#' group  <- MEDseq_nameclusts(names)
#' 
#' # Use the output in plots
#' plot(mod, type="d", soft=FALSE, weighted=FALSE, cluster=FALSE, size=TRUE, border=TRUE)
#' # same as:
#' # seqplot(mvad.seq, type="d", group=group)
#' 
#' # Indeed, this function is invoked by default for certain plot types
#' plot(mod, type="d", soft=TRUE, weighted=TRUE)
#' plot(mod, type="d", soft=TRUE, weighted=TRUE, SPS=FALSE)
#' 
#' # Invoke this function when printing the gating network coefficients
#' print(mod$gating, SPS=FALSE)
#' print(mod$gating, SPS=TRUE)
#' 
#' # Invoke this function in a call to MEDseq_meantime
#' MEDseq_meantime(mod, SPS=TRUE)
#'  
#' # Invoke this function in other plots
#' plot(mod, type="clusters", SPS=TRUE)
#' plot(mod, type="precision", SPS=TRUE)}
MEDseq_clustnames        <- function(x, cluster = TRUE, size = FALSE, weighted = FALSE, ...) {
    UseMethod("MEDseq_clustnames")
}

#' @method MEDseq_clustnames MEDseq
#' @export
MEDseq_clustnames.MEDseq <- function(x, cluster = TRUE, size = FALSE, weighted = FALSE, ...) {
  x               <- if(inherits(x, "MEDseqCompare")) x$optimal else x
  if(length(cluster)     != 1 ||
     !is.logical(cluster))       stop("'cluster' must be a single logical indicator", call.=FALSE)
  if(length(size) != 1   ||
     !is.logical(size))          stop("'size' must be a single logical indicator",    call.=FALSE)
  X               <- utils::capture.output(print(x$params$theta, SPS=TRUE))[-1L]
  X               <- paste0("(", sub("^.*?\\((.*)\\)[^)]*$", "\\1", X), ")")
  noise           <- vapply(X, function(x) grepl("*", x, fixed = TRUE), logical(1L))
  X               <- if(any(noise))      c(X[!noise], "Noise")                                                         else X
  X               <- if(isTRUE(cluster)) paste0("Cluster ", replace(seq_len(x$G), noise, 0L), ": ", X)                 else X
  X               <- if(isTRUE(size))    paste0(X, "~[", round(colMeans2(x$z * if(isTRUE(weighted)) attr(x, "Weights") else 1L) * 100L, 1L), "%]") else X
  attr(X, "MAP")  <- x$MAP
  attr(X, "G")    <- x$G
  class(X)        <- "MEDnames"
    return(X)
}

#' @rdname MEDseq_clustnames
#' @keywords utility
#' @usage MEDseq_nameclusts(names)
#' @export
MEDseq_nameclusts <- function(names) {
    UseMethod("MEDseq_nameclusts")
}

#' @method MEDseq_nameclusts MEDnames
#' @export
MEDseq_nameclusts.MEDnames  <- function(names) {
  MAP             <- attr(names, "MAP")
    factor(replace(MAP, MAP == 0, attr(names, "G")), labels=names)
}

#' Entropy of a fitted MEDseq model
#'
#' Calculates the normalised entropy of a fitted MEDseq model.
#' @param x An object of class \code{"MEDseq"} generated by \code{\link{MEDseq_fit}} or an object of class \code{"MEDseqCompare"} generated by \code{\link{MEDseq_compare}}.
#'
#' @details This function calculates the normalised entropy via \deqn{H=-\frac{1}{n\log(G)}\sum_{i=1}^n\sum_{g=1}^G\hat{z}_{ig}\log(\hat{z}_{ig}),}
#' where \eqn{n} and \eqn{G} are the sample size and number of components, respectively, and \eqn{\hat{z}_{ig}} is the estimated posterior probability at convergence that observation \eqn{i} belongs to component \eqn{g}.
#' @return A single number, given by \eqn{1-H}, in the range [0,1], such that \emph{larger} values indicate clearer separation of the clusters.
#' @note This function will always return a normalised entropy of \code{1} for models fitted using the \code{"CEM"} algorithm (see \code{\link{MEDseq_control}}), or models with only one component.
#' @seealso \code{\link{MEDseq_fit}}, \code{\link{MEDseq_control}}, \code{\link{MEDseq_AvePP}}
#' @references Murphy, K., Murphy, T. B., Piccarreta, R., and Gormley, I. C. (2021). Clustering longitudinal life-course sequences using mixtures of exponential-distance models. \emph{Journal of the Royal Statistical Society: Series A (Statistics in Society)}, 184(4): 1414-1451. <\href{https://rss.onlinelibrary.wiley.com/doi/abs/10.1111/rssa.12712}{doi:10.1111/rssa.12712}>.
#' @author Keefe Murphy - <\email{keefe.murphy@@mu.ie}>
#' @keywords utility
#' @usage
#' MEDseq_entropy(x)
#' @export
#'
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
#' 
#' # Fit a model with weights and a gating covariate
#' # Have the probability of noise-component membership be constant
#' mod           <- MEDseq_fit(mvad.seq, G=11, modtype="UUN", weights=mvad$weights, 
#'                             gating=~ gcse5eq, covars=mvad.cov, noise.gate=FALSE)
#' 
#' # Calculate the normalised entropy
#' MEDseq_entropy(mod)                             
MEDseq_entropy <- function(x) {
    UseMethod("MEDseq_entropy")
}

#' @method MEDseq_entropy MEDseq
#' @export
MEDseq_entropy.MEDseq      <- function(x) {
  x           <- if(inherits(x, "MEDseqCompare")) x$optimal else x
  z           <- x$z
  G           <- ncol(z)
  n           <- nrow(z)
    ifelse(attr(x, "Algo") == "CEM" || G == 1, 1L, pmax(0L, 1 - .entropy(z)/(n * log(G))))
}

#' Average posterior probabilities of a fitted MEDseq model
#'
#' Calculates the per-component average posterior probabilities of a fitted MEDseq model.
#' @param x An object of class \code{"MEDseq"} generated by \code{\link{MEDseq_fit}} or an object of class \code{"MEDseqCompare"} generated by \code{\link{MEDseq_compare}}.
#'
#' @details This function calculates AvePP, the average posterior probability of membership for each component for the observations assigned to that component via MAP probabilities.
#' @return A named vector of numbers, of length equal to the number of components (G), in the range [1/G,1], such that \emph{larger} values indicate clearer separation of the clusters.
#' @note This function will always return values of \code{1} for all components for models fitted using the \code{"CEM"} algorithm (see \code{\link{MEDseq_control}}), or models with only one component.
#' @seealso \code{\link{MEDseq_fit}}, \code{\link{MEDseq_control}}, \code{\link{MEDseq_entropy}}
#' @references Murphy, K., Murphy, T. B., Piccarreta, R., and Gormley, I. C. (2021). Clustering longitudinal life-course sequences using mixtures of exponential-distance models. \emph{Journal of the Royal Statistical Society: Series A (Statistics in Society)}, 184(4): 1414-1451. <\href{https://rss.onlinelibrary.wiley.com/doi/abs/10.1111/rssa.12712}{doi:10.1111/rssa.12712}>.
#' @author Keefe Murphy - <\email{keefe.murphy@@mu.ie}>
#' @keywords utility
#' @usage
#' MEDseq_AvePP(x)
#' @export
#'
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
#' 
#' # Fit a model with weights and a gating covariate
#' # Have the probability of noise-component membership be constant
#' mod           <- MEDseq_fit(mvad.seq, G=11, modtype="UUN", weights=mvad$weights, 
#'                             gating=~ gcse5eq, covars=mvad.cov, noise.gate=FALSE)
#' 
#' # Calculate the AvePP
#' MEDseq_AvePP(mod)
MEDseq_AvePP               <- function(x) {
    UseMethod("MEDseq_AvePP")
}

#' @method MEDseq_AvePP MEDseq
#' @export
MEDseq_AvePP.MEDseq        <- function(x) {
  x           <- if(inherits(x, "MEDseqCompare")) x$optimal else x
  z    <- x$z
  G    <- ncol(z)
  nX   <- attr(x, "Noise")
  gnam <- paste0("Group", seq_len(G))
  MAP  <- replace(x$MAP, x$MAP == 0, G)
    stats::setNames(vapply(seq_len(G), function(g) mean(z[MAP == g,g]), numeric(1L)),
                    if(isTRUE(nX)) replace(gnam, G, "Noise") else gnam)
}
#' Weighted K-Modes Clustering with Tie-Breaking
#'
#' Perform k-modes clustering on categorical data with observation-specific sampling weights and tie-breaking adjustments.
#' @param data A matrix or data frame of categorical data. Objects have to be in rows, variables in columns.
#' @param modes Either the number of modes or a set of initial (distinct) cluster modes (where each mode is a row and \code{modes} has the same number of columns as \code{data}). If a number, a random set of (distinct) rows in \code{data} is chosen as the initial modes. Note, this randomness is always present, and is not governed by \code{random} below.
#' @param weights Optional numeric vector containing non-negative observation-specific case weights.
#' @param iter.max The maximum number of iterations allowed. Defaults to \code{.Machine$integer.max}. The algorithm terminates when \code{iter.max} is reached or when the partition ceases to change between iterations.
#' @param freq.weighted A logical indicating whether the usual simple-matching (Hamming) distance between objects is used, or a frequency weighted version of this distance. Defaults to \code{FALSE}; when \code{TRUE}, the frequency weights are computed within the algorithm and are \emph{not} user-specified. Distinct from the observation-level \code{weights} above, the frequency weights are assigned on a per-feature basis and derived from the categories represented in each column of \code{data}. For convenience, the function \code{dist_freqwH} is provided for calculating the corresponding pairwise dissimilarity matrix for subsequent use.
#' @param fast A logical indicating whether a fast version of the algorithm should be applied. Defaults to \code{TRUE}.
#' @param random A logical indicating whether ties for the modal values &/or assignments are broken at random. Defaults to \code{TRUE}: the implied default had been \code{FALSE} prior to version \code{1.3.2} of this package, as per \code{klaR::kmodes} prior to version \code{1.7-1} (see Note). Note that when \code{modes} is specified as the number of modes, the algorithm is \emph{always} randomly initialised, regardless of the specification of \code{random}.
#' 
#' Regarding the modes, ties are broken at random when \code{TRUE} and the first candidate state is always chosen for the mode when \code{FALSE}. Regarding assignments, tie-breaking is always first biased in favour of the observation's most recent cluster: regarding ties thereafter, these are broken at random when \code{TRUE} or the first other candidate cluster is always chosen when \code{FALSE}.
#' @param ... Catches unused arguments.
#'
#' @details The k-modes algorithm (Huang, 1997) is an extension of the k-means algorithm by MacQueen (1967).
#' 
#' The data given by \code{data} is clustered by the k-modes method (Huang, 1997) which aims to partition the objects into k groups such that the distance from objects to the assigned cluster modes is minimised. 
#' 
#' By default, the simple-matching (Hamming) distance is used to determine the dissimilarity of two objects. It is computed by counting the number of mismatches in all variables. Alternatively, this distance can be weighted by the frequencies of the categories in data, using the \code{freq.weighted} argument (see Huang, 1997, for details). For convenience, the function \code{dist_freqwH} is provided for calculating the corresponding pairwise dissimilarity matrix for subsequent use.
#' 
#' If an initial matrix of modes is supplied, it is possible that no object will be closest to one or more modes. In this case, fewer clusters than the number of supplied modes will be returned and a warning will be printed.
#' 
#' If called using \code{fast = TRUE}, the reassignment of the data to clusters is done for the entire data set before recomputation of the modes is done. For computational reasons, this option should be chosen for all but the most moderate of data sizes.
#' 
#' @note This code is adapted from the \code{kmodes} function in the \pkg{klaR} package. Specifically, modifications were made to allow for random tie-breaking for the modes and assignments (see \code{random} above) and the incorporation of observation-specific sampling \code{weights}, with a view to using this function as a means to initialise the allocations for MEDseq models (see the \code{\link{MEDseq_control}} argument \code{init.z} and the related options \code{"kmodes"} and \code{"kmodes2"}). 
#' 
#' Notably, the \code{wKModes} function, when invoked inside \code{\link{MEDseq_fit}}, is used regardless of whether the weights are true sampling weights, or the weights are merely aggregation weights, or there are no weights at all. Furthermore, the \code{\link{MEDseq_control}} argument \code{random} is \emph{also} passed to \code{wKModes} when it is invoked inside \code{\link{MEDseq_fit}}.
#' 
#' \strong{Update}: as of version \code{1.7-1} of \pkg{klaR}, \code{klaR::kmodes} now breaks assignment ties at random only when \code{fast=TRUE}. It still breaks assignment ties when \code{fast=FALSE} and all ties for modal values in the non-random manner described above. Thus, the old behaviour of \code{klaR::kmodes} can be recovered by specifying \code{random=FALSE} here, but \code{random=TRUE} allows random tie-breaking for both types of ties in all situations.
#' @return An object of class \code{"wKModes"} which is a list with the following components:
#' \describe{
#' \item{\code{cluster}}{A vector of integers indicating the cluster to which each object is allocated.}
#' \item{\code{size}}{The number of objects in each cluster.}
#' \item{\code{modes}}{A matrix of cluster modes.}
#' \item{\code{withindiff}}{The within-cluster (weighted) simple-matching distance for each cluster.}
#' \item{\code{tot.withindiff}}{The total within-cluster (weighted) distance over all clusters. \code{tot.withindiff} can be used to guide the choice of the number of clusters, but beware of inherent randomness in the algorithm, which is liable to yield a jagged elbow plot (see examples).}
#' \item{\code{iterations}}{The number of iterations the algorithm reached.}
#' \item{\code{weighted}}{A logical indicating whether observation-level \code{weights} were used or not throughout the algorithm.}
#' \item{\code{freq.weighted}}{A logical indicating whether feature-level \code{freq.weights} were used or not in the computation of the distances. For convenience, the function \code{dist_freqwH} is provided for calculating the corresponding pairwise dissimilarity matrix for subsequent use.}
#' \item{\code{random}}{A logical indicating whether ties were broken at random or not throughout the algorithm.}}
#' @references Huang, Z. (1997). A fast clustering algorithm to cluster very large categorical data sets in data mining. In H. Lu, H. Motoda, and H. Luu (Eds.), \emph{KDD: Techniques and Applications}, pp. 21-34. Singapore: World Scientific.
#' 
#' MacQueen, J. (1967). Some methods for classification and analysis of multivariate observations. In L. M. L. Cam and J. Neyman (Eds.), \emph{Proceedings of the Fifth Berkeley Symposium on  Mathematical Statistics and Probability}, Volume 1, pp. 281-297. Berkeley, CA, USA: University of California Press.
#' @author Keefe Murphy - <\email{keefe.murphy@@mu.ie}>
#' (adapted from \code{klaR::kmodes})
#' @seealso \code{\link{MEDseq_control}}, \code{\link{MEDseq_fit}}, \code{\link{dist_freqwH}}
#' @keywords utility
#' @importFrom matrixStats "colMeans2"
#' @importFrom TraMineR "seqdef" "seqformat"
#' @importFrom WeightedCluster "wcAggregateCases"
#' @export
#' @usage 
#' wKModes(data,
#'         modes,
#'         weights = NULL,
#'         iter.max = .Machine$integer.max,
#'         freq.weighted = FALSE,
#'         fast = TRUE,
#'         random = TRUE,
#'         ...)
#' @examples
#' suppressMessages(require(WeightedCluster))
#' set.seed(99)
#' # Load the MVAD data & aggregate the state sequences
#' data(mvad)
#' agg      <- wcAggregateCases(mvad[,17:86], weights=mvad$weight)
#' 
#' # Create a state sequence object without the first two (summer) time points
#' states   <- c("EM", "FE", "HE", "JL", "SC", "TR")
#' labels   <- c("Employment", "Further Education", "Higher Education", 
#'               "Joblessness", "School", "Training")
#' mvad.seq <- seqdef(mvad[agg$aggIndex, 17:86], 
#'                    states=states, labels=labels, 
#'                    weights=agg$aggWeights)
#' 
#' # Run k-modes without the weights
#' resX     <- wKModes(mvad.seq, 2)
#' 
#' # Run k-modes with the weights
#' resW     <- wKModes(mvad.seq, 2, weights=agg$aggWeights)
#' 
#' # Examine the modal sequences of both solutions
#' seqformat(seqdef(resX$modes), from="STS", to="SPS", compress=TRUE)
#' seqformat(seqdef(resW$modes), from="STS", to="SPS", compress=TRUE)
#' 
#' # Using tot.withindiff to choose the number of clusters
#' \donttest{
#' TWdiffs   <- sapply(1:5, function(k) wKModes(mvad.seq, k, weights=agg$aggWeights)$tot.withindiff)
#' plot(TWdiffs, type="b", xlab="K")
#' 
#' # Use multiple random starts to account for inherent randomness
#' TWDiff    <- sapply(1:5, function(k) min(replicate(10, 
#'                     wKModes(mvad.seq, k, weights=agg$aggWeights)$tot.withindiff)))
#' plot(TWDiff, type="b", xlab="K")}
wKModes       <- function(data, modes, weights = NULL, iter.max = .Machine$integer.max, 
                          freq.weighted = FALSE, fast = TRUE, random = TRUE, ...) {
  data        <- as.data.frame(data)
  n           <- nrow(data)
  P           <- ncol(data)
  isnumeric   <- vapply(data, is.numeric, logical(1L))
  isfactor    <- vapply(data, is.factor,  logical(1L))
  if(any(isfactor))            {
    levs      <- vector("list", P)
    for(j in which(isfactor))  {
      levsj   <- levels(data[,j])
      data[,j]    <- levsj[data[,j]]
      levs[[j]]   <- levsj
    }
  }
  if(any(isnumeric))           {
    lengths   <- vapply(data[,isnumeric], function(z) length(unique(z)), numeric(1L))
    if(any(lengths > 30))        warning("data has numeric columns with more than 30 different levels!", call.=FALSE, immediate.=TRUE)
  }
  if(wtd      <- !is.null(weights))  {
    if(!is.numeric(weights) ||
       length(weights) != n)     stop(paste0("'weights' must be a numeric vector of length N=", n),      call.=FALSE)
    if(any(weights < 0)     || 
       any(!is.finite(weights))) stop("'weights' must be non-negative and finite",           call.=FALSE)
  }
  if(!is.numeric(iter.max)  ||
     floor(iter.max)   != 
     iter.max ||
     iter.max <= 0)              stop("'itmax' must contain strictly positive integers",     call.=FALSE)
  iter.max    <- pmin(iter.max, .Machine$integer.max)
  if(length(freq.weighted)   > 1    ||
     !is.logical(freq.weighted)) stop("'freq.weighted' must be a single logical indicator",  call.=FALSE)
  if(length(fast)       > 1 ||
     !is.logical(fast))          stop("'fast' must be a single logical indicator",           call.=FALSE)
  if(length(random)     > 1 ||
     !is.logical(random))        stop("'random' must be a single logical indicator",         call.=FALSE)
  nseq        <- seq_len(n)
  Pseq        <- seq_len(P)
  cluster     <- numeric(n)
  names(cluster)  <- nseq
  if(missing(modes))             stop("'modes' must be a number or a data frame",            call.=FALSE)
  if(iter.max  < 1)              stop("'iter.max' must be positive",                         call.=FALSE)
  if(length(modes)     == 1) {
    k         <- modes
    kseq      <- seq_len(k)
    modes     <- unique(data)[sample(nrow(unique(data)), size=k, replace=FALSE),,             drop=FALSE]
    kseq      -> cluster[which(rownames(data) %in% rownames(modes))]
  } else       {
    if(any(duplicated(modes)))   stop("Initial modes are not distinct",                      call.=FALSE)
    if(P      != ncol(modes))    stop("'data' and 'modes' must have same number of columns", call.=FALSE)
    modes     <- as.data.frame(modes)
    if(any(isfactor)) {
      if(!all(vapply(modes[,isfactor], is.factor, 
              logical(1L))))     stop("Types of modes do not match data!",                   call.=FALSE)
      for(j in which(isfactor))      {
        modes[,j] <- levels(modes[,j])[modes[,j]]
      }
    }
    if(freq.weighted) {
      for(j in Pseq)  {
        if(!all(modes[,j] %in% 
            unique(data[,j])))   stop("Values of modes must exist in data when 'freq.weighted' is TRUE", call.=FALSE)
      }
    }
    k         <- nrow(modes)
    kseq      <- seq_len(k)
  }
  if(k > nrow(unique(data)))     stop("More cluster modes than distinct data points",        call.=FALSE)
  frwts       <- if(freq.weighted) lapply(Pseq, function(i) table(data[,i]))
  if(!fast)    {
    for(j in which(cluster  == 0))   {
      dist    <- apply(modes, 1L, .km_dist, data[j,], frwts)
      if(!random)  {
       cluster[j] <- which.min(dist)
      } else       {
       cluster[j] <- .rand_MIN(dist)
      }
      modes[cluster[j],]    <- .update_mode(cluster[j], cluster, data, random, weights)
    }
    for(i in seq_len(iter.max))      {
      continue    <- FALSE
      for(j in nseq)  {
        dist      <- apply(modes, 1L, .km_dist, data[j,], frwts)
        clust_new <- .random_ass(dist, cluster[j], min, random, na.rm=TRUE)
        clust_old <- cluster[j]
        if(clust_new        != clust_old) {
          cluster[j]        <- clust_new
          modes[clust_new,] <- .update_mode(clust_new, cluster, data, random, weights)
          modes[clust_old,] <- .update_mode(clust_old, cluster, data, random, weights)
          continue          <- TRUE
        }
      }
      if(!continue)              break
    }
  }     else   {
    dists     <- matrix(NA, nrow = n, ncol = k)
    if(!freq.weighted)       {
      for(i in kseq)         {
        di    <- vapply(Pseq, function(j) return(data[,j] != rep.int(modes[i,j], n)), logical(n))
        dists[,i] <- if(wtd) rowSums2(di) * weights        else rowSums2(di)
      }
    }   else   {
      n_obj   <- matrix(NA, nrow = n, ncol = P)
      n_mode  <- matrix(NA, nrow = nrow(modes), ncol = P)
      for(j in Pseq)      n_obj[,j] <- frwts[[j]][vapply(as.character(data[,j]),  function(z) return(which(names(frwts[[j]]) == z)), numeric(1L))]
      for(j in Pseq)     n_mode[,j] <- frwts[[j]][vapply(as.character(modes[,j]), function(z) return(which(names(frwts[[j]]) == z)), numeric(1L))]
      for(i in kseq)   {
        di    <- vapply(Pseq, function(j) return(data[,j] != rep.int(modes[i,j], n)), logical(n))
        wts   <- 1/n_mode[rep.int(i, n),] + 1/n_obj
        dists[,i] <- if(wtd) rowSums2(di  * wts) * weights else rowSums2(di * wts)
      }
    }
    cluster   <- apply(dists, 1L, .rand_MIN, random)
    for(j in seq_len(nrow(modes)))   {
      modes[j,]   <- .update_mode(j, cluster, data, random, weights)
    }
    for(i in seq_len(iter.max))      {
      continue    <- FALSE
      dists   <- matrix(NA, nrow = n, ncol = k)
      if(!freq.weighted)     {
       for(i in kseq)        {
        di    <- vapply(Pseq, function(j) return(data[,j] != rep.int(modes[i,j], n)), logical(n))
        dists[,i] <- if(wtd) rowSums2(di) * weights        else rowSums2(di)
       }
      } else   {
       n_mode <- matrix(NA, nrow = nrow(modes), ncol = P)
       for(j in Pseq)    n_mode[,j] <- frwts[[j]][vapply(as.character(modes[,j]), function(z) return(which(names(frwts[[j]]) == z)), numeric(1L))]
       for(i in kseq) {
        di    <- vapply(Pseq, function(j) return(data[,j] != rep.int(modes[i,j], n)), logical(n))
        wts   <- 1/n_mode[rep.int(i, n),] + 1/n_obj
        dists[,i] <- if(wtd) rowSums2(di  * wts) * weights else rowSums2(di * wts)
       }
      }
      old.cluster <- cluster
      cluster <- vapply(nseq, function(i) .random_ass(dists[i,], cluster[i], fun=min, random, na.rm=TRUE), integer(1L))
      for(j in seq_len(nrow(modes))) {
        modes[j,] <- .update_mode(j, cluster, data, random, weights)
      }
      continue    <- any(cluster    != old.cluster)
      if(!continue)              break
    }
  }
  cluster.size    <- table(cluster)
  if(length(cluster.size) < k)   warning("One or more clusters are empty",                               call.=FALSE, immediate.=TRUE)
  dists       <- matrix(NA, nrow = n, ncol = k)
  if(freq.weighted)   {
    n_mode    <- matrix(NA, nrow = nrow(modes), ncol = P)
    for(j in Pseq)       n_mode[,j] <- frwts[[j]][vapply(as.character(modes[,j]), function(z) return(which(names(frwts[[j]]) == z)), numeric(1L))]
  }
  for(i in kseq)      {
    di        <- vapply(Pseq, function(j) return(data[,j] != rep.int(modes[i,j], n)), logical(n))
    if(freq.weighted) {
      if(!fast)    {
        n_obj     <- matrix(NA, nrow = n, ncol = P)
        for(j in Pseq)    n_obj[,j] <- frwts[[j]][vapply(as.character(data[,j]),  function(z) return(which(names(frwts[[j]]) == z)), numeric(1L))]
      }
      wts     <- 1/n_mode[rep.int(i, n),] + 1/n_obj
      di      <- rowSums2(di * wts)
    } else di <- rowSums2(di)
    dists[,i] <- if(wtd)  di * weights else di
  }
  diffs       <- numeric(k)
  for(i in seq_along(cluster.size))  {
    diffs[i]  <- sum(dists[cluster  == i,i])
  }   
  dimnames(modes) <- list(kseq, colnames(data))
  if(any(isfactor))  for(j in which(isfactor))  modes[,j] <- factor(modes[,j], levels = levs[[j]])
  if(any(isnumeric)) for(j in which(isnumeric)) modes[,j] <- as.numeric(modes[,j])
  result      <- list(cluster = cluster, size = cluster.size, modes = modes, 
                      withindiff = diffs, tot.withindiff = sum(diffs), iterations = i, 
                      weighted = wtd, freq.weighted = freq.weighted, random = random)
  class(result)   <- "wKModes"
    return(result)
}

#' Pairwise frequency-Weighted Hamming distance matrix for categorical data
#'
#' Computes the matrix of pairwise distance using a frequency-weighted variant of the Hamming distance often used in k-modes clustering.
#' @param data A matrix or data frame of categorical data. Objects have to be in rows, variables in columns.
#' @param full.matrix Logical. If \code{TRUE} (the default), the full pairwise distance matrix is returned, otherwise an object of class \code{\link[stats]{dist}} is returned, i.e. a vector containing only values from the upper triangle of the distance matrix. Objects of class \code{dist} are smaller and can be passed directly as arguments to most clustering functions.
#'
#' @details As per \code{\link{wKModes}}, the frequency weights are computed within the function and are \emph{not} user-specified. These frequency weights are assigned on a per-feature basis and derived from the categories represented in each column of \code{data}.
#' @return The whole matrix of pairwise distances if \code{full.matrix=TRUE}, otherwise the corresponding \code{\link[stats]{dist}} object.
#' @references Huang, Z. (1997). A fast clustering algorithm to cluster very large categorical data sets in data mining. In H. Lu, H. Motoda, and H. Luu (Eds.), \emph{KDD: Techniques and Applications}, pp. 21-34. Singapore: World Scientific.
#' @author Keefe Murphy - <\email{keefe.murphy@@mu.ie}>
#' @seealso \code{\link{wKModes}}
#' @keywords utility
#' @importFrom matrixStats "rowSums2"
#' @importFrom WeightedCluster "wcAggregateCases" "wcSilhouetteObs"
#' @export
#' @usage 
#' dist_freqwH(data,
#'             full.matrix = TRUE)
#' @examples
#' suppressMessages(require(WeightedCluster))
#' set.seed(99)
#' # Load the MVAD data & aggregate the state sequences
#' data(mvad)
#' agg      <- wcAggregateCases(mvad[,17:86], weights=mvad$weight)
#' 
#' # Create a state sequence object without the first two (summer) time points
#' states   <- c("EM", "FE", "HE", "JL", "SC", "TR")
#' labels   <- c("Employment", "Further Education", "Higher Education", 
#'               "Joblessness", "School", "Training")
#' weights  <- agg$aggWeights
#' mvad.seq <- seqdef(mvad[agg$aggIndex, 17:86], 
#'                    states=states, labels=labels, weights=agg$aggWeights)
#' 
#' # Run k-modes with weights
#' resW     <- wKModes(mvad.seq, 2, weights=agg$aggWeights)
#' 
#' # Run k-modes with additional frequency weights
#' resF     <- wKModes(mvad.seq, 2, weights=agg$aggWeights, freq.weighted=TRUE)
#' 
#' # Examine the average silhouette widths of both weighted solutions
#' weighted.mean(wcSilhouetteObs(seqdist(mvad.seq, method="HAM"), resW$cluster, weights), weights)
#' # weighted.mean(wcSilhouetteObs(seqdist(mvad.seq, method="HAM"), resF$cluster, weights), weights)
#' weighted.mean(wcSilhouetteObs(dist_freqwH(mvad.seq), resF$cluster, weights), weights)
dist_freqwH   <- function(data, full.matrix = TRUE) {
  data        <- as.data.frame(data)
  n           <- nrow(data)
  P           <- ncol(data)
  isnumeric   <- vapply(data, is.numeric, logical(1L))
  isfactor    <- vapply(data, is.factor,  logical(1L))
  if(any(isfactor))             {
    levs      <- vector("list", P)
    for(j in which(isfactor))   {
      levsj   <- levels(data[,j])
      data[,j]    <- levsj[data[,j]]
      levs[[j]]   <- levsj
    }
  }
  if(any(isnumeric))            {
    lengths   <- vapply(data[,isnumeric], function(z) length(unique(z)), numeric(1L))
    if(any(lengths > 30))        warning("data has numeric columns with more than 30 different levels!", call.=FALSE, immediate.=TRUE)
  }
  Pseq        <- seq_len(P)
  frwts       <- lapply(Pseq, function(i) table(data[,i]))
  nX          <- matrix(NA, nrow=n, ncol=P)
  dX          <- matrix(NA, nrow=n, ncol=n)
  for(j in Pseq)       {
    nX[,j]    <- frwts[[j]][vapply(as.character(data[,j]), function(z) which(names(frwts[[j]]) == z), numeric(1L))]
  }    
  for(i in seq_len(n)) {
    di        <- vapply(Pseq, function(j) data[,j] != data[i,j], logical(n))
    wts       <- 1/nX[rep.int(i, n),] + 1/nX
    dX[i,]    <- rowSums2(di * wts)
  }
    if(isTRUE(full.matrix)) dX else stats::as.dist(dX)
}

#' @method print wKModes
#' @export
print.wKModes <- function(x, ...) {
  cat("K-modes clustering with ", length(x$size), " clusters of sizes ", paste(x$size, collapse = ", "), "\n", sep = "")
  cat("\nCluster modes:\n")
  print(x$modes,      ...)
  cat("\nClustering vector:\n")
  print(x$cluster,    ...)
  cat("\nWithin cluster simple-matching distance by cluster:\n")
  print(x$withindiff, ...)
  cat("\nTotal within cluster simple-matching distance over all clusters:\n")
  print(x$tot.withindiff, ...)
  cat("\nAvailable components:\n")
  print(names(x))
    invisible(x)
}

#' Predictions from MEDseq gating networks
#' 
#' Predicts mixing proportions from MEDseq gating networks. Effectively akin to predicting from a multinomial logistic regression via \code{\link[nnet]{multinom}}, although here the noise component (if any) is properly accounted for. So too are models with no gating covariates at all, or models with the equal mixing proportion constraint. Prior probabilities are returned by default.
#' @param object An object of class \code{"MEDgating"} (typically \code{x$gating}, where \code{x} is of class \code{"MEDseq"}).
#' @param newdata A matrix or data frame of test examples. If omitted, the fitted values are used.
#' @param type The type of output desired. The default (\code{"probs"}) returns prior probabilities, while \code{"class"} returns labels indicating the most likely group \emph{a priori}. Note that observations classified assigned the noise component (if any) are given a label of \code{0}.
#' @param keep.noise A logical indicating whether the output should acknowledge the noise component (if any). Defaults to \code{TRUE}; when \code{FALSE}, this column is discarded and the matrix of probabilities is renormalised accordingly.
#' @param droplevels A logical indicating whether unseen factor levels in categorical variables within \code{newdata} should be dropped (with \code{NA} predicted in their place). Defaults to \code{FALSE}.
#' @param ... Catches unused arguments or allows the \code{type} and \code{keep.noise} arguments to be passed through \code{fitted} and the \code{keep.noise} argument to be passed through \code{residuals}.
#' 
#' @return The return value depends on whether \code{newdata} is supplied or not and whether the model includes gating covariates to begin with. When \code{newdata} is not supplied, the fitted values are returned (as matrix if the model contained gating covariates, otherwise as a vector as per \code{x$params$tau}). If \code{newdata} is supplied, the output is always a matrix with the same number of rows as the \code{newdata}.
#' @details This function (unlike the \code{predict} method for \code{\link[nnet]{multinom}} on which \code{predict.MEDgating} is based) accounts for sampling weights and the various ways of treating gating covariates, equal mixing proportion constraints, and noise components, although its \code{type} argument defaults to \code{"probs"} rather than \code{"class"}.
#' 
#' @references Murphy, K., Murphy, T. B., Piccarreta, R., and Gormley, I. C. (2021). Clustering longitudinal life-course sequences using mixtures of exponential-distance models. \emph{Journal of the Royal Statistical Society: Series A (Statistics in Society)}, 184(4): 1414-1451. <\href{https://rss.onlinelibrary.wiley.com/doi/abs/10.1111/rssa.12712}{doi:10.1111/rssa.12712}>.
#' @author Keefe Murphy - <\email{keefe.murphy@@mu.ie}>
#' @seealso \code{\link[nnet]{multinom}}
#' @method predict MEDgating
#' @keywords utility
#' @importFrom nnet "multinom"
#' @export
#' @usage 
#' \method{predict}{MEDgating}(object,
#'         newdata = NULL,
#'         type = c("probs", "class"),
#'         keep.noise = TRUE,
#'         droplevels = FALSE,
#'         ...)
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
#' # Fit a model with weights and a gating covariate
#' # Have the probability of noise-component membership be constant
#' mod    <- MEDseq_fit(mvad.seq, G=11, modtype="UUN", weights=mvad$weights, 
#'                      gating=~ gcse5eq, covars=mvad.cov, noise.gate=FALSE)
#' (preds <- predict(mod$gating, newdata=mvad.cov[1:5,]))
#' 
#' # Note that the predictions are not the same as the multinom predict method
#' # in this instance, owing to the invocation of noise.gate=FALSE above
#' mod2   <- mod
#' class(mod2$gating) <- c("multinom", "nnet")
#' predict(mod2$gating, newdata=mvad.cov[1:5,], type="probs")
#' 
#' # We can make this function behave in the same way by invoking keep.noise=FALSE
#' predict(mod$gating, keep.noise=FALSE, newdata=mvad.cov[1:5,])}
predict.MEDgating <- function(object, newdata = NULL, type = c("probs", "class"), keep.noise = TRUE, droplevels = FALSE, ...) {
  if(!is.null(newdata) &&
     all(!is.data.frame(newdata), 
         !is.matrix(newdata)))   stop("'newdata' must be a matrix or data frame if supplied", call.=FALSE)
  if(!missing(type)    && 
     length(type)  > 1 || 
     !is.character(type))        stop("'type' must be a single character string",             call.=FALSE)
  if(length(keep.noise) > 1   ||
     !is.logical(keep.noise))    stop("'keep.noise' must be a single logical indicator",      call.=FALSE)
  if(length(droplevels) > 1   ||
     !is.logical(droplevels))    stop("'droplevels' must be a single logical indicator",      call.=FALSE)
  class(object)   <- class(object)[-1L]
  gat             <- attr(object, "Formula")
  fits            <- object$fitted.values
  if(attr(object, "DoUni")    &&
     gat          != 1)        {
    fits          <- if(is.matrix(fits)) fits[attr(object, "Disagg"),, drop=FALSE] else fits[attr(object, "Disagg")]
  }
  noise           <- attr(object, "Noise")
  G               <- attr(object, "G")
  gnames          <- paste0("Cluster",      if(noise) replace(seq_len(G), G, "0")                       else seq_len(G))
  if(gat          == "~1")     {
    fits          <- matrix(attr(object, "NoGate"), nrow=ifelse(is.null(newdata), 1L, NROW(newdata)), ncol=G, byrow=TRUE)
  } else if(!is.null(newdata)) { 
    fits          <- stats::predict(object, if(isTRUE(droplevels)) .drop_levels(object, newdata)        else newdata, type="probs")
    fits          <- if(is.matrix(fits)) fits else matrix(fits, nrow=1L, ncol=length(fits), byrow=TRUE)
  }
  if(all(noise, !attr(object, "NoiseGate")))   {
    fits          <- .tau_noise(fits, attr(object, "NoisePro"))
  }
  colnames(fits)  <- NULL
  if(isFALSE(keep.noise)      && noise)        {
    if(attr(object, "G0"))       stop("Nothing to return as the model has only a noise component: use keep.noise=TRUE", call.=FALSE)
    fits          <- .renorm_z(fits[,-G, drop=FALSE])
  }
    switch(EXPR=match.arg(type), 
           probs=provideDimnames(fits, base=list(ifelse(nrow(fits) == 1, "pro", ""), gnames)), {
      if(all(attr(object, "EqualPro"), G > 1)) {
       if(!all(noise, fits[G] == max(fits), !attr(object, "EqualNoise"), 
               keep.noise))      message("class predicted at random due to the equal mixing proportion constraint\n")  
       }
       CLS       <- max.col(fits)
         if(all(noise, isTRUE(keep.noise)))  replace(CLS, CLS == G, 0L) else CLS
    })
}
#' @rdname predict.MEDgating
#' @method fitted MEDgating
#' @keywords prediction utility
#' @importFrom nnet "multinom"
#' @usage 
#' \method{fitted}{MEDgating}(object,
#'        ...)
#' @export
fitted.MEDgating  <- function(object, ...) {
  args            <- c(list(object=object, newdata=NULL), as.list(match.call())[-1L])
  fits            <- do.call(predict.MEDgating, args[unique(names(args))])
  if(!is.null(fits)) return(fits)
}

#' @rdname predict.MEDgating
#' @method residuals MEDgating
#' @keywords prediction utility
#' @importFrom nnet "multinom"
#' @usage 
#' \method{residuals}{MEDgating}(object,
#'           ...)
#' @export
residuals.MEDgating  <- function(object, ...)   {
  dat.z           <- attr(object, "Data")
  args            <- c(list(object=object, type="probs", newdata=dat.z), as.list(match.call())[-1L])
  keep            <- !any(names(list(...))     == "keep.noise") || isTRUE(args$keep.noise)
  fits            <- do.call(fitted.MEDgating, args[unique(names(args))])
  dat.z           <- if(isTRUE(keep) || !attr(object, "Noise")) dat.z else .renorm_z(dat.z[,-ncol(dat.z), drop=FALSE])
  dat.z[is.nan(dat.z)]               <- 0L
  fits            <- if(nrow(fits)   == nrow(dat.z))             fits else matrix(fits, nrow=nrow(dat.z), ncol=ncol(dat.z), byrow=TRUE)
    provideDimnames(tryCatch(dat.z - fits, error=function(e) 1L - fits), base=list(as.character(seq_len(nrow(dat.z))), ""))
}

#' Show the NEWS file
#'
#' Show the \code{NEWS} file of the \code{MEDseq} package.
#' @return The \code{MEDseq} \code{NEWS} file, provided the session is interactive.
#' @export
#' @keywords utility
#'
#' @usage MEDseq_news()
#' @examples
#' MEDseq_news()
MEDseq_news       <- function() {
  newsfile        <- file.path(system.file(package  = "MEDseq"), "NEWS.md")
    if(interactive()) file.show(newsfile) else message("The session is not interactive\n")
}

#' Create a state sequence object
#' 
#' Create a state sequence object with attributes such as alphabet, colour palette, and state labels. This a minimalist copy of \code{\link[TraMineR]{seqdef}} in the \pkg{TraMineR} package.
#' @param data A data frame or matrix containing sequence data.
#' @param ... All other arguments; see \code{\link[TraMineR]{seqdef}} in the \pkg{TraMineR} package.
#'
#' @return An object of class \code{"stslist"}, for which dedicated \code{print} and \code{summary} methods are inherited from \pkg{TraMineR}.
#' @details This function exists only so experienced users of \pkg{MEDseq} and \pkg{TraMineR} can use the former without explicitly requiring the latter to be loaded. In particular, \code{\link{MEDseq_fit}} requires a state-sequence object of class \code{"stslist"} (as created by \code{seqdef}) as input. Users are encouraged to see the documentation at \code{\link[TraMineR:seqdef]{TraMineR::seqdef}} for complete details and further examples.
#' @keywords utility
#' @references Gabadinho, A., Ritschard, G., Mueller, N. S., and Studer, M. (2011). Analyzing and visualizing state sequences in R with \pkg{TraMineR}. \emph{Journal of Statistical Software}, 40(4): 1-37.
#' 
#' Gabadinho, A., Ritschard, G., Studer, M., and Mueller, N. S. (2010). Mining sequence data in \code{R} with the \pkg{TraMineR} package: a user's guide. \emph{Department of Econometrics and Laboratory of Demography, University of Geneva}.
#' @author Alexis Gabadinho and Gilbert Ritschard
#' @seealso \code{\link[TraMineR:seqdef]{TraMineR::seqdef}}, \code{\link{MEDseq_fit}}
#' @importFrom TraMineR "seqdef"
#' @usage seqdef(data, ...)
#' @export
#' @examples
#' data(mvad)
#' # Create a state sequence object with the first two (summer) time points removed
#' states   <- c("EM", "FE", "HE", "JL", "SC", "TR")
#' labels   <- c("Employment", "Further Education", "Higher Education", 
#'               "Joblessness", "School", "Training")
#' mvad.seq <- seqdef(mvad[,17:86], states=states, labels=labels, weights=mvad$weight)
seqdef <- function(data, ...) {
  TraMineR::seqdef(data, ...)
}
#