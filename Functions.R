dbs               <- function(z, tol = log(1E-100), weights = NULL, ...) {
  if(any(!is.matrix(z), !is.numeric(z)) ||
     ncol(z)      <= 1     ||
     nrow(z)      <= 1)          stop("'z' must be a numeric matrix with 2 or more columns & 2 or more rows", call.=FALSE)
  if(!missing(weights))     {
    N             <- nrow(z)
   if(!is.numeric(weights) ||
      length(weights)      != N) stop(paste0("'weights' must be a numeric vector of length N=", N), call.=FALSE)
   if(any(weights < 0)     || 
      any(!is.finite(weights)))  stop("'weights' must be positive and finite", call.=FALSE)
  }
  MAP             <- max.col(z)
  z               <- matrix(z[order(row(z), -z)], nrow(z), byrow=TRUE)
  l2              <- log(z[,2L])
  zz              <- log(z[,1L])     - l2
  zz.inf          <- is.infinite(zz) | l2 < tol
  ds              <- zz/max(abs(zz[!zz.inf]))
  ds[zz.inf]      <- 1L
  ds[is.nan(ds)]  <- 0L
  DS              <- cbind(cluster=MAP, dbs_width=ds)
  class(DS)       <- "MEDsil"
  msw             <- median(ds)
    return(list(silvals = DS, msw = msw, wmsw = ifelse(is.null(weights), msw, matrixStats::weightedMedian(ds, weights))))
}

get_results                   <- function(x, what = c("z", "MAP", "sils"), rank = 1L, criterion = c("bic", "icl", "aic", "cv", "nec", "dbs", "loglik"), G = NULL, modtype = NULL, noise = TRUE, ...) {
    UseMethod("get_results")
}

get_results.MEDseq            <- function(x, what = c("z", "MAP", "sils"), rank = 1L, criterion = c("bic", "icl", "aic", "cv", "nec", "dbs", "loglik"), G = NULL, modtype = NULL, noise = TRUE, ...) {
  x               <- if(inherits(x, "MEDseqCompare")) x$optimal else x
  what            <- match.arg(what)
  minG            <- 1L  + (what == "sils")
  if(!(missing(G) -> m.G)     &&
    (length(G)    != 1        ||
     !is.numeric(G)           ||
     (G < minG    || floor(G) != G))) {
    if(what == "sils")   {       stop("'G' must be a single integer > 1 when 'what=sils'", call.=FALSE)
    } else                       stop("'G' must be a single integer >= 1",                 call.=FALSE)
  }  
  if(!(missing(modtype) ->
       m.M) &&
    (length(modtype) > 1      ||
     !is.character(modtype)))    stop("'modtype' must be a single character string", call.=FALSE)
  if(length(noise) > 1        ||
     !is.logical(noise))         stop("'noise' must be a single logical indicator",  call.=FALSE)
  if(!missing(criterion)      ||
     !missing(rank)           ||
     any(m.G, m.M))            {
    criterion     <- match.arg(criterion)
    if((criterion == "nec"    ||
       criterion  == "dbs")   &&
      !m.G  &&  G == 1)          stop(paste0("Can't select based on the ", toupper(criterion), " criterion when G=1"), call.=FALSE)
    if(criterion  == "cv"     &&
       !attr(x, "CV"))           stop("Can't select based on the CV criterion as cross-validated likelihood wasn't performed", call.=FALSE)
    tmp           <- switch(EXPR=criterion, bic=x$BIC, icl=x$ICL, aic=x$AIC, cv=x$CV, nec=x$NEC, dbs=x$DBS, loglik=x$LOGLIK)
    if(!noise)     {
      tmp         <- tmp[,colnames(tmp) %in% c("CC", "UC", "CU", "UU"), drop=FALSE]
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
  switch(EXPR=what, sils=    {
    SILS          <- x$SILS
    if(!(G  %in%
       as.numeric(names(SILS)))) stop("Invalid 'G' value", call.=FALSE)
    S             <- SILS[[as.character(G)]]
    s.ind         <- if(noise) names(S) else names(S)[names(S) %in% c("CC", "UC", "CU", "UU")]
    if(!(modtype %in% s.ind))    stop("Invalid 'modtype'", call.=FALSE)
    res           <- S[[modtype]]
    if(anyNA(res))               message("Selected model didn't converge: no silhouettes available\n")
  }, {
    ZS            <- x$ZS
    if(!(G  %in%
       as.numeric(names(ZS))))   stop("Invalid 'G' value", call.=FALSE)
    Z             <- ZS[[as.character(G)]]
    z.ind         <- if(noise) names(Z) else names(Z)[names(Z) %in% c("CC", "UC", "CU", "UU")]
    if(!(modtype %in% z.ind))    stop("Invalid 'modtype'", call.=FALSE)
    res           <- Z[[modtype]]
    if(anyNA(res))               message("Selected model didn't converge: no partition available\n")
    if(what == "MAP")       {
      MAP         <- max.col(res)
      res         <- if(noise) replace(MAP, MAP == G, 0L) else MAP
    }
  })
  attr(res, "G")           <- G
  attr(res, "ModelType")   <- modtype
  attr(res, "Noise")       <- is.element(modtype, c("CCN", "UCN", "CUN", "UUN"))
    return(res)
}

MEDseq_compare    <- function(..., criterion = c("bic", "icl", "aic", "cv", "nec", "dbs"), pick = 10L, optimal.only = FALSE) {
  crit.miss       <- missing(criterion)
  if(!missing(criterion)  && (length(criterion) > 1 ||
     !is.character(criterion)))  stop("'criterion' must be a single character string", call.=FALSE)
  criterion       <- match.arg(criterion)
  num.miss        <- missing(pick)
  opt.miss        <- missing(optimal.only)
  if(length(pick) != 1    ||
     !is.numeric(pick))          stop("'pick' must be a single number", call.=FALSE)
  if(floor(pick)  != pick ||
     pick          < 1)          stop("'pick' must be a strictly positive integer", call.=FALSE)
  if(length(optimal.only)  > 1 ||
     !is.logical(optimal.only))  stop("'optimal.only' must be a single logical indicator", call.=FALSE)
  call            <- match.call(expand.dots=TRUE)[-1L]
  call            <- if(crit.miss) call else call[-which(names(call) == "criterion")]
  call            <- if(num.miss)  call else call[-which(names(call) == "pick")]
  call            <- if(opt.miss)  call else call[-which(names(call) == "optimal.only")]
  len.call        <- length(as.list(call))
  if(len.call     == 1    && inherits(..., "list") && !inherits(..., "MEDseq")) {
    mod.names     <- unique(names(...))
    MEDs          <- as.list(...)[mod.names]
    if(is.null(mod.names))       stop("When supplying models as a list, every element of the list must be named", call.=FALSE)
  } else           {
    mod.names     <- vapply(call, deparse, character(1L))
    MEDs          <- stats::setNames(list(...), mod.names)
    mod.names     <- unique(mod.names)
    MEDs          <- MEDs[mod.names]
  }
  Mclass          <- vapply(MEDs, class,          character(1L))
  if(any(Mclass   != "MEDseq"))  stop("All models must be of class 'MEDseq'!", call.=FALSE)
  if(length(unique(lapply(MEDs, "[[", 
     "seqs")))    != 1)          stop("All models being compared must have been fit to the same data set!", call.=FALSE)
  title           <- "Mixtures of Exponential-Distance Models with Covariates"
  dat.name        <- deparse(MEDs[[1L]]$call$seqs)
  gate.x          <- lapply(MEDs, "[[", "gating")
  algo            <- sapply(MEDs,   attr, "Algo")
  equalNoise      <- sapply(MEDs,   attr, "EqualNoise")
  equalPro        <- sapply(MEDs,   attr, "EqualPro")
  noise.gate      <- sapply(MEDs,   attr, "NoiseGate")
  weights         <- sapply(MEDs,   attr, "Weighted")
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
  cvnull          <- vapply(CVs,  is.null, logical(1L))
  necnull         <- vapply(NECs, is.null, logical(1L))
  dbsnull         <- vapply(DBSs, is.null, logical(1L))
  if(all(cvnull)  && 
     criterion    == "cv")       stop("'criterion' cannot be 'cv' when cross-validation was not performed for any of the supplied models", call.=FALSE)
  if(all(necnull) &&             
     criterion    == "nec")      stop("'criterion' cannot be 'nec' when all models being compared contain only 1 component", call.=FALSE)
  if(all(dbsnull) &&             
     criterion    == "dbs")      stop("'criterion' cannot be 'dbs' when all models being compared contain only 1 component", call.=FALSE)
  choice          <- max(lengths(BICs))
  bics            <- lapply(BICs, function(x) .pick_MEDCrit(x, choice)$crits)
  icls            <- lapply(ICLs, function(x) .pick_MEDCrit(x, choice)$crits)
  aics            <- lapply(AICs, function(x) .pick_MEDCrit(x, choice)$crits)
  llxs            <- lapply(LLxs, function(x) .pick_MEDCrit(x, choice)$crits)
  dfxs            <- lapply(DFxs, function(x) .pick_MEDCrit(x, choice)$crits)
  itxs            <- lapply(ITxs, function(x) .pick_MEDCrit(x, choice)$crits)
  cvs             <- lapply(CVs,  function(x) if(!is.null(x)) .pick_MEDCrit(x, choice)$crits)[!cvnull]
  necs            <- lapply(NECs, function(x) if(!is.null(x)) .pick_MEDCrit(x, choice)$crits)[!necnull]
  dbss            <- lapply(DBSs, function(x) if(!is.null(x)) .pick_MEDCrit(x, choice)$crits)[!dbsnull]
  if(optimal.only) {
    opt.names     <- names(.crits_names(lapply(switch(EXPR=criterion, bic=bics, icl=icls, aic=aics, cv=cvs, nec=necs, dbs=dbss), "[", 1L)))
  }
  bics            <- .crits_names(bics)
  icls            <- .crits_names(icls)
  aics            <- .crits_names(aics)
  llxs            <- .crits_names(llxs)
  dfxs            <- .crits_names(dfxs)
  itxs            <- .crits_names(itxs)
  cvs             <- .crits_names(cvs)
  necs            <- .crits_names(necs)
  dbss            <- .crits_names(dbss)
  if(criterion    == "cv"  &&
    (length(cvs)  != 
     length(bics)))              warning("Discarding models for which the CV criterion was not computed\n",  call.=FALSE, immediate.=TRUE)
  if(criterion    == "nec" &&
    (length(necs) != 
     length(bics)))              warning("Discarding models for which the NEC criterion was not computed\n", call.=FALSE, immediate.=TRUE)
  if(criterion    == "dbs" &&
    (length(dbss) != 
     length(bics)))              warning("Discarding models for which the DBS criterion was not computed\n", call.=FALSE, immediate.=TRUE)
  if(optimal.only) {
    bics          <- bics[names(bics) %in% opt.names]
    icls          <- icls[names(icls) %in% opt.names]
    aics          <- aics[names(aics) %in% opt.names]
    llxs          <- llxs[names(llxs) %in% opt.names]
    dfxs          <- dfxs[names(dfxs) %in% opt.names]
    itxs          <- itxs[names(itxs) %in% opt.names]
    cvs           <- cvs[names(cvs)   %in% opt.names]
    necs          <- necs[names(necs) %in% opt.names]
    dbss          <- dbss[names(dbss) %in% opt.names]
  }
  crits           <- switch(EXPR=criterion, bic=bics, icl=icls, aic=aics, cv=cvs, nec=necs, dbs=dbss)
  pick            <- min(pick, length(crits))
  max.crits       <- sort(crits, decreasing=criterion != "nec")[seq_len(pick)]
  if(length(unique(max.crits))  < pick) {
    ties          <- max.crits == max.crits[1L]
    if(any(ties[-1L]))     {     warning(paste0("Ties for the optimal model exist according to the '", criterion, "' criterion: choosing the most parsimonious model\n"), call.=FALSE, immediate.=TRUE)
      df.ties     <- dfxs[names(max.crits)][which(ties)]
      max.crits[ties]     <- max.crits[order(df.ties)]
      if(any((df.ties     == df.ties[1L])[-1L])) {
        max.crits[ties]   <- max.crits[order(as.numeric(gsub(".*,", "", names(max.crits[ties]))))]
      }
    } else                       warning(paste0("Ties exist according to the '", criterion, "' criterion\n"), call.=FALSE, immediate.=TRUE)
  }
  max.names       <- names(max.crits)
  crit.names      <- gsub("\\|.*", "",          max.names)
  G               <- as.numeric(gsub(".*,", "", max.names))
  gating          <- unname(unlist(gating[crit.names]))
  modelNames      <- gsub(",.*", "", gsub(".*\\|", "", max.names))
  best.model      <- MEDs[[crit.names[1L]]]
  if(best.model$modName != modelNames[1L] || best.model$G != G[1L]) {
    cat("Re-fitting optimal model due to mismatched 'criterion'...\n\n")
    old.call    <- best.model$call
    old.call    <- c(as.list(old.call)[1L], list(criterion=criterion), as.list(old.call)[-1L])
    old.call    <- as.call(old.call[!duplicated(names(old.call))])
    best.call   <- c(list(data=best.model$data, l.meth=modelNames[1L], G=G[1L], criterion="bic", verbose=FALSE, do.cv=FALSE, do.nec=FALSE), as.list(old.call[-1L]))
    best.mod    <- try(do.call(MEDseq_fit, best.call[!duplicated(names(best.call))]), silent=TRUE)
    if(!inherits(best.model, "try-error")) {
      best.model$call               <- old.call
      best.model$modName            <- best.mod$modName
      best.model$G                  <- best.mod$G
      best.model$bic                <- best.mod$bic
      best.model$icl                <- best.mod$icl
      best.model$aic                <- best.mod$aic
      best.model$cv                 <- if(attr(best.model, "CV"))  best.model$CV[best.mod$G,best.mod$modName]
      best.model$nec                <- if(attr(best.model, "NEC")) best.model$NEC[which(best.mod$G == as.numeric(rownames(best.model$NEC))),best.mod$modName]
      best.model$dbs                <- if(attr(best.model, "DBS")) best.model$DBS[which(best.mod$G == as.numeric(rownames(best.model$DBS))),best.mod$modName]
      best.model$gating             <- best.mod$gating
      best.model$loglik             <- best.mod$loglik
      best.model$df                 <- best.mod$df
      best.model$iters              <- best.mod$iters
      best.model$params             <- best.mod$params
      best.model$z                  <- best.mod$z
      best.model$MAP                <- best.mod$MAP
      best.model$uncert             <- best.mod$uncert
      attributes(best.model)        <- attributes(best.mod)
    } else best.model               <- paste0("Failed to re-fit the optimal model: ", gsub("\"", "'", deparse(old.call, width.cutoff=500L), fixed=TRUE))
  }
  gating[gating == "~1" | G   == 1] <- "None"
  noise         <- modelNames %in% c("CCN", "UCN", "CUN", "UUN")
  noise.gate    <- ifelse(!noise, NA, noise.gate[crit.names])
  equalPro      <- replace(unname(equalPro[crit.names]), gating != "None" | G == 1, NA)
  equalNoise    <- ifelse(!noise | G == 1, NA, equalNoise[crit.names] & vapply(equalPro, isTRUE, logical(1L)))
  comp          <- list(title = title, data = dat.name, optimal = best.model, pick = pick, MEDNames = crit.names, modelNames = modelNames, G = as.integer(G), 
                        df = as.integer(unname(dfxs[max.names])), iters = as.integer(unname(itxs[max.names])), bic = unname(bics[max.names]), icl = unname(icls[max.names]), 
                        aic = unname(aics[max.names]), cv = unname(cvs[max.names]), nec = replace(unname(necs[max.names]), G == 1, NA), dbs = replace(unname(dbss[max.names]), G == 1, NA), 
                        loglik = unname(llxs[max.names]), gating = gating, algo = unname(algo[crit.names]), weights = unname(weights[crit.names]), equalPro = equalPro, noise = unname(noise), 
                        noise.gate = unname(replace(noise.gate, gating == "None" | G <= 2, NA)), equalNoise = unname(replace(equalNoise, !equalPro | is.na(equalPro), NA)))
  class(comp)   <- c("MEDseqCompare", "MEDseq")
  bic.tmp       <- sapply(BICs, as.vector)
  attr(comp, "Crit")   <- criterion
  attr(comp, "Opt")    <- optimal.only
  attr(comp, "NMods")  <- c(tried = sum(vapply(bic.tmp, function(x) length(x[!is.na(x)]),    numeric(1L))),
                            ran   = sum(vapply(bic.tmp, function(x) length(x[is.finite(x)]), numeric(1L))))
    comp
}

MEDseq_control    <- function(algo = c("EM", "CEM", "cemEM"), init.z = c("kmedoids", "hc", "random"), dist.mat = NULL, unique = FALSE, nstarts = 1L,
                              criterion = c("bic", "icl", "aic", "cv", "nec", "dbs"), do.cv = FALSE, do.nec = FALSE, nfolds = 10L, stopping = c("aitken", "relative"),
                              tau0 = NULL, opti = c("mode", "first", "GA", "medoid"), ordering = c("none", "decreasing", "increasing"), noise.gate = TRUE,
                              equalPro = FALSE, equalNoise = FALSE, tol = c(1E-05, 1E-08), itmax = c(.Machine$integer.max, 100L), nonzero = TRUE, verbose = TRUE, ...) {
  miss.args                <- list(tau0=missing(tau0), unique=missing(unique))
  if(!missing(algo)        &&
    (length(algo)      > 1 ||
     !is.character(algo)))       stop("'algo' must be a character vector of length 1",      call.=FALSE)
  if(!missing(init.z)      &&
    (length(init.z)    > 1 ||
     !is.character(init.z)))     stop("'init.z' must be a character vector of length 1",    call.=FALSE)
  init.z                   <- match.arg(init.z)
  if(!missing(dist.mat))    {
    dist.mat               <- tryCatch(suppressWarnings(as.dist(dist.mat)), error=function(e)   {
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
  if(length(do.cv)     > 1 ||
     !is.logical(do.cv))         stop("'do.cv' must be a single logical indicator",         call.=FALSE)
  if(length(do.nec)    > 1 ||
     !is.logical(do.nec))        stop("'do.nec' must be a single logical indicator",        call.=FALSE)
  if(!missing(stopping)    &&
    (length(stopping)  > 1 ||
     !is.character(stopping)))   stop("'stopping' must be a character vector of length 1",  call.=FALSE)
  if(!miss.args$tau0       &&
    (length(tau0)      > 1 ||
     !is.numeric(tau0)     ||
     tau0 < 0 || tau0 >= 1))     stop("'tau0' must be a scalar in the interval [0, 1)",     call.=FALSE)
  if(!missing(opti)        &&
    (length(opti)      > 1 ||
     !is.character(opti)))       stop("'opti' must be a character vector of length 1",      call.=FALSE)
  if(!missing(ordering)    &&
    (length(ordering)  > 1 ||
     !is.character(ordering)))   stop("'ordering' must be a character vector of length 1",  call.=FALSE)
  if(length(noise.gate)     > 1 ||
     !is.logical(noise.gate))    stop("'noise.gate' must be a single logical indicator",    call.=FALSE)
  if(length(equalPro)  > 1 ||
     !is.logical(equalPro))      stop("'equalPro' must be a single logical indicator",      call.=FALSE)
  if(length(equalNoise)     > 1 ||
     !is.logical(equalNoise))    stop("'equalNoise' must be a single logical indicator",    call.=FALSE)
  if((len.tol     <- 
      length(tol)) > 2     ||
     !is.numeric(tol))           stop("'tol' must be a numeric vector of length at most 2", call.=FALSE)
  if(any(tol   < 0,
         tol  >= 1))             stop("'tol' must be in the interval [0, 1)",               call.=FALSE)
  if(len.tol  == 1)    tol <- rep(tol, 2L)
  if(length(itmax)         == 1) {
    itmax     <- c(itmax, 100L)
  } else if(length(itmax)  != 2) stop("'itmax' must be of length 2",                        call.=FALSE)
  if(!is.numeric(itmax)    ||
     any(floor(itmax) != itmax) ||
     any(itmax    <= 0))         stop("'itmax' must contain strictly positive integers",    call.=FALSE)
  inf         <- is.infinite(itmax)
  if(any(inf))   itmax[inf]     <- .Machine$integer.max
  itmax[1L]   <- ifelse(itmax[1L] == .Machine$integer.max, itmax[1L], itmax[1L] + 2L)
  if(length(nonzero)   > 1 ||
     !is.logical(nonzero))       stop("'nonzero' must be a single logical indicator",       call.=FALSE)
  if(length(verbose)   > 1 ||
     !is.logical(verbose))       stop("'verbose' must be a single logical indicator",       call.=FALSE)
  control                  <- list(algo = match.arg(algo), init.z = init.z, dist.mat = dist.mat, nstarts = nstarts, criterion = match.arg(criterion), nfolds = nfolds, do.cv = do.cv, 
                                   do.nec = do.nec, stopping = match.arg(stopping), tau0 = tau0, opti = match.arg(opti), ordering = match.arg(ordering), noise.gate = noise.gate, unique = unique,
                                   equalPro = equalPro, equalNoise = equalNoise, tol = tol[1L], g.tol = tol[2L], itmax = itmax[1L], g.itmax = itmax[2L], nonzero = nonzero, verbose = verbose)
  attr(control, "missing") <- miss.args
    return(control)
}

MEDseq_fit        <- function(seqs, G = 1L:9L, l.meth = c("CC", "UC", "CU", "UU", "CCN", "UCN", "CUN", "UUN"), 
                              gating = NULL, covars = NULL, weights = NULL, ctrl = MEDseq_control(...), ...) {
  call            <- match.call()
  if(!inherits(seqs, "stslist")) stop("'seqs' must be of class 'stslist'",        call.=FALSE)
  if(any(seqs     ==
         attr(seqs, "nr")))      stop("Missing values in 'seqs' are not allowed", call.=FALSE)
  SEQ             <- apply(.fac_to_num(seqs), 1L, .num_to_char)
  seqX            <- seqs
  levs            <- attr(seqs, "alphabet")
  attr(SEQ, "N")  <- N   <- attr(SEQ, "W") <- nrow(seqs)
  attr(SEQ, "P")  <- P   <- ncol(seqs)
  attr(SEQ, "V")  <- V   <- length(levs)
  attr(SEQ, "V1") <- V1  <- V - 1L
  attr(SEQ, "V1V")       <- V1/V
  attr(SEQ, "logV1")     <- log(V1)
  if(any(c(N, P, V)      <= 1))  stop("The number of sequences, the sequence length, and the sequence vocabulary must all be > 1", call.=FALSE)
  if(!is.character(l.meth))      stop("'l.meth' must be a character vector of length 1",       call.=FALSE)
  l.meth          <- match.arg(l.meth, several.ok=TRUE)
  algo            <- ctrl$algo
  criterion       <- ctrl$criterion
  init.z          <- ctrl$init.z
  dist.mat        <- ctrl$dist.mat
  nstarts         <- switch(EXPR=init.z, random=ctrl$nstarts, 1L)
  startseq        <- seq_len(nstarts)
  equalPro        <- ctrl$equalPro
  equalNoise      <- ctrl$equalNoise
  noise.gate      <- ctrl$noise.gate
  nonzero         <- ctrl$nonzero
  verbose         <- ctrl$verbose
  ctrl$warn       <- TRUE
  x.ctrl          <- list(equalPro=equalPro, noise.gate=noise.gate, equalNoise=equalNoise)
  ctrl$ordering   <- ifelse(ctrl$opti == "first", ctrl$ordering, "none")
  miss.args       <- attr(ctrl, "missing")
  covmiss         <- missing(covars)
  mt1             <- unique(vapply(l.meth, function(lx) switch(EXPR=lx, CC=, UC="CC", CU=, UU="CU", "CCN"), character(1L)))
  mt2             <- unique(vapply(l.meth, function(lx) switch(EXPR=lx, UCN="CCN", UUN="CUN", lx),          character(1L)))
  mtg             <- unique(l.meth)
  n.meths         <- c("CCN", "UCN", "CUN", "UUN")
  n.meth          <- l.meth %in% n.meths
  if(!(tmiss      <- miss.args$tau0) && 
     ctrl$tau0    == 0)        {
    if(all(n.meth))              stop("'tau0' is zero: models with noise component will not be fitted",      call.=FALSE)
    if(any(n.meth))              warning("'tau0' is zero: models with noise component will not be fitted\n", call.=FALSE, immediate.=TRUE)
    l.meths       <- c("CC", "UC", "CU", "UU")
  } else l.meths  <- c("CC", "UC", "CU", "UU", "CCN", "UCN", "CUN", "UUN")
  if(is.null(dist.mat))        { 
    dist.mat      <- suppressMessages(seqdist(seqs, "HAM", full.matrix=FALSE))
  } else if((sqrt(length(dist.mat) * 2L + N) !=
             N))                 stop("Invalid 'dist.mat' dimensions", call.=FALSE)
  
  if((gate.x      <- !missing(gating))) {
    if(!inherits(gating, 
                 "formula"))     stop("'gating' must be a formula", call.=FALSE)
    if(!covmiss   &&
       !is.data.frame(covars))   stop("'covars' must be a data.frame if supplied", call.=FALSE)
    if(inherits(try(stats::terms(gating), silent=TRUE), "try-error")) {
      if(covmiss)                stop("Can't use '.' in 'gating' formula without supplying 'covars' argument", call.=FALSE)
      gating      <- stats::reformulate(attr(stats::terms(gating, data=covars), "term.labels"), response="z")
    }
    gating        <- tryCatch(stats::update.formula(stats::as.formula(gating), z ~ .), error=function(e) {
                                 stop("Invalid 'gating' network formula supplied", call.=FALSE) })
    if(gating[[3L]]      == 1) { 
      if(verbose)                message("Not including gating network covariates with only intercept on gating formula RHS\n")
      gate.x      <- FALSE
    }
    gate.names    <- labels(stats::terms(gating))
  } 
  gate.names      <- if(gate.x)  gate.names[!is.na(gate.names)]
  if(!covmiss)     {
    if(!all(gate.names  %in% 
            colnames(covars)))   stop("Supplied gating covariates not found in supplied 'covars'", call.=FALSE)
    covars        <- if(gate.x)  covars[,gate.names, drop=FALSE] else as.data.frame(matrix(0L, nrow=N, ncol=0L))
  } else {
    if(any(grepl("\\$", 
                 gate.names)))   stop("Don't supply covariates to the gating network using the $ operator: use the 'covars' argument instead", call.=FALSE)
    covars        <- if(gate.x)  stats::model.frame(gating[-2L]) else as.data.frame(matrix(0L, nrow=N, ncol=0L))
  }
  if(nrow(covars) != N)          stop("'gating' covariates must contain the same number of rows as 'seqs'", call.=FALSE)
  glogi           <- vapply(covars, is.logical, logical(1L))
  covars[,glogi]  <- sapply(covars[,glogi], as.factor)
  if(covmiss)      {
    covars        <- data.frame(if(ncol(covars) > 0) covars[,gate.names] else covars, stringsAsFactors=TRUE)
  }
  
  if(ctrl$do.wts  <- do.wts   <- 
     !missing(weights))   {
    if(is.null(weights))         stop("Invalid 'weights' supplied", call.=FALSE)
    if(!is.numeric(weights)   ||
       length(weights)   != N)   stop(paste0("'weights' must be a numeric vector of length N=", N), call.=FALSE)
    if(any(weights < 0)  || 
       any(!is.finite(weights))) stop("'weights' must be positive and finite", call.=FALSE)
    if(ctrl$do.wts       <- 
       do.wts     <- (length(unique(weights)) > 1)) {
      attr(SEQ, "Weights")    <- weights
      attr(SEQ, "W")          <- sum(weights)
      do.uni      <- ifelse(miss.args$unique, TRUE, ctrl$unique)
    }
  } else do.uni   <- ctrl$unique
  DF              <- seqs        
  DF              <- if(gate.x)  cbind(DF, covars)  else DF
  DF              <- if(do.wts)  cbind(DF, weights) else DF
  dup.ind         <- !duplicated(DF)
  sum.dup         <- sum(dup.ind)
  if(sum.dup < N  && !do.uni)   message(paste0("Number of unique observations (", sum.dup, ") is less than N (", N, "): Consider setting 'unique'=TRUE\n"))
  if(do.uni       <- do.uni   &&
     sum.dup < N)  {            message(paste0("Proceeding with ", sum.dup, " unique observations, out of N=", N, "\n"))
    agg.DF        <- merge(data.frame(cbind(id=seq_len(N)), DF), data.frame(aggregate(cbind(DF[0L], count=1L), DF, length)), by=colnames(DF), sort=FALSE)
    agg.id        <- order(agg.DF$id)
    c2            <- agg.DF$count[agg.id]
    agg.DF        <- agg.DF[agg.id[dup.ind],, drop=FALSE]
    counts        <- agg.DF$count
    weights       <- if(do.wts) counts * agg.DF$weights else counts
    if(ctrl$do.wts       <- (length(unique(weights)) > 1)) {
      dist.mat2   <- dist.mat
      w2          <- rep(0L, N)
      w2[dup.ind]        <- weights
      dist.mat    <- as.dist(as.matrix(dist.mat)[dup.ind,dup.ind])
      seqs        <- seqs[dup.ind,,   drop=FALSE]
      atts        <- attributes(SEQ)
      SEQ         <- apply(.fac_to_num(seqs), 1L, .num_to_char)
      attributes(SEQ)         <- atts
      attr(SEQ, "N")     <- N <- nrow(seqs)
      attr(SEQ, "Weights")    <- weights
      attr(SEQ, "W")          <- sum(weights)
      covars      <- covars[dup.ind,, drop=FALSE]
      dis.agg     <- rep(seq_len(N), counts)[agg.id]
    }
  }
  else           {
    dup.ind       <- rep(FALSE, N)
    w2            <- weights
    dist.mat2     <- dist.mat
  }
  
  if(any(G        != floor(G))    &&
     any(G         < 1))         stop("'G' must be strictly positive", call.=FALSE)
  if(any(G        >= N))       {
    G             <- G[G <= N]
    if(length(G)   > 1)        { warning("Removing G values >= the number of observations\n",  call.=FALSE, immediate.=TRUE)
    } else                       stop("G values must be less than the number of observations", call.=FALSE)
  }
  G               <- rG  <- unique(sort(as.integer(G)))
  do.nec          <- ctrl$do.nec
  if(!do.nec      &&
     algo         != "CEM"    &&
     criterion    == "nec"    &&
     !all(G == 1))             { message("Forcing 'do.nec' to TRUE as criterion='nec'\n")
    do.nec        <- TRUE
  }
  if(do.nec &&
     algo         != "CEM"    &&
     !any(G == 1))             { message("Forcing G=1 models to be fitted for NEC criterion computation\n")
    G             <- rG  <- unique(c(1L, G))  
  }
  if(all(G  == 1)) {             message("Density-based silhouettes not computed as only single component models are being fit\n")
    do.dbs        <- FALSE
    if(criterion  == "dbs")      stop("DBS criterion cannot be used to select among only single-component models", call.=FALSE)
    if(do.nec     && 
       !(do.nec   <- FALSE))     message("Forcing 'do.nec' to FALSE as only single component models are being fit\n")
    if(criterion  == "nec")      stop("NEC criterion cannot be used to select among only single-component models", call.=FALSE)
  } else if(algo  == "CEM")    { message("Density-based silhouettes not computed as the CEM algorithm is employed\n")
    do.dbs        <- FALSE
    if(criterion  == "dbs")      stop("DBS criterion cannot be used to select among models fitted via CEM", call.=FALSE)
    if(do.nec     &&
       !(do.nec   <- FALSE))     message("Forcing 'do.nec' to FALSE as models are being fit via CEM\n")
    if(criterion  == "nec")      stop("NEC criterion cannot be used to select among models fitted via CEM", call.=FALSE)
  } else do.dbs   <- TRUE
  if(any(G > 1L))  {
    G1            <- any(G == 1L)
    G2            <- any(G == 2L)
    GG            <- any(G  > 2L)
    all.mod       <- if(all(G1, G2, GG)) unique(c(mtg, mt2, mt1)) else if(all(G1, G2)) unique(c(mt2, mt1)) else if(all(G1, GG)) unique(c(mtg, mt1)) else if(all(G2, GG)) unique(c(mtg, mt2)) else if(G2) mt2 else mtg
    all.mod       <- l.meths[l.meths %in% all.mod]
    if(init.z     == "hc")     {
      hcZ         <- if(do.wts) hclust(dist.mat2, method="ward.D2", members=w2) else cluster::agnes(dist.mat2, diss=TRUE, method="ward")
    }
  } else all.mod  <- l.meths[l.meths %in% mt1]
  len.G           <- length(G)
  multi           <- length(all.mod) > 1
  cvsel           <- ctrl$do.cv
  ctrl$do.cv      <- FALSE
  if(!cvsel       && 
     criterion    == "cv")     { message("Forcing 'do.cv' to TRUE as criterion='cv'\n")
    cvsel         <- TRUE
  }
  if(ctrl$numseq  <- any(c("CU", "UU", "CUN", "UUN") %in% all.mod, ctrl$opti == "mode", ctrl$ordering != "none")) {
    numseq        <- sapply(SEQ, .char_to_num)
    attr(numseq, "P")    <- P
  } else numseq   <- NULL
  
  BICs            <-
  ICLs            <-
  AICs            <-
  NECs            <- 
  DBSs            <- 
  LL.x            <- 
  DF.x            <-
  IT.x            <-
  Nzero.x         <-
  Ninfty.x        <- provideDimnames(matrix(NA, nrow=len.G, ncol=length(all.mod)), base=list(as.character(G), all.mod))
  ZS              <- 
  SILS            <- replicate(len.G, list())
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

  gate.G          <- matrix(ifelse(rG > 1, gate.x, FALSE), nrow=2L, ncol=length(rG), byrow=TRUE)
  if(gate.x)       {
    Gn            <- G - !noise.gate
    if((verbose   && gate.x)       &&
      ((any(Gn    <= 1)  && noise) ||
       any(G      <= 1))) {      message(paste0("Can't include gating network covariates ", ifelse(noise.gate, "in a single component mixture", "where G is less than 3 when 'noise.gate' is FALSE\n")))
     gate.G[2L,Gn <= 1]  <- FALSE
    }
  } else           {
    Gn            <- G
    gating        <- stats::as.formula(z ~ 1)
    environment(gating)  <- environment()
  }
  noise.gate      <- ifelse(gate.G, noise.gate, TRUE)
  if(all(equalPro, gate.x)) { 
    if(verbose)                  message("Can't constrain mixing proportions to be equal when gating covariates are supplied\n")
    equalPro      <- FALSE
  }
  equal.tau       <- rbind(ifelse(G == 1, TRUE, equalPro), ifelse(Gn < 1, TRUE, equalPro)) & !gate.G
  equal.n0        <- (rbind(G == 1, Gn == 1) | equalNoise) & equal.tau
  attr(covars, "Gating") <- gate.names
  colnames(covars)       <- if(gate.x) gate.names
  if(!identical(gating, 
                .drop_constants(covars, 
                gating)))        stop("Constant columns exist in gating formula; remove offending gating covariate(s) and try again", call.=FALSE)
  G.last          <- G[len.G]
  
  for(g  in G)     {
    if(isTRUE(verbose))     {    cat(paste0("\n", g, " cluster model", ifelse(multi, "s", ""), " -\n"))
      last.G      <- g   == G.last
    }
    if(ctrl$numseq)         {
      attr(numseq,   "G")  <- g
    } 
    attr(SEQ, "G")         <- g
    h             <- which(G   == g)
    g0 <- attr(SEQ, "G0")  <- g - noise
    if(isTRUE(noise))   {
      tau0        <- ifelse(tmiss, 1/g, ctrl$tau0)
    }

    if(g  > 1)     {
      algog       <- algo
      if(init.z   == "random"  &&
         nstarts   > 1)     {
        if(isTRUE(nonoise)) {
          zg      <- replicate(nstarts, list(mclust::unmap(sample(seq_len(g),  size=N, replace=TRUE), groups=seq_len(g))))
        }
        if(isTRUE(noise))   {
          if(g0    > 1)     {
            zg0   <- replicate(nstarts, list(mclust::unmap(sample(seq_len(g0), size=N, replace=TRUE), groups=seq_len(g0))))
            zg0   <- lapply(zg0, function(x) cbind(x * (1 - tau0), tau0))
          } else   {
            zg0   <- matrix(tau0, nrow=N, ncol=2L)
          }
        }
      } else       {
        if(isTRUE(nonoise)) {
          zg      <- mclust::unmap(switch(EXPR=init.z, 
                                          random=sample(seq_len(g),  size=N, replace=TRUE),
                                          kmedoids= if(do.wts) {
                                            zz <- WeightedCluster::wcKMedoids(dist.mat, k=g,  weights=weights, cluster.only=TRUE)
                                              as.numeric(factor(zz, labels=seq_along(unique(zz))))
                                            } else cluster::pam(dist.mat2, k=g,  cluster.only=TRUE)[dup.ind], 
                                          hc=cutree(hcZ, k=g)[dup.ind]),  groups=seq_len(g))
        }
        if(isTRUE(noise))   {
          if(g0 > 1)        {
            zg0   <- mclust::unmap(switch(EXPR=init.z, 
                                          random=sample(seq_len(g0), size=N, replace=TRUE),
                                          kmedoids= if(do.wts) {
                                            zz <- WeightedCluster::wcKMedoids(dist.mat, k=g0, weights=weights, cluster.only=TRUE)
                                              as.numeric(factor(zz, labels=seq_along(unique(zz))))
                                            } else cluster::pam(dist.mat2, k=g0, cluster.only=TRUE)[dup.ind],  
                                          hc=cutree(hcZ, k=g0)[dup.ind]), groups=seq_len(g0))
            zg0   <- cbind(zg0 * (1 - tau0), tau0)
          } else   {
            zg0   <- matrix(tau0, nrow=N, ncol=2L)
          }
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
      if(isTRUE(verbose))      { cat(paste0("\n\tModel: ", modtype, "\n"))
        last.T    <- modtype  == T.last
      }
      ctrl$nmeth  <- is.element(modtype, n.meths)
      ctrl$equalNoise    <- equal.n0[ctrl$nmeth    + 1L,h]
      ctrl$equalPro      <- equal.tau[ctrl$nmeth   + 1L,h]
      ctrl$gate.g        <- gate.G[ctrl$nmeth      + 1L,h]
      ctrl$noise.gate    <- ifelse(ctrl$nmeth, noise.gate[2L,h], TRUE)
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
         if(isTRUE(verbose))     cat(paste0("\tRandom Start: #", i, "...\n"))
         switch(EXPR=algog,
               cemEM=          {
            ctrl$algo      <- "CEM"
            EMX[[i]]       <- .EM_algorithm(SEQ=SEQ, numseq=numseq, g=g, modtype=modtype, z=zm[[i]],  ctrl=ctrl, gating=gating, covars=covars)
            if(!EMX[[i]]$ERR)  {
              ctrl$algo    <- "EM"
              tmpEMX       <- EMX[[i]]
              j.i          <- pmax(tmpEMX$j, 2L)
              EMX[[i]]     <- .EM_algorithm(SEQ=SEQ, numseq=numseq, g=g, modtype=modtype, z=tmpEMX$z, ctrl=ctrl, gating=gating, covars=covars, ll=tmpEMX$ll[c(j.i - 1L, j.i)])
              if(EMX[[i]]$ERR) {
                EMX[[i]]   <- tmpEMX
                ctrl$algo  <- "CEM"
              }
            }
         }, EMX[[i]]       <- .EM_algorithm(SEQ=SEQ, numseq=numseq, g=g, modtype=modtype, z=zm[[i]],  ctrl=ctrl, gating=gating, covars=covars))
        }
        EMX       <- EMX[[which.max(vapply(lapply(EMX, "[[", "ll"), max, numeric(1L)))]]
      } else       {
        switch(EXPR=algog,
              cemEM=        {
          ctrl$algo        <- "CEM"
          EMX              <- .EM_algorithm(SEQ=SEQ, numseq=numseq, g=g, modtype=modtype, z=zm,       ctrl=ctrl, gating=gating, covars=covars)
          if(!EMX$ERR)      {
            ctrl$algo      <- "EM"
            tmpEMX         <- EMX
            j.i            <- pmax(EMX$j, 2L)
            EMX            <- .EM_algorithm(SEQ=SEQ, numseq=numseq, g=g, modtype=modtype, z=EMX$z,    ctrl=ctrl, gating=gating, covars=covars, ll=EMX$ll[c(j.i - 1L, j.i)])
            if(EMX$ERR)     {
              EMX          <- tmpEMX
              ctrl$algo    <- "CEM"
            }
          }
        }, EMX             <- .EM_algorithm(SEQ=SEQ, numseq=numseq, g=g, modtype=modtype, z=zm,       ctrl=ctrl, gating=gating, covars=covars))
      }
      ERR         <- EMX$ERR
      j           <- EMX$j
      Mstep       <- EMX$Mstep
      ll          <- EMX$ll
      z           <- EMX$z
      cvsel.X     <- cvsel && !ERR
      j2          <- max(1L, j - switch(EXPR=algog, cemEM=1L, 2L))
      if(isTRUE(verbose))        cat(paste0("\t\t# Iterations: ", ifelse(ERR, "stopped at ", ""), j2, ifelse(last.G && last.T, "\n\n", "\n")))

      if(all((Mstep$lambda -> lambda) == 0) && cvsel.X) {
        CVll      <- -N * P * log(V)
      } else if(cvsel.X)    {
        ctrl$warn <- FALSE
        lCV       <- vector("numeric", nfolds)
        zCV       <- z[cv.ind,,      drop=FALSE]
        gCV       <- covars[cv.ind,, drop=FALSE]
        for(i in foldseq)   {
          testX   <- which(cv.folds == i)
          CVS     <- cv.SEQ[testX]
          SCV     <- cv.SEQ[-testX]
          CVz     <- zCV[-testX,,    drop=FALSE]
          CVg     <- gCV[-testX,,    drop=FALSE]
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
          if(any(modtype %in% c("CU", "UU", "CUN", "UUN"), ctrl$opti == "mode", ctrl$ordering != "none")) {
            nCV            <- cv.numseq[,-testX, drop=FALSE]
            CVn            <- cv.numseq[,testX,  drop=FALSE]
            attr(nCV, "G") <- attr(numseq, "G")
            attr(nCV, "P") <- attr(numseq, "P")
          } else   {
            nCV   <-
            CVn   <- NULL
          }
          EMX     <- .EM_algorithm(SEQ=SCV, numseq=nCV, g=g, modtype=modtype, z=CVz, ctrl=ctrl, gating=gating, covars=CVg)
          MCV     <- EMX$Mstep
          if(is.matrix(MCV$tau)) {
            tau.tmp        <- stats::predict(MCV$fitG, newdata=gCV[testX,, drop=FALSE], type="probs")
            MCV$tau        <- if(ctrl$noise.gate) tau.tmp else cbind(tau.tmp * (1 - MCV$tau[1L,g]), MCV$tau[1L,g])
            rm(tau.tmp)
          }
          MCV$dG  <- NULL
          ctrl$do.cv       <- TRUE
          if(is.infinite(lCV[i] <- 
             .E_step(seqs=CVS, params=MCV, l.meth=modtype, ctrl=ctrl, 
             numseq=CVn)))       break
          ctrl$do.cv       <- FALSE
        }
        CVll      <- sum(lCV)
        ctrl$warn <- TRUE
      }
      
      z           <- if(do.uni) z[dis.agg,] else z
      log.lik     <- ll[j]
      nzero       <- sum(lambda == 0)
      ninfty      <- sum(is.infinite(lambda))
      Gfit        <- if(ctrl$gate.g) Mstep$fitG
      gate.pen    <- ifelse(ctrl$gate.g, length(stats::coef(Gfit)) + !ctrl$noise.gate, 
                     ifelse(ctrl$equalPro, as.integer(ctrl$nmeth && !ctrl$equalNoise), g - 1L))
      choice      <- .choice_crit(ll=log.lik, seqs=SEQ, z=z, l.meth=modtype, nonzero=ifelse(nonzero, sum(lambda != 0), NA), gate.pen=gate.pen)
      bicx        <- choice$bic
      iclx        <- choice$icl
      aicx        <- choice$aic
      dfx         <- choice$df
      if(do.dbs   && g > 1) {
        DBS       <- if(ctrl$do.wts) dbs(z, weights=w2, ...) else dbs(z, ...)
        dbsx      <- DBS$wmsw
        SILS[[h]][[m]]     <- if(ERR)      NA else DBS$silvals
        attr(SILS[[h]][[m]], "G")         <- g
        attr(SILS[[h]][[m]], "ModelType") <- modtype
      } else dbsx <- NA
      necx        <- ifelse(g > 1 && do.nec, -sum(apply(z, 1L, .entropy))/(log.lik - LL.x[1L,switch(EXPR=modtype, CC=, UC="CC", CU=, UU="CU", "CCN")]), NA)
      crit.t      <- switch(EXPR=criterion, cv=CVll, bic=bicx, icl=iclx, aic=aicx, nec=necx, dbs=dbsx)
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
      NECs[h,modtype]      <- ifelse(ERR,  Inf, -necx)
      DBSs[h,modtype]      <- ifelse(ERR, -Inf, dbsx)
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
    ZS[[h]]       <- stats::setNames(ZS[[h]],   modtypes)
    if(do.dbs     && g > 1)    {
      SILS[[h]]   <- stats::setNames(SILS[[h]], modtypes)  
    }
  } # for (g)

  seqs            <- seqX
  if(any(l.warn   <- x.ll  != cummax(x.ll))) {
    if(which.max(l.warn)   != 
       length(x.ll))             warning("Log-likelihoods are not strictly increasing\n", call.=FALSE)
  }
  if(any(IT.x[!is.na(IT.x)]
         == ctrl$itmax))         warning(paste0("One or more models failed to converge in the maximum number of allowed iterations (", itmax, ")\n"), call.=FALSE)
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
  CRITs           <- switch(EXPR=criterion, cv=CV.x, bic=BICs, icl=ICLs, aic=AICs, nec=NECs, dbs=DBSs)
  best.ind        <- which(CRITs == switch(EXPR=criterion, nec=-crit.gx, crit.gx), arr.ind=TRUE)
  if(nrow(best.ind) > 1)    {    warning(paste0("Ties for the optimal model exist according to the '", toupper(criterion), "' criterion: choosing the most parsimonious model\n"), call.=FALSE, immediate.=TRUE)
    best.ind      <- which(DF.x  == min(DF.x[best.ind]), arr.ind=TRUE)
    best.ind      <- best.ind[which.min(best.ind[,1L]),]
  }
  best.G          <- best.ind[1L]
  best.mod        <- colnames(CRITs)[best.ind[2L]]
  G               <- G[best.G]
  x.bic           <- BICs[best.ind]
  x.icl           <- ICLs[best.ind]
  x.aic           <- AICs[best.ind]
  x.nec           <- NECs[best.ind]
  x.dbs           <- DBSs[best.ind]
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
  
  noise           <- best.mod %in% c("CCN", "UCN", "CUN", "UUN")
  attr(x.lambda, "Nzero")       <- Nzero.x[best.ind]
  attr(x.lambda, "Ninfty")      <- Ninfty.x[best.ind]
  attr(DF.x,     "Nzero")       <- Nzero.x
  attr(DF.x,     "Ninfty")      <- Ninfty.x
  attr(DF.x,     "Gate.Pen")    <- x.gp
  if(any(apply(x.lambda == 0, 1L, all))) {
    x.theta                     <- if(G > 1) rbind(do.call(rbind, lapply(x.theta[-G], .char_to_num)), NA) else matrix(NaN, nrow=1L, ncol=P)
  } else x.theta                <- do.call(rbind, lapply(x.theta, .char_to_num))
  storage.mode(x.theta)         <- "integer"
  attr(x.theta, "alphabet")     <- levs
  attr(x.theta, "labels")       <- attr(seqs, "labels")
  attr(x.theta, "lambda")       <- switch(EXPR=best.mod, CCN=, CUN=rbind(matrix(x.lambda[1L,], nrow=G - 1L, ncol=P, byrow=best.mod == "CUN"), 0L), matrix(x.lambda, nrow=G, ncol=P, byrow=best.mod == "CU"))
  class(x.theta)                <- "MEDtheta"
  Gseq            <- seq_len(G)
  colnames(x.z)   <- if(G == 1  && noise) "Cluster0" else paste0("Cluster", if(noise) replace(Gseq, G, 0L) else Gseq)
  MAP             <- MAP2       <- max.col(x.z)
  MAP             <- if(noise) replace(MAP, MAP == G, 0L) else MAP
  equalPro        <- equal.tau[1L   + noise,best.G] 
  equalNoise      <- equal.n0[1L    + noise,best.G] 
  noise.gate      <- noise.gate[1L  + noise,best.G]
  noise.gate      <- ifelse(noise, noise.gate, TRUE)
  covars          <- if(do.uni) covars[dis.agg,] else covars
  rownames(covars)                 <- seq_along(MAP)
  if(!(gate.G[1L   + noise,best.G] -> bG))  {
    if(G > 1)      {
      if(ctrl$do.wts)            {
        z         <- x.z * w2
        z[apply(z == 0, 1L, all),] <- .Machine$double.eps
      } else z    <- x.z
      fitG        <- nnet::multinom(gating, trace=FALSE, data=covars, maxit=ctrl$g.itmax, reltol=ctrl$g.tol)
      if(equalPro && !equalNoise   && !noise) {
        tau0      <- mean(z[,G])
        x.tau     <- c(rep((1 - tau0)/(G - 1L), G  - 1L), tau0)
      } else       {
        x.tau     <- if(equalPro   || equalNoise) rep(1/G, G) else fitG$fitted.values[1L,]
      }
      x.tau       <- stats::setNames(x.tau, paste0("Cluster", if(noise) replace(Gseq, G, 0L) else Gseq))
    }   else       {
      fitG        <- suppressWarnings(stats::glm(z ~ 1, family=stats::binomial()))
    }
  }
  fitG$lab        <- if(noise.gate && G  > 1) replace(Gseq, G, 0L) else Gseq
  attr(fitG, "EqualNoise")      <- equalNoise
  attr(fitG, "EqualPro")        <- equalPro
  attr(fitG, "Formula")         <- Reduce(paste, deparse(gating[-2L]))
  attr(fitG, "Maxit")           <- ctrl$g.itmax
  attr(fitG, "Noise")           <- noise
  attr(fitG, "NoiseGate")       <- noise.gate
  attr(fitG, "Reltol")          <- ctrl$g.tol
  class(fitG)     <- c("MEDgating", class(fitG))
  if(isTRUE(verbose))            cat(paste0("\n\t\tBest Model", ifelse(length(CRITs) > 1, paste0(" (according to ", toupper(criterion), "): "), ": "), best.mod, ", with ",  paste0(G, " component", ifelse(G > 1, "s", "")),
                                     ifelse(bG | x.gcov, paste0(" (incl. ", ifelse(do.wts, "weights and ", ""), "gating network covariates)"), ifelse(do.wts, " (incl. weights)", "")), "\n\t\t",
                                     ifelse(cvsel, paste0("CV = ",  round(x.cv,  2L), " | "), ""),
                                     ifelse(G > 1 && do.nec, paste0("NEC = ", round(x.nec, 2L), " | "), ""),
                                     ifelse(G > 1 && do.dbs, paste0("DBS = ", round(x.dbs, 2L), " | "), ""),
                                     "BIC =", round(x.bic, 2L), " | ICL =", round(x.icl, 2L), " | AIC =", round(x.aic, 2L), "\n\n"))
  params          <- list(theta   = x.theta,
                          lambda  = x.lambda,
                          tau     = x.tau)
  results         <- list(call    = call,
                          data    = seqs,
                          modName = best.mod,
                          G       = G,
                          params  = params,
                          gating  = fitG,
                          z       = x.z,
                          MAP     = MAP,
                          ZS      = stats::setNames(ZS, rG))
  if(cvsel)        {
    results       <- c(results, list(CV = CV.x, cv = x.cv))
  }
  results         <- c(results, list(
                          BIC     = BICs,
                          bic     = x.bic,
                          ICL     = ICLs,
                          icl     = x.icl,
                          AIC     = AICs,
                          aic     = x.aic))
  if(do.dbs)       {
    DBSs          <- if(any(rG  == 1)) DBSs[-1L,, drop=FALSE]            else DBSs
    class(DBSs)   <- "MEDcriterion"
    attr(DBSs, "Criterion")     <- "DBS"
    attr(DBSs, "G")             <- if(any(rG  == 1)) rownames(BICs)[-1L] else rownames(BICs)
    attr(DBSs, "modelNames")    <- colnames(BICs)
    SILS          <- if(any(rG  == 1)) stats::setNames(SILS, rG)[-1L]    else stats::setNames(SILS, rG)
    attr(DBSs, "Weighted")      <- 
    attr(SILS, "Weighted")      <- do.wts
    results       <- c(results, list(DBS = DBSs, SILS = SILS))
    if(G > 1)      {
      x.sils      <- SILS[[as.character(G)]][[best.mod]]
      attr(x.dbs,  "Weighted")  <-
      attr(x.sils, "Weighted")  <- do.wts
      results     <- c(results, list(dbs = x.dbs, silvals = x.sils))
    }
  }
  if(do.nec)       {
    NECs          <- NECs[-1L,, drop=FALSE]
    class(NECs)   <- "MEDcriterion"
    attr(NECs, "Criterion")     <- "NEC"
    attr(NECs, "G")             <- rownames(BICs)[-1L]
    attr(NECs, "modelNames")    <- colnames(BICs)
    results       <- c(results, list(NEC = NECs))
    results       <- if(G > 1) c(results, list(nec = x.nec)) else results
  }
  x.ll            <- x.ll[if(G > 1) switch(EXPR=algo, cemEM=-1L, -seq_len(2L)) else 1L]
  attr(x.ll, "Weighted")        <- do.wts
  results         <- c(results, list(
                          LOGLIK  = LL.x,
                          loglik  = x.ll,
                          uncert  = if(G > 1) 1 - matrixStats::rowMaxs(x.z) else vector("integer", N),
                          covars  = covars,
                          DF      = DF.x,
                          df      = DF.x[best.ind],
                          ITERS   = IT.x,
                          iters   = IT.x[best.ind]))
  dist.mat        <- as.matrix(dist.mat2)
  N               <- length(MAP)
  if(!all((unip   <- unique(MAP)) == 0)) {
    Nseq          <- seq_len(N)
    meds          <- NULL
    set.seed(100)
    for(g in sort(unip[unip > 0]))  {
      srows       <- Nseq[MAP2 == g]
      meds        <- c(meds, srows[which.min(matrixStats::rowSums2(dist.mat[srows,srows]))])
    }
    perm          <- seriation::get_order(seriation::seriate(as.dist(dist.mat[meds,meds]), method="TSP"))
  } else perm     <- NULL
  attr(results, "Algo")         <- algo
  attr(results, "Counts")       <- if(do.uni) c2         else rep(1L, N)
  attr(results, "Criterion")    <- criterion
  attr(results, "CV")           <- cvsel
  attr(results, "DBS")          <- do.dbs
  attr(results, "DistMat")      <- dist.mat
  attr(results, "EqualPro")     <- x.ctrl$equalPro
  attr(results, "EqualNoise")   <- x.ctrl$equalNoise
  attr(results, "Gating")       <- x.gcov | bG
  attr(results, "N")            <- N
  attr(results, "NEC")          <- do.nec
  attr(results, "Noise")        <- noise
  attr(results, "NoiseGate")    <- x.ctrl$noise.gate
  attr(results, "P")            <- P
  attr(results, "Seriate")      <- if(noise)  c(perm, G) else perm
  attr(results, "Unique")       <- do.uni
  attr(results, "V")            <- V
  attr(results, "Weighted")     <- do.wts
  attr(results, "Weights")      <- if(do.wts) weights    else rep(1L, N)
  class(results)  <- "MEDseq"
    return(results)
}

plot.MEDseq       <- function(x, type = c("clusters", "mean", "lambda", "gating", "cv", "bic", "icl", "aic", "nec", "dbs", "LOGLIK", "silhouette", "uncert.bar", 
                              "uncert.profile", "loglik", "d", "f", "Ht", "i", "I"), seriate = TRUE, preczero = TRUE, log.scale = NULL, ...) {
  x               <- if(inherits(x, "MEDseqCompare")) x$optimal else x
  if(!missing(type)           &&
     (length(type)       > 1  ||
      !is.character(type)))      stop("'type' must be a character vector of length 1", call.=FALSE)
  if(length(seriate)     > 1  ||
     !is.logical(seriate))       stop("'seriate' must be a single logical indicator",  call.=FALSE)
  type            <- match.arg(type)
  savepar         <- par(no.readonly=TRUE)
  on.exit(par(savepar))
  G               <- x$G
  N               <- attr(x, "N")
  P               <- attr(x, "P")
  V       <- mV   <- attr(x, "V")
  noise           <- attr(x, "Noise")
  dmat            <- attr(x, "DistMat")
  weighted        <- attr(x, "Weighted")
  dat             <- x$data
  modName         <- x$modName
  Gseq    <- perm <- seq_len(G)
  Nseq            <- seq_len(N)
  Pseq            <- seq_len(P)
  symbols         <- c(17, 2, 16, 10, 13, 18, 15, 7)
  use.col         <- c("gray", "dodgerblue1", "red3", "slateblue", "green3", "violet", "gold", "hotpink")
  alpha.x         <- attr(dat, "alphabet")
  cpal.x          <- attr(dat, "cpal")
  label.x <- lab  <- attr(dat, "labels")
  if(is.element(type,
     c("mean", "lambda")))     {
    theta         <- x$params$theta
    lambda        <- attr(theta, "lambda")
    class(theta)  <- NULL
  }
  dots            <- list(...)
  dots            <- dots[unique(names(dots))]
  switch(EXPR=type, silhouette={
    has.dot       <- length(dots) > 0
  }, {
    if((has.dot   <- length(dots) > 0 && any(names(dots)  == "what")))  {
      MAP         <- dots$MAP
      if(length(MAP)    != N  ||
         !is.numeric(MAP)     ||
        MAP       != floor(MAP)) stop(paste0("'MAP' must be an integer vector of length N=", N),  call.=FALSE)
    }
  })
  if(is.element(type, c("clusters", "d", "f", "Ht", "i", "I"))   ||
     all(isTRUE(seriate), is.element(type, c("mean", "lambda")))) {
    switch(EXPR=type,
           mean=,
           lambda= {
      perm        <- attr(x, "Seriate")
      switch(EXPR=type,
             mean= {
        theta     <- theta[perm,,  drop=FALSE]
      },   lambda= {
        lambda    <- lambda[perm,, drop=FALSE]
      })
    },             {
      if(length(dots) > 0)     {
        if(has.dot)   {
          G       <- length(unique(MAP))
          noise   <- any(MAP  == 0)
        } else     {
          MAP     <- do.call(get_results, c(list(x=x, what="MAP"), dots[!(names(dots) %in% c("x", "what"))]))
          G       <- attr(MAP, "G")
          noise   <- attr(MAP, "Noise")
        }
        Gseq      <-
        perm      <- seq_len(G)
        if(isTRUE(seriate))    {
          unip    <- unique(MAP)
          meds    <- NULL
          set.seed(100)
          for(g in sort(unip[unip > 0])) {
            srows <- Nseq[MAP == g]
            meds  <- c(meds, srows[which.min(matrixStats::rowSums2(dmat[srows,srows]))])
          }
          perm    <- seriation::get_order(seriation::seriate(as.dist(dmat[meds,meds]), method="TSP"))
        }
      } else       {
        MAP       <- x$MAP
        perm      <- if(isTRUE(seriate)) attr(x, "Seriate") else perm
      }
      perm        <- if(isTRUE(noise)) replace(perm, G, 0L) else perm
    })
  }
  perm            <- if(noise && type != "clusters") replace(perm, G, "Noise") else perm

  switch(EXPR=type,
         clusters=     {
    glo.order     <-
    num.cl        <- NULL
    set.seed(200)
    for(g in Gseq) {
      srows       <- Nseq[MAP == perm[g]]
      num.cl      <- c(num.cl, length(srows))
      glo.order   <- c(glo.order, srows[seriation::get_order(seriation::seriate(as.dist(dmat[srows,srows]), method="TSP"))])
    }
    cum.cl        <- cumsum(num.cl)

    layout(rbind(1, 2), heights=c(0.85, 0.15), widths=1)
    OrderedStates <- data.matrix(.fac_to_num(dat))[glo.order,]
    image(x=Pseq, z=t(OrderedStates), y=Nseq, zlim=c(1, length(cpal.x)), xlim=c(1, P), ylim=c(1, N),
          axes=FALSE, xlab="", ylab="Clusters", col=cpal.x, main=ifelse(seriate, "Observations Ordered by Cluster", "Clusters"))
    box(lwd=2)
    axis(side=1, at=Pseq, labels=attr(x$data, "names"), cex.axis=0.75)
    gcl           <- c(0, cum.cl)
    axis(side=2, at=gcl[-length(gcl)] + diff(gcl)/2, labels=if(noise) replace(perm, G, "Noise") else perm, lwd=1, line=-0.5, las=2, tick=FALSE, cex.axis=0.75)
    for(g in Gseq) {
      segments(0, cum.cl[g], P + 0.5, cum.cl[g], lwd=2)
    }

    par(mar=c(1, 1, 0.5, 1) + 0.1, xpd=FALSE)
    plot.new()
    legend("bottom", fill=cpal.x, legend=label.x, ncol=ceiling(V/ifelse(V > 6, 3, 2)), cex=0.75)
    graphics::layout(1)
      invisible(if(isTRUE(noise)) replace(perm, G, G) else perm)
  }, mean=         {
    miss.prec     <- missing(preczero)
    if(length(preczero)  > 1  ||
       !is.logical(preczero))    stop("'preczero' must be a single logical indicator", call.=FALSE)
    if(indmiss    <- preczero) {
      lmiss       <- lambda   == 0
      if(indmiss  <-
         any(lmiss))     {
        lab       <- c(label.x, expression(paste(lambda, " = 0")))
        mV        <- V   + 1L
        if(any(gmiss    <- apply(lmiss, 1L, all))) {
          if(G    == 1L) {       message("The single central sequence is entirely missing\n")
          } else                 message(paste0("One or more central sequences (", paste(shQuote(which(gmiss)), collapse=" & "), ") are entirely missing\n"))
        }
      }     else if(!miss.prec)  message("No missing values to discard\n")
      missind           <- which(tabulate(theta[!lmiss], nbins=V) == 0)
      theta[lmiss]      <- NA
    } else missind      <- which(tabulate(theta,         nbins=V) == 0)
    l.ncol        <- ceiling(mV/ifelse(mV > 6, 3, 2))
    if(vmiss      <-
       any(missind))             message(paste0("One or more sequence categories (", paste(shQuote(label.x[missind]), collapse=" & "), ") are entirely missing\n"))
    dat           <- suppressMessages(seqdef(as.data.frame(theta), states=alpha.x, labels=label.x, cpal=if(vmiss) c(cpal.x[-missind], rep(NA, length(missind))) else cpal.x))
    attr(dat, "names")  <- attr(x$data, "names")
    layout(rbind(1, 2), heights=c(0.85, 0.15), widths=1)
    seqIplot(dat, with.legend=FALSE, main="Central Sequences Plot", border=NA, missing.color=par()$bg, ylab=paste0(G, " seq. (N=", N, ")"), yaxis=FALSE, cex.axis=0.75)
    if(G > 1) graphics::axis(2, at=seq_len(G) - 0.5, labels=as.character(perm), tick=FALSE, las=2, line=-0.5, cex.axis=0.75)
    par(mar=c(1, 1, 0.5, 1) + 0.1, xpd=FALSE)
    plot.new()
    legend("bottom",          fill=if(indmiss) c(cpal.x, par()$bg) else cpal.x, legend=lab, ncol=l.ncol, cex=0.75)
    graphics::layout(1)
      invisible(theta)
  }, lambda=       {
    log.scale     <- ifelse(missing(log.scale), is.element(modName, c("UU", "UUN")), log.scale)
    if(length(log.scale) > 1  ||
       !is.logical(log.scale))   stop("'log.scale' must be a single logical indicator", call.=FALSE)
    layout(rbind(1, 2), heights=c(0.85, 0.15))
    par(mar=c(4.1, 4.1, 4.1, 3.1))
    i.ind         <- is.infinite(lambda)
    num.ind       <- !i.ind  &  lambda >  0
    dat           <- if(log.scale) log(lambda[num.ind]) else lambda[num.ind]
    facs          <- if(length(dat) > 1) cut(dat, 30L, include.lowest=TRUE) else 1L
    cmat          <- matrix("", nrow=G, ncol=P)
    cols          <- rev(heat.colors(30L))
    cmat[i.ind]            <- "grey50"
    cmat[lambda   == 0]    <- "green3"
    cmat[num.ind]          <- cols[as.numeric(facs)]
    levels        <- sort(unique(as.vector(cmat)))
    z             <- matrix(unclass(factor(cmat, levels=levels, labels=seq_along(levels))), nrow=P, ncol=G, byrow=TRUE)
    Gseq          <- seq_len(G)
    image(Pseq, Gseq, z, col=levels, axes=FALSE, xlab="Positions", ylab="Clusters", main=paste0("Precision Parameters Plot", ifelse(log.scale, " (Log Scale)",  "")))
    axis(1, at=Pseq, tick=FALSE, las=1, cex.axis=0.75, labels=Pseq)
    axis(2, at=Gseq, tick=FALSE, las=1, cex.axis=0.75, labels=as.character(perm))
    box(lwd=2)
    bx            <- par("usr")
    xpd           <- par()$xpd
    box.cx        <- c(bx[2L] + (bx[2L]  - bx[1L])/1000, bx[2L] + (bx[2L] - bx[1L])/1000 + (bx[2L] - bx[1L])/50)
    box.cy        <- c(bx[3L],   bx[3L])
    box.sy        <- (bx[4L]  -  bx[3L]) / length(cols)
    xx            <- rep(box.cx, each=2L)
    par(xpd = TRUE)
    for(i in seq_along(cols)) {
      yy          <- c(box.cy[1L] + (box.sy * (i - 1)),
                       box.cy[1L] + (box.sy * (i)),
                       box.cy[1L] + (box.sy * (i)),
                       box.cy[1L] + (box.sy * (i - 1)))
      polygon(xx, yy, col = cols[i], border = cols[i])
    }
    par(new=TRUE)
    plot(0, 0, type="n", ylim=if(length(dat) > 1) range(dat, na.rm=TRUE) else c(0, 1), yaxt="n", ylab="", xaxt="n", xlab="", frame.plot=FALSE)
    axis(side=4, las=2, tick=FALSE, line=0.5)
    par(mar=c(0, 0, 0, 0))
    plot.new()
    legend("center", c(expression(paste(lambda, " = 0")), expression(paste(lambda %->% infinity))), fill=c("green3", "grey50"), ncol=2, text.width=0.1, cex=1.25)
    graphics::layout(1)
      invisible()
  }, gating=       {
    suppressWarnings(par(pty="m"))
    Tau        <- .mat_byrow(x$params$tau, nrow=N, ncol=ncol(x$z))
    miss.x     <- length(dots) > 0 && any(names(dots) == "x.axis")
    x.axis     <- if(miss.x) dots$x.axis else seq_len(N)
    xlab       <- ifelse(miss.x, "", "Observation")
    type       <- ifelse(miss.x, "p", "l")
    col        <- if(noise)  replace(seq_len(G), G, "grey65") else seq_len(G)
    if(length(x.axis) != N)      stop("'x.axis' must be of length N", call.=FALSE)
    if(x.fac   <- is.factor(x.axis)) {
      xlev     <- levels(x.axis)
      x.axis   <- as.integer(x.axis)
      xaxt     <- "n"
    } else      {
      xaxt     <- "s"
    }
    graphics::matplot(x=x.axis, y=Tau, type=type, main="Gating Network", xaxt=xaxt, xlab=xlab, ylab="", col=col, pch=1)
    mtext(expression(widehat(tau)[g]), side=2, las=2, line=3)
    if(x.fac)   {
      graphics::axis(1, at=unique(x.axis), labels=xlev)
    }
  }, cv=,
     bic=,
     icl=,
     aic=,
     nec=,
     dbs=,
     LOGLIK=       {
    if(all(type   == "cv",
       !attr(x, "CV")))          stop("Cross-validated log-likelihood values cannot be plotted as cross-validation didn't take place during model fitting\n", call.=FALSE)
    dat           <- switch(EXPR=type, cv=x$CV, bic=x$BIC, icl=x$ICL, aic=x$AIC, nec=x$NEC, dbs=x$DBS, LOGLIK=x$LOGLIK)
    if(type ==  "nec"     &&
      (!attr(x, "NEC")     || 
       all(type   == "nec",
       is.null(dat))))           stop("NEC values cannot be plotted", call.=FALSE)
    if(type ==  "dbs"      &&
      (!attr(x, "DBS")     ||
       all(type   == "dbs",
       is.null(dat))))           stop("DBS values cannot be plotted", call.=FALSE)
    ms            <- which(c("CC", "UC", "CU", "UU", "CCN", "UCN", "CUN", "UUN") %in% colnames(dat))
    symbols       <- symbols[ms]
    use.col       <- use.col[ms]
    matplot(dat, type="b", xlab="Number of Components", ylab=switch(EXPR=type, cv=, LOGLIK=, dbs="", toupper(type)), col=use.col, pch=symbols, ylim=range(as.vector(dat[!is.na(dat) & is.finite(dat)])), xaxt="n", lty=1)
    if(type == "cv")     mtext(ifelse(weighted, expression("\u2113"["cv"]^"w"), expression("\u2113"["cv"])), side=2, line=3, las=1, cex=1.5)
    if(type == "LOGLIK") mtext(paste0(ifelse(weighted, "Weighted ", ""), "Log-Likelihood"), side=2, line=3, las=3)
    if(type == "dbs")    mtext(paste0(ifelse(weighted, "Weighted ", ""), "Median DBS"), side=2, line=3, las=3)
    axis(1, at=seq_len(nrow(dat)), labels=rownames(dat))
    legend(switch(EXPR=type, nec=, dbs="topright", "bottomright"), ncol=2, cex=1, inset=0.01, legend=colnames(dat), pch=symbols, col=use.col)
      invisible()
  }, silhouette=   {
    if(!attr(x, "DBS")     ||
       all(type   == "dbs",
       is.null(x$SILS)))         stop("DBS values cannot be plotted as only 1-component models were fitted", call.=FALSE)
    object        <- if(has.dot) do.call(get_results, c(list(x=x, what="sils"), dots[!(names(dots) %in% c("x", "what"))])) else x$silvals
    rownames(object)       <- as.character(Nseq)
    cl            <- object[,"cluster"]
    X             <- object[order(cl, -object[,"dbs_width"]),, drop=FALSE]
    dbs           <- X[,"dbs_width"]
    space         <- c(0L, rev(diff(cli <- X[,"cluster"])))
    space[space   != 0]    <- 0.5
    ng            <- table(cl)
    G             <- attr(object, "G")
    dat           <- rev(barplot(rev(dbs), space=space, xlim=c(min(0, min(dbs)), 1), horiz=TRUE, las=1, mgp=c(2.5, 1, 0), col="gray", border=0, cex.names=par("cex.axis"), axisnames=FALSE))
    if(weighted)   {
      weights     <- attr(x, "Weights")
      title(main="(Weighted) Density-based Silhouette Plot", sub=paste0("(Weighted) Median DBS Width : ", round(matrixStats::weightedMedian(dbs, weights), digits=3)), adj=0)
    } else title(main="Density-based Silhouette Plot", sub=paste0("Median DBS Width : ", round(median(dbs), digits=3)), adj=0)
    mtext(paste0("n = ", N),  adj=0)
    modtype       <- attr(object, "ModelType")
    noise         <- modtype %in% c("CCN", "UCN", "CUN", "UUN")
    mtext(substitute(G~modtype~"clusters"~C[g], list(G=G, modtype=modtype)), adj=1)
    if(weighted)   {
      mtext(expression(paste(g, " : ", n[g], " | ", med[i %in% Cg] ~ ~w[i]~s[i])), adj=1.05, line=-1.2)
      dbcl        <- split(dbs,  cli)
      meds        <- sapply(seq_along(dbcl), function(i) matrixStats::weightedMedian(dbcl[[i]], weights[cli == i]))
    } else         {
      mtext(expression(paste(g, " : ", n[g], " | ", med[i %in% Cg] ~ ~s[i])), adj=1.05, line=-1.2)
      meds        <- tapply(dbs, cli, median)
    }
    medy          <- tapply(dat, cli, median)
    for(g in seq_len(G)) {
      text(1, medy[g], paste(ifelse(g == G && noise, 0L, g), ":  ", ng[g], " | ", format(meds[g], digits=1, nsmall=2)), xpd=NA, adj=0.8)
    }
      invisible()
  }, uncert.profile=,
     uncert.bar=   {
    par(pty="m", mar=c(5.1, 4.1, 4.1, 3.1))
    uncX          <- x$uncert
    oneG          <- 1/G
    min1G         <- 1 - oneG
    yx            <- unique(c(0, pretty(c(0, min1G))))
    yx            <- replace(yx, length(yx), min1G)
    cm            <- c("dodgerblue2", "red3", "green3")
    switch(EXPR=type,
           uncert.bar=         {
             cu   <- cm[seq_len(2L)][(uncX >= oneG) + 1L]
             cu[uncX == 0] <- NA
             plot(uncX, type="h", ylim=range(yx), col=cu, yaxt="n", ylab="", xlab="Observations", lend=1)
             lines(x=c(0, N), y=c(oneG, oneG), lty=2, col=cm[3L])
             axis(2, at=yx, labels=replace(yx, length(yx), "1 - 1/G"), las=2, xpd=TRUE)
             axis(2, at=oneG, labels="1/G", las=2, xpd=TRUE, side=4, xpd=TRUE)
           }, uncert.profile=  {
             ord  <- order(uncX, decreasing=FALSE)
             ucO  <- uncX[ord]
             plot(ucO, type="n", ylim=c(-max(uncX)/32, max(yx)), ylab="", xaxt="n", yaxt="n", xlab=paste0("Observations in order of increasing uncertainty"))
             lines(x=c(0, N), y=c(0, 0), lty=3)
             lines(ucO)
             points(ucO, pch=15, cex=0.5, col=1)
             lines(x=c(0, N), y=c(oneG, oneG), lty=2, col=cm[3L])
             axis(2, at=yx,   las=2, xpd=TRUE, labels=replace(yx, length(yx), "1 - 1/G"))
             axis(2, at=oneG, las=2, xpd=TRUE, labels="1/G", side=4)
           })
    mtext("Uncertainty", side=2, line=3)
    title(main=list(paste0("Clustering Uncertainty ", switch(EXPR=type, uncert.bar="Barplot", "Profile Plot"))))
      invisible(uncX)
  }, loglik=       {
    x             <- x$loglik
    if(all(x      != cummax(x))) warning("Log-likelihoods are not strictly increasing\n", call.=FALSE)
    plot(x, type=ifelse(length(x) == 1, "p", "l"), xlab="Iterations", ylab=paste0(ifelse(weighted, "Weighted ", ""), "Log-Likelihood"), xaxt="n")
    seqll         <- seq_along(x)
    llseq         <- pretty(seqll)
    llseq         <- if(any(llseq != floor(llseq))) seqll else llseq
    axis(1, at=llseq, labels=llseq)
      invisible(list(ll = x, converge = x[length(x)]))
  },               {
    MAP           <- factor(replace(MAP, MAP == 0, "Noise"), levels=perm)
    dots          <- c(list(seqdata=dat, with.legend=FALSE, group=MAP, type=type, with.missing=FALSE), dots[!(names(dots) %in% c("G", "modtype", "noise"))])
    dots          <- if(type == "i") dots else c(list(border=NA), dots)
    dots          <- dots[unique(names(dots))]
    if(type       != "Ht")     {
      l.ncol      <- ceiling(V/ifelse(V > 6 | G %% 2 != 0, 3, 2))
      if(G  %% 2  == 0)        {
        par(oma=c(7.5, 0, 0, 0),       xpd=TRUE)
        suppressWarnings(do.call(seqplot, dots))
        par(fig=c(0, 1, 0, 1), oma=c(0.5, 0, 0, 0), mar=c(0.5, 0, 0, 0), new=TRUE)
        plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
        legend("bottom",      fill=attr(dat, "cpal"), legend=lab, ncol=l.ncol)
      } else       {
        suppressWarnings(do.call(seqplot, dots))
        par(mar=c(1, 1, 0.5, 1) + 0.1, xpd=TRUE)
        legend("bottomright", fill=attr(dat, "cpal"), legend=lab, ncol=l.ncol)
      }
    } else         {
      try(suppressWarnings(do.call(seqplot, dots)),   silent=TRUE)
    }
      invisible()
  })
}

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
  weighted        <- attr(x, "Weighted")
  attributes(x)   <- NULL
  attr(x, "dim")         <- dim1
  attr(x, "dimnames")    <- dim2
  cat(switch(EXPR=crit,
             NEC="Normalised Entropy Criterion (NEC):\n",
             DBS=paste0(ifelse(weighted, "(Weighted) ", ""), "Median Density-based Silhouette (DBS):\n"),
             CV=paste0("Cross-Validated ", ifelse(weighted, "(Weighted) ", ""), "Log-Likelihood (CV):\n"),
             BIC="Bayesian Information Criterion (BIC):\n",
             ICL="Integrated Completed Likelihood (ICL):\n",
             AIC="Akaike Information Criterion (AIC):\n",
              DF="Number of Estimated Parameters (Residual DF):\n",
           ITERS=paste0("Number of ", algo, " Iterations:\n"),
          loglik=paste0("Maximal ", ifelse(weighted, "(Weighted) ", ""), "Log-Likelihood:\n")))
  print(unclass(x))
  cat(paste0("\nTop ", ifelse(pick > 1, paste0(pick, " models"), "model"), " based on the ", toupper(crit), " criterion:\n"))
    print(choice$crits)
}

print.MEDgating    <- function(x, ...) {
  noise            <- attr(x, "Noise")
  equalpro         <- attr(x, "EqualPro")
  formula          <- attr(x, "Formula")
  equalNoise       <- noise && equalpro
  gateNoise        <- noise && !equalpro && formula != "~1"
  class(x)         <- class(x)[class(x) != "MEDgating"]
  print(x, ...)
  cat(paste("Formula:", formula, "\n"))
  cat(paste("Noise:",   noise,   "\n"))
  if(gateNoise)                  cat(paste("Noise Component Gating:", attr(x, "NoiseGate"), "\n"))
  cat(paste("EqualPro:", equalpro, ifelse(equalNoise, "\n", "")))
  if(equalNoise)                 cat(paste("Noise Proportion Estimated:", attr(x, "EqualNoise")))
  message("\n\nUsers are cautioned against making inferences about statistical significance from summaries of the coefficients in the gating network\n")
    invisible(x)
}

print.MEDseq      <- function(x, digits = 2L, ...) {
  cat("Call:\t");  print(x$call)
  if(length(digits)  > 1    || !is.numeric(digits) ||
     digits       <= 0)          stop("Invalid 'digits'", call.=FALSE)
  name            <- x$modName
  G               <- x$G
  noise           <- attr(x, "Noise")
  gating          <- attr(x$gating, "Formula")
  gate.x          <- !attr(x, "Gating")
  equalP          <- G == 1 || attr(x$gating, "EqualPro")
  equalN          <- noise  && attr(x$gating, "EqualNoise") && equalP
  no.dbs          <- G == 1 || !attr(x, "DBS")
  no.nec          <- G == 1 || !attr(x, "NEC")
  crit            <- round(unname(c(x$bic, x$icl, x$aic)), digits)
  crit            <- if(no.dbs) crit else c(crit, round(unname(x$dbs), digits))
  crit            <- if(no.nec) crit else c(crit, round(unname(x$nec), digits))
  cat(paste0("\nBest Model", ifelse(length(x$BIC)  > 1, paste0(" (according to ", toupper(attr(x, "Criterion")), "): "), ": "), name, ", with ",
             G, " component",      ifelse(G > 1, "s ", " "),
             paste0("and ", ifelse(attr(x, "Weighted"), "", "no "), "weights"),
             ifelse(gate.x,        " (no covariates)\n", " (incl. gating network covariates)\n"),
             ifelse(!equalP ||
                     G == 1, "",   paste0("Equal Mixing Proportions", ifelse(equalN | G == 1 | !noise, "\n", " (with estimated noise component mixing proportion)\n"))),
             ifelse(attr(x, "CV"), paste0("CV = ", round(unname(x$cv), digits), " | "), ""),
             ifelse(no.nec  && 
                    no.dbs,        paste0(
             "BIC = ",             crit[1L]), 
             ifelse(no.nec,        paste0(
             "DBS = ",             crit[4L],
             " | BIC = ",          crit[1L]),
             ifelse(no.dbs,        paste0(
             "NEC = ",             crit[4L],
             " | BIC = ",          crit[1L]),
             paste0("DBS = ",      crit[5L],
             " | NEC = ",          crit[4L],
             " | BIC = ",          crit[1L])))),
             " | ICL = ",          crit[2L],
             " | AIC = ",          crit[3L],
             ifelse(gate.x,  "",   paste0("\nGating: ", gating, "\n"))))
    invisible()
}

print.MEDseqCompare    <- function(x, index=seq_len(x$pick), digits = 3L, ...) {
  index           <- if(is.logical(index)) which(index) else index
  if(length(index) < 1 || (!is.numeric(index) &&
     (any(index    < 1  | 
          index    > x$pick))))  stop("Invalid 'index'",  call.=FALSE)
  if(length(digits)     > 1 ||
     !is.numeric(digits)    ||
     digits       <= 0)          stop("Invalid 'digits'", call.=FALSE)
  crit            <- attr(x, "Crit")
  opt             <- attr(x, "Opt")
  x$bic           <- round(x$bic,    digits)
  x$icl           <- round(x$icl,    digits)
  x$aic           <- round(x$aic,    digits)
  x$loglik        <- round(x$loglik, digits)
  na.nec          <- is.na(x$nec)
  x$nec[!na.nec]  <- if(!all(na.nec)) round(x$nec[!na.nec], digits)
  nec             <- replace(x$nec, na.nec, "")
  x$nec           <- NULL
  x$nec           <- if(all(na.nec))             NULL else nec
  na.dbs          <- is.na(x$dbs)
  x$dbs[!na.dbs]  <- if(!all(na.dbs)) round(x$dbs[!na.dbs], digits)
  dbs             <- replace(x$dbs, na.dbs, "")
  x$dbs           <- NULL
  x$dbs           <- if(all(na.dbs))             NULL else dbs
  na.cvs          <- is.na(x$cv)
  x$cv[!na.cvs]   <- if(!all(na.cvs)) round(x$cv[!na.cvs],  digits)
  cvs             <- replace(x$cv, na.cvs,  "")
  x$cv            <- NULL
  x$cv            <- if(all(na.cvs))             NULL else cvs
  x               <- c(x[seq_len(which(names(x) == "loglik") - 1L)], list(nec=x$nec, cvs=x$cv, dbs=x$dbs), x[seq(from=which(names(x) == "loglik"), to=length(x), by=1L)])
  x               <- x[unique(names(x))]
  n.all           <- all(x$noise)
  x$noise         <- if(n.all)                   NULL else x$noise
  noise.gate      <- if(n.all)                   NULL else replace(x$noise.gate, is.na(x$noise.gate), "")
  x$noise.gate    <- NULL
  x$noise.gate    <- if(all(x$gating == "None")) NULL else noise.gate
  equalPro        <- if(all(is.na(x$equalPro)))  NULL else replace(x$equalPro,   is.na(x$equalPro),   "")
  x$equalPro      <- NULL
  x$equalPro      <- equalPro
  na.equalNoise   <- is.na(x$equalNoise)
  equalNoise      <- replace(x$equalNoise, na.equalNoise,    "")
  x$equalNoise    <- NULL
  x$equalNoise    <- if(all(na.equalNoise))      NULL else equalNoise
  cat(paste0("------------------------------------------------------------------------------\n", 
             x$title, "\nData: ", x$data, "\nRanking Criterion: ", toupper(crit), "\nOptimal Only: ", opt,
             "\n------------------------------------------------------------------------------\n\n"))
  compX           <- data.frame(do.call(cbind, x[-seq_len(4L)]))[index,]
  compX           <- cbind(rank = rownames(compX), compX)
  rownames(compX) <- NULL
  print(compX, row.names = FALSE)
    invisible()
}

print.MEDtheta    <- function(x, preczero = TRUE, ...) {
  lambda          <- attr(x, "lambda")
  alpha           <- attr(x, "alphabet")
  lab.x           <- attr(x, "labels")
  G               <- nrow(x)
  V               <- length(alpha)
  class(x)        <- NULL
  miss.prec       <- missing(preczero)
  if(length(preczero)  > 1  ||
     !is.logical(preczero))      stop("'preczero' must be a single logical indicator", call.=FALSE)
  if(preczero)     {
    if(any(lam0   <- lambda == 0)) {
      x[lam0]     <- V + 1L
      alpha       <- c(alpha, "*")
      if(any(gm0  <- apply(lam0, 1L, all))) {
        if(G      == 1L)     {   message("The single central sequence is entirely missing\n")
        } else                   message(paste0("One or more central sequences (", paste(shQuote(which(gm0)), collapse=" & "), ") are entirely missing\n"))
      }
    }     else if(!miss.prec)    message("No missing values to discard\n")
    missind       <- which(tabulate(x[!lam0], nbins=V) == 0)
  } else missind  <- which(tabulate(x,        nbins=V) == 0)
  if(any(missind))               message(paste0("One or more sequence categories (", paste(shQuote(lab.x[missind]), collapse=" & "), ") are entirely missing\n"))
  theta           <- as.data.frame(lapply(as.data.frame(x), function(theta) .replace_levels(.num_to_char(theta), alpha)))
    print(theta, ...)
}

tabcluster        <- function(x, norm = FALSE) {
    UseMethod("tabcluster")
}

tabcluster.MEDseq <- function(x, norm = FALSE) {
  x               <- if(inherits(x, "MEDseqCompare")) x$optimal else x
  alph            <- attr(x$params$theta, "alphabet")
  V               <- attr(x, "V")
  P               <- ifelse(isTRUE(norm), attr(x, "P"), 1L)
  G               <- x$G
  noise           <- attr(x, "Noise")
  gnames          <- paste0("Cluster", seq_len(G))
  gnames          <- if(noise) replace(gnames, G, "Noise")   else gnames
  MAP             <- if(noise) replace(x$MAP, x$MAP == 0, G) else x$MAP
  temp            <- do.call(rbind, by(x$data, MAP, function(x) tabulate(do.call(base::c, x), V)))
  tabMAP          <- tabulate(MAP)
  temp            <- temp/(tabMAP * P) * ifelse(isTRUE(norm), 100L, 1L)
  temp            <- cbind(tabMAP, temp)
  rownames(temp)  <- gnames
  colnames(temp)  <- c("Size", alph)
    temp
}
#