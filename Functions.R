aitken            <- function(loglik) {
  if(!is.numeric(loglik) ||
     length(loglik)      != 3)   stop("'loglik' must be a numeric vector of length 3", call.=FALSE)
  l1              <- loglik[1L]
  l2              <- loglik[2L]
  l3              <- loglik[3L]
  if(any(is.infinite(loglik))) {
    linf          <- Inf
    a             <- NA
  } else {
    a             <- ifelse(l2 > l1, (l3 - l2) / (l2 - l1),    0L)
    denom         <- max(1L - a, .Machine$double.eps)
    linf          <- ifelse(a  < 1L,  l2 + (l3 - l2) / denom, Inf)
  }
    return(list(ll = l3, linf  = linf, a = a))
}

char_to_num       <- function(x) {
    as.numeric(strsplit(x, "")[[1L]])
}

choice_crit       <- function(ll, seqs, z, l.meth = c("CC", "UC", "CU", "UU", "CCN", "UCN", "CUN", "UUN"), nonzero = NULL, equalPro = FALSE) {
  G               <- ifelse(noise <- is.element(l.meth, c("CCN", "UCN", "CUN", "UUN")), attr(seqs, "G") - 1L, attr(seqs, "G"))
  P               <- attr(seqs, "P")
  kpar            <- ifelse(is.null(nonzero),
                            switch(EXPR= l.meth,
                                   CC  = G * P   + 1L,
                                   CCN = G * P   + (G > 0),
                                   UC  =,
                                   UCN = G * (P  + 1L),
                                   CU  =,
                                   CUN = P * (G  + 1L),
                                   UU  =,
                                   UUN = P * G   * 2L),
                            switch(EXPR= l.meth,
                                   CC  = nonzero * G  * P   + 1L,
                                   CCN = nonzero * G  * P   + (G > 0),
                                   UC  =,
                                   UCN = nonzero * P  + G,
                                   CU  =,
                                   CUN = nonzero * G  + P,
                                   UU  =,
                                   UUN = nonzero + G  * P)) +
                     ifelse(isFALSE(equalPro), G - !noise, noise)
  ll2             <- ll  * 2
  bic             <- ll2 - kpar * log(attr(seqs, "N"))
    return(list(bic = bic, icl = bic + 2L * sum(log(matrixStats::rowMaxs(z)), na.rm=TRUE), aic = ll2 - kpar * 2, df = kpar))
}

dbar              <- function(seqs, theta) {
    mean(stringdist::stringdistmatrix(seqs, theta, method="hamming"))
}

dseq              <- function(seqs, theta) {
    stringdist::stringdistmatrix(seqs, theta, method="hamming")
}

E_step            <- function(seqs, params, l.meth = c("CC", "UC", "CU", "UU", "CCN", "UCN", "CUN", "UUN"), ctrl = NULL, numseq = NULL) {
  G               <- attr(seqs, "G")
  N               <- attr(seqs, "N")
  P               <- attr(seqs, "P")
  V1              <- attr(seqs, "V1")
  G0              <- ifelse(noise <- attr(seqs, "Noise"), attr(seqs, "G0"), G)
  Gseq            <- seq_len(G0)
  N1              <- N > 1
  theta           <- params$theta
  lambda          <- switch(EXPR=l.meth, CC=, CCN=as.vector(params$lambda), params$lambda)
  dG.X            <- is.null(params$dG)
  if(is.null(numseq)       &&
     isFALSE(ctrl$numseq)  &&
    (G  == 1L     || dG.X) &&
     is.element(l.meth,
     c("CU",   "UU",
       "CUN", "UUN")))      {
    numseq        <- sapply(seqs, char_to_num)
  }

  if(G  == 1L)     {
    return(switch(EXPR=l.meth,
                  CC  =,
                  UC  = -ifelse(lambda == 0, N * P * log1p(V1), sum(lambda * N * dbar(seqs, theta), N * P * log1p(V1 * exp(-lambda)), na.rm=TRUE)),
                  CCN = -N * P * log1p(V1),
                  CU  =,
                  UU  = -sum(sweep(numseq != char_to_num(theta), 1L, lambda, FUN="*", check.margin=FALSE), na.rm=TRUE) - N * sum(log1p(V1 * exp(-lambda)))))
  } else {
    if(is.null(ctrl))            stop("'ctrl' must be supplied if G > 1", call.=FALSE)
    dG  <- if(dG.X)  switch(EXPR=l.meth,
                            CC=, UC=, CCN=, UCN=vapply(Gseq, function(g) dseq(seqs, theta[g]), numeric(N)),
                            CU=, UU=, CUN=, UUN=lapply(Gseq, function(g) unname(apply(numseq, 2L, "!=", char_to_num(theta[g]))))) else params$dG
    log.tau       <- log(params$tau)
    if(noise)      {
      tau0        <- log.tau[G]
      log.tau     <- log.tau[-G]
    }
    if(N1)         {
      numer       <- switch(EXPR=l.meth,
                            CC  =,
                            UC  = sweep(-sweep(dG, 2L, lambda, FUN="*", check.margin=FALSE),
                                        2L, log.tau - P * log1p(V1 * exp(-lambda)),  FUN="+", check.margin=FALSE),
                            CU  = sweep(-vapply(lapply(dG, sweep, 1L, lambda, "*", check.margin=FALSE), FUN=matrixStats::colSums2, na.rm=TRUE, numeric(N)),
                                        2L, log.tau - sum(log1p(V1 * exp(-lambda))), FUN="+", check.margin=FALSE),
                            UU  = sweep(-vapply(Gseq, function(g) matrixStats::colSums2(sweep(dG[[g]], 1L, lambda[g,], FUN="*", check.margin=FALSE), na.rm=TRUE), numeric(N)),
                                        2L, log.tau - matrixStats::rowSums2(log1p(V1 * exp(-lambda))), FUN="+", check.margin=FALSE),
                            CCN = cbind(sweep(-sweep(dG, 2L, lambda[1L],  FUN="*", check.margin=FALSE),
                                              2L, log.tau - P * log1p(V1 * exp(-lambda[1L])),   FUN="+", check.margin=FALSE), tau0 - P * log1p(V1)),
                            UCN = cbind(sweep(-sweep(dG, 2L, lambda[-G,], FUN="*", check.margin=FALSE),
                                              2L, log.tau - P * log1p(V1 * exp(-lambda[-G,])),  FUN="+", check.margin=FALSE), tau0 - P * log1p(V1)),
                            CUN = cbind(sweep(-vapply(lapply(dG, sweep, 1L, lambda[1L,], "*", check.margin=FALSE), FUN=matrixStats::colSums2, na.rm=TRUE, numeric(N)),
                                              2L, log.tau - sum(log1p(V1 * exp(-lambda[1L,]))), FUN="+", check.margin=FALSE), tau0 - P * log1p(V1)),
                            UUN = cbind(sweep(-vapply(Gseq, function(g) matrixStats::colSums2(sweep(dG[[g]], 1L, lambda[g,],  FUN="*", check.margin=FALSE), na.rm=TRUE), numeric(N)),
                                              2L, log.tau - matrixStats::rowSums2(log1p(V1 * exp(-lambda[-G,, drop=FALSE]))), FUN="+", check.margin=FALSE), tau0 - P * log1p(V1)))
    } else         {
      numer       <- switch(EXPR=l.meth,
                            CC  = -dG * lambda + log.tau - P * log1p(V1 * exp(-lambda)),
                            UC  = sweep(-dG * lambda, 2L, log.tau - P * log1p(V1 * exp(-lambda)), FUN="+", check.margin=FALSE),
                            CU  = -vapply(lapply(dG, sweep, 1L, lambda, "*", check.margin=FALSE), FUN=matrixStats::colSums2, na.rm=TRUE, numeric(N)) + log.tau - sum(log1p(V1 * exp(-lambda))),
                            UU  = -vapply(Gseq, function(g) matrixStats::colSums2(sweep(dG[[g]],   1L, lambda[g,], FUN="*", check.margin=FALSE), na.rm=TRUE), numeric(N)) + log.tau - matrixStats::rowSums2(log1p(V1 * exp(-lambda))),
                            CCN = c(-dG * lambda[1L] + log.tau - P * log1p(V1 * exp(-lambda[1L])), tau0 - P * log1p(V1)),
                            UCN = c(sweep(-dG * lambda[-G,, drop=FALSE], 2L, log.tau - P * log1p(V1 * exp(-lambda[-G])), FUN="+", check.margin=FALSE), tau0 - P * log1p(V1)),
                            CUN = c(-vapply(lapply(dG, sweep, 1L, lambda[1L,], "*", check.margin=FALSE), FUN=matrixStats::colSums2, na.rm=TRUE, numeric(N)) + log.tau - sum(log1p(V1 * exp(-lambda[1L,]))), tau0 - P * log1p(V1)),
                            UUN = c(-vapply(Gseq, function(g) matrixStats::colSums2(sweep(dG[[g]], 1L, lambda[g,], FUN="*", check.margin=FALSE), na.rm=TRUE), numeric(N)) + log.tau - matrixStats::rowSums2(log1p(V1 * exp(-lambda[-G,, drop=FALSE]))), tau0 - P * log1p(V1)))
      numer       <- matrix(numer, nrow=1L)
    }

    if(any(iLAM   <- is.infinite(lambda))) {
      switch(EXPR  = l.meth,
             CC    =,
             CCN   = numer[is.na(numer)]  <- matrix(log.tau - P * log1p(V1), nrow=N, ncol=G0, byrow=TRUE),
             UC    =,
             UCN   =         {
             i.x  <- which(dG == 0, arr.ind=TRUE)
             i.x  <- i.x[i.x[,2L] %in% which(switch(EXPR=l.meth, UC=iLAM, UCN=iLAM[-G])),, drop=FALSE]
             numer[i.x]     <- log.tau[i.x[,2L]] - P * log1p(V1)
      })
    }
    if(ctrl$do.cv && !noise &&
       any(i.x    <- apply(is.infinite(numer), 1L, all))) {
      numer[i.x,] <- matrix(log.tau - P * log1p(V1), nrow=sum(i.x), ncol=G, byrow=TRUE)
    }

    switch(EXPR=ctrl$algo,
           EM=     {
      if(N1)       {
        denom     <- matrixStats::rowLogSumExps(numer)
        loglike   <- sum(denom)
        if(ctrl$do.cv)       {
            return(loglike)
        } else     {
            return(list(loglike = loglike, z = exp(sweep(numer, 1L, denom, FUN="-", check.margin=FALSE))))
        }
      }   else     {
            return(matrixStats::logSumExp(numer))
      }
    },    CEM=     {
      if(N1)       {
        z         <- mclust::unmap(max.col(numer), groups=seq_len(G))
        loglike   <- sum(z * numer, na.rm=TRUE)
        if(ctrl$do.cv)       {
            return(loglike)
        } else     {
            return(list(loglike = loglike, z = z))
        }
      }   else     {
            return(max(numer))
      }
    })
  }
}

EM_algorithm      <- function(SEQ, numseq, g, modtype, z, ctrl, ll = NULL) {
  algo            <- ctrl$algo
  itmax           <- ctrl$itmax
  st.ait          <- ctrl$stopping == "aitken"
  tol             <- ctrl$tol
  if(!(runEM      <- g > 1))   {
    Mstep         <- M_step(SEQ, l.meth=modtype, numseq=numseq, ctrl=ctrl)
    ll            <- E_step(SEQ, l.meth=modtype, numseq=numseq, params=Mstep)
    j             <- 1L
    ERR           <- FALSE
  } else           {
    ll            <- if(is.null(ll)) c(-Inf, -sqrt(.Machine$double.xmax)) else ll
    j             <- 2L
    emptywarn     <- TRUE
    while(isTRUE(runEM))       {
      j           <- j + 1L
      Mstep       <- M_step(SEQ, l.meth=modtype, ctrl=ctrl, numseq=numseq, z=z)
      Estep       <- E_step(SEQ, l.meth=modtype, ctrl=ctrl, numseq=numseq, params=Mstep)
      z           <- Estep$z
      ERR         <- any(is.nan(z))
      if(isTRUE(ERR))            break
      if(any(colz <- matrixStats::colSums2(z)  == 0)   &&
         isTRUE(emptywarn))    { warning(paste0("\tThere were empty components: ", modtype, " (G=", g, ")\n"), call.=FALSE, immediate.=TRUE) # REORDER
        emptywarn <- FALSE
      }
      ll          <- c(ll, Estep$loglike)
      if(isTRUE(st.ait))       {
        ait       <- aitken(ll[seq(j - 2L, j, 1L)])
        dX        <- ifelse(is.numeric(ait$a)  && ait$a < 0, 0L, abs(ait$linf - ll[j - 1L]))
        dX[is.nan(dX)]     <- Inf
      } else       {
        dX        <- abs(ll[j] - ll[j - 1L])/(1 + abs(ll[j]))
      }
      runEM       <- dX >= tol && j   < itmax  && !ERR
    } # while (j)
  }
    return(list(ERR = ERR, j = j, Mstep = Mstep, ll = ll, z = z))
}

entropy           <- function(p) {
  p               <- p[p > 0]
    sum(-p * log(p))
}

fac_to_num        <- function(x) {
  Vseq            <- seq_len(nlevels(x[[1L]]))
    as.data.frame(lapply(x, function(y) { levels(y) <- Vseq; y} ))
}

get_partition              <- function(x, rank = 1L, criterion = c("bic", "icl", "aic", "cv", "loglik"), G = NULL, modtype = NULL, MAP = FALSE, noise = TRUE, ...) {
    UseMethod("get_partition")
}

get_partition.MEDseq       <- function(x, rank = 1L, criterion = c("bic", "icl", "aic", "cv", "loglik"), G = NULL, modtype = NULL, MAP = FALSE, noise = TRUE, ...) {
  ZS              <- x$ZS
  if(!(missing(G) -> m.G)  &&
    (length(G)    != 1     ||
     !is.numeric(G)        ||
     (G < 1 || floor(G) != G)))  stop("'G' must be a single integer >= 1", call.=FALSE)
  if(!(missing(modtype) ->
       m.M) &&
    (length(modtype) > 1   ||
     !is.character(modtype)))    stop("'modtype' must be a single character string", call.=FALSE)
  if(length(noise) > 1     ||
     !is.logical(noise))         stop("'noise' must be a single logical indicator", call.=FALSE)
  if(!missing(criterion)   ||
     !missing(rank)        ||
     any(m.G, m.M))         {
    criterion     <- match.arg(criterion)
    if(criterion  == "cv"  &&
       !attr(x, "CV"))           stop("Can't select based on the CV criterion as cross-validated likelihood wasn't performed", call.=FALSE)
    tmp           <- switch(EXPR=criterion, bic=x$BIC, icl=x$ICL, aic=x$AIC, cv=x$CV, loglik=x$LOGLIK)
    if(!noise)     {
      tmp         <- tmp[,colnames(tmp) %in% c("CC", "UC", "CU", "UU"), drop=FALSE]
    }
    if(!m.G)       {
      Gallow      <- as.numeric(rownames(tmp))
      if(!(G %in% Gallow))       stop("Invalid 'G'", call.=FALSE)
      tmp         <- tmp[as.character(G),, drop=FALSE]
    }
    if(!m.M)       {
      Mallow      <- colnames(tmp)
      if(!(modtype %in% Mallow)) stop("Invalid 'modtype'", call.=FALSE)
      tmp         <- tmp[,modtype, drop=FALSE]
    }
    class(tmp)    <- "MEDcriterion"
    if(length(rank) > 1    ||
       !is.numeric(rank)   ||
       rank != floor(rank) ||
       rank <= 0  ||
       rank  > sum(!is.na(tmp))) stop("Invalid 'rank'", call.=FALSE)
    best          <- print(tmp, pick=rank, show=FALSE)
    best          <- strsplit(names(best[rank]), ",")[[1L]]
    modtype       <- best[1L]
    G             <- as.numeric(best[2L])
  }
  if(length(MAP)     > 1   ||
     !is.logical(MAP))           stop("'MAP' must be a single logical indicator", call.=FALSE)
  if(!(G  %in%
       as.numeric(names(ZS))))   stop("Invalid 'G' value", call.=FALSE)
  Z               <- x$ZS[[as.character(G)]]
  z.ind           <- if(noise) names(Z) else names(Z)[names(Z) %in% c("CC", "UC", "CU", "UU")]
  if(!(modtype  %in% z.ind))     stop("Invalid 'modtype'", call.=FALSE)
  z               <- Z[[modtype]]
  if(anyNA(z))                   message("Selected model didn't converge: no partition available\n")
  noise           <- is.element(modtype, c("CCN", "UCN", "CUN", "UUN"))
  if(isTRUE(MAP))  {
    MAP           <- max.col(z)
    z             <- if(noise) replace(MAP, MAP == G, 0L) else MAP
  }
  attr(z, "G")             <- G
  attr(z, "ModelType")     <- modtype
  attr(z, "Noise")         <- noise
    return(z)
}

hamming_control   <- function(algo = c("EM", "CEM", "cemEM"), init.z = c("kmedoids", "hc", "random"), dist.mat = NULL, nstarts = 1L,
                              criterion = c("bic", "icl", "aic", "cv"), do.cv = TRUE, nfolds = 10L, stopping = c("aitken", "relative"),
                              tau0 = NULL, opti = c("mode", "first", "GA", "medoid"), ordering = c("none", "decreasing", "increasing"),
                              equalPro = FALSE, tol = 1E-05, itmax = .Machine$integer.max, nonzero = TRUE, verbose = TRUE, ...) {
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
  if(init.z == "random")    {
   if(length(nstarts) != 1 ||
      !is.numeric(nstarts) ||
      (nstarts         < 1 ||
       floor(nstarts) !=
       nstarts))                 stop(paste0("'nstarts' must be a single integer >= 1 if when 'init.z'=", init.z), call.=FALSE)
  }
  miss.args                <- list(do.cv=missing(do.cv), tau0=missing(tau0))
  if(!missing(criterion)   &&
    (length(criterion) > 1 ||
     !is.character(criterion)))  stop("'criterion' must be a character vector of length 1", call.=FALSE)
  if(length(do.cv)     > 1 ||
     !is.logical(do.cv))         stop("'do.cv' must be a single logical indicator",         call.=FALSE)
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
  if(length(equalPro)  > 1 ||
     !is.logical(equalPro))      stop("'equalPro' must be a single logical indicator",      call.=FALSE)
  if(length(tol)  != 1     ||
     !is.numeric(tol)      ||
     (tol          < 0     ||
      tol         >= 1))         stop("'tol' must be a scalar in the interval [0, 1)",      call.=FALSE)
  if(length(itmax)    != 1 ||
     !is.numeric(itmax)    ||
     itmax        <= 0)          stop("'itmax' must be a strictly positive scalar",         call.=FALSE)
  if(length(nonzero)   > 1 ||
     !is.logical(nonzero))       stop("'nonzero' must be a single logical indicator",       call.=FALSE)
  if(length(verbose)   > 1 ||
     !is.logical(verbose))       stop("'verbose' must be a single logical indicator",       call.=FALSE)
  control                  <- list(algo = match.arg(algo), init.z = init.z, dist.mat = dist.mat, nstarts = nstarts, criterion = match.arg(criterion), nfolds = nfolds, do.cv = do.cv, stopping = match.arg(stopping), tau0 = tau0,
                                   opti = match.arg(opti), ordering = match.arg(ordering), equalPro = equalPro, tol = tol, itmax = ifelse(is.infinite(itmax), .Machine$integer.max, itmax), nonzero = nonzero, verbose = verbose)
  attr(control, "missing") <- miss.args
    return(control)
}

hamming_fit       <- function(seqs, G = 1L:9L, l.meth = c("CC", "UC", "CU", "UU", "CCN", "UCN", "CUN", "UUN"), ctrl = hamming_control(...), ...) {
  call            <- match.call()
  if(!inherits(seqs, "stslist")) stop("'seqs' must be of class 'stslist'",        call.=FALSE)
  if(any(seqs     ==
         attr(seqs, "nr")))      stop("Missing values in 'seqs' are not allowed", call.=FALSE)
  SEQ             <- apply(fac_to_num(seqs), 1L, num_to_char)
  levs            <- attr(seqs, "alphabet")
  attr(SEQ, "N")  <- N   <- nrow(seqs)
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
  nonzero         <- ctrl$nonzero
  verbose         <- ctrl$verbose
  ctrl$ordering   <- ifelse(ctrl$opti == "first", ctrl$ordering, "none")
  miss.args       <- attr(ctrl, "missing")
  if(any(G        != floor(G))    &&
     any(G         < 1))         stop("'G' must be strictly positive", call.=FALSE)
  if(any(G        >= N))       {
    G             <- G[G <= N]
    if(length(G)   > 1)        { warning("Removing G values >= the number of observations\n",  call.=FALSE, immediate.=TRUE)
    } else                       stop("G values must be less than the number of observations", call.=FALSE)
  }
  G               <- rG  <- sort(as.integer(unique(G)))
  len.G           <- length(G)
  mt1             <- unique(vapply(l.meth, function(lx) switch(EXPR=lx, CC=, UC="CC", CU=, UU="CU", "CCN"), character(1L)))
  mt2             <- unique(vapply(l.meth, function(lx) switch(EXPR=lx, UCN="CCN", UUN="CUN", lx),          character(1L)))
  mtg             <- unique(l.meth)
  if(!(tmiss      <- miss.args$tau0) &&
     ctrl$tau0    == 0)        {
    n.meth        <- l.meth %in% c("CCN", "UCN", "CUN", "UUN")
    if(all(n.meth))             stop("'tau0' is zero: models with noise component will not be fitted", call.=FALSE)
    if(any(n.meth))             warning("'tau0' is zero: models with noise component will not be fitted\n", call.=FALSE, immediate.=TRUE)
    l.meths       <- c("CC", "UC", "CU", "UU")
  } else l.meths  <- c("CC", "UC", "CU", "UU", "CCN", "UCN", "CUN", "UUN")
  if(is.null(dist.mat))        {
    dist.mat      <- as.dist(suppressMessages(seqdist(seqs, method="OM", indel=1, sm=suppressMessages(seqsubm(seqs, method="TRATE")))))
  } else if((sqrt(length(dist.mat) * 2L + N) !=
             N))                stop("Invalid 'dist.mat' dimensions", call.=FALSE)
  if(any(G > 1L))  {
    G1            <- any(G == 1L)
    G2            <- any(G == 2L)
    GG            <- any(G  > 2L)
    all.mod       <- if(all(G1, G2, GG)) unique(c(mtg, mt2, mt1)) else if(all(G1, G2)) unique(c(mt2, mt1)) else if(all(G1, GG)) unique(c(mtg, mt1)) else if(all(G2, GG)) unique(c(mtg, mt2)) else if(G2) mt2 else mtg
    all.mod       <- l.meths[l.meths %in% all.mod]
    if(init.z     == "hc")     {
      hcZ         <- cluster::agnes(dist.mat, diss=TRUE, method="ward")
    }
  } else all.mod  <- l.meths[l.meths %in% mt1]
  if(any(len.G     > 1,
    (multi        <- length(all.mod) > 1))) {
    cvsel         <- ctrl$do.cv
  } else           {
    cvsel         <- ifelse(miss.args$do.cv, FALSE, ctrl$do.cv)
  }
  ctrl$do.cv      <- FALSE
  cvsel           <- any(cvsel, criterion == "cv")
  if(ctrl$numseq  <- any(c("CU", "UU", "CUN", "UUN") %in% all.mod, ctrl$opti == "mode", ctrl$ordering != "none")) {
    numseq        <- sapply(SEQ, char_to_num)
    attr(numseq, "P")    <- P
  } else numseq   <- NULL
  BICs            <-
  ICLs            <-
  AICs            <-
  DF.x            <-
  IT.x            <-
  Nzero.x         <-
  Ninfty.x        <- provideDimnames(matrix(NA, nrow=len.G, ncol=length(all.mod)), base=list(as.character(G), all.mod))
  ZS              <- replicate(len.G, list())
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
    cv.folds      <- cut(seq_len(N), breaks=nfolds, labels=FALSE)
    fcount        <- tabulate(cv.folds, nfolds)
    Nfcount       <- N - fcount
    foldseq       <- seq_len(nfolds)
  }
  crit.tx         <-
  crit.gx         <- -sqrt(.Machine$double.xmax)
  noise           <- all.mod %in% c("CCN", "UCN", "CUN", "UUN")
  nonoise         <- all(noise)
  noise           <- any(noise)

  G.last          <- G[len.G]
  for(g  in G)     {
    if(isTRUE(verbose))     {    cat(paste0("\n", g, " cluster model", ifelse(multi, "s", ""), " -\n"))
      last.G      <- g == G.last
    }
    if(ctrl$numseq)         {
      attr(numseq,   "G")  <-
      attr(SEQ,      "G")  <- g
    } else attr(SEQ, "G")  <- g
    h             <- which(G  == g)
    g0 <- attr(SEQ, "G0")  <- g - noise
    if(isTRUE(noise))   {
      tau0        <- ifelse(tmiss, 1/g, ctrl$tau0)
    }

    if(g  > 1)     {
      algog       <- algo
      if(init.z   == "random" &&
         nstarts   > 1)     {
        if(!nonoise)        {
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
        if(!nonoise)        {
          zg      <- mclust::unmap(switch(EXPR=init.z, kmedoids=cluster::pam(dist.mat, k=g,  cluster.only=TRUE), hc=cutree(hcZ, k=g)),  groups=seq_len(g))
        }
        if(isTRUE(noise))   {
          if(g0 > 1)        {
            zg0   <- mclust::unmap(switch(EXPR=init.z, kmedoids=cluster::pam(dist.mat, k=g0, cluster.only=TRUE), hc=cutree(hcZ, k=g0)), groups=seq_len(g0))
            zg0   <- cbind(zg0 * (1 - tau0), tau0)
          } else   {
            zg0   <- matrix(tau0, nrow=N, ncol=2L)
          }
        }
      }
      modtypes    <- if(g == 2L) mt2 else mtg
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

      zm          <- if(attr(SEQ, "Noise") <- gN0 <- is.element(modtype, c("CCN", "UCN", "CUN", "UUN"))) zg0 else zg
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
            EMX[[i]]       <- EM_algorithm(SEQ=SEQ, numseq=numseq, g=g, modtype=modtype, z=zm[[i]],    ctrl=ctrl)
            if(!EMX[[i]]$ERR)  {
              ctrl$algo    <- "EM"
              tmpEMX       <- EMX[[i]]
              EMX[[i]]     <- EM_algorithm(SEQ=SEQ, numseq=numseq, g=g, modtype=modtype, z=EMX[[i]]$z, ctrl=ctrl, ll=tmpEMX$ll[c(tmpEMX$j - 1L, tmpEMX$j)])
              if(EMX[[i]]$ERR) {
                EMX[[i]]   <- tmpEMX
                ctrl$algo  <- "CEM"
              }
            }
         }, EMX[[i]]       <- EM_algorithm(SEQ=SEQ, numseq=numseq, g=g, modtype=modtype, z=zm[[i]],    ctrl=ctrl))
        }
        EMX       <- EMX[[which.max(vapply(lapply(EMX, "[[", "ll"), max, numeric(1L)))]]
      } else       {
        switch(EXPR=algog,
              cemEM=        {
          ctrl$algo        <- "CEM"
          EMX              <- EM_algorithm(SEQ=SEQ, numseq=numseq, g=g, modtype=modtype, z=zm,         ctrl=ctrl)
          if(!EMX$ERR)      {
            ctrl$algo      <- "EM"
            tmpEMX         <- EMX
            EMX            <- EM_algorithm(SEQ=SEQ, numseq=numseq, g=g, modtype=modtype, z=EMX$z,      ctrl=ctrl, ll=tmpEMX$ll[c(tmpEMX$j - 1L, tmpEMX$j)])
            if(EMX$ERR)     {
              EMX          <- tmpEMX
              ctrl$algo    <- "CEM"
            }
          }
        }, EMX             <- EM_algorithm(SEQ=SEQ, numseq=numseq, g=g, modtype=modtype, z=zm,         ctrl=ctrl))
      }
      ERR         <- EMX$ERR
      j           <- EMX$j
      Mstep       <- EMX$Mstep
      ll          <- EMX$ll
      z           <- EMX$z
      cvsel.X     <- cvsel && !ERR
      j2          <- max(1L, j - switch(EXPR=algo, cemEM=1L, 2L))
      if(isTRUE(verbose))        cat(paste0("\t\t# Iterations: ", ifelse(ERR, "stopped at ", ""), j2, ifelse(last.G && last.T, "\n\n", "\n")))

      if(all((Mstep$lambda -> lambda) == 0) && cvsel.X) {
        CVll      <- -N * P * log(V)
      } else if(cvsel.X)    {
        lCV       <- vector("numeric", nfolds)
        zCV       <- z[cv.ind,, drop=FALSE]
        for(i in foldseq)   {
          testX   <- which(cv.folds == i)
          CVS     <- cv.SEQ[testX]
          SCV     <- cv.SEQ[-testX]
          CVz     <- zCV[-testX,, drop=FALSE]
          attributes(CVS)  <-
          attributes(SCV)  <- attributes(SEQ)
          attr(CVS, "N")   <- fcount[i]
          attr(SCV, "N")   <- Nfcount[i]
          if(any(modtype %in% c("CU", "UU", "CUN", "UUN"), ctrl$opti == "mode", ctrl$ordering != "none")) {
            nCV            <- cv.numseq[,-testX, drop=FALSE]
            CVn            <- cv.numseq[,testX,  drop=FALSE]
            attr(nCV, "G") <- attr(numseq, "G")
            attr(nCV, "P") <- attr(numseq, "P")
          } else   {
            nCV   <-
            CVn   <- NULL
          }
          EMX     <- EM_algorithm(SEQ=SCV, numseq=nCV, g=g, modtype=modtype, z=CVz, ctrl=ctrl)
          MCV     <- EMX$Mstep
          MCV$dG  <- NULL
          ctrl$do.cv       <- TRUE
          lCV[i] <- E_step(seqs=CVS, params=MCV, l.meth=modtype, ctrl=ctrl, numseq=CVn)
          ctrl$do.cv       <- FALSE
        }
        CVll      <- sum(lCV)
      }

      nzero       <- sum(lambda == 0)
      ninfty      <- sum(is.infinite(lambda))
      choice      <- choice_crit(ll=ll[j], seqs=SEQ, z=z, l.meth=modtype, nonzero=ifelse(nonzero, sum(lambda != 0), NA), equalPro=equalPro)
      bicx        <- choice$bic
      iclx        <- choice$icl
      aicx        <- choice$aic
      dfx         <- choice$df
      crit.t      <- switch(EXPR=criterion, cv=CVll, bic=bicx, icl=iclx, aic=aicx)
      crit.t      <- ifelse(is.na(crit.t) || ERR, -Inf, crit.t)
      if(crit.t    > crit.tx)  {
        crit.tx   <- crit.t
        theta.x   <- Mstep$theta
        lambda.x  <- lambda
        tau.x     <- Mstep$tau
        z.x       <- z
        ll.x      <- ll
      }
      BICs[h,modtype]         <- ifelse(ERR, -Inf, bicx)
      ICLs[h,modtype]         <- ifelse(ERR, -Inf, iclx)
      DF.x[h,modtype]         <- ifelse(ERR, -Inf, dfx)
      AICs[h,modtype]         <- ifelse(ERR, -Inf, aicx)
      IT.x[h,modtype]         <- ifelse(ERR,  Inf, j2)
      Nzero.x[h,modtype]      <- ifelse(ERR,  NA,  nzero)
      Ninfty.x[h,modtype]     <- ifelse(ERR,  NA,  ninfty)
      ZS[[h]][[m]]            <- if(ERR)      NA   else z
      if(cvsel)    {
        CV.x[h,modtype]       <- ifelse(ERR, -Inf, CVll)
      }
    } # for (modtype)

    if(crit.tx     > crit.gx)  {
      crit.gx     <- crit.tx
      x.theta     <- theta.x
      x.lambda    <- lambda.x
      x.tau       <- tau.x
      x.z         <- z.x
      x.ll        <- ll.x
    }
    ZS[[h]]       <- stats::setNames(ZS[[h]], modtypes)
  } # for (g)

  if(any(x.ll     !=
         cummax(x.ll)))          warning("Log-likelihoods are not strictly increasing\n", call.=FALSE)
  if(any(IT.x[!is.na(IT.x)]
         == ctrl$itmax))         warning(paste0("One or more models failed to converge in the maximum number of allowed iterations (", itmax, ")\n"), call.=FALSE)
  CRITs           <- switch(EXPR=criterion, cv=CV.x, bic=BICs, icl=ICLs, aic=AICs)
  best.ind        <- which(CRITs == crit.gx, arr.ind=TRUE)
  if(nrow(best.ind) > 1)       { warning(paste0("Ties for the optimal model exist according to the '", criterion, "' criterion: choosing the most parsimonious model\n"), call.=FALSE, immediate.=TRUE)
    best.ind      <- which(DF.x  == min(DF.x[best.ind]), arr.ind=TRUE)
    best.ind      <- best.ind[which.min(best.ind[,1L]),]
  }
  best.G          <- best.ind[1L]
  best.mod        <- colnames(CRITs)[best.ind[2L]]
  G               <- G[best.G]
  x.bic           <- BICs[best.ind]
  x.icl           <- ICLs[best.ind]
  x.aic           <- AICs[best.ind]
  LL.x            <- AICs/2 + DF.x
  class(BICs)     <-
  class(ICLs)     <-
  class(AICs)     <-
  class(LL.x)     <- "MEDcriterion"
  attr(BICs, "Criterion")        <- "BIC"
  attr(ICLs, "Criterion")        <- "ICL"
  attr(AICs, "Criterion")        <- "AIC"
  attr(LL.x, "Criterion")        <- "loglik"
  if(cvsel)        {
    CV.x          <- CV.x * 2
    x.cv          <- CV.x[best.ind]
    class(CV.x)   <- "MEDcriterion"
    attr(CV.x, "Criterion")      <- "CV"
  }
  if(len.G > 1    && verbose)     {
    if(G          == min(rG))    message("Best model occurs at the min of the number of components considered\n")
    if(G          == max(rG))    message("Best model occurs at the max of the number of components considered\n")
  }

  if(isTRUE(verbose))            cat(paste0("\n\t\tBest Model", ifelse(length(CRITs) > 1, paste0(" (according to ", toupper(criterion), "): "), ": "), best.mod, ", with ",  paste0(G, " component", ifelse(G > 1, "s", ""))),
                                     ifelse(cvsel, paste0("\n\t\tCV = ", round(x.cv, 2L), " |"), "\n\t       "),
                                     "BIC =", round(x.bic, 2L), "| ICL =", round(x.icl, 2L), "| AIC =", round(x.aic, 2L), "\n\n")
  attr(x.lambda, "Nzero")        <- Nzero.x[best.ind]
  attr(x.lambda, "Ninfty")       <- Ninfty.x[best.ind]
  attr(DF.x,     "Nzero")        <- Nzero.x
  attr(DF.x,     "Ninfty")       <- Ninfty.x
  x.theta                        <- do.call(rbind, lapply(x.theta, char_to_num))
  storage.mode(x.theta)          <- "integer"
  attr(x.theta, "alphabet")      <- levs
  attr(x.theta, "labels")        <- attr(seqs, "labels")
  attr(x.theta, "lambda")        <- switch(EXPR=best.mod, CCN=, CUN=rbind(matrix(x.lambda[1L,], nrow=G - 1L, ncol=P, byrow=best.mod == "CUN"), 0L), matrix(x.lambda, nrow=G, ncol=P, byrow=best.mod == "CU"))
  class(x.theta)                 <- "MEDtheta"
  params          <- list(theta   = x.theta,
                          lambda  = x.lambda,
                          tau     = x.tau)
  MAP             <- max.col(x.z)
  noise           <- best.mod %in% c("CCN", "UCN", "CUN", "UUN")
  MAP             <- if(noise) replace(MAP, MAP == G, 0L) else MAP
  results         <- list(call    = call,
                          data    = seqs,
                          modName = best.mod,
                          G       = G,
                          params  = params,
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
                          aic     = x.aic,
                          LOGLIK  = LL.x,
                          loglik  = x.ll[if(G > 1) switch(EXPR=algo, cemEM=-1L, -c(1L:2L)) else 1L],
                          uncert  = if(G > 1) 1 - matrixStats::rowMaxs(x.z) else vector("integer", N),
                          DF      = DF.x,
                          df      = DF.x[best.ind],
                          ITERS   = IT.x,
                          iters   = IT.x[best.ind]))
  dist.mat        <- as.matrix(dist.mat)
  unip            <- unique(MAP)
  Nseq            <- seq_len(N)
  meds            <- NULL
  set.seed(100)
  for(g in sort(unip[unip > 0]))  {
    srows         <- Nseq[MAP == g]
    meds          <- c(meds, srows[which.min(matrixStats::rowSums2(dist.mat[srows,srows]))])
  }
  perm            <- seriation::get_order(seriation::seriate(as.dist(dist.mat[meds,meds]), method="TSP"))
  attr(results, "Criterion")     <- criterion
  attr(results, "CV")            <- cvsel
  attr(results, "DistMat")       <- dist.mat
  attr(results, "EqualPro")      <- equalPro
  attr(results, "N")             <- N
  attr(results, "Noise")         <- noise
  attr(results, "P")             <- P
  attr(results, "Seriate")       <- if(noise) c(perm, G) else perm
  attr(results, "V")             <- V
  class(results)  <- "MEDseq"
    return(results)
}

lambda_mle        <- function(seqs, params, l.meth = c("CC", "UC", "CU", "UU", "CCN", "UCN", "CUN", "UUN"), numseq = NULL) {
  theta           <- params$theta
  P               <- attr(seqs, "P")
  V1V             <- attr(seqs, "V1V")
  lV1             <- attr(seqs, "logV1")
  l.meth          <- match.arg(l.meth)
  noise           <- attr(seqs, "Noise")
  if(is.null(numseq) && is.element(l.meth, c("CU", "UU", "CUN", "UUN"))) {
    numseq        <- sapply(seqs, char_to_num)
  }
  if((G           <- attr(seqs, "G")) == 1L) {
    if(l.meth == "CCN")                   {
        return(list(lambda = matrix(0L, nrow=1L, ncol=1L)))
    }
    numer         <- switch(EXPR=l.meth, CC=P,                 CU=1)
    denom         <- switch(EXPR=l.meth, CC=dbar(seqs, theta), CU=rowMeans(numseq != char_to_num(theta)))
  } else           {
    N             <- attr(seqs, "N")
    G0            <- ifelse(noise, attr(seqs, "G0"), G)
    if(is.null(z  <- params$z))  stop("'z' must be supplied if 'G'>1",   call.=FALSE)
    if(is.element(l.meth, c("UC", "UU", "UCN", "UUN")) &&
       is.null(tau    <-
                 params$tau))    stop("'tau' must be supplied if 'G'>1", call.=FALSE)
    if(noise)      {
      z           <- z[,-G, drop=FALSE]
      switch(EXPR=l.meth, UCN=, UUN=      {
        tau       <- tau[-G]
      })
    }
    switch(EXPR=l.meth, CU=, UU=, CUN=, UUN= {
      dG          <- lapply(seq_len(G0), function(g) unname(apply(numseq,                 2L, "!=",  char_to_num(theta[g]))))
      dGp         <- vapply(seq_len(G0), function(g) matrixStats::rowSums2(sweep(dG[[g]], 2L, z[,g], FUN="*", check.margin=FALSE)), numeric(P))
    },             {
      pN          <- P * switch(EXPR=l.meth, CCN=sum(z), N)
      dG          <- vapply(seq_len(G0), function(g) dseq(seqs, theta[g]), numeric(N))
    })
    numer         <- switch(EXPR=l.meth, CC=, CCN=pN,       CUN=sum(z),
                                         UC=, UCN=pN * tau, N)
    denom         <- switch(EXPR=l.meth, CC=, CCN=sum(z * dG),
                                         UC=, UCN=matrixStats::colSums2(z * dG),
                                         CU=, CUN=matrixStats::rowSums2(dGp),
                                         UU=, UUN=sweep(dGp, 2L, tau, FUN="/", check.margin=FALSE))
  }
  lambda          <- ifelse(denom < numer * V1V, lV1 + log(numer - denom) - log(denom), 0L)
  lambda          <- matrix(lambda, nrow=switch(EXPR=l.meth, UC=, UU=, UCN=, UUN=G0, 1L), ncol=switch(EXPR=l.meth, CU=, UU=, CUN=, UUN=P, 1L), byrow=is.element(l.meth, c("UU", "UUN")))
  switch(EXPR=l.meth, UC=,  UU=,
                     UCN=, UUN=           {
    lambda[tau == 0,] <- Inf
  })
    return(list(lambda = if(noise) rbind(lambda, 0L) else lambda, dG = if(G > 1) dG))
}

M_step            <- function(seqs, l.meth = c("CC", "UC", "CU", "UU", "CCN", "UCN", "CUN", "UUN"), ctrl, z = NULL, numseq = NULL) {
  G               <- attr(seqs, "G")
  if(is.null(numseq)     && isFALSE(ctrl$numseq)) {
    numseq        <- sapply(seqs, char_to_num)
    attr(numseq,  "G")   <- G
    attr(numseq,  "P")   <- attr(seqs, "P")
  }
  if(G > 1L)       {
    if(is.null(z))               stop("'z' must be supplied when 'G'>1", call.=FALSE)
    tau2          <- colMeans(z)
    tau           <- if(isFALSE(ctrl$equalPro)) tau2 else if(attr(seqs, "Noise")) c(rep((1 - tau2[G])/attr(seqs, "G0"), attr(seqs, "G0")), tau2[G]) else rep(1/G, G)
  } else tau2     <- tau <- 1L
  theta           <- optimise_theta(seqs=seqs, ctrl=ctrl, z=z, numseq=numseq)
  MLE             <- lambda_mle(seqs=seqs, params=list(theta=theta, z=z, tau=tau2), l.meth=l.meth, numseq=numseq)
  param           <- list(theta=theta, lambda=MLE$lambda, dG=if(G > 1) MLE$dG, tau=tau)
  attr(param, "l.meth")  <- l.meth
    return(param)
}

modal             <- function(x) {
  ux              <- unique(x[x != 0])
  tab             <- tabulate(match(x, ux))
    ux[which(tab  == max(tab))]
}

num_to_char       <- function(x) {
    paste(x, sep="", collapse="")
}

optimise_theta    <- function(seqs, ctrl, z = NULL, numseq = NULL) {
  opti            <- ctrl$opti
  ordering        <- ctrl$ordering
  P               <- attr(seqs, "P")
  G               <- attr(seqs, "G")
  Gseq            <- seq_len(G)
  if(G     > 1L   && is.null(z))  stop("'z' must be supplied when 'G'>1", call.=FALSE)
  if(opti         == "mode" || ordering != "none")   {
    if(is.null(numseq)      && isFALSE(ctrl$numseq)) {
      numseq      <- sapply(seqs, char_to_num)
      attr(numseq, "G")     <- G
      attr(numseq, "P")     <- P
    }
    if(opti       == "mode") {
      if(G > 1L)   {
        theta     <- weighted_mode(numseq=numseq, z=z)
        if(is.list(theta)   && any(nonu <- apply(theta,   2L, function(x) any(nchar(x) > 1)))) {
          theta[,nonu]      <- lapply(theta[,nonu], "[[", 1L)
        } else       nonu   <- rep(FALSE, G)
        theta     <- apply(theta,  2L, num_to_char)
      } else       {
        theta     <- apply(numseq, 1L, modal)
        if(nonu   <- is.list(theta)) {
          theta   <- lapply(theta, "[[",  1L)
        }
        theta     <- num_to_char(theta)
      }
      attr(theta, "NonUnique")   <- nonu
        return(theta)
    }
  }

  V               <- attr(seqs, "V")
  theta.opt       <- theta_data(seqs=seqs, z=z)
  if(opti == "medoid")       {
      return(theta.opt$theta)
  }
  N               <- attr(seqs, "N")
  pseq            <- seq_len(P)
  vseq            <- seq_len(V)
  opt             <- theta.opt$dsum
  if(G > 1L)       {
    theta         <- lapply(theta.opt$theta, char_to_num)
  } else           {
    z             <- matrix(1L, nrow=N)
    theta         <- list(char_to_num(theta.opt$theta))
  }
  if(ordering     != "none") {
    stab          <- if(G   == 1) list(apply(numseq, 1L, tabulate, V)) else lapply(Gseq, function(g) apply(sweep(numseq, 2L, z[,g], FUN="*", check.margin=FALSE), 1L, tabulate, V))
    sorder        <- lapply(Gseq, function(g) order(apply(stab[[g]], 2L, entropy), decreasing=ordering == "decreasing"))
    theta         <- lapply(Gseq, function(g) theta[[g]][sorder[[g]]])
    seqs          <- lapply(Gseq, function(g) unname(apply(numseq[sorder[[g]],], 2L, num_to_char)))
  } else seqs     <- replicate(G, list(seqs))

  for(g in Gseq)   {
    p.opt         <- Inf
    switch(EXPR    = opti,
           first   =     {
             opts       <- rep(NA, V)
             while(p.opt > opt[g])       {
               p.opt    <- opt[g]
               for(p    in pseq)         {
                 vdiff  <- setdiff(vseq, theta[[g]][p])
                 for(v in vdiff)         {
                   opts[v]              <- sum(z[,g] * dseq(seqs[[g]], paste(replace(theta[[g]], p, v), collapse="")))
                 }
                 opts[theta[[g]][p]]    <- opt[g]
                 opt[g] <- min(opts)
                 theta[[g]][p]          <- which.min(opts)
               }
             }
           }, GA   =     {
             opts       <- matrix(NA, nrow=P, ncol=V)
             while(p.opt > opt[g])       {
               p.opt    <- opt[g]
               for(p   in pseq)          {
                 vdiff  <- setdiff(vseq, theta[[g]][p])
                 for(v in vdiff)         {
                   opts[p, v]           <- sum(z[,g] * dseq(seqs[[g]], paste(replace(theta[[g]], p, v), collapse="")))
                 }
                 opts[p, theta[[g]][p]] <- opt[g]
               }
               opt[g]   <- min(opts)
               if(p.opt  > opt[g])       {
                 oi     <- which(opts == opt[g], arr.ind=TRUE)[1L,]
                 theta[[g]][oi[1L]]     <- oi[2L]
               }
             }
           })
  }

  theta          <- switch(EXPR=ordering,
                           none=if(G > 1)  {
                             do.call(base::c, lapply(theta, num_to_char))
                           } else num_to_char(theta[[1L]]),
                           if(G > 1)       {
                             do.call(base::c, lapply(Gseq, function(g) num_to_char(theta[[g]][match(pseq, sorder[[g]])])))
                           } else num_to_char(theta[[1L]][match(pseq, sorder[[1L]])]))
}

plot.MEDseq       <- function(x, type = c("clusters", "mean", "lambda", "cv", "bic", "icl", "aic", "LOGLIK", "quality", "uncert.bar", "uncert.profile",
                              "loglik", "d", "f", "Ht", "i", "I"), seriate = TRUE, preczero = TRUE, log.scale = NULL, qual = "ASW", ...) {
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
  dat             <- x$data
  modName         <- x$modName
  Gseq    <- perm <- seq_len(G)
  Nseq            <- seq_len(N)
  Pseq            <- seq_len(P)
  symbols         <- c(17, 2, 16, 10, 13, 18, 15, 7)
  use.col         <- c("gray", "dodgerblue1", "red3", "slateblue", "green3", "skyblue1", "gold", "hotpink")
  alpha.x         <- attr(dat, "alphabet")
  cpal.x          <- attr(dat, "cpal")
  label.x <- lab  <- attr(dat, "labels")
  if(is.element(type,
     c("mean", "lambda")))     {
    theta         <- x$params$theta
    lambda        <- attr(theta, "lambda")
    class(theta)  <- NULL
  }
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
      dots        <- list(...)
      dots        <- dots[names(dots) != "MAP"]
      if(length(dots) > 0)     {
        MAP       <- do.call(get_partition, c(list(x=x, MAP=TRUE), dots))
        G         <- attr(MAP, "G")
        Gseq      <-
        perm      <- seq_len(G)
        noise     <- attr(MAP, "Noise")
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
    OrderedStates <- data.matrix(fac_to_num(dat))[glo.order,]
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
    plot(0, 0, type = "n", ylim=if(length(dat) > 1) range(dat, na.rm=TRUE) else c(0, 1), yaxt="n", ylab="", xaxt="n", xlab="", frame.plot=FALSE)
    axis(side=4, las=2, tick=FALSE, line=0.5)
    par(mar=c(0, 0, 0, 0))
    plot.new()
    legend("center", c(expression(paste(lambda, " = 0")), expression(paste(lambda %->% infinity))), fill=c("green3", "grey50"), ncol=2, text.width=0.1, cex=1.25)
    graphics::layout(1)
      invisible()
  }, cv=,
     bic=,
     icl=,
     aic=,
     LOGLIK=       {
    if(all(type   == "cv",
       !attr(x, "CV")))          stop("Cross-validated log-likelihood values cannot be plotted as cross-validation didn't take place during model fit\n", call.=FALSE)
    dat           <- switch(EXPR=type, cv=x$CV, bic=x$BIC, icl=x$ICL, aic=x$AIC, LOGLIK=x$LOGLIK)
    ms            <- which(c("CC", "UC", "CU", "UU", "CCN", "UCN", "CUN", "UUN") %in% colnames(dat))
    symbols       <- symbols[ms]
    use.col       <- use.col[ms]
    matplot(dat, type="b", xlab="Number of Components", ylab=switch(EXPR=type, cv=, LOGLIK="", toupper(type)), col=use.col, pch=symbols, ylim=range(as.vector(dat[!is.na(dat) & is.finite(dat)])), xaxt="n", lty=1)
    if(type == "cv")     mtext(expression("\u2113"["cv"]), side=2, line=3, las=1, cex=1.5)
    if(type == "LOGLIK") mtext("Log-Likelihood", side=2, line=3, las=3)
    axis(1, at=seq_len(nrow(dat)))
    legend("bottomright", ncol=2, cex=1, inset=0.01, legend=colnames(dat), pch=symbols, col=use.col)
      invisible()
  }, quality=      {
    Qual          <- list()
    ZS            <- x$ZS
    Gseq          <- seq_along(ZS)
    if(all(Gseq   == 1)) {    message("No solutions with G > 1 clusters\n")
        break
    }
    Gseq          <- Gseq[Gseq > 1]
    znames        <- names(ZS[[which.max(lengths(ZS))]])
    for(g   in Gseq)     {
      z           <- ZS[[g]]
      qg          <- provideDimnames(matrix(NA, nrow=length(znames), ncol=10), base=list(znames, c("PBC", "HG", "HGSD", "ASW", "ASWw", "CH", "R2", "CHsq", "R2sq", "HC")))
      for(m in names(z)) {
        qg[m,]    <- WeightedCluster::wcClusterQuality(dmat, max.col(z[[m]]))$stats
      }
      Qual[[g]]   <- qg
    }
    Qual          <- stats::setNames(Qual[-1L], as.character(Gseq))
    if(length(qual)      > 1 ||
       !is.character(qual)   ||
       !(qual %in%
         colnames(Qual[[1L]])))  stop("Invalid 'qual' measure", call.=FALSE)
    ms            <- which(c("CC", "UC", "CU", "UU", "CCN", "UCN", "CUN", "UUN") %in% znames)
    symbols       <- symbols[ms]
    use.col       <- use.col[ms]
    quality       <- sapply(Qual, "[" , TRUE, qual)
    matplot(t(quality), type="b", pch=symbols, col=use.col, lty=1, ylab=qual, xaxt="n",
            xlab="Number of Components", main=paste0("Cluster Quality Index: ", qual))
    axis(1, at=seq_along(Gseq), labels=Gseq)
    legend("bottomright", ncol=2, cex=1, inset=0.01, legend=znames, pch=symbols, col=use.col)
    attr(quality, "quality") <- qual
      invisible(list(all = Qual, qual = quality))
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
             cu   <- cm[1L:2L][(uncX >= oneG) + 1L]
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
    plot(x, type=ifelse(length(x) == 1, "p", "l"), xlab="Iterations", ylab="Log-Likelihood", xaxt="n")
    seqll         <- seq_along(x)
    llseq         <- pretty(seqll)
    llseq         <- if(any(llseq != floor(llseq))) seqll else llseq
    axis(1, at=llseq, labels=llseq)
      invisible(list(ll = x, converge = x[length(x)]))
  },               {
    MAP           <- factor(replace(MAP, MAP == 0, "Noise"), levels=perm)
    dots          <- dots[!(names(dots) %in% c("G", "modtype", "noise"))]
    dots          <- c(dots, list(seqdata=dat, with.legend=FALSE, group=MAP, type=type, with.missing=FALSE))
    dots          <- if(type == "i") dots else c(dots, list(border=NA))
    if(type       != "Ht")     {
      l.ncol      <- ceiling(V/ifelse(V > 6 | G %% 2 != 0, 3, 2))
      if(G  %% 2  == 0)        {
        par(oma=c(7.5, 0, 0, 0),       xpd=TRUE)
        try(suppressWarnings(do.call(seqplot, dots)), silent=TRUE)
        par(fig=c(0, 1, 0, 1), oma=c(0.5, 0, 0, 0), mar=c(0.5, 0, 0, 0), new=TRUE)
        plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
        legend("bottom",      fill=attr(dat, "cpal"), legend=lab, ncol=l.ncol)
      } else       {
        try(suppressWarnings(do.call(seqplot, dots)), silent=TRUE)
        par(mar=c(1, 1, 0.5, 1) + 0.1, xpd=TRUE)
        legend("bottomright", fill=attr(dat, "cpal"), legend=lab, ncol=l.ncol)
      }
    } else         {
      try(suppressWarnings(do.call(seqplot, dots)),   silent=TRUE)
    }
      invisible()
  })
}

print.MEDcriterion       <- function(x, pick = 3L, show = TRUE, ...) {
  if(length(pick)        != 1 ||
     !is.numeric(pick))          stop("'pick' must be a single number", call.=FALSE)
  if(floor(pick)  != pick     ||
     pick          < 1)          stop("'pick' be a strictly positive integer", call.=FALSE)
  if(length(show) != 1        ||
     !is.logical(show))          stop("'show' must be a single logical indicator", call.=FALSE)
  crit            <- attr(x, "Criterion")
  x2              <- replace(x, !is.finite(x), NA)
  pick            <- min(pick,        length(x2[!is.na(x2)]))
  x.crit          <- x2  >= sort(x2,  decreasing=TRUE)[pick]
  x.ind           <- which(x.crit,    arr.ind=TRUE)
  x.val           <- sort(x2[x.ind],  decreasing=TRUE)
  ind.x           <- order(x2[x.ind], decreasing=TRUE)
  x.ind           <- x.ind[ind.x,,    drop=FALSE]
  x.ind[,1L]      <- gsub(".*= ", "", rownames(x)[x.ind[,1L]])
  x.ind[,2L]      <- colnames(x2)[as.numeric(x.ind[,2L])]
  crits           <- stats::setNames(x.val, vapply(seq_len(pick), function(p, b=x.ind[p,]) paste0(b[2L], ",", b[1L]), character(1L)))
  dim1            <- attr(x, "dim")
  dim2            <- attr(x, "dimnames")
  attributes(x)   <- NULL
  attr(x, "dim")       <- dim1
  attr(x, "dimnames")  <- dim2
  if(isTRUE(show)) {
    cat(switch(EXPR= crit,
               CV="Cross-Validated Log-Likelihood (CV):\n",
               BIC="Bayesian Information Criterion (BIC):\n",
               ICL="Integrated Completed Likelihood (ICL):\n",
               AIC="Akaike Information Criterion (AIC):\n"),
            LOGLIK="Maximal Log-Likelihood:\n")
    print(unclass(x))
    cat(paste0("\nTop ", ifelse(pick > 1, paste0(pick, " models"), "model"), " based on the ", crit, " criterion:\n"))
    print(crits)
  }
    invisible(crits)
}

print.MEDseq      <- function(x, digits = 2L, ...) {
  cat("Call:\t");  print(x$call)
  if(length(digits)  > 1   || !is.numeric(digits) ||
     digits       <= 0)          stop("Invalid 'digits'", call.=FALSE)
  name            <- x$modName
  G               <- x$G
  equalP          <- G < 1 || attr(x, "EqualPro")
  crit            <- round(unname(c(x$bic, x$icl, x$aic)), digits)
  cat(paste0("\nBest Model", ifelse(length(x$BIC)  > 1, paste0(" (according to ", toupper(attr(x, "Criterion")), "): "), ": "), name, ", with ",
             G, " component",      ifelse(G > 1, "s\n", "\n"),
             ifelse(!equalP, "",   paste0("Equal Mixing Proportions\n")),
             ifelse(attr(x, "CV"), paste0("CV = ", round(unname(x$cv), digits), " | "), ""),
             "BIC = ",             crit[1L],
             " | ICL = ",          crit[2L],
             " | AIC = ",          crit[3L]))
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
  theta           <- as.data.frame(lapply(as.data.frame(x), function(theta) replace_levels(num_to_char(theta), alpha)))
    print(theta, ...)
}

replace_levels    <- function(seq, levels = NULL) {
  seq             <- as.numeric(strsplit(seq, "")[[1L]])
    if(is.null(levels)) factor(seq) else factor(seq, levels=seq_along(levels) - any(seq == 0), labels=as.character(levels))
}

summary.MEDseq    <- function(object, digits = 2L, ...) {
  x               <- object
  cat("Call:\t");  print(x$call)
  if(length(digits)  > 1   || !is.numeric(digits) ||
     digits       <= 0)          stop("Invalid 'digits'", call.=FALSE)
  name            <- x$modName
  G               <- x$G
  equalP          <- G < 1 || attr(x, "EqualPro")
  crit            <- round(unname(c(x$bic, x$icl, x$aic)), digits)
  cat(paste0("\nBest Model", ifelse(length(x$BIC)  > 1, paste0(" (according to ", toupper(attr(x, "Criterion")), "): "), ": "), name, ", with ",
             G, " component",      ifelse(G > 1, "s\n", "\n"),
             ifelse(!equalP, "",   paste0("Equal Mixing Proportions\n")),
             ifelse(attr(x, "CV"), paste0("CV = ", round(unname(x$cv), digits), " | "), ""),
             "BIC = ",             crit[1L],
             " | ICL = ",          crit[2L],
             " | AIC = ",          crit[3L]))
    invisible()
}

tabcluster        <- function(x, norm = FALSE) {
    UseMethod("tabcluster")
}

tabcluster.MEDseq <- function(x, norm = FALSE) {
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
  temp            <- sweep(temp, 1L, tabMAP, FUN="/")/P * ifelse(isTRUE(norm), 100L, 1L)
  temp            <- cbind(tabMAP, temp)
  rownames(temp)  <- gnames
  colnames(temp)  <- c("Size", alph)
    temp
}

theta_data        <- function(seqs, z = NULL) {
  if((G <- attr(seqs, "G"))   == 1L)  {
    sumdist       <- vapply(seq_len(attr(seqs, "N")), function(i) sum(dseq(seqs, seqs[i])), numeric(1L))
    sumdist       <- sumdist[sumdist  > 0]
      return(list(theta = seqs[which.min(sumdist)], dsum = min(sumdist)))
  } else           {
    if(is.null(z))               stop("'z' must be supplied if 'G'>1", call.=FALSE)
    theta         <- dsum     <- list()
    for(g in seq_len(G))       {
      sumdist     <- vapply(seq_len(attr(seqs, "N")), function(i) sum(dseq(seqs, seqs[i]) * z[,g]), numeric(1L))
      sumdist     <- sumdist[sumdist  > 0]
      theta[[g]]  <- seqs[which.min(sumdist)]
      dsum[[g]]   <- min(sumdist)
    }
      return(list(theta = do.call(base::c, theta), dsum = do.call(base::c, dsum)))
  }
}

weighted_mode     <- function(numseq, z) {
    sapply(seq_len(attr(numseq, "G")), function(g)
    sapply(seq_len(attr(numseq, "P")), function(p, x=tapply(z[,g], numseq[p,], sum)) names(which(x == max(x)))))
}

