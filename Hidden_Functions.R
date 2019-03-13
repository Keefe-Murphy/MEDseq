.aitken           <- function(loglik) {
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

.char_to_num      <- function(x) {
    as.numeric(strsplit(x, "")[[1L]])
}

.choice_crit      <- function(ll, seqs, z, gate.pen, l.meth = c("CC", "UC", "CU", "UU", "CCN", "UCN", "CUN", "UUN"), nonzero = NULL) {
  G               <- attr(seqs, ifelse(noise <- is.element(l.meth, c("CCN", "UCN", "CUN", "UUN")), "G0", "G"))
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
                                   UUN = nonzero + G  * P)) + gate.pen
  ll2             <- ll  * 2
  bic             <- ll2 - kpar * log(attr(seqs, "W"))
    return(list(bic = bic, icl = bic + 2L * sum(log(matrixStats::rowMaxs(z)), na.rm=TRUE), aic = ll2 - kpar * 2, df = kpar))
}

.crits_names      <- function(x) {
    unlist(lapply(seq_along(x), function(i) stats::setNames(x[[i]], paste0(names(x[i]), "|", names(x[[i]])))))
}

.dbar             <- function(seqs, theta) {
    mean(stringdist::stringdistmatrix(seqs, theta, method="hamming"))
}

.dseq             <- function(seqs, theta) {
    stringdist::stringdistmatrix(seqs, theta, method="hamming")
}

.drop_constants   <- function(dat, formula, sub = NULL) {
  if(!is.data.frame(dat))        stop("'dat' must be a data.frame", call.=FALSE)
  Nseq            <- seq_len(nrow(dat))
  sub             <- if(missing(sub)) Nseq else sub
  numsubs         <- all(is.numeric(sub))
  if(!any(numsubs, all(is.logical(sub)) &&
     length(sub)  == nrow(dat))) stop("'sub' must be a numeric vector, or logical vector with length equal to the number of rows in 'dat'", call.=FALSE)
  if(numsubs      &&
     any(match(sub, Nseq, 
     nomatch = 0) == 0))         stop("Numeric 'sub' must correspond to row indices of data", call.=FALSE)
  if(!inherits(formula, 
               "formula"))       stop("'formula' must actually be a formula!", call.=FALSE)
  intercept       <- attr(stats::terms(formula), "intercept")
  dat             <- dat[sub,colnames(dat) %in% all.vars(stats::update.formula(formula, 0 ~ .)), drop=FALSE]
  ind             <- names(which(!apply(dat, 2L, function(x) all(x == x[1L], na.rm=TRUE))))
  fterms          <- attr(stats::terms(formula), "term.labels")
  ind             <- c(ind[ind %in% fterms], fterms[grepl(":", fterms) | grepl("I\\(", fterms)])
  response        <- all.vars(stats::update.formula(formula, . ~ 0))
  form            <- if(length(ind) > 0) stats::reformulate(ind, response=response) else stats::as.formula(paste0(response, " ~ 1"))
  form            <- if(intercept  == 0) stats::update.formula(form, ~ . -1)        else form
  environment(form)        <- environment(formula)
    form
}

.E_step           <- function(seqs, params, l.meth = c("CC", "UC", "CU", "UU", "CCN", "UCN", "CUN", "UUN"), ctrl, numseq = NULL) {
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
     ctrl$nmeth)   {
    numseq        <- sapply(seqs, .char_to_num)
  }
  
  if(G  == 1L)     { 
    return(if(ctrl$do.wts)  {
             switch(EXPR=l.meth, 
                    CC  = -ifelse(lambda == 0, attr(seqs, "W") * P * log1p(V1), lambda * sum(attr(seqs, "Weights") * .dseq(seqs, theta), na.rm=TRUE) + attr(seqs, "W") * P * log1p(V1 * exp(-lambda))),
                    CCN = -attr(seqs, "W") * P * log1p(V1),
                    CU  = -sum(sweep(sweep(numseq != .char_to_num(theta), 1L, lambda, FUN="*", check.margin=FALSE), 2L, attr(seqs, "Weights"), FUN="*", check.margin=FALSE), na.rm=TRUE) - attr(seqs, "W") * sum(log1p(V1 * exp(-lambda))))
           } else  {
             switch(EXPR=l.meth,
                    CC  = -ifelse(lambda == 0, N * P * log1p(V1), sum(lambda * N * .dbar(seqs, theta), N * P * log1p(V1 * exp(-lambda)), na.rm=TRUE)),
                    CCN = -N * P * log1p(V1),
                    CU  = -sum(sweep(numseq != .char_to_num(theta), 1L, lambda, FUN="*", check.margin=FALSE), na.rm=TRUE) - N * sum(log1p(V1 * exp(-lambda))))
           })
  } else {
    dG  <- if(dG.X)  switch(EXPR=l.meth, 
                            CC=, UC=, CCN=, UCN=vapply(Gseq, function(g) .dseq(seqs, theta[g]), numeric(N)),
                            CU=, UU=, CUN=, UUN=lapply(Gseq, function(g) unname(apply(numseq, 2L, "!=", .char_to_num(theta[g]))))) else params$dG
    log.tau       <- .mat_byrow(log(params$tau), nrow=N, ncol=G)
    if(noise)      {
      tau0        <- log.tau[,G]
      log.tau     <- log.tau[,-G, drop=FALSE]
    }
    if(N1)         {
      numer       <- switch(EXPR=l.meth,
                            CC  =-sweep(dG, 2L, lambda, FUN="*", check.margin=FALSE) + log.tau - P * log1p(V1 * exp(-lambda)),
                            UC  = sweep(-sweep(dG, 2L, lambda, FUN="*", check.margin=FALSE), 
                                        2L, P * log1p(V1 * exp(-lambda)),  FUN="-", check.margin=FALSE) + log.tau,
                            CU  = sweep(-vapply(lapply(dG, sweep, 1L, lambda, "*",  check.margin=FALSE), FUN=matrixStats::colSums2, na.rm=TRUE, numeric(N)), 
                                        2L, sum(log1p(V1 * exp(-lambda))), FUN="-", check.margin=FALSE) + log.tau,
                            UU  = sweep(-vapply(Gseq, function(g) matrixStats::colSums2(sweep(dG[[g]], 1L, lambda[g,], FUN="*", check.margin=FALSE), na.rm=TRUE), numeric(N)), 
                                        2L, matrixStats::rowSums2(log1p(V1 * exp(-lambda))), FUN="-", check.margin=FALSE) + log.tau,
                            CCN = cbind(sweep(-sweep(dG, 2L, lambda[1L],  FUN="*",  check.margin=FALSE), 
                                              2L, P * log1p(V1 * exp(-lambda[1L])),   FUN="-", check.margin=FALSE) + log.tau, tau0 - P * log1p(V1)),
                            UCN = cbind(sweep(-sweep(dG, 2L, lambda[-G,], FUN="*",  check.margin=FALSE), 
                                              2L, P * log1p(V1 * exp(-lambda[-G,])),  FUN="-", check.margin=FALSE) + log.tau, tau0 - P * log1p(V1)),
                            CUN = cbind(sweep(-vapply(lapply(dG, sweep, 1L, lambda[1L,], "*",  check.margin=FALSE), FUN=matrixStats::colSums2, na.rm=TRUE, numeric(N)), 
                                              2L, sum(log1p(V1 * exp(-lambda[1L,]))), FUN="-", check.margin=FALSE) + log.tau, tau0 - P * log1p(V1)),
                            UUN = cbind(sweep(-vapply(Gseq, function(g) matrixStats::colSums2(sweep(dG[[g]], 1L, lambda[g,],  FUN="*", check.margin=FALSE), na.rm=TRUE), numeric(N)),
                                              2L, matrixStats::rowSums2(log1p(V1 * exp(-lambda[-G,, drop=FALSE]))), FUN="-", check.margin=FALSE) + log.tau, tau0 - P * log1p(V1)))
    } else         {
      numer       <- switch(EXPR=l.meth,
                            CC  = -dG * lambda - P * log1p(V1 * exp(-lambda)),
                            UC  = sweep(-dG * lambda, 2L, P * log1p(V1 * exp(-lambda)), FUN="-", check.margin=FALSE),
                            CU  = -vapply(lapply(dG, sweep, 1L, lambda, "*", check.margin=FALSE), FUN=matrixStats::colSums2, na.rm=TRUE, numeric(N)) - sum(log1p(V1 * exp(-lambda))),
                            UU  = -vapply(Gseq, function(g) matrixStats::colSums2(sweep(dG[[g]], 1L, lambda[g,], FUN="*", check.margin=FALSE), na.rm=TRUE), numeric(N)) - matrixStats::rowSums2(log1p(V1 * exp(-lambda))),
                            CCN = c(-dG * lambda[1L] - P * log1p(V1 * exp(-lambda[1L])), - P * log1p(V1)),
                            UCN = c(sweep(-dG * lambda[-G,, drop=FALSE], 2L, P * log1p(V1 * exp(-lambda[-G])), FUN="-", check.margin=FALSE), -P * log1p(V1)),
                            CUN = c(-vapply(lapply(dG, sweep, 1L, lambda[1L,], "*", check.margin=FALSE), FUN=matrixStats::colSums2, na.rm=TRUE, numeric(N)) - sum(log1p(V1 * exp(-lambda[1L,]))), - P * log1p(V1)),
                            UUN = c(-vapply(Gseq, function(g) matrixStats::colSums2(sweep(dG[[g]], 1L, lambda[g,], FUN="*", check.margin=FALSE), na.rm=TRUE), numeric(N)) - matrixStats::rowSums2(log1p(V1 * exp(-lambda[-G,, drop=FALSE]))), - P * log1p(V1)))
      numer       <- matrix(numer, nrow=1L) + switch(EXPR=l.meth, CC=, UC=, CU=, UU=log.tau, c(log.tau, tau0))
    }
    
    if(any(iLAM   <- is.infinite(lambda))) {
      switch(EXPR  = l.meth,
             CC    =,
             CCN   = numer[is.na(numer)]  <- log.tau - P * log1p(V1),
             UC    =,
             UCN   =         {
             i.x  <- which(dG == 0, arr.ind=TRUE)
             i.x  <- i.x[i.x[,2L] %in% which(switch(EXPR=l.meth, UC=iLAM, UCN=iLAM[-G])),, drop=FALSE]
             numer[i.x]     <- log.tau[,i.x[,2L]]    - P * log1p(V1)
      })
    }
    if(ctrl$do.cv && !noise &&
       any(i.x    <- apply(is.infinite(numer), 1L, all))) {
      numer[i.x,] <- log.tau[i.x,, drop=FALSE]       - P * log1p(V1)
    }
    
    switch(EXPR=ctrl$algo,
           EM=     {
      if(N1)       {
        denom     <- matrixStats::rowLogSumExps(numer)
        loglike   <- sum(if(ctrl$do.wts) denom * attr(seqs, "Weights") else denom)
        if(ctrl$do.cv)       {
            return(loglike)
        } else     {
            return(list(loglike = loglike, z = exp(sweep(numer, 1L, denom, FUN="-", check.margin=FALSE))))
        }
      }   else     {
            return(matrixStats::logSumExp(numer) * if(ctrl$do.wts) attr(seqs, "Weights") else 1L)
      }
    },    CEM=     {
      if(N1)       {
        z         <- mclust::unmap(max.col(numer), groups=seq_len(G))
        loglike   <- sum(if(ctrl$do.weights) sweep(z * numer, 1L, attr(seqs, "Weights"), FUN="*", check.margin=FALSE) else z * numer, na.rm=TRUE)
        if(ctrl$do.cv)       {
            return(loglike)
        } else     {
            return(list(loglike = loglike, z = z))
        }
      }   else     {
            return(max(numer) * ifelse(ctrl$do.wts, attr(seqs, "Weights"), 1L))
      }
    })
  }
}

.EM_algorithm     <- function(SEQ, numseq, g, modtype, z, ctrl, gating = NULL, covars = NULL, ll = NULL) {
  algo            <- ctrl$algo
  itmax           <- ctrl$itmax
  st.ait          <- ctrl$stopping == "aitken"
  tol             <- ctrl$tol
  if(!(runEM      <- g > 1))   {
    Mstep         <- .M_step(SEQ, l.meth=modtype, ctrl=ctrl, numseq=numseq)
    ll            <- .E_step(SEQ, l.meth=modtype, ctrl=ctrl, numseq=numseq, params=Mstep)
    j             <- 1L
    ERR           <- FALSE
  } else           {
    ll            <- if(is.null(ll)) c(-Inf, -sqrt(.Machine$double.xmax)) else ll
    j             <- 2L
    emptywarn     <- TRUE
    while(isTRUE(runEM))       {
      j           <- j + 1L
      Mstep       <- .M_step(SEQ, l.meth=modtype, ctrl=ctrl, numseq=numseq, gating=gating, covars=covars, z=z)
      Estep       <- .E_step(SEQ, l.meth=modtype, ctrl=ctrl, numseq=numseq, params=Mstep)
      z           <- Estep$z
      ERR         <- any(is.nan(z))
      if(isTRUE(ERR))            break
      if(isTRUE(emptywarn) && ctrl$warn        &&
         any(matrixStats::colSums2(z)
                  == 0))    {    warning(paste0("\tThere were empty components: ", modtype, " (G=", g, ")\n"), call.=FALSE, immediate.=TRUE)
        emptywarn <- FALSE
      }
      ll          <- c(ll, Estep$loglike)
      if(isTRUE(st.ait))       {
        ait       <- .aitken(ll[seq(j - 2L, j, 1L)])
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

.entropy          <- function(p) {
  p               <- p[p > 0]
    sum(-p * log(p))
}

.fac_to_num       <- function(x) {
  Vseq            <- seq_len(nlevels(x[[1L]]))
    as.data.frame(lapply(x, function(y) { levels(y) <- Vseq; y} ))
}

.lambda_mle       <- function(seqs, params, l.meth = c("CC", "UC", "CU", "UU", "CCN", "UCN", "CUN", "UUN"), ctrl, numseq = NULL) {
  theta           <- params$theta
  P               <- attr(seqs, "P")
  V1V             <- attr(seqs, "V1V")
  lV1             <- attr(seqs, "logV1")
  W               <- attr(seqs, "W")
  l.meth          <- match.arg(l.meth)
  noise           <- attr(seqs, "Noise")
  n.meth          <- is.element(l.meth, c("CU", "UU", "CUN", "UUN"))
  p.meth          <- is.element(l.meth, c("UC", "UU", "UCN", "UUN"))
  if(is.null(numseq) && n.meth)        {
    numseq        <- sapply(seqs, .char_to_num)
  }
  if((G           <- attr(seqs, "G")) == 1L) {
    if(l.meth == "CCN")                {
        return(list(lambda = matrix(0L, nrow=1L, ncol=1L)))
    }
    numer         <- switch(EXPR=l.meth, CC=P,                  CU=1L)
    if(ctrl$do.wts)                    {
      ws          <- attr(seqs, "Weights")
      denom       <- switch(EXPR=l.meth, CC=sum(.dseq(seqs, theta) * ws)/W, CU=rowSums(sweep(numseq != .char_to_num(theta), 2L, ws, FUN="*", check.margin=FALSE))/W)
    } else denom  <- switch(EXPR=l.meth, CC=.dbar(seqs, theta),             CU=rowMeans(numseq      != .char_to_num(theta)))  
  } else           {
    N             <- attr(seqs, "N")
    G0            <- ifelse(noise, attr(seqs, "G0"), G)
    if(is.null(z  <- params$z))  stop("'z' must be supplied if 'G'>1",    call.=FALSE)
    if(p.meth     &&
       is.null(prop      <-
               params$prop))     stop("'prop' must be supplied if 'G'>1", call.=FALSE)
    if(noise)      {
      z           <- z[,-G, drop=FALSE]
      switch(EXPR=l.meth, UCN=, UUN=   {
        prop      <- prop[-G]
      })
    }
    switch(EXPR=l.meth, CU=, UU=, CUN=, UUN= {
      dG          <- lapply(seq_len(G0), function(g) unname(apply(numseq,                 2L, "!=",  .char_to_num(theta[g]))))
      dGp         <- vapply(seq_len(G0), function(g) matrixStats::rowSums2(sweep(dG[[g]], 2L, z[,g], FUN="*", check.margin=FALSE)), numeric(P))
    },             {
      pN          <- P * switch(EXPR=l.meth, CCN=sum(z), W)
      dG          <- vapply(seq_len(G0), function(g) .dseq(seqs, theta[g]), numeric(N))
    })
    numer         <- switch(EXPR=l.meth, CC=,  CCN=pN,       
                                         UC=,  UCN=pN * prop,  
                                        CUN=sum(z), W)
    denom         <- switch(EXPR=l.meth, CC=,  CCN=sum(z * dG),
                                         UC=,  UCN=matrixStats::colSums2(z * dG),
                                         CU=,  CUN=matrixStats::rowSums2(dGp),
                                         UU=,  UUN=sweep(dGp, 2L, prop, FUN="/", check.margin=FALSE))
  }
  lambda          <- ifelse(denom < numer * V1V, lV1 + log(numer - denom) - log(denom), 0L)
  lambda          <- matrix(lambda, nrow=switch(EXPR=l.meth, UC=, UU=, UCN=, UUN=G0, 1L), ncol=switch(EXPR=l.meth, CU=, UU=, CUN=, UUN=P, 1L), byrow=is.element(l.meth, c("UU", "UUN")))
  switch(EXPR=l.meth, UC=,  UU=,
                     UCN=, UUN=           {
    lambda[prop   == 0,] <- Inf
  })
    return(list(lambda = if(noise) rbind(lambda, 0L) else lambda, dG = if(G > 1) dG))
}

.M_step           <- function(seqs, l.meth = c("CC", "UC", "CU", "UU", "CCN", "UCN", "CUN", "UUN"), 
                              ctrl, gating = NULL, covars = NULL, z = NULL, numseq = NULL) {
  if(!is.null(gating))    {
    environment(gating)  <- environment()
  }
  G               <- attr(seqs, "G")
  noise           <- attr(seqs, "Noise")
  if(is.null(numseq)     && isFALSE(ctrl$numseq)) {
    numseq        <- sapply(seqs, .char_to_num)
    attr(numseq,  "G")   <- G
    attr(numseq,  "P")   <- attr(seqs, "P")
  }
  if(G > 1L)       {
    if(is.null(z))               stop("'z' must be supplied when 'G'>1", call.=FALSE)
    if(ctrl$do.wts)       {
      z           <- sweep(z, 1L, attr(seqs, "Weights"), FUN="*", check.margin=FALSE)
    }
    if((gate.g    <- ctrl$gate.g))    {
      prop        <- if(is.element(l.meth, c("UC", "UU", "UCN", "UUN"))) { if(ctrl$do.wts) matrixStats::colSums2(z)/attr(seqs, "W") else colMeans(z) }
      if(!noise   || ctrl$noise.gate) {
        fitG      <- nnet::multinom(gating, trace=FALSE, data=covars, maxit=ctrl$g.itmax, reltol=ctrl$g.tol)
        tau       <- fitG$fitted.values
      } else       {
        zN        <- z
        z         <- z[,-G, drop=FALSE]
        rz        <- row(z)
        z         <- z/matrixStats::rowSums2(z)[rz]
        z[is.nan(z)]     <- .Machine$double.eps
        fitG      <- nnet::multinom(gating, trace=FALSE, data=covars, maxit=ctrl$g.itmax, reltol=ctrl$g.tol)
        tau       <- .tau_noise(fitG$fitted.values, zN[,G], rz)
        z         <- zN
      }
    } else         {
      prop        <- if(isFALSE(ctrl$equalPro) || (noise && !ctrl$equalNoise)) { if(ctrl$do.wts) matrixStats::colSums2(z)/attr(seqs, "W") else colMeans(z) }
      tau         <- if(isFALSE(ctrl$equalPro)) prop else if(noise && !ctrl$equalNoise) c(rep((1 - prop[G])/attr(seqs, "G0"), attr(seqs, "G0")), prop[G]) else rep(1/G, G)
    }
  }   else prop   <- tau <- 1L
  theta           <- .optimise_theta(seqs=seqs, ctrl=ctrl, z=z, numseq=numseq)
  MLE             <- .lambda_mle(seqs=seqs, params=list(theta=theta, z=z, prop=prop), l.meth=l.meth, ctrl=ctrl, numseq=numseq)
  param           <- list(theta=theta, lambda=MLE$lambda, dG=if(G > 1) MLE$dG, tau=tau, fitG=if(G > 1 && gate.g) fitG)
  attr(param, "l.meth")  <- l.meth
    return(param)
}

.mat_byrow        <- function(x, nrow, ncol) {
    matrix(x, nrow=nrow, ncol=ncol, byrow=any(dim(as.matrix(x)) == 1))
}

.modal            <- function(x) {
  ux              <- unique(x[x != 0])
  tab             <- tabulate(match(x, ux))
    ux[which(tab  == max(tab))]
}

.num_to_char      <- function(x) {
    paste(x, sep="", collapse="")
}

.optimise_theta   <- function(seqs, ctrl, z = NULL, numseq = NULL) {
  opti            <- ctrl$opti
  ordering        <- ctrl$ordering
  P               <- attr(seqs, "P")
  nmeth           <- ctrl$nmeth
  G               <- attr(seqs, "G") - nmeth
  Gseq            <- seq_len(G)
  if(G == 0)       {
      return(rep(NA, P))
  }
  if(G     > 1L   && is.null(z))  stop("'z' must be supplied when 'G'>1", call.=FALSE)
  if(opti         == "mode" || ordering != "none")   {
    if(is.null(numseq)      && isFALSE(ctrl$numseq)) {
      numseq      <- sapply(seqs, .char_to_num)
      attr(numseq, "P")     <- P
    }
    attr(numseq, "G")       <- G
    if(opti       == "mode") {
      if(G > 1    || ctrl$do.wts)    {
        theta     <- .weighted_mode(numseq=numseq, z=if(G == 1) as.matrix(attr(seqs, "Weights")) else if(nmeth) z[,Gseq, drop=FALSE] else z)
        if(is.list(theta)   && any(nonu <- apply(theta,   2L, function(x) any(nchar(x) > 1)))) {
          theta[,nonu]      <- lapply(theta[,nonu], "[[", 1L)
        } else       nonu   <- rep(FALSE, G)
        theta     <- apply(theta,  2L, .num_to_char)
      }   else     {
        theta     <- apply(numseq, 1L, .modal)
        if(nonu   <- is.list(theta)) {
          theta   <- lapply(theta, "[[",   1L)
        }
        theta     <- .num_to_char(theta)
      }
      theta       <- if(nmeth) c(theta, NA) else theta
      attr(theta, "NonUnique")   <- nonu
        return(theta)
    }
  }
  
  V               <- attr(seqs, "V")
  theta.opt       <- .theta_data(seqs=seqs, z=z, ctrl=ctrl)
  if(opti == "medoid")       {
      return(if(nmeth) c(theta.opt$theta, NA) else theta.opt$theta)
  }
  pseq            <- seq_len(P)
  vseq            <- seq_len(V)
  opt             <- theta.opt$dsum
  if(G > 1L)       {
    theta         <- lapply(theta.opt$theta, .char_to_num)
  } else           {
    z             <- matrix(if(ctrl$do.wts) attr(seqs, "Weights") else 1L, nrow=attr(seqs, "N"))
    theta         <- list(.char_to_num(theta.opt$theta))
  }
  if(ordering     != "none") {
    stab          <- if(G == 1 && !ctrl$do.wts) list(apply(numseq, 1L, tabulate, V)) else lapply(Gseq, function(g) apply(sweep(numseq, 2L, z[,g], FUN="*", check.margin=FALSE), 1L, tabulate, V))
    sorder        <- lapply(Gseq, function(g)   order(apply(stab[[g]], 2L, .entropy), decreasing=ordering == "decreasing"))
    theta         <- lapply(Gseq, function(g)   theta[[g]][sorder[[g]]])
    seqs          <- lapply(Gseq, function(g)   unname(apply(numseq[sorder[[g]],], 2L, .num_to_char)))
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
                   opts[v]              <- sum(z[,g] * .dseq(seqs[[g]], paste(replace(theta[[g]], p, v), collapse="")))
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
                   opts[p, v]           <- sum(z[,g] * .dseq(seqs[[g]], paste(replace(theta[[g]], p, v), collapse="")))
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
  
  theta           <- switch(EXPR=ordering,
                            none=if(G > 1)  {
                              do.call(base::c, lapply(theta, .num_to_char))
                            } else .num_to_char(theta[[1L]]),
                            if(G > 1)       {
                              do.call(base::c, lapply(Gseq, function(g) .num_to_char(theta[[g]][match(pseq, sorder[[g]])])))
                            } else .num_to_char(theta[[1L]][match(pseq, sorder[[1L]])]))
    return(if(nmeth) c(theta, NA) else theta)
} 

.pick_MEDCrit     <- function(x, pick = 3L) {
  if(!inherits(x, 
     "MEDcriterion"))            stop("'x' must be an object of class 'MEDcriterion'", call.=FALSE)
  x               <- replace(x, !is.finite(x), NA)
  pick            <- min(pick,        length(x[!is.na(x)]))
  decrease        <- !is.element(attr(x, "Criterion"), c("DF", "ITERS", "NEC"))
  x.sx            <- sort(x,          decreasing=decrease)[pick]
  x.crit          <- if(decrease)     x   >= x.sx else x <= x.sx
  x.ind           <- which(x.crit,    arr.ind=TRUE)
  x.val           <- sort(x[x.ind],   decreasing=decrease)
  ind.x           <- order(x[x.ind],  decreasing=decrease)
  x.ind           <- x.ind[ind.x,,    drop=FALSE]
  x.ind[,1L]      <- gsub(".*= ", "", rownames(x)[x.ind[,1L]])
  x.ind[,2L]      <- colnames(x)[as.numeric(x.ind[,2L])]
    return(list(crits = stats::setNames(x.val[seq_len(pick)], vapply(seq_len(pick), function(p, b=x.ind[p,]) paste0(b[2L], ",", b[1L]), character(1L))), pick = pick))
}

.replace_levels   <- function(seq, levels = NULL) {
  seq             <- as.numeric(strsplit(seq, "")[[1L]])
    if(is.null(levels)) factor(seq) else factor(seq, levels=seq_along(levels) - any(seq == 0), labels=as.character(levels))
}

.tau_noise        <- function(tau, z0, row.ind = row(tau)) {
  t0              <- mean(z0)
    cbind(tau/(matrixStats::rowSums2(tau)[row.ind]) * (1 - t0), unname(t0))
}

.theta_data       <- function(seqs, z = NULL, ctrl = NULL) {
  if((G <- attr(seqs, "G"))   == 1L)  {
    sumdist       <- if(ctrl$do.wts)  {
                          vapply(seq_len(attr(seqs, "N")), function(i) sum(.dseq(seqs, seqs[i]) * attr(seqs, "Weights")), numeric(1L))
                   } else vapply(seq_len(attr(seqs, "N")), function(i) sum(.dseq(seqs, seqs[i])), numeric(1L))
      return(list(theta = seqs[which.min(sumdist)], dsum = min(sumdist)))
  } else           {
    if(is.null(z))               stop("'z' must be supplied if 'G'>1",     call.=FALSE)
    if(is.null(ctrl$nmeth))      stop("'nmeth' must be supplied if 'G'>1", call.=FALSE)
    theta         <- dsum     <- list()
    for(g in seq_len(G - ctrl$nmeth)) {
      sumdist     <- vapply(seq_len(attr(seqs, "N")), function(i) sum(.dseq(seqs, seqs[i]) * z[,g]), numeric(1L))
      theta[[g]]  <- seqs[which.min(sumdist)]
      dsum[[g]]   <- min(sumdist)
    }
      return(list(theta = do.call(base::c, theta), 
                  dsum  = do.call(base::c, dsum)))
  }
}

.weighted_mode    <- function(numseq, z) {
    sapply(seq_len(attr(numseq, "G")), function(g)
    sapply(seq_len(attr(numseq, "P")), function(p, x=tapply(z[,g], numseq[p,], sum)) names(which(x == max(x)))))
}
#