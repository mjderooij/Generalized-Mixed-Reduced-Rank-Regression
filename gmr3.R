# function for generalized mixed reduced rank regression
#' @param Yn is response variable matrix with numeric responses - size N x Rn
#' @param Yb is response variable matrix with binary responses (coded 0/1)- size N x Rb
#' @param Yo is response variable matrix with ordinal responses - size N x Ro
#' @param X is predictor variable matrix - size N x P
#' @param S dimensionality/rank
#' @param Xscale list with optimal scaling level for predictors. Implemented options:
#'           - "N"(numeric)
#'           - "O" (ordinal)
#'           - "C" (nominal) 
#' @param trace whether convergence information needs to be printed on screen
#' @param maxiter maximum number of iterations
#' @param dcrit convergence criterion
# ----------------------------------------------------

gmr3 = function(Yn = NULL, Yb = NULL, Yo = NULL, X, S = 2, Xscale, trace = TRUE, maxiter = 65536, dcrit = 1e-6){
  # load dependencies
  library(MASS) # polr
  library(nnet) # class.ind
  # library(splines2) # for OS-levels MS and NMS
  library(monotone) # for OS-level O
  
  # ----------------------------------------------------
  eps = 1e-16
  sigma2 = 4
  Xori = X
  
  N = nrow(X)
  P = ncol(X)
  if(!is.null(Yn)) if(nrow(Yn) != N) stop("Number of rows of X and Yn should be equal")
  if(!is.null(Yb)) if(nrow(Yb) != N) stop("Number of rows of X and Yb should be equal")
  if(!is.null(Yo)) if(nrow(Yo) != N) stop("Number of rows of X and Yo should be equal")
  ones.n = rep(1, N)
  
  Rn = Rb = Ro = idd = 0
  Nseq = Bseq = Oseq = NULL
  if(!is.null(Yn)){
    Rn = ncol(Yn)
    Nseq = 1:Rn
    idd = Rn
  }  
  
  if(!is.null(Yb)){
    Rb = ncol(Yb)
    Bseq = (idd + 1):(idd + Rb)
    idd = idd + Rb
    Qb = 2 * as.matrix(Yb) - 1 # assumes Yb is binary 0, 1
  }  
  
  if(!is.null(Yo)){
    Ro = ncol(Yo)
    Oseq = (idd + 1):(idd + Ro)
    idd = idd + Ro
    xi = matrix(0, N, Ro) # needed for computing expected value
  }  
  R = Rn + Rb + Ro
  Y = cbind(Yn, Yb, Yo)
  
  # ----------------------------------------------------
  # initalization optimal scaling for X
  originalCategories = quantifications = DG = D = G = list()
  
  for (p in 1:P){
    if(Xscale[p] != "N"){
      originalCategories[[p]] = sort(unique(X[ , p]))
      quantifications[[p]] = as.matrix(sort(unique(scale(as.numeric(X[ , p])))))
      row.names(quantifications[[p]]) = originalCategories[[p]]
      G[[p]] = class.ind(X[, p])
      D[[p]] = t(G[[p]]) %*% G[[p]]
      DG[[p]] = solve(D[[p]], t(G[[p]]))
    }
  }
  
  PHI = X
  
  a = scale(X[, which(Xscale == "N")])
  mx = attr(a, "scaled:center")
  sdx = attr(a, "scaled:scale")
  
  PHI[, which(Xscale == "N")] = a # standardize numeric predictors
  for(p in 1:P){
    if(Xscale[p] != "N") PHI[ , p] = G[[p]] %*% quantifications[[p]]
  }
  
  # starting values
  eig.out = eigen(t(PHI) %*% PHI)
  iRx = eig.out$vectors %*% diag(1/sqrt(eig.out$values)) %*% t(eig.out$vectors)
  iRxX = iRx %*% t(X)
  udv = svd(iRxX %*% scale(Y, center = TRUE, scale = FALSE))
  B = iRx %*% matrix(udv$u[, 1:S], P, S) %*% diag(udv$d[1:S], nrow = S, ncol = S)
  V = matrix(udv$v[, 1:S], R, S) 
  # B = matrix(0, P, S); V = matrix(0, R, S) # zero starting values
  theta = PHI %*% B %*% t(V) # structural part
  
  mm = rep(0, R) # same as elements of m[[]] for numeric/binary with zeros for ordinal
  m = list() # intercepts for each response variable (where for ordinal there are multiple thresholds)
  
  dev = rep(NA, R)
  rdevs.n  = rdevs.b = rdevs.o = c()
  
  if(!(is.null(Yn))) mm[Nseq] = colMeans(Yn - PHI %*% B %*% t(V[Nseq, , drop = FALSE]))
  if(!(is.null(Yb))) mm[Bseq] = colMeans(Yb) # Zb - PHI %*% B %*% t(V[Bseq, ])
  if(!is.null(Yo)){
    for(r in 1:Ro){
      polr.out = polr(as.factor(Yo[ , r]) ~ 1 + offset(theta[, Oseq[r]]))
      m[[Oseq[r]]] = polr.out$zeta
      dev[Oseq[r]] = polr.out$deviance
    }
  }
  
  eta = outer(ones.n, mm) + theta # "linear predictor"
  
  # ----------------------------------------------------
  dev.n = dev.b = dev.o = 0 #initialisation - stays 0 when either is not available. 
  if(!is.null(Yn)){
    RES = Yn - eta[ , Nseq]
    sigma2 = sd(RES)^2
    # sigma2b = sum(diag(t(RES) %*% RES))/(N*Rn)
    dev.n = sum(RES^2)/(2 * sigma2) + N*Rn/2 * log(sqrt(2 * pi * sigma2)) 
  } # numeric part
  if(!is.null(Yb)) dev.b = -sum(log(plogis(Qb * theta[, Bseq]))) # binary part
  if(!is.null(Yo)) dev.o = sum(dev[Oseq])/2  # ordinal part
  
  loss.old = dev.n + dev.b + dev.o
  
  # ----------------------------------------------------
  # Working responses matrix Z
  Zn = Zb = Zo = NULL
  if(!is.null(Yn)) Zn = Yn
  
  for (iter in 1:maxiter){
    
    # ----------------------------------------------------
    # update majorization constant
    # ----------------------------------------------------
    kappa = max(1/sigma2, 0.25) 
    ikappa = 1/kappa
    
    # ----------------------------------------------------
    # update working responses
    # ----------------------------------------------------
    # continuous
    if(!is.null(Yn)) Zn = eta[ , Nseq] - ikappa * (-RES/sigma2)
    # ------------------------------
    # binary
    if(!is.null(Yb)) Zb = as.matrix(eta[, Bseq] + ikappa * Qb * (1 - plogis(Qb * eta[, Bseq])))
    # if(!is.null(Yb)) Zb = as.matrix(eta[, Bseq] + 4 * Qb * (1 - plogis(Qb * eta[, Bseq])))
    # ------------------------------
    #ordinal
    if(!is.null(Yo)){
      for(r in 1:Ro){
        Ep = expected.p(Yo[ , r], theta[ , Oseq[r] ], m[[ Oseq[r] ]])
        xi[ , r] = 1 - 2 * Ep
      }
      Zo = as.matrix(theta[, Oseq] - ikappa * xi)
      # Zo = as.matrix(theta[, Oseq] - 4 * xi)
    } 
    # ------------------------------
    Z = cbind(Zn, Zb, Zo)
    
    # ----------------------------------------------------
    # update transformations
    # ----------------------------------------------------
    BV = B %*% t(V)
    for(p in 1:P){
      if(Xscale[p] != "N"){
        U = Z - PHI[ , -p, drop = FALSE] %*% BV[-p, ]
        
        # Calculate the nominal (unrestricted) solution to the system
        quantifications[[p]] = DG[[p]] %*% U %*% t(BV[p, , drop = FALSE]) %*% solve(crossprod(BV[p, ]))
        
        # restrictions
        if(Xscale[p] == "O"){ #if predictors are ordinal
          lp = length(quantifications[[p]])
          wp = as.numeric(diag(D[[p]]))
          # monotone increasing
          qq1 = monotone(quantifications[[p]], w = wp)
          l1 = sum(wp * (quantifications[[p]] - qq1)^2)
          # monotone decreasing
          qq2 = monotone(quantifications[[p]][lp:1], w = wp[lp:1])[lp:1]
          l2 = sum(wp * (quantifications[[p]] - qq2)^2)
          # best fitting wins
          if(l1 < l2){
            quantifications[[p]] = qq1
          }
          else{
            quantifications[[p]] = qq2
          }
        }
        
        # Standardize the quantifications
        quantifications[[p]] = matrix(scale(rep(quantifications[[p]], times = diag(D[[p]])))[cumsum(c(1, diag(D[[p]])))[1:(nrow(D[[p]]))]], ncol = 1)
        PHI[, p] = G[[p]] %*% quantifications[[p]]
      }
    }
    
    A = kronecker(diag(S), t(PHI) %*% PHI)
    
    # ----------------------------------------------------
    # update B
    # ----------------------------------------------------
    Zc = (Z - outer(ones.n, mm))
    vz = Vec(Zc)
    b = solve(A, t(kronecker(V, PHI)) %*% vz)
    B = matrix(b, ncol = S)
    
    # ----------------------------------------------------
    # update V
    # ----------------------------------------------------
    udv = svd(t(B) %*% t(PHI) %*% (Zc))
    V = t(udv$u[, 1:S] %*% t(udv$v[, 1:S]))    
    
    theta = PHI %*% B %*% t(V)
    
    # ----------------------------------------------------
    # update m
    if(!(is.null(Yn))) mm[Nseq] = colMeans(Zn - theta[ , Nseq])
    if(!(is.null(Yb))) mm[Bseq] = colMeans(Zb - theta[ , Bseq])
    if(!is.null(Yo)){
      for(r in 1:Ro){
        rr = Oseq[r]
        polr.out = polr(as.factor(Yo[ , r]) ~ 1 + offset(theta[, rr]))
        m[[rr]] = polr.out$zeta
        dev[rr] = polr.out$deviance
      }
    }
    
    eta = outer(ones.n, mm) + theta # PHI %*% B %*% t(V)
    
    # ----------------------------------------------------
    # compute loss and determine convergence
    # ----------------------------------------------------

    # NOTE: The devs computed here are actually negative log-likelihoods!!
    if(!is.null(Yn)){
      RES = Yn - (outer(ones.n, mm[Nseq]) + PHI %*% B %*% t(V[Nseq, , drop = FALSE]))
      sigma2 = sd(RES)^2
      # sigma2 = sum(diag(t(RES) %*% RES))/(N*Rn)
      dev.n = sum(RES^2)/(2 * sigma2) + N*Rn/2 * log(sqrt(2 * pi * sigma2))  # numeric part
      rdevs.n = colSums(RES^2)/(2 * sigma2) + N/2 * log(sqrt(2 * pi * sigma2))
    } 
    if(!is.null(Yb)){
      dev.b =  -sum(log(plogis(Qb * eta[, Bseq]))) # binary part
      rdevs.b = -colSums(log(plogis(Qb * eta[, Bseq])))
    } 
    if(!is.null(Yo)){
      dev.o = sum(dev[Oseq])/2  # ordinal part
      rdevs.o = dev[Oseq]/2
    }
    
    loss.new = dev.n + dev.b + dev.o
    
    if (trace) {cat(iter, loss.old, loss.new, (2*(loss.old - loss.new)/(loss.old + loss.new)), "\n")}
    # if (loss.new > loss.old) stop("Divergence - should not occur!")
    if ( (2*(loss.old - loss.new)/(loss.old + loss.new)) < dcrit ) break
    if (iter == maxiter) warning("Maximum number of iterations reached - not converged (yet)")
    loss.old = loss.new
  }
  
  xnames = colnames(X)
  ynames = c(colnames(Yn), colnames(Yb), colnames(Yo))
  
  # number of parameters
  npar = (P + R - S)*S + Rn + Rb + sum(sapply(m, length)) + sum(unlist(sapply(G, ncol)) - 2)
  npar.change = (P + R - S)*S + sum(unlist(sapply(G, ncol)) - 2)
  
  # make output object
  output = list(
    Xoriginal = Xori,
    Xscale = Xscale,
    xnames = xnames,
    ynames = ynames,
    mx = mx, 
    sdx = sdx,
    Y = Y,
    G = G,
    Yscale = list(Nseq = Nseq, Bseq = Bseq, Oseq = Oseq),
    R = c(Rn, Rb, Ro),
    Yn = Yn,
    Yb = Yb,
    Yo = Yo,
    Zn = Zn,
    Zb = Zb,
    Zo = Zo,
    PHI = PHI, 
    sigma2 = sigma2,
    m = m,
    mm = mm,
    B = B,
    V = V,
    originalcategories = originalCategories,
    quantifications = quantifications,
    iter = iter,
    loss = loss.new,
    rloss = c(rdevs.n, rdevs.b, rdevs.o),
    npar = npar,
    npar.change = npar.change, # the increase in npar compared to intercept only model 
    AIC = 2 * loss.new + 2 * npar,
    BIC = 2 * loss.new + log(N) * npar
  )
  # ----------------------------------------------------
  return(output)
}

