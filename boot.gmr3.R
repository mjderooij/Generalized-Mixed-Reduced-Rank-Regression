# function for balanced bootstrap of a gmr3-objects
# takes Bsample bootstrap samples

boot.gmr3 = function(object, Bsamples = 1000){
  
  # get data from object
  Yn = object$Yn
  Yb = object$Yb
  Yo = object$Yo
  Xscale = object$Xscale
  X = object$Xoriginal
  N = nrow(X)
  P = nrow(object$B)
  S = ncol(object$B)
  R = nrow(object$V)
  
  start = list(B = object$B, 
               V = object$V,
               m = object$m,
               mm = object$mm,
               quantifications = object$quantifications)

  # bootstrap samples given
  if(is.matrix(Bsamples)){f = Bsamples; Bsamples = ncol(f)}
  # else create balanced bootstrap samples 
  else{
    f = matrix(1:N, N, Bsamples)
    ff = matrix(f,prod(dim(f)),1)
    fff = sample(ff)
    f = matrix(fff, N, Bsamples)
  }
  
  # create empty matrices for bootstrap estimates
  BB = matrix(NA, P*S, Bsamples)
  BV = matrix(NA, R*S, Bsamples)
  Bm = matrix(NA, length(unlist(object$m)), Bsamples)
  Bmm = matrix(NA, length(object$mm), Bsamples)
  BBdf = matrix(NA, P*Bsamples, (S + 2))
  BVdf = matrix(NA, R*Bsamples, (S + 2))
  
  # run bootstraps
  for (b in 1:Bsamples) {
    cat("This is analysis", b, "from a total of", Bsamples, "Bootstraps", "\n")
    obs <- f[ , b]
    repout = bgmr3(Yb = Yb[obs, ], Yo = Yo[obs, ], X = X[obs, ], S = S, Xscale = Xscale, start = start)
    # bootstrap estimates weights
    BB[, b] = matrix(repout$B, ncol = 1)
    BBdf[((b-1)*P + 1):(b*P), 1] = b
    BBdf[((b-1)*P + 1):(b*P), 2] = 1:P
    BBdf[((b-1)*P + 1):(b*P), 3:(S+2)] = repout$B
    
    # bootstrap estimates loadings
    BV[, b] = matrix(repout$V, ncol = 1)
    BVdf[((b-1)*R + 1):(b*R), 1] = b
    BVdf[((b-1)*R + 1):(b*R), 2] = 1:R
    BVdf[((b-1)*R + 1):(b*R), 3:(S+2)] = repout$V
    
    # bootstrap estimates intercepts
    Bmm[, b] = repout$mm
    Bm[ , b] = unlist(repout$m)
  }
  
  BBdf = as.data.frame(BBdf)
  colnames(BBdf) = c("Bootstrap", "Predictor", paste0("dim", 1:S))
  BBdf$Predictor = factor(BBdf$Predictor, levels = 1:P, labels = object$xnames)
  
  BVdf = as.data.frame(BVdf)
  colnames(BVdf) = c("Bootstrap", "Response", paste0("dim", 1:S))
  BVdf$Response = factor(BVdf$Response, levels = 1:R, labels = object$ynames)

  output = list(
    gmr3obj = object,
    BB = BB,
    BBdf = BBdf,
    BV = BV,
    BVdf = BVdf,
    Bmm = Bmm,
    Bm = Bm
  )
}



  