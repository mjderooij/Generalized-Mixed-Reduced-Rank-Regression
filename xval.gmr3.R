xval.gmr3 = function(Yn = NULL, Yb = NULL, Yo = NULL, X, Xscale, S = 1:4, K = 10, repeats = 1, mseed = 12){
  # xval function to select the rank/dimensionality
  # have to choose a range of S
  # getting some constants from the data
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
  Ylist = vector(mode = "list", length = ncol(Y))
  for (r in 1:R) {
    if (Rn > 0 && r %in% Nseq) {
      Ylist[[r]] = Y[ , r]      }
    if (Rb > 0 && r %in% Bseq) {
      Ylist[[r]] = Y[, r] 
    }
    if (Ro > 0 && r %in% Oseq) {
      Ylist[[r]] = class.ind(Y[ , r])
    }
  }
  
  # preparation for cross-validation
  set.seed(mseed)
  folds <- cut(seq(1, N), breaks = K, labels = FALSE)
  
  pe.df = data.frame(matrix(nrow = (length(S) * repeats * K), ncol = (3 + R + 1)))
  prederror = rep(NA, R)
  
  colnames(pe.df) = c("S", "repeats", "fold", paste0("PE", 1:R), "PEtotal")
  teller = 0 # this just counts - so that we can put things in the correct rows of pe.df
  
    
  for(rep in 1:repeats){
    folds = sample(folds)
    for(i in 1:length(S)){
      for(k in 1:K){
        teller = teller + 1
        cat("This is analysis", teller, "from a total of", nrow(pe.df), "\n")
        idx = which(folds == k, arr.ind=TRUE) 
        # train
        out = gmr3(Yn = Yn[-idx, ], Yb = Yb[-idx, ], Yo = Yo[-idx, ], X = X[-idx, ], Xscale = Xscale, S = S[i], trace = FALSE)
        # predict
        Yhat =  predict.gmr3(out, newX = X[idx, ])
        # compute deviance of the predictions
        for (r in 1:R) {
          if (Rn > 0 && r %in% Nseq) {
            YY = Ylist[[r]] # vector
            YY = YY[idx] # select test cases
            # prederror[r] = sum((YY - Yhat[[r]])^2)/(2 * out$sigma2) + N/2 * log(sqrt(2 * pi * out$sigma2)) 
            prederror[r] = sum((YY - Yhat[[r]])^2)/(2 * out$sigma2) + N/2 * log(sqrt(2 * pi * out$sigma2))/length(YY) 
          }
          if (Rb > 0 && r %in% Bseq) {
            YY = Ylist[[r]] # vector
            YY = YY[idx] # select test cases
            PI = Yhat[[r]]
            # prederror[r] = -(sum(log(PI[which(YY == 1)])) + sum(log(1 - PI[which(YY == 0)])))
            prederror[r] = -(sum(log(PI[which(YY == 1)])) + sum(log(1 - PI[which(YY == 0)])))/length(YY)
          }
          if (Ro > 0 && r %in% Oseq) {
            PI = Yhat[[r]]
            YY = Ylist[[r]] # matrix
            YY = YY[idx, ] #select test cases
            # prederror[r] = -sum(log(PI[which(YY == 1)]))
            prederror[r] = -mean(log(PI[which(YY == 1)]))
          }
        } #response
        pe.df[teller, ] = c(S[i], rep, k, prederror, sum(prederror))
      } # folds
    } # repeats
  } # ranks/dimensionalities
  return(pe.df)
}
