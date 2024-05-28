predict.gmr3 <- function(object, newX) {
  
  Rn <- Rb <- Ro <- idd <- 0
  
  if (!is.null(object$Yn)) {
    Rn <- ncol(object$Yn)
    Nseq <- 1:Rn
    idd <- Rn
  }  
  
  if (!is.null(object$Yb)) {
    Rb <- ncol(object$Yb)
    Bseq <- (idd + 1):(idd + Rb)
    idd = idd + Rb
  }  
  
  if (!is.null(object$Yo)) {
    Ro <- ncol(object$Yo)
    Oseq <- (idd + 1):(idd + Ro)
    idd = idd + Ro
  }  
  
  
  newPHI = newX
  P = ncol(newX)
  
  newPHI[, which(object$Xscale == "N")] = scale(newX[, which(object$Xscale == "N")], 
                                                center = object$mx, scale = object$sdx)
  
  for (p in 1: P) {
    if (Xscale[p] != "N") {
      newPHI[, p]  = myrecode(newX[, p], object$originalcategories[[p]], object$quantifications[[p]])
    }
  }
  
  
  structural_theta <- newPHI %*% (object$B %*% t(object$V))
  R = ncol(structural_theta)
  Yhat <- vector(mode = "list", length = R)
  
  
  for (r in 1:R) {
    
    if (Rn > 0 && r %in% Nseq) {
      theta <- outer(rep(1, nrow(newX)), object$mm[Nseq]) + structural_theta[, Nseq]
      Yhat[[r]] = theta[,r]
    }
    
    if (Rb > 0 && r %in% Bseq) {
      theta <- outer(rep(1, nrow(newX)), object$mm[Bseq]) + structural_theta[, Bseq]
      Yhat[[r]] =  plogis(theta[,r])
      
    }
    
    if (Ro > 0 && r %in% Oseq) {
      zeta = c(-Inf, object$m[[r]], Inf)
      Cr = length(zeta)
      P = plogis(outer(rep(1,nrow(newX)), zeta) - outer(structural_theta[, r], rep(1, Cr)))
      Yhat[[r]] = t(apply(P, 1, diff)) 
    }
    
  }
  
  return(Yhat) 
}