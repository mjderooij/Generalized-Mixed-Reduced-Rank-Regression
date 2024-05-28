fit.gmr3 = function(object){
  # get response matrices from object
  Yn = object$Yn
  Yb = object$Yb
  Yo = object$Yo

  R = sum(out$R) # number of response variables
  LL0 = rep(NA, R)
  rr = 0

  # Fit NULL models
  if(!is.null(Yn)){
    for(r in 1:ncol(Yn)){
      rr = rr + 1
      LL0[rr] = glm(Yn[, r] ~ 1)$dev
    }
  }
  if(!is.null(Yb)){
    for(r in 1:ncol(Yb)){
      rr = rr + 1
      LL0[rr] = glm(Yb[, r] ~ 1, family = binomial)$dev/2
    }
  }
  if(!is.null(Yo)){
    for(r in 1:ncol(Yo)){
      rr = rr + 1
      LL0[rr] = polr(as.factor(Yo[, r]) ~ 1)$dev/2
    }
  }
  
  # Compare to fitted model
  loss = -object$loss
  LL0sum = -sum(LL0)
  
  McFadden = 1 - loss/LL0sum
  McFadden.adj = 1 - (loss - object$npar.change)/LL0sum
  
  fit = list(
    LL0sum = LL0sum,
    LLmsum = loss,
    LL0 = LL0,
    LLm = object$rloss,
    McFadden = McFadden,
    McFadden.adj = McFadden.adj,
    McFadden.item = 1 - object$rloss/LL0
  )
  return(fit)
}