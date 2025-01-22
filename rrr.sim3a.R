# adaptions from rrpack::rrr.sim3 
# rrr.sim3a => categorization of predictor variables

rrr.sim3a = function (n = 100, p = 30, q.mix = c(5, 20, 5), nrank = 2, intercept = rep(0.5, 30), mis.prop = 0.2) 
{
  q <- sum(q.mix)
  q1 <- q.mix[1]
  q2 <- q.mix[2]
  q3 <- q.mix[3]
  X <- matrix(rnorm(n * p), n, p)
  # ---------------------------------------------------------------------------------------------
  # categorize the predictors
  for(j in 1:p){
    X[ , j] = categorize(X[, j], 5)$xa
  }
  X = scale(X)
  # ---------------------------------------------------------------------------------------------
  u0 <- matrix(nrow = p, ncol = nrank, rnorm(p * nrank, 0, 
                                             1))
  u0 <- qr.Q(qr(u0))
  v0 <- matrix(nrow = q, ncol = nrank, runif(q * nrank, 0.5, 
                                             1) * sample(x = c(1, -1), q * nrank, replace = TRUE))
  C <- u0 %*% t(v0)
  C0 <- rbind(intercept, C)
  X0 <- cbind(1, X)
  MU <- X0 %*% C0
  family <- list(gaussian(), binomial(), poisson())
  familygroup <- c(rep(1, q1), rep(2, q2), rep(3, q3))
  cfamily <- unique(familygroup)
  nfamily <- length(cfamily)
  Y <- matrix(nrow = n, ncol = q, 0)
  sigma <- 1
  if (sum(familygroup == 1) > 0) {
    Y[, familygroup == 1] <- MU[, familygroup == 1] + matrix(nrow = n, 
                                                             ncol = q1, rnorm(n * q1, 0, sigma))
  }
  if (sum(familygroup == 2) > 0) {
    prob <- as.matrix(family[[2]]$linkinv(MU[, familygroup == 2]))
    Y[, familygroup == 2] <- apply(prob, 2, function(a) rbinom(n = n, 
                                                               size = 1, a))
  }
  if (sum(familygroup == 3) > 0) {
    prob <- as.matrix(family[[3]]$linkinv(MU[, familygroup == 3]))
    Y[, familygroup == 3] <- apply(prob, 2, function(a) rpois(n = n, 
                                                              lambda = a))
  }
  N <- n * q
  if (is.numeric(mis.prop) & mis.prop < 1 & mis.prop > 0) {
    N_mis <- round(N * mis.prop)
    index_mis <- sample(1:N, size = N_mis, replace = FALSE)
    Y_mis <- Y
    Y_mis[index_mis] <- NA
  }
  else {
    index_mis <- NULL
    Y_mis <- Y
  }
  list(Y = Y, Y.mis = Y_mis, index.miss = index_mis, X = X, 
       C = C, family = family[q.mix != 0], familygroup = familygroup)
}

