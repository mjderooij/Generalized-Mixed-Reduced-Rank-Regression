# ----------------------------------------------------
# ----------------------------------------------------
# FUNCTIONS, LIBRARIES, SOURCE code
# ----------------------------------------------------
# ----------------------------------------------------
# categorize = function(x, ncat){
#   # categorizes a continuous variable into a 2, 3, or 4 categorical variable
#   # uses random cutpoints
#   if(ncat == 2){
#     breaks = quantile(x, probs = runif(1, 0.3, 0.7))
#     xx = ifelse(x < breaks[1], 1, 2)
#   }
#   else if(ncat == 3){
#     breaks = quantile(x, probs = c(runif(1, 0.1, 0.4), runif(1, 0.6, 0.9)))
#     xx = ifelse(x < breaks[1], 1, 
#                 ifelse(x < breaks[2], 2, 3))
#   }
#   else if(ncat == 4){
#     breaks = quantile(x, probs = c(runif(1, 0.1, 0.3), runif(1, 0.4, 0.6), runif(1, 0.7, 0.9)))
#     xx = ifelse(x < breaks[1], 1, 
#                 ifelse(x < breaks[2], 2, 
#                        ifelse(x < breaks[3], 3, 4)))
#   }
#   return(xx)
# } 

library(ggplot2)

categorize = function(x, ncat){
  mybreaks = quantile(x, probs = seq(0, 1, length = (ncat+1)))
  # xx gives the ordered scores, 1, 2, 3
  xx = cut(x, breaks = mybreaks, labels = FALSE, include.lowest = T)
  # breaks = breaks[-1]
  # xx = rep(NA, length(x))
  # for(c in 2:length(breaks)){
  #   if((breaks[c-1] <= x) & (x < breaks[c])){xx = c-1}
  # }
  # xxx gives the averages of x per category of xx
  xxx = xx
  for(c in 1:ncat) xxx[xx == c] = mean(x[xx == c])
  output = list(xo = xx, xa = xxx)
  return(output)
}

# ----------------------------------------------------
# ----------------------------------------------------
# Constants
# ----------------------------------------------------
# ----------------------------------------------------

N = 500
P = 8 # predictors
Rn = 4 # numeric variables
Rb = 4 # binary variables
Ro = 4 # ordinal variables
S = 2 # rank

#------------------------------------------------------------------
#------------------------------------------------------------------
# 1 simulation that compares all
#------------------------------------------------------------------
#------------------------------------------------------------------

source("~/surfdrive/LogitMDA/gmr4/lorenza/paper1-gmr3/new/gmr3.R")
source("~/surfdrive/LogitMDA/gmr4/lorenza/paper1-gmr3/new/auxiliary.gmr3.R")
source("~/surfdrive/LogitMDA/gmr4/lorenza/paper1-gmr3/rrr.sim3a.R")
source("~/surfdrive/LogitMDA/gmr4/lorenza/paper1-gmr3/rrr.sim3b.R")
set.seed(1234)
replications = 250
measures = data.frame(replication = integer(), 
                      method = integer(), 
                      Rmse = double())
teller = 0
for(r in 1:replications){
  cat("This is replication:", r, "\n")
  teller = teller + 1
  mydata = rrpack::rrr.sim3(n = N, p = P, q.mix = c(Rn,Rb,0), nrank = S, intercept = rep(0, (Rn + Rb)), mis.prop = 0)
  BVpop = mydata$C
  Y = mydata$Y; colnames(Y) = paste0("Y", 1:(Rn + Rb))
  X = scale(mydata$X); colnames(X) = paste0("X", 1:P)
  
  # fit met rrpack::mrrr
  cat("This is replication:", r, "met mrrr", "\n")
  
  family <- mydata$family
  familygroup <- mydata$familygroup
  fit.mrrr <- rrpack::mrrr(Y, X, family = family, familygroup = familygroup,
                           penstr = list(penaltySVD = "rankCon", lambdaSVD = 2))
  BVhat = coef(fit.mrrr)[-1, ]
  measures[teller, 1] = r
  measures[teller, 2] = 1
  measures[teller, 3] = sqrt(mean((BVpop - BVhat)^2))
  
  # fit met gmr3
  cat("This is replication:", r, "met gmr3", "\n")
  teller = teller + 1
  fit.gmr3 = gmr3(Yn = Y[ , 1:Rn], Yb = Y[, (Rn + 1):(Rn + Rb)], X = X, S = 2, Xscale = rep("N", P), trace= F)
  BVhat = fit.gmr3$B %*% t(fit.gmr3$V)
  measures[teller, 1] = r
  measures[teller, 2] = 2
  measures[teller, 3] = sqrt(mean((BVpop - BVhat)^2))
  
  # data generatie met ordinale predictoren
  cat("This is replication:", r, "met gmr3-p", "\n")
  teller = teller + 1
  mydata = rrr.sim3a(n = N, p = P, q.mix = c(Rn,Rb,0), nrank = S, intercept = rep(0, (Rn + Rb)), mis.prop = 0)
  BVpop = mydata$C
  Y = mydata$Y; colnames(Y) = paste0("Y", 1:(Rn + Rb))
  X = scale(mydata$X); colnames(X) = paste0("X", 1:P)
  
  for(p in 1:P){X[, p] = categorize(X[, p], 5)$xo} # 4 categories -- vary??
  fit.gmr3 = gmr3(Yn = Y[ , 1:Rn], Yb = Y[, (Rn + 1):(Rn + Rb)], X = X, S = 2, Xscale = rep("O", P), trace= F)
  BVhat = fit.gmr3$B %*% t(fit.gmr3$V)
  measures[teller, 1] = r
  measures[teller, 2] = 3
  measures[teller, 3] = sqrt(mean((BVpop - BVhat)^2))

  # data generatie met binaire en ordinale responses
  cat("This is replication:", r, "met gmr3-r1", "\n")
  teller = teller + 1
  mydata = rrr.sim3b(n = N, p = P, q.mix = c(0,Rb,Ro), nrank = S, intercept = rep(0, (Rb + Ro)), mis.prop = 0)
  BVpop = mydata$C
  Y = mydata$Y; colnames(Y) = paste0("Y", 1:(Rn + Rb))
  X = scale(mydata$X); colnames(X) = paste0("X", 1:P)
  fit.gmr3 = gmr3(Yb = Y[, 1:Rb], Yo = Y[ , (Rb+1):(Rb + Ro)], X = X, S = 2, Xscale = rep("N", P), trace= F)
  BVhat = fit.gmr3$B %*% t(fit.gmr3$V)
  measures[teller, 1] = r
  measures[teller, 2] = 4
  measures[teller, 3] = sqrt(mean((BVpop - BVhat)^2))

  # data generatie met numerieke en ordinale responses
  cat("This is replication:", r, "met gmr3-r2", "\n")
  teller = teller + 1
  mydata = rrr.sim3b(n = N, p = P, q.mix = c(Rn,0,Ro), nrank = S, intercept = rep(0, (Rn + Ro)), mis.prop = 0)
  BVpop = mydata$C
  Y = mydata$Y; colnames(Y) = paste0("Y", 1:(Rn + Rb))
  X = scale(mydata$X); colnames(X) = paste0("X", 1:P)
  fit.gmr3 = gmr3(Yn = Y[, 1:Rn], Yo = Y[ , (Rn+1):(Rn + Ro)], X = X, S = 2, Xscale = rep("N", P), trace= F)
  BVhat = fit.gmr3$B %*% t(fit.gmr3$V)
  measures[teller, 1] = r
  measures[teller, 2] = 5
  measures[teller, 3] = sqrt(mean((BVpop - BVhat)^2))
}

measures$Method = factor(measures$method, labels = c("mrrr", "gmr3", "gmr3-p", "gmr3-r-1", "gmr3-r-2"))
myplot = ggplot(data = measures, aes(x = Method, y = Rmse)) + geom_boxplot() + theme_bw()

ggsave("~/surfdrive/LogitMDA/gmr4/lorenza/paper1-gmr3/simresults.pdf", plot = myplot, width = 297, height = 210, units = "mm", limitsize = FALSE)


# ----------------------------------------------------
# ----------------------------------------------------
# DATA GENERATION AND ANALYSIS
# ----------------------------------------------------
# ----------------------------------------------------
# source("~/surfdrive/LogitMDA/gmr4/lorenza/paper1-gmr3/new/gmr3.R")
# source("~/surfdrive/LogitMDA/gmr4/lorenza/paper1-gmr3/new/auxiliary.gmr3.R")
# 
# set.seed(1234)
# replications = 500
# measures = data.frame(replication = integer(), 
#                       method = integer(), 
#                       Rmse = double())
# teller = 0
# for(r in 1:replications){
#   cat("This is replication:", r, "\n")
#   teller = teller + 1
#   mydata = rrpack::rrr.sim3(n = N, p = P, q.mix = c(Rn,Rb,0), nrank = S, intercept = rep(0, (Rn + Rb)), mis.prop = 0)
#   BVpop = mydata$C
#   Y = mydata$Y; colnames(Y) = paste0("Y", 1:(Rn + Rb))
#   X = scale(mydata$X); colnames(X) = paste0("X", 1:P)
#   
#   # fir met rrpack::mrrr
#   family <- mydata$family
#   familygroup <- mydata$familygroup
#   fit.mrrr <- rrpack::mrrr(Y, X, family = family, familygroup = familygroup,
#                            penstr = list(penaltySVD = "rankCon", lambdaSVD = 2))
#   BVhat1 = coef(fit.mrrr)[-1, ]
#   measures[teller, 1] = r
#   measures[teller, 2] = 1
#   measures[teller, 3] = sqrt(mean((BVpop - BVhat1)^2))
#   
#   # fir met gmr3
#   teller = teller + 1
#   fit.gmr3 = gmr3(Yn = Y[ , 1:Rn], Yb = Y[, (Rn + 1):(Rn + Rb)], X = X, S = 2, Xscale = rep("N", P), trace= F)
#   BVhat2 = fit.gmr3$B %*% t(fit.gmr3$V)
#   measures[teller, 1] = r
#   measures[teller, 2] = 2
#   measures[teller, 3] = sqrt(mean((BVpop - BVhat2)^2))
# }
# 
# measures$Method = factor(measures$method, labels = c("mrrr", "gmr3"))
# ggplot(data = measures, aes(x = Method, y = Rmse)) + geom_boxplot() + theme_bw()

#------------------------------------------------------------------
#------------------------------------------------------------------
# repeating the experiment above, 
# where we categorize the numeric responses to ordinal
# theoretically this should not influence the B and V - ??
# compare 2 gmr3's
#------------------------------------------------------------------
#------------------------------------------------------------------

# source("~/surfdrive/LogitMDA/gmr4/lorenza/paper1-gmr3/rrr.sim3b.R")
# 
# measures = data.frame(replication = integer(), 
#                       responses = integer(), 
#                       Rmse = double())
# teller = 0
# for(r in 1:replications){
#   cat("This is replication:", r, "\n")
#   teller = teller + 1
#   mydata = rrr.sim3b(n = N, p = P, q.mix = c(0,Rb,Rn), nrank = S, intercept = rep(0, (Rn + Rb)), mis.prop = 0)
#   BVpop = mydata$C
#   Y = mydata$Y; colnames(Y) = paste0("Y", 1:(Rn + Rb))
#   X = scale(mydata$X); colnames(X) = paste0("X", 1:P)
#   
#   fit.gmr3 = gmr3(Yn = Y[ , 1:Rn], Yb = Y[, (Rn + 1):(Rn + Rb)], X = X, S = 2, Xscale = rep("N", P), trace= F)
#   BVhat1 = fit.gmr3$B %*% t(fit.gmr3$V)
#   measures[teller, 1] = r
#   measures[teller, 2] = 1
#   measures[teller, 3] = sqrt(mean((BVpop - BVhat1)^2))
#   
#   # with categorical responses 
#   teller = teller + 1
#   # categorize numeric to ordinal
#   fit.gmr3 = gmr3(Yb = Y[, 1:Rb], Yo = Y[ , (Rb+1):(Rb + Rn)], X = X, S = 2, Xscale = rep("N", P), trace= F)
#   BVhat2 = fit.gmr3$B %*% t(fit.gmr3$V)
#   measures[teller, 1] = r
#   measures[teller, 2] = 2
#   measures[teller, 3] = sqrt(mean((BVpop - BVhat2)^2))
# }
# 
# measures$Responses = factor(measures$responses, labels = c("Numeric", "Ordinal"))
# ggplot(data = measures, aes(x = Responses, y = Rmse)) + geom_boxplot() + theme_bw()


#------------------------------------------------------------------
#------------------------------------------------------------------
# repeating the experiment above, 
# where we categorize the numeric predictors to ordinal
# do not know whether this influences B and V
# compare 2 gmr3's
#------------------------------------------------------------------
#------------------------------------------------------------------

# source("~/surfdrive/LogitMDA/gmr4/lorenza/paper1-gmr3/rrr.sim3a.R")
# 
# measures = data.frame(replication = integer(), 
#                       predictors = integer(), 
#                       Rmse = double())
# teller = 0
# for(r in 1:replications){
#   cat("This is replication:", r, "\n")
#   teller = teller + 1
#   mydata = rrr.sim3a(n = N, p = P, q.mix = c(Rn,Rb,0), nrank = S, intercept = rep(0, (Rn + Rb)), mis.prop = 0)
#   BVpop = mydata$C
#   Y = mydata$Y; colnames(Y) = paste0("Y", 1:(Rn + Rb))
#   X = scale(mydata$X); colnames(X) = paste0("X", 1:P)
#   
#   fit.gmr3 = gmr3(Yn = Y[ , 1:Rn], Yb = Y[, (Rn + 1):(Rn + Rb)], X = X, S = 2, Xscale = rep("N", P), trace= F)
#   BVhat1 = fit.gmr3$B %*% t(fit.gmr3$V)
#   measures[teller, 1] = r
#   measures[teller, 2] = 1
#   measures[teller, 3] = sqrt(mean((BVpop - BVhat1)^2))
#   
#   # with categorical predictors 
#   teller = teller + 1
#   for(p in 1:P){X[, p] = categorize(X[, p], 5)$xo} # 4 categories -- vary??
#   fit.gmr3 = gmr3(Yn = Y[ , 1:Rn], Yb = Y[, (Rn + 1):(Rn + Rb)], X = X, S = 2, Xscale = rep("O", P), trace= F)
#   BVhat2 = fit.gmr3$B %*% t(fit.gmr3$V)
#   measures[teller, 1] = r
#   measures[teller, 2] = 2
#   measures[teller, 3] = sqrt(mean((BVpop - BVhat2)^2))
# }
# 
# measures$Predictors = factor(measures$predictors, labels = c("Numeric", "Ordinal"))
# ggplot(data = measures, aes(x = Predictors, y = Rmse)) + geom_boxplot() + theme_bw()
# 
