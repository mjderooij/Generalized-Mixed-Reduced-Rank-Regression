# READ in DATA
source("~/surfdrive/Lorenza Cotugno/read_data.R")
mydata = read_data()
Yo = mydata$Yo
Yb = mydata$Yb
X = mydata$X
Xscale <- c("N", "O","C", "O", "O")

# SOURCE functions

source("~/surfdrive/Lorenza Cotugno/new/plot.boot.gmr3.R")
source("~/surfdrive/Lorenza Cotugno/new/gmr3.R")
source("~/surfdrive/Lorenza Cotugno/new/boot.gmr3.R")
source("~/surfdrive/Lorenza Cotugno/new/plot.boot.gmr3.R")
source("~/surfdrive/Lorenza Cotugno/new/auxiliary.gmr3.R")
source("~/surfdrive/Lorenza Cotugno/new/fit.gmr3.R")

# recoding Yo en Yb
YYo = matrix(NA, 837, 5)
for(r in 1:5){YYo[, r] = myrecode(Yo[, r], old = c(1,2,3,4), new = c(4,3,2,1))}

YYb = matrix(NA, 837, 2)
for(r in 1:2){YYb[, r] = myrecode(Yb[, r], old = c(0,1), new = c(1, 0))}

rm(Yb, Yo)

colnames(YYo) = c("CI", "MW", "FS", "DI", "RE")
colnames(YYb) = c("T", "FE")
colnames(X) = c("A", "PA", "G", "U", "E")

# ANALYZE the DATA
out1 = gmr3(Yb = YYb, Yo = YYo, X = X, Xscale = Xscale, S = 1)
out2 = gmr3(Yb = YYb, Yo = YYo, X = X, Xscale = Xscale, S = 2)
out3 = gmr3(Yb = YYb, Yo = YYo, X = X, Xscale = Xscale, S = 3)
out4 = gmr3(Yb = YYb, Yo = YYo, X = X, Xscale = Xscale, S = 4)
out5 = gmr3(Yb = YYb, Yo = YYo, X = X, Xscale = Xscale, S = 5)

fit1 = fit.gmr3(out1)
fit2 = fit.gmr3(out2)
fit3 = fit.gmr3(out3)
fit4 = fit.gmr3(out4)
fit5 = fit.gmr3(out5)

source("~/surfdrive/Lorenza Cotugno/new/predict.gmr3.R")
source("~/surfdrive/Lorenza Cotugno/new/xval.gmr3.R")
source("~/surfdrive/Lorenza Cotugno/new/plot.xval.gmr3.R")

PE = xval.gmr3(Yb = YYb, Yo = YYo, X = X, Xscale = Xscale, S = 1:5, repeats = 10, mseed = 12)
PE2 = plot.xval.gmr3(PE)

dimsel = matrix(NA, 5, 7)
dimsel[, 1] = 1:5
dimsel[, 2] = c(out1$loss, out2$loss, out3$loss, out4$loss, out5$loss)
dimsel[, 3] = c(out1$npar, out2$npar, out3$npar, out4$npar, out5$npar)
dimsel[, 4] = c(out1$AIC, out2$AIC, out3$AIC, out4$AIC, out5$AIC)
dimsel[, 5] = c(out1$BIC, out2$BIC, out3$BIC, out4$BIC, out5$BIC)
dimsel[, 6] = c(fit1$McFadden, fit2$McFadden, fit3$McFadden, fit4$McFadden, fit5$McFadden)
dimsel[, 7] = c(fit1$McFadden.adj, fit2$McFadden.adj, fit3$McFadden.adj, fit4$McFadden.adj, fit5$McFadden.adj)
colnames(dimsel) = c("S", "-LL", "k", "AIC", "BIC", "R^2", "R_a^2")

library(xtable)
xtable(dimsel)
out = out2

b.out = boot.gmr3(out)
plts = plot.boot.gmr3(b.out)