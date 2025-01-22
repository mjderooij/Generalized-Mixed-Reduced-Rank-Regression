# READ in DATA
setwd("~/surfdrive/LogitMDA/gmr4/lorenza/paper1-gmr3")
source("~/surfdrive/LogitMDA/gmr4/lorenza/paper1-gmr3/read_data.R")
mydata = read_data()
Yo = mydata$Yo
Yb = mydata$Yb
X = mydata$X
Xscale <- c("N", "O","C", "O", "O")

# SOURCE functions

source("~/surfdrive/LogitMDA/gmr4/lorenza/paper1-gmr3/new/plot.boot.gmr3.R")
source("~/surfdrive/LogitMDA/gmr4/lorenza/paper1-gmr3/new/gmr3.R")
source("~/surfdrive/LogitMDA/gmr4/lorenza/paper1-gmr3/new/boot.gmr3.R")
source("~/surfdrive/LogitMDA/gmr4/lorenza/paper1-gmr3/new/plot.boot.gmr3.R")
source("~/surfdrive/LogitMDA/gmr4/lorenza/paper1-gmr3/new/auxiliary.gmr3.R")
source("~/surfdrive/LogitMDA/gmr4/lorenza/paper1-gmr3/new/fit.gmr3.R")

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

source("~/surfdrive/LogitMDA/gmr4/lorenza/paper1-gmr3/new/predict.gmr3.R")
source("~/surfdrive/LogitMDA/gmr4/lorenza/paper1-gmr3/new/xval.gmr3.R")
source("~/surfdrive/LogitMDA/gmr4/lorenza/paper1-gmr3/new/plot.xval.gmr3.R")

PE = xval.gmr3(Yb = YYb, Yo = YYo, X = X, Xscale = Xscale, S = 1:5, repeats = 10, mseed = 12)
PE2 = plot.xval.gmr3(PE)

# ggsave("~/surfdrive/LogitMDA/gmr4/lorenza/paper1-gmr3/xvalplot.pdf", plot = PE2$plot, width = 11.7, height = 8.3, units = "in", limitsize = FALSE)

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

# ggsave("~/surfdrive/LogitMDA/gmr4/lorenza/paper1-gmr3/bootstrap.pdf", plot = plts$BV, width = 11.7, height = 8.3, units = "in", limitsize = FALSE)

source("~/surfdrive/LogitMDA/gmr4/lorenza/paper1-gmr3/new/plot.quantifications.R")

mylabs2 = c("Left", "Center", "Right")
mylabs3 = c("Male", "Female")
mylabs4 = c("Rural Area/Village", "Small/Middle Town", "Large Town")
mylabs5 = c("Pre-primary", "Primary" , "Low Secondary", "Up Secondary", 
            "Post Secondary", "Tertiary", "Bachelor", "Master", "Doctorate")

plt2 = plot.quantifications(out2, var = 2, labels = mylabs2)
plt3 = plot.quantifications(out2, var = 3, labels = mylabs3)
plt4 = plot.quantifications(out2, var = 4, labels = mylabs4)
plt5 = plot.quantifications(out2, var = 5, labels = mylabs5)

cplot = ggarrange(plt2, plt3, plt4, plt5, nrow = 2, ncol = 2)
# ggsave("~/surfdrive/LogitMDA/gmr4/lorenza/paper1-gmr3/os.pdf", plot = cplot, width = 11.7, height = 11.7, units = "in", limitsize = FALSE)
