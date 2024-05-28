plot.xval.gmr3 = function(object){
  # auxiliary function
  data_summary <- function(data, varname, groupnames){
    require(plyr)
    summary_func <- function(x, col){
      c(mean = mean(x[[col]], na.rm=TRUE),
        sd = sd(x[[col]], na.rm=TRUE))
    }
    data_sum<-ddply(data, groupnames, .fun=summary_func, varname)
    data_sum <- rename(data_sum, c("mean" = varname))
    return(data_sum)
  }
  
  K = max(object$fold)
  nrep = max(object$repeats)
  
  PEsummary = data_summary(object, "PEtotal", "S")
  PEsummary$se = PEsummary$sd/sqrt(nrep * K) # repeats * K
  PEsummary
  
  xvalplt = ggplot(PEsummary, aes(x=S, y=PEtotal)) + 
    geom_line(col = "red") +
    geom_point(size = 3, col = "red") +
    geom_errorbar(aes(ymin=PEtotal-se, ymax=PEtotal+se), col = "red", width= 0, position=position_dodge(0.05)) + 
    labs(x = "Rank/Dimensionality", y = "Prediction Error") + theme_bw()
  
  print(xvalplt)
  
  output = list(
    plot = xvalplt,
    summary = PEsummary
  )
  return(output)
}