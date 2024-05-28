plot.boot.gmr3 = function(object){
  library(ggplot2)
  library(ggpubr)
  
  B = as.data.frame(object$gmr3obj$B)
  colnames(B) = paste0("dim", 1:ncol(B))
  rownames(B) = object$gmr3obj$xnames
  P = nrow(B)
  B$Predictor = c(1:P)
  B$Predictor = factor(B$Predictor, levels = 1:P, labels = rownames(B))
  
  V = as.data.frame(object$gmr3obj$V)
  colnames(V) = paste0("dim", 1:ncol(V))
  rownames(V) = object$gmr3obj$ynames
  R = nrow(V)
  V$Response = c(1:R)
  V$Response = factor(V$Response, levels = 1:R, labels = rownames(V))
  
  Bplot = ggplot(object$BBdf, aes(dim1, dim2, color = Predictor), show.legend = F) +
    geom_point(alpha = 0.05, show.legend = F) + 
    stat_ellipse(type = "norm", show.legend = F) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_vline(xintercept = 0, linetype = 2) +
    coord_fixed() + theme_bw() +
    labs(title = "Weights plot", x = "Dimension 1", y = "Dimension 2")
  
  Bplot = Bplot + 
    geom_text(data = B, aes(x = dim1, y = dim2, label = Predictor), color = "black", show.legend = F)
  
  Vplot = ggplot(object$BVdf, aes(dim1, dim2, color = Response), show.legend = F) + 
    geom_point(alpha = 0.05, show.legend = F) + 
    stat_ellipse(type = "norm", show.legend = F) + 
    geom_hline(yintercept = 0, linetype = 2) + 
    geom_vline(xintercept = 0, linetype = 2) +
    coord_fixed() + theme_bw() + 
    labs(title = "Loadings plot", x = "Dimension 1", y = "Dimension 2")
  
  Vplot = Vplot + 
    geom_text(data = V, aes(x = dim1, y = dim2, label = Response), color = "black", show.legend = F)
  
  # combine the two plots
  BVplot = ggarrange(Bplot, Vplot, ncol = 2)
  print(BVplot)
  
  output = list(
    BVplot = BVplot, 
    Bplt = Bplot,
    Vplt = Vplot
  )
  return(output)
}
