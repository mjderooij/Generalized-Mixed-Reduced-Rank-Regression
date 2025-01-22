plot.quantifications = function(object, var, labels = NULL){
  library(ggplot2)
  
  x = object$originalcategories[[var]]
  nc = length(x)
  if(is.null(labels)){labels = paste0("c", 1:nc)}
  # x <- factor(x) # , levels = x, labels = mlabels, ordered = T)

  y = object$quantifications[[var]] 
  df = as.data.frame(cbind(x,y, labels))
  
  plt = ggplot(df, aes(x = x, y = y, group = 1)) + 
    geom_point() + geom_step() + 
    scale_x_discrete(labels = labels) + 
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 30, vjust = 0.5)) + 
    xlab("Original Categories") +
    ylab("Quantifications") + 
    labs(title = paste("Quantification:", object$xnames[var], sep = " "))
  

  print(plt)
  return(plt)
}
