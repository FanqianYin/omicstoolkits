# Statistical parameters of data

# second column of pdata are group information
calculate_RSD <- function(exp, pdata){
  groups <- unique(pdata$class)
  RSD.df <- data.frame("variable" = colnames(exp))
  for (g in groups) {
    exp.group <- exp[pdata$class==g,]
    RSD.group <- apply(exp.group, 2, function(x) sd(x)/mean(x)*100)
    RSD.df <- cbind(RSD.df, RSD.group)
  }
  colnames(RSD.df) <- c("variable", groups)
  return(RSD.df)
}
