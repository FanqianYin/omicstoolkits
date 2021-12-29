# Statistical parameters of data

# second column of pdata are group information
#' Title
#'
#' @param exp data.frame: row as sample, colum as metabolite(or features/variables)
#' @param pdata data.frame: row-samples, column-phenotype variables or factors.
#'
#' @return data.frame includes RSD of every variable within pre-defined groups.
#' @export
#'
#' @examples
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
