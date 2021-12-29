#regression analysis


#' Regression analysis
#'
#' @param pdata
#' @param exp
#' @param pr.variable
#' @param variable.cofounder character vector, indicating cofounder like age, region, gender etc. in pdata
#' @param resp.variable character vector, indicating 1 or more resp.var (if > 1 resp.var, they will be calculated seperately)
#'
#' @return
#' @export
#'
#' @examples
#'
DATA_lmregression <- function(pdata, exp, pr.variable, variable.cofounder = NULL, resp.variable){
  df <- cbind(pdata, exp)
  if(!all(c(pr.variable, variable.cofounder) %in% colnames(df))) stop("Predict variable defined by pr.varible should exist in pdata or exp")
  if(!all(resp.variable %in% colnames(df))) stop("Responsable variable defined by pr.varible should exist in pdata or exp")
  lm.res <- list()
  fmla.pr.variable <- paste0(variable.cofounder,collapse = "+")

  if (length(resp.variable)==1) {
    for (v in 1:length(pr.variable)) {
	  if(!is.null(variable.cofounder)){
        fmla <- formula(paste0(resp.variable, "~", pr.variable[v], "+", fmla.pr.variable))
	  }else {
	  fmla <- formula(paste0(resp.variable, "~", pr.variable[v]))
	  }
      fit.summary <- summary(lm(fmla, df))

      cofounder.name <- rownames(fit.summary$coefficients)[-1]
      lm.res.colname <- c("varibale",paste(rep(cofounder.name,each=2),c("coef","p"),sep = "_"))
      lm.res[pr.variable[v]] <- c(pr.variable[v])
	  if(nrow(fit.summary$coefficients)>=2) {
      for (vn in 1:(nrow(fit.summary$coefficients)-1)) {
        lm.res[[pr.variable[v]]] <- c(lm.res[[pr.variable[v]]],fit.summary$coefficients[vn+1,1],fit.summary$coefficients[vn+1,4])
       }
	  }
    }
  }else if (length(resp.variable)>1){
    for (v in 1:length(resp.variable)) {
      if(!is.null(variable.cofounder)){
        fmla <- formula(paste0(resp.variable[v], "~", pr.variable, "+", fmla.pr.variable))
      }else {
        fmla <- formula(paste0(resp.variable[v], "~", pr.variable))
      }
      fit.summary <- summary(lm(fmla, df))

      cofounder.name <- rownames(fit.summary$coefficients)[-1]
      lm.res.colname <- c("varibale",paste(rep(cofounder.name,each=2),c("coef","p"),sep = "_"))
      lm.res[resp.variable[v]] <- c(resp.variable[v])
      if(nrow(fit.summary$coefficients)>=2) {
        for (vn in 1:(nrow(fit.summary$coefficients)-1)) {
          lm.res[[resp.variable[v]]] <- c(lm.res[[resp.variable[v]]],fit.summary$coefficients[vn+1,1],fit.summary$coefficients[vn+1,4])
        }
      }
    }
  }
  #unlist the res to df
  res.df <- lm.res[[1]]
  for (v in 2:length(lm.res)) {
    res.df <- rbind(res.df, lm.res[[v]])
  }
  colnames(res.df) <- lm.res.colname
  colnames(res.df)[2:3] <- c("coeffecient","p.value")
  res.df <- as.data.frame(res.df)
  for (v in colnames(res.df)[-1]) {
    res.df[[v]] <- as.numeric(res.df[[v]])
  }

  return(res.df)

}
