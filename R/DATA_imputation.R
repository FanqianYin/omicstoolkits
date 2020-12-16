# data imputation

#method: min0.5 (half of the minimum positive value), min0.2, knn, constant(0,), bpca, svd,
#' Title
#'
#' @param exp
#' @param pdata
#' @param method min0.5 (half of the minimum positive value), min0.2, knn, constant(0,), bpca, svd,
#' @param constant
#' @param na
#'
#' @return
#' @export
#'
#' @examples
fill_NA <- function(exp, pdata, method = "knn", constant = 100, na=NA){
  exp <- as.matrix(exp)
  if(!is.na(na)) exp[exp==na] <- NA
  exp.imptuted <- exp
  if (method == "knn"){
    p <- ncol(exp)
    exp.imptuted <- impute::impute.knn(t(exp.imptuted), maxp = p)
    exp.imptuted <- t(exp.imptuted[["data"]])
  }else if (method == "min0.2"){
    exp.imptuted[is.na(exp.imptuted)] <- min(exp.imptuted, na.rm = T)*0.2
  }else if (method == "min0.5"){
    exp.imptuted[is.na(exp.imptuted)] <- min(exp.imptuted, na.rm = T)*0.5
  }else if (method == "constant"){
    exp.imptuted <- exp
  }else if (method == "bpca"){
    stop("to be developed")
  }else if (method == "svd"){
    stop("to be developed")
  }
  return(exp.imptuted)
}
