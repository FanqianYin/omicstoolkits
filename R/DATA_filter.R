# data filteration


#method: ALL: RSD+missing.value; RSD: all metabolites of RSD > 25% within QC samples are remvoed. missing.value: all metabolites of missing.value ratio > 0.3 within all samples are remvoed.
#na = NA ,or na = 0, or else.
filter_data <- function(exp, pdata, method = "ALL", RSD = 25, missing.value = 0.3, na = NA){
  exp.rsd <- calculate_RSD(exp, pdata)
  exp.missing.ratio <- apply(exp, 2, function(x) sum(x==na)/nrow(exp))
  if (method == "ALL"){
    filter.index <- (exp.rsd["QC"] < RSD) & (exp.missing.ratio < missing.value)
    exp <- as.data.frame(exp)
    exp.filtered <- exp[,filter.index]
    print(paste0(sum(filter.index), " variables are kept, ", sum(!filter.index), " variables are removed."))
  }else if(method == "RSD"){
    filter.index <- (exp.rsd["QC"] < RSD)
    exp <- as.data.frame(exp)
    exp.filtered <- exp[,filter.index]
    print(paste0(sum(filter.index), " variables are kept, ", sum(!filter.index), " variables are removed."))
  }else if(method == "missing.value"){
    filter.index <- (exp.missing.ratio < missing.value)
    exp <- as.data.frame(exp)
    exp.filtered <- exp[,filter.index]
    print(paste0(sum(filter.index), " variables are kept, ", sum(!filter.index), " variables are removed."))
  }
  return(exp.filtered)
}



