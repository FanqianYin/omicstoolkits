# data filteration


#method: ALL: RSD+missing.value; RSD: all metabolites of RSD > 25% within QC samples are remvoed. missing.value: all metabolites of missing.value ratio > 0.3 within all samples are remvoed.
#na = NA ,or na = 0, or else.
#exp, row: samples
#padat, need one "class" colum
filter_data <- function(exp, pdata, method = "ALL", RSD = 25, missing.value = 0.3, na = NA){
  #exp.rsd <- calculate_RSD(exp, pdata)
  #exp.missing.ratio <- apply(exp, 2, function(x) sum(x==na)/nrow(exp))
  exp <- as.matrix(exp)
  if(!is.na(na)) exp[exp==na] <- NA
  if (method == "ALL"){
    exp.rsd <- calculate_RSD(exp, pdata)
    exp.missing.ratio <- apply(exp, 2, function(x) sum(is.na(x))/nrow(exp))
    filter.index <- (exp.rsd["QC"] < RSD) & (exp.missing.ratio < missing.value)
    exp <- as.data.frame(exp)
    exp.filtered <- exp[,filter.index]
    print(paste0(sum(filter.index), " variables are kept, ", sum(!filter.index), " variables are removed."))
  }else if(method == "RSD"){
    exp.rsd <- calculate_RSD(exp, pdata)
    filter.index <- (exp.rsd["QC"] < RSD)
    exp <- as.data.frame(exp)
    exp.filtered <- exp[,filter.index]
    print(paste0(sum(filter.index), " variables are kept, ", sum(!filter.index), " variables are removed."))
  }else if(method == "missing.value"){
    exp.missing.ratio <- apply(exp, 2, function(x) sum(is.na(x))/nrow(exp))
    filter.index <- (exp.missing.ratio < missing.value)
    exp <- as.data.frame(exp)
    exp.filtered <- exp[,filter.index]
    print(paste0(sum(filter.index), " variables are kept, ", sum(!filter.index), " variables are removed."))
  }
  return(exp.filtered)
}



