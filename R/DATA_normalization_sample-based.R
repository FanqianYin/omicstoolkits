#sample-based normalization

# normalize_sample_based methods hypothesize that: most metabolite don't change, so that the distributions of sample data are similar.
# data: sample as columan, metabolites as rows
# chose method as one of: tic, mean, median, range, mad, MSTUS, IS (internal standard).
# design: sample information including case/control, QC samples, run battch, time point, or other classffication.
#' Sample-based normalization
#' Normalise data by sample-based scaling factor, like sum or mean of sample.
#' @param expData data.frame: row as sample, colum as metabolite(or features/variables)
#' @param design data.frame, including 4 columns: "sample", "batch", "class", "order" (run order). QC sample as "QC" in column "class".
#' @param method including: "tic", "median", "mean", "range", "mad", "MSTUS".
#' @param na only take 0 as missing value, not NA.
#'
#' @return
#' @export
#'
#' @examples
normalize_sample_based <- function(expData, design, method, na=0){
  data <- t(expData)
    if (method == "tic") {
      data.norm <- apply(data, 2, function(x) x/sum(x))
    }else if (method == "mean") {
      data.norm <- apply(data, 2, function(x) x/mean(x))
    }else if (method == "median") {
      data.norm <- apply(data, 2, function(x) x/median(x))
    }else if (method == "range") {
      data.norm <- apply(data, 2, function(x) x/range(x))
    }else if (method == "mad") {
      data.norm <- apply(data, 2, function(x) x/mad(x))
    }else if (method == "MSTUS") {
      data <- as.matrix(data)
      index.tus <- apply(1, data, function(x) sum(x != na)==length(x)) #take value=0 as missing value
      data.tus <- data[index.tus,] #metabolites with value in all sample
      normfactor.tus <- colSums(data.tus)
      data <- rbind(data,normfactor.tus)
      data.norm <- apply(data, 2, function(x) x/tail(x))
      data.norm <- data.norm[1:(nrow(data.norm)-1),]
    }

  return(data.norm)
}
