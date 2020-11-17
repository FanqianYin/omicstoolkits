# variable-specific normalization (QC-based signal correction, )

# data as data.frame: row as sample, colum as metabolite(or features/variables)
# method: qc.mean, qc.median, loess (QCRLSC), random forest (QCRFSC), QC-SVR
# details for methods:
#                    BRDG(qc.median): for each metabolite, divide its value by the median level in bridge samples.
#                    MED: for each metabolite, divide its value by the median across the experimental samples (groups).
#                    QCRFSC or QCRLSC from statTarget.
# design: including 4 columns: "sample", "batch", "class", "order" (run order). QC sample as "QC" in column "class".
normalize_metabolite_specific <- function(data, design, method, keep.unit = FALSE){
  qc <- data[design$group == "QC"]

  if (keep.unit == FALSE) {
    if (method == "qc.mean") {
      mean.qc <- rowMeans(qc)
      data.norm <- "to_be_complished"
    } else if (method == "qc.median") {

    }

  }

  if (keep.unit == TRUE){
    if (method == "QCRFSC"){


    }else if(method == "QCRLSC"){

    }
  }


}



#QCRFSC, QCRLSC from statTarget
# data as data.frame: row as sample, colum as metabolite(or features/variables); must contains row.names as feature names
#method.QC are one of: "QCRFSC" or "QCRLSC"
#design including 4 columns: sample, batch, class, order (run order)
normalization_QCRFSC_statTarget <- function(data, design, method.QC, impute.method = "KNN"){
  #require(statTarget)
  #require(tidyverse)
  print("Using statTarget for QC-based normalization")
  #create temporary working dir
  dir.temp <- "./tempDir_statTarget"
  #dir.tempInput <- "./tempDir_statTarget/inputdata"
  dir.create(dir.temp)
  dir.create("./tempDir_statTarget/inputdata")
  #write formated input data to ./tempDir_statTarget/inputdata
  #design$class[] <- stringr::str_replace_all(design$class, "QC", "NA")
  design$class[design$class=="QC"] <- NA
  write.csv(design, file = "./tempDir_statTarget/inputdata/pdata.csv", row.names = FALSE, quote = FALSE)
  #exp data
  data.t <- data
  data.t <- t(data.t)
  data.t <- as.data.frame(data.t)
  data.t <- cbind("name"=rownames(data.t),data.t)
  write.csv(data.t, file = "./tempDir_statTarget/inputdata/exp.csv", row.names = FALSE, quote = FALSE)
  #statTarget::shiftCor
  samPeno <- "./inputdata/pdata.csv"
  samFile <- "./inputdata/exp.csv"
  setwd(dir.temp)
  statTarget::shiftCor(samPeno,samFile, MLmethod = method.QC, imputeM = impute.method, coCV = 30)
  #read result
  #"statTarget/shiftCor/After_shiftCor/shift_all_cor.csv"
  result.data <- read.csv("./statTarget/shiftCor/After_shiftCor/shift_all_cor.csv")
  rownames(result.data) <- result.data$sample
  result.data <- result.data[-1:-2] #remove sample and class columns.
  #delete temp dir
  setwd("..")
  unlink(dir.temp, recursive = TRUE)

  return(result.data)
}


