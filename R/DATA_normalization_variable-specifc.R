# variable-specific normalization (QC-based signal correction)

# data as data.frame: row as sample, colum as metabolite(or features/variables)
# method: qc.mean, qc.median, loess (QCRLSC), random forest (QCRFSC), QC-SVR
# details for methods:
#                    BRDG(qc.median): for each metabolite, divide its value by the median level in bridge samples.
#                    MED: for each metabolite, divide its value by the median across the experimental samples (classs).
#                    QCRFSC or QCRLSC from statTarget.
# design: including 4 columns: "sample", "batch", "class", "order" (run order). QC sample as "QC" in column "class".
#' variable-specific normalization (QC-based signal correction)
#'
#' @param expData data.frame: row as sample, colum as metabolite(or features/variables)
#' @param design data.frame, including 4 columns: "sample", "batch", "class", "order" (run order). QC sample as "QC" in column "class".
#' @param method QC-based signal correction methods, including: qc.mean, qc.median, loess (QCRLSC), random forest (QCRFSC), QC-SVR, IORLSC, IORFSC
#' @param rf.ntree ntree par for randomForest::randomForest()
#' @param loess.span span par for loess()
#' @param loess.degree degree par for loess()
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
#'
normalize_metabolite_specific <- function(expData, design, method,
                                          rf.ntree = 500, loess.span = 0.5, loess.degree = 2,
                                          ...){
  expData <- as.matrix(expData)
  expData[expData==0] <- NA
  #missing value check
  if (sum(is.na(expData))>0){
    missing.number <- sum(is.na(expData))
    cat("This data contians ", missing.number, " missing values.\n")
    minHalf <- min(expData, na.rm = T)/2
    expData[is.na(expData)] <- minHalf
    cat("All missing values are automatically replaced by: ", minHalf,", half minimum positive value of this dataset.\n")
    cat("Recommandation: use knn to impute missing value before data normalization.\n")
  }

  expData <- t(expData)
  expData.qc <- expData[, design$class == "QC"]
  qc.index <- grep("QC",design$class, fixed = TRUE)
  runorder <- design$order

  if (method == "qc.mean") {
    nf <- rowMeans(expData.qc)
    expData.normlized <- expData/nf
  }else if (method == "qc.median") {
    nf <- base::apply(expData.qc, 1, median)
    expData.normlized <- expData/nf
  }else if (method == "QCRFSC"){
    nf.matirx <- matrix(nrow = nrow(expData), ncol = ncol(expData))
    pb <- txtProgressBar(min = 1, max = nrow(expData), style = 3)
    cat("Using QCRFSC for data normalization:\n")
    for (v in 1:nrow(expData.qc)) {
      model.rf <- randomForest::randomForest(data.frame(order=qc.index), expData.qc[v,], ntree=rf.ntree)
      nf <- predict(model.rf, data.frame(order=runorder))
      nf.matirx[v,] <- nf
      setTxtProgressBar(pb, v)
    }
    close(pb)
    expData.normlized <- expData/nf.matirx
  }else if (method == "IORFSC"){
    nf.matirx <- matrix(nrow = nrow(expData), ncol = ncol(expData))
    pb <- txtProgressBar(min = 1, max = nrow(expData), style = 3)
    cat("Using IORFSC for data normalization:\n")
    for (v in 1:nrow(expData)) {
      model.rf <- randomForest::randomForest(data.frame(order=runorder), expData[v,], ntree=rf.ntree)
      nf <- predict(model.rf, data.frame(order=runorder))
      nf.matirx[v,] <- nf
      setTxtProgressBar(pb, v)
    }
    close(pb)
    expData.normlized <- expData/nf.matirx
  }else if(method == "QCRLSC"){
    nf.matirx <- matrix(nrow = nrow(expData), ncol = ncol(expData))
    pb <- txtProgressBar(min = 1, max = nrow(expData), style = 3)
    cat("Using QCRLSC for data normalization:\n")
    for (v in 1:nrow(expData.qc)) {
      model.loess <- stats::loess(expData.qc[v,]~qc.index, span = loess.span, degree = loess.degree)
      nf <- predict(model.loess, runorder)
      nf.matirx[v,] <- nf
      setTxtProgressBar(pb, v)
    }
    close(pb)
    expData.normlized <- expData/nf.matirx
  }else if(method == "IORLSC"){
    nf.matirx <- matrix(nrow = nrow(expData), ncol = ncol(expData))
    pb <- txtProgressBar(min = 1, max = nrow(expData), style = 3)
    cat("Using IORLSC for data normalization:\n")
    for (v in 1:nrow(expData)) {
      model.loess <- stats::loess(expData[v,]~runorder, span = loess.span, degree = loess.degree)
      nf <- predict(model.loess, runorder)
      nf.matirx[v,] <- nf
      setTxtProgressBar(pb, v)
    }
    close(pb)
    expData.normlized <- expData/nf.matirx
  }
  return(t(expData.normlized))
}


#' @export
#QCRFSC, QCRLSC from statTarget
# expData as data.frame: row as sample, colum as metabolite(or features/variables); must contains row.names as feature names
#method.QC are one of: "QCRFSC" or "QCRLSC"
#design including 4 columns: sample, batch, class, order (run order)
normalization_QCRFSC_statTarget <- function(expData, design, method.QC, impute.method = "KNN", ...){
  #require(statTarget)
  #require(tidyverse)
  print("Using statTarget for QC-based normalization")
  #create temporary working dir
  dir.temp <- "./tempDir_statTarget"
  #dir.tempInput <- "./tempDir_statTarget/inputexpData"
  dir.create(dir.temp)
  dir.create("./tempDir_statTarget/inputexpData")
  #write formated input expData to ./tempDir_statTarget/inputexpData
  #design$class[] <- stringr::str_replace_all(design$class, "QC", "NA")
  design$class[design$class=="QC"] <- NA
  write.csv(design, file = "./tempDir_statTarget/inputexpData/pexpData.csv", row.names = FALSE, quote = FALSE)
  #exp expData
  expData.t <- expData
  expData.t <- t(expData.t)
  expData.t <- as.data.frame(expData.t)
  expData.t <- cbind("name"=rownames(expData.t),expData.t)
  write.csv(expData.t, file = "./tempDir_statTarget/inputexpData/exp.csv", row.names = FALSE, quote = FALSE)
  #statTarget::shiftCor
  samPeno <- "./inputexpData/pexpData.csv"
  samFile <- "./inputexpData/exp.csv"
  setwd(dir.temp)
  statTarget::shiftCor(samPeno,samFile, MLmethod = method.QC, imputeM = impute.method, ...)
  #read result
  #"statTarget/shiftCor/After_shiftCor/shift_all_cor.csv"
  result.expData <- read.csv("./statTarget/shiftCor/After_shiftCor/shift_all_cor.csv")
  rownames(result.expData) <- result.expData$sample
  result.expData <- result.expData[-1:-2] #remove sample and class columns.
  #delete temp dir
  setwd("..")
  unlink(dir.temp, recursive = TRUE)

  return(result.expData)
}



#to be dev
#1.find a way to optimize span par in loess.
#1.add qc-svr
#2.dev new method:qc_bin, for every metabolite, take qc value as normlaize factor within the same bin? (bin methods: bin with same interval, or find neighberhood? wihin on feature.)
#3.dev new method:qc_free_rfsc/rlsc (based on the randomed sample class injection order.)

#3.dev new method:qc_free_rfsc/rlsc

#
