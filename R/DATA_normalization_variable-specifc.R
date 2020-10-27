# variable-specific normalization

# data as data.frame
# method: qc.mean, qc.median, loess,
# details for methods:
#                    BRDG(qc.median): for each metabolite, divide its value by the median level in bridge samples.
#                    MED: for each metabolite, divide its value by the median across the experimental samples (groups).
# design: with group colum, QC sample as "QC" in which.
normalize_metabolite_specific <- function(data, design, method, keep.unit = FALSE){
  qc <- data[design$group == "QC"]

  if (keep.unit == FALSE) {
    if (method == "qc.mean") {
      mean.qc <- rowMeans(qc)
      data.norm <-
    } else if (method == "qc.median") {

    }

  }


}
