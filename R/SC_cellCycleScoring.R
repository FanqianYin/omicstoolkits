#cell cycle scoring and normalzation

#methods for cell cycle classification rpkgs: peco, scran, Seurat


#' Score cell cycle phases
#' Cell cycle classification including G2/M, S, G1.
#'
#' @param counts raw counts matrix: row as genes, column as cells/samples. rownames must be gene SYMBOL.
#' @param method one of: "Seurat",
#' @param cc.genes a list, with cc.genes and s.genes vectors. Defult are cc.genes from Seurat package.
#' @param normalization.method Method for counts normalization: "LogNormalize", "CLR", "RC". see \code{\link[Seurat]{NormalizeData}}
#' @param graph.pca whether to show pca graph which only based on 97 cell cycle genes or not.
#'
#' @return factor vector for classfication (G1, G2/M, S)
#' @export
#' @seealso \link{https://satijalab.org/seurat/v3.2/cell_cycle_vignette.html}
#' @examples
SC_cellCycleScoring <- function(counts, method = "Seurat",
                                graph.pca = TRUE, cc.genes = NULL, normalization.method = "LogNormalize"){
  if(is.null(cc.genes)){
    #cc.genes <- c("MCM5","PCNA","TYMS","FEN1","MCM2","MCM4","RRM1","UNG","GINS2","MCM6","CDCA7","DTL","PRIM1","UHRF1","MLF1IP","HELLS","RFC2","RPA2","NASP","RAD51AP1","GMNN","WDR76","SLBP","CCNE2","UBR7","POLD3","MSH2","ATAD2","RAD51","RRM2","CDC45","CDC6","EXO1","TIPIN","DSCC1","BLM","CASP8AP2","USP1","CLSPN","POLA1","CHAF1B","BRIP1","E2F8","HMGB2","CDK1","NUSAP1","UBE2C","BIRC5","TPX2","TOP2A","NDC80","CKS2","NUF2","CKS1B","MKI67","TMPO","CENPF","TACC3","FAM64A","SMC4","CCNB2","CKAP2L","CKAP2","AURKB","BUB1","KIF11","ANP32E","TUBB4B","GTSE1","KIF20B","HJURP","CDCA3","HN1","CDC20","TTK","CDC25C","KIF2C","RANGAP1","NCAPD2","DLGAP5","CDCA2","CDCA8","ECT2","KIF23","HMMR","AURKA","PSRC1","ANLN","LBR","CKAP5","CENPE","CTCF","NEK2","G2E3","GAS2L3","CBX5","CENPA")
    cc.genes <- Seurat::cc.genes
    s.genes <- cc.genes$s.genes
    g2m.genes <- cc.genes$g2m.genes
  } else {
    s.genes <- cc.genes$s.genes
    g2m.genes <- cc.genes$g2m.genes
  }

  if (method == "Seurat") {
    require(Seurat)
    obj <- CreateSeuratObject(counts = counts)
    obj <- NormalizeData(obj, normalization.method = normalization.method)
    obj <- FindVariableFeatures(obj, selection.method = "vst")
    obj <- ScaleData(obj, features = rownames(obj))
    obj <- CellCycleScoring(obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

    if(graph.pca==TRUE){
      obj <- RunPCA(obj, features = c(s.genes, g2m.genes))
      DimPlot(obj)
    }
    cell_cycle_class <- obj@active.ident
    return(cell_cycle_class)

  }else if (method == "peco"){

  }

}
