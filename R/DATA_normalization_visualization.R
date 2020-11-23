# data normalization-visualization

#exp, row as sample, column as varibles
#plot: boxplot (all, QC), boxplot of RSD, PCA, traj: trajectory of single metabolite (all and QC: with or without high RSD.)
#' Data normalization-visualization
#' Using PCA, boxplot, trajectory of single variable to visulize data.
#' @param exp data.frame: row as sample, colum as metabolite(or features/variables)
#' @param pdata data.frame: row as sample, colum as phenotype data(or factors)
#' @param plot type of plot, including: boxplot, RSD, PCA, traj
#' @param order default=TRUE, fix the order of sample in the boxplot.
#' @param type "default": display all samples or "onlyQC": only display QC samples in boxplot.
#' @param variable.name
#' @param expression.unit
#' @param traj.single.metabo.RSD.quantile RSD quantile for drawing the RSD boxplot.
#' @param log whether log2 transform the exp data.
#'
#' @return
#' @export
#'
#' @examples
normalize_data_visualization <- function(exp, pdata, plot, order.fixed = TRUE, type = "default", variable.name = "variable",
                                         expression.unit = "expression", traj.single.metabo.RSD.quantile = 0.1, log = TRUE){
  if (log==TRUE) {
    df <- cbind(pdata, log2(exp))
  }else if (log==FALSE) {
    df <- cbind(pdata, exp)
  }
  df.longer <- pivot_longer(df, -1:-ncol(pdata), names_to = variable.name, values_to = expression.unit)
  if (order.fixed == TRUE) df.longer$sample <- factor(df.longer$sample, levels = pdata$sample)

  if (plot == "boxplot") {
    p1 <- ggplot(df.longer, aes(x=sample, y=expression, fill = class))+geom_boxplot()+
      theme(axis.text.x = element_text(angle = 60, vjust = 0.5, hjust=0.3))
    p2 <- ggplot(filter(df.longer,class=="QC"), aes(x=sample, y=expression, fill = class))+geom_boxplot()+
      theme(axis.text.x = element_text(angle = 60, vjust = 0.5, hjust=0.3))
    if (type == "onlyQC") {p <- p2} else {p <- p1}

  }else if(plot == "RSD"){
    exp.rsd <- calculate_RSD(exp, pdata)
    exp.rsd.longer <- pivot_longer(exp.rsd, -1, names_to = "group", values_to = "RSD")
    p <- ggplot(exp.rsd.longer, aes(group, RSD, fill=group))+geom_boxplot(outlier.shape = NA)+ylim(c(1,100))
  }else if(plot == "PCA"){
    pca <- FactoMineR::PCA(log2(exp+1), graph = F)
    p <- fviz_pca_ind(pca, geom = "point", habillage = as.factor(pdata$class))
  }else if(plot == "traj"){
    #exp <- exp.neg
    #pdata <- pdata.neg
    #traj.single.metabo.RSD.quantile = 0.1
    #rm(exp, pdata, exp.rsd, exp.df, metabolite.index, plot.df, traj.single.metabo.RSD.quantile)

    exp.rsd <- calculate_RSD(exp, pdata) %>% arrange(QC)
    exp.df <- as.data.frame(exp)
    metabolite.index <- exp.rsd[["variable"]][round(nrow(exp.rsd)*traj.single.metabo.RSD.quantile,0)]

    plot.df <- data.frame(pdata, metabolite = exp.df[metabolite.index]) %>% mutate("color.group" = str_detect(class,"QC"))
    colnames(plot.df)[ncol(pdata)+1] <- "metabolite"
    plot.df$color.group[plot.df$color.group == T] <- "QC"
    plot.df$color.group[plot.df$color.group == F] <- "sample"
    #fix the runorder for plot
    plot.df <- arrange(plot.df, runorder)
    #plot.df$runorder <- factor(plot.df$runorder, levels = plot.df$runorder)

    p <- ggplot(plot.df, aes(x=runorder, y=metabolite, color = color.group))+
      geom_point()+
      scale_color_manual(values=c("red", "gray"))+scale_size_manual(values=c(3,1))+
      theme_classic()+ggtitle("Expression signal drift and run order")+geom_smooth(se = F)

  }

  return(p)
}

#example
#

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
