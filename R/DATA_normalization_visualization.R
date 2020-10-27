# data normalization-visualization

#exp, row as sample, column as varibles
#plot: boxplot (all, QC), boxplot of RSD, PCA, traj: trajectory of single metabolite (all and QC: with or without high RSD.)
normalize_data_visualization <- function(exp, pdata, plot, order = "fix", type = "default", variable.name = "variable",
                                         expression.unit = "expression", traj.single.metabo.RSD.quantile = 0.1, log = T){
  if (log==T) {
    df <- cbind(pdata, log2(exp))
  }else if (log==F) {
    df <- cbind(pdata, exp)
  }
  df.longer <- pivot_longer(df, -1:-ncol(pdata), names_to = variable.name, values_to = expression.unit)
  if (order == "fix") df.longer$sample <- factor(df.longer$sample, levels = pdata$sample)

  if (plot == "boxplot") {
    p1 <- ggplot(df.longer, aes(x=sample, y=expression, fill = time))+geom_boxplot()+
      theme(axis.text.x = element_text(angle = 60, vjust = 0.5, hjust=0.3))
    p2 <- ggplot(filter(df.longer,time=="QC"), aes(x=sample, y=expression, fill = time))+geom_boxplot()+
      theme(axis.text.x = element_text(angle = 60, vjust = 0.5, hjust=0.3))
    if (type == "onlyQC") {p <- p2} else {p <- p1}

  }else if(plot == "RSD"){
    exp.rsd <- calculate_RSD(exp, pdata)
    exp.rsd.longer <- pivot_longer(exp.rsd, -1, names_to = "group", values_to = "RSD")
    p <- ggplot(exp.rsd.longer, aes(group, RSD, fill=group))+geom_boxplot()
  }else if(plot == "PCA"){
    pca <- PCA(log2(exp+1), graph = F)
    p <- fviz_pca_ind(pca, geom = "point", habillage = as.factor(pdata$time))
  }else if(plot == "traj"){
    #exp <- exp.neg
    #pdata <- pdata.neg
    #traj.single.metabo.RSD.quantile = 0.1
    #rm(exp, pdata, exp.rsd, exp.df, metabolite.index, plot.df, traj.single.metabo.RSD.quantile)

    exp.rsd <- calculate_RSD(exp, pdata) %>% arrange(QC)
    exp.df <- as.data.frame(exp)
    metabolite.index <- exp.rsd[["variable"]][round(nrow(exp.rsd)*traj.single.metabo.RSD.quantile,0)]

    plot.df <- data.frame(pdata, metabolite = exp.df[metabolite.index]) %>% mutate("color.group" = str_detect(time,"QC"))
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
