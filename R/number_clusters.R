

#number of clusters. Determined by kmean clustering.
number_clusters <- function(Exp, k.max, FUNcluster = c(kmeans, pam, clara, fanny, hcut), scale = TRUE){
  d <- as.matrix(Exp)
  if (scale == TRUE) d <- scale(d)
  require(factoextra)
  require(cluster)
  p1 <- fviz_nbclust(d, FUNcluster, method = "wss")+labs(subtitle = "Elbow method")
  p2 <- fviz_nbclust(d, FUNcluster, method = "silhouette")+labs(subtitle = "Silhouette method")
  require(cowplot)
  plot_grid(p1,p2)
}
