




#' cluster_number
#'
#' Produce one referenced number of clusters which determined by kmean/pam/clara/fanny/hcut clustering.
#'
#' @param Exp data.frame with rownames represent samples, colnames represent features like proteins or genes.
#' @param k.max maximum numbers of cluster.
#' @param FUNcluster Clustering method
#' @param scale Logical. whether input Exp should be scaled.
#'
#' @return Two plots. Indicating number of clusters.
#' @export
#'
#' @examples
#' cluster_number(iris[-5], 5, kmeans)
#'



cluster_number <- function(Exp, k.max, FUNcluster = c(kmeans, pam, clara, fanny, hcut), scale = TRUE){
  d <- as.matrix(Exp)
  if (scale == TRUE) d <- scale(d)
  require(factoextra)
  require(cluster)
  p1 <- fviz_nbclust(d, FUNcluster, method = "wss")+labs(subtitle = "Elbow method")
  p2 <- fviz_nbclust(d, FUNcluster, method = "silhouette")+labs(subtitle = "Silhouette method")
  require(cowplot)
  plot_grid(p1,p2)
}
