#' cluster_number
#'
#' Produce one referenced number of clusters which determined by kmean/pam/clara/fanny/hcut clustering.
#'
#' @param Exp data.frame or matrix with rownames represent samples, colnames represent features like proteins or genes.
#' @param k.max maximum numbers of cluster.
#' @param FUNcluster Clustering method. Type without quote mark. Chose one of following methods: kmean, pam, clara, fanny, hcut.
#' @param scale Logical. whether input Exp should be scaled.
#' @param ... Any parameter passed to fviz_nbclust().
#'
#' @return Two plots. Indicating number of clusters by using Elbow method and Silhouette method.
#'
#' @examples
#' cluster_number(iris[-5], 5, kmeans)
#'
#' @seealso \code{\link[factoextra]{fviz_nbclust}}
#'
#' @export

cluster_number <- function(Exp, k.max, FUNcluster = c(kmeans, pam, clara, fanny, hcut), scale = TRUE, ...){
  d <- as.matrix(Exp)
  if (scale == TRUE) d <- scale(d)
  p1 <- fviz_nbclust(d, FUNcluster, method = "wss", ...)+labs(subtitle = "Elbow method")
  p2 <- fviz_nbclust(d, FUNcluster, method = "silhouette", ...)+labs(subtitle = "Silhouette method")
  cowplot::plot_grid(p1,p2)
}
