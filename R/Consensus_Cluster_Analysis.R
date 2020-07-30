

#' Consensus_Cluster_Analysis
#'
#' Consensus clustering result of multi-clustering algorithms.
#'
#' @param Exp data.frame or matrix with rownames represent samples, colnames represent features like proteins or genes.
#' @param saveplot under developing
#' @param plot.dir.suffix under developing
#' @param k number of clusters
#' @param methods Charactor vector. Default "ALL" will perform all 20 clustering methods.
#' @param dist.method Method for calculating dissimilarity. Default "euclidean", can be any methods supported by dist()
#' @param methods.consensus Method used for the final consensus cluster by hieraichical clustering. Default are set to "average".
#' @param scale Logical value, whether to scale your input data. Default: TRUE.
#' @param plot.cluster.pca Logical value, whether to plot pca of every clustering results. Default: FALSE.
#' @param ConsensusClusterPlus.reps Integer value, used by consensusClusterPlus.
#'
#' @return List contain all cluster results. Also return consensus cluster heatmap and PCA plot indicating the consensus cluster result.
#' @export
#'
#' @examples
#' consensus_clusters <- Consensus_Cluster_Analysis(iris[-5], 3, methods = "ALL")
#' plot_consensus_clusters(consensus_clusters)
#'
#' consensus_clusters <- Consensus_Cluster_Analysis(iris[-5], 3, methods = c("kmeans", "pam"))
#'
Consensus_Cluster_Analysis <- function(Exp, k, methods = "ALL", dist.method = "euclidean", methods.consensus = "average", scale = TRUE,
                                       plot.cluster.pca = FALSE, saveplot = FALSE, plot.dir.suffix = NULL,
                                       ConsensusClusterPlus.reps = 100){
  #Prepare input data
  d <- as.matrix(Exp)
  d <- scale(d)
  res.dist <- dist(d, method = dist.method)
  methods.all <- unlist(clustering.methods)
  if(methods == "ALL") methods <- methods.all #If TRUE, all methods are used
  #Hierarchical Clustering
  if(any(methods == "hc_ward.D")) res_hc_ward.D <- hclust(d = res.dist, method = "ward.D"); res_hc_ward.D$cluster <- cutree(res_hc_ward.D, k)
  if(any(methods == "hc_ward.D2")) res_hc_ward.D2 <- hclust(d = res.dist, method = "ward.D2"); res_hc_ward.D2$cluster <- cutree(res_hc_ward.D2, k)
  if(any(methods == "hc_complete")) res_hc_complete <- hclust(d = res.dist, method = "complete"); res_hc_complete$cluster <- cutree(res_hc_complete, k)
  if(any(methods == "hc_single")) res_hc_single <- hclust(d = res.dist, method = "single"); res_hc_single$cluster <- cutree(res_hc_single, k)
  if(any(methods == "hc_average")) res_hc_average <- hclust(d = res.dist, method = "average"); res_hc_average$cluster <- cutree(res_hc_average, k)
  if(any(methods == "hc_diana")) res_hc_diana <- diana(x = d, stand = F, metric = "euclidean" ); res_hc_diana$cluster <- cutree(res_hc_diana, k)
  if(any(methods == "hc_mcquitty")) res_hc_mcquitty <- hclust(d = res.dist, method = "mcquitty"); res_hc_mcquitty$cluster <- cutree(res_hc_mcquitty, k)
  if(any(methods == "hc_median")) res_hc_median <- hclust(d = res.dist, method = "median"); res_hc_median$cluster <- cutree(res_hc_median, k)
  if(any(methods == "hc_centroid")) res_hc_centroid <- hclust(d = res.dist, method = "centroid"); res_hc_centroid$cluster <- cutree(res_hc_centroid, k)

  #kmean
  if(any(methods == "kmeans")) res_kmeans <- kmeans(d, k, nstart = 10, iter.max = 20)
  if(any(methods == "pam")) res_pam <- pam(d, k)
  if(any(methods == "fanny")) res_fanny <- fanny(d, k)
  if(any(methods == "hkmeans")) res_hkmeans <- hkmeans(d, k)

  #Other methods
  if(any(methods == "Mclust")) res_mc <- Mclust(d, G=1:k)

  #ConsensusClusterPlus
  if(any(methods == "con_kmeans")) res_con_kmeans <- ConsensusClusterPlus(t(d), k+2, reps = ConsensusClusterPlus.reps ,clusterAlg = "km")
  if(any(methods == "con_pam")) res_con_pam <- ConsensusClusterPlus(t(d), k+2, reps = ConsensusClusterPlus.reps ,clusterAlg = "pam")
  if(any(methods == "con_hc_ward.D2")) res_con_hc_ward.D2 <- ConsensusClusterPlus(t(d), k+2, reps = ConsensusClusterPlus.reps ,clusterAlg = "hc", innerLinkage = "ward.D2")
  if(any(methods == "con_hc_complete")) res_con_hc_complete <- ConsensusClusterPlus(t(d), k+2, reps = ConsensusClusterPlus.reps ,clusterAlg = "hc", innerLinkage = "complete")
  if(any(methods == "con_hc_average")) res_con_hc_average <- ConsensusClusterPlus(t(d), k+2, reps = ConsensusClusterPlus.reps ,clusterAlg = "hc", innerLinkage = "average")


  #All result
  consensus_clusters <- list("consensus_cluster" = NULL,
                             "Hierarchical Clustering" = list("res_hc_ward.D" = res_hc_ward.D, "res_hc_ward.D2" = res_hc_ward.D2, "res_hc_complete" = res_hc_complete,
                                                              "res_hc_single" = res_hc_single, "res_hc_average" = res_hc_average,
                                                              "res_hc_mcquitty" = res_hc_mcquitty, "res_hc_median" = res_hc_median,
                                                              "res_hc_centroid" = res_hc_centroid, "res_hc_diana" = res_hc_diana),
                             "Partitioning clustering" = list("res_kmeans" = res_kmeans, "res_pam" = res_pam, "res_fanny" = res_fanny, "res_hkmeans" = res_hkmeans),
                             "Other clustering" = list("Mclust" = res_mc),
                             "ConsensusClusterPlus" = list("con_kmeans" = res_con_kmeans, "con_pam" = res_con_pam, "con_hc_ward.D2" = res_con_hc_ward.D2,
                                                           "con_hc_complete" = res_con_hc_complete, "con_hc_average" = res_con_hc_average))
  #cluster matrix
  cluster_matrix <- data.frame("hc_ward.D" = consensus_clusters[["Hierarchical Clustering"]][["res_hc_ward.D"]][["cluster"]],
                               "hc_ward.D2" = consensus_clusters[["Hierarchical Clustering"]][["res_hc_ward.D2"]][["cluster"]],
                               "hc_complete" = consensus_clusters[["Hierarchical Clustering"]][["res_hc_complete"]][["cluster"]],
                               "hc_average" = consensus_clusters[["Hierarchical Clustering"]][["res_hc_average"]][["cluster"]],
                               "hc_single" = consensus_clusters[["Hierarchical Clustering"]][["res_hc_single"]][["cluster"]],
                               "hc_mcquitty" = consensus_clusters[["Hierarchical Clustering"]][["res_hc_mcquitty"]][["cluster"]],
                               "hc_median" = consensus_clusters[["Hierarchical Clustering"]][["res_hc_median"]][["cluster"]],
                               "hc_centroid" = consensus_clusters[["Hierarchical Clustering"]][["res_hc_centroid"]][["cluster"]],
                               "hc_diana" = consensus_clusters[["Hierarchical Clustering"]][["res_hc_diana"]][["cluster"]],
                               "res_kmeans" = consensus_clusters[["Partitioning clustering"]][["res_kmeans"]][["cluster"]],
                               "res_pam" = consensus_clusters[["Partitioning clustering"]][["res_pam"]][["clustering"]],
                               "res_fanny" = consensus_clusters[["Partitioning clustering"]][["res_fanny"]][["clustering"]],
                               "res_hkmeans" = consensus_clusters[["Partitioning clustering"]][["res_hkmeans"]][["cluster"]],
                               "Mclust" = consensus_clusters[["Other clustering"]][["Mclust"]][["classification"]],
                               "con_kmeans" = consensus_clusters[["ConsensusClusterPlus"]][["con_kmeans"]][[k]][["consensusClass"]],
                               "con_pam" = consensus_clusters[["ConsensusClusterPlus"]][["con_pam"]][[k]][["consensusClass"]],
                               "con_hc_ward.D2" = consensus_clusters[["ConsensusClusterPlus"]][["con_hc_ward.D2"]][[k]][["consensusClass"]],
                               "con_hc_complete" = consensus_clusters[["ConsensusClusterPlus"]][["con_hc_complete"]][[k]][["consensusClass"]],
                               "con_hc_average" = consensus_clusters[["ConsensusClusterPlus"]][["con_hc_average"]][[k]][["consensusClass"]]
                               ) %>% as.matrix()
  #Consensus cluster result:
  consensus_cluster <- dist(cluster_matrix) %>% hclust(method = methods.consensus) %>% cutree(k)
  consensus_clusters$consensus_cluster <- consensus_cluster
  consensus_clusters$cluster_matrix <- cluster_matrix
  consensus_clusters$d <- d
  consensus_clusters$k <- k
  return(consensus_clusters)
}

# all supported cluster methods
clustering.methods <- list("Hierarchical Clustering" = c("hc_ward.D", "hc_ward.D2", "hc_single", "hc_complete", "hc_average", "hc_mcquitty", "hc_median", "hc_centroid", "hc_diana"),
                           "Partitioning clustering" = c("kmeans", "pam", "fanny", "hkmeans"),
                           "Model-based clustering" = c("Mclust"),
                           "Density-based clustering" = c("dbscan"),
                           "ConsensusClusterPlus" = c("con_kmeans", "con_pam", "con_hc_ward.D2", "con_hc_complete", "con_hc_average"))

