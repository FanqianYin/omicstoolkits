

#' Consensus_Cluster_Analysis
#'
#' Consensus clustering result of multi-clustering algorithms.
#'
#' Consensus clustering result of multi-clustering algorithms.
#'
#' @param Exp data.frame or matrix with rownames represent samples, colnames represent features like proteins or genes.
#' @param saveplot under developing
#' @param plot.dir.suffix under developing
#' @param k number of clusters
#' @param methods Charactor vector. Default "ALL" will perform all 20 clustering methods.
#' @param dist.method Method for calculating dissimilarity. Default "euclidean", can be any methods supported by dist()
#' @param consensus.method Method used for the final consensus cluster by hieraichical clustering.
#' Default are set to "average". Alternatively, "complete" may get a diferent view of result, which you should try. It Can be set as any method supported by hclust.
#' @param scale Logical value, whether to scale your input data. Default: TRUE.
#' @param ConsensusClusterPlus.reps Integer value, used by consensusClusterPlus.
#' @param res.con result of \code{Consensus_Cluster_Analysis()}.
#' @return
#' \code{Consensus_Cluster_Analysis} returns res.con--one list contains all cluster results.
#'
#' \code{clustering.methods} returns all supported clutersing algorithms.
#'
#' \code{Consensus_Cluster_conMethod} returns one new res.con with the new consensus.method.
#' @examples
#' #Use ALL supported clustering methods
#' res.con <- Consensus_Cluster_Analysis(iris[-5], 3, methods = "ALL")
#' #Use a subset of clustering methods. Note: not avaliable now.
#' #Under development: res.con <- Consensus_Cluster_Analysis(iris[-5], 3, methods = c("kmeans", "pam"))
#'
#' #plot result:
#' plot_consensus_cluster.heatmap(res.con)
#' plot_consensus_cluster.pca(res.con)
#'
#' #list all supported clustering algorithms:
#' clustering.methods()
#' all.methods <- unlist(clustering.methods)
#'
#' #Change consensus cluster method:
#' res.con <- Consensus_Cluster_conMethod(res.con, "complete")
#' @name Consensus_Cluster_Analysis
NULL

#' @rdname Consensus_Cluster_Analysis
#' @export
Consensus_Cluster_Analysis <- function(Exp, k, methods = "ALL", dist.method = "euclidean", consensus.method = "average", scale = TRUE,
                                       saveplot = FALSE, plot.dir.suffix = NULL,
                                       ConsensusClusterPlus.reps = 100){
  #Prepare input data
  d <- as.matrix(Exp)
  d <- scale(d)
  res.dist <- dist(d, method = dist.method)
  methods.all <- unlist(clustering.methods())
  if(methods == "ALL") methods <- methods.all #If TRUE, all methods are used
  #Hierarchical Clustering
  if(any(methods == "hc_ward.D")) {res_hc_ward.D <- hclust(d = res.dist, method = "ward.D"); res_hc_ward.D$cluster <- cutree(res_hc_ward.D, k)}
  if(any(methods == "hc_ward.D2")) {res_hc_ward.D2 <- hclust(d = res.dist, method = "ward.D2"); res_hc_ward.D2$cluster <- cutree(res_hc_ward.D2, k)}
  if(any(methods == "hc_complete")) {res_hc_complete <- hclust(d = res.dist, method = "complete"); res_hc_complete$cluster <- cutree(res_hc_complete, k)}
  if(any(methods == "hc_single")) {res_hc_single <- hclust(d = res.dist, method = "single"); res_hc_single$cluster <- cutree(res_hc_single, k)}
  if(any(methods == "hc_average")) {res_hc_average <- hclust(d = res.dist, method = "average"); res_hc_average$cluster <- cutree(res_hc_average, k)}
  if(any(methods == "hc_diana")) {res_hc_diana <- diana(x = d, stand = F, metric = "euclidean" ); res_hc_diana$cluster <- cutree(res_hc_diana, k)}
  if(any(methods == "hc_mcquitty")) {res_hc_mcquitty <- hclust(d = res.dist, method = "mcquitty"); res_hc_mcquitty$cluster <- cutree(res_hc_mcquitty, k)}
  if(any(methods == "hc_median")) {res_hc_median <- hclust(d = res.dist, method = "median"); res_hc_median$cluster <- cutree(res_hc_median, k)}
  if(any(methods == "hc_centroid")) {res_hc_centroid <- hclust(d = res.dist, method = "centroid"); res_hc_centroid$cluster <- cutree(res_hc_centroid, k)}

  #kmean
  if(any(methods == "kmeans")) res_kmeans <- kmeans(d, k, nstart = 10, iter.max = 20)
  if(any(methods == "pam")) res_pam <- pam(d, k)
  if(any(methods == "fanny")) res_fanny <- fanny(d, k)
  if(any(methods == "hkmeans")) res_hkmeans <- hkmeans(d, k)

  #Other methods
  if(any(methods == "Mclust")) res_mc <- Mclust(d, G=1:k)

  #ConsensusClusterPlus
  if(any(methods == "con_kmeans")) res_con_kmeans <- ConsensusClusterPlus(t(d), k+1, reps = ConsensusClusterPlus.reps ,clusterAlg = "km")
  if(any(methods == "con_pam")) res_con_pam <- ConsensusClusterPlus(t(d), k+1, reps = ConsensusClusterPlus.reps ,clusterAlg = "pam")
  if(any(methods == "con_hc_ward.D2")) res_con_hc_ward.D2 <- ConsensusClusterPlus(t(d), k+1, reps = ConsensusClusterPlus.reps ,clusterAlg = "hc", innerLinkage = "ward.D2")
  if(any(methods == "con_hc_complete")) res_con_hc_complete <- ConsensusClusterPlus(t(d), k+1, reps = ConsensusClusterPlus.reps ,clusterAlg = "hc", innerLinkage = "complete")
  if(any(methods == "con_hc_average")) res_con_hc_average <- ConsensusClusterPlus(t(d), k+1, reps = ConsensusClusterPlus.reps ,clusterAlg = "hc", innerLinkage = "average")


  #All result
  res.con <- list("consensus_cluster" = NULL,
                             "par" = list("k" = k,
                                          "methods" = ifelse(methods == "ALL", unlist(clustering.methods), methods),
                                          "dist.method" = dist.method,
                                          "consensus.method" = consensus.method,
                                          "scale" = scale,
                                          "saveplot" = saveplot,
                                          "plot.dir.suffix" = plot.dir.suffix,
                                          "ConsensusClusterPlus.reps" = ConsensusClusterPlus.reps
                                          ),
                             "Hierarchical Clustering" = list("hc_ward.D" = res_hc_ward.D, "hc_ward.D2" = res_hc_ward.D2, "hc_complete" = res_hc_complete,
                                                              "hc_single" = res_hc_single, "hc_average" = res_hc_average,
                                                              "hc_mcquitty" = res_hc_mcquitty, "hc_median" = res_hc_median,
                                                              "hc_centroid" = res_hc_centroid, "hc_diana" = res_hc_diana),
                             "Partitioning clustering" = list("kmeans" = res_kmeans, "pam" = res_pam, "fanny" = res_fanny, "hkmeans" = res_hkmeans),
                             "Other clustering" = list("Mclust" = res_mc),
                             "ConsensusClusterPlus" = list("con_kmeans" = res_con_kmeans, "con_pam" = res_con_pam, "con_hc_ward.D2" = res_con_hc_ward.D2,
                                                           "con_hc_complete" = res_con_hc_complete, "con_hc_average" = res_con_hc_average))
  #cluster matrix
  cluster_matrix <- data.frame("hc_ward.D" = res.con[["Hierarchical Clustering"]][["hc_ward.D"]][["cluster"]],
                               "hc_ward.D2" = res.con[["Hierarchical Clustering"]][["hc_ward.D2"]][["cluster"]],
                               "hc_complete" = res.con[["Hierarchical Clustering"]][["hc_complete"]][["cluster"]],
                               "hc_average" = res.con[["Hierarchical Clustering"]][["hc_average"]][["cluster"]],
                               "hc_single" = res.con[["Hierarchical Clustering"]][["hc_single"]][["cluster"]],
                               "hc_mcquitty" = res.con[["Hierarchical Clustering"]][["hc_mcquitty"]][["cluster"]],
                               "hc_median" = res.con[["Hierarchical Clustering"]][["hc_median"]][["cluster"]],
                               "hc_centroid" = res.con[["Hierarchical Clustering"]][["hc_centroid"]][["cluster"]],
                               "hc_diana" = res.con[["Hierarchical Clustering"]][["hc_diana"]][["cluster"]],
                               "kmeans" = res.con[["Partitioning clustering"]][["kmeans"]][["cluster"]],
                               "pam" = res.con[["Partitioning clustering"]][["pam"]][["clustering"]],
                               "fanny" = res.con[["Partitioning clustering"]][["fanny"]][["clustering"]],
                               "hkmeans" = res.con[["Partitioning clustering"]][["hkmeans"]][["cluster"]],
                               "Mclust" = res.con[["Other clustering"]][["Mclust"]][["classification"]],
                               "con_kmeans" = res.con[["ConsensusClusterPlus"]][["con_kmeans"]][[k]][["consensusClass"]],
                               "con_pam" = res.con[["ConsensusClusterPlus"]][["con_pam"]][[k]][["consensusClass"]],
                               "con_hc_ward.D2" = res.con[["ConsensusClusterPlus"]][["con_hc_ward.D2"]][[k]][["consensusClass"]],
                               "con_hc_complete" = res.con[["ConsensusClusterPlus"]][["con_hc_complete"]][[k]][["consensusClass"]],
                               "con_hc_average" = res.con[["ConsensusClusterPlus"]][["con_hc_average"]][[k]][["consensusClass"]]
                               ) %>% as.matrix()
  #Consensus cluster result:
  res.con$consensus_result <- dist(cluster_matrix) %>% hclust(method = consensus.method)
  consensus_cluster <- res.con$consensus_result %>% cutree(k)
  res.con$consensus_cluster <- consensus_cluster
  res.con$cluster_matrix <- cluster_matrix
  res.con$d <- d
  res.con$k <- k
  return(res.con)
}

#' @rdname Consensus_Cluster_Analysis
#' @export
clustering.methods <- function(){
  list("Hierarchical Clustering" = c("hc_ward.D", "hc_ward.D2", "hc_single", "hc_complete", "hc_average", "hc_mcquitty", "hc_median", "hc_centroid", "hc_diana"),
                           "Partitioning clustering" = c("kmeans", "pam", "fanny", "hkmeans"),
                           "Model-based clustering" = c("Mclust"),
                           "Density-based clustering" = c("dbscan"),
                           "ConsensusClusterPlus" = c("con_kmeans", "con_pam", "con_hc_ward.D2", "con_hc_complete", "con_hc_average"))
}

#' @rdname Consensus_Cluster_Analysis
#' @export
# Change consensus.method and return new result.
Consensus_Cluster_conMethod <- function(res.con, consensus.method){
  res.con[["par"]][["consensus.method"]] <- consensus.method
  #Consensus cluster result:
  res.con$consensus_result <- dist(res.con$cluster_matrix) %>% hclust(method = consensus.method)
  consensus_cluster <- res.con$consensus_result %>% cutree(res.con$k)
  res.con$consensus_cluster <- consensus_cluster
  return(res.con)
}

#Future dev
#cluster funtion
#fastcluster::hclust, NMF,
