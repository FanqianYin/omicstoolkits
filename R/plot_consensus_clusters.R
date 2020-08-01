#' @name plot_consensus_clusters
#' @aliases plot_consensus_clusters.heatmap
#' @aliases plot_consensus_clusters.pca
#' @title plot consensus result
#' @param consensus_clusters Result of Consensus_Cluster_Analysis().
#'
#' @rdname plot_consensus_clusters.heatmap
#' @export
#'
#' @example
#' plot_consensus_clusters.heatmap(consensus_clusters)
#'
plot_consensus_clusters.heatmap <- function(consensus_clusters, clustering_method_columns = "ward.D2", clustering_method_rows = "complete"){
  #consenesus cluster heatmap
  if (consensus_clusters$k <= 9) { #k <=9, use RColorBrewer::brewer.pal
    Heatmap(consensus_clusters[["cluster_matrix"]],
            col = structure(RColorBrewer::brewer.pal(n = consensus_clusters$k, name = "Set1"), names = as.character(1:consensus_clusters$k)),
            cluster_rows = T, cluster_columns = T,
            clustering_method_columns = clustering_method_columns,
            clustering_method_rows = clustering_method_rows, name="clusters", row_title = "samples", column_title = "clutering algorithms",
            row_split = consensus_clusters$k)+
      Heatmap(as.character(consensus_clusters$consensus_cluster), col = structure(RColorBrewer::brewer.pal(n = consensus_clusters$k, name = "Set1"), names = as.character(1:consensus_clusters$k)))
  }else {#k > 9, use random colors()
    Heatmap(consensus_clusters[["cluster_matrix"]],
            col = 1:consensus_clusters$k,
            cluster_rows = T, cluster_columns = T,
            clustering_distance_rows = "pearson", clustering_distance_columns = "pearson", clustering_method_columns = "ward.D2",
            clustering_method_rows = "average", name="clusters", row_title = "samples", column_title = "clutering algorithms",
            row_split = consensus_clusters$k)+
      rowAnnotation(df = data.frame("consensus_cluster" = as.character(consensus_clusters$consensus_cluster)))
  }

}
#rowAnnotation(df = data.frame("consensus_cluster" = as.character(consensus_clusters$consensus_cluster)), col = list(consensus_cluster = structure(RColorBrewer::brewer.pal(n = consensus_clusters$k, name = "Set1"), names = as.character(1:consensus_clusters$k))))
#rowAnnotation(df = data.frame("consensus_cluster" = as.character(consensus_clusters$consensus_cluster)), col = list(consensus_cluster = structure(RColorBrewer::brewer.pal(n = consensus_clusters$k, name = "Set1"), names = as.character(1:consensus_clusters$k))))
#' @rdname plot_consensus_clusters.pca
#' @export
#'
#' @example
#' plot_consensus_clusters.pca(consensus_clusters)
plot_consensus_clusters.pca <- function(consensus_clusters){
  #pca plot
  fviz_cluster(list(data = consensus_clusters$d, cluster = consensus_clusters$consensus_cluster),geom = "point",
               palette = c("#2E9FDF", "#00AFBB", "#E7B800", "#FC4E07"),
               ellipse.type = "convex", # Concentration ellipse
               repel = TRUE, # Avoid label overplotting (slow)
               show.clust.cent = FALSE, ggtheme = theme_minimal(),main = "Consensus cluster")
}

#' @rdname plot_cluster_algorithms.pca
#' @export
#'
#' @example
#' plot_cluster_algorithms.pca(consensus_clusters)
plot_cluster_algorithms.pca <- function(consensus_clusters){
  #pca plot
  fviz_cluster(list(data = consensus_clusters$cluster_matrix, cluster = consensus_clusters$consensus_cluster),geom = "point",
               palette = c("#2E9FDF", "#00AFBB", "#E7B800", "#FC4E07"),
               ellipse.type = "convex", # Concentration ellipse
               repel = TRUE, # Avoid label overplotting (slow)
               show.clust.cent = FALSE, ggtheme = theme_minimal(),main = "consistency of clustering algorithms")
}

#to be developed functions:
#plot_all_clusters.cpa

#' @export
#plot_cluster.pca
#plot pca of single clustering algorithm result
plot_cluster.pca <- function(consensus_clusters, method){
  if (method == "ALL") {

  }
  fviz_cluster(list(data = consensus_clusters$d, cluster = as.data.frame(consensus_clusters$cluster_matrix)[[method]]),geom = "point",
               palette = c("#2E9FDF", "#00AFBB", "#E7B800", "#FC4E07"),
               ellipse.type = "convex", # Concentration ellipse
               repel = TRUE, # Avoid label overplotting (slow)
               show.clust.cent = FALSE, ggtheme = theme_minimal(),main = method)
}

#plot_exp_heatmap

