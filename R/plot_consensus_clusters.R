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
plot_consensus_clusters.heatmap <- function(consensus_clusters){
  #consenesus cluster heatmap
  Heatmap(consensus_clusters[["cluster_matrix"]], cluster_rows = T, cluster_columns = T, clustering_method_columns = "ward.D2",
          clustering_method_rows = "average", name="clusters", row_title = "samples", column_title = "clutering algorithms",
          row_split = consensus_clusters$k)+
    rowAnnotation(df = data.frame("consensus_cluster" = as.character(consensus_clusters$consensus_cluster)))
}

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
