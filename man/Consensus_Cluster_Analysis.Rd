% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CA_Consensus_Cluster_Analysis.R
\name{Consensus_Cluster_Analysis}
\alias{Consensus_Cluster_Analysis}
\alias{clustering.methods}
\alias{Consensus_Cluster_conMethod}
\title{Consensus_Cluster_Analysis}
\usage{
Consensus_Cluster_Analysis(
  Exp,
  k,
  methods = "ALL",
  dist.method = "euclidean",
  consensus.method = "average",
  scale = TRUE,
  saveplot = FALSE,
  plot.dir.suffix = NULL,
  ConsensusClusterPlus.reps = 100
)

clustering.methods()

Consensus_Cluster_conMethod(res.con, consensus.method)
}
\arguments{
\item{Exp}{data.frame or matrix with rownames represent samples, colnames represent features like proteins or genes.}

\item{k}{number of clusters}

\item{methods}{Charactor vector. Default "ALL" will perform all 20 clustering methods.}

\item{dist.method}{Method for calculating dissimilarity. Default "euclidean", can be any methods supported by dist()}

\item{consensus.method}{Method used for the final consensus cluster by hieraichical clustering.
Default are set to "average". Alternatively, "complete" may get a diferent view of result, which you should try. It Can be set as any method supported by hclust.}

\item{scale}{Logical value, whether to scale your input data. Default: TRUE.}

\item{saveplot}{under developing}

\item{plot.dir.suffix}{under developing}

\item{ConsensusClusterPlus.reps}{Integer value, used by consensusClusterPlus.}

\item{res.con}{result of \code{Consensus_Cluster_Analysis()}.}
}
\value{
\code{Consensus_Cluster_Analysis} returns res.con--one list contains all cluster results.

\code{clustering.methods} returns all supported clutersing algorithms.

\code{Consensus_Cluster_conMethod} returns one new res.con with the new consensus.method.
}
\description{
Consensus clustering result of multi-clustering algorithms.
}
\details{
Consensus clustering result of multi-clustering algorithms.
}
\examples{
#Use ALL supported clustering methods
res.con <- Consensus_Cluster_Analysis(iris[-5], 3, methods = "ALL")
#Use a subset of clustering methods. Note: not avaliable now.
#Under development: res.con <- Consensus_Cluster_Analysis(iris[-5], 3, methods = c("kmeans", "pam"))

#plot result:
plot_consensus_cluster.heatmap(res.con)
plot_consensus_cluster.pca(res.con)

#list all supported clustering algorithms:
clustering.methods()
all.methods <- unlist(clustering.methods)

#Change consensus cluster method:
res.con <- Consensus_Cluster_conMethod(res.con, "complete")
}
