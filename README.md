# omicstoolkits
## Integrated bioinformatics toolkits for quantitative omics data analysis

This is an integrated toolkits focusing on data analysis of omics data like transcriptomic, proteomic, metabolic and any other quantitative omics dataset.

Under development...

## Main toolkits

### Consensus_Cluster_Analysis() 
An toolkit for consensus clustering with multi-clustering algorithms.

Currently supported clustering methods: 

hc_ward.D, hc_ward.D2, hc_single, hc_complete, hc_average, hc_mcquitty, hc_median, hc_centroid, hc_diana, kmeans, pam, fanny, hkmeans, Mclust, dbscan, con_kmeans, con_pam, con_hc_ward.D2, con_hc_complete, con_hc_average.
## Under development

1. Consensus differential expression analysis

2. Network analysis

3. Regression analysis

4. Time serie anaylsis

5. TBD

Details for future development plans can be viewed at https://github.com/FanqianYin/omicstoolkits/Features_under_developing.html.

## Install package

```{r}
install.packages("devtools")
devtools::install_github("FanqianYin/omicstoolkits")
```
## Contact

If you encounter any problem, don't hesitate to cantact me: fanqianyin@gmail.com or yinfanqian@mail.kiz.ac.cn.

Or report bugs at: https://github.com/FanqianYin/omicstoolkits/issues.

Suggestion: If you have   any adivce on these tools, please let me know;

Collaboration: If you are interested in integrating these omicstoolkits, welcome to cooperate.

## Development history

7/30/2020  Consensus_Cluster_Analysis: an toolkit focus on sample-based subtyping by using consensus clustering result of multi-clustering algorithms.

## Version history

**Current version: v0.1.1**

v0.1.1  8/2/2020

Consensus_Cluster_Analysis toolkit: Add plot_consensus_clusters.heatmap(), plot_consensus_clusters.pca(), plot_cluster_algorithms.pca(), plot_cluster.pca(), for visulaztion result of Consensus_Cluster_Analysis().
    
v0.1.0  7/30/2020

Consensus_Cluster_Analysis toolkit: First release.
