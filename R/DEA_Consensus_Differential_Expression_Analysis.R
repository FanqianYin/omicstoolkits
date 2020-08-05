
#Differential expression analysis(DEA)

#Consensus Differential Expression analysis
#' Consensus_Differential_Expression
#'
#' @param Exp data.frame or matrix with rownames represent samples, colnames represent features like proteins or genes.
#' @param groups Charactor or factor vector indicates group information.
#' @param adj.method
#'
#' @return One list including all results of Differential Expression analysis
#'
#'
#' @examples
#'
Consensus_Differential_Expression <- function(Exp, groups, adj.method = "BH"){

}

#t.test
DEA_t.test <- function(Exp, groups, adj.method = "BH"){

}

#one-way analysis of variance (ANOVA)
DEA_ANOVA <- function(){

}


# To be dev: nonparametric alternatives:
# Mann–Whitney U.test for unpaired Student’s t‑test; Wilcoxon signed‑rank testforpaired Student’s t‑test ,Kruskal–Wallis test as the nonparametric equivalent
# of ANOVA; the Friedman’s test as the counterpart of repeated measures ANOVA.
