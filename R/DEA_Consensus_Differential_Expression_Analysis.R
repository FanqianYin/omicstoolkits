
#Differential expression analysis(DEA)

#Consensus Differential Expression analysis
#' Consensus_Differential_Expression
#' @aliases DEA_ANOVA
#' @param Exp data.frame or matrix with rownames represent samples, colnames represent features like proteins or genes.
#' @param groups Charactor or factor vector indicates group information. In two-way anova, groups are a data.frame, each column indicate one group factor.
#' @param adj.method
#' @param type "one-way" or "two-way" anova.
#' @return
#' \code{Consensus_Differential_Expression}
#' One list including all results of Differential Expression analysis
#' \code{DEA_ANOVA}
#' Anova result with pairwised comparation.
#'
#' @examples
#' res.anova <- DEA_ANOVA(iris[-5], groups = iris$Species)
#' res.ttest <- DEA_t.test(exp, groups)
Consensus_Differential_Expression <- function(Exp, groups, methods = "ALL", adj.method = "BH", fearture.name = "gene", center.fun = "mean"){

}

#' @export
#t.test
DEA_t.test <- function(Exp, groups, adj.method = "BH", fearture.name = "gene", center.fun = "mean"){
  Exp <- as.matrix(Exp)
  res.df <- data.frame("variale" = colnames(Exp))

  p.value <- apply(Exp, 2, function(x) {
    #fo <- as.formula("x", " ~ ", "group")
    t.test.result <- t.test(x[groups==levels(groups)[1]],x[groups==levels(groups)[2]])
    p <- t.test.result$p.value
    return(p)
  })
  p.value <- as.vector(p.value)
  res.df <- cbind(res.df,p.value)
  colnames(res.df) <- c(fearture.name, "p")
  #p.adjust
  res.df <- res.df %>% mutate("p.adj" = p.adjust(p, method = adj.method))
  #add other descriptive statistics


  res.fc <- calculate_FC(Exp, groups, group.number = 2, fun = center.fun, fearture.name = fearture.name)
  res.df <- dplyr::left_join(res.df, res.fc) %>% dplyr::arrange(p.adj)
  return(res.df)

}

#' @export
#one-way analysis of variance (ANOVA), parametric, >2 groups.
DEA_ANOVA <- function(Exp, groups, type = "one-way", group.number = 3, adj.method = "BH", fearture.name = "gene", center.fun = "mean"){
  Exp <- as.data.frame(Exp)
  res.df <- data.frame()
  if (type == "one-way") {
    cn <- colnames(Exp)
    for (n in 1:ncol(Exp)) {
      fo <- as.formula(paste0(cn[n], " ~ ", "groups"))
      res.aov <- aov(fo, Exp)
      p <- summary(res.aov)[[1]][["Pr(>F)"]][1]

      res.tukeyHSD <- TukeyHSD(res.aov)
      if(n == 1) {rnames <- rownames(res.tukeyHSD[[1]]) #Only one run for creating comparation names list
        rnames <- paste0("p.adj", "_", rnames)
        }
      p.multi.compare <- res.tukeyHSD[[1]][,4] #pairwise comparation p value of groups. length are decided by the comparation.

      res.df <- rbind(res.df, c(p, p.multi.compare))
    }
    res.df <- data.frame(fearture.name = cn, res.df)
    colnames(res.df) <- c(fearture.name, "p", rnames) #finish calcluating multi-p values
    res.df <- res.df %>% mutate("p.adj" = p.adjust(p, method = adj.method)) %>% select(1,2,ncol(.), everything())
    #add other descriptive statistics
    res.fc <- calculate_FC(Exp, groups, fun = center.fun, fearture.name = fearture.name)
    res.df <- dplyr::left_join(res.df, res.fc)
    return(res.df)
  }
}

#limma imbeded test method
DEA_limma.test <- function(Exp, groups, adj.method = "BH", fearture.name = "gene", center.fun = "mean"){

}

#Kruskal-Wallis test, non-parametric, >2 groups.



# To be dev: nonparametric alternatives:
# Mann–Whitney U.test for unpaired Student’s t‑test; Wilcoxon signed‑rank testforpaired Student’s t‑test ,Kruskal–Wallis test as the nonparametric equivalent
# of ANOVA; the Friedman’s test as the counterpart of repeated measures ANOVA.

#todo before perform any test, the normalty and homogeneity of variance should be checked. If not pass, no test will be performed until it pass.


#supported functions
#calculate descriptive statistcs of multi-groups; return df with n*m, n is number of features, and m is groups*funs
#fun = c("mean", "median"); groups, charactor or logical vector
calculate_mean <- function(data, groups, fun = "mean", fearture.name = "gene"){
  data <- as.data.frame(data)
  res.means <- data.frame(fearture.name = colnames(data))
  #fun = "mean"
  if (fun == "mean") {
    res.descri <- data.frame("groups" = groups, data) %>% group_by(groups) %>% summarise_all(mean)
    name.groups <- as.character(res.descri$groups);name.groups <- paste0(fun, ".", name.groups)
    res.descri <- t(res.descri[-1]) #matrix, row as varibles(genes), column as result of group means

    res.means <- cbind(res.means, res.descri)
    colnames(res.means) <- c(fearture.name, name.groups)
    return(res.means)
  }
  if (fun == "median") {
    res.descri <- data.frame("groups" = groups, data) %>% group_by(groups) %>% summarise_all(median)
    name.groups <- as.character(res.descri$groups);name.groups <- paste0(fun, ".", name.groups)
    res.descri <- t(res.descri[-1]) #matrix, row as varibles(genes), column as result of group means

    res.medians <- cbind(res.medians, res.descri)
    colnames(res.medians) <- c(fearture.name, name.groups)
    return(res.medians)
  }
}
#return groups of mean, FC, log2(FC)
calculate_FC <- function(data, groups, group.number, fun = "mean", fearture.name = "gene"){
  if (fun == "mean") {
    res.means <- calculate_mean(data, groups, fun, fearture.name)
    name.groups <- colnames(res.means)[-1];name.groups <- stringr::str_sub(name.groups, 6) #get group names
  }else if (fun == "median"){
    res.means <- calculate_mean(data, groups, fun, fearture.name)
    name.groups <- colnames(res.means)[-1];name.groups <- stringr::str_sub(name.groups, 8) #get group names
  }else{
      stop("fun should be either mean or median")
  }

  if(group.number == 2){
    res.FCs <- as.data.frame(res.means)
    res.FC <- res.FCs[[2]]/res.FCs[[3]]
    res.FC.log <- log2(res.FC)
    res <- cbind(res.FCs, res.FC, res.FC.log)
    colnames(res) <- c(fearture.name, name.groups, "FC", "log2FC")

  }else if(group.number > 2){
    res.means.l <- res.means[-1]
    levs <- levels(as.factor(groups))
    res.index <- pairwise.compare.index(groups)

    res.FCs <- data.frame(fearture.name = res.means[1])
    for (c in 1:length(res.index)) {
      i <- res.index[[c]][1];j <- res.index[[c]][2]
      res.FC <- res.means.l[i]/res.means.l[j] %>% as.data.frame()
      colnames(res.FC) <- paste0("FC.", levs[i], "_", levs[j])
      res.FC.log <- log2(res.FC) %>% as.data.frame()
      colnames(res.FC.log) <- paste0("log2(FC).", levs[i], "_", levs[j])
      res.FCs <- cbind(res.FCs, res.FC, res.FC.log)
    }
    res <- dplyr::left_join(res.means, res.FCs)
  }

  return(res)
}

pairwise.compare.index <- function(groups){#n number of levels
  lev <- levels(as.factor(groups))
  n <- length(lev)
  #produce compare index
  t <- 1
  res.index <- list()
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      res.index[[t]] <- c(i,j)
      t <- t+1
    }
  }
  return(res.index)
}
