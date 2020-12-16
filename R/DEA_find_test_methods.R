#pre-test step for determining which statistic test methods should be used to compare means of groups (1, 2, >2 groups)

#
parametric_test <- function(data){

}


# example: ckeck_normality(iris$Sepal.Length), ckeck_normality(iris$Sepal.Length, method = "shapiro.test")
#check whether data is normally distributed by visual inspection [normal plots (histogram), Q-Q plot (quantile-quantile plot)] or significance tests [Shapiro-Wilk test].
ckeck_normality <- function(value, method = "qqplot", p.cutoff = 0.01){
  if(method == "densityplot"){
    ggpubr::ggdensity(value, fill = "lightgray")

  }else if(method == "qqplot"){
    ggpubr::ggqqplot(value)

  }else if(method == "shapiro.test"){
    res <- shapiro.test(value)
    if (res[["p.value"]] > p.cutoff){
      cat("p=", res[["p.value"]], ", data is normally distributed")
    }else cat("p=", res[["p.value"]], ", data isn't normally distributed")

  }
}

#Test of Homogeneity of Variances: bartlett.test (parametric), fligner.test (non-parametric)
#bartlett.test()
#fligner.test()
