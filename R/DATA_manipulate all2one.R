#DATA_manipulate all2one


DATA_manipulate_all2one <- function(df, by = ";"){
  id <- unique(df[[1]])
  df.res <- data.frame("id"=id)
  counter <- 1
  for (i in id) {
    df.id <- df[[2]][df[[1]]==i]
    df.res1 <- paste(df.id,collapse = by)
    df.res$multi[counter] <- df.res1
      counter <- counter+1
  }
  return(df.res)
}
