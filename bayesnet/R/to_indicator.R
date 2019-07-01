to_indicator <- function(X, text, directed=FALSE)
{
  for (k in 1:ncol(X)){
      X[,k] <- paste0(text, X[,k])
      factor.levels <- unique(unlist(X))
      X[,k] <- factor(X[,k], levels=factor.levels)
    }
  
  X
}
