#' @param X
#' @param text

to_indicator <- function(X, text, directed=FALSE)
{
  if (isTRUE(directed)) {
    stop("Directed networks not yet supported.")
  } else {
    for (k in 1:ncol(X)){
      X[,k] <- paste0(text, X[,k])
      factor.levels <- unique(unlist(X))
      X[,k] <- factor(X[,k], levels=factor.levels)
    }
  }
  X
}
