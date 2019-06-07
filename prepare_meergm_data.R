#' @importFrom dplyr
#' @importFrom ergm
#' @importFrom fergm
#' @importFrom rstan
#'
#' @param net the outcome network
#' @param group.data the group level data for mixed effect estimation
#' @param form first stage model

prepare_meergm_data <- function(net, group.data, form, verbose=FALSE)
{
  ## Temporary until directed code is added.
  if (network::is.directed(net)) stop("Directed networks not yet supported.")

  if (isTRUE(verbose)) cat("\n## Preparing FERGM dataset...")
  nodes <- nrow(as.matrix(net))
  ndyads <- network::network.dyadcount(net)
  form <- stats::as.formula(paste("net ~", form))

  if (isTRUE(verbose)) cat("\n##   building array...")
  dta.array <- ergm::ergmMPLE(form, output="array", maxMPLEsamplesize=+Inf,
                              control=ergm::control.ergm(MPLE.max.dyad.types=ndyads*10))

  if (isTRUE(verbose)) cat("\n##   building data.frame...")
  ncoef <- length(dta.array$predictor[1,2,])
  dta <- matrix(0, nrow=ndyads, ncol=5+ncoef)
  group <- ergm::get.node.attr(net,"Group")

  idx <- 1
  for (tail in 1:(nodes-1)) {
    for (head in (tail+1):nodes) {
      dta[idx,] <- c(dta.array$response[tail, head],
                     dta.array$predictor[tail, head, ],
                     group[tail],
                     group[head],
                     tail,
                     head)
      idx <- idx+1
    }
  }

  dta <- data.frame(dta)
  nm <- c("Y", names(dta.array$predictor[tail, head, ]), "Group1", "Group2",
          "Sociality1", "Sociality2")
  names(dta) <- nm

  if (isTRUE(verbose)) cat("\n##   setting random effects indicators...\n")
  if (network::is.directed(net)) {
    stop("Directed networks not yet supported.")
  } else {
    Soc <- to_indicator(dta[,c("Sociality1", "Sociality2")], "Node")
    Grp <- to_indicator(dta[,c("Group1", "Group2")], "Group")
    dta[, "Sociality1"] <- Soc[,1]
    dta[, "Sociality2"] <- Soc[,2]
    dta[, "Group1"] <- Grp[,1]
    dta[, "Group2"] <- Grp[,2]
  }

  dta
}
