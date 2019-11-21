prepare_meergm_data <- function(net, group.data, form, verbose=FALSE)
{
  if (network::is.directed(net) == F) 
  {
    if (isTRUE(verbose)) cat("\n## Preparing MEERGM dataset...")
    nodes <- nrow(as.matrix(net))
    ndyads <- network::network.dyadcount(net)
    form_temp = form
    form <- stats::as.formula(paste("net ~", form))
    
    if (isTRUE(verbose)) cat("\n##   building array...")
    dta.array <- ergm::ergmMPLE(form, output="array", maxMPLEsamplesize=+Inf,
                                control=ergm::control.ergm(MPLE.max.dyad.types=ndyads*10))
    
    if (verbose > 0) cat("\n##   building data.frame...")
    ncoef <- length(dta.array$predictor[1,2,])
    dta <- matrix(0, nrow=ndyads, ncol=5+ncoef)
    group_var <- names(net$val[[1]])[grepl('group', names(net$val[[1]]), ignore.case = T)]
    group <- ergm::get.node.attr(net,group_var)
    
    
    
    
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
    
    if (sum(is.null(names(dta.array$predictor[tail, head, ]))) > 0)
    {
      nm <- c("Y", form_temp, "Group1", "Group2",
              "Sociality1", "Sociality2")
    }
    names(dta) <- nm
    
    if (isTRUE(verbose)) cat("\n##   setting group effects indicators...\n")
    

      #Soc <- to_indicator(dta[,c("Sociality1", "Sociality2")], "Node")
      #Grp <- to_indicator(dta[,c("Group1", "Group2")], "Group")
      #dta[, "Sociality1"] <- Soc[,1]
      #dta[, "Sociality2"] <- Soc[,2]
      #dta[, "Group1"] <- Grp[,1]
      #dta[, "Group2"] <- Grp[,2]
  
  } else {
    if (verbose > 0) cat("\n## Preparing MEERGM dataset...")
    nodes <- nrow(as.matrix(net))
    ndyads <- network::network.dyadcount(net)
    form <- stats::as.formula(paste("net ~", form))
    
    if (isTRUE(verbose)) cat("\n##   building array...")
    dta.array <- ergm::ergmMPLE(form, output="array", maxMPLEsamplesize=+Inf,
                                control=ergm::control.ergm(MPLE.max.dyad.types=ndyads*10))
    
    if (isTRUE(verbose)) cat("\n##   building data.frame...")
    ncoef <- length(dta.array$predictor[1,2,])
    group_var <- names(net$val[[1]])[grepl('group', names(net$val[[1]]), ignore.case = T)]
    group <- ergm::get.node.attr(net,group_var)
    
    
    response <- expand.grid(colnames(dta.array$response), colnames(dta.array$response))
    response <- subset(response, Var1 != Var2)
    
    colnames(response) <- c("head", "tail")
    
    dta <- matrix(nrow = nrow(response), ncol = 5+ncoef)
    
    for (directed_dyads in 1:nrow(response))
    {
      temp_response <- dta.array$response[rownames(dta.array$response) == response[directed_dyads,1],
                                             colnames(dta.array$response) == response[directed_dyads,2]]
      
      temp_predictors <- dta.array$predictor[rownames(dta.array$predictor) == response[directed_dyads,1],
                                             colnames(dta.array$predictor) == response[directed_dyads,2], ]
      
      group_1 <- group[response[directed_dyads,1]]
      group_2 <- group[response[directed_dyads,2]]
      
      temp_response <- as.vector(temp_response)
      temp_predictors <- as.vector(temp_predictors)
      group_1 <- as.vector(group_1)
      group_2 <- as.vector(group_2)
      head_tail <- as.vector(response[directed_dyads,])
      dta[directed_dyads,] <- c(temp_response, temp_predictors, group_1, group_2, response[directed_dyads,1], response[directed_dyads,2])
    }
    
    
    
    dta <- data.frame(dta)
    nm <- c("Y", names(dta.array$predictor[1, 2, ]), "Group1", "Group2",
            "Sociality1", "Sociality2")
    
    if (sum(is.null(names(dta.array$predictor[tail, head, ]))) > 0)
    {
      nm <- c("Y", form_temp, "Group1", "Group2",
              "Sociality1", "Sociality2")
    }
    names(dta) <- nm
    
    names(dta) <- nm
    
    for (dta_cols in 1:ncol(dta))
    {
      dta[,dta_cols] <- as.numeric(as.character(dta[,dta_cols]))
    }
    dta <- arrange(dta, Sociality1, Sociality2, Group1, Group2)
    if (isTRUE(verbose)) cat("\n##   setting random effects indicators...\n")
    #Soc <- to_indicator(dta[,c("Sociality1", "Sociality2")], "Node")
    #Grp <- to_indicator(dta[,c("Group1", "Group2")], "Group")
    #dta[, "Sociality1"] <- Soc[,1]
    #dta[, "Sociality2"] <- Soc[,2]
    #dta[, "Group1"] <- Grp[,1]
    #dta[, "Group2"] <- Grp[,2]
    
  }
  
  dta
}
