create_group_data <- function(group.data, formula.stage.1)
{
  require('ergm')
  require('dplyr')
  ## Turn the group data into inter and intra group data
  group.nodes <- nrow(as.matrix(group.data))
  
  ## Create an adjacency matrix for group data
  adj.mat <- matrix(0, nrow = group.nodes, ncol = group.nodes)
  
  ## Names of all the groups
  colnames(adj.mat) <- rownames(adj.mat) <- group.data$Group_ID
  
  ## For all groups
  group.data.for.analysis <- network::as.network(adj.mat, directed = F)
  
  ## Turn data into a network
  group.dyads <- network::network.dyadcount(group.data.for.analysis)
  group.var.names <- colnames(group.data)
  group.var.names <- group.var.names[!(group.var.names %in% c("ingroup"))]
  
  ## Add in the network attributes
  for (vert.atts in 1:length(group.var.names))
  {
    group.data.for.analysis %v% group.var.names[vert.atts] <- group.data[,group.var.names[vert.atts]]
  }
  
  if (
    ## Variables that cannot work for group data
    grepl("absdiffcat", formula.stage.1) == T |
    grepl("transitive", formula.stage.1) == T |
    grepl("triad", formula.stage.1) == T |
    grepl("atleast", formula.stage.1) == T |
    grepl("atmost", formula.stage.1) == T |
    grepl("degree", formula.stage.1) == T |
    grepl("dgwe", formula.stage.1) == T |
    grepl("edgecov", formula.stage.1) == T |
    grepl("kstar", formula.stage.1) == T |
    grepl("absdiffcat", formula.stage.1) == T)
  {
    cat("The program does not support one or more stage 1 coefficients.")
  } else {
    
    ## Formula for group data
    
    formula.stage.1 = gsub(" ", "", formula.stage.1)
    group_form <- stats::as.formula(paste0("group.data.for.analysis ~", formula.stage.1))
      dta.array.group <- ergm::ergmMPLE(group_form, output="array", maxMPLEsamplesize=+Inf,
                                        control=ergm::control.ergm(MPLE.max.dyad.types=group.dyads*10))
      
      ncoef.group <- length(dta.array.group$predictor[1,2,])
      dta.group <- matrix(0, nrow=group.dyads, ncol=3 + ncoef.group)
      
      ## Creating group data
      idx.group <- 1
      for (tail.group in 1:(group.nodes-1)) {
        for (head.group in (tail.group+1):group.nodes) {
          dta.group[idx.group,] <- c(dta.array.group$response[tail.group, head.group],
                                     dta.array.group$predictor[tail.group, head.group, ],
                                     tail.group,
                                     head.group)
          idx.group <- idx.group+1
        }
      }
      ## Name the group data
      dta.group <- data.frame(dta.group)
      if (ncoef.group != 1)
      {
        nm.group <- c("Y", names(dta.array.group$predictor[tail.group, head.group, ]),
                      "Group1", "Group2")
        colnames(dta.group) <- nm.group  
      } else 
      {
        var_name <- formula.stage.1
        var_name <- gsub("\\(", "", var_name)
        var_name <- gsub("\\)", "", var_name)
        var_name <- gsub("'", "\\.", var_name)
        var_name <- substr(var_name, 1, nchar(var_name) - 1)
        nm.group <- c("Y", var_name,
                      "Group1", "Group2")
        colnames(dta.group) <- nm.group
      }
      
      
      ## Creating group data
      all_groups <- unique(c(dta.group$Group1,dta.group$Group2))
      self_loop_data <- data.frame(matrix(NA, ncol = ncol(dta.group), nrow = length(all_groups)))
      colnames(self_loop_data) <- nm.group
      self_loop_data$Y <- 0
      self_loop_data$Group1 <- all_groups
      self_loop_data$Group2 <- all_groups
      
      ## Self loop data for groups
      for (self_loop in 2:ncol(self_loop_data))
      {
        if (grepl("nodematch.", colnames(self_loop_data)[self_loop]) == T)
        {
          self_loop_data[,colnames(self_loop_data)[self_loop]] <- 1
        } else if (grepl("absdiff", colnames(self_loop_data)[self_loop]) == T)
        {
          self_loop_data[,colnames(self_loop_data)[self_loop]] <- 0
        }
      }
      
      ## Add in the self loop data
      dta.group <- dplyr::bind_rows(dta.group,
                                    self_loop_data)
      dta.group <- dplyr::arrange(dta.group, Group1, Group2)
      dta.group <- dplyr::select(dta.group, -Y)
      dta.group <- dplyr::mutate(dta.group, intercept = 1)
      dta.group$Group_ID <- 1:nrow(dta.group)
      
      #Grp <- to_indicator(dta.group[,c("Group1", "Group2", "Group_ID")], "Group")
      #dta.group[, "Group1"] <- Grp[,1]
      #dta.group[, "Group2"] <- Grp[,2]
      #dta.group[, "Group_ID"] <- Grp[,3]
      
      
      if (sum(grepl("nodefactor", formula.stage.1)) > 0)
      {
        dta.group$ID_temp <- 1:nrow(dta.group)
        dta.group.factor <- dta.group[,grepl(paste(c("Group_ID", "Group1", "Group2", "ID_temp", "nodefactor"), collapse = "|"), colnames(dta.group))]
        dta.group.factor <- dta.group.factor[,!(colnames(dta.group.factor) %in% "nodematch.Group_ID")]
        
        factor_list_output <- list()
        for (fact_vars in 1:(ncol(dta.group.factor) - 4))
        {
          dta.group.factor.temp <- dta.group.factor[,colnames(dta.group.factor) %in% c(colnames(dta.group.factor)[fact_vars], "Group_ID", "Group1", "Group2", "ID_temp")]
          dta.group.factor.temp.missing <- subset(dta.group.factor.temp, is.na(dta.group.factor.temp[,1]))
          dta.group.factor.temp.missing$G1 <- as.numeric(gsub("[^\\d]+", "", dta.group.factor.temp.missing$Group1, perl=TRUE))
          dta.group.factor.temp.missing$G2 <- as.numeric(gsub("[^\\d]+", "", dta.group.factor.temp.missing$Group2, perl=TRUE))
          
          factor_level <- colnames(dta.group.factor)[fact_vars]
          factor_level <- gsub("nodefactor", "", factor_level)
          factor_level <- gsub("\\.", " ", factor_level)
          factor_level <- trimws(factor_level)
          factor_level <- gsub(".*\\s","",factor_level)
          factor_variable <- colnames(dta.group.factor)[fact_vars]
          factor_variable <- gsub("nodefactor", "", factor_variable)
          factor_variable <- gsub("\\.", " ", factor_variable)
          factor_variable <- trimws(factor_variable)
          factor_variable <- gsub("\\s.*","",factor_variable)
          
          group.temp <- data.frame(group.data[,colnames(group.data) %in% c("Group_ID", factor_variable)])
          
          if (factor_variable == "Group_ID")
          {
            colnames(group.temp)[1] <- "Group_ID"
          }
          group.temp[,factor_variable] <- as.character(group.temp[,factor_variable])
          for (missing_data in 1:nrow(dta.group.factor.temp.missing))
          {
            dta.group.factor.temp.missing[missing_data,1] <- 2*sum(group.temp[group.temp$Group_ID == dta.group.factor.temp.missing$G1[missing_data],factor_variable] == factor_level)
          }
          factor_list_output[[fact_vars]] <- dta.group.factor.temp.missing
        }
        
        if (length(factor_list_output) > 1)
        {
          for (loop_fact in 2:length(factor_list_output))
          {
            factor_list_output[[1]] <- merge(factor_list_output[[1]],
                                             factor_list_output[[loop_fact]],
                                             by = c("Group1", "Group2", "Group_ID", "ID_temp", "G1", "G2"))
          }
        }
        
        dta.group.merge.data <- dta.group[,grepl("nodefactor", colnames(dta.group)) == F]
        factor_list_output[[1]] <- merge(factor_list_output[[1]], 
                                         dta.group.merge.data, 
                                         by = c("Group_ID", "Group1", "Group2", "ID_temp"),
                                         all.x = T)
        dta.group <- subset(dta.group, !(ID_temp %in%  factor_list_output[[1]]$ID_temp))
        dta.group <- bind_rows(dta.group, 
                               factor_list_output[[1]])
        dta.group <- dta.group[,!(colnames(dta.group) %in% c("ID_temp", "G1", "G2"))]
      }
     
        
        dta.group <- arrange(dta.group,
                             Group1, Group2)
      }
    
    return(dta.group)
}

