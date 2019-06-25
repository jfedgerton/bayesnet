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
  group.var.names <- group.var.names[!(group.var.names %in% c("ingroup", "Group_ID"))]
  
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
    if (grepl("\\+ingroup", formula.stage.1) | substr(formula.stage.1, 1, 7) == "ingroup")
    {
      if (grepl("\\+ingroup", formula.stage.1))
      {
        formula.stage.1 = gsub("\\+ingroup", "", formula.stage.1)  
      } else
      {
        formula.stage.1 <- substr(formula.stage.1, 9, nchar(formula.stage.1))
      }
      
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
      nm.group <- c("Y", names(dta.array.group$predictor[tail.group, head.group, ]),
                    "Group1", "Group2")
      colnames(dta.group) <- nm.group
      
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
      
      Grp <- to_indicator(dta.group[,c("Group1", "Group2", "Group_ID")], "Group")
      dta.group[, "Group1"] <- Grp[,1]
      dta.group[, "Group2"] <- Grp[,2]
      dta.group[, "Group_ID"] <- Grp[,3]
      
      ## Fix factor variables
      if (sum(grepl("nodefactor", formula.stage.1)) > 0)
      {
        fix_factor_vars <- data.frame(dta.group[,grepl("nodefactor", colnames(dta.group))])
        colnames(fix_factor_vars) <- colnames(dta.group)[grepl("nodefactor", colnames(dta.group))]
        factor_vars_form1 <- colnames(fix_factor_vars)
        factor_vars_form1 <- gsub("nodefactor", "", factor_vars_form1)
        factor_vars_form1 <- gsub("[[:digit:]]+", "", factor_vars_form1)
        factor_vars_form1 <- unique(gsub("\\.", "", factor_vars_form1))
        
        new_factor_variable_list <- list()
        for (col_check_factor in 1:length(factor_vars_form1))
        {
          factor_permutations <- expand.grid(unique(group.data$Group_ID), unique(group.data$Group_ID))
          factor_permutations$Var1 <- as.numeric(as.character(factor_permutations$Var1))
          factor_permutations$Var2 <- as.numeric(as.character(factor_permutations$Var2))
          factor_permutations$Keep <- factor_permutations$Var1 <= factor_permutations$Var2
          factor_permutations <- subset(factor_permutations, Keep == T)
          factor_permutations <- factor_permutations[,colnames(factor_permutations) != "Keep"]
          colnames(factor_permutations) <- c("Group1", "Group2")
          original_data_to_factorize <- group.data[,colnames(group.data) %in% c("Group_ID", factor_vars_form1[col_check_factor])]
          factor_permutations <- merge(factor_permutations,
                                       original_data_to_factorize,
                                       by.x = "Group1",
                                       by.y = "Group_ID",
                                       all = T)
          
          colnames(factor_permutations)[3] <- "Group1_Factor"
          
          factor_permutations <- merge(factor_permutations,
                                       original_data_to_factorize,
                                       by.x = "Group2",
                                       by.y = "Group_ID",
                                       all = T)
          
          colnames(factor_permutations)[4] <- "Group2_Factor"
          
          factor_levels_to_create <- sort(unique(c(factor_permutations$Group1_Factor, factor_permutations$Group2_Factor)))
          factor_levels_to_create <- factor_levels_to_create[2:length(factor_levels_to_create)]
          
          new_col_names <- paste0("nodefactor.", factor_vars_form1[col_check_factor], ".",  factor_levels_to_create)
          
          factorial_levels_list <- list()
          for (check_rows in 1:length(factor_levels_to_create))
          {
            var_1_fac <- factor_permutations$Group1_Factor == factor_levels_to_create[check_rows]
            var_2_fac <- factor_permutations$Group2_Factor == factor_levels_to_create[check_rows]
            
            factorial_levels_list[[check_rows]] <- rowSums(cbind(var_1_fac, var_2_fac))
          }
          
          factorial_output <- do.call(cbind, factorial_levels_list)
          colnames(factorial_output) <- new_col_names
          factor_permutations <- cbind(factor_permutations, factorial_output)
          factor_permutations <- factor_permutations[,colnames(factor_permutations) %in% c("Group1", "Group2", new_col_names)]
          
          new_factor_variable_list[[col_check_factor]] <- factor_permutations
        }
        factor_permutations <- factor_permutations[,colnames(factor_permutations) %in% c("Group1", "Group2")]
        for (list_length in 1:length(new_factor_variable_list))
        {
          factor_permutations <- merge(factor_permutations,
                                       new_factor_variable_list[[list_length]],
                                       by = c("Group1", "Group2"))
        }
        
        dta.group <- dta.group[,!(colnames(dta.group) %in% colnames(fix_factor_vars))]
        factor_permutations$Group1 <- paste0("Group", factor_permutations$Group1)
        factor_permutations$Group2 <- paste0("Group", factor_permutations$Group2)
        dta.group <- merge(dta.group,
                           factor_permutations,
                           by = c("Group1", "Group2"))
        
        dta.group <- arrange(dta.group,
                             Group1, Group2)
      }
      
      ## Creats and in and outgroup variable.
      dta.group$ingroup <- 0
      dta.group$ingroup[dta.group$Group1 == dta.group$Group2] <- 1
      
    } else 
    {
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
      nm.group <- c("Y", names(dta.array.group$predictor[tail.group, head.group, ]),
                    "Group1", "Group2")
      colnames(dta.group) <- nm.group
      
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
      
      Grp <- to_indicator(dta.group[,c("Group1", "Group2", "Group_ID")], "Group")
      dta.group[, "Group1"] <- Grp[,1]
      dta.group[, "Group2"] <- Grp[,2]
      dta.group[, "Group_ID"] <- Grp[,3]
      
      ## Fix factor variables
      if (sum(grepl("nodefactor", formula.stage.1)) > 0)
      {
        fix_factor_vars <- data.frame(dta.group[,grepl("nodefactor", colnames(dta.group))])
        colnames(fix_factor_vars) <- colnames(dta.group)[grepl("nodefactor", colnames(dta.group))]
        factor_vars_form1 <- colnames(fix_factor_vars)
        factor_vars_form1 <- gsub("nodefactor", "", factor_vars_form1)
        factor_vars_form1 <- gsub("[[:digit:]]+", "", factor_vars_form1)
        factor_vars_form1 <- unique(gsub("\\.", "", factor_vars_form1))
        
        new_factor_variable_list <- list()
        for (col_check_factor in 1:length(factor_vars_form1))
        {
          factor_permutations <- expand.grid(unique(group.data$Group_ID), unique(group.data$Group_ID))
          factor_permutations$Var1 <- as.numeric(as.character(factor_permutations$Var1))
          factor_permutations$Var2 <- as.numeric(as.character(factor_permutations$Var2))
          factor_permutations$Keep <- factor_permutations$Var1 <= factor_permutations$Var2
          factor_permutations <- subset(factor_permutations, Keep == T)
          factor_permutations <- factor_permutations[,colnames(factor_permutations) != "Keep"]
          colnames(factor_permutations) <- c("Group1", "Group2")
          original_data_to_factorize <- group.data[,colnames(group.data) %in% c("Group_ID", factor_vars_form1[col_check_factor])]
          factor_permutations <- merge(factor_permutations,
                                       original_data_to_factorize,
                                       by.x = "Group1",
                                       by.y = "Group_ID",
                                       all = T)
          
          colnames(factor_permutations)[3] <- "Group1_Factor"
          
          factor_permutations <- merge(factor_permutations,
                                       original_data_to_factorize,
                                       by.x = "Group2",
                                       by.y = "Group_ID",
                                       all = T)
          
          colnames(factor_permutations)[4] <- "Group2_Factor"
          
          factor_levels_to_create <- sort(unique(c(factor_permutations$Group1_Factor, factor_permutations$Group2_Factor)))
          factor_levels_to_create <- factor_levels_to_create[2:length(factor_levels_to_create)]
          
          new_col_names <- paste0("nodefactor.", factor_vars_form1[col_check_factor], ".",  factor_levels_to_create)
          
          factorial_levels_list <- list()
          for (check_rows in 1:length(factor_levels_to_create))
          {
            var_1_fac <- factor_permutations$Group1_Factor == factor_levels_to_create[check_rows]
            var_2_fac <- factor_permutations$Group2_Factor == factor_levels_to_create[check_rows]
            
            factorial_levels_list[[check_rows]] <- rowSums(cbind(var_1_fac, var_2_fac))
          }
          
          factorial_output <- do.call(cbind, factorial_levels_list)
          colnames(factorial_output) <- new_col_names
          factor_permutations <- cbind(factor_permutations, factorial_output)
          factor_permutations <- factor_permutations[,colnames(factor_permutations) %in% c("Group1", "Group2", new_col_names)]
          
          new_factor_variable_list[[col_check_factor]] <- factor_permutations
        }
        factor_permutations <- factor_permutations[,colnames(factor_permutations) %in% c("Group1", "Group2")]
        for (list_length in 1:length(new_factor_variable_list))
        {
          factor_permutations <- merge(factor_permutations,
                                       new_factor_variable_list[[list_length]],
                                       by = c("Group1", "Group2"))
        }
        
        dta.group <- dta.group[,!(colnames(dta.group) %in% colnames(fix_factor_vars))]
        factor_permutations$Group1 <- paste0("Group", factor_permutations$Group1)
        factor_permutations$Group2 <- paste0("Group", factor_permutations$Group2)
        dta.group <- merge(dta.group,
                           factor_permutations,
                           by = c("Group1", "Group2"))
        
        dta.group <- arrange(dta.group,
                             Group1, Group2)
      }
    }
    return(dta.group)
  }
}

