MCMC_sample_par <- function(obj) {
  
  
  sim_stats <- list()
  chains = obj$chains
  for (num_chains in 1:chains)
  {
    sim_stats[[num_chains]] <- par_sim_fun(obj = obj, chain_var = num_chains)  
  }
  
  f.test.data <- list()
  for (num_chains in 1:chains)
  {
    f.data.temp <- data.frame(sim_stats[[num_chains]]$stat_matrix)
    f.data.temp <- f.data.temp[(floor(0.5*nrow(f.data.temp)) + 1):nrow(f.data.temp),]
    f.test.data[[num_chains]] <- f.data.temp
  }
  
  f.test.data <- do.call(bind_rows, f.test.data)
  chain_test <- sort(rep(1:chains, ceiling(0.5*obj$sim$interval)))
  f.test.check <<- apply(f.test.data, 2, function(x, chain_mix = chain_test)
  {
    sum_model <- aov(x ~ chain_mix)
    p.value <- summary(sum_model)[[1]]$'Pr(>F)'[1]
    return(p.value)
  })
  
  if (sum(is.na(f.test.check)) != length(f.test.check))
  {
    if (obj$verbose > 0)
    {
      if (min(f.test.check) < 0.05)
      {
        cat("\n\n")
        cat("Chains did not mix. Changing the MPLE prior or increasing the interval/burnin may improve the fit.")
      }  
    }
    stat_matrix_all <- list()
    stat_random_all  <- list()
    for (check_data in 1:length(sim_stats))
    {
      stat_matrix_all[[check_data]] <- data.frame(sim_stats[[check_data]]$stat_matrix)
      stat_random_all[[check_data]]  <- data.frame(sim_stats[[check_data]]$sim_random)
    }
    
  }
  
  
  
  obj$sim$stats <- do.call(bind_rows, stat_matrix_all)
  obj$sim$re_stats <- do.call(bind_rows, stat_random_all)
  names(obj$sim$re_stats) <- names(summary(obj$form_re))
  obj$net$obs_stats <- sim_stats[[1]]$obs_stats
  return(obj)
}


par_sim_fun <- function(obj, chain_var) {
  require(tidygraph)
  
  # Setup net and formula for ergm::simulate simulation
  cur_theta <- obj$est$chat
  vertex_df <- obj$vertex_data
  #stat_matrix <- matrix(0, nrow = obj$sim$num_obs,
  #                     ncol = obj$net$num_terms)
  
  
  if (sum(colnames(vertex_df) == "names") > 1)
  {
    name_data <- vertex_df[,colnames(vertex_df) == "names"]
    name_data <- data.frame(names = name_data[,1])
    vertex_df <- vertex_df[,colnames(vertex_df) != "names"]
    vertex_df <- data.frame(name_data, vertex_df)
  }
  
  sim_data = obj$hierarchical_data 
  sim_output <- matrix(NA, ncol = length(fixef(obj$est$chat)), nrow = (obj$sim$burnin + obj$sim$interval))  
  sim_random <- matrix(NA, ncol = (length(fixef(obj$est$chat)) + choose(obj$group.count, 2) + obj$group.count - 1), nrow = (obj$sim$burnin + obj$sim$interval))  
  # Simulate sufficient statistics
  
  if (obj$verbose > 0)
  {
    cat("\n\n")
    cat(paste0("Chain ", chain_var, ":\n"))
  }
  
  sing.check = F
  pct_complete <- round(seq(0.02, 1, 0.02)*(obj$sim$burnin + obj$sim$interval))
  sim_data_loop  <- 1
  repeat
  {
    if (obj$verbose > 0)
    {
      if (sim_data_loop == 1)
      {
        cat(paste0(cat("\n",rep("*", round(sim_data_loop/(obj$sim$burnin + obj$sim$interval), 2)*100), sep = ""), "|", round(sim_data_loop/(obj$sim$burnin + obj$sim$interval), 2)*100, "%"),  "\r")
        flush.console()
      } else if (sim_data_loop %in% pct_complete){
        cat(paste0(cat(rep("*", round(sim_data_loop/(obj$sim$burnin + obj$sim$interval), 2)*50), sep = ""), "|", round(sim_data_loop/(obj$sim$burnin + obj$sim$interval), 2)*100, "%"), "\r")
        flush.console()
      }
    }
    
    sing.check = F
    iter = 0
    if (sim_data_loop == (obj$sim$burnin + obj$sim$interval + 1))
    {
      break
    }
    while(!sing.check & iter != 10)
    {
      iter = iter + 1
      ## predict ties
      tie_pred <- predict(cur_theta, sim_data, type = "response")
      
      ## randomly pick ties 
      tie_pred <- sapply(tie_pred, rbinom, n = 1, size = 1)
      
      ## turn the nodes into numeric demarcation
      node_names <- data.frame(names = unique(c(sim_data$Sociality1, sim_data$Sociality2)))
      
      ## turn into edgelist 
      node_ties <- subset(sim_data, tie_pred == 1)
      
      ## subset data
      node_ties <- node_ties[,colnames(node_ties) %in% c("Sociality1", "Sociality2")]
      node_ties <- node_ties[order(node_ties$Sociality1, node_ties$Sociality2),]
      
      colnames(node_ties)[colnames(node_ties) == "Sociality1"] <- "from"
      colnames(node_ties)[colnames(node_ties) == "Sociality2"] <- "to"
      
      simulated_network <- tbl_graph(nodes = node_names, edges = node_ties, directed = is.directed(obj$net$net)) 
      
      suppressMessages(simulated_network %>%
                         activate(nodes) %>%
                         inner_join(vertex_df) %>%
                         intergraph::asNetwork(.) ->> simulated_network)
      
      sim_output[sim_data_loop,] <- summary(obj$form_mcmc) 
      sim_random[sim_data_loop,] <- summary(obj$form_re) 
      
      sim_data$Y <- tie_pred 
      
      suppressWarnings(cur_theta_temp <<- suppressMessages(try({
        bglmer(obj$form_sim,
               data = sim_data,
               family = binomial,
               fixef.prior = normal(sd = c(obj$mcmc.prior, obj$mcmc.prior)))
      }, silent = T)))
      
      if(class(cur_theta_temp) != "try-error")
      {
        if (length(coef(cur_theta_temp)$Group_ID) == ncol(sim_output))
        {
          cur_theta <<- cur_theta_temp
          sim_data_loop <- sim_data_loop + 1
          btw.group.theta <- data.frame(VarCorr(cur_theta))
          if (!isSingular(cur_theta))
          {
            sing.check <- T
          }  
        }
      }
    }
    if (iter == 20)
    {
      break
    }
  }
  rm(cur_theta)
  colnames(sim_random) <- names(summary(obj$form_re))
  rm(simulated_network)
  if (iter == 20)
  {
    stop(paste0("The MCMC failed to mix"))
  } else {
    sim_output <- sim_output[-c(1:obj$sim$burnin),]
    colnames(sim_output) <- names(summary(obj$form_mcmc))
    stat_matrix = sim_output
    if (is.curved(obj$form_net)) {
      num_curved <- sum(obj$net$model$etamap$canonical == 0) / 2
      mod_temp <- ergm_model(form, cur_net)
      cur_curved_ind <- numeric(0)
      curved_ind <- numeric(0)
      for (cur_t in 1:num_curved) { 
        curved_ind_ <- obj$net$model$etamap$curved[[cur_t]]$to
        curved_ind <- c(curved_ind, curved_ind_)
        
        cur_curved_ind_ <- mod_temp$etamap$curved[[cur_t]]$to
        cur_curved_ind <- c(cur_curved_ind, cur_curved_ind_)
        
        cur_len <- length(cur_curved_ind_)
        
        stat_matrix[ , curved_ind_[1:cur_len]] <- cur_sim_stat[ , cur_curved_ind_] 
      }
      stat_matrix[ , -curved_ind] <- cur_sim_stat[ , -cur_curved_ind]
    } 
    
    # Impute observed sufficient statistics if missing data
    obs_ <- summary(obj$form_net)
    
    if (is.curved(obj$form_net)) {
      for (cur_t in 1:num_curved) { 
        curved_ind_ <- obj$net$model$etamap$curved[[cur_t]]$to
        cur_curved_ind_ <- mod_temp$etamap$curved[[cur_t]]$to
        
        cur_len <- length(cur_curved_ind_)
        
        obs_stats[curved_ind_[1:cur_len]] <- obs_[cur_curved_ind_]
      }
      obs_stats[-curved_ind] <- obs_[-cur_curved_ind]
    } else {
      obs_stats <-  obs_
    }
  }
  
  stat_list <- list(stat_matrix = as.matrix(stat_matrix), sim_random = sim_random, obs_stats = obs_)
  return(stat_list)
}