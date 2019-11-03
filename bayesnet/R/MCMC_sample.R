MCMC_sample_par <- function(obj) {
  
  
  
  
  sim_stats <- par_sim_fun(obj = obj, 
                           data = obj$hierarchical_data, 
                           check_curve_form = obj$check_curve_form, 
                           form_mcmc = obj$form_mcmc, 
                           form_sim = obj$form_sim,
                           form_net = obj$form_net,
                           net = net)
  obj$sim$stats <- sim_stats$stat_matrix
  obj$net$obs_stats <- sim_stats$obs_stats
  return(obj)
}


par_sim_fun <- function(obj, data, check_curve_form, form_mcmc, form_sim, form_net, net) {
  require(tidygraph)
  
  # Setup net and formula for ergm::simulate simulation
  cur_theta <- obj$est$chat
    
    #stat_matrix <- matrix(0, nrow = obj$sim$num_obs,
     #                     ncol = obj$net$num_terms)
    vertex_df <- node_data(net = net)
    
    if (sum(colnames(vertex_df) == "names") > 1)
    {
      name_data <- vertex_df[,colnames(vertex_df) == "names"]
      name_data <- data.frame(names = name_data[,1])
      vertex_df <- vertex_df[,colnames(vertex_df) != "names"]
      vertex_df <- data.frame(name_data, vertex_df)
    }
    
    sim_data = data 
    sim_output <- matrix(NA, ncol = ncol(obj$sim$stats), nrow = (obj$sim$burnin + obj$sim$interval))  
    
    # Simulate sufficient statistics
    
    sing.check = F
    pct_complete <- round(seq(0.01, 0.99, 0.01)*(obj$sim$burnin + obj$sim$interval))
    sim_data_loop  <- 1
    repeat
    {
      if (sim_data_loop == 1)
      {
        cat(paste0(cat("\n",rep("*", round(sim_data_loop/(obj$sim$burnin + obj$sim$interval), 2)*100), sep = ""), "|", round(sim_data_loop/(obj$sim$burnin + obj$sim$interval), 2)*100, "%"), "\r")
        flush.console()
      } else if (sim_data_loop %in% pct_complete){
        cat(paste0(cat(rep("*", round(sim_data_loop/(obj$sim$burnin + obj$sim$interval), 2)*100), sep = ""), "|", round(sim_data_loop/(obj$sim$burnin + obj$sim$interval), 2)*100, "%"), "\r")
        flush.console()
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
          
          simulated_network <- tbl_graph(nodes = node_names, edges = node_ties, directed = is.directed(net)) 
          
          suppressMessages(simulated_network %>%
                             activate(nodes) %>%
                             inner_join(vertex_df) %>%
                             intergraph::asNetwork(.) ->> simulated_network)
          
          sim_output[sim_data_loop,] <- summary(form_mcmc) 
          
          sim_data$Y <- tie_pred 
          
          suppressWarnings(cur_theta_temp <<- suppressMessages(try({
            bglmer(form_sim,
                   data = sim_data,
                   family = binomial,
                   cov.prior = NULL,
                   fixef.prior = t())
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
      if (iter == 10)
      {
        break
      }
    }
    rm(cur_theta)
    rm(simulated_network)
    if (iter == 10)
    {
      stop(paste0("The MCMC failed to mix"))
    } else {
      sim_output <- sim_output[-c(1:obj$sim$burnin),]
      colnames(sim_output) <- names(summary(form_mcmc))
      stat_matrix = sim_output
      if (is.curved(form_net)) {
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
      obs_ <- summary(form_net)
        
      if (is.curved(form_net)) {
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
    
    stat_list <- list(stat_matrix = as.matrix(stat_matrix), obs_stats = obs_)
    return(stat_list)
}
