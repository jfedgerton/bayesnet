mple_boot <- function(obj, group, form)
{
  require('tidygraph')
  cur_theta <- obj$est$chat
  vertex_df <- obj$vertex_data
  
  if (sum(colnames(vertex_df) == "names") > 1)
  {
    name_data <- vertex_df[,colnames(vertex_df) == "names"]
    name_data <- data.frame(names = name_data[,1])
    vertex_df <- vertex_df[,colnames(vertex_df) != "names"]
    vertex_df <- data.frame(name_data, vertex_df)
  }
  
  sim_data = obj$hierarchical_data 
  par_simulations <- obj$sim$interval*obj$chains
  pct_complete <- round(seq(0.02, 1, 0.02) * par_simulations)
  sim_boot <- matrix(NA, ncol = (length(fixef(obj$est$chat)) + choose(obj$group.count, 2) + obj$group.count - 1), nrow  = par_simulations)  
  # Simulate sufficient statistics
  cat("\n\n")
  
  if (obj$verbose > 0)
  {
    cat(paste0("Starting parametric bootstrap:\n")) 
  }
  
  sing.check = F
  pct_complete <- round(seq(0.02, 1, 0.02) * par_simulations)
  sim_data_loop  <- 1
  boot_fe <- matrix(ncol = length(names(coef(obj$est$chat)$Group_ID)),
                    nrow = par_simulations)
  boot_re <- matrix(ncol = nrow(coef(obj$est$chat)$Group_ID),
                    nrow = par_simulations)
  
  colnames(boot_fe) <- names(coef(obj$est$chat)$Group_ID)
  colnames(boot_fe)[colnames(boot_fe) == "(Intercept)"] <- "grand mean"
  colnames(boot_re) <- paste0("group", 1:ncol(boot_re))
  group.var <- c()
  iter = 0
  while (iter < par_simulations){
    if (obj$verbose > 0)
    {
      if (iter == 1)
      {
        cat(paste0(cat("\n",rep("*", round(iter/par_simulations, 2)*100), sep = ""), "|", round(iter/par_simulations, 2)*100, "%"),  "\r")
        flush.console()
      } else if (iter %in% pct_complete){
        cat(paste0(cat(rep("*", round(iter/par_simulations, 2)*50), sep = ""), "|", round(iter/par_simulations, 2)*100, "%"), "\r")
        flush.console()
      }
    }
    
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
    
    simulated_network <- tbl_graph(nodes = node_names, edges = node_ties, directed = obj$directed) 
    
    suppressMessages(simulated_network %>%
                       activate(nodes) %>%
                       inner_join(vertex_df) %>%
                       intergraph::asNetwork(.) -> simulated_network)
    par_node <- prepare_meergm_data(net = simulated_network, group.data = group, form = form)
    
    boot_data <- merge(group,
                       par_node,
                       by = c("Group1", "Group2"))
    
    
    suppressWarnings(theta_boot <- suppressMessages(try({
        bglmer(obj$form_sim,
               data = boot_data,
               family = binomial,
               fixef.prior = normal(sd = c(obj$mcmc.prior, obj$mcmc.prior)))
      }, silent = T)))
      
    if (class(theta_boot) != "try-error")
    {
      boot_est_fe <- fixef(theta_boot)
      boot_est_re <- coef(theta_boot)$Group_ID[,colnames(coef(theta_boot)$Group_ID) == "(Intercept)"]
      names(boot_est_fe)[names(boot_est_fe) == "(Intercentp)"] <- "grand mean"
      names(boot_est_re) <- paste0("group", 1:length(boot_est_re))
      iter = iter + 1
      boot_fe[iter,] <- boot_est_fe
      boot_re[iter,] <- boot_est_re
      group.var[iter] <- data.frame(VarCorr(theta_boot))$sdcor
    }
  }
  
  if (obj$verbose > 0)
  {
    cat(paste0(cat(rep("*", round(iter/par_simulations, 2)*50), sep = ""), "|", round(iter/par_simulations, 2)*100, "%"), "\r")
    flush.console()
  }
  
  boot_fe_means <- apply(boot_fe, 2, mean)
  boot_fe_sd    <- apply(boot_fe, 2, sd)
  boot_re_mean  <- apply(boot_re, 2, mean)
  group.var <- mean(group.var)
  names(group.var) <- "between variance"
  theta    <- c(boot_fe_means, group.var)
  theta.sd <- c(boot_fe_sd)
  theta.group <- c(boot_re_mean)
  return.items <- list(theta = theta, se = theta.sd, group.re = theta.group, boot_fe = boot_fe, boot_re = boot_re)
  return(return.items)
}
