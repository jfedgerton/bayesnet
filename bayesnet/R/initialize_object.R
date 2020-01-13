initialize_object <- function(net, 
                              theta_init,
                              sim_param,
                              est_param,
                              verbose,
                              parameterization, 
                              form_mple,
                              check_curve_form, 
                              group.count,
                              hierarchical_data,
                              form_mcmc,
                              form_net,
                              form_sim,
                              form_re,
                              chains,
                              mcmc.prior) {
  
  ##### Create object top structure
  obj <- list(net = NULL,
              sim = NULL,
              est = NULL)
  obj$verbose <- verbose; 
  
  ##### Construct obj$net
  obj$net <- list(model = form_mple,
                  net = net,
                  directed_flag = is.directed(net), 
                  na_flag = FALSE,
                  obs_stats = NULL,
                  obs_stats_step = NULL,
                  theta_names = NULL)
  
  # Fill in the number of nodes in each cluster
  
 
  
  
  term_data <- unlist(strsplit(as.character(form_mple), split = "\\+"))
  theta_names <- term_data[!(term_data %in% c("Y", "~"))]
  
  
  # Set the theta_names for printing 
  theta_names <- data.frame(num_terms = length(theta_names), curved = !is.curved(check_curve_form))
  
  
  obj$net$theta_names <- term_data  
  obj$net$num_terms   <- theta_names$num_terms
  obj$net$is_curved   <- theta_names$curved
  
  ##### Construct obj$sim
  obj$sim <- sim_param
  obj$sim$seq <- TRUE
  obj$sim$stats <- matrix(0, nrow = obj$sim$num_obs,
                          ncol = obj$net$num_terms)
  if (obj$net$na_flag) {
    obj$sim$cond_stats <- obj$sim$stats
  }
  
  
  
  ##### Construct obj$est 
  obj$est           <- est_param
  obj$est$theta     <- theta_init
  obj$est$theta_0   <- theta_init
  obj$est$gamma     <- 1
  obj$est$score_val <- NULL
  
  obj$est$parameterization <- parameterization;
  
  obj$group.count <- group.count
  obj$hierarchical_data <- hierarchical_data
  obj$check_curve_form <- check_curve_form
  obj$form_mcmc <- form_mcmc
  obj$form_sim <- form_sim
  obj$form_net <- form_net
  obj$form_re <- form_re
  obj$chains <- chains
  obj$vertex_data <- node_data(net = obj$net$net)
  obj$directed <- is.directed(obj$net$net)
  obj$mcmc.prior <- mcmc.prior
  
  return(obj)
}
