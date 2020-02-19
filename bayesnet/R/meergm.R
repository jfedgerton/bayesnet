meergm <- function(net, 
                   form.stage.1 = NULL,
                   form.stage.2 = NULL,
                   group.data = NULL,
                   parameterization = "standard",
                   options = set_options(),  
                   theta_init = NULL,
                   verbose = 0,
                   chains = 2,
                   eval_loglik = TRUE,
                   seed = 123, 
                   estimation = "MCMC-MLE",
                   prior.var = 1e4, 
                   mcmc.prior = 1, 
                   n_prior = 1e4,
                   df_prior = 1e4) {
  
  require('plyr')
  require('blme')
  require('stringr')
  require('parallel')
  
  net <<- net
  if (grepl("gwesp|gwidegree|gwodegreegwesp", form.stage.2))
  {
    cat("Fixed effects are only supported at this time for curved space parameters.\n")
  }
  
  group_data <- create_group_data(group.data = group.data, formula.stage.1 = form.stage.1)
  group_data2 <- group_data
  group_data2$Group1 <- group_data$Group2
  group_data2$Group2 <- group_data$Group1
  group_data <- rbind(group_data, 
                      group_data2)
  group_data <- distinct(group_data, Group1, Group2, Group_ID, .keep_all = T)
  model_data <- prepare_meergm_data(net = net, group.data = group_data, form = form.stage.2)
  hierarchical_data <- merge(group_data,
                             model_data, 
                             by = c("Group1", "Group2"))
  
  hierarchical_data <- arrange(hierarchical_data, 
                               Sociality1, Sociality2)
  # If a seed is provided, set it
  if (!is.null(seed)) { 
    check_integer(seed, "seed")
    set.seed(seed, "L'Ecuyer")
  } 
  group_var <- names(net$val[[1]])[grepl("group", names(net$val[[1]]), ignore.case = T)]
  node_memb <- get.vertex.attribute(net, group_var)
  group.count <- length(unique(node_memb))  
  
  
  # Adjust formula if necessary 
  form_mple        <<- adjust_formula_mple(form = form.stage.2, form.stage.1 = form.stage.1) 
  form_sim         <<- adjust_formula_sim(data = hierarchical_data, form_ref = form_mple)
  form_net         <<- adjust_formula_net(form = form.stage.2, form.stage.1 = form.stage.1) 
  form_mcmc        <<- adjust_formula_sim_net(form = form.stage.2, form.stage.1 = form.stage.1)
  form_re          <<- adjust_formula_sim_net_re(form = form.stage.2, form.stage.1 = form.stage.1, net = net)
  check_curve_form <<- adjust_formula_curve(form = form.stage.2, form.stage.1 = form.stage.1)
  # Parse formula to get network and model
  
 
  # Initialize object
  obj <<- initialize_object(net = net, 
                           theta_init = theta_init,
                           sim_param = options$sim_param,
                           est_param = options$est_param,
                           verbose = verbose,
                           parameterization = parameterization,
                           form_mple = form_mple,
                           check_curve_form = check_curve_form,
                           group.count = group.count,
                           hierarchical_data = hierarchical_data,
                           form_mcmc = form_mcmc,
                           form_net = form_net,
                           form_sim = form_sim,
                           form_re = form_re,
                           chains = chains,
                           mcmc.prior = mcmc.prior)
  
  # Remove objects that are no longer needed 
  
  
  # Initialize estimate if an initial estimate is not provided 
  cd_flag <- TRUE
  if (is.null(obj$est$theta)) {
    if (verbose > 0) { 
      cat("\n\nComputing initial estimate.")  
    }
    
    
      obj$est$theta <- numeric(length(obj$net$theta_names)) 
      obj <- compute_initial_estimate(obj)
    
    if (verbose > 0) {
      cat("\n    Initial fixed effect estimate:")
      cat(paste("\n     ", names(obj$est$theta), " = ", formatC(obj$est$theta, digits = 4, format = "f")))
    }
  } else { 
    obj$est$theta <- theta_init
  }
  
  
  if (estimation == "MCMC-MLE")
  {
    if (verbose > 0) { 
      cat("\n\n\nBegining Monte-Carlo maximum likelihood estimation\n")
      cat("===================================================")
      cat("\n")
    }
    
    obj <- MCMLE(obj)
    
    # Estimate the between block model (if possible)
    
    # Evalaute results of MCMLE procedure 
    if (verbose > 0 & !obj$est$ML_status_fail) {
      cat(paste0("\n\nComputing approximate loglikelihood at estimate using ", 
                 obj$est$bridge_num, " bridges."))
      cat("\n\n")
      
    } else if (verbose > 0 & obj$est$ML_status_fail) {
      cat("\n\nEstimation procedure stopping. Estimation unsuccesful.\n\n", call. = FALSE)
    }
    
    
    # Make structure to be returned
    if (!obj$est$ML_status_fail) {
      if (obj$est$inCH_counter > 0) {
        ## Start here
        if (eval_loglik) { 
          obj$likval <- lik_fun(form = form_net, theta = obj$est$theta, 
                                bridge_num = obj$est$bridge_num, ncores = obj$est$par_n_cores,
                                form_net = obj$form_net,
                                offset = obj$est$parameterization == "offset",
                                burnin = obj$sim$bridge_burnin, 
                                interval = obj$sim$bridge_interval, 
                                sample_size = obj$sim$bridge_sample_size) 
        } else { 
          obj$likval <- NULL
          obj$bic <- NULL
        }
        mcmc_path <- obj$sim$stats
        #names(obj$se) <- colnames(obj$sim$stats)
        obj$est$theta <- as.vector(obj$est$theta)
        names(obj$est$theta) <- colnames(obj$sim$stats)
        obj <- compute_se(obj)
        names(obj$se) <- names(obj$est$theta)
        nodes.per.group <- table(get.node.attr(net, group_var))
        dyad.count <- count.edges(nodes.per.group, net = net)
        sum1 <- apply(obj$sim$re_stats[,grepl("mix\\.group", names(obj$sim$re_stats))], 2, mean) #can be any model including nodemix
        
        grp.efct <- sum1[startsWith(names(sum1),"mix")]/dyad.count
        suppressWarnings(grp.efct <- logit(grp.efct))
        mean.grp.efct <- mean(grp.efct)
        center.mean.grp.efct <- obj$est$theta[names(obj$est$theta) == "edges"] - mean.grp.efct
        grp.efct <- grp.efct + center.mean.grp.efct
        group.var <- c(`grand mean` = mean(grp.efct), `between variance` = var(grp.efct))
        names(obj$est$theta)[names(obj$est$theta) == "edges"] <- "grand mean"
        names(obj$se)[names(obj$se) == "edges"] <- "grand mean"
        
        pvalue = compute_pvalue(obj)
        estimates <- list(theta = obj$est$theta,
                          se = obj$se,
                          dyad.group.intercepts = grp.efct,
                          btw.var = group.var[names(group.var) == "between variance"],
                          pvalue = obj$pvalue, 
                          logLikval = obj$likval,  
                          mcmc_chain = mcmc_path,
                          estimation_status = "success",
                          parameterization = obj$est$parameterization,
                          formula = form_net,
                          network = net, 
                          mple.estimate = obj$est$chat,
                          type = "MCMC-MLE")
        class(estimates) <- "meergm" 
        rm(mcmc_path); clean_mem()
      } else if (obj$est$inCH_counter == 0) { 
        cat("\n\nWarning: Maximum number of iterations reached without the observation lying in the")
        cat(" interior of the simulated convex hull. Parameters not estimated.\n\n")
        estimates <- list(theta = NA,
                          se = NA,
                          formula = form_net, 
                          network = net,
                          mcmc_chain = NULL,
                          estimation_status = "failed")
        class(estimates) <- "meergm"
      }
    } else {
      estimates <- list(posterior.theta = NA, 
                        posterior.se = NA,
                        formula = form_net, 
                        network = net, 
                        mcmc_chain = NULL,
                        estimation_status = "failed")
      class(estimates) <- "meergm" 
    }
    
  } else if (estimation == "BOOTSTRAP-MPLE")
  {
    boot_est <- mple_boot(obj = obj, group = group_data, form = form.stage.2)
    obj$theta_est <- boot_est$theta
    obj$se        <- boot_est$se
    obj$boot_fe   <- boot_est$boot_fe
    obj <- compute_pvalue_mple(obj)
    
    estimates <- list(theta = obj$theta_est,
                      se = obj$se,
                      pvalue = obj$pvalue, 
                      estimation_status = "success",
                      parameterization = obj$est$parameterization,
                      formula = form_net,
                      network = net, 
                      mple.estmate = obj$est$chat,
                      btw.var = boot_est$group.re,
                      type = "BOOTSTRAP-MPLE")
    
  } else {
    warning("Improper estimation process.")
  }
  
  # Call MCMLE to perform estimation
  
  return(estimates)
}

