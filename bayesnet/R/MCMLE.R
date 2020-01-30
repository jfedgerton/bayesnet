MCMLE <- function(obj) {
  
  # Run the estimation procedure until termination criterions are met
  conv_flag_1 <- FALSE
  conv_flag_2 <- FALSE
  obj$est$inCH_counter <- 0
  while (!obj$est$MCMLE_status & (obj$est$MCMLE_iter <= obj$est$MCMLE_max_iter) & !obj$est$ML_status_fail) {
    
    if (obj$verbose > 0) {
      if (obj$est$MCMLE_iter == 1) {  
        cat("\n")
      } else { 
        cat("\n\n")
      }
      cat(paste0("MCMLE Iteration: ", obj$est$MCMLE_iter))
    }
    
    # Get MCMC sample
    if (obj$verbose > 0) { 
      cat("\n")
      cat("  Sampling networks:")
    } 
    
    
    obj <- MCMC_sample_par(obj)
    
    # Check that the distribution is not degenerate
    degen_flag <- FALSE 
    sd_ <- apply(obj$sim$stats, 2, sd)
    
    sd_check_sum <- diag(length(obj$est$theta)) %*% sd_
    if (sum(sd_check_sum < 0.05) > 0) { 
      msg <- "Standard deviation of some simulated statistics are near zero;" 
      msg <- paste(msg, "estimation cannot continue.\n\n")
      msg <- paste(msg, "Possible reasons include:\n") 
      msg <- paste(msg, "    - The proposed model may be nearly degenerate or poorly specified.\n")
      msg <- paste(msg, "    - The  Markov chain has failed to mix (adjusting simulation parameters may help).")
      cat("\n\n")
      stop(msg, call. = FALSE)
      degen_flag <- TRUE
    }
    
    if (!degen_flag) { 
      # Use stepping algorithm to put observed vector within the convex hull of the sim sufficient statistics 
      obj <- step_to_chull(obj)
      
      if (obj$est$gamma == 1) { 
        obj$est$inCH_counter <- obj$est$inCH_counter + 1
      }
      
      # Optimize log-likelihood approximation
      if (obj$verbose > 0) { 
        cat("\n\n  Optimizing likelihood:")
      }
      obj <- opt_fisher(obj)
      
      # If two consecutive samples have been within the convex hull, accept convergence 
      if (obj$est$inCH_counter == 2) { 
        obj$est$MCMLE_status <- TRUE
      }
      
      if (obj$est$inCH_counter < 2) { 
        obj$est$MCMLE_iter <- obj$est$MCMLE_iter + 1
      }
      
      if ((obj$est$MCMLE_iter <= obj$est$MCMLE_max_iter) & !obj$est$ML_status_fail &
          !obj$est$max_iter_flag) { 
        if (obj$verbose > 0) { 
          cat("\n         - Current estimate: ")
          cat(paste0("\n            ",
                     names(obj$net$obs_stats),
                     " = ",
                     round(obj$est$theta, digits = 4)))
          
          score_val <- sum(abs(obj$est$score_val^2))
          print_1 <- ifelse(score_val < 0.000001, "<.000001", as.character(round(score_val, digits = 6)))
          cat("\n\n         - L2-norm of score function at estimate: ", print_1)
        }
      } else if (obj$est$MCMLE_status) {
        # Bypass this 
      } 
    } else {
      obj$est$ML_status_fail <- TRUE  
    }
  }
  return(obj)
}

