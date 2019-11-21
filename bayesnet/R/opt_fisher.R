opt_fisher <- function(obj) {
  outer_fun <- function(v) {
    base::outer(v, v)
  }
  # Reset iter counters and status
  obj$est$NR_iter <- 1
  obj$est$NR_status <- FALSE
  obj$est$max_iter_flag <- FALSE
  obj$est$step_err <- 1
  obj$est$NR_step_len_ <- obj$est$NR_step_len
  
  # Compute iterations while termination criterion is not met
  if (obj$verbose > 1) {
    cat("\n")
  }
  
  if (any(is.nan(obj$est$theta))) {
    obj$est$ML_status_fail <- TRUE
  }
  
  while (!obj$est$NR_status & (obj$est$NR_iter <= obj$est$NR_max_iter) & !obj$est$ML_status_fail) {
    
    # Compute weights
    if (obj$verbose > 1) {
      cat("\n    Computing step.")
    }
    weights  <- comp_weights(obj)
    if (sum(weights == "broken") > 0) { 
      obj$est$ML_status_fail <- TRUE
    }
    
    # Compute information matrix
    if (!obj$est$par_flag & !obj$est$ML_status_fail) {
      info_mat <- comp_info_mat(obj, weights)
    } else if (obj$est$par_flag & !obj$est$ML_status_fail) {
      info_mat <- comp_info_mat_par(obj, weights)
    }
    if (obj$est$ML_status_fail == FALSE) { 
      obj$est$info_mat <- info_mat
    }
    
    # Compute the step
    if (!obj$est$ML_status_fail) { 
      step_list <- comp_step(obj, weights, info_mat)
      step <- step_list$step
      if (!is.null(step_list$ML_status_fail)) {
        obj$est$ML_status_fail <- step_list$ML_status_fail
      }
      if (obj$est$adaptive_step_len) {
        obj$est$NR_step_len_ <- 1 / (1 + sum(abs(step)^2))
      }
      if (any(is.nan((obj$est$theta + step * obj$est$NR_step_len_)))) {
        obj$est$ML_status_fail <- TRUE
      }
    }
    
    if (!obj$est$ML_status_fail) {
      if (obj$verbose > 1) {
        val <- sum(abs(step))
        print_1 <- val
        print_1[val < 0.0001] <- "<.0001"
        print_1[val >= 0.0001] <- formatC(val[val >= 0.0001], digits = 4, format = "f")
        cat(paste0("\n      L1-norm of increment of parameter vector: ", print_1))
      }
      obj$est$score_val <- step_list$score_val
      
      # Check if lowest point
      if (obj$est$NR_iter > 1) { 
        if (norm(obj$est$score_val, type = "2") < norm(obj$est$score_min, type = "2")) { 
          obj$est$theta_min <- obj$est$theta
          obj$est$score_min <- obj$est$score_val
        }
      } else { 
        obj$est$theta_min <- obj$est$theta
        obj$est$score_min <- obj$est$score_val
      }
      
      # Iterate the next step
      obj$est$theta <- obj$est$theta + step * obj$est$NR_step_len_
      if (any(is.nan(obj$est$theta))) {
        obj$est$ML_status_fail <- TRUE
      }
      
      if (obj$verbose > 1) {
        cat(paste0("\n      Iteration: ", 
                   obj$est$NR_iter, 
                   ",  ", 
                   names(obj$net$obs_stats), 
                   " = ", 
                   round(obj$est$theta, digits = 4)))
        cat("\n")
      }
      
      # Update status and iterations
      if (!obj$est$ML_status_fail) {
        obj$est$NR_iter <- obj$est$NR_iter + 1
        obj <- check_convergence(obj, step)
      }
      
    } else {
      if (obj$est$step_err < 3 & !obj$est$adaptive_step_len) {
        if (obj$verbose > 1) {
          cat("\n    Optimization did not converge. Decreasing step length and restarting.\n")
        }
        obj$est$step_err <- obj$est$step_err + 1
        obj$est$NR_step_len_ <- obj$est$NR_step_len_ * obj$est$NR_step_len_multiplier 
        obj$est$theta <- obj$est$theta_0
        obj$est$NR_iter <- 1
        obj$est$ML_status_fail <- FALSE
      } else { 
        cat("\n    Optimization failed to converge.")
        cat("\n    Proposed model may be near-degenerate or the MCMLE may be unstable or not exist.") 
        cat("\n    Decreasing the step length (argument 'NR_step_len' in set_options())") 
        cat(" or increasing the MCMC sample-size may help.")
      }
    }
  }
  
  if (obj$est$NR_iter >= obj$est$NR_max_iter) {
    if (obj$verbose > 0) { 
      cat("\n    NOTE: Optimization reached maximum number of allowed iterations.")
      cat(paste0("\n\n      - Minimum theta found:"))
      cat(paste0("\n        ",
                 names(obj$net$obs_stats),
                 " = ",
                 round(obj$est$theta_min, digits = 4)))
      
      cat("\n\n")
      cat("      - L2-norm of score at this theta: ")
      cat(paste(round(sum(abs(obj$est$score_min^2)), digits = 4)))
    }
    obj$est$max_iter_flag <- TRUE
  }
  obj$est$theta_0 <- obj$est$theta
  obj$est$NR_conv_thresh <- sqrt(sum(step^2))
  return(obj)
}





#### comp_weights
comp_weights <- function(obj) {
  
  weights <- list(rep(numeric(obj$sim$num_obs)))
  obj$sim$stats <- as.matrix(obj$sim$stats)
  theta_diff <- obj$est$theta - obj$est$theta_0
  adjust_flag <- FALSE
  if (any(is.na(obj$sim$stats %*% theta_diff))) { 
      return("broken")
    } else {
    norm_const_full <- sum(exp_fun(obj$sim$stats %*% theta_diff))
    if (norm_const_full == Inf) {
      norm_const_full <- .Machine$double.xmax
      adjust_flag <- TRUE
    }
    if (norm_const_full == 0) {
      norm_const_full <- .Machine$double.xmin
      adjust_flag <- TRUE
    }
    weights$weight_full <- exp_fun(obj$sim$stats %*% theta_diff) /
      norm_const_full
    if (adjust_flag) {
      weights$weight_full <- weights$weight_full / sum(weights$weight_full)
    }
    if (obj$net$na_flag) {
      if (any(is.na(obj$sim$cond_stats[[i]] %*% theta_diff))) { 
        return("broken")
      }
      norm_const_cond <- sum(exp_fun(obj$sim$cond_stats[[i]] %*% theta_diff))
      weights$weight_cond[[i]] <- exp_fun(obj$sim$cond_stats[[i]] %*% theta_diff) /
        norm_const_cond
    }

  return(weights)
 }
}




#### comp_info_mat
comp_info_mat <- function(obj, weights) {
  info_mat <- list(NULL)
  theta_grad_val <- t(diag(length(obj$est$theta)))
  term_1_full <- numeric(obj$net$num_terms)
  term_2_full <- numeric(obj$net$num_terms)
    
  if (obj$net$na_flag) {
    term_1_cond <- numeric(obj$net$num_terms)
    term_2_cond <- numeric(obj$net$num_terms)
  }
    
  list_ <- as.list(data.frame(t(obj$sim$stats)))
  outers_ <- lapply(list_, outer_fun)
  term_1_full <- Reduce("+", Map("*", outers_, weights$weight_full))
  term_2_full <-  as.vector(t(obj$sim$stats) %*% weights$weight_full)
    
  J_full <- t(theta_grad_val) %*%
      (term_1_full - outer(term_2_full, term_2_full)) %*%
      theta_grad_val
    
  if (obj$net$na_flag) {
    list_ <- as.list(data.frame(t(obj$sim$cond_stats[[i]])))
    outers_ <- lapply(list_, outer_fun)
    term_1_cond <- Reduce("+", Map("*", outers_, weights$weight_cond[[i]]))
    term_2_cond <-  as.vector(t(obj$sim$cond_stats[[i]]) %*%
                                weights$weight_cond[[i]])
    
    J_cond <- t(theta_grad_val) %*%
      (term_1_cond - outer(term_2_cond, term_2_cond)) %*%
      theta_grad_val
  }
  
  if (obj$net$na_flag) {
    info_mat[[i]] <- list(J_full - J_cond)
  } else {
    info_mat <- list(J_full)
  }
  
  return(info_mat)
}


comp_info_mat_par <- function(obj, weights) {
  
  theta_grad_val <- t(diag(length(obj$est$theta)))
  
  outer_fun <- function(v) {
      base::outer(v, v)
  }
    info_mat_list <- rep(list(NULL), 1)
    term_1_cond <- numeric(obj$net$num_terms)
    term_2_cond <- numeric(obj$net$num_terms)
    list_ <- as.list(data.frame(t(obj$sim$stats)))
    outers_ <- lapply(list_, outer_fun)
    term_1_full <- Reduce("+", Map("*", outers_, weights$weight_full))
    term_2_full <-  as.vector(t(obj$sim$stats) %*% weights$weight_full)
      
      J_full <- t(theta_grad_val) %*%
        (term_1_full - outer(term_2_full, term_2_full)) %*%
        theta_grad_val
      
      if (obj$net$na_flag) {
        list_ <- as.list(data.frame(t(obj$sim$cond_stats)))
        outers_ <- lapply(list_, outer_fun)
        term_1_cond <- Reduce("+", Map("*", outers_, weights$weight_cond))
        term_2_cond <-  as.vector(t(obj$sim$cond_stats) %*%
                                    weights$weight_cond)
        
        J_cond <- t(theta_grad_val) %*%
          (term_1_cond - outer(term_2_cond, term_2_cond)) %*%
          theta_grad_val
      }
      
      if (obj$net$na_flag) {
        info_mat <- J_full - J_cond
      } else {
        info_mat <- J_full
      }
      info_mat_list <- list(info_mat)
    
    return(info_mat_list)
}




### comp_step
comp_step <- function(obj, weights, info_mat) {
  obj$sim$stats <- as.matrix(obj$sim$stats)
  if (is.nan(norm(info_mat[[1]]))) {
    ML_status_fail <- TRUE
    cat("\n    Information matrix is near-singular.")
  } else {
    ML_status_fail <- FALSE
  }
  if (!ML_status_fail) {
    info_mat_inv <- tryCatch({ solve(info_mat[[1]]) },
                             error = function(e) {
                               return("error")
                             })
    
    if (is.character(info_mat_inv)) {
      ML_status_fail <- TRUE
      cat("\n    Information matrix is near-singular.")
    } 
  }## Start working from here
  if (!ML_status_fail) {
    theta_grad_val <- t(diag(length(obj$est$theta)))
    exp_approx <- t(weights$weight_full) %*% obj$sim$stats
    if (obj$net$na_flag) {
      obs_approx <- t(weights$weight_cond) %*% obj$sim$cond_stats
    } else {
      obs_approx <- obj$net$obs_stats_step
    }
    
    if ((norm(theta_grad_val) == Inf) | is.nan(norm(theta_grad_val))) {
      ML_status_fail <- TRUE
      cat("\n    Gradient is not finite. Stopping optimization.")
    }
    
    if (!ML_status_fail) {
      score_val <- t(theta_grad_val) %*% (as.vector(obs_approx) - as.vector(exp_approx))
      step <- info_mat_inv %*% score_val
    }
    
  }
  if (ML_status_fail) {
    step <- NA
    score_val <- NA
  }
  return(list(step = step, score_val = score_val, ML_status_fail = ML_status_fail))
}



### check_convergence
check_convergence <- function(obj, step) {
  
  step_norm <- sqrt(sum(step^2))
  
  if (!(step_norm == Inf)) {
    
    conv_check <- sum(abs(step) > obj$est$adj_NR_tol) == 0
    
    # Check if step is smaller than some tolerance
    if (conv_check) {
      obj$est$NR_status <- TRUE
    }
    
  } else {
    obj$est$step_err <- obj$est$step_err + 1
    if (obj$est$step_err < 3) {
      obj$est$NR_step_len_ <- obj$est$NR_step_len_ * obj$est$NR_step_len_multiplier
      obj$est$NR_iter <- 1
    } else {
      obj$est$ML_status_fail <- TRUE
    }
  }
  return(obj)
}

