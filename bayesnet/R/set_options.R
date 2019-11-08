set_options <- function(burnin = 100,
                        interval = 1500,
                        sample_size = 2000,
                        NR_tol = 1e-4,
                        NR_max_iter = 50,
                        MCMLE_max_iter = 10,
                        do_parallel = TRUE,
                        number_cores = detectCores(all.tests = FALSE, logical = TRUE) - 1,
                        adaptive_step_len = TRUE,
                        step_len_multiplier = 0.5,
                        step_len = 1,
                        bridge_num = 10,
                        bridge_burnin = 1e+4,
                        bridge_interval = 500,
                        bridge_sample_size = 5000) {
  
  
  sim_param <- list(burnin = burnin,
                    interval = interval,
                    num_obs = sample_size,
                    stats = NULL,
                    cond_stats = NULL,
                    bridge_burnin = bridge_burnin,
                    bridge_interval = bridge_interval,
                    bridge_sample_size = bridge_sample_size)
  
  
  est_param <- list(eta = NULL,
                    eta_0 = NULL,
                    eta_grad = NULL,
                    eta_fun = NULL,
                    score_val = NULL,
                    NR_tol = NR_tol,
                    NR_iter = 1,
                    NR_max_iter = NR_max_iter,
                    NR_status = FALSE,
                    step_err = 0,
                    MCMLE_iter = 1,
                    MCMLE_max_iter = MCMLE_max_iter,
                    MCMLE_status = FALSE,
                    info_mat = NULL,
                    bridge_num = bridge_num,
                    adaptive_step_len = adaptive_step_len,
                    NR_step_len = step_len,
                    NR_step_len_multiplier = step_len_multiplier,
                    NR_conv_thresh = NULL,
                    MCMLE_conv_thresh = NULL,
                    par_flag = do_parallel,
                    par_n_cores = number_cores,
                    ML_status_fail = FALSE)
  
  return(list(sim_param = sim_param,
              est_param = est_param))
}
