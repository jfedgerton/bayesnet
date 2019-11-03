compute_initial_estimate <- function(obj, form_sim, data) {
  
  # Find the initial point for a non-curved ERGM 
  
    if (summary(net ~ edges) == 0) { 
      stop("Network provided is empty and has no edges. Maximum likelihood estimator will not exist.", 
           call. = FALSE) 
    }
    init <- suppressMessages(
      bglmer(form_sim,  data,
            family = binomial,
            #cov.prior = NULL,
            fixef.prior = t())
    )
    check_fe <- apply(coef(init)$Group_ID, 2, sd) == 0
    re_var <- names(check_fe)[check_fe == F]
    theta <- fixef(init)
    theta <- theta[names(theta) != re_var]
    extract_re <- t(coef(init)$Group_ID[names(coef(init)$Group_ID) == re_var])
    row.names(extract_re) <- NULL
    names(extract_re) <- row.names(coef(init)$Group_ID)
    obj$est$theta_wg <- c(extract_re, theta)
    obj$est$theta <- fixef(init)
    names(obj$est$theta )[names(obj$est$theta) == "(Intercept)"] <- "Grand mean"
    obj$est$theta_0 <- obj$est$theta 
    btw.var <- data.frame(VarCorr(init))
    obj$est$btw.var <- btw.var$vcov
    names(obj$est$btw.var) <- "Between group variance"
  obj$est$chat <- init
  return(obj) 
}
