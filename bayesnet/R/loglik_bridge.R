bridge_fun <- function(net, form, theta, offset, burnin, interval, num_bridges, sample_size) { 
  form <- as.formula(paste0("net ~ ", as.character(form)[3]))
  model <- ergm_model(form, net)
  etamap <- model$etamap
  coef_names <- get_coef_names(model, !is.curved(model))
  if (offset == TRUE) { 
    if ("edges" %in% coef_names) {
      edge_loc <- which(coef_names == "edges")
      theta[edge_loc] <- theta[edge_loc] - log(network.size(net))
    }
    if ("mutual" %in% coef_names) {
      mutual_loc <- which(coef_names == "mutual")
      theta[mutual_loc] <- theta[mutual_loc] + log(network.size(net))
    }
  }
  bridge_val <- suppressMessages(
    ergm.bridge.llr(form, 
                    to = theta, 
                    from = rep(0, length(theta)), 
                    llronly = TRUE,
                    control = control.ergm.bridge(MCMC.samplesize = sample_size, 
                                                  MCMC.interval = interval,
                                                  MCMC.burnin = burnin, 
                                                  nsteps = num_bridges)))
  
  return(bridge_val) 
}

lik_fun <- function(form, memb, theta, bridge_num = 10, ncores = 3, form_net, offset = FALSE, 
                    burnin = NULL, interval = NULL, sample_size = NULL) {
  
  
  # Make net_list + compute obs
  network <- net
  obs <- summary(form_net)
  
  terms <- as.character(form)[3]
  # Simulate bridges
  bridges <- bridge_fun(net = network, theta = theta, offset = offset, form = form_net, 
                        num_bridges = bridge_num, burnin = burnin, interval = interval, 
                        sample_size = sample_size)
  
  null_bridge <- null_bridge_function(net)
  
  lik_val <- bridges - null_bridge
  
  return(lik_val)
}

null_bridge_function <- function(net) {
  if (is.directed(net)) {
    return.item <- 2 * choose(network.size(net), 2) * log(2)
    } else { 
      return.item <- choose(network.size(net), 2) * log(2)
    }
  return(return.item)
}

get_coef_names <- function(model_obj, is_canonical) { 
  if(is_canonical) {
    model_obj$coef.names
  } else { 
    unlist(lapply(model_obj$terms,
                  function(term) {
                    find_first_non_null(names(term$params),  term$coef.names)
                  }))
  }
}


find_first_non_null <- function(...) { 
  for (x in list(...)) {
    if (!is.null(x)) {
      break
    }
  }
  x
}
