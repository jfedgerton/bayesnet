meergm <- function(formula.stage.1, formula.stage.2, group.data, net, chains = 4, warmup = 100, iter = 5000, cores = 2, seed = NULL, control = NULL, algorithm = "NUTS")
{
  require('ergm')
  require('dplyr')
  require('fergm')
  require('rstan')
  
  
  if (grepl("gwesp", formula.stage.2) == T){
    warning("Curved spaced parameters may be biased.")
  }
  
  
  group.data <- create_group_data(group.data = group.data, formula.stage.1 = formula.stage.1)
  
  if ("Group1" %in% colnames(group.data))
  {
    ties <- prepare_meergm_data(net = net, group.data = group.data, form = formula.stage.2)
    ties$Group1 <- factor(ties$Group1)
    ties$Group2 <- factor(ties$Group2)
    group.data$Group1   <- factor(group.data$Group1)
    group.data$Group2   <- factor(group.data$Group2)
    group.data$Group_ID <- factor(group.data$Group_ID)
    
    ties$Group1 <- as.character(ties$Group1)
    ties$Group2 <- as.character(ties$Group2)
    group.data$Group1 <- as.character(group.data$Group1)
    group.data$Group2 <- as.character(group.data$Group2)
    group.data$Group_ID <- as.character(group.data$Group_ID)
    
    temp_gr <- group.data
    temp_gr$Group1 <- group.data$Group2
    temp_gr$Group2 <- group.data$Group1
    
    ties_temp <- dplyr::bind_rows(group.data, temp_gr)
    ties_temp <-   dplyr::distinct(ties_temp,
                                   Group1, Group2, Group_ID, .keep_all = T)
    ties_temp <- right_join(ties_temp, ties, by = c("Group1", "Group2"))
    ties_temp <- dplyr::arrange(ties_temp, Group_ID)
    
    if (nrow(ties_temp) == nrow(ties))
    {
      ties = ties_temp
    } else
    {
      cat("Error in group identification.")
    }
    cat("Starting data preparation")
    ## Pull out the first stage variables
    first_stage_variables <- names(group.data)
    first_stage_variables <- first_stage_variables[!(first_stage_variables %in% c("Group1", "Group2", "intercept", "Group_ID"))]
    u <- ties[,names(ties) %in% first_stage_variables]
    
    ## Pull out the second stage variables
    second_stage_variables <- names(ties)
    second_stage_variables <- second_stage_variables[!(second_stage_variables %in% c(c("Y", "Sociality1", "Sociality2"), names(group.data)))]
    x <- ties[,names(ties) %in% second_stage_variables]
    
    for (convert_to_numeric in 1:ncol(x))
    {
      x[,convert_to_numeric] <- as.numeric(as.character(x[,convert_to_numeric]))
    }
    
    for (u_convert_to_numeric in 1:ncol(u))
    {
      u[,u_convert_to_numeric] <- as.numeric(as.character(u[,u_convert_to_numeric]))
    }
    
    g.idx <- as.numeric(as.factor(ties$Group_ID))
    y <- as.numeric(as.character(ties[,"Y"]))
    group.id <- unique(cbind(g.idx, ties$Group_ID))
    stan.dta  <- list(L1 = ncol(u),
                      L2 = ncol(x),
                      GD = nrow(group.data),
                      D = nrow(x),
                      group = g.idx,
                      u = u,
                      x = x,
                      y = y)
    stan.dta$x <- as.matrix(stan.dta$x)
    
    rstan_options(auto_write = TRUE)
    options(mc.cores = cores)
    
    # fit the model using STAN
    scode <- "data {
    int<lower=1> L1;		// number of level 1 predictors
    int<lower=1> L2;		// number of level 2 predictors
    int<lower=1> GD;		// number of node group pairs
    int<lower=0> D;		// number of potential ties
    int<lower=1, upper=GD> group[D]; // numerical indicator of edge group (group_ID)
    row_vector[L1] u[D];		// level 1 predictors
    row_vector[L2] x[D];		// level 2 predictors
    int<lower=0,upper=1> y[D];	// ties (outcome)
  }
    parameters {
    real<lower=0> tau;		// sociality dispersion
    vector[GD] a;       // separate edge group intercepts
    real mu;            // group population intercept
    vector[L1] beta1;		// L1 predictor coefficients
    vector[L2] beta2;		// L2 predictor coefficients
    }
    transformed parameters{
    vector[D] phi; // container for linear predictor
    vector[D] alpha; // container for linear predictor
    for (i in 1:D) {
    alpha[i] = a[group[i]] + dot_product(beta1, u[i]);
    phi[i] = dot_product(beta2, x[i]) + 0.5 * (alpha[i]);
    }
    }
    model {
    mu ~ normal(0,1);
    tau ~ gamma(0.0001, 0.0001);
    a ~ normal(mu,tau);
    beta1 ~ normal(0, 0.1);
    beta2 ~ normal(0, 1);
    
    y ~ bernoulli_logit(phi);
    }
    generated quantities {
    vector[D] predictions;
    for (i in 1:D) {
    predictions[i] = bernoulli_rng(
    inv_logit(
    dot_product(beta2, x[i]) + 0.5 * normal_rng( a[group[i]] +
    dot_product(beta1, u[i]),tau)
    )
    );
    }
    }"
    if(is.null(seed)){
      seed = 12345
      set.seed(seed)
      warning("Note: This function relies on simulation so seed was auto set to 12345.  Consider specifying a seed to set to ensure replicability.")
    }
    stan.fit <- stan(model_code = scode,
                     data=stan.dta, chains=chains, warmup=warmup, iter=iter, algorithm = algorithm)
    
    stan.output = rstan::extract(stan.fit)
    var.names.stage1 <- names(u)
    colnames(stan.output$beta1) <- var.names.stage1
    var.names.stage2 <- names(x)
    colnames(stan.output$beta2) <- var.names.stage2
    return(list(stan.fit = stan.fit, stan.output = stan.output, network.group.data = group.data, stan.dta = stan.dta, group.id = group.id, form.stage.1 = formula.stage.1, form.stage.2 = formula.stage.2))
    } else{
      cat('Error in the group data.')
    }
  
}
