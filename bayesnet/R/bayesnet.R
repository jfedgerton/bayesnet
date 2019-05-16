#' @importFrom dplyr
#' @importFrom ergm
#' @importFrom fergm
#' @importFrom rstan
#'
#' @param create_group_data
#' @param group.data the group data for first stage model
#' @param formula.stage.1 the formula for the first stage model

create_group_data <- function(group.data, formula.stage.1)
{
  require('ergm')
  require('dplyr')
  ## Turn the group data into inter and intra group data
  group.nodes <- nrow(as.matrix(group.data))

  ## Create an adjacency matrix for output
  adj.mat <- matrix(0, nrow = group.nodes, ncol = group.nodes)
  colnames(adj.mat) <- rownames(adj.mat) <- group.data$Group_ID
  group.data.for.analysis <- stat.net::as.network(adj.mat, directed = F)
  group.dyads <- network::network.dyadcount(group.data.for.analysis)
  group.var.names <- colnames(group.data)
  group.var.names <- group.var.names[group.var.names != "Group_ID"]

  for (vert.atts in 1:length(group.var.names))
  {
    group.data.for.analysis %v% group.var.names[vert.atts] <- group.data[,group.var.names[vert.atts]]
  }

  if (grepl("absdiffcat", formula.stage.1) == T |
      grepl("transitive", formula.stage.1) == T |
      grepl("triad", formula.stage.1) == T |
      grepl("atleast", formula.stage.1) == T |
      grepl("atmost", formula.stage.1) == T |
      grepl("degree", formula.stage.1) == T |
      grepl("dgwe", formula.stage.1) == T |
      grepl("edgecov", formula.stage.1) == T |
      grepl("kstar", formula.stage.1) == T |
      grepl("absdiffcat", formula.stage.1) == T)
  {
    cat("The program does not support one or more stage 1 coefficients.")
  } else {
    group_form <- stats::as.formula(paste("group.data.for.analysis ~", formula.stage.1))
    dta.array.group <- ergm::ergmMPLE(group_form, output="array", maxMPLEsamplesize=+Inf,
                                      control=ergm::control.ergm(MPLE.max.dyad.types=group.dyads*10))

    ncoef.group <- length(dta.array.group$predictor[1,2,])
    dta.group <- matrix(0, nrow=group.dyads, ncol=3 + ncoef.group)

    idx.group <- 1
    for (tail.group in 1:(group.nodes-1)) {
      for (head.group in (tail.group+1):group.nodes) {
        dta.group[idx.group,] <- c(dta.array.group$response[tail.group, head.group],
                                   dta.array.group$predictor[tail.group, head.group, ],
                                   tail.group,
                                   head.group)
        idx.group <- idx.group+1
      }
    }

    dta.group <- data.frame(dta.group)
    nm.group <- c("Y", names(dta.array.group$predictor[tail.group, head.group, ]),
                  "Group1", "Group2")
    colnames(dta.group) <- nm.group

    all_groups <- unique(c(dta.group$Group1,dta.group$Group2))
    self_loop_data <- data.frame(matrix(NA, ncol = ncol(dta.group), nrow = length(all_groups)))
    colnames(self_loop_data) <- nm.group
    self_loop_data$Y <- 0
    self_loop_data$Group1 <- all_groups
    self_loop_data$Group2 <- all_groups

    for (self_loop in 2:ncol(self_loop_data))
    {
      if (grepl("nodematch.", colnames(self_loop_data)[self_loop]) == T)
      {
        self_loop_data[,colnames(self_loop_data)[self_loop]] <- 1
      } else if (grepl("absdiff.", colnames(self_loop_data)[self_loop]) == T)
      {
        self_loop_data[,colnames(self_loop_data)[self_loop]] <- 0
      }
    }
    dta.group <- dplyr::bind_rows(dta.group,
                                  self_loop_data)
    dta.group <- dplyr::arrange(dta.group, Group1, Group2)
    dta.group <- dplyr::select(dta.group, -Y)
    dta.group <- dplyr::mutate(dta.group, intercept = 1)
    dta.group$Group_ID <- 1:nrow(dta.group)

    Grp <- to_indicator(dta.group[,c("Group1", "Group2", "Group_ID")], "Group")
    dta.group[, "Group1"] <- Grp[,1]
    dta.group[, "Group2"] <- Grp[,2]
    dta.group[, "Group_ID"] <- Grp[,3]

    return(dta.group)
  }
}

to_indicator <- function(X, text, directed=FALSE)
{
  if (isTRUE(directed)) {
    stop("Directed networks not yet supported.")
  } else {
    for (k in 1:ncol(X)){
      X[,k] <- paste0(text, X[,k])
      factor.levels <- unique(unlist(X))
      X[,k] <- factor(X[,k], levels=factor.levels)
    }
  }
  X
}

prepare_meergm_data <- function(net, group.data, form, verbose=FALSE)
{
  ## Temporary until directed code is added.
  if (network::is.directed(net)) stop("Directed networks not yet supported.")

  if (isTRUE(verbose)) cat("\n## Preparing FERGM dataset...")
  nodes <- nrow(as.matrix(net))
  ndyads <- network::network.dyadcount(net)
  form <- stats::as.formula(paste("net ~", form))

  if (isTRUE(verbose)) cat("\n##   building array...")
  dta.array <- ergm::ergmMPLE(form, output="array", maxMPLEsamplesize=+Inf,
                              control=ergm::control.ergm(MPLE.max.dyad.types=ndyads*10))

  if (isTRUE(verbose)) cat("\n##   building data.frame...")
  ncoef <- length(dta.array$predictor[1,2,])
  dta <- matrix(0, nrow=ndyads, ncol=5+ncoef)
  group <- ergm::get.node.attr(group.data,"Group")

  idx <- 1
  for (tail in 1:(nodes-1)) {
    for (head in (tail+1):nodes) {
      dta[idx,] <- c(dta.array$response[tail, head],
                     dta.array$predictor[tail, head, ],
                     group[tail],
                     group[head],
                     tail,
                     head)
      idx <- idx+1
    }
  }

  dta <- data.frame(dta)
  nm <- c("Y", names(dta.array$predictor[tail, head, ]), "Group1", "Group2",
          "Sociality1", "Sociality2")
  names(dta) <- nm

  if (isTRUE(verbose)) cat("\n##   setting random effects indicators...\n")
  if (network::is.directed(net)) {
    stop("Directed networks not yet supported.")
  } else {
    Soc <- to_indicator(dta[,c("Sociality1", "Sociality2")], "Node")
    Grp <- to_indicator(dta[,c("Group1", "Group2")], "Group")
    dta[, "Sociality1"] <- Soc[,1]
    dta[, "Sociality2"] <- Soc[,2]
    dta[, "Group1"] <- Grp[,1]
    dta[, "Group2"] <- Grp[,2]
  }

  dta
}

meergm <- function(formula.stage.1, formula.stage.2, group.data, net, chains=4, warmup=100, iter=1000, cores = 1, seed = 12345)
{
  require('ergm')
  require('dplyr')
  require('fergm')
  if (exists('group.data') && is.data.frame(get('group.data')))
  {
    if (grepl("gwesp", formula.stage.2) == T){
      warning("Curved spaced parameters may be biased.")
    }

    if(is.directed(net) == T)
    {

    } else
    {




      ties <- prepare_meergm_data(net, formula.stage.2)
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

      g.idx <- as.numeric(as.factor(ties$Group_ID))
      y <- ties[,"Y"]

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
      tau ~ gamma(0.001, 0.001);
      a ~ normal(mu,tau);
      beta1 ~ normal(0, 1);
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
    if(!is.null(seed)){
      set.seed(seed)
      cat("Setting seed at the default value of 12345 for the seed argument.")
  } else {
    warning("Note: This function relies simulation.  Consider specifying a seed to set to ensure replicability.")
  }
      stan.fit <- stan(model_code = scode,
                       data=stan.dta, chains=chains, warmup=warmup, iter=iter)

      stan.output = rstan::extract(stan.fit)
      var.names.stage1 <- names(u)
      colnames(stan.output$beta1) <- var.names.stage1
      var.names.stage2 <- names(x)
      colnames(stan.output$beta2) <- var.names.stage2
      return(list(stan.fit = stan.fit, stan.output = stan.output, stan.dta = stan.dta, form.stage.1 = formula.stage.1, form.stage.2 = formula.stage.2))
  }

  } else {
    cat('Compile group data before running meergm.')
  }

}

