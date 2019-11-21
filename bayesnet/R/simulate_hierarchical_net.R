#rm(list=ls())
library('ergm')
library('boot')
library('lme4')
library('sna')
## Number of simulations
formula.stage.1 = "absdiff('SES')"
## Initial function parameters
simulate_network_hierarchy <- function(Nodes = 30,  
                                       Groups = 3, 
                                       directed = F, 
                                       b_sex_match = 1, 
                                       b_absdiff_SES = -1,
                                       wth_group = 0.1,
                                       btw_group = 1,
                                       in_group = 1,
                                       alpha = -3,
                                       endog_var = "triangle",
                                       b_endog_var = 1)
{
  require('tidyverse')
  require('ergm')
  require('boot')
  require('MASS')
  require('tidygraph')
 
  n  <- Nodes
  k <- rep(1:Groups, n/Groups)
  position <- runif(n)
  k <- k[order(position)]
  sex <- rep(0:1, n/2)
  position <- runif(n)
  sex <- sex[order(position)]
  ses <- sample(Groups)
  
  original.data <- data.frame(Nodes  = 1:n, "Group_ID" = k, "sex" = sex)
  group.data <- data.frame("Group_ID" = 1:Groups, "SES" = ses)
  
  original.data <- left_join(original.data, 
                             group.data, 
                             by = "Group_ID")
  ## Find all dyadic comparisons 
  simulate_nodes <- expand.grid(1:Nodes, 1:Nodes)
  
  ## Remove self loops
  simulate_nodes <- simulate_nodes[simulate_nodes$Var1 != simulate_nodes$Var2,]
  
  ## If the network is undirected then get rid of self loops 
  if (directed == F)
  {
    simulate_nodes <- simulate_nodes[simulate_nodes$Var1 < simulate_nodes$Var2,]
  }
  
  ## Merge in additional data
  colnames(simulate_nodes) <- c("from", "to")
  simulate_nodes <- arrange(simulate_nodes, 
                            from, to)
  
  ## Merge in additional data
  simulate_nodes.merge.node.1 <- original.data
  
  ## What the the senders' data
  colnames(simulate_nodes.merge.node.1) <- paste0("from_", colnames(simulate_nodes.merge.node.1))
  colnames(simulate_nodes.merge.node.1)[1] <- "from"
  
  ## What the the recievers' data
  simulate_nodes.merge.node.2 <- original.data
  colnames(simulate_nodes.merge.node.2) <- paste0("to_", colnames(simulate_nodes.merge.node.2))
  colnames(simulate_nodes.merge.node.2)[1] <- "to"
  
  ## Merge in sender data
  simulate_nodes <- left_join(simulate_nodes,
                           simulate_nodes.merge.node.1, 
                           by = "from")
  
  ## Merge in reciever data
  simulate_nodes <- left_join(simulate_nodes,
                           simulate_nodes.merge.node.2, 
                           by = "to")
  
  ## Keep this data for later export
  original.data.group <- group.data
  
  colnames(original.data.group)[1] <- "Group_ID"
  
  ## Turn the group data into dyadic data
  group.net.data <- create_group_data(original.data.group, formula.stage.1 = formula.stage.1)
  
  ## Add in a random group effect
  group.net.data$group.effect <- rnorm(n = nrow(group.net.data), mean = 0, sd = btw_group)
  
  ## Create unique ID for distinct observations
  group.net.data.temp <- group.net.data
  group.net.data.temp$Group1 <- group.net.data$Group2
  group.net.data.temp$Group2 <- group.net.data$Group1
  
  group.net.data <- group.net.data.temp %>%
    bind_rows(group.net.data) %>%
    distinct(Group1, Group2, Group_ID, .keep_all = T)
  
  group.net.data %>%
    mutate(from_Group_ID = Group1,
           to_Group_ID = Group2) %>%
    dplyr::select(absdiff.SES, group.effect, from_Group_ID, to_Group_ID, Group_ID) %>%
    inner_join(simulate_nodes, 
               by = c("from_Group_ID", "to_Group_ID")) -> simulate_nodes
  
  simulate_nodes$Sex_Match <- 0
  simulate_nodes$Sex_Match[simulate_nodes$to_Sex == simulate_nodes$from_Sex] <- 1
  simulate_nodes$Group_Match <- 0
  simulate_nodes$to_Group_ID[simulate_nodes$from_Group_ID == simulate_nodes$to_Group_ID] <- 1
  
  simulate_nodes$indv_eff <- rnorm(nrow(simulate_nodes), mean = 0, wth_group)
  
  simulate_nodes$phi = 
    with(simulate_nodes, alpha + b_absdiff_SES *  absdiff.SES + group.effect + 
         b_sex_match * Sex_Match + in_group * to_Group_ID + indv_eff) 
  simulate_nodes$tie <- rbinom(nrow(simulate_nodes),1,inv.logit(simulate_nodes$phi))
  
  node_names <- data.frame(names = 1:n)
  
  ## turn into edgelist 
  node_ties <- subset(simulate_nodes, tie == 1)
  
  ## subset data
  node_ties <- node_ties[,colnames(node_ties) %in% c("from", "to")]
  
  simulated_network <- tbl_graph(nodes = node_names, edges = node_ties, directed = directed) 
  
  vertex_df <- original.data
  colnames(vertex_df)[colnames(vertex_df) == "Group_ID"] <- "Group"
  colnames(vertex_df)[colnames(vertex_df) == "sex"] <- "Sex"
  colnames(vertex_df)[colnames(vertex_df) == "Nodes"] <- "names"
  suppressMessages(simulated_network %>%
                     activate(nodes) %>%
                     inner_join(vertex_df) %>%
                     intergraph::asNetwork(.) -> simulated_network)
  
  #simulate_nodes$Sociality1 <- paste0("Node", simulate_nodes$from)
  #simulate_nodes$Sociality2 <- paste0("Node", simulate_nodes$to)
  sim_names <- names(simulate_nodes)
  ## Add in the endogenous variable
  net_obj_dyadic <- prepare_meergm_data(net = simulated_network, 
                                                    group.data = original.data.group, 
                                                    form = paste0("nodematch('Sex')+nodematch('Group')+", endog_var))
    
    #net_obj_dyadic$Sociality1 <- as.character(net_obj_dyadic$Sociality1)
    #net_obj_dyadic$Sociality2 <- as.character(net_obj_dyadic$Sociality2)
    #simulate_nodes$Sociality1 <- paste0("Node", simulate_nodes$from)
    #simulate_nodes$Sociality2 <- paste0("Node", simulate_nodes$to)
    simulate_nodes <- simulate_nodes %>%
      mutate(Sociality1 = from,
             Sociality2 = to) %>%
      left_join(net_obj_dyadic, by = c("Sociality1", "Sociality2"))
    
    ## Estimate the model again with endogenous attribute
    simulate_nodes$phi = 
      with(simulate_nodes, alpha + b_absdiff_SES *  absdiff.SES + group.effect + 
           b_sex_match * Sex_Match + b_sex_match * nodematch.Group + b_endog_var * simulate_nodes[,colnames(simulate_nodes) %in% gsub("[[:punct:]]", "", endog_var)] + indv_eff ) 
    simulate_nodes$tie <- rbinom(nrow(simulate_nodes),1,inv.logit(simulate_nodes$phi))
    #test_output <- length(list.files("meergm/New simulation output/Test freq/")) + 1
    #save(simulate_nodes, file = paste0("meergm/New simulation output/Test freq/sim_test_", test_output, ".Rda"))
    
    node_names <- data.frame(names = 1:n)
    
    ## turn into edgelist 
    node_ties <- subset(simulate_nodes, tie == 1)
    
    ## subset data
    node_ties <- node_ties[,colnames(node_ties) %in% c("Sociality1", "Sociality2")]
    node_ties <- node_ties[order(node_ties$Sociality1, node_ties$Sociality2),]
    
    colnames(node_ties) <- c("from", "to")
    simulated_network <- tbl_graph(nodes = node_names, edges = node_ties, directed = directed) 
    
    
    suppressMessages(simulated_network %>%
                       activate(nodes) %>%
                       inner_join(vertex_df) %>%
                       intergraph::asNetwork(.) -> simulated_network)
    
    group_data <- create_group_data(group.data = original.data.group, formula.stage.1 = form.stage.1)
    group_data2 <- group_data
    group_data2$Group1 <- group_data$Group2
    group_data2$Group2 <- group_data$Group1
    group_data <- rbind(group_data, 
                        group_data2)
    group_data <- distinct(group_data, Group1, Group2, Group_ID, .keep_all = T)
    simulate_nodes <- prepare_meergm_data(net = simulated_network,
                                          group.data = original.data.group,
                                          form = form.stage.2)
    
    simulate_nodes <- merge(group_data,
                               simulate_nodes, 
                               by = c("Group1", "Group2"))
    
    simulate_nodes <- arrange(simulate_nodes, 
                                 Sociality1, Sociality2)
    
    endo_var_temp <- gsub("\\(", "", endog_var)
    endo_var_temp <- gsub("\\)", "", endo_var_temp)
    formula_hierarhcical <- as.formula(paste0("Y ~ nodematch.Sex + nodematch.Group + absdiff.SES + ", endo_var_temp, "+ (1  | Group_ID)")) 
    
    #merge_id <- group.net.data %>%
    #  dplyr::select(Group1, Group2, Group_ID)
    #simulate_nodes %>%
    #  dplyr::select(-Group_ID) %>%
    #  left_join(merge_id, by = c("Group1", "Group2")) -> simulate_nodes
    
    
    init <- suppressMessages(
      bglmer(formula_hierarhcical,  simulate_nodes,
             family = binomial,
             cov.prior = wishart,
             fixef.prior = normal(sd = c(10, 10)))
    )
    theta <- apply(coef(init)$Group_ID, 2, mean)
    theta <- theta 
    ## Create a new network with the edogenous variables in it
    
    
    
    
  
  
  plot(simulated_network,  vertex.col = simulated_network%v%"Group", vertex.cex = 2)
  return.items <- list(net = simulated_network, group.data = original.data.group, theta = theta)
  return(return.items)
}

seed.sim <- 5000
return.list <- list()
for (i in 1:seed.sim)
{
  set.seed(i)
  return.list[[i]] <- simulate_network_hierarchy()
}

check_triangle <- list()
for (i in 1:length(return.list))
{
  check_triangle[[i]] <- return.list[[i]]$theta
}

all_models <- do.call(bind_rows, check_triangle)
all_models$seed <- 1:nrow(all_models)
all_models$`(Intercept)` <- all_models$`(Intercept)` - -3
all_models$nodematch.Group <- all_models$nodematch.Group - 1
all_models$absdiff.SES <- all_models$absdiff.SES - -1
all_models$nodematch.Sex <- all_models$nodematch.Sex - 1
all_models$triangle <- all_models$triangle - 1

all_models$min <- rowSums(abs(all_models[,1:(ncol(all_models) - 1)]))
all_models <- arrange(all_models, min)

set.seed(5)
net.sim <- simulate_network_hierarchy()
## First simulation coefficients
coefficients_1 <- list(Nodes = 100,  
                       Groups = 5, 
                       directed = F, 
                       b_diff_ses = -1,
                       b_diff_gpa = -0.20,
                       b_sex_match = 0.20, 
                       b_match_rural = 0.5, 
                       sch_sig = 0.1,
                       within_group_sig = 1.5,
                       alpha = -3.5,
                       endog_var = "triangle",
                       b_endog_var = 0.3)

simulate_triangle_75_nodes <- list()
simulate_triangle_100_nodes <- list()
simulate_triangle_125_nodes <- list()
simulate_triangle_150_nodes <- list()
for (i in 1:simulations)
{
  try({simulate_triangle_75_nodes[[i]]   <- simulate_network_hierarchy(Nodes = 75,  directed = F, endog_var = "triangle", alpha = -3.5)})  
  try({simulate_triangle_100_nodes[[i]]  <- simulate_network_hierarchy(Nodes = 100, directed = F, endog_var = "triangle", alpha = -3.5)})  
  try({simulate_triangle_125_nodes[[i]]  <- simulate_network_hierarchy(Nodes = 125, directed = F, endog_var = "triangle", alpha = -3.5)})  
  try({simulate_triangle_150_nodes[[i]]  <- simulate_network_hierarchy(Nodes = 150, directed = F, endog_var = "triangle", alpha = -3.5)})  
}


save(simulate_triangle_150_nodes, file = "meergm/New simulation output/Test data/simulate_triangle_150_nodes.RData")
save(simulate_triangle_125_nodes, file = "meergm/New simulation output/Test data/simulate_triangle_125_nodes.RData")
save(simulate_triangle_100_nodes, file = "meergm/New simulation output/Test data/simulate_triangle_100_nodes.RData")
save(simulate_triangle_75_nodes, file = "meergm/New simulation output/Test data/simulate_triangle_75_nodes.RData")
save(coefficients_1, file = "meergm/New simulation output/Test data/coefficients_1.Rda")

set.seed(54980)
simulate_degree7_75_nodes <- list()
simulate_degree7_100_nodes <- list()
simulate_degree7_125_nodes <- list()
simulate_degree7_150_nodes <- list()
for (i in 1:simulations)
{
  try({simulate_degree7_75_nodes[[i]]   <- simulate_network_hierarchy(Nodes = 75, directed = F, endog_var = "degree(7)", alpha = -3.5)})
  try({simulate_degree7_125_nodes[[i]]  <- simulate_network_hierarchy(Nodes = 125, directed = F, endog_var = "degree(7)", alpha = -3.5)})  
  try({simulate_degree7_100_nodes[[i]]  <- simulate_network_hierarchy(Nodes = 100, directed = F, endog_var = "degree(7)", alpha = -3.5)})  
  try({simulate_degree7_150_nodes[[i]]  <- simulate_network_hierarchy(Nodes = 150, directed = F, endog_var = "degree(7)", alpha = -3.5)})  
}

save(simulate_degree7_75_nodes, file = "meergm/New simulation output/Test data/simulate_degree7_75_nodes.RData")
save(simulate_degree7_125_nodes, file = "meergm/New simulation output/Test data/simulate_degree7_125_nodes.RData")
save(simulate_degree7_100_nodes, file = "meergm/New simulation output/Test data/simulate_degree7_100_nodes.RData")
save(simulate_degree7_150_nodes, file = "meergm/New simulation output/Test data/simulate_degree7_150_nodes.RData")

set.seed(3769)
simulate_mutual_75_nodes <- list()
simulate_mutual_100_nodes <- list()
simulate_mutual_125_nodes <- list()
simulate_mutual_150_nodes <- list()
for (i in 1:simulations)
{
  try({simulate_mutual_75_nodes[[i]]   <- simulate_network_hierarchy(Nodes = 75,  directed = T, endog_var = "mutual", alpha = -3.5)})  
  try({simulate_mutual_100_nodes[[i]]  <- simulate_network_hierarchy(Nodes = 100, directed = T, endog_var = "mutual", alpha = -3.5)})  
  try({simulate_mutual_125_nodes[[i]]  <- simulate_network_hierarchy(Nodes = 125, directed = T, endog_var = "mutual", alpha = -3.5)})  
  try({simulate_mutual_150_nodes[[i]]  <- simulate_network_hierarchy(Nodes = 150, directed = T, endog_var = "mutual", alpha = -3.5)})  
}


save(simulate_mutual_150_nodes, file = "meergm/New simulation output/Test data/simulate_mutual_150_nodes.RData")
save(simulate_mutual_125_nodes, file = "meergm/New simulation output/Test data/simulate_mutual_125_nodes.RData")
save(simulate_mutual_100_nodes, file = "meergm/New simulation output/Test data/simulate_mutual_100_nodes.RData")
save(simulate_mutual_75_nodes, file = "meergm/New simulation output/Test data/simulate_mutual_75_nodes.RData")
