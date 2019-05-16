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
  group.data.for.analysis <- network::as.network(adj.mat, directed = F)
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
