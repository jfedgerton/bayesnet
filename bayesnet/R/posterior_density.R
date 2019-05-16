#' @importFrom ggplot
#' @importFrom ergm
#' @importFrom RColorBrewer
#'
#' @param meergm.fit model output
#' @param custom_var_names overwrite automatic variable names

model_check <- function(meergm.fit = NULL, custom_var_names = NULL){
  require(ggplot2)
  require(RColorBrewer)
  stage.1 <- data.frame(meergm.fit$stan.output$beta1)
  stage.2 <- data.frame(meergm.fit$stan.output$beta2)
  group.estimates <- data.frame(meergm.fit$stan.output$a)
  group.estimates.tau <- data.frame(meergm.fit$stan.output$tau)
  colnames(group.estimates) <- paste0("Group", 1:ncol(group.estimates))
  colnames(group.estimates.tau) <- "Tau"
  all_var_to_plot <- c(colnames(stage.1), colnames(stage.2), colnames(group.estimates), colnames(group.estimates.tau))

  data_to_plot <- dplyr::bind_cols(stage.1, stage.2,
                                   group.estimates, group.estimates.tau)

  for (coef_loop in 1:length(all_var_to_plot))
  {
    temp_data <- data.frame(data_to_plot[,names(data_to_plot) %in% all_var_to_plot[coef_loop]])
    colnames(temp_data) <- "Var_Plot"
    iteration_check <- nrow(temp_data)/ncol(meergm.fit$stan.fit)
    temp_data$iteration <- rep(1:iteration_check, ncol(meergm.fit$stan.fit))
    temp_data$chains <- as.character(sort(rep(1:ncol(meergm.fit$stan.fit), iteration_check)))
    dens_plot <- ggplot(temp_data) +
      geom_density(aes(Var_Plot), fill = "lightblue", adjust = 3) +
      labs(x = paste0("n = ", nrow(data_to_plot)),
           title = paste0("Posterior Distribution ", all_var_to_plot[coef_loop])) +
      theme_bw()

    trace_plot <- ggplot(temp_data) +
      geom_line(aes(x = iteration, y = Var_Plot, col = chains)) +
      labs(x = paste0("iterations = ", iteration_check), y = all_var_to_plot[coef_loop],
           title = paste0("Trace Plot ", all_var_to_plot[coef_loop])) +
      theme_bw() +
      scale_color_brewer(palette = "Set1") +
      theme(legend.position="none")
    gridExtra::grid.arrange(dens_plot, trace_plot, ncol=2)
  }
}
