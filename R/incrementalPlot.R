#' Plot the difference in the non-group component of inequality at increasing
#' bandwidth sizes. 
#' 
#' @title Incremental Local Indicators of Dispersion plot
#' @param inc The list output of the \code{\link[lbmech]{incrementalLID}} function
#' This does not work when both \code{bw} and \code{def.neigh} were permuted. 
#' @return A ggplot object
#' @examples 
#' 
#' # Generate dummy observations
#' x <- runif(10, 1, 100)
#' 
#' # Get distance matrix
#' dists <- dist(x)
#' 
#' # Bandwidth sizes from 3 to 5
#' bws <- 3:5
#' 
#' inc <- incrementalLID(x, dist = dists, bws = bws, index = 'gini', type = 'local',
#'                       weighting = 'distance', FUN = function(x) 1/x^2, minval = 1)
#' 
#' # Plot the results
#' incrementalPlot(inc)
#' @export
incrementalPlot<- function(inc){
  plotme=NULL
  # This bit to silence CRAN check warnings
  bw=delta_J_NG=NG_Class=NULL
  
  plotVar <- unlist(inc$stats[, lapply(.SD, function(x) max(x) != min(x))][,1:2])
  plotVar <- names(plotVar)[plotVar]
  
  inc$stats$plotme <- inc$stats[,.SD,.SDcols = plotVar]
  inc$plotme <- inc[[2]]
  
  # Line plot with vertical lines at each extrema
  gg <- ggplot2::ggplot(inc$stats, ggplot2::aes(x=plotme,y=delta_J_NG, color = NG_Class)) + 
    ggplot2::geom_vline(data = data.table(bw = inc$plotme), 
                        ggplot2::aes(xintercept = inc$plotme), 
                        linetype = 'dashed', 
                        color = 'yellow', 
                        linewidth = 1) +
    ggplot2::geom_point() + 
    ggplot2::scale_color_manual(values = c(Significant = "red",
                                           `Not Significant` = 'blue')) + 
    ggplot2::geom_line(ggplot2::aes(color = 'black')) +
    ggplot2::scale_y_continuous(expression(paste(delta, J[NG]))) + 
    ggplot2::scale_x_continuous("Bandwidth") + 
    ggplot2::ggtitle(paste0("Out-group ",
                            stringr::str_extract(names(inc$index),
                                                 "[A-Za-z\\-]+"),
                            " Inequality\n Difference at Incremental Bandwidths")) + 
    ggplot2::labs(subtitle = paste0("Extrema at Bandwidth(s) = ",
                                    paste(inc$plotme,collapse = ', '))) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                   plot.subtitle = ggplot2::element_text(hjust = 0.5),
                   legend.title = ggplot2::element_blank())
  return(gg)
}


