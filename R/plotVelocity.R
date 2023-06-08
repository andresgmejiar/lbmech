#' Plot the log-probability of the observed velocity data points versus the 
#' regressed nonlinear quantile regression
#'
#' @title Plot GPS Velocities
#' @param velocity The list output to \code{\link[lbmech]{getVelocity}}.
#' @param v_lim The maximum velocities to plot (y-axis limit). Default is 3 m/s.
#' @param v_min The minimum velocities to plot (y-axis limit). Default is 0.
#' @param slope_lim The maximum slopes to plot (x-axis limits). Default is 1.
#' @param bins Into how many bins are the axes divided? 
#' @param x.bins Into how many bins are the axes divided? 
#' @param y.bins Into how many bins are the axes divided? 
#' @importFrom data.table data.table
#' @return A ggplot object
#' @examples
#' # Note that the output results should be senseless since they
#' # are computed on random data
#' 
#' # If the data contains an 'elevation' or 'z' column
#' data <- data.table(x = runif(10000,10000,20000),
#'                    y = runif(10000,30000,40000),
#'                    elevation = runif(10000,0,200),
#'                    dt = 120,
#'                    ID = rep(1:10,each=1000))
#' velocity <- getVelocity(data = data, z = 'elevation')
#' plotVelocity(velocity)
#' @export
plotVelocity <- function(velocity, v_lim = 3, v_min = 0, slope_lim = 1,
                         bins = 100, x.bins = bins, y.bins = bins){
  # This bit to silence CRAN check warnings
  dl_dt=dz_dl=.N=N_x=N_xy=x_plot=N_plot=y_plot=NULL
  
  # Extract the relevant variables
  v_max <- velocity$vmax
  s <- velocity$s
  k <- velocity$k
  
  # Crop the data.table to the values that will be plotted
  d1 <- velocity$data[dl_dt < v_lim & dl_dt > v_min & abs(dz_dl) < slope_lim]
  
  # Bin the values to the chosen resolution
  d1$y_plot <- round(d1$dl_dt*y.bins) / y.bins
  d1$x_plot <- round(d1$dz_dl*x.bins) / x.bins
  
  # Calculate the densities
  d1_plot <- d1[,.(N_xy = .N) ,by=c("x_plot","y_plot")
  ][,N_x := sum(N_xy), by=x_plot][
    ,N_plot := N_xy/N_x ][]
  
  # Get the points for the regressed function
  regressed <- data.table(dz_dl = seq(-slope_lim,slope_lim,length.out = 1001))
  regressed[, dl_dt := v_max * exp(-k * abs(dz_dl - s))]
  
  # Value is the log of the number of observations per x
  ggplot2::ggplot(d1_plot, ggplot2::aes(x = x_plot, y = y_plot, fill = log(N_plot,10))) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_viridis_c(bquote("Log"[10]*" Prob")) +
    ggplot2::geom_line(data = regressed, 
                       mapping = ggplot2::aes(x = dz_dl, 
                                              y = dl_dt,
                                              color = "Regressed\nFunction", 
                                              fill = NULL),
              linewidth = 1) + 
    ggplot2::scale_color_discrete("") + 
    ggplot2::ggtitle("Slope vs. Planimetric Speed\nHeatmap and Regressed Function",
          subtitle=paste0("(N = ",length(unique(velocity$data$ID)),
                   " tracks, ",nrow(velocity$data)," points)")) + 
    ggplot2::xlab("Slope (dz/dl, dimensionless)") + 
    ggplot2::ylab("Planimetric Speed (dl/dt, in m/s)") + 
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                   plot.subtitle = ggplot2::element_text(hjust = 0.5))
}
