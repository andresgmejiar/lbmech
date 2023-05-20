#' Plot the local group and non-group components of a local indicator of dispersion,
#' colored by their inference-based class.
#' 
#' colorLID() acts as a function converting class names to the hex codes corresponding
#' to the colors used by scatterLID when \code{table = FALSE} (the default), and
#' returns the color table itself when \code{table = FALSE}.
#' 
#' @title Local Indicators of Dispersion Scatterplot
#' @name scatterLID
#' @aliases colorLID
#' @param lid The list output of the \code{\link[lbmech]{LID}} function.
#' @param inference The list output of the \code{\link[lbmech]{inferLID}} function.
#' @param log.scale Logical. Should the axes be log-transformed? Default is \code{FALSE}.
#' If \code{TRUE}, log transformation is \code{log(1+x,10)}.
#' @param x.lim One of \code{NULL} to determine the x-range automatically (the default), 
#' a numeric vector of length two providing the x boundaries, or a function that accepts
#' the automatic boundaries and returns new limits (see \code{\link[ggplot2]{scale_x_continuous}}).
#' @param y.lim One of \code{NULL} to determine the y-range automatically (the default), 
#' a numeric vector of length two providing the y boundaries, or a function that accepts
#' the automatic boundaries and returns new limits (see \code{\link[ggplot2]{scale_y_continuous}}).
#' @return A ggplot object with two elements---the LID Scatter plot and its scale.
#' @examples 
#' 
#' # Generate dummy observations
#' x <- runif(10, 0, 100)
#' 
#' # Get distance matrix
#' dists <- dist(x)
#' 
#' # Get fuzzy weights considering 5 nearest neighbors based on 
#' # inverse square distance
#' weights <- makeWeights(dists, bw = 5, 
#'                        mode = 'adaptive', weighting = 'distance',
#'                        FUN = function(x) 1/x^2, minval = 0.1,
#'                        row.stand = 'fuzzy')
#'                        
#' # Obtain the 'local gini' value
#' lid <- LID(x, w = weights, index = 'gini', type = 'local')
#' 
#' # Infer whether values are significant relative to the spatial distribution
#' # of the neighbots
#' inference <- inferLID(lid, w = weights, ntrials = 100)
#' 
#' # Plot the inferences
#' scatterLID(lid, inference)
#' @export
scatterLID <- function(lid, inference, log.scale = FALSE, x.lim = NULL, y.lim = NULL){
# lid.cols is an internal variable

# Create a color assignment function
colvect <- lid.cols$Color
names(colvect) <- lid.cols$Class

# Perform a log transformation if needed
if (!log.scale){
xlab <- bquote("Inequality Within Groups (J "["G"]^.(lid$index) * " )")
ylab <- bquote("Inequality Across Groups (J "["NG"]^.(lid$index) *" )")
FUN <- function(x) x
} else if (log.scale){
  xlab <- bquote("Inequality Within Groups (log"["10"] * "1 + J "["G"]^.(lid$index) * " )")
  ylab <- bquote("Inequality Across Groups (log"["10"] *"1 + J "["NG"]^.(lid$index) * " )")
FUN <- function(x) log(1 + x, 10)
}

# LID plot, transforming the values if needed
gglid <- ggplot2::ggplot(cbind(lid$local,inference$local), ggplot2::aes(x = FUN(G_Gi),
                                                      y = FUN(G_NGi),
                                                      color = Class)) + 
  ggplot2::geom_point() + 
  ggplot2::scale_color_manual(values = colvect) + 
  ggplot2::ggtitle(paste0("Local Indicators of Dispersion\n",names(lid$index)," Scatterplot")) +
  ggplot2::theme(legend.position = 'none',
        plot.title = ggplot2::element_text(hjust = 0.5)) + 
  ggplot2::scale_x_continuous(xlab, limits = x.lim) +
  ggplot2::scale_y_continuous(ylab, limits = y.lim)

# Scale is just a tile form of the lid.cols variable
ggcolors <- ggplot2::ggplot(lid.cols, ggplot2::aes(fill = Class, x = x, y = y)) + 
  geom_tile() + 
  ggplot2::theme_bw() + 
  ggplot2::theme(legend.position = 'none',
        plot.title = ggplot2::element_text(hjust = 0.5, size = 12),
        axis.text = ggplot2::element_text(angle=45,size = 8),
        axis.title = ggplot2::element_text(size = 10)) + 
  ggplot2::scale_fill_manual(values = colvect) +
  ggplot2::scale_x_continuous("Local Inequality",
                     expand=c(0,0),
                     breaks = c(-1,0,1), 
                     labels = c('-1' = 'Low', '0' = 'Avg.', '1' = 'High')) +
  ggplot2::scale_y_continuous("Global Inequality", 
                     expand=c(0,0),
                     breaks = c(-1,0,1), 
                     labels = c('-1' = 'Low', '0' = 'Avg.', '1' = 'High'),
                     position = 'right') + 
  ggplot2::ggtitle("Significance") + 
  ggplot2::coord_equal()

# Use cowplot to arrange the grid
return(cowplot::plot_grid(gglid,ggcolors,ncol=2,rel_widths = c(0.8,0.2)))
}

#' @rdname scatterLID
#' @param x A character string or vector containing a LID significance class. 
#' Ignored if \code{table = TRUE}.
#' @param table Logical. Should the function convert character strings of classes
#' to hex codes of colors (\code{table = FALSE}, the default),
#' or should it return the conversion table itself?
#' @export
colorLID <- function(x = NULL, table = FALSE){
  if (table){
    return(lid.cols[,.(Class,Color)])
  } else if (!table){
    return(
      leaflet::colorFactor(lid.cols$Color,levels = lid.cols$Class,ordered=FALSE)(x))
  }
}


