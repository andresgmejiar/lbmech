#' Plot the local group and non-group components of a local indicator of dispersion,
#' colored by their inference-based class.
#' 
#' \code{colorLID()} acts as a function converting class names to the hex codes corresponding
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
#' @param key Logical. Should a color key be printed? Default is \code{Key = TRUE}. Set as 
#' \code{FALSE} if you wish to further manipulate the ggplot object
#' @param only.key Logical. If \code{only.key = TRUE}, only a color key will be output and all
#' other parameters will be ignored. 
#' @param arrows Logical. Should the points be displayed as arrows pointing to the appropriate
#' quadrant according to the key? Default is \code{arrows = FALSE}. Ignored if \code{table = TRUE} in 
#' colorLID (Cairo may be needed to export images and PDFs). 
#' @return A ggplot object with two elements---the LID Scatter plot and its scale.
#' @examples 
#' 
#' # Generate dummy observations
#' x <- runif(10, 1, 100)
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
#' inference <- inferLID(lid, w = weights, ntrials = 100, pb = FALSE)
#' 
#' # Plot the inferences
#' scatterLID(lid, inference)
#' @export
scatterLID <- function(lid = NULL, inference = NULL, log.scale = FALSE, key = TRUE, 
                       only.key = FALSE, arrows = FALSE, x.lim = NULL, y.lim = NULL){
  # CRAN Check silencing bit
  J_Gi=J_NGi=Class=x=y=NULL
  
  # lid.cols is an internal variable
  
  
  # Create a color assignment function
  colvect <- lid.cols$Color
  names(colvect) <- lid.cols$Class
  if (!only.key){
    
    # Perform a log transformation if needed
    if (!log.scale){
      xlab <- bquote("In-group Inequality (J "["G"]^.(lid$index) * " )")
      ylab <- bquote("Out-group Inequality (J "["NG"]^.(lid$index) *" )")
      FUN <- function(x) x
    } else if (log.scale){
      xlab <- bquote("In-group Inequality (log"["10"] * "1 + J "["G"]^.(lid$index) * " )")
      ylab <- bquote("Out-group Inequality (log"["10"] *"1 + J "["NG"]^.(lid$index) * " )")
      FUN <- function(x) log(1 + x, 10)
    }
    
    if (arrows) {
      windowsFonts <- grDevices::windowsFonts
      pdfFonts <- grDevices::pdfFonts
      postscriptFonts <- grDevices::postscriptFonts
      extrafont::loadfonts(quiet = TRUE)
    }
    # LID plot, transforming the values if needed
    gglid <- ggplot2::ggplot(cbind(lid$local,inference$local), ggplot2::aes(x = FUN(J_Gi),
                                                                            y = FUN(J_NGi),
                                                                            color = Class)) + 
      {if (!arrows) ggplot2::geom_point()} + 
      {if (arrows) ggplot2::geom_text(label = colorLID(inference$local$Class, 
                                                       arrows = TRUE))} + 
      ggplot2::scale_color_manual(values = colvect) + 
      ggplot2::ggtitle(paste0("Local Indicators of Dispersion\n",names(lid$index)," Scatterplot")) +
      ggplot2::theme_dark() + 
      ggplot2::theme(legend.position = 'none',
                     plot.title = ggplot2::element_text(hjust = 0.5)) + 
      ggplot2::scale_x_continuous(xlab, limits = x.lim) +
      ggplot2::scale_y_continuous(ylab, limits = y.lim)
  } else key <- TRUE
  
  if (key){
    # Scale is just a tile form of the lid.cols variable
    ggcolors <- ggplot2::ggplot(lid.cols, ggplot2::aes(fill = Class, x = x, y = y)) + 
      ggplot2::geom_tile() + 
      ggplot2::theme_bw() + 
      ggplot2::theme(legend.position = 'none',
                     plot.title = ggplot2::element_text(hjust = 0.5, size = 12),
                     axis.text = ggplot2::element_text(angle=45,size = 8),
                     axis.title = ggplot2::element_text(size = 10)) + 
      ggplot2::scale_fill_manual(values = colvect) +
      ggplot2::scale_x_continuous("In-group Inequality",
                                  expand=c(0,0),
                                  breaks = c(-1,0,1), 
                                  labels = c('-1' = 'Low', '0' = 'Avg.', '1' = 'High')) +
      ggplot2::scale_y_continuous("Out-group Inequality", 
                                  expand=c(0,0),
                                  breaks = c(-1,0,1), 
                                  labels = c('-1' = 'Low', '0' = 'Avg.', '1' = 'High'),
                                  position = 'right') + 
      ggplot2::ggtitle("Significance") + 
      ggplot2::coord_equal()
    
    # Use cowplot to arrange the grid
    if (!only.key) {
      return(cowplot::plot_grid(gglid,ggcolors,ncol=2,rel_widths = c(0.8,0.2)))
    } else return(ggcolors)
  } else return(gglid)  
}

#' @rdname scatterLID
#' @param x A character string or vector containing a LID significance class. 
#' Ignored if \code{table = TRUE}.
#' @param table Logical. Should the function convert character strings of classes
#' to hex codes of colors (\code{table = FALSE}, the default) or arrows (when \code{arrows = TRUE})
#' or should it return the conversion table itself?
#' @param arrows Logical. Should the points be displayed as arrows pointing to the appropriate
#' quadrant according to the key? Default is \code{arrows = FALSE}. Ignored if \code{table = TRUE} in 
#' colorLID.
#' @export
colorLID <- function(x = NULL, table = FALSE, arrows = FALSE){
  # This bit to silence CRAN warnings
  Class=Color=Arrow=NULL
  if (table){
    return(lid.cols[,.(Class,Color,Arrow)])
  } else if (!table){
    if (!arrows){
      return(lid.cols$Color[match(x, lid.cols$Class)])
    } else{
      return(lid.cols$Arrow[match(x,lid.cols$Class)])
    }
  }
}
