

#' Define a function that generates an approximation function for a given 
#' subset of data. This is designed to be used to predict values in one data.table
#' using values in another via the j slot. See examples.
#'
#' @title Create an approximate mapping function
#' @param dt A data.table containing the observed data
#' @param id.val The value to look up in \code{dt[[id.col]]}. When used in a data.table \code{target}'s
#' \code{j} slot, this would be the name of the column shared between \code{dt} and \code{target} without
#' quotes.
#' @param x A character string representing the name of the independent variable 
#' column in \code{dt}. Default is \code{x = 'x'}.
#' @param y A character string representing the name of the dependent variable 
#' column in \code{dt}. Default is \code{y = 'y'}.
#' @param id.col A character string representing the name of the group column in
#' \code{dt}. Default is \code{id.col = 'ID'}
#' @param rand.val A character vector with values present in \code{dt[[id.col]]}
#' for which to return a random value. Used for psuedo-observations where the 
#' outcome is known regardless of the value of the predicted variable.
#' @param rand.range  A vector of two numbers repretenting t he range of
#'  potential values for \code{rand.range}. Default is \code{rand.range = c(0,1)}.
#' @param FUN The approximation function. Default is \code{\link[stats]{approxfun}},
#' but any function that accepts a two-column data.frame is acceptable. 
#' @param ... Additional parameters to pass to \code{FUN}.
#' @importFrom data.table copy
#' @importFrom data.table :=
#' @examples 
#' # Generate a data.table with different observations for different categories
#' data <- data.table(ID = rep(c("A","B","C"), each = 10), 
#'                    x  = runif(30),
#'                    y  = runif(30, min = 10, max = 20))
#' 
#' # Create a target data.table
#' target <- data.table(ID = rep(c("A","B","C"), each = 101),
#'                      x = rep(seq(0,1,length.out=101), each = 3))
#' 
#' target[, y := approxMap(data, id.val = ID, rule = 2, na.rm=TRUE)(x),by='ID']
#' @export
approxMap <- function(dt, id.val, x = 'x', y = 'y', id.col = 'ID',
                      rand.val = NULL, rand.range = c(0,1), 
                      FUN = stats::approxfun, ...){
  # This bit to silence CRAN warnings
  ..x=..y=..id.col=id=i=NULL
  
  dt <- copy(dt[,.(x = get(..x), y = get(..y), id = get(..id.col))
  ][id == id.val])
  dt <- unique(stats::na.omit(dt))
  
  if (!is.null(rand.val)){
    if (id.val %in% rand.val){
      # Return a random value from 0 to 1 if it's one of the
      # pseudo-observations
      return (function(x) stats::runif(1, min = rand.range[1], max = rand.range[2]))
    } 
  } else if (nrow(dt) == 1){
    # If there's only one observation, return the original value
    return (function(x) mean(unlist(dt$x),na.rm=TRUE))
  } else if (nrow(dt) == 0){
    return (function(x) NA)
  }
  # Otherwise, return a linear map 
  dt <- dt[, .(x,y)]
  if (nrow(stats::na.omit(dt)) == 0){
    return (function(x) NA)
  } else if (nrow(stats::na.omit(dt)) == 1){
    return (function(x) mean(unlist(dt$x),na.rm=TRUE))
  }
  return(FUN(dt, ...))
}
