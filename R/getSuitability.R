#' This function measures how suitable an observation is given set of 
#' metrics that tend to increase the favorability of a location. Z-scores
#' are first calculated for each metric (with an optional threshold). 
#' These are then standardized between [0,1], and a weighted mean between
#' all metrics is performed to find the global suitability. Values of zero
#' indicate that all of the least-favorable values are coincident at that
#' observation; values of one indicate that the location sees all of the most
#' favorable conditions for each metric.
#' 
#' @title Calculate the suitability of a set of observations
#' @param x a numeric vector, data.frame, or data.table with metrics representing
#' values that tend to increase (or decrease, if negative \code{weights} are provided)
#' the favorability of an observation
#' @param group.var (Optional) A character string representing the name of the 
#' column with group IDs, such that z-scores and suitabilities are calculated
#' only within groups. 
#' @param weights (Optional) A named list of numerics or vector of class numeric
#' by which to weigh favorability  observations in the final suitability estimate. 
#' If \code{NULL}, all weights are given equal (positive) importance. 
#' If a negative weight is provided, the assumption is that higher values
#' *decrease* favorability. If vector, it must be be of the 
#' same length as the number of variable columns in \code{x}. If named list,
#' the names should correspond with the favorability metrics to consider, and the
#' values to the strength of the weights. Missing values are assumed to have a
#' zero weight.
#' @param suit.name The name of the output net suitability column. Default is 
#' \code{suit.name = 'Suit'}. If \code{NULL}, only thresheld z-scores are calculated
#' @param stdev The threshold value for z-scores (see \code{\link[lbmech]{rescaleSD}}).
#' Default is \code{stdev = 2}
#' @param keep.vars Should the thresheld z-scores for each variable be kept? Default
#' is \code{keep.vars = FALSE}
#' @importFrom data.table copy
#' @examples 
#' # Creata dummy data
#' x <- data.table(A = c(1:10),
#'                 B = c(11:20),
#'                 C = c(21:30),
#'                 D = c(31:40),
#'                 ID = rep(c("a","b"),5))
#' 
#' # Create weights, don't consider 'C'
#' weights <- list(A = 1,
#'                 B = 2,
#'                 D = -1)
#'                 
#' # Get suitabilities
#' getSuitability(x, group.var = 'ID', weights = weights, keep.vars = TRUE)
#' @export
getSuitability <- function(x, group.var = NULL, weights = NULL,
                           suit.name = "Suit", stdev = 2, keep.vars = FALSE){
  # This bit to silence CRAN warnings
  ..stdev=..weights=NULL
  
  if (methods::is(x,'data.frame') | methods::is(x,'numeric')){
    x <- as.data.table(x)
  }
  
  # Standardize the weights entry into a named list of positive weights,
  x <- copy(x)
  if (is.null(weights)){
    weights <- as.list(rep(1,(ncol(x) - !is.null(group.var) )))
    names(weights) <- names(x)[which(names(x) != suit.name)]
  }
  
  if (methods::is(weights,'numeric')){
    weights <- as.list(weights)
    names(weights) <- names(x)
  } 
  
  if (methods::is(weights,'data.frame')){
    weights <- as.list(weights)
  }
  
  # Negate all variables with negative weights
  weights <- weights[which(weights != 0)]
  negweights <- names(weights)[which(weights < 0)]
  weights <- lapply(weights, abs)
  
  if (length(negweights) != 0){
    x[, (negweights) := .SD * -1, .SDcols = negweights]
  }
  
  # Perform a standard deviational stretch
  x[, (names(weights)) := lapply(.SD, as.double), .SDcols = names(weights)
  ][, (names(weights)) := lapply(.SD, rescaleSD, stdev = stdev),
    by = group.var, .SDcols = names(weights)
  ][, (names(weights)) := lapply(.SD, function(x) (x + ..stdev)/(2*..stdev)),
    by = group.var,
    .SDcols = names(weights)]
  
  if (!is.null(suit.name)){
    x[, (suit.name) := rowSums(.SD * t(unlist(..weights))/sum(unlist(..weights)),
                               na.rm = TRUE), 
      .SDcols = names(weights)]
  }
  
  if (keep.vars == TRUE){
  return(x[,.SD,.SDcols = c(names(weights),suit.name)])
  } else {
    return(x[,(.SD),.SDcols = suit.name])
  }
}
