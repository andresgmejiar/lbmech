#' Calculate the weighted average that maximizes the maximum likelihood estimator

#' @title Compute weighted means and medians.
#' @param x A numeric vector with the values to be averaged.
#' @param w A vector of weights of the same length as \code{x}.
#' @param type A character string, either: 'mean' or 'L2' (the default),
#' such that a weighted mean is computed; or 'median' or 'L1', such
#' that a weighted median is computed using \code{\link[DescTools]{Median}}.
#' @param na.rm Should \code{NA} values be ignored? Default is \code{na.rm = FALSE}.
#' @return A numeric value, representing the weighted mean or median
#' @examples 
#' 
#' # Create dummy values
#' x = runif(10,1,100)
#' w = runif(10,1,10)
#' 
#' # Compute weighted mean
#' average(x, w, type = 'mean')
#' 
#' # Compute weighted median
#' average(x, w, type = 'median')
#' @export
average <- function(x, w = rep(1,length(x)), type = 'mean',na.rm = FALSE){
  if (type == 'mean' | type == 'L2'){
    return(sum(x * w/sum(w,na.rm=na.rm), na.rm=na.rm))
  } else if (type == 'median' | type == 'L1'){
  return(DescTools::Median(x, weights = w, na.rm = na.rm))
  }
}