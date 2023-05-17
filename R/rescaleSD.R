#' Perform a standard deviational stretch on any numeric coercible to a 
#' vector object. Z-scores are first calculated, and values that exceed the 
#' threshold are reassigned the theshold value
#' 
#' @title Standard deviational stretch
#' @param x The input value, vector, list, or anything coercible to a vector with 
#' \code{\link[base]{unlist}}. Values must be numeric
#' @param stdev The threshold number of standard deviations. Input values \code{n}
#' whith a z score whose absolute value is  beyond this threshold will be 
#' reassigned \code{x[n] <- stdev * z/abs(z)}
#' @return A numeric vector with threshold-limited z-scores
#' @examples 
#' # Create a vector
#' x <- rnorm(10, mean = 10, sd = 10)
#' 
#' y <- rescaleSD(x, stdev = 1)
#' @export
rescaleSD <- function(x,stdev=2) {
  if (length(stats::na.omit(x)) == 1){
    return(0.5)
  } else {
    x <- unlist(x)
    out <- (x - mean(x,na.rm=TRUE))/(stats::sd(x,na.rm=TRUE))
    out[out < -1 * stdev] <- -1 * stdev
    out[out > stdev] <- stdev
    out <- (out + stdev) / (2 * stdev)
    if (all(is.na(out))){
      return(0.5)
    } else {
      return(out)
    }
  }
}