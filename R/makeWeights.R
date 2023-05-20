#' Given a distance matrix (symmetrical or not) and weighting function, 
#' calculate a spatial weights matrix with various ways of handling row-sums
#' 
#' @title Calculate spatial weights matrices.
#' @param x A matrix or distance object representing pairwise distances. The
#' distances need not be symmetrical.
#' @param bw A number representing the bandwidth within neighbors are considered.
#' If \code{mode = 'adaptive'}, \code{bw} is the number of nearest neighbors. If \code{mode = 'fixed'},
#' \code{bw} is the radius of the window in the map units.
#' @param mode One of \code{'adaptive'}, which considers a \code{bw} number of nearest 
#' neighbors; or \code{'fixed'}, which considers a fixed bandwidth window of radius \code{bw}.
#' @param weighting One of \code{'membership'}, which considers binary membership 
#' such that neighbors are weighted 1 and non-neighbors 0; \code{'distance'} which 
#' weighs neighbors according to \code{FUN} with the raw distance matrix providing the
#' distance; or \code{'rank'} which uses the rank-distance (i.e. 1 for nearest neighbor,
#' 2 for second nearest...) as the distance variable.
#' @param FUN The distance function. Default is \code{NULL} for \code{'membership'}, and
#' \code{function(x) 1/(offset + x)} otherwise. 
#' @param inf.val When singularities arise, (i.e. whenever the value is 1/0), by what value are
#' they replaced? Default is the \code{FUN} of the lowest non-\code{minval} value.
#' @param offset What value is added to the denominator to prevent singularities from arising
#' (e.g. whenever the value is 1/0)? Larger values imply smaller distance-decay
#' \code{Inf}, what is the value by which they are replaced? Default \code{NULL} uses the
#' value of the smallest neighbor pair from the entire dataset. 
#' @param minval When distances are raw, what is the minimum allowable distance?
#' Default is 50. 
#' @param row.stand Logical or \code{'fuzzy'}. If \code{TRUE} (the default), rows are standardized such 
#' that they sum to one. If \code{'fuzzy'}, rows are standardized as a proportion of the 
#' largest value. 
#' @param clear.mem Logical. Should \code{\link[base]{gc}} be run in the middle of the 
#' calculation? Default is \code{clear.mem} but set as \code{TRUE} if memory limits are a concern. 
#' @importFrom data.table as.data.table
#' @importFrom data.table set
#' @return A matrix of dimensions \code{length(x) x length(x)}.
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
#' @export                      
makeWeights <- function(x, bw, mode = 'adaptive', weighting = 'membership', 
                        FUN = NULL, offset = 0, inf.val = NA, minval = -Inf, 
                        row.stand = FALSE, clear.mem = FALSE) {
  if (is.null(FUN)){
    FUN <- function(x) 1/(offset + x)
  }
  
  # Coerce inputs to matrix
  x <- as.matrix(x)
  if (minval != -Inf) x[x <= minval] <- minval
  
  # x <- as.data.table(x)
  # names(x) <- as.character(1:nrow(x))
  
  # Deal with group membership
  if (mode == 'adaptive'){
    # Convert distances to ranks of distances; need for both
    y <- matrixStats::colRanks(x)
    # Remove observations with rank lower than bw; keep only obs within bw
    if (weighting == 'rank'){
      x <- y
      rm(y)
      x[x >= bw+1] <- NA 
    } else {
      y[y >= bw+1] <- NA 
      x <- x * y/y
      rm(y)
    }
  } else if (mode == 'fixed'){
    if (weighting == 'rank'){
      # Convert distances to ranks of distances; only need for rank
      y <- matrixStats::colRanks(x)
      x[x >= bw] <- NA
      x <- y * x/x
      rm(y)
    } else {
      x[x >= bw] <- NA
    }
  }
  if (clear.mem) gc()
  # Deal with distance weighting
  if (weighting == 'membership'){
    # Binary membership, equal weights
    x[!is.na(x)] <- 1

  } else if (weighting %in% c('distance','rank')){
    # Calculate distances based on raw distances
    x <- FUN(x)
    if (is.null(inf.val)){
      # Replace infinites
      inf.val <-  FUN(min(x[x > minval],na.rm=TRUE))
    }
    x[is.infinite(x)] <- inf.val

  } 
  if (clear.mem) gc()
  if (row.stand == TRUE){
    # Traditional Row standardization, where all rows sum to one
    x <- x/rowSums(weights,na.rm=TRUE)
  }
  
  if (row.stand == 'fuzzy'){
    # Each row has at least one neighbor, with the remainder being a fraction
    # of the maximum
    x[is.na(x)] <- 0
    x <- t(apply(x,1,scales::rescale,to=c(0,1)))
  }
  x[is.na(x)] <- 0
  gc()
  return(x)
}
