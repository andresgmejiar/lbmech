#' Given a a description of relationships,  
#' calculate a spatial weights matrix with various ways of handling row-sums
#' 
#' The description of relationships can either be a vector describing to which
#' group an observation belongs, or a distance matrix and associated decay function. 
#' @title Calculate weights matrices.
#' @param x A vector representing group membership, or a matrix/distance object 
#' representing pairwise distances. The distances need not be symmetrical.
#' @param ID A vector given the unique observation ID names. Default is 
#' \code{ID = 1:length(x)}. Ignored if \code{x} is a matrix or distance object.
#' @param bw A number representing the bandwidth within neighbors are considered.
#' If \code{mode = 'adaptive'}, \code{bw} is the number of nearest neighbors. 
#' If \code{mode = 'fixed'}, \code{bw} is the radius of the window in the map units.
#' Ignored if \code{x} is a vector, required if \code{weighting = 'membership'}, 
#' and assumed to be \code{bw = Inf} otherwise.
#' @param mode One of \code{'adaptive'}, which considers a \code{bw} number of nearest 
#' neighbors; or \code{'fixed'}, which considers a fixed bandwidth window of radius \code{bw}.
#' Ignored if \code{x} is a vector.
#' @param weighting One of \code{'membership'}, which considers binary membership 
#' such that neighbors are weighted 1 and non-neighbors 0; \code{'distance'} which 
#' weighs neighbors according to \code{FUN} with the raw distance matrix providing the
#' distance; or \code{'rank'} which uses the rank-distance (i.e. 1 for nearest neighbor,
#' 2 for second nearest...) as the distance variable. Ignored if \code{x} is a vector.
#' @param FUN The distance function. Default is \code{NULL} for \code{'membership'}, and
#' \code{function(x) 1/(offset + x)} otherwise. Ignored if \code{x} is a vector.
#' @param inf.val When singularities arise, (i.e. whenever the value is 1/0), by what value are
#' they replaced? Default is the \code{FUN} of the lowest non-\code{minval} value.
#' Ignored if \code{x} is a vector.
#' @param offset What value is added to the denominator to prevent singularities from arising
#' (e.g. whenever the value is 1/0)? Larger values imply smaller distance-decay. 
#' Default is \code{offset = 0}. Ignored if \code{x} is a vector.
#' @param minval When distances are raw, what is the minimum allowable distance?
#' Default is \code{0}. Ignored if \code{x}. Use this if you don't want to offset values otherwise. 
#' @param def.neigh Numeric. At what distance (in the map units) are observations definitely neighbors?
#' All distances are subtracted by this value, and all resulting distances less than zero are reassigned
#' to \code{minval}.
#' @param row.stand Logical or \code{'fuzzy'}. If \code{TRUE} (the default), rows are standardized such 
#' that they sum to one. If \code{'fuzzy'}, rows are standardized as a proportion of the 
#' largest value. If \code{x} is a vector, \code{row.stand} must be logical. 
#' @param clear.mem Logical. Should \code{\link[base]{gc}} be run in the middle of the 
#' calculation? Default is \code{clear.mem == FALSE} but set as \code{TRUE} if 
#' memory limits are a concern. Ignored if \code{x} is a vector.
#' @importFrom data.table as.data.table
#' @importFrom data.table set
#' @importFrom data.table CJ
#' @return A matrix of dimensions \code{length(x) x length(x)}.
#' @examples 
#' 
#' # Example 1: Calculate group-based weights
#' # Generate dummy group names
#' groups <- 1:4
#' 
#' # Create 10 dummy values assigned to those groups
#' x <- sample(groups, 10, replace = TRUE)
#' 
#' # Create group membership weights matrix
#' weights <- makeWeights(x)
#'
#'
#' # Example 2: Calculate distance-based weights 
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
makeWeights <- function(x, ID = NULL, bw = NULL,
                        mode = 'adaptive', weighting = 'membership', 
                        FUN = NULL, offset = 0, inf.val = NA, minval = 0, 
                        def.neigh = 0, row.stand = FALSE, clear.mem = FALSE) {
  
  # First do group membership vectors, then do distance matrices 
  if (is.vector(x)){
    # Deal with IDs. If not provided create them. If they are, check to make
    # sure they're the same length as x
    if (is.null(ID)) ID <- 1:length(x)
    if (length(x) != length(ID)){
      stop("'x' and 'index' must be the same length")
    }
    if (length(x) != length(unique(ID))){
      stop("'ID' vector most contain all unique IDs")
    }
    # To assign group values, cross join values and indexes on themselves,
    # then wherever i == j, assign 1 otherwise 0
    cj <- cbind(CJ(i = x,
                   j = x,
                   sorted = FALSE, unique = FALSE),
                CJ(i_id = ID,
                   j_id = ID,
                   sorted = FALSE, unique = FALSE)
    )[,Weight := fifelse(i == j, 1, 0)][]
    
    # Perform row standardization if parameter is set
    if (row.stand) cj[, Weight := Weight/sum(Weight), by = 'i_id']
    
    # Change ID column names to i & j
    cj[, `:=`(i = i_id, j = j_id)]
    
    # Convert table to matrix
    cj <- as.matrix(stats::xtabs(Weight ~ i + j, data = cj))
    
    # Return matrix, sorted in the original order
    return(cj[as.character(ID),
              as.character(ID)])
  } else {
    if (!is.null(bw)){
      if (is.infinite(bw)) bw <- NULL
    } 
    
    if (is.null(bw) & weighting == 'membership'){
      stop("argument 'bw' is missing, with no default")
    }
    if (is.null(FUN)){
      FUN <- function(x) 1/(offset + x)
    }
    
    # Coerce inputs to matrix
    x <- as.matrix(x)
    if (minval != 0) x[x <= minval] <- minval
    
    # Adjust based on definite neighbor value
    if (def.neigh > 0) {
      x <- x - def.neigh
      x[x < 0] <- minval
    }
    
    # x <- as.data.table(x)
    # names(x) <- as.character(1:nrow(x))
    
    # Deal with group membership
    if (mode == 'adaptive'){
      # Convert distances to ranks of distances; need for both
      y <- matrixStats::colRanks(x, ties.method = "min")
      # Remove observations with rank lower than bw; keep only obs within bw
      if (weighting == 'rank'){
        x <- y
        rm(y)
        if (!is.null(bw)) x[x >= bw+1] <- NA 
      } else {
        if (!is.null(bw)) y[y >= bw+1] <- NA 
        x <- x * y/y
        rm(y)
      }
    } else if (mode == 'fixed'){
      if (weighting == 'rank'){
        # Convert distances to ranks of distances; only need for rank
        y <- matrixStats::colRanks(x, ties.method = "min")
        if (!is.null(bw)) x[x >= bw] <- NA
        x <- y * x/x
        rm(y)
      } else if (!is.null(bw)) {
        x[x >= bw] <- NA
      }
    }
    if (clear.mem) gc()
    # Deal with distance weighting
    if (weighting == 'membership'){
      # Binary membership, equal weights
      x[!is.na(x)] <- 1
      
    } else if (weighting %in% c('distance','rank')){
      # Calculate transformed distances
      x <- FUN(x)
      if (is.null(inf.val)){
        # Replace infinites
        if (def.neigh > 0){
          inf.val <- 1/def.neigh
        } else {
          inf.val <-  FUN(min(x[x > minval],na.rm=TRUE))
        }
        x[is.infinite(x)] <- inf.val
        
      } 
      if (clear.mem) gc()
      if (row.stand == TRUE){
        # Traditional Row standardization, where all rows sum to one
        x <- x/rowSums(x,na.rm=TRUE)
      }
      
      if (row.stand == 'fuzzy'){
        # Each row has at least one neighbor, with the remainder being a fraction
        # of the maximum
        x[is.na(x)] <- 0
        x <- t(apply(x,1,scales::rescale,to=c(0,1)))
      }
    }
    x[is.na(x)] <- 0
    gc()
    return(x)
  }
}