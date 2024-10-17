#' A function to compare implied 2D topologies based on different distance
#' matrices for the same data points
#' 
#' The function relies on \code{\link[MASS]{isoMDS}} from the MASS package to 
#' generate pseudo-coordinates for distance matrices using a 2D non-metric 
#' Multi Dimensional Scaling. Since the output orientation can be unpredictable,
#' normMDS aligns all topologies based on a common user-defined axis 
#' (the ends of which will appear equidistant and equioriented on an output plot),
#' and selects a chiral reference point to ensure that plots are oriented in a 
#' common orientation. 
#' 
#' @title Normalized Non-Metric Multidimensional Scaling
#' @param dists A list (named, ideally) of distance matrices with the same points
#' @param xy (Optional) The points in geographic space. This can be a 
#' data.frame, data.table, SpatVect, or SpatialPointsDataFrame object. If provided, 
#' axes orientations and chirality will be taken with respect to \code{xy}.
#' @param id Character. The name of the column in \code{xy} containing the unique
#' name or id for each observation. Default is \code{xy = 'ID'}.
#' @param ax The index (if numeric) or Unique ID (if character) of the two points
#' defining the reference axis. The first point will be placed at the origin.
#' If \code{xy = NULL}, the second point will be placed at v = 0. Default is 
#' \code{ax = 1:2}.
#' @param chir The index (if numeric) or Unique ID (if character) of the point
#' defining the reference handendess/chirality. This must be a point not in \code{ax}.
#' Default is \code{ax = 3}
#' @param x Character. The name of the column in \code{xy} containing the x 
#' coordinate for each observation. Default is \code{xy = 'x'}. Ignored if 
#' \code{xy} is a SpatVect or Spatial object. 
#' @param y Character. The name of the column in \code{xy} containing the y 
#' coordinate for each observation. Default is \code{xy = 'y'}. Ignored if 
#' \code{xy} is a SpatVect or Spatial object. 
#' @param ... Additional parameters to pass to \code{\link[MASS]{isoMDS}}
#' @return A named list, with each entry containing a data.table of estimated
#' pseudo-coordinates in \emph{(u,y)} space. One unit is equal to the distance
#' between both points in the \code{ax} parameter. If \code{xy} was provided,
#' the first entry \code{'ref'} will be the geographic data. 
#' @importFrom terra geomtype
#' @importFrom terra vect
#' @importFrom terra geom
#' @importFrom data.table data.table
#' @importFrom data.table as.data.table
#' @examples 
#' 
#' # Create fake data
#' xy <- data.table(ID = c('A','B','C','D','E'),
#'                  x = runif(5, -100, 100),
#'                  y = runif(5, -100, 100))
#'                  
#' # Calculate a distance matrix
#' dists <- as.matrix(stats::dist(xy[,.(x,y)], upper = TRUE, diag = TRUE) )
#' rownames(dists) <- xy$ID
#' colnames(dists) <- xy$ID               
#'                  
#'                  
#' # Make a list of five fake distance matrices from five fake coordinates
#' dists <- list()
#' for (i in 1:5){
#'   fake <- data.table(ID = c('A','B','C','D','E'),
#'                  x = runif(5, -100, 100),
#'                  y = runif(5, -100, 100))
#'                  
#'   fake <- as.matrix(stats::dist(fake[,.(x,y)], upper = TRUE, diag = TRUE) )
#'   rownames(fake) <- xy$ID
#'   colnames(fake) <- xy$ID 
#'   dists[[i]] <- fake
#' }
#' 
#' # Get pseudo-coordinates
#' normMDS(dists, xy = xy)
#' @export
normMDS <- function(dists, xy = NULL, id = 'ID', ax = 1:2, chir = 3, 
                    x = 'x', y = 'y', ...){
  # Deal with missing names
  if (is.null(names(dists))){
    names(dists) <- 1:length(dists)
  }
  iternames <- names(dists)
  
  
  # If no reference data is provided, use the output of the first MDS
  if (is.null(xy)){
    invisible({
    xy <- as.data.table(MASS::isoMDS(dists[[1]], k = 2, ...)$points,
                        keep.rownames = TRUE)
    })
    names(xy) <- c(id,x,y) 
    iternames <- iternames[2:length(iternames)]
  } else if (is(xy, 'SpatialPointsDataFrame') | is(xy,'SpatVector')){
    # If a spatial object is provided, coerce to data.table
    if (!is(xy,'SpatVector')) xy <- vect(xy)
    if (geomtype(xy) != 'points') stop("'xy' can only be point locations.")
    
    xy <- data.table('id' = xy[[id]], geom(xy))
    names(xy)[1] <- id
  } else xy = as.data.table(xy)
  
  # Get consistent naming
  xy <- xy[, `:=`(x = get(..x), y = get(..y))][]  
  
  # Get the names of the reference axes if indexes were provided
  if (is(ax,'numeric')){
    ax <- xy[ax, get(id)]
  }
  
  # And likewise for chiral reference
  if (is(chir,'numeric')){
    chir <- xy[chir, get(id)]
  }
  
  # Check to make sure chiral reference is different from axes points
  if (chir %in% ax){
    stop('Chiral refernce point cannot be in reference axis')
  }
  
  # Get the original coordinates of the reference standard
  p <- xy[match(ax, get(id))]
  
  # Shift original data to chosen centroid
  xy[, `:=`(x = x - ..p[1]$x, 
            y = y - ..p[1]$y)]
  p <- xy[match(ax, get(id))]
  
  # Scale to unity
  S <- sqrt(p[2]$x^2 + p[2]$y^2)
  xy[, `:=`(x = x / S, y = y / S)]
  p <- xy[match(ax, get(id))]
  
  # Get the angle of orientation of chosen axis
  angle <- atan(p[2]$y/p[2]$x)
  
  # If it's just mds coordinates, the original orientation has no bearing, so
  # rotate so the chosen axis is algined with v = 0
  if (is.null(xy)){
    xy[, `:=`(theta = atan(y/x), l = sqrt(x^2 + y^2))
    ][, `:=`(x = l * cos(theta - angle), y =  l * sin(theta - angle))]
    xy[is.na(pts)] <- 0
    angle <- 0
    p <- xy[match(ax, get(id))]
  }
  
  # Instantiate list where we'll store the transformations
  mds <- list()
  firstname <- ifelse(length(iternames) == length(dists), 'ref', names(dists)[1])
  
  mds[[firstname]] <- xy[,.(get(id), u = x, v = y)]
  names(mds[[firstname]])[1] <- id
  
  for (i in iternames){
    # Perform 2D Multidimensional Scaling, convert to data.table
    invisible({
      pts <- as.data.table(MASS::isoMDS(dists[[i]], k = 2, ...)$points,
                         keep.rownames = TRUE)
    })
    names(pts) <- c(id,'u','v')
    
    # Get points for chosen reference axis
    q <- pts[match(ax, get(id))]
    
    # Shift coordinates so the chosen centroid is aligned
    dx <- q[1]$u - p[1]$x 
    dy <- q[1]$v - p[1]$y 
    
    pts[, `:=`(u = u - dx,
               v = v - dy)]
    q <- pts[match(ax, get(id))]
    
    # Calculate the difference between the angle of the chosen axis in this one
    # and the reference one
    dtheta <- atan(q[2]$v/q[2]$u) - angle
    
    # Rotate data so the chosen axes are aligned
    pts[, `:=`(theta = atan(v/u), l = sqrt(u^2 + v^2))
    ][, `:=`(u = l * cos(theta - ..dtheta), v =  l * sin(theta - ..dtheta))]
    pts[is.na(pts)] <- 0
    
    # Scale data so the chosen reference distances are equivalent
    S <- sqrt(q[2]$u^2 + q[2]$v^2) / sqrt(p[2]$x^2 + p[2]$y^2)
    
    pts[, `:=`(u = u / ..S,
               v = v / ..S)]
    
    # We need to check for chirality. 
    if (q[2]$v == 0){
      # If the chosen axis is a vertical line see if the chiral reference point is 
      # on the same side of the reference axis.
      flip <- xy[get(id) == chir, x > p[2]$x] != pts[get(id) == chir, u > p[2]$x]
      asympt <- TRUE # We'll have to deal with verticality again; might as well tag it
      
    } else {
      # Likewise for all other types of lines, but here calculate the slope too
      m <- p[2]$y / p[2]$x
      flip <- xy[get(id) == chir,  x * ..m > y] != pts[get(id) == chir,  u * ..m > v] 
      asympt <- FALSE
    }
    
    # If the points are not on the same side, we'll have to flip this one
    if (flip == TRUE & asympt == FALSE){
      if (m == 0){
        # If the slope is zero, then just negate the values
        pts[, `:=`(v = -v)]
      } else {
        # If the slope is not, rotate the data about the chosen axis
        pts[, `:=`(a = abs(v - u * m), b = abs(u - v / m))]
        left <- ifelse(m > 0, 1, -1)
        pts[,up := ifelse(b > v, 1, -1)
        ][, `:=`(u = u - a * up * left,
                 v = v + b * up)]
      }
    } else if (flip == TRUE & asympt == TRUE){
      # If we flip AND it's a vertical line, just negate the x coordinate
      pts[, `:=`(u = -u)]
    }
    
    # Add to stored list
    mds[[i]] <- pts[,.(get(id), u, v)]
    names(mds[[i]])[1] <- id
  }
  
  return(mds)
}

