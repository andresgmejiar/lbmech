#' A function to automatically perform the raster arithmetic
#' necessary to calculate the cost-of-travel for paths with multiple
#' waypoints, and the predicted cost of taking a detour to any
#' arbitrary point in the landscape (a 'corridor').
#' \code{\link[lbmech]{getCosts}} must have been run before this tool can be used.
#'
#' @title Calculate cost corridors for a path
#' @param rasters One of either a character string or multilayer SpatRaster. 
#' If character string, it represents the filepath
#' to the workspace used as \code{dir} for the previous functions.
#' Default is \code{tempdir()} but unless you are not following best
#' practices you will have to change it to your output directory. If multilayer
#' SpatRaster, it should be the output (or identical in form) to
#' the \code{\link[lbmech]{getCosts}} function with \code{"object" \%in\% output}.
#' @param order A character vector containing the desired path in
#' order of visited nodes. For example, to visit "A" then "B" then "C" then "A"
#' the vector would be \code{c("A","B","C","A")}. Note that these MUST correspond
#' to the ID names for the \code{from} features used in the \code{\link[lbmech]{getCosts}}
#' function and must have previously been calculated
#' @param costs A character vector containing any combination of pre-calculated
#' cost names (e.g. dt for time, dW_l for work using \code{\link[lbmech]{energyCosts}})
#' if the input world data.table is the output of of the \code{\link[lbmech]{calculateCosts}}
#' function. This selects which types of costs will be calculated. 
#' \code{costs = 'all'} is shorthand for \code{costs = c("dt","dW_l","dE_l")}
#' while \code{costs = 'energetics'} is shorthand for \code{c("dW_l","dE_l")}.
#' Default is \code{'all'}. Note that these must have previously been calculated.
#' @param name A character vector representing the \code{outname} in
#' \code{\link[lbmech]{getCosts}}. If none is provided, it will be the name of the
#' variable passed to the function in the \code{from} slot.
#' If \code{length(costs) == 1}, a Raster* If \code{length(costs) > 1}
#' a list of Raster* with one slot for each \code{cost}.
#' @return Rasters representing cost corridors.
#' @importFrom terra rast
#' @importFrom terra global
#' @examples 
#' # Generate a DEM
#' n <- 5
#' dem <- expand.grid(list(x = 1:(n * 100),
#'                         y = 1:(n * 100))) / 100
#' dem <- as.data.table(dem)
#' dem[, z := 250 * exp(-(x - n/2)^2) + 
#'       250 * exp(-(y - n/2)^2)]
#' dem <- rast(dem)
#' ext(dem) <- c(10000, 20000, 30000, 40000)
#' crs(dem) <- "+proj=lcc +lat_1=48 +lat_2=33 +lon_0=-100 +datum=WGS84"
#' 
#' # Export it so it doesn't just exist on the memory
#' dir <- tempdir()
#' writeRaster(dem, paste0(dir,"/DEM.tif"),overwrite=TRUE)
#' 
#'
#' # Import raster, get the grid
#' dem <- rast(paste0(dir,"/DEM.tif"))
#' grid <- makeGrid(dem = dem, nx = n, ny = n, sources = TRUE)
#' 
#' # Import the data lying between x = (12000,16000) and y = (32000,36000)
#' region <- ext(c(12000,16000,32000,36000))
#' region <- as.polygons(region)
#' crs(region) <- crs(grid)
#' 
#' # Generate five random points that fall within the region
#' points <- data.table(ID = 1:5,
#'                      x = runif(5, ext(region)[1], ext(region)[2]),
#'                      y = runif(5, ext(region)[3], ext(region)[4]))
#'                      
#' # Make a world but limit it to the DEM grid size
#' defineWorld(source = grid, cut_slope = 0.5, 
#'             res = res(dem), dir = dir, overwrite=TRUE)
#'             
#' # Calculate cost rasters
#' costRasters  <- getCosts(region, from = points, dir = dir,
#'                          destination = 'all',
#'                          polygons = 'center',
#'                          costs = 'all', costFUN = energyCosts,
#'                          output = c('object','file'),
#'                          m = 70, v_max = 1.5, BMR = 76, k = 3.5, s = 0.05, l_s = 1,
#'                          L = 0.8)
#'                         
#' # Calculating the corridors from a list of RasterStacks,
#' # with path 1 -> 2 -> 4 -> 1 -> 5
#' corridors <- makeCorridor(rasters = costRasters, order = c(1,2,5,1,4),)
#' 
#' #### Example 2:
#' # Calculating the corridors from a workspace directory
#' # with path 1 -> 2 -> 4 -> 1 -> 5
#' corridors <- makeCorridor(rasters = dir, name = 'points', 
#'                           order = c(1,2,5,1,4))
#' @export
makeCorridor <- function(rasters = tempdir(), order, costs = "all",
                         name = NULL){
  # Silence CRAN
  FID=Vector=NULL
  
  # Unlike getCosts, makeCorridor only works if the previous steps have
  # already been completed. It needs either the output list of cost
  # rasters from getCosts OR the directory containing the output
  # cost raster files for all desired nodes.
  if (all(costs == 'all')){
    costs <- c("dt","dW_l","dE_l")
  } else if (all(costs == "energetics")){
    costs <- c("dW_l","dE_l")
  }
  
  
  # Get the names of files/rasters that will be needed for each leg
  starts <- paste0("_from_",order[seq(1,(length(order)-1))])
  stops <-  paste0("_to_",order[seq(2,length(order))])
  
  # Empty output list, as in getCosts
  outList <- rast()
  for (cost in costs){
    # Iterate over costs, again as in getCosts
    
    if (methods::is(rasters, "SpatRaster")){
      # If the input is a SpatRaster, then we can just
      # get the appropriate rasters from the keys
      from <- rasters[[paste0(cost,starts)]]
      to <- rasters[[paste0(cost,stops)]]
      
    } else if (methods::is(rasters,'character')){
      # If the input is a directory with rasters, then we need to import
      # each individual raster and stack them
      rd <- normalizePath(paste0(rasters,"/World/CostRasters"),mustWork=FALSE)
      
      from <- rast(normalizePath(paste0(rd,"/",name,'.tif'),mustWork=FALSE),
                        lyrs = paste0(cost,starts))
      
      to <- rast(normalizePath(paste0(rd,"/",name,'.tif'),mustWork=FALSE),
                      lyrs = paste0(cost,stops))
    }
    
    # The cost for each leg is From Origin + To Destination
    corridor <- from + to
    mins <- unlist(global(corridor,'min',na.rm=TRUE))
    names(mins) <- NULL
    
    # To ensure that the path from the first origin to the last destination
    # contains the same minimum travel cost, subtract the minimum cost per
    # leg from each leg's raster
    corridor <- corridor - mins
    
    # The total corridor is the minimum value at each cell from each
    # leg's raster. Add the sum of the minimums to represent the total
    # minimum cost, not just detour.
    corridor <- min(corridor,na.rm=TRUE) + sum(mins)
    names(corridor) <- cost
    
    # Add to the output list
    outList <- suppressWarnings(c(outList, corridor))
  }
  return(outList)
}