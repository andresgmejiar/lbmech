#' A function to perform the necessary raster arithmetic to determine the 
#' catchment of a location based on a movement cost function and an optional
#' weight. 
#' \code{\link[lbmech]{getCosts}} must have been run before this tool can be used.
#' 
#' See also \code{\link[lbmech]{mwVoronoi}}
#' 
#' @title Generate catchment polygons for groups of locations
#' @param rasters One of either a character string or multilayer SpatRaster. 
#' If character string, it represents the filepath
#' to the workspace used as \code{dir} for the previous functions.
#' Default is \code{tempdir()} but unless you are not following best
#' practices you will have to change it to your output directory. If multilayer
#' SpatRaster, it should be the output (or identical in form) to
#' the \code{\link[lbmech]{getCosts}} function with \code{"object" \%in\% output}.
#' @param cost A character with the desired pre-calculated
#' cost names (e.g. \code{'dt'} for time, \code{'dW_l'} for work using
#' \code{\link[lbmech]{energyCosts}}).
#' @param direction A character vector containing one or both of \code{c("in","out")}
#' or the singular string 'both'. This determines whether costs to or from the nodes
#' are calculated.
#' @param ids A character vector containing the IDs of the locations to consider
#' in the catchment calculation. Default is \code{ids = NULL}, which considers all
#' of the locations in the \code{rasters} multilayer SpatRaster or the \code{name}
#' shapefile if \code{rasters} is a directory string. 
#' @param name A character string representing the \code{outname} in
#' \code{\link[lbmech]{getCosts}} containing all of the locations of interest. 
#' Ignored if \code{ids} is provided or if \code{rasters} is a multilayer SpatRaster
#' @param w A named numreic vector with multiplicative weights for a given location.
#' Names must correspond with ID strings. Missing locations are assumed to have a weight
#' of 1. 
#' @return A SpatVector with catchment polygons
#' @importFrom terra vect
#' @importFrom terra rast
#' @importFrom terra which.min
#' @importFrom terra as.polygons
#' @importFrom data.table data.table
#' @importFrom data.table rbindlist
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
#' # Calculate catchments for work
#' catchments <- makeCatchment(costRasters,
#'                             cost = 'dW_l')

makeCatchment <- function(rasters = tempdir(), cost, direction = "both",
                          ids = NULL, name = NULL, w = NULL){
  
  # Check to make sure direction is properly called
  if (!(direction %in% c('in','out','both'))){
    stop("'direction' must be one of 'in', 'out', or 'both")
  }
  
  # If no ID list is provided, assume it's all the items in a shapefile or
  # raster stack
  if (is.null(ids)){
    if (methods::is(rasters, 'character')){
      if (is.null(name)){
        stop("List of locations or 'outname' must be provided if input is a directory")
      } else {
        ids <- normalizePath(paste0(rasters,"/World/CostRasters/",name,'.gpkg'),mustWork=FALSE)
        ids <- unlist(vect(ids)$ID)
      }
    } else if (methods::is(rasters, 'SpatRaster')){
      ids <- unique(stringr::str_extract(names(rasters),'(?<=_(to|from)_).+$'))
    }
  }
  
  # Calculate the cost to move to, from, or roundtrip depending on direction type
  to <- 0
  from <- 0
  if (methods::is(rasters,'SpatRaster')){
    # If costs are provided as raster it's easy; just go by the layer names
    if (direction %in% c('in','both')){
      from <- rasters[[paste0(cost,'_from_',ids)]]
    }
    if (direction %in% c('out','both')){
      to <- rasters[[paste0(cost,'_to_',ids)]]
    }
  
  } else if (methods::is(rasters, 'character')){
    # Otherwise, we need to look through all of the IDs in the shapefiles
    rd <- normalizePath(paste0(rasters,"/World/CostRasters"),mustWork=FALSE)
    
    # Get a list of all the IDs per shapefile
    f <- list.files(rd, pattern = '.gpkg$',full.names = TRUE)
    f <- lapply(f, function(x) data.table(file = stringr::str_replace(f,'\\.gpkg','.tif'),
                                          ID = vect(f)$ID))
    f <- rbindlist(f)
    
    # Extract the rows we need
    f <- f[match(ids, f$ID)]
    
    # Read in only those layers depending on the direction of travel
    if (direction %in% c('in','both')){
      to <- rast(lapply(1:nrow(f), 
                        function(x) rast(f$file[x],
                                         lyrs = paste0(cost,'_to_',f$ID[x])))
      )
    }
    if (direction %in% c('out','both')){
      from <- rast(lapply(1:nrow(f), 
                          function(x) rast(f$file[x],
                                           lyrs = paste0(cost,'_from_',f$ID[x])))
      )
    }
  }
  
  # Get total cost, zero is added for one-way travel
  costs <- from + to
  names(costs) <- stringr::str_extract(names(costs),'(?<=_(from|to)_).+$')
  
  # Deal with weights; assume 1 if one is not provided
  ws <- rep(1, length(ids))
  names(ws) <- ids
  if (!is.null(w)){
    ws[match(names(w), names(ws))] <- w
  }
  
  # Find catchments, convert raster to polygon
  zone <- which.min(costs / ws)
  zone <- as.polygons(zone)

  # Change index numbers to ID names
  zone$ID <- ids[zone$which.min]
  zone$which.min <- NULL
  
  return(zone)
}

