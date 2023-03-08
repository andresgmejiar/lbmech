#' A function that for a given region imports all cells from the
#' transition \code{.fst} files. If such files have not yet been generated,
#' they can be created by passing along the necessary parameters to this
#' function as with \code{\link[lbmech]{calculateCosts}}.
#' 
#' The default parameters are sufficient for a workflow involving calculating
#' costs with the \code{\link[lbmech]{energyCosts}} function. However, if
#' non-energetic analyses are desired, the user must define their own.
#' 
#' @title Import a world where movement is possible
#' @param region An object of class SpatRraster, Raster* or SpatialPolygons*
#' representing the total area where movement is possible.
#' @param banned An object of class Raster* or SpatialPolygons*
#' representing the total area where movement is \emph{prohibited}. Must lie within
#' the area defined by \code{polys}
#' @param dir A filepath to the directory being used as the workspace, the same
#' one instantiated with \code{\link[lbmech]{defineWorld}}. 
#' Default is \code{tempdir()} but unless the analyses will only be performed a few
#' times it is highly recommended to define a permanent workspace.
#' @param vars The variable names to import.
#' @param costFUN A cost function such as (\code{\link[lbmech]{timeCosts}} or
#' \code{\link[lbmech]{energyCosts}}). The input to such a function should be 
#' a \code{data.table} object with column names present in the makeWorld file
#' (initially \code{c('x_i','x_f','y_i','y_f','z_i','z_f','dz','dl','dr')}), and
#' the output should be an data.table object with the same number of rows and the
#' desired output variables as the only column names. Constants can be passed in the 
#' \code{...} slot. Default is 
#' \code{costFUN = \link[lbmech]{energyCosts}}.
#' @param ... Additional arguments to pass to \code{\link[lbmech]{calculateCosts}}
#' and \code{costFUN}. 
#' @return An object of class data.table containing at least three columns:
#'
#' (1) \code{$from}, a character string of all possible origin cells in format "x,y",
#'
#' (2) \code{$to},  a character string of all possible destination cells in format "x,y"
#'
#' (3+) a numeric representing the imported costs(s)
#' 
#' @importFrom terra vect
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
#' # Select all tiles that exist between x = (12000,16000) and y = (32000,36000)
#' tiles <- ext(c(12000,16000,32000,36000))
#' tiles <- as.polygons(tiles)
#' crs(tiles) <- crs(grid)
#' tiles <- whichTiles(region = tiles, polys = grid)
#' 
#' #' # Make a world but limit it to the DEM grid size
#' defineWorld(source = grid, cut_slope = 0.5, 
#'             res = res(dem), dir = dir, overwrite=TRUE)
#' 
#' # Make a world but limit it to the DEM grid size
#' world <- importWorld(grid[8,], dir = dir, vars = 'dE_l', costFUN = energyCosts,
#' m = 70, v_max = 1.5, BMR = 76, k = 3.5, s = 0.05, l_s = 1,
#' L = 0.8)
#' @export
importWorld <- function(region, banned = NULL, dir = tempdir(), 
                        vars = NULL,  costFUN = NULL, ...){
  # This bit is to silence the CRAN check warnings for literal column names
  from=to=dz=x_f=x_i=y_f=y_i=z_f=z_i=x=y=long_f=lat_f=long_i=lat_i=NULL
  
  
  
  # Set up the local environment from directory
  dir <- normalizePath(paste0(dir,"/World"),mustWork=FALSE)
  subdirs <- c("/Raw","/Local","/Diff")
  subdirs <-  normalizePath(paste0(dir,subdirs),mustWork=FALSE)
  
  callVars <- readRDS(normalizePath(paste0(dir,"/callVars.gz"),mustWork=FALSE))
  list2env(lapply(as.list(callVars),unlist),environment())
  
  source <- vect(normalizePath(paste0(dir,"/z_sources.gpkg"),mustWork=FALSE))
  grid <- vect(normalizePath(paste0(dir,"/z_grid.gpkg"),mustWork=FALSE))
  z_fix <- importRST(normalizePath(paste0(dir,"/z_fix"),mustWork=FALSE))
  
  
  
  vars <- c('from','to',vars)
  
  # Make sure the relevant sector and cost tiles have been created
  tiles <- whichTiles(region,grid,tile_id = 'id')
  if (is.null(costFUN)){
    makeWorld(tiles,
              dir=normalizePath(stringr::str_remove(dir,'World$'),mustWork=FALSE))
  } else {
    calculateCosts(tiles,
                   costFUN = costFUN, 
                   dir = normalizePath(stringr::str_remove(dir,'World$'),mustWork=FALSE),
                   ...)
  }
  
  # Convert the regions into a list of the allowable cells
  region <- regionMask(region, z_fix =  z_fix)
  if (!is.null(banned)){
    banned <- regionMask(banned, z_fix = z_fix)
  } else{
    banned <- c("NULL")
  }
  
  # Create an empty data.table to which to add the imported tensors
  Edges <- data.table()
  
  pb <- utils::txtProgressBar(max = length(tiles), style = 3)
  
  # Import tensors one-by-one, keeping only the allowable cells
  
  for (i in seq(1,length(tiles))) {
    import <- fst::read_fst(normalizePath(paste0(subdirs[3],'/',tiles[i],'.fst')),
                            columns = vars,
                            as.data.table = TRUE)
    import <- import[(!(from %in% banned) | !(to %in% banned)) &
                       ((from %in% region) | (to %in% region)),]
    Edges <- rbind(Edges, import)
    utils::setTxtProgressBar(pb,i)
  } 
  rm(import)
  Edges <- unique(Edges)
  
  return(Edges[])
}