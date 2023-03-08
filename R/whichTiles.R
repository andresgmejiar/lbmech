
#' Get the names of tiles that would be needed to perform an analysis
#' over a given region-of-interest within the maximum possible extent defined
#' by a source grid.
#'
#' @title Identify necessary tiles/sectors
#' @param region An object of class SpatRaster, SpatVector Raster*, Spatial*, 
#' data.frame, or data.table indicating the region-of-interest. If input is of class
#' SpatialPoints*, data.table, or data.frame the \code{region} represents sectors containing the
#' individual points.
#' @param polys A polygon of class SpatVector representing
#' the partitioning grid for the maximum possible area, in the same format as the
#' output of the \code{\link[lbmech]{makeGrid}} function.
#' @param tile_id a character string representing the name of the column
#' in the \code{polys} polygon containing the unique Tile IDs. Default is \code{tile_id = 'TILEID'}
#' @param x a character string representing column name containing the "x"
#' coordinates. Required for SpatialPoints*, data.frame, and data.table \code{region},
#' otherwise ignored.
#' @param y a character string representing column name containing the "y"
#' coordinates. Required for SpatialPoints*, data.frame, and data.table \code{region},
#' otherwise ignored.
#' @return A character vector containing the TILEIDs overlapping with the \code{region}
#' @importFrom terra ext
#' @importFrom terra crs
#' @importFrom terra crs<-
#' @importFrom terra intersect
#' @importFrom terra vect
#' @importFrom terra rast
#' @importFrom terra as.polygons
#' @importFrom raster SpPolygons
#' @examples 
#' 
#' #### Example 1:
#' # If the grid is the product if the makeGrid function
#' # Make the grid
#' n <- 6
#' dem <- rast(ncol = n * 600, nrow = n * 600, vals = 1)
#' ext(dem) <- c(1000, 2000, 3000, 4000)
#' crs(dem) <- "+proj=lcc +lat_1=48 +lat_2=33 +lon_0=-100 +datum=WGS84"
#' # Export it so it doesn't just exist on the memory
#' dir <- tempdir()
#' writeRaster(dem, paste0(dir,"/DEM.tif"),overwrite=TRUE)
#' 
#' # Import raster, get the grid
#' dem <- rast(paste0(dir,"/DEM.tif"))
#' grid <- makeGrid(dem = dem, nx = n, ny = n, sources = TRUE)
#' 
#' # Select five random points that fall within the grid
#' points <- data.table(x = runif(5, ext(dem)[1], ext(dem)[2]),
#'                      y = runif(5, ext(dem)[3], ext(dem)[4]))
#' 
#' tile_list <- whichTiles(region = points, polys = grid)
#' 
#' #### Example 2 (Do not execute):
#' ## If it is a custom polygon "polys", where Tile IDs are stored in 
#' ## a "NAME" column, and coordinates in "Easting" and "Northing"
#' # tile_list <- whichTiles(region = points, grid = polys, 
#' #                         tile_id = "NAME", x = "Easting", y = "Northing")
#' @export
whichTiles <- function(region, polys, tile_id = "TILEID", x = "x", y = "y"){
  # This bit is to silence the CRAN check warnings for literal column names
  ..x=..y=NULL
  #
  
  if ("data.frame" %in% class(region)){
    region <- vect(region, geom = c(x,y), crs = crs(polys), keepgeom=TRUE)
  }   
  if (methods::is(region, "Raster")){
    region <- rast(region)
  }
  if (methods::is(region,"SpatRaster")){
    proj <- crs(region)
    region <- ext(region)
    region <- as.polygons(region)
    crs(region) <- proj
  }
  if (methods::is(region,'Spatial')){
    region <- vect(region)
  } 
  if (methods::is(region,'SpatVector')){
    region <- region[,NA]
  }
  region <- as.data.table(intersect(region,polys))
  tiles <- unique(region[,get(tile_id)])
  return(tiles)
}