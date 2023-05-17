#' Calculate the geodesic distortion on a raster
#' 
#' @title Get geodesic distortion for rasters
#' @param z An input SpatRaster with a known projection
#' @importFrom data.table as.data.table
#' @importFrom terra res
#' @importFrom terra crs
#' @importFrom terra project
#' @importFrom terra rast
#' @return A SpatRaster with three layers:
#' (1) \code{'lx'} with the horizontal pixel size
#' 
#' (2) \code{'ly'} with the vertical pixel size
#' 
#' (3) \code{'A'} with the pixel area
#' @examples 
#' # Generate dummy dem, assign it a projection
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
#' # Get the geodesic distorsion
#' distort <- getDistortion(dem)
#' @export
getDistortion <- function(z){
  # This bit to silence CRAN warnings
  x=ymin=ymax=xmin=y=xmax=lx=ly=A=NULL
  
  # Use rastToTable to get x-y coordinates, keep projection and resolution
  z <- rastToTable(z)
  l_p <- res(z)
  proj <- crs(z)
  
  # Convert to data.table, calculate the coordinates of pixel borders
  z <- as.data.table(z)
  z$xmin <- z$x - l_p[1]/2
  z$xmax <- z$x + l_p[1]/2
  z$ymin <- z$y - l_p[2]/2
  z$ymax <- z$y + l_p[2]/2
  
  # Get the arclength along the center of the pixel in the x and y directions
  z$ly <-geosphere::distGeo(project(as.matrix(z[,.(x,ymin)]),proj,crs('+proj=longlat')),
                            project(as.matrix(z[,.(x,ymax)]),proj,crs('+proj=longlat')))
  
  z$lx <- geosphere::distGeo(project(as.matrix(z[,.(xmin,y)]),proj,crs('+proj=longlat')),
                             project(as.matrix(z[,.(xmax,y)]),proj,crs('+proj=longlat')))
  
  # Calculate pixel size
  z$A <- z$lx + z$ly
  
  # Return distortion raster
  return(rast(z[,.(x,y,lx,ly,A)],crs=proj))
}