#' Get a sample of XY locations from a single raster meeting a certain
#' condition
#' 
#' @title Get a conditional sample from a raster
#' @param z A SpatRaster or RasterLayer object
#' @param n  Numeric. How many samples to draw
#' @param minval Numeric. Select all cells larger than this value
#' @param maxval Numeric. Select all cells lesser than thsi value
#' @param replace Logical. Should samples be taken with replacement?
#' @importFrom data.table as.data.table
#' @importFrom terra rast
#' @importFrom terra cells
#' @importFrom terra xyFromCell
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
#' # Take a sample of 10 values between 150 and 200
#' points <- rangeSample(dem, 10, minval = 150, maxval = 200)
#' @export
rangeSample <-  function(z, n, minval = NULL, maxval = NULL, replace = FALSE){
  if (methods::is(z,'Raster')){
    z <- rast(z)
  }
  if (!is.null(maxval)){
  z[z > maxval] <- NA
  }
  
  if (!is.null(minval)){
    z[z < minval] <- NA
  }
  out <- sample(cells(z),n, replace = replace)
  out <- as.data.table(xyFromCell(z,out))
  return(out)
}
