#' Append x,y data to a raster (wrapper for \code{\link[terra]{init}}).
#' 
#' @title Append x,y data to a raster 
#' @param z An object of class Raster* or SpatRaster
#' @return A data.table containing the raster values in column names after
#' the raster layers, and 'x' and 'y' columns containing the cell locations
#' @importFrom terra init
#' @importFrom terra rast
#' @importFrom data.table as.data.table
#' @examples 
#' n <- 5
#' dem <- expand.grid(list(x = 1:(n * 100),
#'                         y = 1:(n * 100))) / 100
#' dem <- as.data.table(dem)
#' dem[, z := 250 * exp(-(x - n/2)^2) + 
#'       250 * exp(-(y - n/2)^2)]
#' dem <- rast(dem)
#' ext(dem) <- c(10000, 20000, 30000, 40000)
#' crs(dem) <- "+proj=lcc +lat_1=48 +lat_2=33 +lon_0=-100 +datum=WGS84"
#' dem <- rastToTable(dem)
#' @export
rastToTable <- function(z){
  if (methods::is(z,"Raster")){
    z <- rast(z)
  }
  if ("x" %in% names(z)){
    warning("Object contains a layer named 'x'. Did not calculate easting")
  } else {
    z$x <- init(z,'x')
  }
  if ("y" %in% names(z)){
    warning("Object contains a layer named 'y'. Did not calculate easting")
  } else {
    z$y <- init(z,'y')
  }
  noms <- names(z)
  stats::na.omit(as.data.table(z))
  if (nrow(z) == 0){
  z <- rbind(z, t(rep(NA,length(z))),use.names=FALSE)
  }
  return(z)
}
