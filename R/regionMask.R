#' Function that converts raster and SpatialPolygon* objects
#' to a list of cells that fall within such a region.
#'
#' @title Convert SpatRaster, Raster*, SpatVector, or SpatialPolygon* to "x,y"
#' @param region An object of class SpatRaster, Raster*, SpatVector, 
#' or SpatialPolygon* to convert to "x,y"
#' @param proj A crs object or character string representing the output
#' projection. Default is \code{proj = crs(region)}
#' unless \code{proj} or \code{z_fix} are provided in which case 
#' the latter takes precedence. 
#' @param id A character string indicating which column in a Spatial* or Spat*
#' contains each feature's unique ID. Otherwise ignored
#' @param z_fix A SpatRaster with the same origin and resolution as the
#' \code{z_fix} used to generate the 'world' with \code{\link[lbmech]{makeWorld}}.
#' Do not modify this parameter if you didn't modify it when you ran
#' \code{\link[lbmech]{makeWorld}}.
#' @param precision An integer representing the number of decimals to retain
#' in the x and y directions. For grid sizes with nice, round numbers precisions
#' can be low. This factor is controled by \code{\link[terra]{rast}} and
#' must be the same as the one used to generate the 
#' 'world' with \code{\link[lbmech]{makeWorld}}. Default is 2.
#' @param ... Additional arguments to pass to \code{\link[lbmech]{fix_z}}.
#' @return A character vector containing all cells that fall in the same
#' location as the input 'region'. If \code{id} is provided, a data.table.
#' @importFrom terra resample
#' @importFrom terra cells
#' @importFrom terra xFromCell
#' @importFrom terra yFromCell
#' @importFrom terra rast
#' @importFrom terra crs
#' @importFrom terra crs<-
#' @importFrom terra ext
#' @importFrom terra res
#' @importFrom terra project
#' @importFrom terra rasterize
#' @importFrom data.table data.table
#' @importFrom data.table :=
#' @importFrom data.table setnames
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
#' # Generate a polygon that falls within the DEM's extent
#' region <- ext(c(12500,12600,32500,32700))
#' region <- as.polygons(region)
#' crs(region) <- crs(dem)
#' 
#' maskedCells <- regionMask(region = region, res = res(dem))
#' @export
regionMask <- function(region, proj = crs(region), id = NULL,
                       z_fix = NULL, precision = 2, ...){
  # This bit is to silence the CRAN check warnings for literal column names
  x=y=coord=..regionVals=..regionCells=NULL
  #
  if (is.null(z_fix)){
    z_fix <- fix_z(proj = proj, ...)
  }
  
  z_temp <- as.data.table(
    expand.grid(x = seq(from = ext(region)[1],
                        to = ext(region)[2],
                        by = res(z_fix)[1]),
                y = seq(from = ext(region)[3],
                        to = ext(region)[4],
                        by = res(z_fix)[2])))
  z_temp <- rast(z_temp[,.(x,y,z=1)])
  
  crs(z_temp) <- crs(z_fix)
  z_fix <- suppressWarnings(project(z_temp,z_fix,
                                    align = TRUE))
  
  if (methods::is(region,"Spatial")){
    region <- vect(region)
  }
  if (methods::is(region,"SpatVector")){ 
    if (!is.null(id)){
      region <- rasterize(region,z_fix,field = id)
    } else {
      region <- rasterize(region,z_fix)  
    }
  } else if (methods::is(region, "Raster")){
    region <- rast(region)
  }
  if (methods::is(region, "SpatRaster")){
    region <- resample(region, z_fix)
    regionVals <- region
    names(regionVals) <- id
  } else {
    stop("Inappropriate Region Object. Only Spatial*, Spat*, or Raster allowed")
  }
  regionCells <- cells(region)
  region <- data.table(x = round(xFromCell(z_fix, regionCells),precision),
                       y = round(yFromCell(z_fix, regionCells),precision))
  
  region[, `:=`(coord = paste(format(x, scientific = FALSE),
                              format(y, scientific = FALSE), sep=','), x = NULL, y = NULL)
  ][, coord := stringr::str_remove_all(coord," ")]
  
  if (is.null(id)){
    return(region$coord)
  } else {
    region[,(id) := ..regionVals[..regionCells]]
    region <- region[,.(get(id),coord)]
    setnames(region, c(id,"coord"))
    return(region[])
  }
}