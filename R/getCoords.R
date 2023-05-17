#' Function to get the coordinates in "x,y" format for a given set of points
#'
#' @title Get "x,y" coordinates in appropriate format
#' @param data An object of class data.table or something coercible to it
#' containing the coordinates needing conversion, or a SpatialPointsDataFrame.
#' @param proj A crs object or character string representing the output
#' projection. Required unless \code{z_fix} is provided in which case 
#' \code{proj} is ignored. 
#' @param x A character vector representing the column containing the 'x' coordinates.
#' Required if \code{data} is not SpatialPointsDataFrame.
#' @param y A character vector representing the column containing the 'y' coordinates.
#' Required if \code{data} is not SpatialPointsDataFrame.
#' @param z_fix A SpatRaster with the same origin and resolution as the
#' \code{z_fix} used to generate the 'world' with \code{\link[lbmech]{makeWorld}}.
#' @param precision An integer representing the number of decimals to retain
#' in the x and y directions. For grid sizes with nice, round numbers precisions
#' can be low. This factor is controled by \code{\link[terra]{rast}} and
#' must be the same as the one used to generate the 
#' 'world' with \code{\link[lbmech]{makeWorld}}. Default is 2.
#' @param ... Additional arguments to pass to \code{\link[lbmech]{fix_z}}.
#' @return A vector containing the requested coordinates in appropriate format
#' in the same order as the input data.
#' @importFrom data.table :=
#' @importFrom data.table as.data.table
#' @importFrom terra cellFromXY
#' @importFrom terra xyFromCell
#' @importFrom terra res
#' @importFrom terra rast
#' @importFrom terra crs
#' @importFrom terra crs<-
#' @importFrom terra project
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
#' # Generate five random points that fall within the DEM
#' points <- data.table(x = runif(5, ext(dem)[1], ext(dem)[2]),
#'                      y = runif(5, ext(dem)[3], ext(dem)[4]))
#' 
#' # Get the coordinates
#' points$Cell <- getCoords(points, z_fix = dem)
#' @export
getCoords <- function(data, proj = NULL, x = "x", y = "y", 
                      z_fix = NULL, precision = 2,...){
  # This bit is to silence the CRAN check warnings for literal column names
  ..x=..y=Cell=..poi=NULL
  rm(..x,..y)
  # End
  
  if (is.null(z_fix)){
    z_fix <- fix_z(proj = proj, ...)
  }
  
  data <- as.data.table(data)
  
  if (nrow(data) == 0){
    return("NA")
  }
  
  z_temp <- as.data.table(expand.grid(x = seq(from = min(data[,..x]) - 2 * res(z_fix)[1],
                                              to = max(data[,..x]) + 2 * res(z_fix)[1],
                                              by = res(z_fix)[1]),
                                      y = seq(from = min(data[,..y]) - 2 * res(z_fix)[2],
                                              to = max(data[,..y]) + 2 * res(z_fix)[2],
                                              by = res(z_fix)[2])))
  z_temp <- rast(z_temp[,.(x,y,z=1)])
  crs(z_temp) <- crs(z_fix)
  z_fix <- suppressWarnings(project(z_temp,z_fix, align= TRUE))
  data <- data[,.(x = get(..x), y = get(..y))]
  data <- cellFromXY(z_fix, data)
  data <- as.data.table(format(round(as.data.table(xyFromCell(z_fix, data)),precision), 
                              scientific=FALSE))
  data[,Cell := paste(x,y,sep=',')
      ][, Cell := stringr::str_remove(Cell," ")]
  
  
  
  return(as.character(data$Cell))
}
