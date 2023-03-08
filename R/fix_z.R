#' Create a raster that can be used to define
#' the resolution, origin, and projection to be 
#' employed for all least-cost analyses. If a source
#' DEM has such properties you may use that.
#' 
#' @title Define the sampling grid
#' @param proj A \code{\link[raster]{crs}} object or character string containing
#' projection information. Should be conformal and in meters.
#' @param res A numeric of length one or two nrepresenting the spatial resolution.
#' Default is 5. 
#' @param dx The horizontal offset from the origin (see \code{\link[raster]{origin}}).
#' Default is 0 (this does not correspond to an origin of zero however).
#' @param dy The vertical offset from the origin (see \code{\link[raster]{origin}}).
#' Default is 0 (this does not correspond to an origin of zero however).
#' @return A SpatRaster object consisting of four cells, with resolution \code{res} and
#'  the origin at \code{x = nx} and \code{y = ny}.
#' @importFrom terra rast
#' @importFrom data.table data.table
#' @examples 
#' projection <- "+proj=lcc +lat_1=48 +lat_2=33 +lon_0=-100 +datum=WGS84"
#' z_fix <- fix_z(res = 2, proj = projection)
#' @export
fix_z <- function(proj, res = 5, dx = 0, dy = 0){
  if (length(res) == 1){
    res <- c(res,res)
  }
  z <- data.table(x = c(0,0,res[1],res[1]) + dx,
                  y = c(0,res[2],0,res[2]) + dy,
                  z = c(1,1,1,1))
  z <- rast(z, crs=proj)
  return(z)
}