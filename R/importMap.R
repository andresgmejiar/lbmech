#' Import a raster for a specific region from a multisource environment, 
#' such as the outputs of the getMap function.
#' 
#' @title Import a contiguous raster map
#' @param region A SpatRaster, Raster*, SpatVector, Spatial* or character
#' string representing the area of interest
#' @param polys A polygon of class SpatVector representing the partitioning grid
#' for the maximum possible area, in the same format as the output of the 
#' \code{\link[lbmech]{makeGrid}} function. 
#' @param tile_id A character string representing the name of the column in the 
#' \code{polys} polygon pontaining the unique Tile IDs. Default is 
#' \code{tile_id = 'TILEID'}
#' @param z_fix A SpatRaster or Raster* object with the desired output projection, 
#' resolution, and origin. Required if \code{tiles} is of classes SpatVector,
#' Spatial*, or character, unless \code{res} is provided. 
#' @param neighbor_distance An integer representing the number of cells that 
#' tiles are buffered. In other words, to ensure that there are no gaps between tiles,
#' neighboring tiles within \code{neighborhood_distance} cells are also considered as
#' potential sources. Default is 5 cells.
#' @param FUN Function to deal with overlapping values for overlapping tiles.
#' Default is \code{NA}, which uses \code{\link[terra]{merge}}. To use 
#' \code{\link[terra]{mosaic}}, provide a compatible function
#' @param mask If \code{FALSE} (the default), the output map will contain all cells
#' falling within the extent of \code{region}. If \code{TRUE}, places with \code{NA}
#' (if \code{region} is SpatRaster or Raster*) or no coverage (if \code{region} is
#' SpatVector or Spatial*) will be assigned a value of NA. 
#' @param vals A character string or a Raster* object. Required only if the
#' \code{z} parameter is a polygon NOT the output of the 
#' \code{\link[lbmech]{makeGrid}} function as the default is
#' the character string \code{'location'}. If not, the \code{vals} parameter should be
#' set to the column name containing the URL or file path to the DEM for that
#' sector.
#' @param filt Numeric. Size of moving window to apply a low-pass filter. Default 
#' is \code{filt = 0}. Ignored unless the tiles need to be generated from
#' the raw source files. 
#' @param dir A filepath to the directory being used as the workspace. Default
#' is \code{tempdir()}, but unless the analyses will only be performed a few times 
#' it is highly recommended to define a permanent workspace.
#' @param ... Additional agruments to pass to \code{\link[lbmech]{fix_z}}
#' @importFrom terra rast
#' @importFrom terra vect
#' @importFrom terra crs
#' @importFrom terra ext
#' @importFrom terra buffer
#' @importFrom terra intersect
#' @importFrom terra merge
#' @importFrom terra mosaic
#' @importFrom terra as.polygons
#' @importFrom terra project
#' @importFrom terra crs<-
#' @importFrom terra crop
#' @importFrom terra mask
#' @importFrom terra aggregate
#' @importFrom terra res
#' @examples 
#' # Generate a DEM, export it
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
#' dir <- tempdir()
#' writeRaster(dem, paste0(dir,"/DEM.tif"),overwrite=TRUE)
#' 
#' 
#' # Import raster, get the grid
#' dem <- rast(paste0(dir,"/DEM.tif"))
#' grid <- makeGrid(dem = dem, nx = n, ny = n, sources = TRUE)
#' 
#' # Import the map for the center tile resampled to a different resolution
#' dem2 <- importMap('SECTOR_13', grid, res = 20)
#' @export
importMap <- function(region, polys,
                      tile_id = 'TILEID', z_fix = NULL,
                      neighbor_distance = 5, FUN = NULL, mask = FALSE,
                      vals = 'location', filt = 0, 
                      dir = tempdir(), ...){
  dir <- normalizePath(dir,mustWork=FALSE)
  
  # From the region we need the tiles that it covers, and it should be coerced
  # to something of SpatRaster or SpatVector class
  
  #Try to coerce
  if (methods::is(region,'character')){
    region <- polys[which(unlist(polys[[tile_id]]) %in% region),]
  } else if (methods::is(region,'Raster')){
    region <- rast(region)
  } else if (methods::is(region,'Spatial')){
    region <- vect(region)
  }
  if (methods::is(z_fix,'Raster')){
    z_fix <- rast(z_fix)
  }
  
  # Deal with terra objects
  if (is.null(z_fix) & !methods::is(region,'SpatRaster')){
    if (!methods::is(region,'character')){
      z_fix <- fix_z(crs(region), ...)
    } else {
      z_fix <- fix_z(crs(polys), ...)
    }
  }
  if (methods::is(region,'SpatRaster')){
    if (is.null(z_fix)){
      z_fix <- region
    }
  } 
  if ('makeGrid' %in% names(polys)){
    polys <- aggregate(polys,by=c(tile_id,vals,'makeGrid'))
  } else {
    polys <- aggregate(polys,by=c(tile_id,vals))
  }
  neighbor_distance <- neighbor_distance*max(res(z_fix))
  
  # Get the tile names, make sure maps have been downlaoded
  tiles <- buffer(region,width=neighbor_distance)
  tiles <- intersect(region[,NA],polys)[[tile_id]]
  tiles <- unlist(tiles)
  getMap(tiles, polys, tile_id = tile_id, dir = dir, filt = filt)
  
  # Import all DEMs
  dem <- normalizePath(paste0(dir,'/',tiles),mustWork=FALSE)
  dem <- lapply(dem,importRST)
  dem <- suppressWarnings(lapply(dem,
                                 project,
                                 y = dem[[1]],
                                 align = TRUE))
  
  # Merge or mosaic them, depending on whether a FUN is provided 
  if (length(dem) > 1){
    if (is.null(FUN)){
      dem <- do.call(merge,dem)
    } else {
      dem$fun <- FUN
      dem <- do.call(mosaic,dem)
    }
  } else {
    dem <- dem[[1]]
  }
  
  # Create a cropping polygon. This is equivalent to the extent of the
  # region shape, expanded in the X and Y directions by the
  # pixel size times the connectivity contiguity. This has to be
  # reprojected back to the original coordinate system since only it
  # preserves orthogonality. This step saves us from having to initially project
  # a rather large SpatRaster
  if (methods::is(region,'SpatRaster')){
    poly <- region
    poly[poly == poly] <- 1
    poly <- as.polygons(poly)
  } else {
    poly <- region
  }
  
  poly <- project(poly,dem)
  poly <- unlist(ext(poly))
  poly <- c(poly[1] - (neighbor_distance + 1) * res(z_fix)[1], poly[2] + 
              (neighbor_distance + 1) * res(z_fix)[1],
            poly[3] - (neighbor_distance + 1) * res(z_fix)[2], poly[4] + 
              (neighbor_distance + 1) * res(z_fix)[2])
  poly <- ext(poly)
  poly <- as.polygons(poly)
  crs(poly) <- crs(dem)
  
  # Crop the mosaiced raster by the cropping polygon, project it
  dem <- crop(dem,poly)
  dem <- suppressWarnings(project(dem,z_fix,align = TRUE))
  
  # Should the output raster have the same extent, or the same coverage?
  if (mask == TRUE){
    dem <- mask(dem,region)
  }
  
  if (crs(z_fix) != crs(region)){
    warning("Different projections provided for 'z_fix' and 'region'. Output is the same as 'z_fix'")
  }
  return(dem)
}