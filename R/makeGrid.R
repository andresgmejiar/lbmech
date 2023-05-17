#' Generate a partitioning grid for a single raster source representing regional
#' elevations. Smaller partitioning grids (i.e. a greater value of \code{nx * ny})
#' results in a greater number of saved files and a greater number of
#' read-write operations in future operations, but reduces the amount of
#' memory employed.
#'
#' @title Make partitioning grid
#' @param dem One of either a single character string, SpatVector (polygon),
#' SpatialPolygons*, SpatRaster, or Raster
#' 
#' 
#' If character, SpatRaster, Raster, such an object or (filepath to one) 
#' containing the elevations for the maximum possible extent imaginable for a study. 
#' Note that SpatRaster and Raster only work
#' rasters that have been read in, not those that exist exclusively in the memory. 
#' If you  have just generated the raster and it is in memory, export it first with
#' \code{\link[terra]{writeRaster}} then use the filepath string as \code{dem} or
#' re-import it with \code{\link[terra]{rast}} before using the SoatRaster object.
#' 
#' If SpatVector or SpatialPolygons*, the extent of possible movement.
#' @param nx The integer-number of columns in the output grid
#' @param ny The integer-number of rows in the output grid
#' @param path (Optional) The filepath or URL of the source DEM. Ignored if 
#' \code{dem} is of class raster or SpatRaster. If a SpatVector or SpatialPolygons
#' is provided but no path, \code{\link[lbmech]{getMap}} will use
#' \code{\link[elevatr]{get_elev_raster}} to download topographic data.
#' @param sources Logical. Should source information be saved as attributes
#' to the grid for use in \code{\link[lbmech]{getMap}} and 
#' \code{\link[lbmech]{defineWorld}}? Default is \code{sources = FALSE}.
#' @param proj A \code{\link[terra]{crs}} or something coercible to it representing 
#' the desired output projection. Default is the input raster's projection.
#' @param prefix A character string containing the prefix to name individual sectors.
#' Default is \code{prefix = "SECTOR_"}
#' @param crop Logical. If TRUE (the default), the output polygons will be cropped
#' by the original \code{dem} (if SpatVector or SpatialPolygons*), or by
#' the area covered by non-NA cells (if raster or SpatRaster). 
#' @param overlap How much should adjacent polygons overlap to ensure there's
#' contiguity between different tiles? Default is \code{overlap = 0.005}.
#' @param zoom Considered only if \code{var = 'z'} and no data source is set. 
#' The zoom level to be downloaded. See documentation for the \code{z} parameter
#' in \code{\link[elevatr]{get_elev_raster}} for further information.
#' Default is 13, but see documentation for the \code{z} parameter in 
#' \code{\link[elevatr]{get_elev_raster}}.
#' @param var If the polygons point to a data source, what will be the variable
#' name in the internal GIS? Default is 'z' for elevation. 
#' @param extension A character vector representing the extension of the source
#' path. Required only if \code{sources = TRUE} and the extension is not apparent
#' from the URL stored in the \code{var} column.
#' @return Polygons of class SpatVector representing the individual sectors ('tiles'),
#' with a dataframe containing three columns: the "TILEID", the raster's filepath,
#' and a dummy column indicating that the grid was made using the makeGrid function.
#' This will be necessary for future functions. The object MUST be stored 
#' on the disk, it should not be stored in the memory
#' @importFrom terra ext
#' @importFrom terra rast
#' @importFrom terra vect
#' @importFrom terra crs
#' @importFrom terra crs<-
#' @importFrom terra intersect
#' @importFrom terra as.polygons
#' @importFrom terra sources
#' @importFrom data.table as.data.table
#' @importFrom data.table :=
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
#' @export
makeGrid <- function(dem, nx, ny, path = NA, sources = FALSE, extension = NULL,
                     proj = NULL, prefix = "SECTOR_", crop = TRUE,
                     zoom = 13, var = 'z', overlap = 0.005){
  #Silences CRAN check warnings
  x=..dx=y=..dy=left=right=bottom=top=NULL
  # Convert raster to SpatRast if raster is provided
  if (methods::is(dem, 'Raster') | methods::is(dem, 'character')){
    dem <- rast(dem)
  } else if (methods::is(dem,"Spatial")){
    dem <- vect(dem)
  } 
  
  # Save the filepath if it's a raster object
  if (methods::is(dem,'SpatRaster')){
    path <- normalizePath(sources(dem))
  }
  
  # Save original covered area if the polygons will be cropped
  if (crop){
    if (methods::is(dem,'SpatRaster')){
      cropMask <- dem
      cropMask[!is.na(cropMask)] <- 1
      cropMask <- as.polygons(cropMask)}
    else if (methods::is(dem,'SpatVector')){
      cropMask <- dem
    }
  }
  
  # If a different projection is desired, transform it
  if (!is.null(proj)){
    dem <- ext(dem)
    dem <- as.polygons(dem)
    dem <- project(dem,proj)
  }
  
  dx     <- (ext(dem)[2]- ext(dem)[1])/ nx  ## extent of one tile in x
  dy     <- (ext(dem)[4]- ext(dem)[3])/ ny  ## extent of one tile in y
  xs     <- seq(ext(dem)[1], by= dx, length= nx) ## lower left x-coordinate
  ys     <- seq(ext(dem)[3], by= dy, length= ny) ## lower left y-coordinate
  cS     <- expand.grid(x= xs, y= ys)
  
  # What percent overlap should we consider to ensure everything is within the 
  # area even if there's a funky projection
  dx <- dx * (1 + overlap)
  dy <- dy * (1 + overlap)
  
  # Get the eastings/northings of the extent limits
  cS <- as.data.table(cS)
  cS[, `:=`(left = x,
            right = x + ..dx,
            bottom = y,
            top = y + ..dy)]
  
  # Convert the data.table to polygons and assign it a unique TILEID
  polys <- apply(cS[,.(left,right,bottom,top)],
                 1, FUN = function(x) as.polygons(ext(x)))
  if (nx * ny > 1){
    polys <- do.call("rbind",polys)
  } else {
    polys <- polys[[1]]
  }
  polys$TILEID <- paste0(prefix,stringr::str_pad(1:length(polys),
                                                 nchar(length(polys)),
                                                 side='left',
                                                 pad='0'))
  
  # Apply the source projection
  crs(polys) <- crs(dem)
  
  if (crop){
    polys <- intersect(polys,
                       project(cropMask[,NA],polys))
    polys <- aggregate(polys,by='TILEID')
  }
  
  if (sources){
    # Set the location of the data to the filepath of the original raster
    # or the desired location. If it's NA, put the zoom value
    if (is.na(path)){
      path <- zoom
    }
    polys$location <- path
    if (!is.null(extension)){
      polys$extension <- extension
    }
  }
  # If the output vector is input into getMap, this tells getMap that the source
  # DEM is a single file and it should crop; not download.
  polys$makeGrid <- TRUE
  
  return(polys)
}
