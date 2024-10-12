#' Function that defines the grid that can be traversed--the "world"--as well as the
#' cells that can be accessed from each individual cell. This is the most
#' time-intensive function. 
#' 
#' It first checks to see if the required 
#' files contained within a \code{'/World/'} directory have already been created 
#' in the \code{dir} workspace. If not it creates them; if the \code{'/World/'}
#' directory exists though, then an error is thrown (default) OR the user has the option
#' to overwrite the directory. 
#' 
#' The default parameters are sufficient for a workflow involving calculating
#' costs with the \code{\link[lbmech]{calculateCosts}} function.
#' @title Define the topographic landscape
#' @param source An object of class SpatVector, or potential inputs to
#' \code{\link[lbmech]{makeGrid}} that define the paths to the original source files 
#' of the desired rasters in the same format as the
#' output of the \code{\link[lbmech]{makeGrid}(sources = TRUE)}. 
#' If the object is NOT a already a polygon of the appropriate format, set 
#' \code{source_id = NULL} and add the necessary \code{\link[lbmech]{makeGrid}}
#' parameter
#' @param source_id A character string representing the name of the column in the
#' \code{source} polygon containing the unique Tile IDs. Default is \code{source_id = 'TILEID'}
#' @param grid An object of class SpatVector representing
#' the partitioning grid for the maximum possible travel area. Smaller grids increase 
#' the amount of area that can be read into the memory, but require more I/O operations.
#' Default is \code{grid = source}.
#' @param grid_id A character string representing the name of the column
#' in the \code{grid} polygon containing the unique Tile IDs. Default is \code{tile_id = 'TILEID'}
#' @param cut_slope A number representing the dimensionless maximum slope
#' of ascent/descent. To ignore, set \code{cut_slope = Inf}.
#' @param proj A crs object or character string representing the output projection.
#' Default projection is \code{proj = crs(polys)} unless a `z_fix` or `proj` is 
#' provided, in which case the latter is ignored. Great care should be 
#' employed to ensure that the projection is conformal and in meters. 
#' @param directions One of the integers \code{c(4, 8, 16)},
#' the character string \code{'bishop'},or a neighborhood matrix.
#' Default is \code{directions = 16}, implying that all 'knight and one-cell
#' queen moves' are permissible movements on the grid. See \code{\link[raster]{adjacent}}.
#' @param neighbor_distance An integer representing the distance in meters
#' that tiles are buffered. In other words, to ensure that all transitions in the
#' 'world' are recorded, files for each tile will contain a number of observations
#' that fall outside of the tile in other ones. Default is 100 m, but adjust
#' on raster size.
#' @param z_fix A SpatRaster or Raster* that will define the resolution, origin, and
#' projection information for the entire "world" of possible movement. Note that
#' it does NOT need the same extent. Default resolution is 5, and offset is 0.
#' Default projection is \code{proj = crs(polys)} unless a `z_fix` or `proj` is 
#' provided, in which case the latter is ignored. Great care should be 
#' employed to ensure that the projection is conformal and in meters. 
#' @param z_min The minimum allowable elevation. Useful if DEM source includes
#' ocean bathymetry as does the SRTM data from AWS. Default is \code{z_min = NULL},
#' but set to \code{0} for SRTM data.
#' @param unit One of \code{c("m", "km", "ft", "mi")}, representing the unit of the DEM.
#' All will be converted to meters, which is the default.
#' @param vals A character string or a SpatRaster or Raster* object. Ignored unless the
#' \code{polys} parameter is a polygon NOT the output of the \code{\link[lbmech]{makeGrid}}
#' function as the default is the character string \code{'location'},
#' AND the appropriate world \code{.gz} file is NOT
#' present in the workspace directory. In which case it must represent either the
#' original DEM or a character string with the column representing the DEM
#' filepath or URL.
#' @param precision An integer representing the number of decimals to retain
#' in the x and y directions. For grid sizes with nice, round numbers precisions
#' can be low. This factor is controled by \code{\link[terra]{rast}}.
#' Default is 2.
#' @param dist A character string representing the way distances should be
#' calculated. Default is \code{dist = 'proj'} for the default projection units,
#' but the following geodesic methods are also available: 
#' \code{c('haversine','cosine','karney','meeus','vincentyEllipsoid','vincentySphere')}.
#' The developer recommends \code{'karney'}.
#' @param r The earth's radius. Employed only if one of the geodesic methods is used.
#' Default is \code{r = 6378137} for WGS1984.
#' @param f The earth's ellipsoidal flattening. Employed only if 
#' \code{dist \%in\% c('karney','meeus','vincentyEllipsoid')}. Default is 
#' \code{f=1/298.257223563} for WGS1984
#' @param b The earth's semiminor axis. Employed only if \code{dist = 'vincentyEllipsoid'}.
#' Default is \code{b=6356752.3142} for WGS1984.
#' @param FUN Function to deal with overlapping values for overlapping sectors.
#' Default is NA, which uses \code{\link[terra]{merge}}. 
#' To use \code{\link[terra]{mosaic}}, provide a compatible function.
#' @param sampling How to resample rasters. Default is \code{'bilinear'} interpolation,
#' although \code{'ngb'} nearest neighbor is available for categorical rasters. 
#' @param dir A filepath to the directory being used as the workspace.
#' Default is \code{tempdir()} but unless the analyses will only be performed a few
#' times it is highly recommended to define a permanent workspace.
#' @param overwrite If a directory with a \code{World} subdirectory already exists,
#' should the latter be overwritten? Default is \code{overwrite = FALSE}.
#' @param filt Numeric. Size of moving window to apply a low-pass filter. Default 
#' is \code{filt = 0}. Ignored unless the tiles need to be generated from
#' the raw source files. 
#' @param water Optional. One of (1) SpatialPolygon* or SpatVector polygon representing 
#' the area covered by water, in which case \code{water_speed} must also be provided;
#' (2) A RasterLayer or single-layer SpatRaster representing the speed of water at a
#' given location. Flow direction will be calculated from the world's DEM; or (3) A 
#' two layer SpatRaster or Raster* object, representing either the horizontal/vertical
#' components of water velocity (in m/s) or the absolute water speed and flow direction
#' (in that order; see \code{uv}).
#' @param water_speed A character representing the column name in 
#' \code{water} that contains the average speed of the water in m/s. Required
#' if \code{water} is a polygon.
#' @param uv Logical. If \code{TRUE} (the default), a two-layer raster input for
#' \code{water} is taken to be the horizontal and vertical components of water velocity. 
#' If false, a two-layer raster input for \code{water} is taken to be the water speed
#' and flow direction.
#' @param priority One of \code{'water'} (the default), or \code{'land'}, indicating
#' whether land values or water values should take precedence when values for a given
#' location exist in both datasets.  
#' @param cols A character vector. Default is c("x_i","y_i","dz","dl","dr"), and
#' dictates which columns are returned but final values for x and y and final/initial
#' z values are also available
#' @param ... Additional arguments to pass to \code{\link[lbmech]{fix_z}}.
#' @return A \code{/World/} directory, containing \code{/Loc/}, \code{/Diff/},
#' and \code{/Raw/}, directories where cropped and transformed files will be stored,
#' and \code{/callVars.gz/}, \code{/z_sources/}, \code{/z_grid/}, 
#' \code{/z_fix.fst/}, and \code{/z_fix.fstproj/} files defining the world.
#' @importFrom terra crs
#' @importFrom terra values
#' @importFrom terra intersect
#' @importFrom terra values<-
#' @importFrom terra writeVector
#' @importFrom terra nlyr
#' @importFrom data.table as.data.table
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
#' # Make a world but limit it to the DEM grid size
#' defineWorld(source = grid, cut_slope = 0.5, 
#'             res = res(dem), dir = dir, overwrite=TRUE)
#' @export 
defineWorld <- function(source, source_id = 'TILEID', z_min = NULL,
                        grid = NULL, proj = NULL, grid_id = 'TILEID', cut_slope = Inf,
                        directions = 16, neighbor_distance = 10, 
                        z_fix = NULL,unit = "m", vals = 'location', precision = 2, 
                        dist = 'proj', r = 6378137, f = 1/298.257223563, b = 6356752.3142,
                        FUN = NULL, sampling = 'bilinear', 
                        water = FALSE, water_speed = 'speed', uv = TRUE,
                        priority = 'water',
                        overwrite = FALSE,
                        filt = 0,
                        dir = tempdir(),
                        cols = c("x_i","y_i","dz","dl","dr"), ...){
  ..source_id=location=..vals=..grid_id=NULL
  # Source needs to be in the form of makeGrid-type object
  if (is.null(source_id)){
    source <- makeGrid(source,sources = TRUE,...)
    grid <- source
    source_id <- 'TILEID'
    proj <- crs(grid)
  } 
  if (is.null(grid)){
    grid <- source
  } 
  
  if (is.null(proj)){
    proj <- crs(grid)
  }
  
  
  # Normalize input grid names
  if (!("makeGrid" %in% names(source))){
    values(source) <- as.data.table(values(source))[,.(source_id = get(..source_id),
                                                       location = get(..vals))]
  } else {
    values(source) <- as.data.table(values(source))[,.(source_id = get(..source_id),
                                                       location = get(..vals),
                                                       makeGrid = TRUE)] 
  }
  
  values(grid) <- as.data.table(values(grid))[,.(id = get(..grid_id))]
  
  
  # Create the subdirectories
  dir <- normalizePath(paste0(dir,"/World"), mustWork = FALSE)
  subdirs <- c("/Raw","/Local","/Diff")
  if (overwrite != TRUE & dir.exists(dir)){
    warning(paste(dir,"already exists, inputs ignored.\nTo overwrite, use 'overwrite=TRUE'. Carefuly ensure that this is the desired behavior."))
    return()
  } else if (overwrite == TRUE & dir.exists(dir)){
    unlink(dir,recursive=TRUE)
  }
  dir.create(dir)
  writeVector(source,normalizePath(paste0(dir,"/z_sources.gpkg"),mustWork = FALSE))
  writeVector(grid,normalizePath(paste0(dir,"/z_grid.gpkg"),mustWork = FALSE))
  
  subdirs <- normalizePath(paste0(dir,subdirs),mustWork=FALSE)
  for (i in subdirs){
    dir.create(i)
  }
  
  units <- c("m","km","ft","mi")
  if (!(unit %in% units)){
    stop("Error: Unknown units. Only m, km, ft, and mi allowed.
         If your data source is in another unit, please convert units
         before running the getMap function.")
  }
  
  if (directions == 4 | directions == 8 | directions == 'bishop'){
    contiguity <- 1
  } else if (directions == 16){
    contiguity <- 2
  } else{
    contiguity <- base::ceiling(max(nrow(directions),ncol(directions)))
  }

  # If no z_fix is provided, make one
  if (is.null(z_fix)){
    z_fix <- fix_z(proj = proj, ...)
  } else if (methods::is(z_fix, "Raster")){
    z_fix <- rast(z_fix)
  }
  writeRST(z_fix,normalizePath(paste0(dir,"/z_fix"),mustWork = FALSE))
  neighbor_distance <- neighbor_distance*max(res(z_fix))
  
  callVars <- data.table(cut_slope = cut_slope,
                         source_id = source_id,
                         z_min = z_min,
                         grid_id =  grid_id,
                         directions = list(directions),
                         neighbor_distance = neighbor_distance,
                         l_p = list(res(z_fix)),
                         unit = unit,
                         vals = vals,
                         precision = precision,
                         FUN = list(FUN),
                         sampling = sampling,
                         dist = dist,
                         filt = filt,
                         r = r,
                         f = f,
                         b = b,
                         water_speed = water_speed,
                         priority = priority,
                         uv = uv)
  
  if (!methods::is(water,'logical')){
  if (methods::is(water,'Spatial')){
    water <- vect(water)
  } else if (methods::is(water,"Raster")){
    water <- rast(water)
  }
  if (methods::is(water,'SpatVector')){
    if (methods::is(water,'numeric')){
      water <- water[,NA]
      water$speed <- water_speed
    } else if (methods::is(water,'character')){
      water$speed <- water[,water_speed]
      water <- water[,'speed']
    }
    writeVector(project(water[,'speed'],crs(z_fix)),
                normalizePath(paste0(dir,"/water.gpkg"),mustWork = FALSE))
  } else if (methods::is(water,'SpatRaster')){
    if (nlyr(water) > 2){
      stop("Water may not have more than two layers.")
    }
    writeRaster(project(water, crs(z_fix)),
                normalizePath(paste0(dir,"/water.tif"),mustWork = FALSE))
  }
  }
  saveRDS(callVars,normalizePath(paste0(dir,"/callVars.gz"),mustWork = FALSE))
}
