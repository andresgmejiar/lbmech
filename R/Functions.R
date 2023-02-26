#' Convert a raster to a data.table keeping the geometry properties. Note that
#' the raster may not have any layers named 'x' or 'y' (lowercase). This does
#' not preserve projection information.
#' 
#' @title Convert raster to data.table
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
  return(na.omit(as.data.table(z)))
}

#' Read and write rasters using the fst library
#' 
#' @title Fast read and write rasters
#' @name writeRST
#' @aliases importRST
#' @param x An object of class SpatRaster or Raster. It may not contain
#' layers named 'x' or 'y'
#' @param filename Character. Output filename. Do not use extensions.
#' @param ... Additional parameters to pass on to \code{\link[fst]{read_fst}} or
#' \code{\link[fst]{write_fst}}
#' @importFrom terra crs
#' @importFrom terra crs<-
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
#' writeRST(dem, 'DEM.fst')
#' importRST('DEM.fst')
NULL

#' @rdname writeRST
#' @export
writeRST <- function(x, filename, ...){
  filename <- normalizePath(filename, mustWork = FALSE)
  if ('x' %in% names(x) | 'y' %in% names(x)){
    stop("Raster may not contain layers named 'x' or 'y'. Operation aborted")
  }
  xres <- res(x)[1]
  yres <- res(x)[2]
  xoff <- round(as.numeric(xFromCell(x,1))/xres) * 
    xres - as.numeric(xFromCell(x,1))
  yoff <- round(as.numeric(yFromCell(x,1))/yres) * 
    yres - as.numeric(yFromCell(x,1))
  
  fst::write_fst(rbind(data.table(x = c(xres,xoff),
                                  y = c(yres,yoff)),
                       rastToTable(x), fill = TRUE), 
                 paste0(filename,'.fst'), compress = 0, ...)
  write(crs(x), paste0(filename,'.fstproj'))
}

#' @rdname writeRST
#' @export
importRST <- function(filename, ...){
  filename <- normalizePath(filename, mustWork = FALSE)
  x <- fst::read_fst(paste0(filename,'.fst'),as.data.table = TRUE)
  xres <- x[1]$x
  yres <- x[1]$y
  xoff <- x[2]$x
  yoff <- x[2]$y
  x[, `:=`(x = round(x/..xres)*..xres + xoff,
           y = round(y/..yres)*..yres + yoff)]
  x <- rast(x[!c(1,2)],
            crs = as.character(
              noquote(
                paste(
                  readLines(
                    paste0(filename,'.fstproj')),collapse='\n')))
  )
}

#' Generate a partitioning grid for a single raster source representing regional
#' elevations. Smaller partitioning grids (i.e. a greater value of \code{nx * ny})
#' results in a greater number of saved files and a greater number of
#' read-write operations in future operations, but reduces the amount of
#' memory employed.
#'
#' @title Make partitioning grid
#' @param dem One of either a SpatRaster, Raster, SpatVector (polygon),
#' or SpatialPolygons*. 
#' 
#' If SpatRaster or Raster, such an object containing the elevations for the 
#' maximum possible extent imaginable for a study. Note that this only works for 
#' rasters that have been read in, not those that exist exclusively in the memory. 
#' If you  have just generated the raster and it is in memory, export it first with
#' \code{\link[terra]{writeRaster}} then re-import it with \code{\link[terra]{rast}}.
#' 
#' If SpatVector or SpatialPolygons*, the extent of possible movement.
#' @param nx The integer-number of columns in the output grid
#' @param ny The integer-number of rows in the output grid
#' @param path (Optional) The filepath or URL of the source DEM. Ignored if 
#' \code{dem} is of class raster or SpatRaster. If a SpatVector or SpatialPolygons
#' is provided but no path, \code{\link[lbmech]{getMap}} will use
#' \code{\link[elevatr]{get_elev_raster}} to download topographic data. 
#' @param proj A \code{\link[terra]{crs}} or something coercible to it representing 
#' the desired output projection. Default is the input raster's projection.
#' @param prefix A character string containing the prefix to name individual sectors.
#' Default is \code{prefix = "SECTOR_"}
#' @param crop Logical. If TRUE (the default), the output polygons will be cropped
#' by the original \code{dem} (if SpatVector or SpatialPolygons*), or by
#' the area covered by non-NA cells (if raster or SpatRaster). 
#' @param overlap How much should adjacent polygons overlap to ensure there's
#' contiguity between different tiles? Default is \code{overlap = 0.005}.
#' @param zoom Considered only no data source is set. The zoom level to be downloaded. 
#' Default is 13, but see documentation for the \code{z} parameter in 
#' \code{\link[elevatr]{get_elev_raster}}.
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
#' grid <- makeGrid(dem = dem, nx = n, ny = n)
#' 
#' @export
makeGrid <- function(dem, nx, ny, path = NA, 
                     proj = NULL, prefix = "SECTOR_", crop = TRUE,
                     zoom = 13, overlap = 0.005){
  #Silences CRAN check warnings
  x=..dx=y=..dy=left=right=bottom=top=NULL
  # Convert raster to SpatRast if raster is provided
  if (methods::is(dem, 'Raster')){
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
  
  # Set the location of the data to the filepath of the original raster
  # or the desired location. If it's NA, put the zoom value
  if (is.na(path)){
    path <- zoom
  }
  polys$location <- path
  
  # If the output vector is input into getMap, this tells getMap that the source
  # DEM is a single file and it should crop; not download.
  polys$makeGrid <- TRUE
  
  if (crop){
    polys <- intersect(polys,
                       project(cropMask[,NA],polys))
  }
  
  return(polys)
}

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
#' grid <- makeGrid(dem = dem, nx = n, ny = n)
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

#' A function that checks if the DEMs for a given set of sectors exist
#' in the workspace, and if not downloads them or crops them from a larger file
#'
#' @title Download or crop necessary DEMs
#' @param tiles A character vector--such as the output to
#' \code{\link[lbmech]{whichTiles}}---containing the unique tile IDs for sectors that
#' should be in the workspace.
#' @param polys A polygon of class SpatVector representing
#' the partitioning grid for the maximum possible area, in the same format as the
#' output of the \code{\link[lbmech]{makeGrid}} function.
#' @param tile_id a character string representing the name of the column
#' in the \code{polys} polygon containing the unique Tile IDs. Default is \code{tile_id = 'TILEID'}
#' @param vals A character string or a SpatRast or Raster* object. Optional if the
#' \code{polys} polygon is the output of the \code{\link[lbmech]{makeGrid}} function as the default is
#' the character string \code{'location'}. If no DEM was provided when \code{\link[lbmech]{makeGrid}} 
#' was initially run (i.e. polys$location == NA), then the function will 
#' use \code{\link[elevatr]{get_elev_raster}} to download data. 
#' If not from \code{\link[lbmech]{makeGrid}}, the \code{vals} parameter should be
#' set to the column name containing the URL or filepath to the DEM for that
#' sector. 
#' @param z_min The minimum allowable elevation. Useful if DEM source includes
#' ocean bathymetry as does the SRTM data from AWS. Default is \code{z_min = 0},
#' but set to \code{-Inf} to disable.
#' @param filt Should a lowpass filter be applied? Default is a moving window
#' size \code{filt = 3}, but set to \code{NA} or 0 to disable. 
#' @param dir A filepath to the directory being used as the workspace.
#' Default is \code{tempdir()} but unless the analyses will only be performed a few
#' times it is highly recommended to define a permanent workspace.
#' @return Function does not return any objects, but sets up the workspace
#' such that the necessary DEM files are downloaded/cropped and accessible.
#' @importFrom terra rast
#' @importFrom terra crop
#' @importFrom terra focal
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
#' grid <- makeGrid(dem = dem, nx = n, ny = n)
#' 
#' 
#' # Generate five random points that fall within the grid
#' points <- data.table(x = runif(5, ext(dem)[1], ext(dem)[2]),
#'                      y = runif(5, ext(dem)[3], ext(dem)[4]))
#'                
#'                            
#' # Run whichTiles and getMap to prepare appropriate sector files
#' tile_list <- whichTiles(region = points, polys = grid) 
#' getMap(tiles = tile_list, polys = grid, dir = dir)
#' @export
getMap <- function(tiles, polys, tile_id = "TILEID", vals = "location", 
                   z_min = 0, filt = 3,
                   dir = tempdir()){
  # The function first checks to see if the needed "tiles" sector DEMs exist
  # If not, it downloads the correct sectors or crops them from the master
  # DEM
  dir <- normalizePath(dir)
  # Check to see if elevations directory exists. If not, make it
  rd <- normalizePath(dir,mustWork=FALSE)
  
  if(!dir.exists(rd)){
    dir.create(rd)
  }
  
  # Check to see if the file exists; if not add to "down" vector
  down <- c()
  
  for (tile_name in tiles){
    if (!file.exists(normalizePath(paste0(rd,"/",tile_name,".fst"),
                                   mustWork=FALSE))){
      down <- c(down,which(polys[[tile_id]] == tile_name))
    }
  }
  
  # Skip the rest if we have all the needed tiles
  if (!is.null(down)){
    
    # If we don't provide a raster layer and polys doesn't come from makeGrid
    if (!((methods::is(vals,"Raster")) | (methods::is(vals,"SpatRaster"))) &
        is.null(polys$makeGrid)){
      # If what we provide is a polygon with URLs to the source,
      # download the file
      for (i in seq(1,length(down))){
        tile_name <- as.character(polys[down[i],][[tile_id]])
        print(paste0("Downloading Tile ",tile_name," (",
                     i," of ",length(down),")"))
        extension <- stringr::str_split_fixed(polys[down[i],][[vals]],"/.")[2]
        file_path <- normalizePath(paste0(rd,"/",tile_name),mustWork = FALSE)
        utils::download.file(url=polys[down[i],][[vals]],
                             destfile=paste0(file_path,".",extension))
        
        
        # If downloaded file is a compressed file
        exts <- c(".(grd)|(asc)|(sdat)|(rst)|(nc)|(tif)|(envi)|(bil)|(img)|(adf)$")
        if (extension == 'zip'){
          utils::unzip(paste0(file_path,".zip"), 
                       exdir = normalizePath(paste0(rd,"/",tile_name),
                                             mustWork = FALSE),
                       junkpaths = FALSE)        
          files <- unlist(list.files(paste0(rd,"/",tile_name), pattern = exts,
                                     recursive = TRUE,
                                     full.names = TRUE))
          if (length(files) > 1){
            files <- lapply(files,rast)
            files <- do.call(merge,files) * 1.0
          } else {
            files <- rast(files) * 1.0
          }
          files[files < z_min] <- NA
          if (!is.null(filt) | filt != 0){
            files <- focal(files,w=filt,fun=mean,na.policy='omit')
          }
          names(files) <- 'z'
          
          writeRST(rastToTable(files), normalizePath(file_path))
          unlink(paste0(file_path,".zip"),recursive=TRUE)
        } 
        
        if (extension == 'gz'){
          R.utils::decompressFile(normalizePath(paste0(file_path,".gz"),mustWork=FALSE), 
                                  destname = normalizePath(paste0(rd,"/",tile_name),
                                                           mustWork = FALSE),
                                  remove = TRUE,
                                  FUN = gzfile,
                                  ext = 'gz')        
          files <- unlist(list.files(paste0(rd,"/",tile_name), pattern = exts,
                                     recursive = TRUE,
                                     full.names = TRUE))
          if (length(files) > 1){
            files <- lapply(files,rast)
            files <- do.call(merge,files) * 1.0
          } else {
            files <- rast(files) * 1.0
          }
          files[files < z_min] <- NA
          if (!is.null(filt) | filt != 0){
            files <- focal(files,w=filt,fun=mean,na.policy='omit')
          }
          writeRST(files, normalizePath(file_path))
          unlink(normalizePath(paste0(rd,"/",tile_name),
                               mustWork = FALSE),recursive=TRUE)
        } 
        
        # If downloaded file is an uncompressed file
        exts <- unlist(stringr::str_extract_all(exts,pattern="[a-z]+"))
        if (extension %in% exts){
          files <- rast(paste0(file_path,".",extension)) * 1.0
          files[files < z_min] <- NA
          if (!is.null(filt) | filt != 0){
            files <- focal(files,w=filt,fun=mean,na.policy='omit')
          }
          writeRST(files, normalizePath(file_path))
          unlink(paste0(file_path,".",extension),recursive=TRUE)
        }
      }
    } else if (vals == 'location' & methods::is(unique(polys$location),'numeric')){
      # If what we provide is a polygon from makeGrid but no URLs to the source,
      # download from AWS
      zoom <- unique(polys$location)
      for (i in seq(1,length(down))){
        tile_name <- as.character(polys[down[i],][[tile_id]])
        print(paste0("Downloading Tile ",tile_name," (",
                     i," of ",length(down),")"))
        file_path <- normalizePath(paste0(rd,"/",tile_name),mustWork = FALSE)
        clip <- suppressWarnings(rast(elevatr::get_elev_raster(sf::st_as_sf(as(polys[down[i],],'Spatial')),
                                              z = zoom,
                                              src = 'aws')))
        poly <- which(polys[[tile_id]] == tile_name)

        clip[clip < z_min] <- NA
        if (!is.null(filt) | filt != 0){
        clip <- focal(clip,w=filt,fun=mean,na.policy='omit')
        }
        clip <- crop(clip,polys[poly,], snap = 'out')
        writeRST(clip, file_path)
      }
    } else {
      # If what we provide is a polygon and a raster or a path to a singular
      # DEM as the source, crop the necessary DEM
      if (!((methods::is(vals,"Raster")) | (methods::is(vals,"SpatRaster")))){
        # If dem is not a raster, it's because the source raster is stored
        # elsewhere. Import it
        dem <- rast(unique(polys$location))
      } else {
        if (methods::is(vals,"Raster")){
          dem <- rast(vals)
        } else if (methods::is(vals,'SpatRaster')){
          dem <- vals
        }
      }
      
      dem[dem < z_min] <- NA
      if (!is.null(filt) | filt != 0){
        dem <- focal(dem,w=filt,fun=mean,na.policy='omit')
      }
      # For every tile that needs to be acquired...
      for (i in seq(1,length(down))){
        
        # Get the unique tile id, and define the output filepath
        tile_name <- as.character(polys[down[i],][[tile_id]])
        file_path <- normalizePath(paste0(rd,"/",tile_name),mustWork=FALSE)
        print(paste0("Cropping Tile ",tile_name," (",
                     i," of ",length(down),")"))
        # Select the singular tile, and use it to crop the dem.
        # Save it to the above filepath, zip it, and delete the tiff
        poly <- which(polys[[tile_id]] == tile_name)
        clip <- crop(dem,polys[poly,], snap = 'out') * 1.0
        writeRST(clip,file_path)
      }
    }
  }
}


#' Import \code{GPX} tracks from commercial GPS equipment as a data.table
#' ready for velocity function estimation. Note that the output coordinates
#' must still be converted to a conformal projection in meters if they are
#' to be used in functions other than \code{\link[lbmech]{getVelocity}},
#' \code{\link[lbmech]{dtVelocity}}, or \code{\link[lbmech]{downsample}}.
#' 
#' @title Import GPX tracks
#' @param tracks A character string or vector with filepaths pointing to
#' the location of the GPX files
#' @return A data.table with five or six columns:
#' 
#' (1) \code{$TrackID} The name of the input GPX track
#' 
#' (2) \code{$PID} The order that point appeared in the GPX file
#' 
#' (3) \code{$t} The timestamp
#' 
#' (4) \code{$long} The longitude in decimal degrees
#' 
#' (5) \code{$lat} The latitude in decimal degrees
#' 
#' (6) \code{$z} The elevation (if present)
#' @importFrom data.table data.table
#' @examples 
#' # Get a list of GPX tracks in a directory
#' gpx <-  list.files(pattern = ".gpx$")
#' 
#' # Convert to data.table
#' gpx <- importGPX(gpx)
#' @export
importGPX <- function(tracks){
  tracks_out <- data.table()
  j = 0
  
  pb <- utils::txtProgressBar(min=0,max=length(tracks)+0.00001,style=3)
  for (i in tracks){
    try({
      gpx <- tmaptools::read_GPX(i,
                                 layers='track_points')$track_points
      gpx <- sf::as_Spatial(gpx)
      
      gpx <- data.table(TrackID = stringr::str_extract(i,
                                                       pattern="(?<=/)[0-9A-Za-z_\\-]+(?=.gpx)"),
                        PID = seq(1,nrow(gpx)),
                        t = gpx@data$time,
                        long = gpx@coords[,"coords.x1"],
                        lat = gpx@coords[,"coords.x2"],
                        z = gpx@data$ele)
      
      gpx <- gpx[order(t)]
      
      tracks_out <- rbind(tracks_out,gpx)[]
    },silent=TRUE)
    j <- j + 1
    utils::setTxtProgressBar(pb,j)
    return(tracks_out)
  }
}

#' Downsample high-resolution *(x,y,t)* or *(x,y,z,t)* data
#' to match a minimum spatial resolution. 
#' 
#' @title Downsample high-resolution GPS data
#' @param data A data.frame or something coercible to a data.table containing
#' all observations
#' @param t_step The smallest allowable median time interval between observations.
#' Any track with a median value below this will be resampled to a track of 
#' equally-spaced observations with time difference \code{t_step}. This should
#' be selected based on the parameters that will be used to generate a
#' \code{world} object such that each successive step is likely to fall outside
#' of the range of possible raster transitions. For example, for a given source
#' \code{dem} in \code{\link[lbmech]{makeGrid}}, a given \code{distances = 16} in
#' \code{\link[lbmech]{makeWorld}} (implying that \code{contiguity = 2}), and 
#' an estimated animal maximum velocity of \code{v_max = 1.5} m/s, t_step
#' should be at least \code{t_step = res(dem) / v_max * (contiguity + 1) * sqrt(2)}
#' @param t_cut A numeric. Any gap exceeding this value in seconds in a given
#' track will be treated as the start of a new segment. Default is \code{t_step * 10}.
#' @param x A character string representing the data column containing the 'x'
#' coordinates in any projection.
#' @param y A character string representing the data column containing the 'y'
#' coordinates in any projection.
#' @param z Optional. A character string representing the data column containing 
#' the 'z' elevations. 
#' @param t A character string representing the data column containing the time
#' values, either as a numeric of successive seconds or as a date time string
#' @param ID A character string representing the data column containing the
#' unique ID for each observed trajectory. In other words, each set of points
#' for each continuous observation for each observed individual would merit a
#' unique id.
#' @return A data.table, containing the original input data but with all 
#' overly-high resolution tracks downsampled to an acceptable rate of observations
#' and column names prepared for the \code{\link[lbmech]{getVelocity}} function. 
#' @importFrom data.table data.table
#' @examples 
#' # Generate fake data with x,y coordinates, z elevation, and a t
#' # column representing the number of seconds into the observation
#' data <- data.table(x = runif(10000,10000,20000),
#'                    y = runif(10000,30000,40000),
#'                    z = runif(10000,0,200),
#'                    t = 1:1000,
#'                    ID = rep(1:10,each=1000))
#'                    
#' # Set the minimum value at 3 seconds
#' data <- downsample(data = data, t_step = 3, z = 'z')
#' @export
downsample <- function(data, t_step, t_cut = t_step * 10,
                       x = 'x', y = 'y', z = NULL,
                       t = 't', ID = "ID"){
  # This bit is to silence the CRAN check warnings for literal column names
  ..ID=..x=..y=..t=..z=dt=V1=LegID=NULL
  #
  
  if (is.null(z)){
    data <- data[, .(ID = get(..ID), x = get(..x), y = get(..y), t = get(..t))]
  } else {
    data <- data[, .(ID = get(..ID), x = get(..x), y = get(..y), 
                     z = get(..z), t = get(..t))] 
  }
  
  # Start by calculating dt, make sure it's a numeric
  data[, dt := t - data.table::shift(t), by = 'ID'
  ][, dt := as.numeric(dt,units='secs')]
  
  # Get a list of tracks that exceed the minimum allowable time interval
  upsamp <- data[,stats::median(dt,na.rm=TRUE),by='ID'
  ][V1 < t_step, ID]
  
  # We need a new Unique ID that separates individual tracks
  # into sub-tracks where if a track is high-resolution AND
  # any two sub-points are more than `t_Cut` in time from each
  # other, they are two different tracks
  data[ID %in% upsamp & (dt > t_cut), dt := NA
  ][is.na(dt), LegID := 1
  ][!is.na(dt), LegID := 0
  ][, LegID := cumsum(LegID)] 
  
  legs <- unique(data[,.(ID,LegID)])
  
  # Any leg that belongs to a high-res track must be downsampled
  whichlegs <- legs[ID %in% upsamp, LegID]
  
  # Add modified data to this blank data.table
  newData <- data.table()
  pb <- utils::txtProgressBar(min = 0 , max = max(whichlegs),style=3)
  i = 0
  
  # Iterate over every 'leg' subtrack, performing a smooth
  # spline to join the points and sampling at t_step
  for (track in whichlegs){
    
    # Select points corresponding to the leg
    points <- data[LegID == track]
    tryCatch({
      # Save initial time; temporarily strip units if present
      t_i <- points[1,t]
      points[1,dt := 0][, t := cumsum(dt)]
      
      # Get smoothing splines for xy
      x_sp <- stats::smooth.spline(points[,.(x = t, y = x)])
      y_sp <- stats::smooth.spline(points[,.(x = t, y = y)])
      
      # Get a list of times where we'll sample from
      t_new <- seq(from = 0, to = max(points$t), by = t_step)
      if (is.null(z)){
        # If no z is provided
        points <- data.table(
          ID = legs[LegID == track, ID],
          t = t_i + t_new,
          x = stats::predict(x_sp,t_new)$y,
          y = stats::predict(y_sp,t_new)$y,
          dt = t_step,
          LegID = track
        )
      } else {
        # Smooth spline if z was provided
        z_sp <- stats::smooth.spline(points[,.(x = t, y = z)])
        points <- data.table(
          ID = legs[LegID == track, ID],
          t = t_i + t_new,
          x = stats::predict(x_sp,t_new)$y,
          y = stats::predict(y_sp,t_new)$y,
          z = stats::predict(z_sp,t_new)$y,
          dt = t_step,
          LegID = track
        )
      }
      newData <- rbind(newData,points)
    }, error = function(x) {
      1
    })
    i <- i + 1
    utils::setTxtProgressBar(pb,val = i)
  }
  changed <- unique(newData$ID)
  return(
    if (is.null(z)){
      rbind(data[!(ID %in% changed)], newData)[order(t)
      ][order(ID),.(ID,LegID,x,y,t,dt)]
    } else {
    rbind(data[!(ID %in% changed)], newData)[order(t)
      ][order(ID),.(ID,LegID,x,y,z,t,dt)]
    }
  )
}

#' Calculate the velocity function for an animal from \emph{(x,y,z,t)} data
#' such as from GPS collars, assuming a function of form Tobler (see ####)
#'
#' @title Calculate a velocity function from data
#' @param data A data.frame or something coercible to a data.table containing
#' all observations
#' @param x A character string representing the data column containing the 'x'
#' coordinates in meters or degrees. Ignored for distance calculations if \code{dl} is,
#' provided but required to extract elevations if \code{z} is of class
#' Raster* or SpatialPolygonsDataFrame.
#' @param y A character string representing the data column containing the 'y'
#' coordinates in meters or degrees. Ignored for distance calculations if \code{dl} is,
#' provided but required to extract elevations if \code{z} is of class SpatRaster,
#' Raster*, or SpatVector, SpatialPolygonsDataFrame.
#' @param degs Logical. If FALSE, the \code{getVelocity} proceeds as if the
#' input coordinates are meters in a projected coordinate system. If 
#' TRUE, assumes the input coordinates are decimal degrees in lat/long, with
#' \code{x} giving the longitude column name and \code{y} the latitude. See
#' \code{\link[geosphere]{distGeo}}.
#' @param dl A character string representing the data column containing the
#' net displacement between this observation and the previous in meters. Optional.
#' @param z Either a character string, a SpatRaster or Raster*, or a SpatVector
#' object. If character string, it represents the data column containing the 'z'
#' coordinates/elevations in meters. If a SpatRaster or Raster*, a DEM containing
#' the elevation in meters. If SpatVector, it must represent the
#' sectors and filepaths/URLs corresponding to the region's elevations
#' (see the output of \code{\link[lbmech]{makeGrid}}).
#' @param dt A character string representing the data column containing the
#' time difference from the previous observation in seconds.
#' @param ID A character string representing the data column containing the
#' unique ID for each observed trajectory. In other words, each set of points
#' for each continuous observation for each observed individual would merit a
#' unique id.
#' @param tau A number between 0 and 1, representing a global cutoff
#' value for \code{tau_vmax} and \code{tau_nlrq}. Ignored if the latter two are provided.
#' @param tau_vmax A number between 0 and 1, representing the percentile at which
#' the maximum velocity is calculated. In other words, if \code{tau_vmax = 0.995} (the default),
#' the maximum velocity will be at the 99.5th percentile of all observations.
#' @param tau_nlrq A number between 0 and 1, representing the percentile at which
#' the nonlinear regression is calculated. In other words, if \code{tau_nlrq = 0.95} (the default),
#' the total curve will attempt to have at each interval 5\% of the observations above the
#' regression and 95\% below.
#' @param k_init A number representing the value for the topographic sensitivity
#' at which to initiate the nonlinear regression. Default is \code{k_init = 3.5}.
#' @param alpha_init A number representing the value for dimensionless slope of
#' maximum velocity at which to initiate the nonlinear regression. Default is
#' \code{alpha_init = -0.05}.
#' @param v_lim The maximum velocity that will be considered. Any value above
#' this will be excluded from the regression. Default is \code{v_lim = Inf},
#' but it should be set to an animal's maximum possible velocity.
#' @param slope_lim the maximum angle that will be considered. Any value
#' above this will be excluded from the regression. Default is \code{slope_lim = 1}.
#' @param tile_id a character string representing the name of the column
#' in the \code{z} polygon containing the unique Tile IDs. Ignored if elevations are
#' provided as a column or Raster*. Otherwise default is \code{tile_id = 'TILEID'}.
#' @param vals A character string or a Raster* object. Required only if the
#' \code{z} parameter is a polygon NOT the output of the 
#' \code{\link[lbmech]{makeGrid}} function as the default is
#' the character string \code{'location'}. If not, the \code{vals} parameter should be
#' set to the column name containing the URL or file path to the DEM for that
#' sector.
#' @param dir A filepath to the directory being used as the workspace.
#' Default is \code{tempdir()} but unless the analyses will only be performed a few
#' times it is highly recommended to define a permanent workspace.
#' @return A list, containing the following entries:
#'
#' (1) \code{$model}, containing an object of class \code{\link[quantreg]{nlrq}} containing the output
#' model from the nonlinear quantitative regression.
#'
#' (2) \code{$vmax}, containing the identified maximum velocity.
#'
#' (3) \code{$alpha}, containing the identified angle of maximum velocity.
#'
#' (4) \code{$k}, containing the identified topographic sensitivity factor.
#'
#' (5) \code{$tau_max}, containing the employed tau_max.
#'
#' (6) \code{$tau_nlrq}, containing the employed tau_nlrq.
#'
#' (7) \code{$data}, containing a data.table with the original data in a standardized format
#' @importFrom data.table as.data.table
#' @importFrom data.table :=
#' @importFrom data.table data.table
#' @importFrom terra cellFromXY
#' @importFrom terra crs
#' @importFrom terra rast
#' @importFrom terra vect
#' @examples 
#' # Note that the output results should be senseless since they
#' # are computed on random data
#' 
#' #### Example 1:
#' # If the data contains an 'elevation' or 'z' column
#' data <- data.table(x = runif(10000,10000,20000),
#'                    y = runif(10000,30000,40000),
#'                    elevation = runif(10000,0,200),
#'                    dt = 120,
#'                    ID = rep(1:10,each=1000))
#' velocity <- getVelocity(data = data, z = 'elevation')
#' 
#' 
#' #### Example 2:
#' # If the data do not contain elevation, and 'z' is a raster
#' data$z <- NULL
#' 
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
#' velocity <- getVelocity(data = data, z = dem)
#' 
#' 
#' ### Example 3:
#' # If the data do not contain elevation, and 'z' is a sector grid pointing
#' # to file locations
#' 
#' # Export the DEM so it's not just stored on the memory
#' dir <- tempdir()
#' writeRaster(dem, paste0(dir,"/DEM.tif"),overwrite=TRUE)
#' 
#' # Import raster, get the grid
#' dem <- rast(paste0(dir,"/DEM.tif"))
#' grid <- makeGrid(dem = dem, nx = n, ny = n)
#' 
#' velocity <- getVelocity(data = data, z = grid, dir = dir)
#' 
#' @export
getVelocity <- function(data, x = 'x', y ='y', degs = FALSE, dl = NULL, z = 'z', 
                        dt = 'dt', ID = 'ID', tau = NULL, tau_vmax = 0.995,
                        tau_nlrq = 0.95, k_init = 3.5, alpha_init = -0.05,
                        v_lim = Inf, slope_lim = 1,
                        tile_id = "TILEID", vals = "location",
                        dir = tempdir()){
  # This bit is to silence the CRAN check warnings for literal column names
  ..x=..y=..z=..ID=..dt=..z_vals=Order=dx=dy=dz=dl_dt=..dl=dz_dl=NULL
  #
  
  # Standardize filepath names
  dir <- normalizePath(dir)
  rd <- normalizePath(paste0(dir,"/Elevations"),mustWork = FALSE)
  
  data <- as.data.table(data)
  
  if (!(ID %in% names(data))){
    stop(paste0("'",ID,"' not found among data columns"))
  }
  
  # If z is a raster its simply extracting the z points from xy locations
  if (methods::is(z,"Raster")){
    z <- rast(z)
  } 
  if (methods::is(z, "SpatRaster")){
    dem <- z
    z <- dem[cellFromXY(dem,data[,.(get(..x),get(..y))])]
    data[,z := ..z]
    if (is.null(dl)){
    data <- data[,.(ID = get(..ID),x = get(..x),y = get(..y),z,dt = get(..dt))]
    } else {
      data <- data[,.(ID = get(..ID),x = get(..x),y = get(..y),z,dt = get(..dt))]
    }
    
    # If z is a character, do nothing; just rename columns
  } else if ("character" %in% class(z)){
    
    if (is.null(dl)){
    data <- data[,.(ID = get(..ID), x = get(..x), y = get(..y),
                    z = get(..z),dt = get(..dt))]
    } else if ((x %in% names(data)) & (y %in% names(data))){
      data <- data[,.(ID = get(..ID), x = get(..x), y = get(..y), dl = get(..dl),
                      z = get(..z), dt = get(..dt))]  
    } else{
      data <- data[,.(ID = get(..ID),dl = get(..dl), z = get(..z), dt = get(..dt))]  
    }
    
    
    # If z is a polygon defining sector locations...
  } else if (methods::is(z,"SpatVector")){
    
    # Create an FID to keep the observations in order since the process will
    # shuffle them
    data$Order <- seq(1,nrow(data))
    
    # If z locations have already been provided but the input contains a DEM,
    # warn the user
    if (!is.null(data$z)){
      warning("Data contains a 'z' column which will be ignored")
      data[, z:= NULL]
    }
    
    # Create a SPDF object to detect which tiles will be employed
    data_points <- vect(data, crs = crs(z),
                        geom=c(x,y),
                        keepgeom=TRUE)
    data_points <- as.data.table(intersect(data_points,z))
    tiles <- unique(data_points[,get(tile_id)])
    
    # Make sure that the necessary files exist
    getMap(tiles=tiles,polys=z,tile_id=tile_id,vals=vals,
           dir= normalizePath(paste0(dir,"/Elevations/")))
    
    # Create an empty data.table in which to store extracted values. There's a
    # small chance some points will land in two sectors, so we'll have
    # to collect all and consolidate
    data_new <- data.table()
    for (tile in tiles){
      elevs <- importRST(normalizePath(paste0(rd,"/",tile),mustWork=FALSE))
      
      d <- data_points[get(tile_id) == tile]
      z_vals <- elevs[cellFromXY(elevs,d[,.(get(..x),get(..y))])]
      d[,z := ..z_vals]
      
      # Add it to the empty data.table, and remove the uncompressed tiff.
      data_new <- rbind(data_new,d)
    }
    
    # Group by all columns that ARE NOT z, setting z equal to the mean of all
    # non-NA values. This makes sure we have one row per observation and lets
    # us handle edge cases where multiple z values were assigned to points on
    # the border between two tiles.
    cols <- names(data_new)
    cols <- cols[cols != "z"]
    data_new <- data_new[,.(z = mean(z,na.rm=TRUE)),by=cols][]
    
    if (is.null(dl)){
    data <- data_new[order(Order)
    ][,.(ID = get(..ID),x = get(..x),y = get(..y),z,dt = get(..dt))
    ][]
    } else {
      data <- data_new[order(Order)
      ][,.(ID = get(..ID),x = get(..x),y = get(..y),dl = get(..dl),
           z,dt = get(..dt))
      ][]
    }
  }
  
  # Calculate displacement, then speed and slope
  if (is.null(dl)){
    if (degs == FALSE) {
  data[, `:=`(dx = x - data.table::shift(x),
              dy = y - data.table::shift(y),
              dz = z - data.table::shift(z)),
       by = 'ID'
  ][, dl := sqrt(dx^2 + dy^2)]
    } else if (degs == TRUE){
      # If what we provided was latitude and longitude
      data[, `:=`(dl = data.table::shift(
        geosphere::distGeo(
          matrix(c(x,y),
                 nrow=length(x)))),
        dz = z - data.table::shift(z)),
        by = 'ID']
    }
  } else {
    data[, `:=`(dl = get(..dl),
                dz = z - data.table::shift(z)), 
         by = 'ID']
  }
  
  data[, `:=`(dl_dt = dl / dt,
           dz_dl = dz / dl)]
  
  # The user can set a global tau, or a tau for the maximum velocity
  # and the non-linear quantile regression independently. This section
  # lets us deal with the event that the user sets all three or a global
  # and only one of the independents.
  if (!is.null(tau) & (!is.null(tau_vmax) & !is.null(tau_nlrq))){
    warning("Tau is specified for both vmax and nlrq. Ignoring global tau")
  }
  if (!is.null(tau) & (is.null(tau_vmax) | is.null(tau_nlrq))){
    warning("Tau not specified for vmax and/or nlrq. Ignoring tau_vmax/tau_nlrq
            and employing global tau")
    tau_nlrq <- tau
    tau_vmax <- tau
  }
  
  # Get the maximum velocity as the tauth_vmax quantile
  v_max <- as.numeric(stats::quantile(data[(dl_dt <= v_lim) &
                                             abs(dz_dl) <= slope_lim,dl_dt],
                                      tau_vmax,na.rm=TRUE))
  data$v_max <- v_max
  
  
  # And obtain the other coefficients through an nlrq of the form proposed
  # by Tobler (exponential decay from an optimal angle)
  velocity <- quantreg::nlrq(dl_dt ~ v_max * exp(-k * abs(dz_dl - alpha)),
                             data = data[(dl_dt <= v_lim) & 
                                           abs(dz_dl) <= slope_lim], 
                             tau = tau_nlrq, 
                             start=list(k=k_init,alpha=alpha_init))
  data$v_max <- NULL
  
  
  # We'll store the output in a list. Best practices would usually suggest a
  # unique class, but it's not necessary and would require the user to load
  # nlrq even if they don't employ this function.
  
  # Setting the velocity always throws an unnecessary warning, so quiet it temporarily
  defaultW <- getOption("warn")
  options(warn = -1)
  
  out <- list(
    model = velocity,
    vmax = as.numeric(v_max),
    alpha = as.numeric(stats::coefficients(velocity)[['alpha']]),
    k = as.numeric(stats::coefficients(velocity)[['k']]),
    tau_vmax = tau_vmax,
    tau_nlrq = tau_nlrq,
    data = data
  )
  # Return the environment to what it was
  on.exit(options(warn = defaultW))
  return(out)
  
}

#' A wrapper for \code{\link[lbmech]{getVelocity}} for use inside
#' of a data.table with observations from multiple sources
#' (be it different individual animals and/or different instances
#' from the same animal). In a \code{\link[data.table]{data.table}},
#' use this function in the \code{j} slot, passing it along \code{.SD}.
#' 
#' @title Estimate multiple velocity functions in data.tables
#' @param data A data.frame or something coercible to a data.table 
#' containing all observations.
#' @param ... Additional arguments to pass along to
#' \code{\link[lbmech]{getVelocity}}. Set the \code{ID} parameter
#' parameter to the column name representing each individual track
#' segment. Otherwise, employ the same singular column name as in \code{by=}
#' (data.table syntax).
#' @return A list with five entries: \code{$k}, the sample's
#' topographic sensitivity factor and its associated
#' p-value \code{$k_p}; \code{$alpha}, the sample's slope of fastest
#' travel and its associated p-value \code{$alpha_p}; and the
#' estimated maximum walking speed \code{$v_max}. 
#' @examples 
#' # Generate fake data with x,y coordinates, z elevation, and a t
#' # column representing the number of seconds into the observation
#' data <- data.table(x = runif(10000,10000,20000),
#'                    y = runif(10000,30000,40000),
#'                    z = runif(10000,0,200),
#'                    dt = 15,
#'                    ID = rep(1:10,each=1000))
#'                    
#' # To get the velocity function for all observations together
#' v1 <- getVelocity(data)
#' 
#' # This is the same as above, but it only returns a list with the 
#' # coefficients and p-values
#' v2 <- dtVelocity(data)
#' 
#' # Instead this function is best to get the coefficients for 
#' # each individual track un a data.table using .SD
#' v3 <- data[, dtVelocity(.SD), by = 'ID', .SDcols = names(data)]
#' @export
dtVelocity <- function(data, ...){
  j <- list()
  tryCatch(
    {
      n <- getVelocity(data = data, ...)
      j <- as.list(stats::coef(n$model))
      j[['v_max']] <- as.numeric(n$vmax)
      n <- summary(n$model)$coefficients
      j[['k_p']] <- as.numeric(n[,"Pr(>|t|)"][1])
      j[['alpha_p']] <- as.numeric(n[,"Pr(>|t|)"][2])
      return(j)
    }, 
    error = function(e) {
      j[['k']] <- as.numeric(NA)
      j[['alpha']] <- as.numeric(NA)
      j[['v_max']] <- as.numeric(NA)
      j[['k_p']] <- as.numeric(NA)
      j[['alpha_p']] <- as.numeric(NA)
      return(j)
    })
}



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



#' Function that defines the grid that can be traversed--the "world"--as well as the
#' cells that can be accessed from each individual cell. This is the most
#' time-intensive function. 
#' 
#' It first checks to see if the required transition \code{.gz}
#' files have already been created in the \code{dir} workspace. If not, it checks to see if the
#' DEMs required to generate the transition \code{.gz} files have been downloaded/cropped,
#' and generates the latter if so. If not, it downloads/crops them and proceeds.
#' 
#' The default parameters are sufficient for a workflow involving calculating
#' costs with the \code{\link[lbmech]{calculateCosts}} function. However, if
#' non-energetic analyses are desired, or if the user wishes to employ functions
#' relying on the raw values at each cell (not just the *difference*), \code{keep_z}
#' should be modified accordingly. 
#'
#' @title Define which cells are adjacent
#' @param tiles A character vector--such as the output to
#' \code{\link[lbmech]{whichTiles}}---containing the unique tile IDs for sectors that
#' should be in the workspace.
#' @param polys An object of class SpatVector representing
#' the partitioning grid for the maximum possible area, in the same format as the
#' output of the \code{\link[lbmech]{makeGrid}} function.
#' @param tile_id A character string representing the name of the column
#' in the 'polys' polygon containing the unique Tile IDs. Default is \code{tile_id = 'TILEID'}
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
#' @param keep_z A \code{NULL}, character string, or character vector consisting of 
#' one, two, or three of (respectively) \code{c('z_i','z_f','dz')} indicating which
#' of the initial or final cell's values should be stored in the output tensor.
#' \code{keep_z = 'all'} is shorthand for \code{c('z_i','z_f','dz')}.
#' Since the energetic calculations only depend on a *change* in elevation, 
#' the default is \code{'dz'}. However, if the specific origin or destination
#' values are needed this should me modified accordingly. 
#' @param z_fix A SpatRaster or Raster* that will define the resolution, origin, and
#' projection information for the entire "world" of possible movement. Note that
#' it does NOT need the same extent. Default resolution is 5, and offset is 0.
#' Default projection is \code{proj = crs(polys)} unless a `z_fix` or `proj` is 
#' provided, in which case the latter is ignored. Great care should be 
#' employed to ensure that the projection is conformal and in meters. 
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
#' @param FUN Function to deal with overlapping values for overlapping sectors.
#' Default is NA, which uses \code{\link[terra]{merge}}. 
#' To use \code{\link[terra]{mosaic}}, provide a compatible function.
#' @param sampling How to resample rasters. Default is \code{'bilinear'} interpolation,
#' although \code{'ngb'} nearest neighbor is available for categorical rasters. 
#' @param dir A filepath to the directory being used as the workspace.
#' Default is \code{tempdir()} but unless the analyses will only be performed a few
#' times it is highly recommended to define a permanent workspace.
#' @param ... Additional arguments to pass to \code{\link[lbmech]{fix_z}}.
#' @return A \code{.gz} file for each sector named after its sector id,
#' containing a data.table with three columns in the default setting:
#'
#' (1) \code{$from}, a character string of all possible origin cells in format "x,y",
#' rounded to the next-lowest integer
#'
#' (2) \code{$to},  a character string of all possible destination cells in format "x,y"
#' rounded to the next-lowest integer
#'
#' (3) \code{$dz}, an integer representing the change in elevation for each origin-destination pair
#' 
#' Additionally, \code{$z_i} and \code{$z_f} columns containing the values of the
#' initial and final cells may appear depending on the setting of the 
#' \code{keep_z} parameter.
#' @importFrom terra res
#' @importFrom terra project
#' @importFrom terra mosaic
#' @importFrom terra merge
#' @importFrom terra crs
#' @importFrom terra crs<-
#' @importFrom terra xFromCell
#' @importFrom terra yFromCell
#' @importFrom terra ext
#' @importFrom terra crop
#' @importFrom terra cells
#' @importFrom terra adjacent
#' @importFrom terra ncell
#' @importFrom terra intersect
#' @importFrom terra buffer
#' @importFrom data.table as.data.table
#' @importFrom data.table :=
#' @importFrom data.table .SD
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
#' grid <- makeGrid(dem = dem, nx = n, ny = n)
#' 
#' # Select all tiles that exist between x = (12000,16000) and y = (32000,36000)
#' tiles <- ext(c(12000,16000,32000,36000))
#' tiles <- as.polygons(tiles)
#' crs(tiles) <- crs(grid)
#' tiles <- whichTiles(region = tiles, polys = grid)
#' 
#' # Make a world but limit it to the DEM grid size
#' makeWorld(tiles = tiles, polys = grid, cut_slope = 0.5, 
#'           res = res(dem), dir = dir)
#' @export 
makeWorld <- function(tiles,polys,cut_slope,tile_id = 'TILEID', proj = crs(polys),
                      directions = 16, neighbor_distance = 100, keep_z = 'dz',
                      z_fix = NULL, unit = "m", vals = 'location', precision = 2, 
                      FUN = NULL, sampling = 'bilinear', dir = tempdir(), ...){
  # This bit is to silence the CRAN check warnings for literal column names
  from=to=..dem=z_f=z_i=x_f=x_i=y_f=y_i=dz=..cols=NULL
  rm(..cols)
  
  dir <- normalizePath(dir)
  units <- c("m","km","ft","mi")
  if (!(unit %in% units)){
    stop("Error: Unknown units. Only m, km, ft, and mi allowed.
         If your data source is in another unit, please convert units
         before running the getMap function.")
  }
  if (directions == 4 | directions == 8){
    contiguity <- 1
  } else if (directions == 16){
    contiguity <- 2
  } else{
    contiguity <- base::ceiling(max(nrow(directions),ncol(directions)))
  }
  
  if (all(keep_z == 'all')){
    keep_z <- c("dz","z_i","z_f")
  }
  
  # If no z_fix is provided, make one
  if (is.null(z_fix)){
    z_fix <- fix_z(proj = proj, ...)
  } else if (methods::is(z_fix, "Raster")){
    z_fix <- rast(z_fix)
  }
  l_p <- res(z_fix)
  
  rd <- normalizePath(paste0(dir,"/Tensors"),mustWork=FALSE)
  if (!dir.exists(rd)){
    dir.create(rd)
  }
  for (i in tiles){
    if (!file.exists(normalizePath(paste0(rd,"/",i,".fst"),mustWork=FALSE))){
      # Get the appropriate polygon and name; create a temporary folder
      tensors <- normalizePath(paste0(rd,"/",i),mustWork=FALSE)
      poly <- which(polys[[tile_id]] == i)
      poly <- polys[poly,]
      
      # Get the neighboring sectors and make sure they're all downloaded
      poly_buff <- buffer(poly,width=neighbor_distance)
      mzip <- intersect(poly_buff[,NA],polys)[[tile_id]]
      mzip <- unlist(mzip)
      getMap(mzip, polys = polys, tile_id = tile_id, vals = vals,dir = 
               normalizePath(paste0(dir,"/Elevations/")))
      
      # Import all DEMs
      dem <- normalizePath(paste0(dir,"/Elevations/",mzip),mustWork=FALSE)
      dem <- lapply(dem,importRST)
      dem <- suppressWarnings(lapply(dem,
                                     project,
                                     y = dem[[1]],
                                     align = TRUE))
      
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
      
      # Fix the units if they weren't provided in meters
      if (unit == "km"){
        dem <- dem / 1000
      } else if (unit == "ft"){
        dem <- dem*12*2.54/100
      } else if (unit == "mi"){
        dem <- dem*5280*12*2.54/100
      }
      
      # Create a cropping polygon. This is equivalent to the extent of the
      # tile shape, expanded in the X and Y directions by the
      # pixel size times the connectivity contiguity. This has to be
      # reprojected back to the original coordinate system since only it
      # preserves orthogonality
      poly <- project(poly,dem)
      poly <- unlist(ext(poly))
      poly <- c(poly[1] - (contiguity + 1) * l_p[1], poly[2] + (contiguity + 1) * l_p[1],
                poly[3] - (contiguity + 1) * l_p[2], poly[4] + (contiguity + 1) * l_p[2])
      poly <- ext(poly)
      poly <- as.polygons(poly)
      crs(poly) <- crs(dem)
      
      # Crop the mosaiced raster by the cropping polygon, project it
      dem <- crop(dem,poly)
      dem <- suppressWarnings(project(dem,z_fix,align = TRUE))
      
      name <- i
      
      nas <- cells(dem)
      
      # Get pairs of adjacent all adjacent cells; drop those that
      # correspond to NA
      adj <- adjacent(dem,nas,directions=directions,pairs=TRUE)
      adj <- as.data.table(adj)
      
      # Calculate the change in elevation between every accessible cell pairs,
      # then drop all values that would require movement over the
      # cut slope.
      adj[, `:=`(z_i = unlist(..dem[from]), z_f = unlist(..dem[to]),
                 x_i = xFromCell(dem,from), y_i = yFromCell(dem,from),
                 x_f = xFromCell(dem,to), y_f = yFromCell(dem,to))
      ][, (c("x_i","y_i","x_f","y_f")) :=
          lapply(.SD,round,precision),
        .SDcols = c("x_i","y_i","x_f","y_f")
      ][, `:=`(from = paste(format(x_i, scientific = FALSE),
                            format(y_i, scientific = FALSE), sep=","),
               to = paste(format(x_f, scientific = FALSE),
                          format(y_f, scientific = FALSE), sep=","))
      ][, (c("from","to")) := lapply(.SD,stringr::str_remove_all," "),
        .SDcols = c("from","to")]
      
      if ('dz' %in% keep_z){
        adj[, `:=`(dz = z_f - z_i)]
      }
      
      adj <- stats::na.omit(adj)
      
      # This exports the cost tensor so that it can be called later by
      # the lbm.distance tool to calculate paths and catchments
      cols <- c("from","to",keep_z)
      fst::write_fst(adj[,..cols],
              path = normalizePath(paste0(rd,"/",i,".fst"),mustWork=FALSE))
    }
  }
}


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


#' A function that for a given region imports all cells from the
#' transition \code{.gz} files. If such files have not yet been generated,
#' they can be created by passing along the necessary parameters to this
#' function as with \code{\link[lbmech]{makeWorld}}.
#' 
#' The default parameters are sufficient for a workflow involving calculating
#' costs with the \code{\link[lbmech]{calculateCosts}} function. However, if
#' non-energetic analyses are desired, or if the user wishes to employ functions
#' relying on the coordinates of each cell or the distance traveled in each
#' segment, \code{xy1}, \code{xy2} and \code{dl} should be modified accordingly.
#' 
#' After running this function, the user may plug in the output world into
#' the \code{\link[lbmech]{calculateCosts}} function. Alternatively, the user
#' may perform algebraic operations using any of the 
#' \code{c('x_i','x_f','y_i','y_f','z_i','z_f','dl','dz')} columns.
#' @title Import a world where movement is possible
#' @param region An object of class Raster* or SpatialPolygons*
#' representing the total area where movement is possible. Must lie within
#' the area defined by \code{polys}
#' @param polys An object of class SpatialPolygonsDataFrame representing
#' the partitioning grid for the maximum possible area, in the same format as the
#' output of the \code{\link[lbmech]{makeGrid}} function.
#' @param proj A crs object or character string representing the output
#' projection. Default is \code{proj = crs(polys)}
#' unless \code{proj} or \code{z_fix} are provided in which case 
#' the latter takes precedence. 
#' @param banned An object of class Raster* or SpatialPolygons*
#' representing the total area where movement is \emph{prohibited}. Must lie within
#' the area defined by \code{polys}
#' @param tile_id A character string representing the name of the column
#' in the \code{polys} polygon containing the unique Tile IDs. Default is \code{tile_id = 'TILEID'}
#' @param dir A filepath to the directory being used as the workspace.
#' Default is \code{tempdir()} but unless the analyses will only be performed a few
#' times it is highly recommended to define a permanent workspace.
#' @param xy1 A logical \code{TRUE} or \code{FALSE} indicating whether to generate
#' unique numeric columns for the x and y coordinates of the starting points.
#' Default is \code{FALSE}, and is not needed for the 
#' \code{\link[lbmech]{calculateCosts}}.
#' If you intend on using a custom cost function, set this as TRUE since 
#' \code{\link[lbmech]{getCosts}} will need them. 
#' @param xy2 A logical \code{TRUE} or \code{FALSE} indicating whether to generate
#' unique numeric columns for the x and y coordinates of the ending points.
#' Default is \code{FALSE}, and is not needed for the 
#' \code{\link[lbmech]{calculateCosts}} function. 
#' @param dl A logical \code{TRUE} or \code{FALSE} indicating whether to 
#' Default is \code{FALSE}, and is not needed for the 
#' \code{\link[lbmech]{calculateCosts}} function, which generates it on its own.
#' @param dist A character string specifying the method by which to calculate
#' linear distances. One of \code{proj} (for source projection distance units; the fastest),
#' \code{karney} (the second most accurate and second fastest; the default), 
#' \code{cosine}, \code{haversine}, \code{meeus}, 
#' \code{vincentyEllipsoid} or \code{vincentySphere} (see documentation for the
#' \code{\link[geosphere]{geosphere}} package for further documentation). Ignored
#' unless \code{dl = TRUE}.
#' @param r The earth's radius, or semimajor axis (a) for ellipsoidal models. 
#' Default is 6378137 for WGS1984. Ignored unless \code{dl == TRUE} and 
#' \code{dist != 'proj'}.
#' @param f The ellipsoidal flattening of the earth for ellipsoidal models. 
#' Default is \code{f=1/298.257223563} for WGS1984.  Ignored unless \code{dl == TRUE} and 
#' \code{dist == c('karney','meeus','vincentyEllipsoid')}.
#' @param b The earth's semiminor axis. 
#' Default is \code{b=6356752.3142} for WGS1984.  Ignored unless \code{dl == TRUE} and 
#' \code{dist == 'vincentyEllipsoid'}.
#' @param z_fix A SpatRaster with the same origin and resolution as the
#' \code{z_fix} used to generate the 'world' with \code{\link[lbmech]{makeWorld}}.
#' Do not modify this parameter if you didn't modify it when you ran
#' \code{\link[lbmech]{makeWorld}}.
#' @param res Passed on to \code{\link[lbmech]{fix_z}}, default is 5. Ignored
#' if z_fix is provided.
#' @param dx Passed on to \code{\link[lbmech]{fix_z}}, default is 0. Ignored
#' if z_fix is provided.
#' @param dy Passed on to \code{\link[lbmech]{fix_z}}, default is 0. Ignored
#' if z_fix is provided.
#' @param ... Additional arguments to pass to \code{\link[lbmech]{makeWorld}}
#' @return An object of class data.table containing three  under the default settings:
#'
#' (1) \code{$from}, a character string of all possible origin cells in format "x,y",
#' rounded to the next-lowest integer
#'
#' (2) \code{$to},  a character string of all possible destination cells in format "x,y"
#' rounded to the next-lowest integer
#'
#' (3) \code{$dz}, an integer representing the change in elevation for each origin-destination pair
#' 
#' If \code{xy1} and/or \code{xy2} are \code{TRUE}, then an additional columns
#' for the start and/or end 'x' and 'y' coordinates are included (e.g. \code{$x_i}
#' for the initial 'x' and \code{$y_f} for the final 'y'). If \code{dl = TRUE},
#' then a final \code{$dl} column is added with the total displacement based on the
#' chosen \code{dist} method. 
#' @importFrom data.table data.table
#' @importFrom data.table fread
#' @importFrom terra ext
#' @importFrom terra res
#' @importFrom terra rast
#' @importFrom terra crs
#' @importFrom terra crs<-
#' @importFrom terra project
#' @examples 
#' 
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
#' grid <- makeGrid(dem = dem, nx = n, ny = n)
#' 
#' # Import the data lying between x = (12000,16000) and y = (32000,36000)
#' region <- ext(c(12000,16000,32000,36000))
#' region <- as.polygons(region)
#' crs(region) <- crs(grid)
#' 
#' world <- importWorld(region = region, polys = grid, cut_slope = 0.5, 
#'                      res = res(dem), dir = dir)
#' @export
importWorld <- function(region, polys, proj = crs(polys), banned = NULL,
                        tile_id = "TILEID", xy1 = FALSE, xy2 = FALSE, 
                        res = 5, dx = 0, dy = 0, dl = FALSE, dist = 'karney',
                        r = 6378137, f = 1/298.257223563, b = 6356752.3142,
                        z_fix = NULL, dir = tempdir(), ...){
  # This bit is to silence the CRAN check warnings for literal column names
  from=to=dz=x_f=x_i=y_f=y_i=z_f=z_i=x=y=long_f=lat_f=long_i=lat_i=NULL
  #
  if (is.null(z_fix)){
    z_fix <- fix_z(proj = proj, res = res, dx = dx, dy = dy)
  }
  
  z_temp <- data.table(
    expand.grid(x = seq(from = ext(region)[1] - 2 * res(z_fix)[1],
                        to = ext(region)[2]  + 2 * res(z_fix)[1],
                        by = res(z_fix)[1]),
                y = seq(from = ext(region)[3] - 2 * res(z_fix)[2],
                        to = ext(region)[4] + 2 * res(z_fix)[2],
                        by = res(z_fix)[2])))
  z_temp <- rast(z_temp[,.(x,y,z=1)])
  crs(z_temp) <- crs(z_fix)
  z_fix <- suppressWarnings(project(z_temp, z_fix, align = TRUE))
  
  # First determine which tiles to import,
  # download maps and make tensors if needed
  dir <- normalizePath(dir)
  
  tiles <- whichTiles(region = region, polys = polys,
                      tile_id = tile_id)
  makeWorld(tiles = tiles, polys = polys, tile_id = tile_id, dir = dir,
            z_fix = z_fix, ...)
  
  tiles <- normalizePath(paste0(dir,"/Tensors/",tiles,".fst"),
                         mustWork=FALSE)
  
  # Convert the regions into a list of the allowable cells
  region <- regionMask(region, z_fix =  z_fix, res = res, dx = dx, dy = dy)
  if (!is.null(banned)){
    banned <- regionMask(banned, z_fix = z_fix, res = res, dx = dx, dy = dy)
  } else{
    banned <- c("NULL")
  }
  
  # Create an empty data.table to which to add the imported tensors
  Edges <- data.table()
  
  pb <- utils::txtProgressBar(max = length(tiles), style = 3)
  
  # Import tensors one-by-one, keeping only the allowable cells
  for (i in seq(1,length(tiles))) {
    import <- fst::read_fst(tiles[i], as.data.table = TRUE)
    import <- import[(!(from %in% banned) | !(to %in% banned)) &
                       ((from %in% region) | (to %in% region)),]
    Edges <- rbind(Edges, import)
    utils::setTxtProgressBar(pb,i)
  }
  
  # Convert to numeric any 'z' columns to numeric'
  if ("dz" %in% names(Edges)){
    Edges[, `:=`(dz = as.numeric(dz))]
  }
  if ("z_i" %in% names(Edges)){
    Edges[, `:=`(z_i = as.numeric(z_i))]
  }
  if ("z_f" %in% names(Edges)){
    Edges[, `:=`(z_f = as.numeric(z_f))]
  }
  
  # If the user requests the tensor prepared such that origin
  # XY coordinates are reported as numeric columns
  if ((xy1 == TRUE) | (dl == TRUE)){
    Edges[, c("x_i","y_i") := transpose(stringr::str_split(from,","))
    ][, c("x_i","y_i") := lapply(.SD,as.numeric),
      .SDcols = c("x_i","y_i")]
  }
  
  # Destination XY coordinates as numeric columns
  if ((xy2 == TRUE | dl == TRUE)){
    Edges[, c("x_f","y_f") := transpose(stringr::str_split(to,","))
    ][, c("x_f","y_f") := lapply(.SD,as.numeric),
      .SDcols = c("x_f","y_f")]
  }
  
  # If user wants the distance between cells
  if (dl == TRUE){
    if (dist == 'proj'){
      Edges[, dl := sqrt((x_f - x_i)^2 + (y_f - y_i)^2)]
    } else {
      Edges[, c("long_i","lat_i") :=
              as.data.table((project(as.matrix(data.table(x=x_i,y=y_i)),
                                     from = crs(z_fix), to = "+proj=longlat")))
      ][, c("long_f","lat_f") :=
          as.data.table((project(as.matrix(data.table(x=x_f,y=y_f)),
                                 from = crs(z_fix), to = "+proj=longlat")))]
      if (dist == 'karney'){
        Edges$dl <- geosphere::distGeo(Edges[,.(long_f,lat_f)],Edges[,.(long_i,lat_i)], 
                                       a = r, f = f)
      } else if (dist == 'cosine'){
        Edges$dl <- geosphere::distCosine(Edges[,.(long_f,lat_f)],Edges[,.(long_i,lat_i)],
                                          r = r) 
      } else if (dist == 'haversine'){
        Edges$dl <- geosphere::distHaversine(Edges[,.(long_f,lat_f)],Edges[,.(long_i,lat_i)],
                                             r = r) 
      } else if (dist == 'meeus'){
        Edges$dl <- geosphere::distMeeus(Edges[,.(long_f,lat_f)],Edges[,.(long_i,lat_i)],
                                         a = r, f = f) 
      } else if (dist == 'vincentyEllipsoid'){
        Edges$dl <- geosphere::distVincentyEllipsoid(Edges[,.(long_f,lat_f)],Edges[,.(long_i,lat_i)],
                                                     a = r, b = b, f = f)
      } else if (dist == 'vincentySphere'){
        Edges$dl <- geosphere::distVincentySphere(Edges[,.(long_f,lat_f)],Edges[,.(long_i,lat_i)],
                                                  r = r)
      } else{
        stop("Unknown distance calculation method")
      }
      Edges[, c("long_i","lat_i","long_f","lat_f") := NULL]
    }
    
    if (xy1 == FALSE){
      Edges$x_i <- NULL
      Edges$y_i <- NULL
    }
    if (xy2 == FALSE){
      Edges$x_f <- NULL
      Edges$y_f <- NULL
    }
  }
  
  return(Edges[])
}



#' A function that for a given world of possible movement calculates
#' the transition cost for each in terms of caloric work, caloric energy,
#' and time.
#'
#' @title Calculate movement costs
#' @param world Either the output of the \code{\link[lbmech]{importWorld}} function,
#' or a an object of class SpatialPolygonsDataFrame representing
#' the partitioning grid for the maximum possible area, in the same format as the
#' output of the \code{\link[lbmech]{makeGrid}} function.
#' @param method one of either \code{'kuo'}, \code{'heglund'}, or \code{'oscillator'} defining
#' the method by which to calculate work per stride; see details below.
#' @param m The mass of the animal moving across the landscape, in kilograms.
#' @param v_max The maximum velocity of the animal moving across the landscape,
#' in meters per second; see \code{\link[lbmech]{getVelocity}}.
#' @param epsilon The biomechanical efficiency factor for an animal moving across
#' the landscape. Default is \code{epsilon = 0.2}.
#' @param BMR The base metabolic rate of the object moving across the landscape
#' in Joules per second.
#' @param k The topographic sensitivity factor; see \code{\link[lbmech]{getVelocity}}.
#' @param alpha The dimensionless slope of maximum velocity;
#' see \code{\link[lbmech]{getVelocity}}
#' @param g The acceleration due to gravity, in meters per second per second.
#' Default is \code{g = 9.81} m/s^2, as for the surface of planet Earth.
#' @param l_s The average stride length, in meters. Required for
#' \code{method =  'kuo'} or \code{'oscillator'}, ignored for \code{'heglund'}
#' @param L The average leg length. Required for \code{method =  'kuo'},
#' ignored for \code{'heglund'} and \code{'oscillator'}.
#' @param gamma The fractional maximal deviation from average velocity per stride.
#' Required for \code{method = 'oscillator'}, ignored for \code{'kuo'} and \code{'heglund'}.
#' @param dist A character string specifying the method by which to calculate
#' linear distances. One of \code{proj} (for source projection distance units; the fastest and default),
#' \code{karney} (the second most accurate and second fastest but by far best choice), 
#' \code{cosine}, \code{haversine}, \code{meeus}, 
#' \code{vincentyEllipsoid} or \code{vincentySphere} (see documentation for the
#' \code{\link[geosphere]{geosphere}} package for further documentation). If
#' \code{dist != 'proj'}, one of \code{z_fix}, or \code{proj} and \code{res}
#' must be provided. Recall default \code{res = 5} in other functions but here
#' it is \code{NULL}.
#' @param z_fix A SpatRaster with the same origin and resolution as the
#' \code{z_fix} used to generate the 'world' with \code{\link[lbmech]{makeWorld}}.
#' @param proj The projection used to generate the world. Required if 
#' \code{dist != 'proj'} and \code{z_fix == NULL}, otherwise ignored.
#' @param res The res used to generate the world. Required if 
#' \code{dist != 'proj'} and \code{z_fix == NULL}, otherwise ignored.
#' @param r The earth's radius, or semimajor axis (a) for ellipsoidal models. 
#' Default is 6378137 for WGS1984. Ignored unless \code{dist != 'proj'}.
#' @param f The ellipsoidal flattening of the earth for ellipsoidal models.
#' Default is \code{f=1/298.257223563} for WGS1984.  Ignored unless
#' \code{dist == c('karney','meeus','vincentyEllipsoid')}.
#' @param b The earth's semiminor axis. 
#' Default is \code{b=6356752.3142} for WGS1984.  Ignored unless  
#' \code{dist == 'vincentyEllipsoid'}.
#' @param ... Additional arguments to pass to \code{\link[lbmech]{importWorld}}.
#' @return A data.table with thirteen columns, representing:
#'
#' (1) \code{$from} The "x,y"-format coordinates of each possible origin cell
#'
#' (2) \code{$to} The "x,y"-format coordinates of each possible destination cell
#'
#' (3) \code{$dz} The change in elevation between the \code{from} and \code{to}
#' cells
#'
#' (4) \code{$x_i} The numeric x coordinate of the origin cell
#'
#' (5) \code{$y_i} The numeric y coordinate of the origin cell
#'
#' (6) \code{$dl} The planimetric distance between the \code{from} and \code{to} cells
#' 
#' (7) \code{$dr} The 3D distance between the \code{from} and \code{to} cells
#'
#' (8) \code{$dl_t} The predicted walking speed in meters per second
#' when walking between the \code{from} and \code{to} cells
#'
#' (9) \code{$dt} The predicted amount of time spent walking between
#' the \code{from} and \code{to} cells
#'
#' (10) \code{$dU_l} The predicted work against gravitational potential energy
#' in Joules when walking between the \code{from} and \code{to} cells
#'
#' (11) \code{$dK_l} The predicted kinematic work in Joules when walking
#' between the \code{from} and \code{to} cells
#'
#' (12) \code{$dW_l} The total predicted energy lost due to biomechanical
#' work when walking between the \code{from} and \code{to} cells.
#'
#' (13) \code{$dE_l} The net metabolic expenditure exerted when walking
#' between the \code{from} and \code{to} cells.
#' @importFrom data.table transpose
#' @importFrom data.table :=
#' @examples 
#' 
#' #### Example 1:
#' # Running importWorld before calculating costs
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
#' grid <- makeGrid(dem = dem, nx = n, ny = n)
#' 
#' # Import the data lying between x = (12000,16000) and y = (32000,36000)
#' region <- ext(c(12000,16000,32000,36000))
#' region <- as.polygons(region)
#' crs(region) <- crs(grid)
#' 
#' world <- importWorld(region = region, polys = grid,
#'                      cut_slope = 0.5, res = res(dem), dir = dir)
#'
#' # We'll run calculateCosts using the canonical parameters for Tobler's 
#' # function on a human being. Note dist = 'proj' since the coordinates
#' # are somewhat nonsensical. If real-world scenario, employ a geodesic method
#' world <- calculateCosts(world = world, m = 60, v_max = 1.5,
#'                         BMR = 93, k = 3, alpha = -0.05, l_s = 1.6, L = 0.8,
#'                         dist = 'proj')
#'                         
#'                         
#' #### Example 2:
#' # Running all of the previous steps within calculateCosts
#' # We'll first get velocity parameters from the getVelocity function
#' data <- data.table(x = runif(10000,10000,20000),
#'                    y = runif(10000,30000,40000),
#'                    elevation = runif(10000,0,200),
#'                    dt = 120,
#'                    ID = rep(1:10,each=1000))
#' velocity <- getVelocity(data = data, z = 'elevation')
#' 
#' # Run calculateCosts
#' world <- calculateCosts(world = grid, m = 60, v_max = velocity$vmax,
#'                         k = velocity$k, alpha = velocity$alpha, BMR = 93,
#'                         l_s = 1.6, L = 0.8, region = region,
#'                         proj = crs(dem), res = res(dem),
#'                         cut_slope = 0.5, dir = dir, dist = 'proj')
#' 
#' @export
calculateCosts <- function(world, method = 'kuo', m = NULL, v_max = NULL,
                           epsilon = 0.2, BMR = NULL, k = NULL, alpha = NULL,
                           g = 9.81, l_s = NULL, L = NULL, gamma = NULL, 
                           dist = 'proj', z_fix = NULL, proj = NULL, res = NULL,
                           r = 6378137, f = 1/298.257223563, b = 6356752.3142,
                           ...){
  # This bit is to silence the CRAN check warnings for literal column names
  dz=from=to=dl=x_f=x_i=y_f=y_i=dl_t=..v_max=..k=..alpha=dt=dU_l=..m=..g=NULL
  dK_l=..l_s=..L=..gamma=dW_l=..epsilon=dE_l=..BMR=long_f=lat_f=long_i=lat_i=dr=NULL
  
  # Create z_fix if missing
  if (is.null(z_fix) & dist != 'proj'){
    z_fix <- fix_z(proj = proj, ...)
  }
  
  
  
  # If what we provided was a sector tile grid, make sure the appropriate
  # files have been imported/created
  if("SpatialPolygonsDataFrame" %in% class(world)){
    world <- vect(world)
  }
  if("SpatVector" %in% class(world)){
    world <- importWorld(polys = world, ...)
  }
  
  # Calculate the transition costs. This is broken up into three
  # principal transformations.
  ## (1) Convert the cell names into coordinates, calculate distances,
  ##     calculate velocity, calculate time spent in cell, and calculate
  ##     work against gravity
  world[, `:=`(dz = as.numeric(dz))
  ][, c("x_i","y_i") := transpose(stringr::str_split(from,","))
  ][, c("x_f","y_f") := transpose(stringr::str_split(to,","))
  ][, c("x_i","y_i","x_f","y_f") := lapply(.SD,as.numeric),
    .SDcols = c("x_i","y_i","x_f","y_f")
  ]
  if (dist == 'proj'){
    world[, dl := sqrt((x_f - x_i)^2 + (y_f - y_i)^2)]
  } else {
    world[, c("long_i","lat_i") :=
            as.data.table((project(as.matrix(data.table(x=x_i,y=y_i)),
                                   from = crs(z_fix), to = "+proj=longlat")))
    ][, c("long_f","lat_f") :=
        as.data.table((project(as.matrix(data.table(x=x_f,y=y_f)),
                               from = crs(z_fix), to = "+proj=longlat")))]
    if (dist == 'karney'){
      world$dl <- geosphere::distGeo(world[,.(long_f,lat_f)],world[,.(long_i,lat_i)], 
                                     a = r, f = f)
    } else if (dist == 'cosine'){
      world$dl <- geosphere::distCosine(world[,.(long_f,lat_f)],world[,.(long_i,lat_i)],
                                        r = r) 
    } else if (dist == 'haversine'){
      world$dl <- geosphere::distHaversine(world[,.(long_f,lat_f)],world[,.(long_i,lat_i)],
                                           r = r) 
    } else if (dist == 'meeus'){
      world$dl <- geosphere::distMeeus(world[,.(long_f,lat_f)],world[,.(long_i,lat_i)],
                                       a = r, f = f) 
    } else if (dist == 'vincentyEllipsoid'){
      world$dl <- geosphere::distVincentyEllipsoid(world[,.(long_f,lat_f)],world[,.(long_i,lat_i)],
                                                   a = r, b = b, f = f)
    } else if (dist == 'vincentySphere'){
      world$dl <- geosphere::distVincentySphere(world[,.(long_f,lat_f)],world[,.(long_i,lat_i)],
                                                r = r)
    } else{
      stop("Unknown distance calculation method")
    }
    world[, c("long_i","lat_i","long_f","lat_f") := NULL]
  }
  world[, dr := sqrt(dl^2 + dz^2)
  ][, `:=`(x_f = NULL, y_f = NULL)
  ][, dl_t := ..v_max * exp(-(..k) * abs(dz/dl - ..alpha))
  ][, dt := dl_t ^ -1 * dl #* ..l_p / ..l_s
  ][, dU_l := ..m * ..g * dz
  ][dU_l < 0, dU_l := 0]
  
  ## (2) Calculate the work based on a user-selected function
  if (method == 'kuo'){
    # Kuo's function for human movement
    world[, dK_l := 1 / 4 * ..m * dl_t^2 * ..l_s * dl / ..L^2]
  } else if (method == 'heglund'){
    # Heglund et al.'s function for arbitrary quadripeds
    world[, dK_l := (0.478 * dl_t ^ 1.53 + 0.685 * dl_t + 0.072) * dl_t ^ -1 * dl * ..m]
  } else if (method == 'oscillator'){
    world[, dK_l := 2 * ..m * dl_t ^2 * ..gamma * dl / ..l_s]
  }
  
  ## (3) Finally, calculate the total work and energy
  world[, dW_l := (dU_l + dK_l) / ..epsilon
  ][, dE_l := dW_l + ..BMR * dt][]
  
  return(world)
}


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
  
  z_temp <- as.data.table(expand.grid(x = seq(from = min(data[,..x]) - 2 * res(z_fix)[1],
                                to = max(data[,..x]) + 2 * res(z_fix)[1],
                                by = res(z_fix)[1]),
                        y = seq(from = min(data[,..y]) - 2 * res(z_fix)[2],
                                to = max(data[,..y]) + 2 * res(z_fix)[2],
                                by = res(z_fix)[2])))
  z_temp <- rast(z_temp[,.(x,y,z=1)])
  crs(z_temp) <- crs(z_fix)
  z_fix <- suppressWarnings(project(z_temp,z_fix, align= TRUE))
  
  for (i in seq(1,nrow(data))){
    centroid <- data[i, .(x = get(..x),y = get(..y))]
    poi <- cellFromXY(z_fix, centroid[,.(x,y)])
    poi <- paste(format(round(xyFromCell(z_fix, poi),precision), 
                        scientific=FALSE), 
                 collapse = ",")
    data[i, Cell := stringr::str_remove_all(..poi," ")]
  }
  
  
  return(data$Cell)
}

#' A function that calculates the cost to travel between two sets of points.
#' This can be between two circumscribed sets of points, or one circumscribed one
#' ('nodes') and all other points on the landscape.
#'
#' There are four possible workflows:
#'
#' (1) If you simply desire the distance between two sets of points, provide
#' entries for \code{from} and \code{to} (or just \code{from} if the interest is
#' in all distances between locations in that object). Output is a distance matrix.
#' The computational time for this operation is comparable to generating a raster
#' for the distance to \emph{all} cells in the world (unless all of the locations
#' in the object are close to each other). So unless the operation is to be done
#' multiple times, it is highly recommended to generate the rasters as below and extract
#' values
#'
#' (2) If you wish to generate a RasterStack of costs from and/or to all nodes
#' in the \code{from} object, set the \code{output = 'object'} and 
#' \code{destination = 'all'}.
#'
#' (3) You may also save the rasters as a series of \code{.tif} files in the same workspace
#' directory as the transition \code{.gz} tensor files and the cropped/downloaded
#' DEMs. This allows us to use \code{getCosts} within a loop for large numbers of origin
#' nodes without running into random access memory limitations. Do this by
#' setting \code{output = 'file'} and \code{destination = 'all'}.
#'
#' (4) You may perform (2) and (3) simultaneously by setting
#' \code{output == c('file','object')} and \code{destination = 'all'}.
#'
#' @title Get cost of travel
#' @param world The data.table output of the \code{\link[lbmech]{calculateCosts}} 
#' function.
#' @param from One of data.frame, data.table, SpatVector, SpatialPointsDataFrame, or
#' SpatialPolygonsDataFrame representing the origin locations. If \code{to = NULL}
#' and \code{destination = 'pairs'}, this will also be used as the \code{to} parameter.
#' If object is of class SpatialPolygonsDataFrame, the location will be taken at the
#' centroid.
#' @param proj A crs object or character string representing the workspace 
#' projection. Required, unless \code{z_fix} is provided in which case
#' it's ignored.
#' @param to One of data.frame, data.table, SpatVector SpatialPointsDataFrame, or
#' SpatialPolygonsDataFrame representing the origin locations. Optional if
#' \code{destination =  'pairs'}, ignored if \code{destination =  'all'}.
#' @param id Character string representing the column containing the unique
#' IDs for for each location in the \code{from} (and \code{to}) objects.
#' @param x A character vector representing the column containing the 'x' coordinates.
#' Required if \code{data} is not Spatial*.
#' @param y A character vector representing the column containing the 'y' coordinates.
#' Required if \code{data} is not Spatial*.
#' @param destination One of \code{'pairs'} or \code{'all'}. If \code{'pairs'},
#' a distance matrix will be generated between every pair of locations in
#' \code{from}, or every pair of locations between \code{from} and \code{to}.
#' If \code{'all'}, rasters will be generated for each node representing the cost
#' to travel to every cell in the given \code{world}.
#' @param region An object of class SpatRaster, SpatVector, Raster* or SpatialPolygons*
#' representing the total area where movement is possible. Optional; used
#' to constrain the \code{world} more than it may already have been.
#' @param costs A character vector containing any combination of the strings
#' \code{c("time","work","energy")} if the input world data.table is the
#' output of of the \code{\link[lbmech]{calculateCosts}} function. Otherwise
#' a character string or vector containing the names of the columns where the 
#' cost values are stored in the custom-produced data.table.
#' This selects which types of costs will be calculated.
#' \code{costs = 'all'} is shorthand for \code{costs = c("time","work","energy")}
#' while \code{costs = 'energetics'} is shorthand for \code{c("work","energy")}.
#' Default is \code{'all'}.
#' @param direction A character vector containing one or both of \code{c("in","out")}
#' or the singular string 'both'. This determines whether costs to or from the nodes
#' are calculated. Ignored for \code{destination = 'pairs'}.
#' @param output A character vector containing one or both of \code{c("object","file")}.
#' If \code{"object"} is included, then a list of RasterStacks will be returned.
#' If \code{"file"} is included, then the appropriate cost rasters will be saved
#' to the workspace directory \code{dir}. Ignored if \code{destination = 'pairs'}.
#' @param z_fix A SpatVector with the same origin and resolution as the
#' \code{z_fix} used to generate the 'world' with \code{\link[lbmech]{makeWorld}}.
#' Do not modify this parameter if you didn't modify it when you ran
#' \code{\link[lbmech]{makeWorld}}.
#' @param res The res used to generate the world. Required if 
#' \code{dist != 'proj'} and \code{z_fix == NULL}, otherwise ignored.
#' @param dir A filepath to the directory being used as the workspace.
#' Default is \code{tempdir()} but unless the analyses will only be performed a few
#' times it is highly recommended to define a permanent workspace.
#' @param precision An integer representing the number of decimals to retain
#' in the x and y directions. For grid sizes with nice, round numbers precisions
#' can be low. This factor is controled by \code{\link[terra]{rast}} and
#' must be the same as the one used to generate the 
#' 'world' with \code{\link[lbmech]{makeWorld}}. Default is 2.
#' @param ... Additional arguments to pass to \code{\link[lbmech]{fix_z}}.
#' @return A list, with entries for each combination of selected \code{costs}
#' and \code{directions}. The object class of the list entries depends on the
#' \code{destination} and \code{output} parameters:
#'
#' (1) If \code{destination = 'pairs'}, entries are distance matrices.
#' (2) If \code{destination = 'all'} and \code{'object' \%in\% output},
#' entries are RasterStacks with each Raster* representing the cost to/from
#' each node.
#' (3) If \code{destination = 'all'} and \code{output == 'file'}, no list is output.
#'
#' Moreover, if \code{file \%in\% output}, then the cost rasters are saved in the
#' workspace directory.
#' @importFrom data.table as.data.table
#' @importFrom data.table :=
#' @importFrom data.table data.table
#' @importFrom data.table tstrsplit
#' @importFrom terra writeRaster
#' @importFrom terra rast
#' @importFrom terra crs
#' @importFrom terra crs<-
#' @importFrom terra res
#' @importFrom terra project
#' @importFrom terra geomtype
#' @importFrom terra centroids
#' @importFrom terra geom
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
#' @examples 
#' #### Example 1:
#' # Calculate the 'time' and 'work' cost between a set of points
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
#' grid <- makeGrid(dem = dem, nx = n, ny = n)
#' 
#' # Import the data lying between x = (12000,16000) and y = (32000,36000)
#' region <- ext(c(12000,16000,32000,36000))
#' region <- as.polygons(region)
#' crs(region) <- crs(grid)
#' 
#' # Use the canonical parameters for a human in Tobler's Function
#' world <- calculateCosts(world = grid, m = 60, v_max = 1.5,
#'                         k = 3.5, alpha = -0.05, BMR = 93,
#'                         l_s = 1.6, L = 0.8, region = region, 
#'                         proj = crs(dem), res = res(dem),
#'                         cut_slope = 0.5, dir = dir, dist = 'proj')
#'                         
#' # Generate five random points that fall within the region
#' points <- data.table(ID = 1:5,
#'                      x = runif(5, ext(region)[1], ext(region)[2]),
#'                      y = runif(5, ext(region)[3], ext(region)[4]))
#'                      
#' # Get 'time' and 'work' costs
#' costMatrix <- getCosts(world = world, from = points, 
#'                        costs = c("time","work"),
#'                        proj = crs(dem), res = res(dem), dir = dir)
#'                      
#' #### Example 2:
#' # Calculate the cost rasters to travel to and from a set of points
#' costRasters <- getCosts(world = world, from = points,
#'                         destination = 'all', costs = c('time','work'),
#'                         output = c("object","file"), 
#'                         proj = crs(dem), res = res(dem),
#'                         dir = dir)
#' 
#' @export
getCosts <- function(world, from, to = NULL, proj = NULL, id = 'ID', 
                     x = "x", y = "y", destination = 'pairs', region = NULL, 
                     costs = 'all', direction = 'both', output = 'object', 
                     precision = 2, z_fix = NULL, res = NULL, 
                     dir = tempdir(), ...){
  # This bit is to silence the CRAN check warnings for literal column names
  ..id=..x=..y=dt=dW_l=dE_l=Var2=value=Var1=coord=Replace=Masked=x_n=y_n=NULL
  Cell=From_ID=ID=..cost=..d=Value=x_i=y_i=NULL
  #
  
  dir <- normalizePath(dir)
  # First, "all" and "energetics" are just shorthand for a vector
  # of multiple cost types. Fix that.
  if (all(costs == 'all')){
    costs <- c("time","work","energy")
  } else if (all(costs == "energetics")){
    costs <- c("work","energy")
  }
  
  id2 <- 1
  # If no Unique ID column is given, make one
  if (is.null(id)){
    from$Node_ID <- seq(1,nrow(from))
    id <- "Node_ID"
    id2 <- NULL
  }
  
  if (is.null(z_fix)){
    z_fix <- fix_z(proj = proj, res = res, ...)
  }
  
  if(is.null(world$x_i) | is.null(world$y_i)){
    world[, c("x_i","y_i") := tstrsplit(from,",")
          ][, c("x_i","y_i") := lapply(.SD,as.numeric),.SDcols = c("x_i","y_i")]
  }
  
  z_temp <- as.data.table(expand.grid(x = seq(from = min(world[,x_i]) - 2 * res(z_fix)[1],
                                to = max(world[,x_i]) + 2 * res(z_fix)[1],
                                by = res(z_fix)[1]),
                        y = seq(from = min(world[,y_i]) - 2 * res(z_fix)[2],
                                to = max(world[,y_i]) + 2 * res(z_fix)[2],
                                by = res(z_fix)[2])))
  z_temp <- rast(z_temp[,.(x,y,z=1)])
  crs(z_temp) <- crs(z_fix)
  z_fix <- suppressWarnings(project(z_temp,z_fix,align = TRUE))
  
  # We need to coerce everything to an (x,y) list stored as a data.table
  # First, if the input are polygons we'll just find the centroid
  fromMask <- NULL
  toMask <- NULL
  polysmask <- NULL
  if (methods::is(from,'Spatial')){ 
    from <- vect(from)
  }
  if (methods::is(from, "SpatVector")){
    if (geomtype(from) == 'polygons'){
    fromMask <- regionMask(from, z_fix = z_fix, id = id)
    }
    from <- centroids(from,inside = TRUE)
    from[[x]] <- as.data.table(geom(from))[,x]
    from[[y]] <- as.data.table(geom(from))[,y]
  }
  
  
  # Coerce to data.table
  from <- as.data.table(from)[,.(ID = get(..id), x = get(..x), y = get(..y))]
  from$Cell <- getCoords(from, z_fix = z_fix, precision = precision)
  
  # Do the same for the destination UNLESS it's set as "all".
  # If it's empty, then from is the same as to.
  if (destination != "all"){
    if (!is.null(to)){
      if (methods::is(to,'Spatial')){
      to <- vect(to)
        }
      if (methods::is(to,"SpatVector")){
        # If no column ID is provided, make one as above
        if (is.null(id2)){
          to$Node_ID <- seq(1,nrow(to))
        }
        if (geomtype(from) == 'polygons'){
        toMask <- regionMask(to, z_fix = z_fix, id = id)
        }
        to <- centroids(to,inside = TRUE)
        to[[x]] <- as.data.table(geom(to))[,x]
        to[[y]] <- as.data.table(geom(to))[,y]
      }
      
      to <- as.data.table(to)[,.(ID = get(..id), x = get(..x), y = get(..y))]
      to$Cell <- getCoords(to, z_fix = z_fix, precision = precision)
    } else {
      to <- from
      toMask <- fromMask
    }
  } else{
    # If destination is set as "all", then to is "NULL" so that igraph
    # calculaltes the distance to all cells.
    if (!is.null(to)){
      warning("Destinations have been provided, but method selected is 'all'.
            Ignoring input for 'to'.")
      Sys.sleep(1)
      to <- NULL
    }
  }
  
  if (!is.null(fromMask)){
  fromMask <- merge(fromMask,from,by='ID'
  )[,`:=`(Masked = coord,
          Replace = Cell,
          coord = NULL,
          Cell = NULL)][]
  polysmask <- fromMask
  } 
  if (!is.null(toMask)){
  toMask <- merge(toMask,to,by='ID'
  )[,`:=`(Masked = coord,
          Replace = Cell,
          coord = NULL,
          Cell = NULL)][]
  polysmask <- rbind(polysmask,toMask)
  }
  if (!is.null(polysmask)){
  polysmask[, (id) := NULL]
  polysmask <- unique(polysmask)
  
  
  polysmask[, c("x_n","y_n") := tstrsplit(Replace,",")
  ][, c("x_n","y_n") := lapply(.SD,as.numeric),.SDcols = c("x_n","y_n")
  ]
  
  w <- world[from %in% polysmask$Masked | to %in% polysmask$Masked]
  w[from %in% polysmask$Masked, 
    c("from","x_i","y_i") := polysmask[Masked %in% from, 
                                      .(Replace, x_n, y_n)][1]
  ][to %in% polysmask$Masked, 
    to := polysmask[Masked %in% to, .(Replace)][1]
  ][,(names(world)[!(names(world) %in% c("from","to"))]) :=
      lapply(.SD,min),by=c("from","to")]
  w <- unique(w)
  
  world <- rbind(world[!(from %in% polysmask$Masked | to %in% polysmask$Masked)],
                 w)
  rm(w)
  }
  # If a region further constraining the world more than what was previosly
  # provided is input, then drop those observations from the dataframe
  if (!is.null(region)){
    region <- regionMask(region = region, z_fix = z_fix)
    world <- world[(from %in% region) & (to %in% region),]
  }
  
  # If direction is set to both, we'll need to create a "to"
  # and "from" category
  if (all(direction == "both")){
    direction <- c("in","out")
  }
  
  # Create a list object to which we'll append either the matrices or
  # rasterstacks representing the costs
  outList <- list()
  
  pb <- txtProgressBar(min = 0, max = length(costs) * length(direction), style=3)
  iter = 0
  
  # We first iterate over the desired costs. Each desired cost will
  # create one or two entries (depending on the number of directions)
  # in the output list.
  for (cost in costs){
    # For each loop, we'll generate the igraph object using the
    # appropriate cost variable
    if (all(cost == "time")){
      Graph <- igraph::graph_from_data_frame(world[,.(from,to,weight = dt)])
    } else if (all(cost == "work")){
      Graph <- igraph::graph_from_data_frame(world[,.(from,to,weight = dW_l)])
    } else if (all(cost == "energy")){
      Graph <- igraph::graph_from_data_frame(world[,.(from,to,weight = dE_l)])
    } else {
      Graph <- igraph::graph_from_data_frame(world[,.(from,to,weight = get(cost))])
    }
    
    costVals <- data.table()
    
    # Now we also iterate over each desired direction
    for (d in direction){
      
      # This section is if the desired output is a *matrix*
      # and we are not storing the rasters
      if (all(destination == 'pairs')){
        # Calculate the distances
        Distances <- igraph::distances(Graph,
                                       v = from$Cell,
                                       to = to$Cell,
                                       mode = d)
        
        # Convert to data.table
        Distances <- as.data.table(reshape2::melt(Distances))
        
        # Append the appropriate names/ID
        Distances <- merge(Distances[,.(Var2,value,Cell = Var1)],
                           from[,.(Cell,From_ID = get(..id))],
                           on="Cell")[, Cell := NULL][]
        
        Distances <- merge(Distances[,.(From_ID,value,Cell = Var2)],
                           to[,.(Cell,To_ID = get(..id))],
                           on="Cell")[, Cell := NULL][]
        
        # Convert to distance matrix with xtabs,
        # Add to the output list
        outList[[paste0(cost,"_",d)]] <- stats::xtabs(value~.,Distances)
        
      }
      
      # This section is if the desired costs are from a set of
      # cells to every location on the landscape (a raster of costs)
      if (all(destination == 'all')){
        # Calculate the distances
        Distances <- igraph::distances(Graph,
                                       v = from$Cell,
                                       mode = d)
        
        # Convert to data.table
        Distances <- as.data.table(reshape2::melt(Distances))
        
        # Append the appropriate IDs to each origin/destination node
        # then get the x,y coordinates from the cell name
        Distances <- merge(Distances[,.(Var2,value,Cell = Var1)],
                           from[,.(Cell,ID = get(..id))],
                           on="Cell")[, Cell := NULL
                           ][,.(ID,Cell = Var2, Cost = ..cost,
                                Direction = ..d,Value = value)]
        Distances[, c("x","y") := tstrsplit(Cell,",")
        ][, c("x","y") := lapply(.SD,as.numeric), .SDcols = c("x","y")]
        
        if (!is.null(polysmask)){
          Distances <- rbind(Distances,
            merge(Distances, polysmask[,.(Masked,Replace)],
                             by.x='Cell', by.y='Replace',
                             allow.cartesian=TRUE
                             )[!is.na(Masked), Cell := Masked
                               ][, Masked := NULL
                                 ][, c("x","y") := tstrsplit(Cell,",")
                                 ][, c("x","y") := lapply(.SD,as.numeric), .SDcols = c("x","y")][])
        }
        
        # We'll add all of the generated rasters to a stack
        outStack <- rast()
        for (i in unique(Distances$ID)){
          # Convert the data.table to raster and add to the stack
          outStack <- suppressWarnings(c(outStack,
                               rast(Distances[ID == i,.(x,y,Value)],
                                             crs=crs(z_fix))))
        }
        
        # Select the appropriate prefix names for the filepaths/layer names
        cost <-  stringr::str_to_title(cost)
        if (all(d == "out")){
            names(outStack) <- paste0(cost,"_From_",unique(Distances$ID))
        } else if (all(d == "in")){
            names(outStack) <- paste0(cost,"_To_",unique(Distances$ID))
        }
       
        
        # If we want to output a rasterstack, add it to the list
        if ("object" %in% output){
          outList[[paste0(cost,"_",d)]] <- outStack
        }
        
        # If we want to export the file, save it
        if ("file" %in% output){
          dir.create(normalizePath(paste0(dir,"/CostRasters/"),mustWork=FALSE),
                     showWarnings=FALSE)
          writeRaster(outStack,
                      filename=normalizePath(paste0(dir,"/CostRasters/",
                                                    names(outStack),".tif"),
                                             mustWork=FALSE),
                      overwrite=TRUE)
        }
        
        
      }
      iter <- iter + 1
      setTxtProgressBar(pb,iter)
    }
    
  }
  
  if (methods::is(outList[[1]], 'SpatRaster')){
    if (length(outList) > 1){
    names(outList) <- NULL
    }
    outList <- rast(outList)
    outList <- outList[[gtools::mixedsort(names(outList))]]
  } else {
    if (length(outList) == 1){
      outList <- unlist(outList)
    }
  }
  return(outList)
}

#' A function to automatically perform the raster arithmetic
#' necessary to calculate the cost-of-travel for paths with multiple
#' waypoints, and the predicted cost of taking a detour to any
#' arbitrary point in the landscape (a 'corridor').
#' \code{\link[lbmech]{getCosts}} must have been run before this tool can be used.
#'
#' @title Calculate cost corridors for a path
#' @param rasters One of either a character string or multilayer SpatRaster. 
#' If character string, it represents the filepath
#' to the workspace used as \code{dir} for the previous functions.
#' Default is \code{tempdir()} but unless you are not following best
#' practices you will have to change it to your output directory. If multilayer
#' SpatRaster, it should be the output (or identical in form) to
#' the \code{\link[lbmech]{getCosts}} function with \code{"object" \%in\% output}.
#' @param order A character vector containing the desired path in
#' order of visited nodes. For example, to visit "A" then "B" then "C" then "A"
#' the vector would be \code{c("A","B","C","A")}. Note that these MUST correspond
#' to the ID names for the \code{from} features used in the \code{\link[lbmech]{getCosts}}
#' function and must have previously been calculated
#' @param costs A character vector containing any combination of the strings
#' \code{c("time","work","energy")} if the input world data.table is the
#' output of of the \code{\link[lbmech]{calculateCosts}} function. Otherwise
#' a character string or vector containing the names of the columns where the 
#' cost values are stored in the custom-produced data.table.
#' This selects which types of costs will be calculated. 
#' \code{costs = 'all'} is shorthand for \code{costs = c("time","work","energy")}
#' while \code{costs = 'energetics'} is shorthand for \code{c("work","energy")}.
#' Default is \code{'all'}. Note that these must have previously been calculated.
#' @return Rasters representing cost corridors.
#' If \code{length(costs) == 1}, a Raster* If \code{length(costs) > 1}
#' a list of Raster* with one slot for each \code{cost}.
#' @importFrom terra rast
#' @importFrom terra global
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
#' grid <- makeGrid(dem = dem, nx = n, ny = n)
#' 
#' # Import the data lying between x = (12000,16000) and y = (32000,36000)
#' region <- ext(c(12000,16000,32000,36000))
#' region <- as.polygons(region)
#' crs(region) <- crs(grid)
#' 
#' 
#' # Use the canonical parameters for a human in Tobler's Function
#' world <- calculateCosts(world = grid, m = 60, v_max = 1.5,
#'                         k = 3.5, alpha = -0.05, BMR = 93,
#'                         l_s = 1.6, L = 0.8, region = region, 
#'                         proj = crs(dem), res = res(dem),
#'                         cut_slope = 0.5, dir = dir)
#'                         
#' # Generate five random points that fall within the region
#' points <- data.table(ID = 1:5,
#'                      x = runif(5, ext(region)[1], ext(region)[2]),
#'                      y = runif(5, ext(region)[3], ext(region)[4]))
#'                      
#' # Calculate cost rasters
#' costRasters <- getCosts(world = world, from = points, 
#'                         destination = 'all', costs = 'all',
#'                         output = c("object","file"), 
#'                         proj = crs(dem), res = res(dem), dir = dir)
#'                         
#' #### Example 1:
#' # Calculating the corridors from a list of RasterStacks,
#' # with path 1 -> 2 -> 4 -> 1 -> 5
#' corridors <- makeCorridor(rasters = costRasters, order = c(1,2,5,1,4))
#' 
#' #### Example 2:
#' # Calculating the corridors from a workspace directory
#' # with path 1 -> 2 -> 4 -> 1 -> 5
#' corridors <- makeCorridor(rasters = dir, order = c(1,2,5,1,4))
#' @export
makeCorridor <- function(rasters = tempdir(), order, costs = "all"){
  # Unlike getCosts, makeCorridor only works if the previous steps have
  # already been completed. It needs either the output list of cost
  # rasters from getCosts OR the directory containing the output
  # cost raster files for all desired nodes.
  if (all(costs == 'all')){
    costs <- c("time","work","energy")
  } else if (all(costs == "energetics")){
    costs <- c("work","energy")
  }
  
  costs <- stringr::str_to_sentence(costs)
  # Get the names of files/rasters that will be needed for each leg
  starts <- paste0("_From_",order[seq(1,(length(order)-1))])
  stops <-  paste0("_To_",order[seq(2,length(order))])
  
  # Empty output list, as in getCosts
  outList <- rast()
  for (cost in costs){
    # Iterate over costs, again as in getCosts
    
    if (methods::is(rasters, "SpatRaster")){
      # If the input is a SpatRaster, then we can just
      # get the appropriate rasters from the keys
      from <- rasters[[paste0(cost,starts)]]
      to <- rasters[[paste0(cost,stops)]]
      
    } else if (methods::is(rasters,'character')){
      # If the input is a directory with rasters, then we need to import
      # each individual raster and stack them
      rd <- normalizePath(paste0(rasters,"/CostRasters"),mustWork=FALSE)
      
      from <- normalizePath(paste0(rd,"/",
                                   stringr::str_to_sentence(cost),
                                   starts,".tif"),mustWork = FALSE)
      to <- normalizePath(paste0(rd,"/",
                                 stringr::str_to_sentence(cost),
                                 stops,".tif"),mustWork =FALSE)
      
      from <- rast(from)
      to <- rast(to)
    }
    
    # The cost for each leg is From Origin + To Destination
    corridor <- from + to
    mins <- unlist(global(corridor,'min',na.rm=TRUE))
    names(mins) <- NULL
    
    # To ensure that the path from the first origin to the last destination
    # contains the same minimum travel cost, subtract the minimum cost per
    # leg from each leg's raster
    corridor <- corridor - mins
    
    # The total corridor is the minimum value at each cell from each
    # leg's raster. Add the sum of the minimums to represent the total
    # minimum cost, not just detour.
    corridor <- min(corridor,na.rm=TRUE) + sum(mins)
    names(corridor) <- cost
    
    # Add to the output list
    outList <- suppressWarnings(c(outList, corridor))
  }
  return(outList)
}


#' Get the shortest path for a given trip that requires travel through a
#' set of nodes. Use is like \code{\link[lbmech]{getCosts}}, but with 
#' \code{nodes} and {order} parameters and no \code{from} or \code{to}.
#'
#' @title Get least-cost paths
#' @param world The data.table output of the \code{\link[lbmech]{calculateCosts}} function.
#' @param nodes An object of class data.frame, data.table, SpatialPointsDataFrame, or
#' SpatialPolygonsDataFrame representing the locations.
#' If object is of class SpatialPolygonsDataFrame, the location will be taken at the
#' centroid.
#' @param proj A crs object or character string representing the workspace 
#' projection. Required, unless \code{z_fix} is provided in which case
#' it's ignored.
#' @param id A character string representing the column containing each \code{node}
#' location's unique ID.
#' @param order A character vector containing the desired path in
#' order of visited nodes by ID. For example, to visit "A" then "B" then "C" then "A"
#' the vector would be \code{c("A","B","C","A")}. If this is not provided,
#' the function assumes that \code{nodes} is already sorted in the desired
#' order.
#' @param z_fix A SpatRaster with the same origin and resolution as the
#' \code{z_fix} used to generate the 'world' with \code{\link[lbmech]{makeWorld}}.
#' @param x A character vector representing the column containing the 'x' coordinates.
#' Required if \code{data} is not Spatial*.
#' @param y A character vector representing the column containing the 'y' coordinates.
#' Required if \code{data} is not Spatial*.
#' @param region An object of class Raster* or SpatialPolygons*
#' representing the total area where movement is possible. Optional; used
#' to constrain the \code{world} more than it may already have been.
#' @param costs A character vector containing any combination of the strings
#' \code{c("time","work","energy")} if the input world data.table is the
#' output of of the \code{\link[lbmech]{calculateCosts}} function. Otherwise
#' a character string or vector containing the names of the columns where the 
#' cost values are stored in the custom-produced data.table.
#' This selects which types of costs will be calculated.
#' \code{costs = 'all'} is shorthand for \code{costs = c("time","work","energy")}
#' while \code{costs = 'energetics'} is shorthand for \code{c("work","energy")}.
#' Default is \code{'all'}.
#' @param precision An integer representing the number of decimals to retain
#' in the x and y directions. For grid sizes with nice, round numbers precisions
#' can be low. This factor is controled by \code{\link[terra]{rast}} and
#' must be the same as the one used to generate the 
#' 'world' with \code{\link[lbmech]{makeWorld}}. Default is 2.
#' @param z_fix A SpatRaster with the same origin and resolution as the
#' \code{z_fix} used to generate the 'world' with \code{\link[lbmech]{makeWorld}}.
#' Do not modify this parameter if you didn't modify it when you ran
#' \code{\link[lbmech]{makeWorld}}.
#' @param ... Additional arguments to pass to \code{\link[lbmech]{fix_z}}.
#' @return SpatialLinesDataFrames representing least-cost paths. For each
#' cost, each entry within the SpatialLinesDataFrame object represents
#' a single leg of the journey, sorted in the original path order.
#' If \code{length(costs) == 1}, only a SpatialLinesDataFrame is returned.
#' If \code{length(costs) > 1}
#' a list of SpatialLinesDataFrame with one slot for each \code{cost} is returned.
#' @importFrom data.table as.data.table
#' @importFrom data.table tstrsplit
#' @importFrom data.table data.table
#' @importFrom data.table rbindlist
#' @importFrom data.table :=
#' @importFrom terra vect
#' @importFrom terra as.lines
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
#' # Import raster, get the grid
#' dem <- rast(paste0(dir,"/DEM.tif"))
#' grid <- makeGrid(dem = dem, nx = n, ny = n)
#' 
#' # Import the data lying between x = (12000,16000) and y = (32000,36000)
#' region <- ext(c(12000,16000,32000,36000))
#' region <- as.polygons(region)
#' crs(region) <- crs(grid)
#' 
#' 
#' # Use the canonical parameters for a human in Tobler's Function
#' world <- calculateCosts(world = grid, m = 60, v_max = 1.5,
#'                         k = 3.5, alpha = -0.05, BMR = 93,
#'                         l_s = 1.6, L = 0.8, region = region,
#'                         proj = crs(dem), res = res(dem),
#'                         cut_slope = 0.5, dir = dir)
#'                         
#' # Generate five random points that fall within the region
#' points <- data.table(ID = 1:5,
#'                      x = runif(5, ext(region)[1], ext(region)[2]),
#'                      y = runif(5, ext(region)[3], ext(region)[4]))
#'                      
#'                      
#' # Calculate the path from 1 -> 2 -> 5 -> 1 -> 4 
#' pathOrder <- c(1,2,5,1,4)
#' 
#' paths <- getPaths(world = world, nodes = points, order = pathOrder,
#'                           proj = crs(dem), res = res(dem),
#'                           costs = c('time','work','energy'))
#'                         
#' #getCosts(world = world, from = points, proj = crs(dem), res = res(dem),
#' #                        destination = 'all', costs = 'all',
#' #                        output = 'file', dir = dir)                       
#' #corridors <- makeCorridor(rasters = dir, order = pathOrder)
#' #plot(corridors$time)
#' #plot(paths$time,add=TRUE)
#' @export
getPaths <- function(world, nodes, proj = NULL, id = "ID", order = NULL, x = "x",
                     y = "y", region = NULL, costs = 'all', precision = 2, 
                     z_fix = NULL, ...){
  # This bit is to silence the CRAN check warnings for literal column names
  from=to=..id=..x=..y=ID=Cell=V1=V2=TempID=dt=dW_l=dE_l=x_i=y_i=NULL
  coord=Replace=Masked=x_n=y_n=NULL
  #
  
  # First, "all" and "energetics" are just shorthand for a vector
  # of multiple cost types. Fix that.
  if (all(costs == 'all')){
    costs <- c("time","work","energy")
  } else if (all(costs == "energetics")){
    costs <- c("work","energy")
  }
  
  # If no z_fix is provided, make one. Here we require crs
  if (is.null(z_fix)){
    z_fix <- fix_z(proj = proj, ...)
  }
  
  z_temp <- as.data.table(expand.grid(x = seq(from = min(world[,x_i]) - 2 * res(z_fix)[1],
                                to = max(world[,x_i]) + 2 * res(z_fix)[1],
                                by = res(z_fix)[1]),
                        y = seq(from = min(world[,y_i]) - 2 * res(z_fix)[2],
                                to = max(world[,y_i]) + 2 * res(z_fix)[2],
                                by = res(z_fix)[2])))
  z_temp <- rast(z_temp[,.(x,y,z=1)])
  crs(z_temp) <- crs(z_fix)
  z_fix <- suppressWarnings(project(z_temp, z_fix, align = TRUE))
  
  # If region is provided, filter out those values
  if (!is.null(region)){
    region <- regionMask(region, z_fix = z_fix, ...)
    world <- world[(from %in% region) & (to %in% region),]
  }

  # Deal with polygons such that all values within are considered a single cell
  nodeMask <- NULL
  if (methods::is(nodes,'Spatial')){ 
    nodes <- vect(nodes)
  }
  if (methods::is(nodes,"SpatVector")){
    if (geomtype(nodes) == 'polygons'){
    nodeMask <- regionMask(nodes, z_fix = z_fix, id = id)
    }
    nodes <- centroids(nodes,inside = TRUE)
    nodes[[x]] <- as.data.table(geom(nodes))[,x]
    nodes[[y]] <- as.data.table(geom(nodes))[,y]
  }
  
  # Coerce to data.table
  nodes <- as.data.table(nodes)[,.(ID = get(..id), x = get(..x), y = get(..y))]
  nodes$Cell <- getCoords(nodes, z_fix = z_fix, precision = precision, ...)
  
  # Consolidate all cells within the world for any polygon to their centroid
  if (!is.null(nodeMask)){
    nodeMask <- merge(nodeMask,nodes,by='ID'
    )[,`:=`(Masked = coord,
            Replace = Cell,
            coord = NULL,
            Cell = NULL)][]
  
  nodeMask[, (id) := NULL]
  nodeMask <- unique(nodeMask)
  nodeMask[, c("x_n","y_n") := tstrsplit(Replace,",")
  ][, c("x_n","y_n") := lapply(.SD,as.numeric),.SDcols = c("x_n","y_n")
  ]
  
  w <- world[from %in% nodeMask$Masked | to %in% nodeMask$Masked]
  w[from %in% nodeMask$Masked, 
        c("from","x_i","y_i") := nodeMask[Masked %in% from, 
                                          .(Replace, x_n, y_n)][1]
  ][to %in% nodeMask$Masked, 
    to := nodeMask[Masked %in% to, .(Replace)][1]
  ][,(names(world)[!(names(world) %in% c("from","to"))]) :=
      lapply(.SD,min),by=c("from","to")]
  w <- unique(w)
  
  world <- rbind(world[!(from %in% nodeMask$Masked | to %in% nodeMask$Masked)],
                 w)
  rm(w)
  
  }
  
  # If no order is given, assume the node object is already in the desired order
  if (is.null(order)){
    order <- unlist(nodes[[id]])
  }
  
  # Get the IDs of nodes that will be needed for each leg
  starts <- order[seq(1,(length(order)-1))]
  stops <-  order[seq(2,length(order))]
  order <- paste0(starts,"_to_",stops)
  
  # We only need to do the calculations once for each origin node,
  # but we need to keep the order of the above vectors
  origins <- unique(starts)
  
  # We'll store each path inside a list
  pathList <- list()
  
  # Iterate over every cost, adding each result as an entry in pathList
  for (cost in costs){
    if (cost == "time"){
      Graph <- igraph::graph_from_data_frame(world[,.(from,to,weight = dt)])
    } else if (cost == "work"){
      Graph <- igraph::graph_from_data_frame(world[,.(from,to,weight = dW_l)])
    } else if (cost == "energy"){
      Graph <- igraph::graph_from_data_frame(world[,.(from,to,weight = dE_l)])
    } else {
      Graph <- igraph::graph_from_data_frame(world[,.(from,to,weight = get(cost))])
      
    }
    
    # Store each *segment* (leg) of each path in a list,
    # Since they're out-of-order, document the name of the iterated leg
    legList <- list()
    nameList <- c()
    
    # Calculate the shortest path ONCE for each origin node...
    for (i in origins){
      # ...destination node pair storing it in the above list/vector
      dests <- unique(stops[which(starts == i)])
      legs <- igraph::shortest_paths(Graph,from=nodes[ID == i,Cell],
                                     to = nodes[ID %in% dests,Cell],
                                     mode='out', output='vpath')$vpath
      legList <- append(legList, legs)
      nameList <- c(nameList,paste0(i,"_to_",dests))
    }
    
    # Convert the output igraph object to the cell names, convert to
    # a numeric data.table containing XY coordinates for each traversed cell
    legList <- lapply(legList, igraph::as_ids)
    legList <- lapply(legList, tstrsplit, ",")
    legList <- rbindlist(legList,idcol="TempID")
    legList[,`:=`(x = as.numeric(V1), y = as.numeric(V2), V1 = NULL, V2 = NULL)]
    
    # Convert the XY list to SpatialPoints, then SpatialLines
    lineList <- list()
    for (i in unique(legList$TempID)){
      lineList[[i]] <- as.lines(vect(legList,geom=c(x,y),crs=crs(z_fix)))
    }
    lineList <- vect(lineList)
    
    
    # Add the information about what segment each leg represents
    lineList$segment <- nameList
    lineList$cost <- cost
      
    # And place it in the originally-requested order
    lineList <- lineList[match(order, nameList),]
    
    # Save to the output pathList
    pathList[[cost]] <- lineList
    
  }
  if (length(pathList) > 0) {
    if (length(pathList) >= 1 ){
      names(pathList) <- NULL
      pathList <- vect(pathList)
      crs(pathList) <- crs(z_fix)
    }
    return(pathList)
  }
}
