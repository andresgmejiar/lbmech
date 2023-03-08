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
#' ocean bathymetry as does the SRTM data from AWS. Default is \code{z_min = NULL},
#' but set to \code{0} for SRTM data.
#' @param filt Numeric. Size of moving window to apply a low-pass filter.
#'  Default is \code{filt = 0}.
#' @param verbose Should the number of remaining sectors be printed? Default
#' is \code{FALSE}
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
#' grid <- makeGrid(dem = dem, nx = n, ny = n, sources = TRUE)
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
                   z_min = NULL, filt = 0, verbose = FALSE,
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
        (!("makeGrid" %in% names(polys)))){
      # If what we provide is a polygon with URLs to the source,
      # download the file
      for (i in seq(1,length(down))){
        tile_name <- as.character(polys[down[i],][[tile_id]])
        if (!verbose){
          print(paste0("Downloading Tile ",tile_name," (",
                       i," of ",length(down),")"))
        }
        if(!('extension' %in% names(polys))){
          extension <- unlist(stringr::str_split(polys[down[i],][[vals]],"\\."))
          extension <- extension[length(extension)]
        }
        file_path <- normalizePath(paste0(rd,"/",tile_name),mustWork = FALSE)
        utils::download.file(url=unlist(polys[down[i],][[vals]]),
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
            files <- do.call(merge,files) 
          } else {
            files <- rast(files) 
          }
          if (!is.null(z_min)){
            files[files < z_min] <- NA
          }
          if (filt != 0){
            files <- focal(files,w=filt,fun=mean,na.policy='omit')
          }
          names(files) <- 'z'
          
          writeRST(files, normalizePath(file_path))
          unlink(paste0(file_path),recursive=TRUE)
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
          if (!is.null(z_min)){
            files[files < z_min] <- NA
          }
          if (filt != 0){
            files <- focal(files,w=filt,fun=mean,na.policy='omit')
          }
          writeRST(files, normalizePath(file_path))
          unlink(normalizePath(paste0(rd,"/",tile_name),
                               mustWork = FALSE),recursive=TRUE)
        } 
        
        # If downloaded file is an uncompressed file
        exts <- unlist(stringr::str_extract_all(exts,pattern="[a-z]+"))
        if (extension %in% exts){
          files <- rast(paste0(file_path,".",extension)) 
          if (!is.null(z_min)){
            files[files < z_min] <- NA
          }
          if (filt != 0){
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
        clip <- suppressWarnings(rast(elevatr::get_elev_raster(sf::st_as_sf(
          methods::as(polys[down[i],],'Spatial')),z = zoom,src = 'aws')))
        poly <- which(polys[[tile_id]] == tile_name)
        
        if (!is.null(z_min)){
          clip[clip < z_min] <- NA
        }
        if (filt != 0){
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
      
      if (!is.null(z_min)){
        dem[dem < z_min] <- NA
      }
      if (filt != 0){
        dem <- focal(dem,w=filt,fun=mean,na.policy='omit')
      }
      # For every tile that needs to be acquired...
      for (i in seq(1,length(down))){
        
        # Get the unique tile id, and define the output filepath
        tile_name <- as.character(polys[down[i],][[tile_id]])
        file_path <- normalizePath(paste0(rd,"/",tile_name),mustWork=FALSE)
        if (!verbose){
          print(paste0("Cropping Tile ",tile_name," (",
                       i," of ",length(down),")"))
        }
        # Select the singular tile, and use it to crop the dem.
        # Save it to the above filepath, zip it, and delete the tiff
        poly <- which(polys[[tile_id]] == tile_name)
        clip <- crop(dem,polys[poly,], snap = 'out') * 1.0
        writeRST(clip,file_path)
      }
    }
  }
}
