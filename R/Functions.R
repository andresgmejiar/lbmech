#' Generate a partitioning grid for a single raster source representing regional
#' elevations. Smaller partitioning grids (i.e. a greater value of \code{nx * ny})
#' results in a greater number of saved files and a greater number of
#' read-write operations in future operations, but reduces the amount of
#' memory employed.
#'
#' @title Make partitioning grid
#' @param dem A raster object containing the elevations for the maximum possible
#' extent imaginable for a study. Note that this only works for rasters that
#' have been read in, not those that exist exclusively in the memory. If you
#' have just generated the raster and it is in memory, export it first with
#' \code{\link[raster]{writeRaster}} then re-import it with \code{\link[raster]{raster}}.
#' @param nx The integer-number of columns in the output grid
#' @param ny The integer-number of rows in the output grid
#' @param crs A crs object representing the desired output projection.
#' Default is the input raster's projection.
#' @param prefix A character string containing the prefix to name individual sectors.
#' Default is \code{prefix = "SECTOR_"}
#' @return An object of class SpatialPolygonsDataFrame representing the individual sectors ('tiles'),
#' with a dataframe containing three columns: the "TILEID", the raster's filepath,
#' and a dummy column indicating that the grid was made using the makeGrid function.
#' This will be necessary for future functions. The object MUST be stored 
#' on the disk, it should not be stored in the memory
#' @importFrom raster extent
#' @importFrom data.table as.data.table
#' @importFrom data.table :=
#' @importFrom sp SpatialPolygonsDataFrame
#' @importFrom sp proj4string
#' @importFrom sp proj4string<-
#' @examples
#' # Generate a DEM, export it
#' n <- 5
#' dem <- expand.grid(list(x = 1:(n * 100),
#'                         y = 1:(n * 100))) / 100
#' dem <- as.data.table(dem)
#' dem[, z := 250 * exp(-(x - n/2)^2) + 
#'       250 * exp(-(y - n/2)^2)]
#' dem <- rasterFromXYZ(dem)
#' extent(dem) <- c(10000, 20000, 30000, 40000)
#' crs(dem) <- "+proj=lcc +lat_1=48 +lat_2=33 +lon_0=-100 +datum=WGS84"
#' 
#' dir <- tempdir()
#' writeRaster(dem, paste0(dir,"/DEM.tif"),format="GTiff",overwrite=TRUE)
#' 
#' 
#' # Import raster, get the grid
#' dem <- raster(paste0(dir,"/DEM.tif"))
#' grid <- makeGrid(dem = dem, nx = n, ny = n)
#' 
#' @export
makeGrid <- function(dem, nx, ny, crs = NULL, prefix = "SECTOR_"){
  # This bit is to silence the CRAN check warnings for literal column names
  x=..dx=y=..dy=left=right=bottom=top=NULL
  #
  
  # If the DEM to be used is a single large DEM,
  # it's necessary to make a grid that will define
  # the extent of each sector.
  
  # The code snippet within was taken from  Jelmer (2016)'s
  # response in the following stackoverflow thread, accessed 27 January 2021:
  # https://stackoverflow.com/questions/29784829/
  #  r-raster-package-split-image-into-multiples
  # It has been modified such that each raster has two rows of overlap (l_p)
  # at the top and bottom to calculate cost values for cells that would be
  # otherwise inaccessible were there no buffer. It is also used to generate
  # a polygon grid and not directly clip rasters.
  # That'll happen in the getMap function.
  
  filepath <- normalizePath(dem@file@name)
  if (!is.null(crs)){
    dem <- extent(dem)
    dem <- methods::as(dem,"SpatialPolygons")
    dem <- sp::spTransform(dem,crs)
  }
  
  
  dx     <- (extent(dem)[2]- extent(dem)[1])/ nx  ## extent of one tile in x
  dy     <- (extent(dem)[4]- extent(dem)[3])/ ny  ## extent of one tile in y
  xs     <- seq(extent(dem)[1], by= dx, length= nx) ## lower left x-coordinate
  ys     <- seq(extent(dem)[3], by= dy, length= ny) ## lower left y-coordinate
  cS     <- expand.grid(x= xs, y= ys)
  
  # Get the eastings/northings of the extent limits
  cS <- as.data.table(cS)
  cS[, `:=`(left = x,
            right = x + ..dx,
            bottom = y,
            top = y + ..dy)]
  ### END CODE SNIPPET ###
  
  # Convert the data.table to an SPDF and assign it a unique TILEID
  polys <- apply(cS[,.(left,right,bottom,top)],
                 1, FUN = function(x) methods::as(extent(x), "SpatialPolygons"))
  polys <- do.call("bind",polys)
  polys <- SpatialPolygonsDataFrame(polys,
                                    data.frame("TILEID" =
                                                 paste0(prefix,1:length(polys))))
  
  # Apply the source projection
  proj4string(polys) <- crs(dem)
  
  # Set the location of the data to the filepath of the original raster
  polys$location <- filepath
  
  # If the output SPDF is input into getMap, this tells getMap that the source
  # DEM is a single file and it should crop; not download.
  polys$makeGrid <- TRUE
  return(polys)
}

#' Get the names of tiles that would be needed to perform an analysis
#' over a given region-of-interest within the maximum possible extent defined
#' by a source grid.
#'
#' @title Identify necessary tiles/sectors
#' @param region An object of class RasterLayer, SpatialPolygons,
#' SpatialPolygonsDataFrame, SpatialPoints, SpatialPointsDataFrame, data.frame,
#' or data.table indicating the region-of-interest. If input is of class
#' SpatialPoints*, data.table, or data.frame the \code{region} represents sectors containing the
#' individual points.
#' @param polys An object of class SpatialPolygonsDataFrame representing
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
#' @importFrom raster extent
#' @importFrom raster crs
#' @importFrom raster crs<-
#' @importFrom raster intersect
#' @importFrom sp SpatialPointsDataFrame
#' @examples 
#' 
#' #### Example 1:
#' # If the grid is the product if the makeGrid function
#' # Make the grid
#' n <- 6
#' dem <- raster(ncol = n * 6000, nrow = n * 6000)
#' extent(dem) <- c(10000, 20000, 30000, 40000)
#' crs(dem) <- "+proj=lcc +lat_1=48 +lat_2=33 +lon_0=-100 +datum=WGS84"
#' grid <- makeGrid(dem = dem, nx = n, ny = n)
#' 
#' # Select five random points that fall within the grid
#' points <- data.table(x = runif(5, extent(dem)[1], extent(dem)[2]),
#'                      y = runif(5, extent(dem)[3], extent(dem)[4]))
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
    region <- SpatialPointsDataFrame(region[,.(get(..x),get(..y))], data=region,
                                     proj4string = crs(polys))
  }
  if (class(region) == "RasterLayer"){
    proj <- crs(region)
    region <- extent(region)
    region <- methods::as(region,"SpatialPolygons")
    crs(region) <- proj
  }
  if ("SpatialPolygonsDataFrame" %in% class(region)){
    region <- methods::as(region,"SpatialPolygons")
  }
  region <- as.data.table(intersect(region,polys)@data)
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
#' @param polys An object of class SpatialPolygonsDataFrame representing
#' the partitioning grid for the maximum possible area, in the same format as the
#' output of the \code{\link[lbmech]{makeGrid}} function.
#' @param tile_id a character string representing the name of the column
#' in the 'polys' polygon containing the unique Tile IDs. Default is \code{tile_id = 'TILEID'}
#' @param vals A character string or a RasterLayer object. Optional if the
#' \code{polys} polygon is the output of the \code{\link[lbmech]{makeGrid}} function as the default is
#' the character string \code{'location'}. If not, the \code{vals} parameter should be
#' set to the column name containing the URL or filepath to the DEM for that
#' sector.
#' @param dir A filepath to the directory being used as the workspace.
#' Default is \code{tempdir()} but unless the analyses will only be performed a few
#' times it is highly recommended to define a permanent workspace.
#' @return Function does not return any objects, but sets up the workspace
#' such that the necessary DEM files are downloaded/cropped and accessible.
#' @importFrom raster raster
#' @importFrom raster crop
#' @importFrom raster writeRaster
#' @examples 
#' # Generate a DEM, export it
#' n <- 5
#' dem <- expand.grid(list(x = 1:(n * 100),
#'                         y = 1:(n * 100))) / 100
#' dem <- as.data.table(dem)
#' dem[, z := 250 * exp(-(x - n/2)^2) + 
#'       250 * exp(-(y - n/2)^2)]
#' dem <- rasterFromXYZ(dem)
#' extent(dem) <- c(10000, 20000, 30000, 40000)
#' crs(dem) <- "+proj=lcc +lat_1=48 +lat_2=33 +lon_0=-100 +datum=WGS84"
#' 
#' dir <- tempdir()
#' writeRaster(dem, paste0(dir,"/DEM.tif"),format="GTiff",overwrite=TRUE)
#' 
#' 
#' # Import raster, get the grid
#' dem <- raster(paste0(dir,"/DEM.tif"))
#' grid <- makeGrid(dem = dem, nx = n, ny = n)
#' 
#' 
#' # Generate five random points that fall within the grid
#' points <- data.table(x = runif(5, extent(dem)[1], extent(dem)[2]),
#'                      y = runif(5, extent(dem)[3], extent(dem)[4]))
#'                
#'                            
#' # Run whichTiles and getMap to prepare appropriate sector files
#' tile_list <- whichTiles(region = points, polys = grid) 
#' getMap(tiles = tile_list, polys = grid, dir = dir)
#' @export
getMap <- function(tiles, polys, tile_id = "TILEID", vals = "location",
                   dir = tempdir()){
  # The function first checks to see if the needed "tiles" sector DEMs exist
  # If not, it downloads the correct sectors or crops them from the master
  # DEM
  dir <- normalizePath(dir)
  # Check to see if elevations directory exists. If not, make it
  rd <- normalizePath(paste0(dir,"/Elevations"),mustWork=FALSE)
  if(!dir.exists(rd)){
    dir.create(rd)
  }
  
  # Check to see if the file exists; if not add to "down" vector
  down <- c()
  
  for (tile_name in tiles){
    if (!file.exists(normalizePath(paste0(rd,"/",tile_name,".gz"),
                                   mustWork=FALSE))){
      down <- c(down,which(polys@data[,tile_id] == tile_name))
    }
  }
  
  # Skip the rest if we have all the needed tiles
  if (!is.null(down)){
    
    # If what we provide is a polygon with URLs to the source,
    # download the file
    if ((class(vals) != "RasterLayer") & is.null(polys$makeGrid)){
      for (i in 1:length(down)){
        tile_name <- tile_id[down[i]]
        print(paste0("Downloading Tile ",tile_name," (",
                     i+1," of ",length(down),")"))
        extension <- stringr::str_split_fixed(polys@data[down[i],vals],"/.")[2]
        file_path <- normalizePath(paste0(rd,tile_name),mustWork = FALSE)
        utils::download.file(url=polys@data[down[i],vals],
                             destfile=paste0(file_path,".",extension))
        
        # For compatibility with the PA vignette--which downloads tiles
        # as zip files, if the download is a zip, unzip and convert to .gz.
        # This will likely fail in Linux systems
        if (extension == 'zip'){
          utils::unzip(paste0(file_path,".zip"), 
                       exdir = paste0(rd,"/",tile_name),
                       junkpaths = TRUE)        
          newFile <- unlist(list.files(paste0(rd,tile_name), full.names = TRUE))
          R.utils::gzip(filename = newFile,
                        destname = paste0(file_path,".gz"),
                        remove = FALSE,
                        overwrite = TRUE)
          unlink(newFile,recursive=TRUE)
        } 
      }
      # If what we provide is a polygon and a RasterLayer or a path to a singular
      # DEM as the source, crop the necessary DEM
    } else{
      if (class(vals) != "RasterLayer"){
        # If dem is not a raster, it's because the source raster is stored
        # elsewhere. Import it
        dem <- raster(unique(polys$location))
      } else {
        if (class(vals == "RasterLayer")){
          dem <- vals
        }
      }
      
      # For every tile that needs to be acquired...
      for (i in 1:length(down)){
        
        # Get the unique tile id, and define the output filepath
        tile_name <- polys@data[down[i],tile_id]
        file_path <- normalizePath(paste0(rd,"/",tile_name),mustWork=FALSE)
        print(paste0("Cropping Tile ",tile_name," (",
                     i," of ",length(down),")"))
        
        # Select the singular tile, and use it to crop the dem.
        # Save it to the above filepath, zip it, and delete the tiff.
        
        
        poly <- which(polys@data[,tile_id] == tile_name)
        clip <- crop(dem,polys[poly,], snap = 'in')
        writeRaster(clip, filename = paste0(file_path,".tif"))
        R.utils::gzip(filename = paste0(file_path,".tif"),
                      destname = paste0(file_path,".gz"),
                      remove = FALSE,
                      overwrite = TRUE)
        unlink(paste0(file_path,".tif"),recursive=TRUE)
      }
    }
  }
}

#' Calculate the velocity function for an animal from \emph{(x,y,z,t)} data
#' such as from GPS collars, assuming a function of form Tobler (see ####)
#'
#' @title Calculate a velocity function from data
#' @param data A data.frame or something coercible to a data.table containing
#' all observations
#' @param x A character string representing the data column containing the 'x'
#' coordinates in meters.
#' @param y A character string representing the data column containing the 'y'
#' coordinates in meters.
#' @param z Either a character string, a RasterLayer, or a SpatialPolygonsDataFrame
#' object. If character string, it represents the data column containing the 'z'
#' coordinates/elevations in meters. If RasterLayer, a DEM containing
#' the elevation in meters. If SpatialPolygonsDataFrame, it must represent the
#' sectors and filepaths/URLs (see the output of \code{\link[lbmech]{makeGrid}})
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
#' \code{alpha_init = -0.1}.
#' @param tile_id a character string representing the name of the column
#' in the \code{z} polygon containing the unique Tile IDs. Ignored if elevations are
#' provided as a column or RasterLayer. Otherwise default is \code{tile_id = 'TILEID'}.
#' @param vals A character string or a RasterLayer object. Required only if the
#' \code{z} parameter is a polygon NOT the output of the \code{\link[lbmech]{makeGrid}} function as the default is
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
#' @importFrom raster cellFromXY
#' @importFrom raster crs
#' @importFrom raster raster
#' @importFrom sp SpatialPointsDataFrame
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
#' # If the data do not contain elevation, and 'z' is a RasterLayer
#' data$z <- NULL
#' 
#' # Generate a DEM
#' n <- 5
#' dem <- expand.grid(list(x = 1:(n * 100),
#'                         y = 1:(n * 100))) / 100
#' dem <- as.data.table(dem)
#' dem[, z := 250 * exp(-(x - n/2)^2) + 
#'       250 * exp(-(y - n/2)^2)]
#' dem <- rasterFromXYZ(dem)
#' extent(dem) <- c(10000, 20000, 30000, 40000)
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
#' writeRaster(dem, paste0(dir,"/DEM.tif"),format="GTiff",overwrite=TRUE)
#' 
#' # Import raster, get the grid
#' dem <- raster(paste0(dir,"/DEM.tif"))
#' grid <- makeGrid(dem = dem, nx = n, ny = n)
#' 
#' velocity <- getVelocity(data = data, z = grid, dir = dir)
#' 
#' @export
getVelocity <- function(data, x = 'x', y ='y', z = 'z', 
                        dt = 'dt', ID = 'ID', tau = NULL, tau_vmax = 0.995,
                        tau_nlrq = 0.95,k_init = 3.5, alpha_init = -0.1,
                        tile_id = "TILEID", vals = "location",
                        dir = tempdir()){
  # This bit is to silence the CRAN check warnings for literal column names
  ..x=..y=..z=..ID=..dt=..z_vals=Order=dl=dx=dy=dz=dl_dt=NULL
  #
  
  # Standardize filepath names
  dir <- normalizePath(dir)
  rd <- normalizePath(paste0(dir,"/Elevations"),mustWork = FALSE)
  
  data <- as.data.table(data)
  
  # If z is a raster layer, its simply extracting the z points from xy locations
  if ("RasterLayer" %in% class(z)){
    dem <- z
    z <- dem[cellFromXY(dem,data[,.(get(..x),get(..y))])]
    data[,z := ..z]
    data <- data[,.(ID = get(..ID),x = get(..x),y = get(..y),z,dt = get(..dt))]
    
    # If z is a character, do nothing; just rename columns
  } else if ("character" %in% class(z)){
    
    data <- data[,.(ID = get(..ID),x = get(..x),y = get(..y),z = get(..z),dt = get(..dt))]
    
    # If z is a polygon defining sector locations...
  } else if ("SpatialPolygonsDataFrame" %in% class(z)){
    
    # Create an FID to keep the observations in order since the process will
    # shuffle them
    data$Order <- 1:nrow(data)
    
    # If z locations have already been provided but the input contains a DEM,
    # warn the user
    if (!is.null(data$z)){
      warning("Data contains a 'z' column which will be ignored")
      data[, z:= NULL]
    }
    
    # Create a SPDF object to detect which tiles will be employed
    data_points <- SpatialPointsDataFrame(data[,.(get(..x),get(..y))], data=data,
                                          proj4string = crs(z))
    data_points <- as.data.table(intersect(data_points,z)@data)
    tiles <- unique(data_points[,get(tile_id)])
    
    # Make sure that the necessary files exist
    getMap(tiles=tiles,polys=z,tile_id=tile_id,vals=vals,dir=dir)
    
    # Create an empty data.table in which to store extracted values. There's a
    # small chance some points will land in two sectors, so we'll have
    # to collect all and consolidate
    data_new <- data.table()
    for (tile in tiles){
      # Unzip the source dem tile-by-tile, import it, and extract the elevation
      R.utils::decompressFile(filename = normalizePath(paste0(rd,"/",tile,".gz"),
                                                       mustWork = FALSE), 
                              destname = normalizePath(paste0(rd,"/",tile,'.tif'),
                                                       mustWork = FALSE),
                              overwrite = TRUE,
                              remove = FALSE,
                              FUN = gzfile,
                              ext = 'gz')
      elevs <- raster(normalizePath(paste0(rd,"/",tile,".tif"),mustWork=FALSE))
      
      d <- data_points[get(tile_id) == tile]
      z_vals <- elevs[cellFromXY(elevs,d[,.(get(..x),get(..y))])]
      d[,z := ..z_vals]
      
      # Add it to the empty data.table, and remove the uncompressed tiff.
      data_new <- rbind(data_new,d)
      unlink(normalizePath(paste0(rd,"/",tile,".tif")),recursive=TRUE)
    }
    
    # Group by all columns that ARE NOT z, setting z equal to the mean of all
    # non-NA values. This makes sure we have one row per observation and lets
    # us handle edge cases where multiple z values were assigned to points on
    # the border between two tiles.
    cols <- names(data_new)
    cols <- cols[cols != "z"]
    data_new <- data_new[,.(z = mean(z,na.rm=TRUE)),by=cols][]
    data <- data_new[order(Order)
    ][,.(ID = get(..ID),x = get(..x),y = get(..y),z,dt = get(..dt))
    ][]
    
  }
  
  # Calculate displacement, then speed and slope
  data[, `:=`(dx = x - data.table::shift(x),
              dy = y - data.table::shift(y),
              dz = z - data.table::shift(z)),
       by = 'ID'
  ][, dl := sqrt(dx^2 + dy^2)
  ][, `:=`(dl_dt = dl / dt,
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
  v_max <- as.numeric(stats::quantile(data[,dl_dt],tau_vmax,na.rm=TRUE))
  data$v_max <- v_max
  
  
  # And obtain the other coefficients through an nlrq of the form proposed
  # by Tobler (exponential decay from an optimal angle)
  velocity <- quantreg::nlrq(dl_dt ~ v_max * exp(-k * abs(dz_dl - alpha)),
                             data = data, tau = tau_vmax, start=list(k=k_init,alpha=-alpha_init))
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
  options(warn = defaultW)
  return(out)
  
}


#' Function that defines the grid that can be traversed--the "world"--as well as the
#' cells that can be accessed from each individual cell. This is the most
#' time-intensive function. It first checks to see if the required transition \code{.gz}
#' files have already been created in the \code{dir} workspace. If not, it checks to see if the
#' DEMs required to generate the transition \code{.gz} files have been downloaded/cropped,
#' and generates the latter if so. If not, it downloads/crops them and proceeds.
#'
#' @title Define which cells are adjacent
#' @param tiles A character vector--such as the output to
#' \code{\link[lbmech]{whichTiles}}---containing the unique tile IDs for sectors that
#' should be in the workspace.
#' @param polys An object of class SpatialPolygonsDataFrame representing
#' the partitioning grid for the maximum possible area, in the same format as the
#' output of the \code{\link[lbmech]{makeGrid}} function.
#' @param tile_id A character string representing the name of the column
#' in the 'polys' polygon containing the unique Tile IDs. Default is \code{tile_id = 'TILEID'}
#' @param cut_slope A number representing the dimensionless maximum slope
#' of ascent/descent.
#' @param z_fix A raster that will define the resolution, origin, and
#' projection information for the entire "world" of possible movement. Note that
#' it does NOT need the same extent.
#' @param directions One of the integers \code{c(4, 8, 16)},
#' the character string \code{'bishop'},or a neighborhood matrix.
#' Default is \code{directions = 16}, implying that all 'knight and one-cell
#' queen moves' are permissible movements on the grid. See \code{\link[raster]{adjacent}}.
#' @param neighbor_distance An integer representing the distance in meters
#' that tiles are buffered. In other words, to ensure that all transitions in the
#' 'world' are recorded, files for each tile will contain a number of observations
#' that fall outside of the tile in other ones. Default is 100 m, but adjust
#' on raster size.
#' @param unit One of \code{c("m", "km", "ft", "mi")}, representing the unit of the DEM.
#' All will be converted to meters, which is the default.
#' @param vals A character string or a RasterLayer object. Ignored unless the
#' \code{polys} parameter is a polygon NOT the output of the \code{\link[lbmech]{makeGrid}}
#' function as the default is the character string \code{'location'},
#' AND the appropriate world \code{.gz} file is NOT
#' present in the workspace directory. In which case it must represent either the
#' original DEM or a character string with the column representing the DEM
#' filepath or URL.
#' @param dir A filepath to the directory being used as the workspace.
#' Default is \code{tempdir()} but unless the analyses will only be performed a few
#' times it is highly recommended to define a permanent workspace.
#' @return A \code{.gz} file for each sector named after its sector id,
#' containing a data.table with three columns:
#'
#' (1) \code{$from}, a character string of all possible origin cells in format "x,y",
#' rounded to the next-lowest integer
#'
#' (2) \code{$to},  a character string of all possible destination cells in format "x,y"
#' rounded to the next-lowest integer
#'
#' (3) \code{$dz}, an integer representing the change in elevation for each origin-destination pair
#' @importFrom raster res
#' @importFrom raster projectRaster
#' @importFrom raster mosaic
#' @importFrom raster crs
#' @importFrom raster crs<-
#' @importFrom raster xFromCell
#' @importFrom raster yFromCell
#' @importFrom raster extent
#' @importFrom raster crop
#' @importFrom raster getValues
#' @importFrom raster adjacent
#' @importFrom raster ncell
#' @importFrom raster intersect
#' @importFrom data.table fwrite
#' @importFrom data.table as.data.table
#' @importFrom data.table :=
#' @importFrom data.table .SD
#' @importFrom sp spTransform
#' @examples 
#' # Generate a DEM
#' n <- 5
#' dem <- expand.grid(list(x = 1:(n * 100),
#'                         y = 1:(n * 100))) / 100
#' dem <- as.data.table(dem)
#' dem[, z := 250 * exp(-(x - n/2)^2) + 
#'       250 * exp(-(y - n/2)^2)]
#' dem <- rasterFromXYZ(dem)
#' extent(dem) <- c(10000, 20000, 30000, 40000)
#' crs(dem) <- "+proj=lcc +lat_1=48 +lat_2=33 +lon_0=-100 +datum=WGS84"
#' 
#' # Export it so it doesn't just exist on the memory
#' dir <- tempdir()
#' writeRaster(dem, paste0(dir,"/DEM.tif"),format="GTiff",overwrite=TRUE)
#' 
#'
#' # Import raster, get the grid
#' dem <- raster(paste0(dir,"/DEM.tif"))
#' grid <- makeGrid(dem = dem, nx = n, ny = n)
#' 
#' # Select all tiles that exist between x = (12000,16000) and y = (32000,36000)
#' tiles <- extent(c(12000,16000,32000,36000))
#' tiles <- as(tiles,"SpatialPolygons")
#' crs(tiles) <- crs(grid)
#' tiles <- whichTiles(region = tiles, polys = grid)
#' 
#' makeWorld(tiles = tiles, polys = grid,
#'           cut_slope = 0.5, z_fix = dem, dir = dir)
#' @export 
makeWorld <- function(tiles,polys,tile_id = 'TILEID',cut_slope,z_fix,
                      directions = 16, neighbor_distance = 100,
                      unit = "m", vals = 'location',
                      dir = tempdir()){
  # This bit is to silence the CRAN check warnings for literal column names
  from=to=..dem=z_f=z_i=x_f=x_i=y_f=y_i=dz=NULL
  #
  
  
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
  
  
  rd <- normalizePath(paste0(dir,"/Tensors"),mustWork=FALSE)
  for (i in tiles){
    if (!file.exists(normalizePath(paste0(rd,"/",i,"_Tensor.gz"),mustWork=FALSE))){
      # Get the appropriate polygon and name; create a temporary folder
      tensors <- normalizePath(paste0(rd,"/",i),mustWork=FALSE)
      poly <- which(polys@data[,tile_id] == i)
      poly <- polys[poly,]
      dir.create(tensors,recursive=TRUE)
      l_p <- res(z_fix)[[1]]
      
      # Get the neighboring sectors and unzip them into temp folder
      poly_buff <- rgeos::gBuffer(poly,byid=TRUE,width=neighbor_distance)
      mzip <- intersect(poly_buff,polys)@data[,paste0(tile_id,".2")]
      
      getMap(mzip, polys = polys, tile_id = tile_id, vals = vals,dir = dir)
      
      mzip <- as.list(mzip)
      
      for (j in mzip){
        R.utils::decompressFile(
          filename = normalizePath(paste0(dir,"/Elevations/",j,".gz"),
                                   mustWork = FALSE), 
          destname = normalizePath(paste0(tensors,"/",j,".tif"),
                                   mustWork=FALSE),
          ext = 'gz',
          FUN = gzfile,
          remove = FALSE,
          overwrite = TRUE
        )
      }
      
      
      # Get list of rasters; put them into the same projection,
      # resample them, then mosaic them
      dem <- list.files(path=tensors,pattern=".tif$",full.names=TRUE)
      dem <- lapply(dem,raster)
      dem_temp <- suppressWarnings(lapply(dem,
                                          projectRaster,
                                          to=dem[[1]],
                                          alignOnly=TRUE))
      dem_temp$fun <- mean
      dem_temp <- do.call(mosaic,dem_temp)
      
      dem <- suppressWarnings(lapply(dem,projectRaster,to=dem_temp))
      dem$fun <- mean
      dem <- do.call(mosaic,dem)
      
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
      poly <- spTransform(poly,crs(dem))
      poly <- extent(poly)
      error <- (res(dem) - res(z_fix))[1] * contiguity
      poly <- c(poly[1], poly[2] + contiguity * l_p / 12 / 2.54 * 100,
                poly[3], poly[4] + contiguity * l_p / 12 / 2.54 * 100)
      poly <- extent(poly)
      poly <- methods::as(poly,"SpatialPolygons")
      crs(poly) <- crs(dem)
      
      # Crop the mosaiced raster by the cropping polygon, project it
      dem <- crop(dem,poly)
      dem_temp <- suppressWarnings(projectRaster(dem,to=z_fix,alignOnly = TRUE))
      dem <- suppressWarnings(projectRaster(dem,to=dem_temp,method='bilinear'))
      
      name <- i
      
      nas <- which(!is.na(getValues(dem)))
      
      # Get pairs of adjacent all adjacent cells; drop those that
      # correspond to NA
      adj <- adjacent(dem,1:ncell(dem),directions=directions,pairs=TRUE,sort=TRUE)
      adj <- as.data.table(adj)
      adj <- adj[from %in% nas | to %in% nas,]
      
      # Calculate the change in elevation between every accessible cell pairs,
      # then drop all values that would require movement over the
      # cut slope.
      adj[, `:=`(z_i = ..dem[from], z_f = ..dem[to],
                 x_i = xFromCell(dem,from), y_i = yFromCell(dem,from),
                 x_f = xFromCell(dem,to), y_f = yFromCell(dem,to))
      ][, `:=`(dz = z_f - z_i)
      ][, (c("x_i","y_i","x_f","y_f")) :=
          lapply(.SD,round,2),
        .SDcols = c("x_i","y_i","x_f","y_f")
      ][, `:=`(from = paste(x_i,y_i,sep=","),
               to = paste(x_f,y_f,sep=","))
      ]
      
      adj <- stats::na.omit(adj)
      
      # This exports the cost tensor so that it can be called later by
      # the lbm.distance tool to calculate paths and catchments
      fwrite(adj[,.(from,to,dz)],
             file = normalizePath(paste0(rd,"/",i,"_Tensor.gz"),mustWork=FALSE))
      unlink(tensors,recursive=TRUE)
    }
  }
}

#' Function that converts raster and SpatialPolygon* objects
#' to a list of cells that fall within such a region.
#'
#' @title Convert RasterLayer or SpatialPolygon* to "x,y"
#' @param region An object of class SpatialPolygons* to convert to "x,y"
#' @param z_fix A RasterLayer with the same origin and resolution as the
#' \code{z_fix} used to generate the 'world' with \code{\link[lbmech]{makeWorld}}.
#' @return A character vector containing all cells that fall in the same
#' location as the input 'region'.
#' @importFrom raster resample
#' @importFrom raster getValues
#' @importFrom raster xFromCell
#' @importFrom raster yFromCell
#' @importFrom data.table data.table
#' @importFrom data.table :=
#' @examples 
#' # Generate a DEM
#' n <- 5
#' dem <- expand.grid(list(x = 1:(n * 100),
#'                         y = 1:(n * 100))) / 100
#' dem <- as.data.table(dem)
#' dem[, z := 250 * exp(-(x - n/2)^2) + 
#'       250 * exp(-(y - n/2)^2)]
#' dem <- rasterFromXYZ(dem)
#' extent(dem) <- c(10000, 20000, 30000, 40000)
#' crs(dem) <- "+proj=lcc +lat_1=48 +lat_2=33 +lon_0=-100 +datum=WGS84"
#' 
#' # Generate a polygon that falls within the DEM's extent
#' region <- extent(c(12500,12600,32500,32700))
#' region <- as(region,"SpatialPolygons")
#' crs(region) <- crs(dem)
#' 
#' maskedCells <- regionMask(region = region, z_fix = dem)
#' 
#' @export
regionMask <- function(region, z_fix){
  # This bit is to silence the CRAN check warnings for literal column names
  x=y=NULL
  #
  
  if (class(region) %in% c("SpatialPolygons","SpatialPolygonsDataFrame")){
    region <- fasterize::fasterize(sf::st_as_sf(region), z_fix)
  } else if (class(region) == "RasterLayer"){
    region <- resample(region, z_fix)
  } else {
    stop("Inappropriate Region Object. Only SpatialPolygons* or Raster allowed")
  }
  region <- which(!is.na(getValues(region)))
  region <- data.table(x = round(xFromCell(z_fix, region),2),
                       y = round(yFromCell(z_fix, region),2))
  region[, `:=`(coord = paste(x,y,sep=','), x = NULL, y = NULL)]
  return(region$coord)
}

#' A function that for a given region imports all cells from the
#' transition \code{.gz} files. If such files have not yet been generated,
#' they can be created by passing along the necessary parameters to this
#' function as with \code{\link[lbmech]{makeWorld}}.
#'
#' @title Import a world where movement is possible
#' @param region An object of class RasterLayer or SpatialPolygons*
#' representing the total area where movement is possible. Must lie within
#' the area defined by \code{polys}
#' @param z_fix A RasterLayer with the same origin and resolution as the
#' \code{z_fix} used to generate the 'world' with \code{\link[lbmech]{makeWorld}}.
#' @param polys An object of class SpatialPolygonsDataFrame representing
#' the partitioning grid for the maximum possible area, in the same format as the
#' output of the \code{\link[lbmech]{makeGrid}} function.
#' @param banned An object of class RasterLayer or SpatialPolygons*
#' representing the total area where movement is \emph{prohibited}. Must lie within
#' the area defined by \code{polys}
#' @param tile_id A character string representing the name of the column
#' in the \code{polys} polygon containing the unique Tile IDs. Default is \code{tile_id = 'TILEID'}
#' @param dir A filepath to the directory being used as the workspace.
#' Default is \code{tempdir()} but unless the analyses will only be performed a few
#' times it is highly recommended to define a permanent workspace.
#' @param ... Additional arguments to pass to \code{\link[lbmech]{makeWorld}}
#' @return An object of class data.table containing three columns:
#'
#' (1) \code{$from}, a character string of all possible origin cells in format "x,y",
#' rounded to the next-lowest integer
#'
#' (2) \code{$to},  a character string of all possible destination cells in format "x,y"
#' rounded to the next-lowest integer
#'
#' (3) \code{$dz}, an integer representing the change in elevation for each origin-destination pair
#' @importFrom data.table data.table
#' @importFrom data.table fread
#' @examples 
#' 
#' # Generate a DEM
#' n <- 5
#' dem <- expand.grid(list(x = 1:(n * 100),
#'                         y = 1:(n * 100))) / 100
#' dem <- as.data.table(dem)
#' dem[, z := 250 * exp(-(x - n/2)^2) + 
#'       250 * exp(-(y - n/2)^2)]
#' dem <- rasterFromXYZ(dem)
#' extent(dem) <- c(10000, 20000, 30000, 40000)
#' crs(dem) <- "+proj=lcc +lat_1=48 +lat_2=33 +lon_0=-100 +datum=WGS84"
#' 
#' # Export it so it doesn't just exist on the memory
#' dir <- tempdir()
#' writeRaster(dem, paste0(dir,"/DEM.tif"),format="GTiff",overwrite=TRUE)
#' 
#'
#' # Import raster, get the grid
#' dem <- raster(paste0(dir,"/DEM.tif"))
#' grid <- makeGrid(dem = dem, nx = n, ny = n)
#' 
#' # Import the data lying between x = (12000,16000) and y = (32000,36000)
#' region <- extent(c(12000,16000,32000,36000))
#' region <- as(region,"SpatialPolygons")
#' crs(region) <- crs(grid)
#' 
#' world <- importWorld(region = region, polys = grid,
#'                      cut_slope = 0.5, z_fix = dem, dir = dir)
#' @export
importWorld <- function(region, z_fix, polys, banned = NULL,
                        tile_id = "TILEID", dir = tempdir(), ...){
  # This bit is to silence the CRAN check warnings for literal column names
  from=to=NULL
  #
  
  # First determine which tiles to import,
  # download maps and make tensors if needed
  dir <- normalizePath(dir)
  
  tiles <- whichTiles(region = region, polys = polys,
                      tile_id = tile_id)
  makeWorld(tiles = tiles, polys = polys, tile_id = tile_id, dir = dir,
            z_fix = z_fix, ...)
  
  tiles <- normalizePath(paste0(dir,"/Tensors/",tiles,"_Tensor.gz"),
                         mustWork=FALSE)
  
  # Convert the regions into a list of the allowable cells
  region <- regionMask(region, z_fix)
  if (!is.null(banned)){
    banned <- regionMask(banned, z_fix)
  } else{
    banned <- c("NULL")
  }
  
  # Create an empty data.table to which to add the imported tensors
  Edges <- data.table()
  
  pb <- utils::txtProgressBar(max = length(tiles), style = 3)
  
  # Import tensors one-by-one, keeping only the allowable cells
  for (i in 1:length(tiles)) {
    import <- fread(tiles[i])
    import <- import[(!(from %in% banned) | !(to %in% banned)) &
                       ((from %in% region) | (to %in% region)),]
    Edges <- rbind(Edges, import)
    utils::setTxtProgressBar(pb,i)
  }
  return(Edges)
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
#' @param ... Additional arguments to pass to \code{\link[lbmech]{importWorld}}.
#' @return A data.table with twelve columns, representing:
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
#' (6) \code{$dl} The distance between the \code{from} and \code{to} cells
#'
#' (7) \code{$dl_t} The predicted walking speed in meters per second
#' when walking between the \code{from} and \code{to} cells
#'
#' (8) \code{$dt} The predicted amount of time spent walking between
#' the \code{from} and \code{to} cells
#'
#' (9) \code{$dU_l} The predicted work against gravitational potential energy
#' in Joules when walking between the \code{from} and \code{to} cells
#'
#' (10) \code{$dK_l} The predicted kinematic work in Joules when walking
#' between the \code{from} and \code{to} cells
#'
#' (11) \code{$dW_l} The total predicted energy lost due to biomechanical
#' work when walking between the \code{from} and \code{to} cells.
#'
#' (12) \code{$dE_l} The net metabolic expenditure exerted when walking
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
#' dem <- rasterFromXYZ(dem)
#' extent(dem) <- c(10000, 20000, 30000, 40000)
#' crs(dem) <- "+proj=lcc +lat_1=48 +lat_2=33 +lon_0=-100 +datum=WGS84"
#' 
#' # Export it so it doesn't just exist on the memory
#' dir <- tempdir()
#' writeRaster(dem, paste0(dir,"/DEM.tif"),format="GTiff",overwrite=TRUE)
#' 
#'
#' # Import raster, get the grid
#' dem <- raster(paste0(dir,"/DEM.tif"))
#' grid <- makeGrid(dem = dem, nx = n, ny = n)
#' 
#' # Import the data lying between x = (12000,16000) and y = (32000,36000)
#' region <- extent(c(12000,16000,32000,36000))
#' region <- as(region,"SpatialPolygons")
#' crs(region) <- crs(grid)
#' 
#' world <- importWorld(region = region, polys = grid,
#'                      cut_slope = 0.5, z_fix = dem, dir = dir)
#'
#' # We'll run calculateCosts using the canonical parameters for Tobler's 
#' # function on a human being. 
#' world <- calculateCosts(world = world, m = 60, v_max = 1.5,
#'                         BMR = 93, k = 3, alpha = -0.05, l_s = 1.6, L = 0.8)
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
#'                         l_s = 1.6, L = 0.8, region = region, z_fix = dem,
#'                         cut_slope = 0.5, dir = dir)
#' 
#' @export
calculateCosts <- function(world, method = 'kuo', m = NULL, v_max = NULL,
                           epsilon = 0.2, BMR = NULL, k = NULL, alpha = NULL,
                           g = 9.81, l_s = NULL, L = NULL, gamma = NULL,...){
  # This bit is to silence the CRAN check warnings for literal column names
  dz=from=to=dl=x_f=x_i=y_f=y_i=dl_t=..v_max=..k=..alpha=dt=dU_l=..m=..g=NULL
  dK_l=..l_s=..L=..gamma=dW_l=..epsilon=dE_l=..BMR=NULL
  #
  
  # If what we provided was a sector tile grid, make sure the appropriate
  # files have been imported/created
  if("SpatialPolygonsDataFrame" %in% class(world)){
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
  ][, dl := sqrt((x_f - x_i)^2 + (y_f - y_i)^2)
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
#' containing the coordinates needing conversion, or a SpatialPointsDataFrame
#' @param x A character vector representing the column containing the 'x' coordinates.
#' Required if \code{data} is not SpatialPointsDataFrame.
#' @param y A character vector representing the column containing the 'y' coordinates.
#' Required if \code{data} is not SpatialPointsDataFrame.
#' @param z_fix A RasterLayer with the same origin and resolution as the
#' \code{z_fix} used to generate the 'world' with \code{\link[lbmech]{makeWorld}}.
#' @return A vector containing the requested coordinates in appropriate format
#' in the same order as the input data.
#' @importFrom data.table :=
#' @importFrom data.table as.data.table
#' @importFrom raster cellFromXY
#' @importFrom raster xyFromCell
#' @examples
#' # Generate a DEM
#' n <- 5
#' dem <- expand.grid(list(x = 1:(n * 100),
#'                         y = 1:(n * 100))) / 100
#' dem <- as.data.table(dem)
#' dem[, z := 250 * exp(-(x - n/2)^2) + 
#'       250 * exp(-(y - n/2)^2)]
#' dem <- rasterFromXYZ(dem)
#' extent(dem) <- c(10000, 20000, 30000, 40000)
#' crs(dem) <- "+proj=lcc +lat_1=48 +lat_2=33 +lon_0=-100 +datum=WGS84"
#' 
#' # Generate five random points that fall within the DEM
#' points <- data.table(x = runif(5, extent(dem)[1], extent(dem)[2]),
#'                      y = runif(5, extent(dem)[3], extent(dem)[4]))
#' 
#' # Get the coordinates
#' points$Cell <- getCoords(points, z_fix = dem)
#' @export
getCoords <- function(data, x = "x", y = "y", z_fix){
  # This bit is to silence the CRAN check warnings for literal column names
  ..x=..y=Cell=..poi=NULL
  # End
  
  data <- as.data.table(data)
  for (i in 1:nrow(data)){
    centroid <- data[i, .(x = get(..x),y = get(..y))]
    poi <- cellFromXY(z_fix, centroid[,.(x,y)])
    poi <- paste(round(xyFromCell(z_fix, poi),2), collapse = ",")
    data[i, Cell := ..poi]
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
#' @param world The data.table output of the \code{\link[lbmech]{calculateCosts}} function.
#' @param from One of data.frame, data.table, SpatialPointsDataFrame, or
#' SpatialPolygonsDataFrame representing the origin locations. If \code{to = NULL}
#' and \code{destination = 'pairs'}, this will also be used as the \code{to} parameter.
#' If object is of class SpatialPolygonsDataFrame, the location will be taken at the
#' centroid.
#' @param to One of data.frame, data.table, SpatialPointsDataFrame, or
#' SpatialPolygonsDataFrame representing the origin locations. Optional if
#' \code{destination =  'pairs'}, ignored if \code{destination =  'all'}.
#' @param id Character string representing the column containing the unique
#' IDs for for each location in the \code{from} (and \code{to}) objects.
#' @param z_fix A RasterLayer with the same origin and resolution as the
#' \code{z_fix} used to generate the 'world' with \code{\link[lbmech]{makeWorld}}.
#' @param x A character vector representing the column containing the 'x' coordinates.
#' Required if \code{data} is not Spatial*.
#' @param y A character vector representing the column containing the 'y' coordinates.
#' Required if \code{data} is not Spatial*.
#' @param destination One of \code{'pairs'} or \code{'all'}. If \code{'pairs'},
#' a distance matrix will be generated between every pair of locations in
#' \code{from}, or every pair of locations between \code{from} and \code{to}.
#' If \code{'all'}, rasters will be generated for each node representing the cost
#' to travel to every cell in the given \code{world}.
#' @param region An object of class RasterLayer or SpatialPolygons*
#' representing the total area where movement is possible. Optional; used
#' to constrain the \code{world} more than it may already have been.
#' @param costs A character vector containing any combination of the strings
#' \code{c("time","work","energy")}. This selects which types of costs will be calculated.
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
#' @param dir A filepath to the directory being used as the workspace.
#' Default is \code{tempdir()} but unless the analyses will only be performed a few
#' times it is highly recommended to define a permanent workspace.
#' @return A list, with entries for each combination of selected \code{costs}
#' and \code{directions}. The object class of the list entries depends on the
#' \code{destination} and \code{output} parameters:
#'
#' (1) If \code{destination = 'pairs'}, entries are distance matrices.
#' (2) If \code{destination = 'all'} and \code{'object' \%in\% output},
#' entries are RasterStacks with each RasterLayer representing the cost to/from
#' each node.
#' (3) If \code{destination = 'all'} and \code{output == 'file'}, no list is output.
#'
#' Moreover, if \code{file \%in\% output}, then the cost rasters are saved in the
#' workspace directory.
#' @importFrom data.table as.data.table
#' @importFrom data.table :=
#' @importFrom data.table data.table
#' @importFrom data.table tstrsplit
#' @importFrom raster writeRaster
#' @importFrom raster rasterFromXYZ
#' @importFrom raster stack
#' @importFrom raster addLayer
#' @importFrom raster crs
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
#' dem <- rasterFromXYZ(dem)
#' extent(dem) <- c(10000, 20000, 30000, 40000)
#' crs(dem) <- "+proj=lcc +lat_1=48 +lat_2=33 +lon_0=-100 +datum=WGS84"
#' 
#' # Export it so it doesn't just exist on the memory
#' dir <- tempdir()
#' writeRaster(dem, paste0(dir,"/DEM.tif"),format="GTiff",overwrite=TRUE)
#' 
#'
#' # Import raster, get the grid
#' dem <- raster(paste0(dir,"/DEM.tif"))
#' grid <- makeGrid(dem = dem, nx = n, ny = n)
#' 
#' # Import the data lying between x = (12000,16000) and y = (32000,36000)
#' region <- extent(c(12000,16000,32000,36000))
#' region <- as(region,"SpatialPolygons")
#' crs(region) <- crs(grid)
#' 
#' 
#' # Use the canonical parameters for a human in Tobler's Function
#' world <- calculateCosts(world = grid, m = 60, v_max = 1.5,
#'                         k = 3.5, alpha = -0.05, BMR = 93,
#'                         l_s = 1.6, L = 0.8, region = region, z_fix = dem,
#'                         cut_slope = 0.5, dir = dir)
#'                         
#' # Generate five random points that fall within the region
#' points <- data.table(ID = 1:5,
#'                      x = runif(5, extent(region)[1], extent(region)[2]),
#'                      y = runif(5, extent(region)[3], extent(region)[4]))
#'                      
#' # Get 'time' and 'work' costs
#' costMatrix <- getCosts(world = world, from = points, z_fix = dem,
#'                        costs = c("time","work"), dir = dir)
#'                      
#' #### Example 2:
#' # Calculate the cost rasters to travel to and from a set of points
#' costRasters <- getCosts(world = world, from = points, z_fix = dem,
#'                         destination = 'all', costs = c('time','work'),
#'                         output = c("object","file"), dir = dir)
#' 
#' @export
getCosts <- function(world, from, to = NULL, id = 'ID', z_fix = NULL, x = "x", y = "y",
                     destination = 'pairs', region = NULL, costs = 'all',
                     direction = 'both', output = 'object',
                     dir = tempdir()){
  # This bit is to silence the CRAN check warnings for literal column names
  ..id=..x=..y=dt=dW_l=dE_l=Var2=value=Var1=Cell=From_ID=ID=..cost=..d=Value=NULL
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
    from$Node_ID <- 1:nrow(from)
    id <- "Node_ID"
    id2 <- NULL
  }
  
  # We need to coerce everything to an (x,y) list stored as a data.table
  # First, if the input are polygons we'll just find the centroid
  if (sum(class(from) %in% 
          c("SpatialPolygons","SpatialPolygonsDataFrame")) > 1){
    from <- rgeos::gCentroid(from, byid = TRUE, id = id)
  }
  # Coerce to data.table
  from <- as.data.table(from)[,.(ID = get(..id), x = get(..x), y = get(..y))]
  from$Cell <- getCoords(from, z_fix = z_fix)
  
  # Do the same for the destination UNLESS it's set as "all".
  # If it's empty, then from is the same as to.
  if (destination != "all"){
    if (!is.null(to)){
      if (class(to) %in% c("SpatialPolygons","SpatialPolygonsDataFrame")){
        # If no column ID is provided, make one as above
        if (is.null(id2)){
          to$Node_ID <- 1:nrow(to)
        }
        to <- rgeos::gCentroid(to, byid = TRUE, id = id)
      }
      
      to <- as.data.table(to)[,.(ID = get(..id), x = get(..x), y = get(..y))]
      to$Cell <- getCoords(to, z_fix = z_fix)
    } else {
      to <- from
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
  
  # If a region further constraining the world more than what was previosly
  # provided is input, then drop those observations from the dataframe
  if (!is.null(region)){
    region <- regionMask(region, z_fix = z_fix)
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
    }
    if (all(cost == "work")){
      Graph <- igraph::graph_from_data_frame(world[,.(from,to,weight = dW_l)])
    }
    if (all(cost == "energy")){
      Graph <- igraph::graph_from_data_frame(world[,.(from,to,weight = dE_l)])
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
        
        # We'll add all of the generated rasters to a stack
        outStack <- stack()
        for (i in unique(Distances$ID)){
          # Convert the data.table to raster and add to the stack
          outStack <- addLayer(outStack,
                               rasterFromXYZ(Distances[ID == i,.(x,y,z=Value)],
                                             crs=crs(z_fix)))
        }
        
        # Select the appropriate prefix names for the filepaths/layer names
        if (all(d == "out")){
          names(outStack) <- paste0("From_",unique(Distances$ID))
        } else if (all(d == "in")){
          names(outStack) <- paste0("To_",unique(Distances$ID))
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
                                                    stringr::str_to_sentence(cost),"_",
                                                    names(outStack),".tif"),mustWork=FALSE),
                      bylayer=TRUE,format="GTiff",overwrite=TRUE)
        }
        
        
      }
      iter <- iter + 1
      setTxtProgressBar(pb,iter)
    }
    
  }
  # Return the output list if it exists
  if (length(outList) > 0) {
    if (length(outList) == 1 ){
      outList <- unlist(outList)
    }
    return(outList)
  }
  
}

#' A function to automatically perform the raster arithmetic
#' necessary to calculate the cost-of-travel for paths with multiple
#' waypoints, and the predicted cost of taking a detour to any
#' arbitrary point in the landscape (a 'corridor').
#' \code{\link[lbmech]{getCosts}} must have been run before this tool can be used.
#'
#' @title Calculate cost corridors for a path
#' @param rasters One of either a character string or list of
#' RasterStacks. If character string, it represents the filepath
#' to the workspace used as \code{dir} for the previous functions.
#' Default is \code{tempdir()} but unless you are not following best
#' practices you will have to change it to your output directory. If list
#' of RasterStacks, it should be the output (or identical in form) to
#' the \code{\link[lbmech]{getCosts}} function with \code{"object" \%in\% output}.
#' @param order A character vector containing the desired path in
#' order of visited nodes. For example, to visit "A" then "B" then "C" then "A"
#' the vector would be \code{c("A","B","C","A")}. Note that these MUST correspond
#' to the ID names for the \code{from} features used in the \code{\link[lbmech]{getCosts}}
#' function and must have previously been calculated
#' @param costs A character vector containing any combination of the strings
#' \code{c("time","work","energy")}. This selects which types of costs will be calculated.
#' \code{costs = 'all'} is shorthand for \code{costs = c("time","work","energy")}
#' while \code{costs = 'energetics'} is shorthand for \code{c("work","energy")}.
#' Default is \code{'all'}. Note that these must have previously been calculated.
#' @return Rasters representing cost corridors.
#' If \code{length(costs) == 1}, a RasterLayer. If \code{length(costs) > 1}
#' a list of RasterLayers with one slot for each \code{cost}.
#' @importFrom raster raster
#' @importFrom raster stack
#' @importFrom raster minValue
#' @examples 
#' # Generate a DEM
#' n <- 5
#' dem <- expand.grid(list(x = 1:(n * 100),
#'                         y = 1:(n * 100))) / 100
#' dem <- as.data.table(dem)
#' dem[, z := 250 * exp(-(x - n/2)^2) + 
#'       250 * exp(-(y - n/2)^2)]
#' dem <- rasterFromXYZ(dem)
#' extent(dem) <- c(10000, 20000, 30000, 40000)
#' crs(dem) <- "+proj=lcc +lat_1=48 +lat_2=33 +lon_0=-100 +datum=WGS84"
#' 
#' # Export it so it doesn't just exist on the memory
#' dir <- tempdir()
#' writeRaster(dem, paste0(dir,"/DEM.tif"),format="GTiff",overwrite=TRUE)
#' 
#'
#' # Import raster, get the grid
#' dem <- raster(paste0(dir,"/DEM.tif"))
#' grid <- makeGrid(dem = dem, nx = n, ny = n)
#' 
#' # Import the data lying between x = (12000,16000) and y = (32000,36000)
#' region <- extent(c(12000,16000,32000,36000))
#' region <- as(region,"SpatialPolygons")
#' crs(region) <- crs(grid)
#' 
#' 
#' # Use the canonical parameters for a human in Tobler's Function
#' world <- calculateCosts(world = grid, m = 60, v_max = 1.5,
#'                         k = 3.5, alpha = -0.05, BMR = 93,
#'                         l_s = 1.6, L = 0.8, region = region, z_fix = dem,
#'                         cut_slope = 0.5, dir = dir)
#'                         
#' # Generate five random points that fall within the region
#' points <- data.table(ID = 1:5,
#'                      x = runif(5, extent(region)[1], extent(region)[2]),
#'                      y = runif(5, extent(region)[3], extent(region)[4]))
#'                      
#' # Calculate cost rasters
#' costRasters <- getCosts(world = world, from = points, z_fix = dem,
#'                         destination = 'all', costs = 'all',
#'                         output = c("object","file"), dir = dir)
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
  
  # Get the names of files/RasterLayers that will be needed for each leg
  starts <- paste0("From_",order[1:(length(order)-1)])
  stops <-  paste0("To_",order[2:length(order)])
  
  # Empty output list, as in getCosts
  outList <- list()
  for (cost in costs){
    # Iterate over costs, again as in getCosts
    
    if (class(rasters) == "list"){
      # If the input is a list of RasterStacks, then we can just
      # get the appropriate rasters from the keys
      from <- rasters[[paste0(cost,"_out")]][[starts]]
      to <- rasters[[paste0(cost,"_in")]][[stops]]
    } else if (class(rasters) == 'character'){
      # If the input is a directory with RasterLayers, then we need to import
      # each individual raster and stack them
      rd <- normalizePath(paste0(rasters,"/CostRasters"),mustWork=FALSE)
      
      from <- normalizePath(paste0(rd,"/",
                                   stringr::str_to_sentence(cost),"_",
                                   starts,".tif"))
      to <- normalizePath(paste0(rd,"/",
                                 stringr::str_to_sentence(cost),"_",
                                 stops,".tif"))
      
      from <- lapply(from,raster)
      from <- stack(unlist(from))
      
      to <- lapply(to,raster)
      to <- stack(unlist(to))
      
    }
    
    # The cost for each leg is From Origin + To Destination
    corridor <- from + to
    mins <- minValue(corridor)
    
    # To ensure that the path from the first origin to the last destination
    # contains the same minimum travel cost, subtract the minimum cost per
    # leg from each leg's raster
    corridor <- corridor - mins
    
    # The total corridor is the minimum value at each cell from each
    # leg's raster. Add the sum of the minimums to represent the total
    # minimum cost, not just detour.
    corridor <- min(corridor) + sum(mins)
    
    # Add to the output list
    outList[[cost]] <- corridor
  }
  if (length(outList) > 0) {
    if (length(outList) == 1 ){
      outList <- outList[[cost]]
    }
    return(outList)
  }
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
#' @param id A character string representing the column containing each \code{node}
#' location's unique ID.
#' @param order A character vector containing the desired path in
#' order of visited nodes by ID. For example, to visit "A" then "B" then "C" then "A"
#' the vector would be \code{c("A","B","C","A")}. If this is not provided,
#' the function assumes that \code{nodes} is already sorted in the desired
#' order.
#' @param z_fix A RasterLayer with the same origin and resolution as the
#' \code{z_fix} used to generate the 'world' with \code{\link[lbmech]{makeWorld}}.
#' @param x A character vector representing the column containing the 'x' coordinates.
#' Required if \code{data} is not Spatial*.
#' @param y A character vector representing the column containing the 'y' coordinates.
#' Required if \code{data} is not Spatial*.
#' @param region An object of class RasterLayer or SpatialPolygons*
#' representing the total area where movement is possible. Optional; used
#' to constrain the \code{world} more than it may already have been.
#' @param costs A character vector containing any combination of the strings
#' \code{c("time","work","energy")}. This selects which types of costs will be calculated.
#' \code{costs = 'all'} is shorthand for \code{costs = c("time","work","energy")}
#' while \code{costs = 'energetics'} is shorthand for \code{c("work","energy")}.
#' Default is \code{'all'}.
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
#' @importFrom sp SpatialPoints
#' @importFrom sp SpatialLinesDataFrame
#' @examples 
#' # Generate a DEM
#' n <- 5
#' dem <- expand.grid(list(x = 1:(n * 100),
#'                         y = 1:(n * 100))) / 100
#' dem <- as.data.table(dem)
#' dem[, z := 250 * exp(-(x - n/2)^2) + 
#'       250 * exp(-(y - n/2)^2)]
#' dem <- rasterFromXYZ(dem)
#' extent(dem) <- c(10000, 20000, 30000, 40000)
#' crs(dem) <- "+proj=lcc +lat_1=48 +lat_2=33 +lon_0=-100 +datum=WGS84"
#' 
#' # Export it so it doesn't just exist on the memory
#' dir <- tempdir()
#' writeRaster(dem, paste0(dir,"/DEM.tif"),format="GTiff",overwrite=TRUE)
#' 
#'
#' # Import raster, get the grid
#' dem <- raster(paste0(dir,"/DEM.tif"))
#' grid <- makeGrid(dem = dem, nx = n, ny = n)
#' 
#' # Import the data lying between x = (12000,16000) and y = (32000,36000)
#' region <- extent(c(12000,16000,32000,36000))
#' region <- as(region,"SpatialPolygons")
#' crs(region) <- crs(grid)
#' 
#' 
#' # Use the canonical parameters for a human in Tobler's Function
#' world <- calculateCosts(world = grid, m = 60, v_max = 1.5,
#'                         k = 3.5, alpha = -0.05, BMR = 93,
#'                         l_s = 1.6, L = 0.8, region = region, z_fix = dem,
#'                         cut_slope = 0.5, dir = dir)
#'                         
#' # Generate five random points that fall within the region
#' points <- data.table(ID = 1:5,
#'                      x = runif(5, extent(region)[1], extent(region)[2]),
#'                      y = runif(5, extent(region)[3], extent(region)[4]))
#'                      
#'                      
#' # Calculate the path from 1 -> 2 -> 5 -> 1 -> 4 
#' pathOrder <- c(1,2,5,1,4)
#' 
#' paths <- getPaths(world = world, nodes = points, z_fix = dem,
#'                           order = pathOrder,
#'                           costs = c('time','work','energy'))
#'                         
#' #getCosts(world = world, from = points, z_fix = dem,
#' #                        destination = 'all', costs = 'all',
#' #                        output = 'file', dir = dir)                       
#' #corridors <- makeCorridor(rasters = dir, order = pathOrder)
#' #plot(corridors$time)
#' #plot(paths$time,add=TRUE)
#' @export
getPaths <- function(world, nodes, z_fix, id = "ID", order = NULL, x = "x",
                     y = "y", region = NULL, costs = 'all'){
  # This bit is to silence the CRAN check warnings for literal column names
  from=to=..id=..x=..y=ID=Cell=V1=V2=TempID=dt=dW_l=dE_l=NULL
  #
  
  # First, "all" and "energetics" are just shorthand for a vector
  # of multiple cost types. Fix that.
  if (all(costs == 'all')){
    costs <- c("time","work","energy")
  } else if (all(costs == "energetics")){
    costs <- c("work","energy")
  }
  
  
  # If region is provided, filter out those values
  if (!is.null(region)){
    region <- regionMask(region, z_fix = z_fix)
    world <- world[(from %in% region) & (to %in% region),]
  }
  
  # Put everything in XY format, calculating centroids if polygons are provided
  if (sum(class(nodes) %in% 
          c("SpatialPolygons","SpatialPolygonsDataFrame")) > 1){
    nodes <- rgeos::gCentroid(nodes, byid = TRUE, id = id)
  }
  # Coerce to data.table
  nodes <- as.data.table(nodes)[,.(ID = get(..id), x = get(..x), y = get(..y))]
  nodes$Cell <- getCoords(nodes, z_fix = z_fix)
  
  # If no order is given, assume the node object is already in the desired order
  if (is.null(order)){
    order <- unlist(nodes[[id]])
  }
  
  # Get the IDs of nodes that will be needed for each leg
  starts <- order[1:(length(order)-1)]
  stops <-  order[2:length(order)]
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
    }
    if (cost == "work"){
      Graph <- igraph::graph_from_data_frame(world[,.(from,to,weight = dW_l)])
    }
    if (cost == "energy"){
      Graph <- igraph::graph_from_data_frame(world[,.(from,to,weight = dE_l)])
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
    # a numeric datatable containing XY coordinates for each traversed cell
    legList <- lapply(legList, igraph::as_ids)
    legList <- lapply(legList, tstrsplit, ",")
    legList <- rbindlist(legList,idcol="TempID")
    legList[,`:=`(x = as.numeric(V1), y = as.numeric(V2), V1 = NULL, V2 = NULL)]
    
    # Convert the XY list to SpatialPoints, then SpatialLines
    lineList <- list()
    for (i in unique(legList$TempID)){
      lineList <- append(lineList,
                         SpatialPoints(legList[TempID==i,.(x,y)],
                                       proj4string = crs(z_fix)))
      
    }
    lineList <- lapply(lineList,methods::as,"SpatialLines")
    lineList <- do.call(rbind,lineList)
    
    # Add the information about what segment each leg represents
    lineList <- SpatialLinesDataFrame(lineList,
                                      data = data.table(segment = nameList),
                                      match.ID = FALSE)
    # And place it in the originally-requested order
    lineList <- lineList[match(order, nameList),]
    
    # Save to the output pathList
    pathList[[cost]] <- lineList
    
  }
  if (length(pathList) > 0) {
    if (length(pathList) == 1 ){
      pathList <- pathList[[cost]]
    }
    return(pathList)
  }
}

#' Create a raster that can be used to define
#' the resolution, origin, and projection to be 
#' employed for all least-cost analyses. If a source
#' DEM has such properties you may use that.
#' 
#' @title Define the sampling grid
#' @param res A numeric of length one or two nrepresenting the spatial resolution
#' @param crs A \code{\link[raster]{crs}} object or character string containing
#' projection information.
#' @param dx The horizontal offset from the origin (see \code{\link[raster]{origin}}).
#' Default is 0.
#' @param dy The vertical offset from the origin (see \code{\link[raster]{origin}}).
#' Default is 0.
#' @return A RasterLayer object consisting of four cells, with resolution \code{res} and
#'  the origin at \code{x = nx} and \code{y = ny}.
#' @importFrom raster rasterFromXYZ
#' @importFrom raster origin<-
#' @importFrom data.table data.table
#' @examples 
#' projection <- "+proj=lcc +lat_1=48 +lat_2=33 +lon_0=-100 +datum=WGS84"
#' z_fix <- fix_z(res = 2, crs = projection)
#' @export
fix_z <- function(res, crs, dx = 0, dy = 0){
  if (length(res) == 1){
    res <- c(res,res)
  }
  z <- data.table(x = c(0,0,res[1],res[1]),
                  y = c(0,res[2],0,res[2]),
                  z = c(1,1,1,1))
  z <- rasterFromXYZ(z, res = res, crs=crs)
  origin(z) <- c(dx,dy)
  return(z)
}
