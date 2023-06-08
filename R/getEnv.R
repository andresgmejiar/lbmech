#' Get environmental variables from trusted, global datasets,
#' and optionally extract it for point locations
#' 
#' @title Download and extract environmental data for known locations.
#' @param region A SpatRaster, Raster*, SpatVector, Spatial* or character
#' string representing the area of interest. Required unless \code{locs} are
#' provided. If \code{region = NULL} but \code{locs} are provided, then 
#' \code{region = ext(locs)}.
#' @param locs (Optional) A data.frame or something coercible to a data.table containing
#' all observations to which data should be extracted. 
#' @param x A character string representing the data column containing the 'x'
#' coordinates in meters or degrees. Ignored if data is of class SpatVector
#' or Spatial.
#' @param y A character string representing the data column containing the 'y'
#' coordinates in meters or degrees. Ignored if data is of class SpatVector
#' or Spatial, and ignored for distance calculations if \code{dl} is
#' provided.
#' @param dem (Optional) A SpatRaster or Raster object to use for slope and
#' topographic wetness calculations. Default is \code{dem = NULL}, which
#' will download the appropriate SRTM data using 
#' \code{\link[elevatr]{get_elev_raster}(z = zoom)}
#' @param z_fix A raster with the origin, projection, and resolution of the 
#' desired output rasters. 
#' @param proj A crs object or character string representing the projection
#' information for x,y coordinates. If z_fix is provided, default is for 
#' \code{proj = z_fix}. Ignored if locs if \code{locs} is 
#' SpatVector or Spatial. If \code{locs} is a data.table or data.frame and 
#' \code{proj = NULL}, the function will attempt \code{proj = crs(dem)} as
#' a last resort.
#' @param zoom  Considered only if \code{dem = NULL}  
#' The zoom level to be downloaded. See documentation for the \code{z} parameter
#' in \code{\link[elevatr]{get_elev_raster}} for further information.
#' @param z_min The minimum allowable elevation. Useful if DEM source includes
#' ocean bathymetry as does the SRTM data from AWS. Default is \code{z_min = 0},
#' but set to \code{NULL} to disable.
#' @param filt Numeric. Size of moving window to apply a low-pass filter to
#' terrain-based metrics (\code{slope} VIA DEM, \code{twi} directly)
#' @param slope Logical. Should the slope at each input location be calculated?
#' Default uses SRTM data from \code{\link[elevatr]{get_elev_raster}}, but
#' this can be overwritten by providing a \code{dem}
#' @param twi Logical. Should the topgraphic wetness index be calculated?
#' Default uses the provided \code{dem} as a slope source, and 
#' \href{https://www.hydrosheds.org/hydrosheds-core-downloads#w-tabs-0-data-w-pane-3}{HydroSHEDS 15s arcsecond}
#' data as a flow accumulation source, however this can
#' be altered by modifying the \code{acc} parameter. All data are corrected
#' for geodesic distorsion, which can be time-consuming if \code{acc = 3}
#' @param depth Logical. Should the depth to bedrock be obtained? If TRUE,
#' SoilGrids v. 1.0 data are downloaded from 
#' \href{https://files.isric.org/soilgrids/former/2017-03-10/data/}{ISRIC SoilGrids}. 
#' Unfortunately, these data are only provided globally and therefore the 
#' download is 1.34 GB. 
#' @param soils A character vector containing zero, all or, some of 
#' \code{c('sand','silt','clay','cec','soc')} representing soil's parts per million
#' sand, silt, and clay content, the cation exchange capacity, and the organic carbon
#' (respectively). Data are downloaded from \href{https://maps.isric.org/}{ISRIC SoilGrids}
#' using their OGC Web Service API allowing the download to be fairly efficient.
#' @param climate A character vector containing zero, all or, some of 
#' c("tmin","tavg","tmax","prec","srad") representing the minimum average temperature,
#' average temperature, maximum average temperature, precipitation, and solar radiation.
#' Data are downloaded from \href{http://www.worldclim.com/version2}{WorldClim v. 2}
#' using \code{\link[geodata]{worldclim_country}}.
#' @param acc One of (1) A Raster or SpatRaster representing the flow accumulation
#' in units of incoming per cell, or (2) an integer equal to one of 
#' \code{acc = c(3, 15, 30)}, representing the resolution in arcseconds to download
#' from \href{https://www.hydrosheds.org/hydrosheds-core-downloads#w-tabs-0-data-w-pane-3}{HydroSHEDS}.
#' Default is \code{acc = 15}. Note that the 3s data is 2.29 GB.
#' @param vars A SpatRaster object (or something coercible to it using the
#' \code{\link[terra]{rast}} function) of one or multiple layers containing
#' custom environmental data. Data with names equivalent to a parameter that has
#' been activated will be ignored. 
#' @param overwrite Should the 
#' @param dir A filepath to the directory being used as the workspace.
#' Default is \code{tempdir()} but unless the analyses will only be performed a few
#' times it is highly recommended to define a permanent workspace so that 
#' files only need to be downloaded once.
#' @param ... Additional parameters to pass to \code{\link[lbmech]{fix_z}}.
#' @importFrom terra extract
#' @importFrom terra rast
#' @importFrom terra vect
#' @importFrom terra project
#' @importFrom terra values
#' @importFrom terra geom
#' @importFrom terra ext
#' @importFrom terra as.polygons
#' @importFrom terra writeVector
#' @importFrom terra resample
#' @importFrom terra res
#' @importFrom terra writeRaster
#' @importFrom terra focal
#' @importFrom terra terrain
#' @importFrom terra extract
#' @importFrom terra crop
#' @importFrom data.table as.data.table
#' @return A folder in the `dir` with the requested datasets.
#' @examples 
#' 
#' # Get coordinates for Mexico City in longlat
#' data <- data.table(x = -99.1332,
#'                    y = 19.4326)
#'
#' ## Not run due to size of downloads
#' # getEnv(data, filt = 3)
#' @export
getEnv <- function(region = NULL, locs = NULL,
                   x = 'x', y = 'y', dem = NULL, proj = NULL,
                   z_fix = NULL, zoom = 10, z_min = 0, filt = 0, 
                   slope = FALSE, twi = FALSE, depth = FALSE, 
                   soils = NULL,
                   climate = NULL,
                   acc = 15, 
                   vars = NULL, 
                   overwrite = FALSE,
                   dir = tempdir(),...){
  # This bit to silence CRAN warnings
  a=lx=ly=acc_deg=..twi=..depth=..s=NULL
  
  dir <- normalizePath(paste0(dir,'/Environment/'),mustWork=FALSE)
  subdirs <- c("/Raw","/Local","/Diff")
  subdirs <-  normalizePath(paste0(dir,subdirs),mustWork=FALSE)
  
  if (overwrite == TRUE){
    if (dir.exists(dir)){
      unlink(dir, recursive = TRUE)
    }
  }
  if (!dir.exists(dir)){
    dir.create(dir)
  }
  
  for (i in subdirs){
    if (!dir.exists(i)){
      dir.create(i)
    }
  }
  
  rd <- normalizePath(dir)
  outNames <- c()
  
  if (methods::is(region,'Raster')){
    region <- rast(region)
  } else if (methods::is(region,'Spatial')){
    region <- vect(region)
  }
  
  if (!is.null(z_fix)){
    proj <- crs(z_fix)
    region <- project(region,z_fix)
  }
  
  if (is.null(proj)){
    proj <- crs(region)
  }
  
  if (!is.null(locs)){
    if (methods::is(locs,"data.frame")){
      locs <- vect(locs, geom = c(x,y), crs=proj)
      locs <- project(locs,proj)
      locs <- cbind(as.data.table(
        values(locs))[,.SD,.SDcols =
                        names(locs)[!(names(locs) %in% c('x','y'))]],
        geom(locs)[,c('x','y')])
    } else if (methods::is(locs,"Spatial")){
      locs <- vect(locs)
    }
    if (methods::is(locs,"SpatVector")){
      locs <- project(locs,proj)
      locs[,c('x','y')] <- geom(locs)[,c('x','y')]
      locs <- as.data.table(locs)
    } 
  }
  
  extentPath <- normalizePath(paste0(dir,'/extent.gpkg'),mustWork = FALSE)
  if (file.exists(extentPath)){
    if (!is.null(z_fix) | (length(list(...)) != 0)){
      warning("Region provided, but extent.gpkg already exists in the /Environment/ dir. Ignoring input. Use overwrite = TRUE to clear directory\n")
    }
    extent <- vect(extentPath)
  } else if (is.null(region)){
    # Create a bounding polygon and raster for input locations
    if (!is.null(locs)) {
      if(methods::is(locs,'data.table')){
        extent <- ext(vect(locs,geom=c(x,y),crs=proj))
        extent <- as.polygons(extent,crs=proj)
        writeVector(extent, extentPath)
      }
      else {
        extent <- ext(locs)
        extent <- as.polygons(extent,crs=proj)
        writeVector(extent, extentPath)
      }
    }
  } else {
    extent <- ext(region)
    extent <- as.polygons(extent,crs=proj)
    writeVector(extent, extentPath)
  } 
  
  extent <- vect(extentPath)
  
  zfixPath <- normalizePath(paste0(dir,'/z_fix.tif'),mustWork = FALSE)
  if (file.exists(zfixPath)){
    if (!is.null(z_fix) | (length(list(...)) != 0)){
      warning("z_fix provided, but z_fix.tif already exists in the /Environment/ dir. Ignoring input. Use overwrite = TRUE to clear directory\n")
    }
  } else if (is.null(z_fix)){
    # If no z_fix is provided, make one
    z_fix <- fix_z(proj = proj, ...)
    z_fix <-  resample(z_fix, rast(
      xmin=ext(extent)[1]%/% res(z_fix)[1] * res(z_fix)[1],
      xmax=ext(extent)[2] %/% res(z_fix)[1] * res(z_fix)[1] +res(z_fix)[1],
      ymin=ext(extent)[3]%/% res(z_fix)[2] * res(z_fix)[2] ,
      ymax=ext(extent)[4] %/% res(z_fix)[2] * res(z_fix)[2] +res(z_fix)[2],
      res = res(z_fix), crs=proj,
    ))
    z_fix[cells(z_fix) & !cells(z_fix)] <- 1
    writeRaster(z_fix, zfixPath)
  } else if (methods::is(z_fix, "Raster")){
    z_fix <- rast(z_fix)
    writeRaster(z_fix, zfixPath)
  } else if (methods::is(z_fix,"SpatRaster")){
    writeRaster(z_fix, zfixPath)
  }
  z_fix <- rast(zfixPath)
  
  res <- res(z_fix)
  
  # If DEM is missing, download elevations. If raster, convert to slope
  slopePath <- normalizePath(paste0(subdirs[2],"/slope.tif"),mustWork = FALSE)
  do_slope <- !file.exists(slopePath)
  twiPath <- normalizePath(paste0(subdirs[2],"/twi.tif"),mustWork = FALSE)
  do_twi <- !file.exists(twiPath)
  if (slope | twi) {
    if (is.null(dem)){
      demPath <- normalizePath(paste0(subdirs[1],"/dem.tif"),mustWork = FALSE)
      if (!file.exists(demPath) & (do_twi | do_slope)){
        dem <- rast(elevatr::get_elev_raster(methods::as(extent,'Spatial'),
                                             z=zoom,src='aws',prj=proj))
        # Apply elevation pre-processing
        if (!is.null(z_min)){
          dem[dem < 0] <- NA
        } 
        dem <- resample(dem,z_fix)
        if (filt != 0){
          dem <- focal(dem,w=filt,fun=mean,na.policy='omit')
        }
        writeRaster(dem,demPath)
      }
      dem <- rast(demPath)
    } else if (methods::is(dem,'Raster')){
      dem <- rast(dem)
    } 
    if (methods::is(dem,'SpatRaster') & !file.exists(demPath)){
      writeRaster(dem, demPath)
    }
    dem <- rast(demPath)
    # Calculate slope
    do_slope <- slope
    slope <- tan(terrain(dem,v='slope',unit='radians'))
    if (do_slope & !is.null(locs)){
      locs[, slope := extract(slope,data.table(x = x, y = y))$slope]
      outNames <- c(outNames,'slope')
    }
    if (do_slope & !file.exists(slopePath)){
      writeRaster(slope,slopePath)
    }
    
    if (twi){
      # If no topographic wetness index raster, make it. 
      if (!file.exists(twiPath)){
        # Get the flow accumulation if not provided. This is an even bigger download
        # from a legacy dataset
        accPath <- normalizePath(paste0(subdirs[1],"/acc.tif"),mustWork=FALSE)
        if (!file.exists(accPath) & methods::is(acc,'numeric')){
          utils::download.file(paste0("https://data.hydrosheds.org/file/hydrosheds-v1-acc/na_acc_",acc,"s.zip"),
                        normalizePath(paste0(rd,'/Temp.zip'),mustWork = FALSE))
          utils::unzip(normalizePath(paste0(rd,"/Temp.zip"), mustWork = FALSE),
                       exdir = normalizePath(paste0(rd,"/Temp"),
                                             mustWork = FALSE),
                       junkpaths = FALSE)        
          files <- unlist(list.files(paste0(rd,"/Temp"), pattern = ".tif$",
                                     recursive = TRUE,
                                     full.names = TRUE))
          files <- rast(files)
          unlink(normalizePath(paste0(rd,"/Temp.zip")),recursive = TRUE)
          unlink(normalizePath(paste0(rd,"/Temp")),recursive = TRUE)
          writeRaster(files,accPath)
        }
        # Import flow accumupation raster, get resolution and crop to aoi
        acc <- rast(accPath)
        accRes <- res(acc)
        acc <- crop(acc,project(extent,crs(acc)))
        names(acc) <- "acc_deg"
        acc <- project(acc,'+proj=longlat')
        
        # Get the area of each pixel in meters, since they're provided in units of 
        # accumulated cells, and the cells would have a resolution of 3 arcseconds = 
        # 9 square arcseconds. Calculate specific flow area for TWI
        accXY <- rastToTable(acc)
        accXY$lx <- geosphere::distGeo(accXY[, .(x = x - accRes[1]/2,y)], 
                                       accXY[, .(x = x + accRes[1]/2,y)])  
        accXY$ly <- geosphere::distGeo(accXY[, .(x,y - accRes[2]/2)], 
                                       accXY[, .(x,y + accRes[2]/2)])  
        accXY[, a := lx * ly * acc_deg / (lx/4 + ly/4 + sqrt(lx^2 + ly^2)/2)]
        acc$a <- acc$acc_deg
        acc$a[cells(acc$acc_deg)] <- accXY$a
        rm(accXY)
        
        # Calculate TWI
        acc <- log(acc$a/project(slope,acc))
        names(acc) <- 'twi'
        
        acc[is.infinite(acc)] <- max(unlist(acc[!is.infinite(acc)]),
                                     na.rm=TRUE)
        acc[!cells(acc)] <- as.numeric(global(acc,max,na.rm=TRUE))
        
        if (filt != 0){
          acc <- focal(acc,w=filt,fun=mean,na.rm=TRUE)
        }
        names(acc) <- 'twi'
        acc <- project(acc,z_fix)
        writeRaster(acc,twiPath)
        rm(acc)
      }
      
      # Import Topographic Wetness Index, project and extract values
      twi <- rast(twiPath)
      names(twi) <- 'twi'
      
      if (!is.null(locs)){
        locs[, twi := extract(..twi,data.table(x=x,y=y))$twi]
      }
      rm(twi)
      outNames <- c(outNames,'twi')
      
    }
    rm(slope)
  }
  
  
  if (depth){
    # Get the depth to bedrock; sadly this is a very large download from a legacy
    # dataset. 
    depth_path_down <- normalizePath(paste0(subdirs[1],"/depth.tif"),mustWork=FALSE)
    depth_path <- normalizePath(paste0(subdirs[2],"/depth.tif"),mustWork=FALSE)
    
    if (!file.exists(depth_path)){
      if (!file.exists(depth_path_down)){
        utils::download.file("https://files.isric.org/soilgrids/former/2017-03-10/data/BDRICM_M_250m_ll.tif",
                      depth_path_down)
        utils::download.file("https://files.isric.org/soilgrids/former/2017-03-10/data/BDRICM_M_250m_ll.tif.xml",
                      normalizePath(paste0(depth_path_down,".xml")))
      }
      depth <- rast(depth_path_down)
      depth <- crop(depth,project(extent,crs(depth)))
      depth <- project(depth,z_fix)
      writeRaster(depth,depth_path)
    }
    depth <- rast(depth_path)
    
    if (!is.null(locs)){
      locs[, depth := extract(..depth,data.table(x=x,y=y))$depth]
    }
    
    rm(depth)
    outNames <- c(outNames,'depth')
  }
  
  soilExt <- project(z_fix,'epsg:4326')
  # Download soil coverage rasters
  for (s in soils){
    soil <- paste0("https://maps.isric.org/mapserv?map=/map/",s,".map&",
                   "SERVICE=WCS&",
                   "VERSION=2.0.1&",
                   "REQUEST=GetCoverage&",
                   "COVERAGEID=",s,"_5-15cm_mean&",
                   "FORMAT=GEOTIFF_INT16&",
                   "SUBSET=X(",paste(ext(soilExt)[1:2],collapse=','),")&",
                   "SUBSET=Y(",paste(ext(soilExt)[3:4],collapse=','),")&",
                   "SUBSETTINGCRS=http://www.opengis.net/def/crs/EPSG/0/4326&",
                   "OUTPUTCRS=http://www.opengis.net/def/crs/EPSG/0/4326")
    
    soilPath_down <- normalizePath(paste0(subdirs[1],"/",s,".tif"),
                                   mustWork = FALSE)
    soilPath <-  normalizePath(paste0(subdirs[2],"/",s,".tif"),
                               mustWork = FALSE)
    if (!file.exists(soilPath)){
      # If soil raster is missing, download it
      if (!file.exists(soilPath_down)){
        utils::download.file(soil,soilPath_down, mode = 'wb')
      }
      soil <- rast(soilPath_down)
      soil <- crop(soil,project(extent,crs(soil)))
      soil <- project(soil,z_fix)
      writeRaster(soil, soilPath)
    } 
    soil <- rast(soilPath)
    if (!is.null(locs)){
      locs[, (s) := extract(soil,data.table(x = x, y = y))[[..s]]]
    }
    outNames <- c(outNames,s)
    rm(soil)
  }
  
  # Import the monthly climate data
  for (var in climate){
    climPath <- normalizePath(paste0(subdirs[2],"/",var,".tif"),
                              mustWork = FALSE)
    if (!file.exists(climPath)){
      countries <- crop(
        buffer(
          project(
            vect(rworldmap::getMap(resolution='low')),proj),5000),z_fix)
      countries <- as.character(unlist(countries[['ISO_A3']]))
      countries <- lapply(countries,geodata::worldclim_country,
                          var = var, path = subdirs[1])
      clim <- do.call(merge,countries)
      names(clim) <- paste0(var,"_",1:12)
      
      clim <- crop(clim,project(extent,crs(clim)))
      clim <- terra::project(clim,z_fix)
      writeRaster(clim,climPath)
    }
    clim <- rast(climPath)
    if (!is.null(locs)){
      locs[, (names(clim)) :=
             extract(clim,data.table(x=x,y=y))[names(clim)]]
    }
    outNames <- c(outNames,names(clim))
    rm(clim)
  }
  
  
  if (!is.null(vars)){
    
    if (!methods::is(vars,'SpatRaster')){
      vars <- rast(vars)
    }
    customNames <- names(vars)
    customPath_down <- normalizePath(paste0(subdirs[1],"/",deparse(substitute(vars)),".tif"),
                                     mustWork = FALSE)
    customPath <- normalizePath(paste0(subdirs[2],"/",deparse(substitute(vars)),".tif"),
                                mustWork = FALSE)
    
    if (!is.null(locs)){
      locs[, (names(vars)) :=
             extract(vars,data.table(x=x,y=y))[names(vars)]]
    }
    outNames <- c(outNames,customNames)
  }
  if (!is.null(locs)){
    return(locs[,.SD,.SDcols = outNames])
  }
}  
  