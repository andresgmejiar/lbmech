#' Calculate the velocity function for an animal from \emph{(x,y,z,t)} data
#' such as from GPS collars, assuming a function of form Tobler (see ####).
#' 
#' \code{dtVelocity} is a wrapper for \code{\link[lbmech]{getVelocity}} for use 
#' inside of a data.table with observations from multiple sources
#' (be it different individual animals and/or different instances
#' from the same animal). In a \code{\link[data.table]{data.table}},
#' use this function in the \code{j} slot, passing it along \code{.SD}.
#'
#' @name getVelocity
#' @aliases dtVelocity
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
#' @param s_init A number representing the value for dimensionless slope of
#' maximum velocity at which to initiate the nonlinear regression. Default is
#' \code{s_init = -0.05}.
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
#' @param ... Additional parameters to pass to \code{\link[lbmech]{importMap}}.
#' Ignored unless \code{z} is a polygon pointing to source dem files.
#' @return A list, containing the following entries:
#'
#' (1) \code{$model}, containing an object of class \code{\link[quantreg]{nlrq}} containing the output
#' model from the nonlinear quantitative regression.
#'
#' (2) \code{$vmax}, containing the identified maximum velocity.
#'
#' (3) \code{$s}, containing the identified angle of maximum velocity.
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
#' suppressWarnings( data[, z := NULL])
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
#' grid <- makeGrid(dem = dem, nx = n, ny = n, sources = TRUE)
#' 
#' velocity <- getVelocity(data = data, z = grid, dir = dir, res = res(dem))
#' 
#' #' # Example 4:
#' # If you want to get the values by group within a data.table
#' 
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
getVelocity <- function(data, x = 'x', y ='y', degs = FALSE, dl = NULL, z = 'z', 
                        dt = 'dt', ID = 'ID', tau = NULL, tau_vmax = 0.995,
                        tau_nlrq = 0.95, k_init = 3.5, s_init = -0.05,
                        v_lim = Inf, slope_lim = 1,
                        tile_id = "TILEID", vals = "location",
                        dir = tempdir(), ...){
  # This bit is to silence the CRAN check warnings for literal column names
  ..x=..y=..z=..ID=..dt=..z_vals=Order=dx=dy=dz=dl_dt=..dl=dz_dl=NULL
  #
  
  # Standardize filepath names
  dir <- normalizePath(dir)
  rd <- dir
  #rd <- normalizePath(paste0(dir),mustWork = FALSE)
  if (!dir.exists(rd)){
    dir.create(rd)
  }
  
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
    
    # Create an empty data.table in which to store extracted values. There's a
    # small chance some points will land in two sectors, so we'll have
    # to collect all and consolidate
    data_new <- data.table()
    for (tile in tiles){
      elevs <- importMap(tile, polys = z, tile_id = tile_id, vals = vals,
                         dir = normalizePath(dir), ...)
      
      d <- data_points[get(tile_id) == tile]
      z_vals <- lapply(cellFromXY(elevs,d[,.(get(..x),get(..y))]),
                       function(x) elevs[x])
      z_vals <- lapply(z_vals, function(x) if (length(unlist(x)) == 0) NA else x)
      d[,z := unlist(..z_vals)]
      
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
  velocity <- quantreg::nlrq(dl_dt ~ v_max * exp(-k * abs(dz_dl - s)),
                             data = data[(dl_dt <= v_lim) & 
                                           abs(dz_dl) <= slope_lim], 
                             tau = tau_nlrq, 
                             start=list(k=k_init,s=s_init))
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
    s = as.numeric(stats::coefficients(velocity)[['s']]),
    k = as.numeric(stats::coefficients(velocity)[['k']]),
    tau_vmax = tau_vmax,
    tau_nlrq = tau_nlrq,
    data = data
  )
  # Return the environment to what it was
  on.exit(options(warn = defaultW))
  return(out)
  
}

#' @rdname getVelocity
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
      j[['s_p']] <- as.numeric(n[,"Pr(>|t|)"][2])
      return(j)
    }, 
    error = function(e) {
      j[['k']] <- as.numeric(NA)
      j[['s']] <- as.numeric(NA)
      j[['v_max']] <- as.numeric(NA)
      j[['k_p']] <- as.numeric(NA)
      j[['s_p']] <- as.numeric(NA)
      return(j)
    })
}