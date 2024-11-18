#' Import \code{GPX} tracks from commercial GPS equipment as a data.table
#' ready for velocity function estimation. Note that the output coordinates
#' must still be converted to a conformal projection in meters if they are
#' to be used in functions other than \code{\link[lbmech]{getVelocity}},
#' \code{\link[lbmech]{dtVelocity}}, or \code{\link[lbmech]{downsampleXYZ}}.
#' 
#' @title Import GPX tracks
#' @param track A character string the filepaths pointing to
#' the location of the GPX file
#' @return A data.table with five or six columns:
#' 
#' (1) \code{$TrackID} The name of the input GPX track
#' 
#' (2) \code{$PID} The order that point appeared in the GPX file
#' 
#' (3) \code{$t} The time stamp, with a starting date of January 1st, 1900 if no
#' day is listed
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
#' ## Convert the first file to data.table (not run)
#' # gpx <- importGPX(gpx[1])
#' # Convert all the files to data.tables
#' gpx <- lapply(gpx, importGPX)
#' @export
importGPX <- function(track){
  tryCatch({
    # Silence CRAN warnings
    hour=time=hdiff=PID=long=lat=z=NULL
    
    suppressWarnings(gpx <- readLines(track))

    gpx <- paste(gpx,collapse = '')
    gpx <- unlist(stringr::str_extract_all(gpx,'(?<=<trk>).*(?!</trk>).*(?=</trk>)'))
    if (length(gpx) == 0) {
      return(data.table())
    }
    gpx <- unlist(stringr::str_split(gpx,'(?<=<trkpt)'))[-1]
  
    # Extract the appropriate strings; we suppress warnings in case the file lacks
    # elevation or date/time information
    suppressWarnings({
    vals <- data.table(lat = stringr::str_extract(gpx, "(?<= lat=\")[\\d\\.]+"),
                       long = stringr::str_extract(gpx, "(?<= lon=\")[\\d\\.]+"),
                       z = stringr::str_extract(gpx,'(?<=<ele>)[\\d\\.]+(?=</ele>)'),
                       time = stringr::str_extract(gpx,'(?<=T)[\\d\\.\\:]+(?=Z)'),
                       date = stringr::str_extract(gpx,'(?<=<time>)[\\d\\-]*(?=T)'))
    })
    # Convert to numeric
    
    vals[, (names(vals)[1:3]) := lapply(.SD, as.double), .SDcols = names(vals)[1:3]]
 
    # If a date isn't present, we need to assign a dummy date way in the past,
    # and add a day everytime we go past midnight. Meanwhile we convert to datetime
    if (!all(is.na(unique(vals$date)))){
    if (all(unique(vals$date) == '')){
      # Get the hour number
      vals[, hour := as.numeric(stringr::str_sub(time, 1,2))]
      
      # If it's over midnight, we mark a new day. We sum the total number of marks
      # to get the number of days we need to add to a given dummy value
      vals[, hdiff := fifelse(abs(hour - data.table::shift(hour)) >= 23, 1, 0) 
      ][1, hdiff := 0][, hdiff := cumsum(hdiff)]
      
      # Assign dummy start date
      date = '1900-01-01'
      
      # Convert to datetime object, adding a day's worth of seconds for every
      # new day mark
      vals[, t := as.POSIXct(paste(date, time),tz = 'GMT') + hdiff * 3600*24]
    } else {
      vals[!is.na(date), t := as.POSIXct(paste(date, time),tz = 'GMT')]
    }
    } else{
      vals[, t := NA]
    }
    
    vals[, PID := 1:nrow(vals)]
    
    vals <- vals[order(t, na.last = FALSE)]
    
    if ('z' %in% names(vals)){
      return(vals[,.(PID,t,long,lat,z)])
    } else {
      return(vals[,.(PID,t,long,lat)])
    }
  }, error = function(e) data.table() ,silent=TRUE)
}
