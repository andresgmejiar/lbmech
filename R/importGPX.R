#' Import \code{GPX} tracks from commercial GPS equipment as a data.table
#' ready for velocity function estimation. Note that the output coordinates
#' must still be converted to a conformal projection in meters if they are
#' to be used in functions other than \code{\link[lbmech]{getVelocity}},
#' \code{\link[lbmech]{dtVelocity}}, or \code{\link[lbmech]{downsampleXYZ}}.
#' 
#' @title Import GPX tracks
#' @param tracks A character string or vector with filepaths pointing to
#' the location of the GPX files
#' @param verbose Should a progress bar be printed? Default is \code{verbose = FALSE},
#' recommended particularly when used inside \code{\link[base]{lapply}}
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
#' @importFrom terra vect
#' @importFrom terra values
#' @importFrom terra geom
#' @examples 
#' # Get a list of GPX tracks in a directory
#' gpx <-  list.files(pattern = ".gpx$")
#' 
#' # Convert to data.table
#' gpx <- importGPX(gpx)
#' @export
importGPX <- function(tracks, verbose = FALSE){
  tracks_out <- data.table()
  j = 0
  if(verbose){
  pb <- utils::txtProgressBar(min=0,max=length(tracks)+0.00001,style=3)
  }
  for (i in tracks){
    try({
      gpx <- tmaptools::read_GPX(i,
                                 layers='track_points')$track_points
      gpx <- vect(gpx)
      
      gpx <- data.table(TrackID = stringr::str_extract(i,
                                                       pattern="[0-9A-Za-z_\\-]+(?=.gpx$)"),
                        PID = seq(1,nrow(gpx)),
                        t = values(gpx)$time,
                        long = geom(gpx)[,'x'],
                        lat = geom(gpx)[,'y'],
                        z = values(gpx)$ele)
      
      gpx <- gpx[order(t)]
      
      tracks_out <- rbind(tracks_out,gpx)
    },silent=TRUE)
    j <- j + 1
    if (verbose){
    utils::setTxtProgressBar(pb,j)
    }
  }
  return(tracks_out)
}