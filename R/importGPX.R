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