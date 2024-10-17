#' Get the circumference of the earth given a known point and initial bearing
#' from that point
#' 
#' @title Great Circle Circumference
#' @param lon Numeric. A known longitude along a Great Circle.
#' @param lat Numeric. A known latitude along a Great circle
#' @param b Numeric. The initial bearing defining the Great Circle
#' @param n How many points along the great circle should be modeled? More points
#' results in higher precision. Default is \code{n = 72}
#' @param a Equatorial radius. Default is for WGS84
#' @param f Ellipsoidal flattening. Default is for WGS84
#' @importFrom data.table data.table
#' @examples 
#' 
#' circum <- GreatCircleCircum(-98.844, 19.693, , b = 175.8)
#' @export
GreatCircleCircum <- function(lon,lat,b,n = 72, 
                              a = 6378137, f = 1/298.257223563){
  # To silence CRAN warnings
  first=NULL
  # The Geosphere package returns NAs if the bearing is zero, since it relies
  # on interpolating the coordinates along longitudes. Filter out vertical cases
  if (!(b %% 360 == 0)){
    # Get a circumference from a known great circle
    g <- geosphere::greatCircleBearing(data.table(lon,lat),
                                       b = b,
                                       n = n)
    # Circumference the summed distance of all thse points (plus a repeat of the)
    # first to close the circle
    return(sum(geosphere::distGeo(rbind(g,first(g)),
                                  a = a,
                                  f = f), na.rm=TRUE))
  } else {
    # In the vertical case, it's just four times the distance from the equator
    # to a pole. 
    return(geosphere::distGeo(data.table(0,0),
                   data.table(0,90), a = a, f = f) * 4)
  } 
}