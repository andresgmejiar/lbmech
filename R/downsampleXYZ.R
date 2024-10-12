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
#' data <- downsampleXYZ(data = data, t_step = 3, z = 'z')
#' @export
downsampleXYZ <- function(data, t_step, t_cut = t_step * 10,
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
  
  # Add modified data to this blank list
  newData <- list()
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
      newData[[track]] <- points
    }, error = function(x) {
      1
    })
    i <- i + 1
    utils::setTxtProgressBar(pb,val = i)
  }
  utils::setTxtProgressBar(pb,val = max(whichlegs))
  newData <- rbindlist(newData)
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
