#' A function that for a given world of possible movement calculates
#' the transition cost for each in terms of a pre-defined time, work, and energy 
#' cost functions. \code{energyCosts} calls \code{timeCosts} if columns
#' named \code{'dt'} and \code{'dl_t'} are not present in the input data.table
#' 
#' @name energyCosts
#' @aliases timeCosts
#' @title Calculate time and energy costs 
#' @param DT A data.table containing at minimum columns 'dz' representing
#' the change in elevation and 'dl' representing planimetric distance
#' @param v_max The maximum velocity of the animal moving across the landscape,
#' in meters per second; see \code{\link[lbmech]{getVelocity}}.
#' @param k The topographic sensitivity factor; see \code{\link[lbmech]{getVelocity}}.
#' @param s The dimensionless slope of maximum velocity;
#' see \code{\link[lbmech]{getVelocity}}.
#' @param v_max The maximum velocity of the animal moving across the landscape,
#' in meters per second; see \code{\link[lbmech]{getVelocity}}.
#' @param method A character string for the method that energy costs per
#' unit stride should be calculated. One of \code{method \%in\% 
#' c('kuo','heglund','pontzer','oscillator')}; see references.
#' @param time The method by which time costs should be calculated by \code{energyCosts}
#' should \code{c('dt','dl_t')} not be column names in the input data.table. 
#' Default is \code{time = timeCosts}.
#' @param m The mass of the animal moving across the landscape, in kilograms.
#' @param epsilon The biomechanical efficiency factor for an animal moving across
#' the landscape. Default is \code{epsilon = 0.2}.
#' @param BMR The base metabolic rate of the object moving across the landscape
#' in Joules per second.
#' @param g The acceleration due to gravity, in meters per second per second.
#' Default is \code{g = 9.81} m/s^2, as for the surface of planet Earth.
#' @param l_s The average stride length, in meters. Required for
#' \code{method =  'kuo'}, \code{'pontzer'} or \code{'oscillator'}, 
#' ignored for \code{'heglund'}
#' @param L The average leg length. Required for \code{method =  'kuo'},
#' ignored for \code{'heglund'}, \code{'pontzer'} and \code{'oscillator'}.
#' @param gamma The fractional maximal deviation from average velocity per stride.
#' Required for \code{method = 'oscillator'}, ignored otherwise
#' @param water Logical. If \code{FALSE} (the default), movement costs are calculated
#' as if over land. If \code{water = TRUE}, movement costs are calculated considering
#' moving water. 
#' @param row_speed How fast can a person move over water? Default is \code{row_speed = NULL},
#' but required if \code{water} are provided. 
#' @param row_work How much work in joules per second does a person use to move over water?
#' Default is \code{row_work = NULL}, but required if \code{water} is provided. 
#' @param ... Additional parameters to pass to \code{timeCosts}
#' @references {
#'  Heglund, N. C., Cavagna, G. A., and Taylor, C. R. (1982). "Energetics and 
#'  mechanics of terrestrial locomotion. III. Energy changes of the centre of
#'  mass as a function of speed and body size in birds and mammals." 
#'  \emph{Journal of Experimental Biology} 97(1):41-56. 
#'  \url{https://doi.org/10.1242/jeb.97.1.41}.
#'  
#'  Kuo, Arthur D. (2007). "The six determinants of gait and the inverted 
#'  pendulum analogy: A dynamic walking perspective." \emph{Human Movement Science}
#'  26(4):617-656. 
#'  \url{https://doi.org/10.1016/j.humov.2007.04.003}.
#'  
#'  Pontzer, Herman (2016). "A unified theory for the energy cost of legged locomotion"
#'  \emph{Biology Letters} 12(2):20150935.
#'  \url{https://doi.org/10.1098/rsbl.2015.0935}
#'  }
#' @return For \code{timeCosts}, A data.table object with two columns:
#'
#' (1) \code{$dl_t} The predicted walking speed in meters per second
#' when walking between the \code{from} and \code{to} cells
#'
#' (2) \code{$dt} The predicted amount of time spent walking between
#' the \code{from} and \code{to} cells
#' 
#' For \code{energyCosts}, a data.table object with five columns:
#' 
#' (1) \code{$dt} The predicted amount of time spent walking between
#' the \code{from} and \code{to} cells
#'
#' (2) \code{$dU_l} The predicted work against gravitational potential energy
#' in Joules when walking between the \code{from} and \code{to} cells
#'
#' (3) \code{$dK_l} The predicted kinematic work in Joules when walking
#' between the \code{from} and \code{to} cells
#'
#' (4) \code{$dW_l} The total predicted energy lost due to biomechanical
#' work when walking between the \code{from} and \code{to} cells.
#'
#' (5) \code{$dE_l} The net metabolic expenditure exerted when walking
#' between the \code{from} and \code{to} cells.
#' 
#' @importFrom data.table :=
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
#' grid <- makeGrid(dem = dem, nx = n, ny = n, sources = TRUE)
#' 
#' # Select all tiles that exist between x = (12000,16000) and y = (32000,36000)
#' tiles <- ext(c(12000,16000,32000,36000))
#' tiles <- as.polygons(tiles)
#' crs(tiles) <- crs(grid)
#' tiles <- whichTiles(region = tiles, polys = grid)
#' 
#' # Make a world but limit it to the DEM grid size
#' defineWorld(source = grid, cut_slope = 0.5, 
#'             res = res(dem), dir = dir, overwrite=TRUE)
#' 
#' # Calculate the energetic and temporal costs
#' calculateCosts(costFUN = energyCosts, 
#' tiles = tiles, dir = dir,
#' m = 70, v_max = 1.5, BMR = 76, k = 3.5, s = 0.05, l_s = 1,
#' L = 0.8)
#' @export
timeCosts <- function(DT, v_max, k, s, row_speed = NULL, water = FALSE){
  # This bit to silence CRAN warnings
  DT <- copy(DT)
  dl_t=..v_max=..k=dz=dl=..s=dt=..WaterSpeed=..row_speed=WaterSpeed=NULL
  if (water == FALSE){
    DT[, dl_t := ..v_max * exp(-(..k) * abs(dz/dl - ..s))
    ][, dt := dl_t ^ -1 * dl] #* ..l_p / ..l_s
    return(DT[,.(dl_t, dt)])
  } else {
    DT[, dl_t := WaterSpeed + ..row_speed
    ][, dt := dl * dl_t ^ -1]
    DT <- DT[dl > 0]
    DT[dl_t <= 0, `:=`(dt = Inf)]
    return(DT[,.(dl_t,dt)])
  }
}

#' @rdname energyCosts
#' @export
energyCosts <- function(DT, method = 'kuo', m = NULL, BMR = NULL, g = 9.81, 
                        epsilon = 0.2, l_s = NULL, L = NULL, gamma = NULL,
                        time = timeCosts, water = FALSE, row_work = NULL, ...){
  DT = copy(DT)
  #This bit to silence the CRAN warnings
  dU_l=..m=..g=dz=dK_l=dl_t=..l_s=dl=..L=..gamma=dW_l=..epsilon=dE_l=..BMR=NULL
  dt=..row_work=NULL
  
  if ((sum(c("dl_t","dt") %in% names(DT)) != 2)){
    DT <- cbind(DT, timeCosts(DT, water = water, ...))
  }
  
  if (water == FALSE){
    DT[, dU_l := ..m * ..g * dz
    ][dU_l < 0, dU_l := 0]
    
    ## (2) Calculate the work based on a user-selected function
    if (method == 'kuo'){
      # Kuo's function for human movement
      DT[, dK_l := 1 / 4 * ..m * dl_t^2 * ..l_s * dl / ..L^2]
    } else if (method == 'heglund'){
      # Heglund et al.'s function for arbitrary quadripeds
      DT[, dK_l := (0.478 * dl_t ^ 1.53 + 0.685 * dl_t + 0.072) * dl_t ^ -1 * dl * ..m]
    } else if (method == 'oscillator'){
      DT[, dK_l := 2 * ..m * dl_t ^2 * ..gamma * dl / ..l_s]
    } else if (method == 'pontzer'){
      # Pontzer's cross-bridge method 
      DT[, dK_l := (8 * ..m^(-0.34) + 
           100 * (1 + sin(2 * atan(dz/dl) - 74 * 3.14159/180)) * ..m^(-0.12)) * 
           dl / ..l_s]
    }
    
    ## (3) Finally, calculate the total work and energy
    DT[, dW_l := (dU_l + dK_l) / ..epsilon
    ][, dE_l := dW_l + ..BMR * dt][]
    
    return(DT[,.(dt,dU_l, dK_l, dW_l, dE_l)])
  } else {
    DT[, dW_l := ..row_work * dt
    ][, dE_l := dW_l/..epsilon + ..BMR * dt]
    
    DT[dl_t <= 0, `:=`(dt = Inf, dU_l = 0, dK_l = Inf,
                       dW_l = Inf, dE_l = Inf)]
    return(DT[,.(dt, dU_l = 0, dK_l = dW_l, dW_l, dE_l)])
  }
  
}
