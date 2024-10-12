#' A function that for a given world of possible movement calculates
#' the transition cost for each in terms of a user-defined cost function
#'
#' @title Calculate movement costs according to a cost function
#' @param tiles A character vector--such as the output to
#' \code{\link[lbmech]{whichTiles}}---containing the unique tile IDs for sectors that
#' should be in the workspace. Default is \code{NULL}.
#' @param costFUN A cost function such as (\code{\link[lbmech]{timeCosts}} or
#' \code{\link[lbmech]{energyCosts}}). The input to such a function should be 
#' a \code{data.table} object with column names present in the makeWorld file
#' (initially \code{c('x_i','x_f','y_i','y_f','z_i','z_f','dz','dl','dr')}), and
#' the output should be an data.table object with the same number of rows and the
#' desired output variables as the only column names. Constants can be passed in the 
#' \code{...} slot. Default is 
#' \code{costFUN = \link[lbmech]{energyCosts}}
#' @param dir A filepath to the directory being used as the workspace, the same
#' one instantiated with \code{\link[lbmech]{defineWorld}}. 
#' Default is \code{tempdir()} but unless the analyses will only be performed a few
#' times it is highly recommended to define a permanent workspace.
#' \code{\link[lbmech]{whichTiles}}---containing the unique tile IDs for sectors that
#' should be in the workspace. Default is \code{NULL}.
#' @param ... Additional parameters to pass to \code{costFUN}
#' @param costname A name to save the cost call parametrs. Default is the name 
#' of the costFUN variable.
#' @return An \code{.fst} file for each sector named after its sector id
#' stored in the \code{/World/Diff} directory, and/or a data.table object (depending
#' on the output parameter) containing a data.table with at least eight columns
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
#' (6) \code{$dl} The planimetric distance between the \code{from} and \code{to} cells
#' 
#' (7) \code{$dr} The 3D distance between the \code{from} and \code{to} cells
#' 
#' (8+) The output cost variables
#' @importFrom terra vect
#' @importFrom data.table transpose
#' @importFrom data.table :=
#' @importFrom data.table melt
#' @importFrom data.table fifelse
#' @importFrom data.table as.data.table
#' @importFrom data.table data.table
#' @importFrom terra terrain
#' @importFrom terra atan2
#' @importFrom terra nlyr
#' @importFrom terra classify
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
#' # Calculate the costs within the world
#' calculateCosts(tiles = tiles, dir = dir,
#' m = 70, v_max = 1.5, BMR = 76, k = 3.5, s = 0.05, l_s = 1,
#' L = 0.8)
#' @export
calculateCosts <- function(tiles = NULL, costFUN = energyCosts, dir = tempdir(),
                           costname = deparse(substitute(costFUN)), ...){
  # This bit is to silence the CRAN check warnings for literal column names
  R0=R45=..dirConvert=R90=R135=R180=R235=R270=R315=cell=flow_theta=WaterSpeed=NULL
  y=x=..z_fix=water_speed=FlowSpeed=from=z_i=z_f=x_i=y_i=x_f=y_f=dz=dl=dr=NULL
  neighbor_distance=uv=z_min=dist=long_i=lat_i=long_f=lat_f=r=f=b=priority=.I=NULL
  
  # Set up the local environment from directory
  dir <- normalizePath(paste0(dir,"/World"),mustWork=FALSE)
  subdirs <- c("/Raw","/Local","/Diff")
  subdirs <-  normalizePath(paste0(dir,subdirs),mustWork=FALSE)
  
  callVars <- readRDS(normalizePath(paste0(dir,"/callVars.gz"),mustWork=FALSE))
  list2env(lapply(as.list(callVars),unlist),environment())
  
  source <- vect(normalizePath(paste0(dir,"/z_sources.gpkg"),mustWork=FALSE))
  grid <- vect(normalizePath(paste0(dir,"/z_grid.gpkg"),mustWork=FALSE))
  z_fix <- importRST(normalizePath(paste0(dir,"/z_fix"),mustWork=FALSE))
  
  do_water <- FALSE
  if (file.exists(normalizePath(paste0(dir,"/water.gpkg"),mustWork = FALSE))){
    waterAll <- vect(normalizePath(paste0(dir,"/water.gpkg")))
    do_water <- TRUE
  } else if (file.exists(normalizePath(paste0(dir,"/water.tif"),mustWork=FALSE))){
    waterAll <- rast(normalizePath(paste0(dir,"/Water.tif")))
    do_water <- TRUE
  }
  # Save parameters as a file with the cost function name
  funName <- costname
  
  costPath <- paste0(dir,"/",funName,".gz")
  costPath <- normalizePath(costPath, mustWork = FALSE)
  
  if (!file.exists(costPath)){
    params <- list(...)
    saveRDS(params,costPath)
  } else {
    params <- readRDS(costPath)
    if (length(list(...) )> 0){
      warning(paste0(costPath," already defined. Ignoring input parameters."))
    }
  }
  
  for (i in tiles){
    water <- FALSE
    tilePath <- normalizePath(paste0(subdirs[3],"/",i,".fst"),
                              mustWork = FALSE)
    doTile <- FALSE
    # Create sector differential if it doesn't exist and calculate
    # cost columns before exporting
    if (!file.exists(tilePath)){
      DT <- makeWorld(tiles = i, 
                      dir = normalizePath(stringr::str_remove(dir,'World$'),
                                          mustWork=FALSE), 
                      output = 'object')
      
      params$water <- FALSE
      preNames <- names(DT)
      DT <- cbind(DT,
                  do.call(costFUN,c(list(DT),params)))
      DT <- stats::na.omit(DT)
      doTile <- TRUE
      
    } else if (sum(c("dt","dU_l","dK_l","dW_l","dE_l") %in% 
                   fst::metadata_fst(tilePath)$columnNames) != 5){
      # If they do exist, check to make sure that the relevant column names
      # are present. If not, import, calculate columns, and overwrite
      
      DT <- fst::read_fst(tilePath,
                          as.data.table = TRUE)
      preNames <- names(DT)
      DT <- cbind(DT,
                  do.call(costFUN,c(list(DT),params)))
      DT <- stats::na.omit(DT)
      doTile <- TRUE
    }
    
    if (do_water && doTile){
      # If the provided water file is a polygon vector, use flow direction
      # from elevation raster and speed from attribute
      waterHere <- FALSE
      if (methods::is(waterAll,"SpatVector")){
        if (length(intersect(grid[which(grid[['id']] == i), ], waterAll)) > 0){
          waterHere <- 1
          # Import terrain, calculate direction of flow
          flowdir <- suppressWarnings(importRST(normalizePath(paste0(subdirs[2],"/",i),mustWork = FALSE)))
          z_fix <- flowdir
          flowdir <- terrain(flowdir,'flowdir')
          
          flowdir <- mask(flowdir,waterAll)
          water <- rasterize(waterAll,flowdir,field=water_speed)
          names(water) <- 'WaterSpeed'
        } 
      } else if (methods::is(waterAll,"SpatRaster")){
        # If water is a raster and it intersects the current grid extent
        
        water <- tryCatch(crop(waterAll,buffer(grid[which(grid[['id']] == i)],
                                               neighbor_distance/max(res(waterAll))),
                               mask=TRUE, snap = 'out'), 
                          error = function(x) FALSE)
        if (methods::is(water,'SpatRaster')){
          waterHere <- nlyr(water)
          water <- project(water,z_fix,align=TRUE)
          if (waterHere == 1){
            # If there's only one layer, assume it's speed
            flowdir <- suppressWarnings(importRST(normalizePath(paste0(subdirs[2],"/",i),
                                                                mustWork = FALSE)))
            z_fix <- flowdir
            flowdir <- terrain(flowdir,'flowdir')
            
            flowdir <- tryCatch(mask(flowdir,water), error =function(x) flowdir)
          } else if (waterHere == 2 & uv){
            # If there's two layers and they're uv form, calculate direction and
            # speed rasters
            
            
            flowdir <- atan2(water[[2]],water[[1]])
            flowdir <- round(flowdir / (pi/4)) * (pi/4)
            radianConvert <- data.table(
              rad = seq(-pi,pi,length.out = 9),
              esri = c(16,8,4,4,1,128,64,32,16)
            )
            flowdir <- classify(flowdir, radianConvert)
            names(flowdir) <- 'flowdir'
            water <- sqrt(sum(water^2))
            names(water) <- 'speed'
            waterHere <- TRUE
          } else if (waterHere == 2 & !uv){
            # If there's two layers and they're speed/flowdir form, just reassign
            flowdir <- water[[2]]
            names(flowdir) <- 'flowdir'
            water <- water[[1]]
            names(water) <- 'speed'
          }
        }
      }
      
      if (!methods::is(water,'logical')){
        # Get adjacency list, convert list into the cell name of the 
        # direction of flow
        flowadj <- as.data.table(adjacent(flowdir,cells(flowdir),
                                          pairs=FALSE,directions = 8),
                                 keep.rownames = TRUE)
        if (nrow(flowadj) != 0){
          
          names(flowadj) <- c("cell",1:8)
          flowadj[, cell := cells(flowdir)]
          flowadj$cell <- getCoords(xyFromCell(flowdir,as.numeric(flowadj$cell)),
                                    z_fix=z_fix)
          flowadj[, cell := as.character(cell)]
          
          # Lookup table for esri notation vs. different orientations
          dirConvert <- data.table(esri = c(8,4,2,16,1,32,64,128),
                                   R0 =   c(1,2,3,4, 5, 6, 7,  8),
                                   R45=   c(2,3,4,5,6,7,8,1),
                                   R90=   c(3,4,5,6,7,8,1,2),
                                   R135=  c(4,5,6,7,8,1,2,3),
                                   R180=  c(5,6,7,8,1,2,3,4),
                                   R235=  c(6,7,8,1,2,3,4,5),
                                   R270=  c(7,8,1,2,3,4,5,6),
                                   R315=  c(8,1,2,3,4,5,6,7))
          
          # Convert from ESRI to matrix index
          flowdir_rast <- copy(flowdir)
          flowdir[cells(flowdir)] <- dirConvert$R0[match(
            flowdir[cells(flowdir)]$flowdir,dirConvert$esri)]
          dirConvert[,names(dirConvert) := lapply(.SD,as.character)]
          
          # Convert flow direction to data.table
          flowdir <- data.table(cells = cells(flowdir),
                                flowdir = unlist(flowdir[cells(flowdir)]))
          names(flowdir) <- c('cell','flowdir')
          flowdir$cell <- getCoords(xyFromCell(flowdir_rast,
                                               as.numeric(flowdir$cell)),
                                    z_fix=z_fix)
          flowdir[, cell := as.character(cell)]
          
          # Combine the adjacency and direction into one table
          flowadj <- merge(flowadj,flowdir,by='cell')
          rm(flowdir)
          
          # Assign cell names to adjacency list (finally)
          flowadj[, R0 := .SD[[flowdir]],by='cell'
          ][, R45 := .SD[[..dirConvert$R45[as.numeric(flowdir)]]],by='cell'
          ][, R90 := .SD[[..dirConvert$R90[as.numeric(flowdir)]]],by='cell'
          ][, R135 := .SD[[..dirConvert$R135[as.numeric(flowdir)]]],by='cell'
          ][, R180 := .SD[[..dirConvert$R180[as.numeric(flowdir)]]],by='cell'
          ][, R235 := .SD[[..dirConvert$R235[as.numeric(flowdir)]]],by='cell'
          ][, R270 := .SD[[..dirConvert$R270[as.numeric(flowdir)]]],by='cell'
          ][, R315 := .SD[[..dirConvert$R315[as.numeric(flowdir)]]],by='cell']
          
          # Melt the table from wide to long
          flowadj <- flowadj[,.(from = cell,R0,R45,R90,R135,R180,R235,R270,R315)]
          flow <- melt(flowadj,id.vars='from',variable.name='flow_theta',
                       value.name ='to')
          
          # Assign geometric costs
          flow[flow_theta == "R0",   WaterSpeed := 1
          ][flow_theta == "R45",  WaterSpeed := sqrt(2)/2
          ][flow_theta == "R90",  WaterSpeed := 0
          ][flow_theta == "R135", WaterSpeed := -sqrt(2)/2
          ][flow_theta == "R180", WaterSpeed := -1
          ][flow_theta == "R235", WaterSpeed := -sqrt(2)/2
          ][flow_theta == "R270", WaterSpeed := 0
          ][flow_theta == "R315", WaterSpeed := sqrt(2)/2
          ]
          
          # Convert cells to coordinate names
          to <- as.data.table(xyFromCell(flowdir_rast,as.numeric(flow$to)))
          to[!is.na(y) & !is.na(x), to:= getCoords(.SD,z_fix=..z_fix)]
          
          flow$to <-   to$to
          rm(to)
          
          flow[,`:=`(flow_theta=NULL)]
          
          flow <- merge(flow,
                        data.table(from=getCoords(xyFromCell(water,cells(water)),z_fix = z_fix),
                                   FlowSpeed = unlist(water[cells(water)])), by = 'from')
          flow[,`:=`(WaterSpeed = WaterSpeed * FlowSpeed, FlowSpeed = NULL)]
          
          
          flow[, (c("x_i","y_i")) := lapply(tstrsplit(from,","),as.numeric)
          ][, (c("x_f","y_f")) := lapply(tstrsplit(to,","),as.numeric)
          ][, (c("z_i","z_f","dz"))  := list(z_min,z_min,0)]
          
          flow <- stats::na.omit(flow)
          
          if (dist == 'proj'){
            flow[, dl := sqrt((x_f - x_i)^2 + (y_f - y_i)^2)]
          } else {
            flow[, c("long_i","lat_i") :=
                   as.data.table((project(as.matrix(data.table(x=x_i,y=y_i)),
                                          from = crs(z_fix), to = "+proj=longlat")))
            ][, c("long_f","lat_f") :=
                as.data.table((project(as.matrix(data.table(x=x_f,y=y_f)),
                                       from = crs(z_fix), to = "+proj=longlat")))]
            
            if (dist == 'karney'){
              flow[, dl := geosphere::distGeo(data.table(x = long_i, y = lat_i),
                                              data.table(x = long_f, y = lat_f), 
                                              a = r, f = f)]
              
            } else if (dist == 'cosine'){
              flow[, dl := geosphere::distCosine(data.table(x = long_i, y = lat_i),
                                                 data.table(x = long_f, y = lat_f), 
                                                 r = r)]
            } else if (dist == 'haversine'){
              flow[, dl := geosphere::distHaversine(data.table(x = long_i, y = lat_i),
                                                    data.table(x = long_f, y = lat_f), 
                                                    r = r)]
            } else if (dist == 'meeus'){
              flow[, dl := geosphere::distMeeus(data.table(x = long_i, y = lat_i),
                                                data.table(x = long_f, y = lat_f), 
                                                a = r, f = f)]
            } else if (dist == 'vincentyEllipsoid'){
              flow[, dl := geosphere::distVincentyEllipsoid(data.table(x = long_i, y = lat_i),
                                                            data.table(x = long_f, y = lat_f), 
                                                            a = r, b = b, f = f)]
            } else if (dist == 'vincentySphere'){
              flow[, dl := geosphere::distVincentySphere(data.table(x = long_i, y = lat_i),
                                                         data.table(x = long_f, y = lat_f), 
                                                         r = r)]
            } else{
              stop("Unknown distance calculation method")
            }
            flow[, c("long_i","lat_i","long_f","lat_f") := NULL
            ]
          }
          
          flow[, dr := dl]
          
          params$water <- TRUE
          flow <- cbind(flow,
                        do.call(costFUN,c(list(flow),params)))
          flow <- stats::na.omit(flow)
          
          if (priority == 'water'){
            DT <- rbind(flow[,.SD,.SDcols = names(DT)], DT)
          } else if (priority == 'land'){
            DT <- rbind(DT, flow[,.SD,.SDcols = names(DT)])
          }
        }
        
        DT <- merge(
          DT[,.(x_i = unique(x_i),
                y_i = unique(y_i),
                dz = fifelse(max(dz,na.rm = TRUE) == z_min, 
                             min(dz,na.rm = TRUE), 
                             max(dz,na.rm = TRUE))), 
             by = c('from','to')],
          DT[DT[,.I[1],by=c('from','to')]$V1, .SD,
             .SDcols = names(DT)[!(names(DT) %in% c("x_i","y_i","dz"))]],
          by = c('from','to')
        )
      }
    } 
    
    if (doTile){
      fst::write_fst(stats::na.omit(DT),tilePath)
    }
    
  } 
} 
