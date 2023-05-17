#' Get the shortest path for a given trip that requires travel through a
#' set of nodes. Use is like \code{\link[lbmech]{getCosts}}, but with 
#' \code{nodes} and {order} parameters and no \code{from} or \code{to}.
#'
#' @title Get least-cost paths
#' @param region A SpatVector, Spatial* or Raster* representing the area of 
#' maximum movement
#' @param nodes One of data.frame, data.table, SpatVector, SpatialPointsDataFrame, or
#' SpatialPolygonsDataFrame representing the node locations. 
#' If it is a polygon, the location will be controlled by the \code{polygons} parameter.
#' @param id A character string representing the column containing each \code{node}
#' location's unique ID.
#' @param order A character vector containing the desired path in
#' order of visited nodes by ID. Required if  For example, to visit "A" then "B" then "C" then "A"
#' the vector would be \code{c("A","B","C","A")}. If this is not provided, 
#' the function assumes that \code{nodes} is  already sorted in the desired order.
#' @param x A character vector representing the column containing the 'x' coordinates.
#' Required if \code{data} is not Spatial*.
#' @param y A character vector representing the column containing the 'y' coordinates.
#' Required if \code{data} is not Spatial*.
#' @param costs A character vector containing any combination of the strings
#' of differential values present in the environment (see 
#' \code{\link[lbmech]{calculateCosts}} function; default is 
#' \code{costs = c("dt","dW_l","dE_l")} anticipating the use of 
#' \code{\link[lbmech]{calculateCosts}(costFUN = \link[lbmech]{energyCosts})}.
#' @param polygons One of \code{c('polygon','centroid','center')}. Ignored unless
#' \code{nodes} are polygons. If \code{polygons = 'centroid'} (the default),
#' the destinations are calculated to the centroid of the polygon whether or not
#' it lies within the polygon itself. If \code{polygons = 'center'}, distances
#' are calculated to the point within the polygon closest to the centroid. If 
#' \code{polygons = 'polygons'}, distances are calculated to *any* point within the 
#' polygon---in essence, the polygon acts as a giant node permitting costless
#' movement within its bounds. This is generally not consistent with real-world
#' scenarios and for that reason is not the default. 
#' @param dir A filepath to the directory being used as the workspace.
#' Default is \code{tempdir()} but unless the analyses will only be performed a few
#' times it is highly recommended to define a permanent workspace.
#' @param ... Additional parameters to pass to \code{\link[lbmech]{importWorld}}.
#' @return SpatialVector lines representing least-cost paths. For each
#' cost, each entry within the SpatialLinesDataFrame object represents
#' a single leg of the journey, sorted in the original path order.
#' If \code{length(costs) == 1}, only a SpatVector is returned.
#' If \code{length(costs) > 1}
#' a list of SpatVectors with one slot for each \code{cost} is returned.
#' @importFrom data.table as.data.table
#' @importFrom data.table tstrsplit
#' @importFrom data.table data.table
#' @importFrom data.table rbindlist
#' @importFrom data.table :=
#' @importFrom terra vect
#' @importFrom terra as.lines
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
#' region <- grid[8, ]
#' 
#' # Select all tiles that exist between x = (12000,16000) and y = (32000,36000)
#' tiles <- ext(c(12000,16000,32000,36000))
#' tiles <- as.polygons(tiles)
#' crs(tiles) <- crs(grid)
#' tiles <- whichTiles(region = tiles, polys = grid)
#' 
#' # Make a world but limit it to the DEM grid size
#' calculateCosts(tiles = tiles, dir = dir,
#' m = 70, v_max = 1.5, BMR = 76, k = 3.5, s = 0.05, l_s = 1,
#' L = 0.8)
#' 
#' # Generate five random points that fall within the region
#' points <- data.table(ID = 1:5,
#'                      x = runif(5, ext(region)[1], ext(region)[2]),
#'                      y = runif(5, ext(region)[3], ext(region)[4]))
#'                      
#'                      
#' # Calculate the path from 1 -> 2 -> 5 -> 1 -> 4 
#' pathOrder <- c(1,2,5,1,4)
#' 
#' # Make a world but limit it to the DEM grid size
#' defineWorld(source = grid, cut_slope = 0.5, 
#'             res = res(dem), dir = dir, overwrite=TRUE)
#'             
#' paths <- getPaths(region = region, nodes = points, order = pathOrder, 
#'                   costs = 'all', costFUN = energyCosts,
#'                        m = 70, v_max = 1.5, BMR = 76, k = 3.5, s = 0.05, l_s = 1,
#'                        L = 0.8)
#'                               
#' ## Plot against corridors (not run)                         
#' #getCosts(region = region, from = points, proj = crs(dem), res = res(dem),
#' #                        destination = 'all', costs = 'all',
#' #                        output = 'file', dir = dir)                       
#' #corridors <- makeCorridor(rasters = dir, order = pathOrder)
#' #plot(corridors$time)
#' #plot(paths$time,add=TRUE)
#' @export
getPaths <- function(region, nodes, id = "ID", order = NULL, x = "x",
                     y = "y", costs = 'all', polygons = 'centroid',
                     dir = tempdir(),...){
  # This bit is to silence the CRAN check warnings for literal column names
  from=to=..id=..x=..y=ID=Cell=V1=V2=TempID=dt=dW_l=dE_l=x_i=y_i=filt=NULL
  coord=Replace=Masked=x_n=y_n=precision=..cost=NULL
  #
  
  # First, "all" and "energetics" are just shorthand for a vector
  # of multiple cost types. Fix that.
  if (all(costs == 'all')){
    costs <- c("dt","dW_l","dE_l")
  } else if (all(costs == "energetics")){
    costs <- c("dW_l","dE_l")
  }
  
  # Set up the local environment from directory
  dir <- normalizePath(paste0(dir,"/World"),mustWork=FALSE)
  subdirs <- c("/Raw","/Local","/Diff")
  subdirs <-  normalizePath(paste0(dir,subdirs),mustWork=FALSE)
  
  callVars <- readRDS(normalizePath(paste0(dir,"/callVars.gz"),mustWork=FALSE))
  list2env(lapply(as.list(callVars),unlist),environment())
  
  source <- vect(normalizePath(paste0(dir,"/z_sources.gpkg"),mustWork=FALSE))
  grid <- vect(normalizePath(paste0(dir,"/z_grid.gpkg"),mustWork=FALSE))
  z_fix <- importRST(normalizePath(paste0(dir,"/z_fix"),mustWork=FALSE))
  
  region_shp <- region
  region <- importWorld(region,vars = c("x_i","y_i"), 
                        dir = normalizePath(stringr::str_remove(dir,'World$'),
                                            mustWork=FALSE),...)
  
  z_temp <- as.data.table(expand.grid(x = seq(from = min(region[,x_i]) - 2 * res(z_fix)[1],
                                              to = max(region[,x_i]) + 2 * res(z_fix)[1],
                                              by = res(z_fix)[1]),
                                      y = seq(from = min(region[,y_i]) - 2 * res(z_fix)[2],
                                              to = max(region[,y_i]) + 2 * res(z_fix)[2],
                                              by = res(z_fix)[2])))
  z_temp <- rast(z_temp[,.(x,y,z=1)])
  crs(z_temp) <- crs(z_fix)
  z_fix <- suppressWarnings(project(z_temp,z_fix,align = TRUE))
  
  # Deal with polygons such that all values within are considered a single cell
  nodeMask <- NULL
  if (methods::is(nodes,'Spatial')){ 
    nodes <- vect(nodes)
  }
  if (methods::is(nodes,"SpatVector")){
    if (geomtype(nodes) == 'polygons'){
      if (polygons == 'border'){
        fromMask <- regionMask(nodes, z_fix = z_fix, id = id)
        nodes <- centroids(nodes,inside = TRUE)
      } else if (polygons == 'centroid'){
        nodes <- centroids(nodes,inside = FALSE)
      } else if (polygons == 'center'){
        nodes <- centroids(nodes,inside = TRUE)
      }
    }
    
    nodes[[x]] <- as.data.table(geom(nodes))[,x]
    nodes[[y]] <- as.data.table(geom(nodes))[,y]
  }
  
  # Coerce to data.table
  nodes <- as.data.table(nodes)[,.(ID = get(..id), x = get(..x), y = get(..y))]
  nodes$Cell <- getCoords(nodes, z_fix = z_fix, precision = precision, ...)
  
  # Consolidate all cells within the world for any polygon to their centroid
  if (!is.null(nodeMask)){
    nodeMask <- merge(nodeMask,nodes,by='ID'
    )[,`:=`(Masked = coord,
            Replace = Cell,
            coord = NULL,
            Cell = NULL)][]
    
    nodeMask[, (id) := NULL]
    nodeMask <- unique(nodeMask)
    nodeMask[, c("x_n","y_n") := tstrsplit(Replace,",")
    ][, c("x_n","y_n") := lapply(.SD,as.numeric),.SDcols = c("x_n","y_n")
    ]
    
    w <- world[from %in% nodeMask$Masked | to %in% nodeMask$Masked]
    w[from %in% nodeMask$Masked, 
      c("from","x_i","y_i") := nodeMask[Masked %in% from, 
                                        .(Replace, x_n, y_n)][1]
    ][to %in% nodeMask$Masked, 
      to := nodeMask[Masked %in% to, .(Replace)][1]
    ][,(names(world)[!(names(world) %in% c("from","to"))]) :=
        lapply(.SD,min),by=c("from","to")]
    w <- unique(w)
    
    world <- rbind(world[!(from %in% nodeMask$Masked | to %in% nodeMask$Masked)],
                   w)
    rm(w)
    
  }
  
  # If no order is given, assume the node object is already in the desired order
  if (is.null(order)){
    order <- unlist(nodes[[id]])
  }
  
  # Get the IDs of nodes that will be needed for each leg
  starts <- order[seq(1,(length(order)-1))]
  stops <-  order[seq(2,length(order))]
  order <- paste0(starts,"_to_",stops)
  
  # We only need to do the calculations once for each origin node,
  # but we need to keep the order of the above vectors
  origins <- unique(starts)
  
  # We'll store each path inside a list
  pathList <- list()
  
  region[,`:=`(x_i = NULL, y_i =NULL,dumval = 1)]
  # Iterate over every cost, adding each result as an entry in pathList
  for (cost in costs){
    # For each loop, we'll generate the igraph object using the
    # appropriate cost variable
    world <- merge(region,
                   importWorld(region_shp,vars = cost, 
                               dir = normalizePath(stringr::str_remove(dir,'World$'),
                                                   mustWork=FALSE),
                               filt = filt),
                   by = c("from","to"), all = FALSE)
    world <- stats::na.omit(world)
    world$dumval <- NULL
    world <- world[, min(get(..cost),na.rm=TRUE), by=c('from','to')]
    Graph <- igraph::graph_from_data_frame(world[,.(from,to,weight = V1)])
    
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
    # a numeric data.table containing XY coordinates for each traversed cell
    legList <- lapply(legList, igraph::as_ids)
    legList <- lapply(legList, tstrsplit, ",")
    names(legList) <- nameList
    legList <- rbindlist(legList,idcol="Leg")
    legList[,`:=`(x = as.numeric(V1), y = as.numeric(V2), V1 = NULL, V2 = NULL)]
    
    # Convert the XY list to points, then polylines
    lineList <- list()
    for (i in unique(legList$Leg)){
      lineList[[i]] <- as.lines(vect(legList[legList$Leg == i,]
                                     ,geom=c(x,y),crs=crs(z_fix)))
      lineList[[i]]$segment <- i
    }
    
    lineList <- vect(lineList,crs=crs(z_fix))
    
    
    # Add the information about what segment each leg represents
    lineList$cost <- cost
    
    # And place it in the originally-requested order
    lineList <- lineList[match(order, nameList),]
    
    # Save to the output pathList
    pathList[[cost]] <- lineList
    
  }
  if (length(pathList) > 0) {
    if (length(pathList) >= 1 ){
      names(pathList) <- NULL
      pathList <- vect(pathList)
      crs(pathList) <- crs(z_fix)
    }
    return(pathList)
  }
}
