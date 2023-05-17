#' A function that calculates the cost to travel between two sets of points.
#' This can be between two circumscribed sets of points, or one circumscribed one
#' ('nodes') and all other points on the landscape.
#'
#' There are four possible workflows:
#'
#' (1) If you simply desire the distance between two sets of points, provide
#' entries for \code{from} and \code{to} (or just \code{from} if the interest is
#' in all distances between locations in that object). Output is a distance matrix.
#' The computational time for this operation is comparable to generating a raster
#' for the distance to \emph{all} cells in the world (unless all of the locations
#' in the object are close to each other). So unless the operation is to be done
#' multiple times, it is highly recommended to generate the rasters as below and extract
#' values
#'
#' (2) If you wish to generate a RasterStack of costs from and/or to all nodes
#' in the \code{from} object, set the \code{output = 'object'} and 
#' \code{destination = 'all'}.
#'
#' (3) You may also save the rasters as a series of \code{.tif} files in the same workspace
#' directory as the transition \code{.gz} tensor files and the cropped/downloaded
#' DEMs. This allows us to use \code{getCosts} within a loop for large numbers of origin
#' nodes without running into random access memory limitations. Do this by
#' setting \code{output = 'file'} and \code{destination = 'all'}.
#'
#' (4) You may perform (2) and (3) simultaneously by setting
#' \code{output == c('file','object')} and \code{destination = 'all'}.
#'
#' @title Get cost of travel
#' @param region A SpatVector, Spatial* or Raster* representing the area of 
#' maximum movement
#' @param from One of data.frame, data.table, SpatVector, SpatialPointsDataFrame, or
#' SpatialPolygonsDataFrame representing the origin locations. If \code{to = NULL}
#' and \code{destination = 'pairs'}, this will also be used as the \code{to} parameter.
#' If it is a polygon, the location will be controlled by the \code{polygons} parameter.
#' @param to One of data.frame, data.table, SpatVector SpatialPointsDataFrame, or
#' SpatialPolygonsDataFrame representing the origin locations. Optional if
#' \code{destination =  'pairs'}, ignored if \code{destination =  'all'}.
#' @param id Character string representing the column containing the unique
#' IDs for for each location in the \code{from} (and \code{to}) objects.
#' @param x A character vector representing the column containing the 'x' coordinates.
#' Required if \code{data} is not Spatial*.
#' @param y A character vector representing the column containing the 'y' coordinates.
#' Required if \code{data} is not Spatial*.
#' @param destination One of \code{'pairs'} or \code{'all'}. If \code{'pairs'},
#' a distance matrix will be generated between every pair of locations in
#' \code{from}, or every pair of locations between \code{from} and \code{to}.
#' If \code{'all'}, rasters will be generated for each node representing the cost
#' to travel to every cell in the given \code{world}.
#' @param polygons One of \code{c('polygon','centroid','center')}. Ignored unless
#' \code{from} and/or \code{to} are polygons. If \code{polygons = 'centroid'} (the default),
#' the destinations are calculated to the centroid of the polygon whether or not
#' it lies within the polygon itself. If \code{polygons = 'center'}, distances
#' are calculated to the point within the polygon closest to the centroid. If 
#' \code{polygons = 'polygons'}, distances are calculated to any point within the 
#' polygon---in essence, the polygon acts as a giant node permitting costless
#' movement within its bounds. This is generally not consistent with real-world
#' scenarios and for that reason is not the default. 
#' @param costs A character vector containing any combination of the strings
#' of differential values present in the environment (see 
#' \code{\link[lbmech]{calculateCosts}} function; default is 
#' \code{costs = c("dt","dW_l","dE_l")} anticipating the use of 
#' \code{\link[lbmech]{calculateCosts}(costFUN = \link[lbmech]{energyCosts})}.
#' @param direction A character vector containing one or both of \code{c("in","out")}
#' or the singular string 'both'. This determines whether costs to or from the nodes
#' are calculated. Ignored for \code{destination = 'pairs'}.
#' @param output A character vector containing one or both of \code{c("object","file")}.
#' If \code{"object"} is included, then a list of RasterStacks will be returned.
#' If \code{"file"} is included, then the appropriate cost rasters will be saved
#' to the workspace directory \code{dir}. Ignored if \code{destination = 'pairs'}.
#' @param outname A character vector describing the name of the set of input
#' nodes for raster analysis. Ignored unless \code{'file' \%in\% output}, in
#' which case the files will be stored in a file named as such. Default is the
#' name of the \code{from} input, which can lead to files being 
#' overwritten if vector sources are arbitrarily changed.
#' @param overwrite Should any output files with the same outname be overwritten?
#' Default is \code{overwrite = FALSE}. 
#' @param dir A filepath to the directory being used as the workspace.
#' Default is \code{tempdir()} but unless the analyses will only be performed a few
#' times it is highly recommended to define a permanent workspace.
#' @param ... Additional arguments to pass to \code{\link[lbmech]{importWorld}}.
#' @return A list, with entries for each combination of selected \code{costs}
#' and \code{directions}. The object class of the list entries depends on the
#' \code{destination} and \code{output} parameters:
#'
#' (1) If \code{destination = 'pairs'}, entries are distance matrices.
#' (2) If \code{destination = 'all'} and \code{'object' \%in\% output},
#' entries are RasterStacks with each Raster* representing the cost to/from
#' each node.
#' (3) If \code{destination = 'all'} and \code{output == 'file'}, no list is output.
#'
#' Moreover, if \code{file \%in\% output}, then the cost rasters are saved in the
#' workspace directory.
#' @importFrom data.table as.data.table
#' @importFrom data.table :=
#' @importFrom data.table data.table
#' @importFrom data.table tstrsplit
#' @importFrom terra writeRaster
#' @importFrom terra writeVector
#' @importFrom terra rast
#' @importFrom terra crs
#' @importFrom terra crs<-
#' @importFrom terra res
#' @importFrom terra project
#' @importFrom terra geomtype
#' @importFrom terra centroids
#' @importFrom terra geom
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
#' @examples 
#' #### Example 1:
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
#' region <- grid[8,]
#' # Generate five random points that fall within the region
#' points <- data.table(ID = 1:5,
#'                      x = runif(5, ext(region)[1], ext(region)[2]),
#'                      y = runif(5, ext(region)[3], ext(region)[4]))
#'                      
#' # Make a world but limit it to the DEM grid size
#' defineWorld(source = grid, cut_slope = 0.5, 
#'             res = res(dem), dir = dir, overwrite=TRUE)
#'             
#' # Get time and work costs between points
#' costMatrix <- getCosts(grid[8,], from = points, dir = dir, 
#'                        costs = c('dt','dW_l'), costFUN = energyCosts,
#'                        m = 70, v_max = 1.5, BMR = 76, k = 3.5, s = 0.05, l_s = 1,
#'                        L = 0.8)
#'                    
#' #### Example 2:
#' # Calculate the cost rasters to travel to ad from the center of a polygon
#' costRasters <- getCosts(grid[8,], from = grid[8,], dir = dir, destination = 'all',
#'                        polygons = 'center',
#'                        costs = c('dt','dW_l'), costFUN = energyCosts,
#'                        m = 70, v_max = 1.5, BMR = 76, k = 3.5, s = 0.05, l_s = 1,
#'                        L = 0.8)
#' 
#' @export
getCosts <- function(region, from, to = NULL, id = 'ID', dir = tempdir(),
                     x = "x", y = "y", destination = 'pairs', polygons = 'centroid',
                     costs = c("dt","dW_l","dE_l"), direction = 'both', output = 'object',
                     outname = deparse(substitute(from)), overwrite = FALSE,
                     ...){
  # This bit is to silence the CRAN check warnings for literal column names
  ..id=..x=..y=dt=dW_l=dE_l=Var2=value=Var1=coord=Replace=Masked=x_n=y_n=NULL
  Cell=From_ID=ID=..cost=..d=Value=x_i=y_i=Vector=proj=precision=V1=filt=NULL
  #
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
  
  if(file.exists(normalizePath(paste0(dir,"/CostRasters/",
                                      outname,'.gpkg'),
                               mustWork=FALSE)) & !overwrite){
    stop(paste(outname,"already exists. To overwrite, use 'overwrite=TRUE'. 
             Alternatively, set output = 'object' to ignore file creation.
             Carefuly ensure that this is the desired behavior."))
  } 
  
  id2 <- 1
  # If no Unique ID column is given, make one
  if (is.null(id)){
    from$Node_ID <- seq(1,nrow(from))
    id <- "Node_ID"
    id2 <- NULL
  }

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
  
  # We need to coerce everything to an (x,y) list stored as a data.table
  # First, if the input are polygons we'll just find the centroid
  fromMask <- NULL
  toMask <- NULL
  polysmask <- NULL
  if (methods::is(from,'Spatial')){ 
    from <- vect(from)
  }
  if (methods::is(from, "SpatVector")){
    if (geomtype(from) == 'polygons'){
      if (polygons == 'border'){
        fromMask <- regionMask(from, z_fix = z_fix, id = id)
        from <- centroids(from,inside = TRUE)
      } else if (polygons == 'centroid'){
        from <- centroids(from,inside = FALSE)
      } else if (polygons == 'center'){
        from <- centroids(from,inside = TRUE)
      }
    }
    
    from[[x]] <- as.data.table(geom(from))[,x]
    from[[y]] <- as.data.table(geom(from))[,y]
  }
  
  # Coerce to data,.table
  from <- as.data.table(from)[,.(ID = get(..id), x = get(..x), y = get(..y))]
  from$Cell <- getCoords(from, z_fix = z_fix, precision = precision)
  if ('file' %in% output){
    fromVect <- vect(from, geom=c('x','y'), crs = crs(z_fix))
  }
  
  # Do the same for the destination UNLESS it's set as "all".
  # If it's empty, then from is the same as to.
  if (destination != "all"){
    if (!is.null(to)){
      if (methods::is(to,'Spatial')){
        to <- vect(to)
      }
      if (methods::is(to,"SpatVector")){
        # If no column ID is provided, make one as above
        if (is.null(id2)){
          to$Node_ID <- seq(1,nrow(to))
        }
        if (geomtype(from) == 'polygons'){
          toMask <- regionMask(to, z_fix = z_fix, id = id)
        }
        to <- centroids(to,inside = TRUE)
        to[[x]] <- as.data.table(geom(to))[,x]
        to[[y]] <- as.data.table(geom(to))[,y]
      }
      
      to <- as.data.table(to)[,.(ID = get(..id), x = get(..x), y = get(..y))]
      to$Cell <- getCoords(to, z_fix = z_fix, precision = precision)
    } else {
      to <- from
      toMask <- fromMask
    }
  } else{
    # If destination is set as "all", then to is "NULL" so that igraph
    # calculaltes the distance to all cells.
    if (!is.null(to)){
      warning("Destinations have been provided, but method selected is 'all'.
            Ignoring input for 'to'.")
      Sys.sleep(1)
      to <- NULL
    }
  }
  
  # Join the fromMask and toMask (if present) into a polysMask that keeps
  # the relations between the masked nodes and their neighbors
  if (!is.null(fromMask)){
    fromMask <- merge(fromMask,from,by='ID'
    )[,`:=`(Masked = coord,
            Replace = Cell,
            coord = NULL,
            Cell = NULL)][]
    polysmask <- fromMask
  } 
  if (!is.null(toMask)){
    toMask <- merge(toMask,to,by='ID'
    )[,`:=`(Masked = coord,
            Replace = Cell,
            coord = NULL,
            Cell = NULL)][]
    polysmask <- rbind(polysmask,toMask)
  }
  if (!is.null(polysmask)){
    polysmask[, (id) := NULL]
    polysmask <- unique(polysmask)
    
    # Get the location of the nodes-to-be-masked
    polysmask[, c("x_n","y_n") := tstrsplit(Replace,",")
    ][, c("x_n","y_n") := lapply(.SD,as.numeric),.SDcols = c("x_n","y_n")
    ]
    
    # Condense the observation for all nodes within a polygon to a single node
    w <- region[from %in% polysmask$Masked | to %in% polysmask$Masked]
    w[from %in% polysmask$Masked, 
      c("from","x_i","y_i") := polysmask[Masked %in% from, 
                                         .(Replace, x_n, y_n)][1]
    ][to %in% polysmask$Masked, 
      to := polysmask[Masked %in% to, .(Replace)][1]
    ][,(names(world)[!(names(world) %in% c("from","to"))]) :=
        lapply(.SD,min),by=c("from","to")]
    w <- unique(w)
    
    region <- rbind(region[!(from %in% polysmask$Masked | to %in% polysmask$Masked)],
                    w)
    rm(w)
  }
  
  
  # If direction is set to both, we'll need to create a "to"
  # and "from" category
  if (all(direction == "both")){
    direction <- c("in","out")
  }
  
  # Create a list object to which we'll append either the matrices or
  # rasterstacks representing the costs
  outList <- list()
  
  pb <- txtProgressBar(min = 0, max = length(costs) * length(direction), style=3)
  iter = 0
  
  region[,`:=`(x_i = NULL, y_i =NULL,dumval = 1)]
  # We first iterate over the desired costs. Each desired cost will
  # create one or two entries (depending on the number of directions)
  # in the output list.
  for (cost in costs){
    
    # For each loop, we'll generate the igraph object using the
    # appropriate cost variable
    world <- merge(region,
                   importWorld(region_shp,vars = cost, 
                               dir = normalizePath(stringr::str_remove(dir,'World$'),
                                                   mustWork=FALSE),
                               filt = filt),
                   by = c("from","to"), all = FALSE, allow.cartesian=TRUE)
    world <- stats::na.omit(world)
    world$dumval <- NULL
    world <- world[, min(get(..cost),na.rm=TRUE), by=c('from','to')]
    Graph <- igraph::graph_from_data_frame(world[,.(from,to,weight = V1)])
    
    
    costVals <- data.table()
    
    # Now we also iterate over each desired direction
    for (d in direction){
      
      # This section is if the desired output is a *matrix*
      # and we are not storing the rasters
      if (all(destination == 'pairs')){
        # Calculate the distances
        Distances <- igraph::distances(Graph,
                                       v = from$Cell,
                                       to = to$Cell,
                                       mode = d)
        
        # Convert to data.table
        Distances <- as.data.table(reshape2::melt(Distances))
        
        # Append the appropriate names/ID
        Distances <- merge(Distances[,.(Var2,value,Cell = Var1)],
                           from[,.(Cell,From_ID = ID)],
                           on="Cell",allow.cartesian=TRUE)[, Cell := NULL][]
        
        Distances <- merge(Distances[,.(From_ID,value,Cell = Var2)],
                           to[,.(Cell,To_ID = ID)],
                           on="Cell",allow.cartesian=TRUE)[, Cell := NULL][]
        
        # Convert to distance matrix with xtabs,
        # Add to the output list
        outList[[paste0(cost,"_",d)]] <- stats::xtabs(value~.,Distances)
        
      }
      
      # This section is if the desired costs are from a set of
      # cells to every location on the landscape (a raster of costs)
      if (all(destination == 'all')){
        # Calculate the distances
        Distances <- igraph::distances(Graph,
                                       v = from$Cell,
                                       mode = d)
        
        # Convert to data.table
        Distances <- as.data.table(reshape2::melt(Distances))
        
        # Append the appropriate IDs to each origin/destination node
        # then get the x,y coordinates from the cell name
        Distances <- merge(Distances[,.(Var2,value,Cell = Var1)],
                           from[,.(Cell,ID = ID)],
                           on="Cell", allow.cartesian=TRUE)[, Cell := NULL
                           ][,.(ID,Cell = Var2, Cost = ..cost,
                                Direction = ..d,Value = value)]
        Distances[, c("x","y") := tstrsplit(Cell,",")
        ][, c("x","y") := lapply(.SD,as.numeric), .SDcols = c("x","y")
        ]
        
        if (!is.null(polysmask)){
          Distances <- rbind(Distances,
                             merge(Distances, polysmask[,.(Masked,Replace)],
                                   by.x='Cell', by.y='Replace',
                                   allow.cartesian=TRUE
                             )[!is.na(Masked), Cell := Masked
                             ][, Masked := NULL
                             ][, c("x","y") := tstrsplit(Cell,",")
                             ][, c("x","y") := lapply(.SD,as.numeric), .SDcols = c("x","y")][])
        }
        
        # We'll add all of the generated rasters to a stack
        outStack <- rast()
        for (i in unique(Distances$ID)){
          # Convert the data.table to raster and add to the stack
          outStack <- suppressWarnings(c(outStack,
                                         rast(Distances[ID == i & 
                                                          !is.infinite(Value)
                                                        ,.(x,y,Value)],
                                              crs=crs(z_fix))))
        }
        
        # Select the appropriate prefix names for the filepaths/layer names
        if (all(d == "out")){
          names(outStack) <- paste0(cost,"_from_",unique(Distances$ID))
        } else if (all(d == "in")){
          names(outStack) <- paste0(cost,"_to_",unique(Distances$ID))
        }
        
        # Add to raster stack
        
        outList[[paste0(cost,"_",d)]] <- outStack
        
        
        
      }
      iter <- iter + 1
      setTxtProgressBar(pb,iter)
    }
    
  }
  
  if (methods::is(outList[[1]], 'SpatRaster')){
    if (length(outList) > 1){
      names(outList) <- NULL
    }
    outList <- rast(outList)
    outList <- outList[[gtools::mixedsort(names(outList))]]
  } else {
    if (length(outList) == 1){
      outList <- unlist(outList)
    }
  }
  
  if ("file" %in% output){
    outPath <- normalizePath(paste0(dir,"/CostRasters/"),mustWork=FALSE)
    if (!dir.exists(outPath)){
      dir.create(outPath)
    }
    writeRaster(outList,
             filename=normalizePath(paste0(outPath,"/",
                                           outname,'.tif'),
                                    mustWork=FALSE))
    
    writeVector(fromVect, normalizePath(paste0(outPath,"/",
                                               outname,'.gpkg'),
                                        mustWork=FALSE), overwrite = TRUE)
  
  }
  if ("object" %in% output){
    return(outList)
  }
  
}