#' Function that defines the grid that can be traversed--the "world"--as well as the
#' cells that can be accessed from each individual cell. This is the most
#' time-intensive function. 
#' 
#' It first checks to see if the required elevation models have been downloaded
#' for each source file for the requested tiles grid, and then if  transition \code{.gz}
#' then if they have been converted to a sector's local
#' files have already been created in the \code{dir} workspace. If not,
#' it generates each at each step
#'
#' @title Build the world from a defined environment
#' @param tiles A character vector--such as the output to
#' \code{\link[lbmech]{whichTiles}}---containing the unique tile IDs for sectors that
#' should be in the workspace. Default is \code{NULL}.
#' @param dir A filepath to the directory being used as the workspace, the same
#' one instantiated with \code{\link[lbmech]{defineWorld}}. 
#' Default is \code{tempdir()} but unless the analyses will only be performed a few
#' times it is highly recommended to define a permanent workspace.
#' @param output A character string or vector, consisting of one or both of 
#' \code{c('file','object')}, representing whether a file should be written and/or
#' whether an object should be returned. Default is \code{output = file}. 
#' @param cols A character vector containing the name of the important 
#' spatial variables to be retained. Default is 
#' \code{cols = c("x_i","y_i","z_i","z_f","dz","dl","dr")},
#' but \code{c("x_f","y_f")} are also available. 
#' @return An \code{.fst} file for each sector named after its sector id
#' stored in the \code{/World/Diff} directory, and/or a data.table object (depending
#' on the output parameter) containing a data.table with five columns
#'
#' (1) \code{$from}, a character string of all possible origin cells in format "x,y",
#' rounded to the next-lowest integer
#'
#' (2) \code{$to},  a character string of all possible destination cells in format "x,y"
#' rounded to the next-lowest integer
#'
#' (3) \code{$dz}, a numeric representing the change in elevation for each origin-destination pair
#' 
#' (4) \code{$dl}, a numeric representing the change in planimetric distance (x,y)
#' 
#' (5) \code{$dr}, a numeric representing the change in displacement (x,y,z)
#' 
#' Likewise, in the \code{/World/Loc} directory the local \code{z} elevation 
#' values projected to the locally defined grid as a \code{\link[lbmech]{writeRST}}  file. 
#' @importFrom terra res
#' @importFrom terra project
#' @importFrom terra mosaic
#' @importFrom terra merge
#' @importFrom terra crs
#' @importFrom terra crs<-
#' @importFrom terra xFromCell
#' @importFrom terra yFromCell
#' @importFrom terra ext
#' @importFrom terra crop
#' @importFrom terra cells
#' @importFrom terra adjacent
#' @importFrom terra ncell
#' @importFrom terra intersect
#' @importFrom terra buffer
#' @importFrom terra vect
#' @importFrom data.table as.data.table
#' @importFrom data.table :=
#' @importFrom data.table .SD
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
#' makeWorld(tiles = tiles, dir = dir)
#' @export 
makeWorld <- function(tiles = NULL, 
                      dir = tempdir(),
                      cols = c("x_i","y_i","dz","dl","dr"),
                      output = 'file'){
  if (length(tiles) == 0){
    tiles <- NULL
  }
  # This bit is to silence the CRAN check warnings for literal column names
  from=to=..dem=z_f=z_i=x_f=x_i=y_f=y_i=dz=..cols=NULL
  cut_slope=source_id=grid_id=directions=neighbor_distance=keepz=unit=vals=NULL
  precision=FUN=sampling=l_p=dist=dl=long_f=lat_f=long_i=lat_i=r=f=b=dr=filt=z_min=NULL
  rm(..cols)
  
  # Import the variable names from the saved files
  dir <- normalizePath(paste0(dir,"/World"),mustWork=FALSE)
  subdirs <- c("/Raw","/Local","/Diff")
  subdirs <-  normalizePath(paste0(dir,subdirs),mustWork=FALSE)
  
  callVars <- readRDS(normalizePath(paste0(dir,"/callVars.gz"),mustWork=FALSE))
  list2env(lapply(as.list(callVars),unlist),environment())
  
  source <- vect(normalizePath(paste0(dir,"/z_sources.gpkg"),mustWork=FALSE))
  grid <- vect(normalizePath(paste0(dir,"/z_grid.gpkg"),mustWork=FALSE))
  z_fix <- suppressWarnings(importRST(normalizePath(paste0(dir,"/z_fix"),mustWork=FALSE)))
  
  loc_tiles <- which(unlist(grid$id) %in% tiles)
  loc_tiles <- grid[loc_tiles,]
  loc_tiles <- buffer(loc_tiles,width = (neighbor_distance * max(l_p) + max(l_p)))
  loc_tiles <- intersect(loc_tiles[,NA],grid)[['id']]
  loc_tiles <- unique(as.character(unlist(loc_tiles)))
  
  for (i in loc_tiles){
    if (!file.exists(normalizePath(paste0(subdirs[2],"/",i,".fst"),mustWork=FALSE))){
      # Use importMap to build a projected DEM, then export it to loc_tiles
      i_poly <- grid[grid$id == i,]
      dem <- importMap(i_poly, polys = source, tile_id = 'source_id', vals = vals,
                       z_fix = z_fix, z_min = z_min, mask = FALSE, filt = filt,
                       dir = normalizePath(paste0(dir,"/Raw/")))
      if (methods::is(dem,'SpatRaster')){
      names(dem) <- 'z'
      }
      writeRST(dem, normalizePath(paste0(subdirs[2],"/",i),mustWork=FALSE))
      
      # Fix the units if they weren't provided in meters
      if (unit == "km"){
        dem <- dem / 1000
      } else if (unit == "ft"){
        dem <- dem*12*2.54/100
      } else if (unit == "mi"){
        dem <- dem*5280*12*2.54/100
      }
    }
  }
  
  
  for (i in tiles){
    if (!file.exists(normalizePath(paste0(subdirs[3],"/",i,".fst"),mustWork=FALSE))){
      # Get the appropriate polygon and name; create a temporary folder
      tensors <- normalizePath(paste0(subdirs[3],"/",i),mustWork=FALSE)
      poly <- which(grid$id == i)
      poly <- grid[poly,]
      
      # import RST
      dem <- suppressWarnings(importRST(normalizePath(paste0(subdirs[2],"/",i),mustWork=FALSE)))
      
      name <- i
      
      nas <- cells(dem)
      tryCatch({
      # Get pairs of adjacent all adjacent cells; drop those that
      # correspond to NA
      adj <- adjacent(dem,nas,directions=directions,pairs=TRUE)
      adj <- as.data.table(adj)
      
      # Calculate the change in elevation between every accessible cell pairs,
      # then drop all values that would require movement over the
      # cut slope.
      adj[, `:=`(z_i = unlist(..dem[from]), z_f = unlist(..dem[to]))
          ][, (c('x_i','y_i')) := as.data.table(xyFromCell(..dem,from))
            ][, (c('x_f','y_f')) := as.data.table(xyFromCell(..dem,to))
      ][, (c("x_i","y_i","x_f","y_f")) :=
          lapply(.SD,round,precision),
        .SDcols = c("x_i","y_i","x_f","y_f")
      ][, `:=`(from = paste(format(x_i, scientific = FALSE),
                            format(y_i, scientific = FALSE), sep=","),
               to = paste(format(x_f, scientific = FALSE),
                          format(y_f, scientific = FALSE), sep=","))
      ][, (c("from","to")) := lapply(.SD,stringr::str_remove_all," "),
        .SDcols = c("from","to")
      ][, `:=`(dz = z_f - z_i)]
      
      if (dist == 'proj'){
        adj[, dl := sqrt((x_f - x_i)^2 + (y_f - y_i)^2)]
      } else {
        adj[, c("long_i","lat_i") :=
              as.data.table((project(as.matrix(data.table(x=x_i,y=y_i)),
                                     from = crs(z_fix), to = "+proj=longlat")))
        ][, c("long_f","lat_f") :=
            as.data.table((project(as.matrix(data.table(x=x_f,y=y_f)),
                                   from = crs(z_fix), to = "+proj=longlat")))]
        if (dist == 'karney'){
          adj$dl <- geosphere::distGeo(adj[,.(long_f,lat_f)],adj[,.(long_i,lat_i)], 
                                       a = r, f = f)
        } else if (dist == 'cosine'){
          adj$dl <- geosphere::distCosine(adj[,.(long_f,lat_f)],adj[,.(long_i,lat_i)],
                                          r = r) 
        } else if (dist == 'haversine'){
          adj$dl <- geosphere::distHaversine(adj[,.(long_f,lat_f)],adj[,.(long_i,lat_i)],
                                             r = r) 
        } else if (dist == 'meeus'){
          adj$dl <- geosphere::distMeeus(adj[,.(long_f,lat_f)],adj[,.(long_i,lat_i)],
                                         a = r, f = f) 
        } else if (dist == 'vincentyEllipsoid'){
          adj$dl <- geosphere::distVincentyEllipsoid(adj[,.(long_f,lat_f)],adj[,.(long_i,lat_i)],
                                                     a = r, b = b, f = f)
        } else if (dist == 'vincentySphere'){
          adj$dl <- geosphere::distVincentySphere(adj[,.(long_f,lat_f)],adj[,.(long_i,lat_i)],
                                                  r = r)
        } else{
          stop("Unknown distance calculation method")
        }
        adj[, c("long_i","lat_i","long_f","lat_f") := NULL
        ]
      }
      
      adj <- stats::na.omit(adj)[, dr := sqrt(dl^2 + dz^2)]
      
      # This exports the cost tensor so that it can be called later by
      # the lbm.distance tool to calculate paths and catchments
      cols <- c("from","to",cols)
      cols <- unique(cols)
      if ('file' %in% output){
        fst::write_fst(adj[,.SD,.SDcols=cols],
                       path = normalizePath(paste0(subdirs[3],"/",i,".fst"),mustWork=FALSE))
      }
      if ('object' %in% output){
        adj <- adj[,.SD,.SDcols=cols]
        if (length(adj) * nrow(adj) == 0){
          rbind(adj, t(rep(NA,length(adj))),use.names=FALSE)
        }
        if (is.null(adj)){
          cols <- c("from","to",cols)
          cols <- unique(cols)
          adj <- 1:length(cols)
          names(adj) <- cols
          adj <- as.data.table(t(adj))
          adj[, (cols) := NA]
        }
        return(adj)
      }
      }, error = function(x){
        cols <- c("from","to",cols)
        cols <- unique(cols)
        adj <- 1:length(cols)
        names(adj) <- cols
        adj <- as.data.table(t(adj))
        adj[, (cols) := NA]
        
        if ('file' %in% output){
          fst::write_fst(adj,
                         path = normalizePath(paste0(subdirs[3],"/",i,".fst"),mustWork=FALSE))
        }
        if ('object' %in% output){
          return(adj)
        }
        
      })
    }
  }
}