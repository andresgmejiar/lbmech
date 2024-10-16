#' Calculate a multiplicatively-weighed Voronoi diagram, where distances are
#' divided by some per-observation weight. On planar an spherical topologies
#' all segments are circles and segments thereof (or straight lines if observations
#' with equal weights exist), with the radii of the circles corresponding to
#' \eqn{\frac{r_e}{4}\left(\frac{w_p}{w_p + w_q}\right)}, with \eqn{r_e} 
#' corresponding to earth's radius. The circle is positioned such that (1) it 
#' encompasses the point of interest and (2) the point immediately between \eqn{p} 
#' and \eqn{q} lies at \eqn{\ell_{pq} \left(\frac{w_p}{w_p + w_q}\right)} units away 
#' from \eqn{p} with \eqn{\ell_{pq}} being the distance between \eqn{p} and {q}, and 
#' (3) oriented with this point normal to the geodesic between \eqn{p} and {q}.
#' 
#' At present, only calculations on the Earth geoid
#' according to WGS 1984 are supported using the \code{\link[geosphere]{geosphere}}
#' package.
#' 
#' @title Multiplicatively-weighted Voronoi Polygons
#' @param xy Something coercible to a SpatVect points object, with xy units in
#' lonlat and at least two observations. 
#' The point locations of interest. In the future lines and polygons will be 
#' supported
#' @param w A vector of equal length to xy. The weights corresponding to each 
#' point
#' @param tolerance How many digits of the lonlat coordinates are kept (to avoid
#' floating point errors?). Default is \code{tolerance = 7}
#' @param prec How many segments does each ellipsoid contain? Default is 
#' \code{prec = 72}, or one every five degrees  
#' @param clip An object of class null, logical, or coercible to a SpatExtent.
#' Should the output polygons encompass the whole world (the default, 
#' \code{clip = NULL})? If not, \code{clip = TRUE} crops them to the extent of the
#' input points, while passing an object with a SpatExtent will crop them by such
#' an extent. 
#' @param topology One of \code{'geoid'}, \code{'spherical'}, or \code{'planar'},
#' corresponding to the underlying topology. At present, only \code{'geoid'} is 
#' supported assuming the WGS 1984 geoid. 
#' @param pb Logical. Should a progress bar be displayed? 
#' Default is \code{pb = FALSE}
#' @importFrom terra vect
#' @importFrom terra crs
#' @importFrom terra erase
#' @importFrom terra intersect
#' @importFrom terra union
#' @importFrom terra crop
#' @importFrom terra relate
#' @importFrom terra aggregate
#' @importFrom terra disagg
#' @importFrom terra geom
#' @importFrom terra values
#' @importFrom terra ext
#' @importFrom terra makeValid
#' @importFrom data.table data.table
#' @importFrom data.table as.data.table
#' @importFrom data.table melt
#' @importFrom data.table fifelse
#' @importFrom data.table shift
#' @importFrom data.table first
#' @return A SpatVector containing the Voronoi polygons for each input observation
#' @examples 
#' 
#' set.seed(3755)
#' # Create dummy observations
#' obs <- data.table(lon = runif(25, -180, 180),
#'                   lat = runif(25, -90, 90),
#'                   N = runif(25, 0, 100))
#'                   
#' mwv <- mwVoronoi(obs[,1:2], w = obs$N)
#' @export
mwVoronoi <- function(xy, w, tolerance = 7, prec = 72, clip = NULL,
                      topology = 'geoid', pb = FALSE){
  # Silence CRAN warnings
  N=Rank=.I=d=..p=lon=lat=ell=theta=r_e=w_ratio=r=p1_lon=p1_lat=x0=y0=..segs=NULL
  ..seg=V1=V2=dLon=LonTest=.N=SegID=Hemi_Init=.GRP=y=x=NULL
  spat <- c('SpatVector','SpatialPoints','SpatialPointsDataFrame')
  if (any(class(xy) %in% spat)){
    xy <- geom(vect(xy))
  }
  
  # Start by ordering them based on weight
  obs <- as.data.table(xy)
  obs$N <- w
  obs$Orig_Order <- 1:nrow(obs)
  
  obs <- obs[order(N)]
  obs[, Rank := .I]
  orig_order <- obs$Orig_Order
  
  shps <- list()
  taken <- vect()
  iter = 0
  
  if (pb) pb1 <- txtProgressBar(min = 0, max = (nrow(obs)^2 / 2 - 1), style = 3)
  
  for (i in 1:(nrow(obs)-1)){
    # Select lowest-ranking point. We need to find the point immediately between that
    # point and all higher-ranking ones (q) that lies along the Voroni border (Point 1),
    # and the radius of the Voroni sphere
    p <- obs[Rank == i]
    cells <- obs[Rank > i]
    
    # Find the distance between lowest-ranking point and those above it
    cells[,d := geosphere::distGeo(..p,data.table(lon,lat))]
    
    # Distance to Point 1 is the distance between p and q times the weights ratio
    cells[,ell := d * ..p$N / (N + ..p$N)]
    
    # Get the initial bearings between the points
    cells[,theta := geosphere::bearing(..p, data.table(lon,lat))]
    
    # And now find the location of Point 1
    nms <- c('p1_lon','p1_lat')
    cells[, (nms) := as.data.table(geosphere::destPoint(..p, 
                                                        b = theta, 
                                                        d = ell))]
    
    # When weights are equal, radius of pairwise voronoi is equal to one-fourth
    # earth's circumference. Linear decay to zero when one's weight equals zero
    # Start by finding the radius of the ellipsoid along the axis defined by the 
    # Great Circle connecting the two points
    cells[, r_e := GreatCircleCircum(p[1],p[2], theta[1], n = prec) / 4,
          by = 'Rank']
    cells[, w_ratio := ..p$N /(..p$N + N)
    ][, r := r_e * w_ratio]
    
    # Find origin
    cells[, (c('x0','y0')) := as.data.table(
      geosphere::destPoint(data.table(p1_lon, p1_lat),
                b = geosphere::bearing(data.table(p1_lon, p1_lat),
                            ..p),
                d = abs(r))),
      by = 'Rank']
    
    # At what bearings from the Voroni ellipsoid centroid are we going to calculate
    # the location of the border?
    segs <- seq(from = 0, to = 360, length.out = prec)
    segNames <- paste0('Seg',1:(prec))
    
    # For each of those sections, calculate the radius along that bearing, and 
    # find the coordinates of the destination point
    for (seg in 1:length(segs)){
      snames <- paste0(segNames[seg],c('_lon','_lat'))
      cells[, (snames) := as.data.table(geosphere::destPoint(
        data.table(x0,y0),
        b = ..segs[..seg],
        d = w_ratio  * GreatCircleCircum(x0,y0,b = ..segs[..seg], n = prec) / 4
      )), by = 'Rank']
    }
    
    lat_cols = names(cells)[stringr::str_detect(names(cells), 'Seg\\d+_lat')]
    lon_cols = names(cells)[stringr::str_detect(names(cells), 'Seg\\d+_lon')]
    
    # Pivot table to long format
    cells <- melt(cells[,-1], measure.vars = lon_cols, variable.name = 'V1',
                  value.name = 'lon', variable.factor = FALSE
    )[, V1 := stringr::str_remove(V1,'_\\w+')][]
    cells <- melt(cells[,-1], measure.vars = lat_cols, variable.name = 'V2',
                  value.name = 'lat', variable.factor = FALSE
    )[, V2 := stringr::str_remove(V2,'_\\w+')
    ][V1 == V2, -1][, `:=`(SegID = paste0(..p$Rank,'-',Rank,'_',V1),
                           V1 = NULL,
                           V2 = NULL)][order(Rank)]
    
    # If the longitude changes signs and the difference is more than 90, then it
    # crossed the international date line
    cells[, dLon := shift(lon) - lon, by = 'Rank']
    
    # Find the shapes that have points meeting the above conditions
    lonchange <- cells[abs(dLon) > 90
    ][, LonTest := unlist(fifelse(.N != 0, list(1:.N), list(integer(0)))),
      by = 'Rank'
    ]
    
    # If it changes signs, to ensure that the shape is drawn as a single part
    # and contiguity is followed add (or subtract, depending on which hemisphere
    # bearing of 0 gets you) 360 to all values after the sign change, unless
    # the sign changes again
    lplus <- lonchange[LonTest == 1]$SegID
    lminus <- lonchange[LonTest == 2]$SegID
    
    # What's the initial hemisphere
    cells[SegID %in% lplus, Hemi_Init := lon/abs(lon)
    ][, Hemi_Init := unique(stats::na.omit(Hemi_Init)),by = 'Rank'
    ][is.na(Hemi_Init), Hemi_Init := first(lon)/abs(first(lon)),by = 'Rank']
    
    # Switch those signs
    cells[, dLon := fifelse(SegID %in% lplus, 
                            TRUE, 
                            NA)][, dLon := zoo::na.locf(dLon, 
                                                        na.rm = FALSE),
                                 by = 'Rank'
                            ][dLon == TRUE ,
                              lon := lon - 360 * Hemi_Init 
                            ][,dLon := NULL
                            ]
    
    # And switch back those jf there were two sign changes
    cells[, dLon := fifelse(SegID %in% lminus, 
                            TRUE, 
                            NA)][, dLon := zoo::na.locf(dLon, 
                                                        na.rm = FALSE),
                                 by = 'Rank'
                            ][dLon == TRUE,
                              lon := lon + 360 * Hemi_Init
                            ][,dLon := NULL
                            ]
    
    # Convert back to SpatVector
    cells_vect <- vect(as.matrix(cells[,.(.GRP,.GRP,lon,lat,FALSE), 
                                       by = 'Rank'][,-1]),
                       type = 'polygons',
                       crs = '+proj=lonlat')
    cells_vect$Index <- 1:nrow(cells_vect)
    p_vect <- vect(p, crs = '+proj=lonlat')
    cells_vect <- makeValid(cells_vect)
    
    
    # If the circle encompasses the north or south poles, Gaussian-like
    # shapes are produced. These need to be extended to encompass 
    # their respective poles
    
    # To identify the Gaussian-like shapes versus closed ellipsoids, 
    # find those whose maximum or minimum Y coordinates are dispersed about a 
    # range of 360 degrees (closed ellipsoids will have a range close to zero).
    tonorth <- as.data.table(geom(cells_vect)
    )[order(-y)
    ][,max(x[1:6]
    ) -min(x[1:6]), by = c('geom','part')
    ][V1 == 360]$geom
    
    tosouth <- as.data.table(geom(cells_vect)
    )[order(y)
    ][,max(x[1:6]
    ) -min(x[1:6]), by = c('geom','part')
    ][V1 == 360]$geom
    
    # Before we actualyl add the extra latitudes we need to make sure that
    # everything is within a range of -180 to 180 longitude. Extract
    # parts that need to be moved left
    cells_vect_move_left <- crop(cells_vect, ext(c(180,180+360,-90,90)))
    cells_attr_left <- values(cells_vect_move_left)
    
    # And those that need to be moved right
    cells_vect_move_right <- crop(cells_vect, ext(c(-180-360,-180,-90,90)))
    cells_attr_right <- values(cells_vect_move_right)
    
    # And those that are already OK
    cells_vect_keep <- crop(cells_vect, ext(c(-180,180,-90,90)))
    crs(cells_vect_keep) <- '+proj=lonlat'
    
    # Move left
    cells_vect_move_left <- vect(as.matrix(as.data.table(geom(
      cells_vect_move_left))[
        , x := x - 360
      ][]), type = 'polygon',
      crs = '+proj=lonlat')
    values(cells_vect_move_left) <- cells_attr_left
    
    # Move right
    cells_vect_move_right <- vect(as.matrix(as.data.table(geom(
      cells_vect_move_right))[
        , x := x + 360
      ][]), type = 'polygon',
      crs = '+proj=lonlat')
    values(cells_vect_move_right) <- cells_attr_right
    
    # And combine
    cells_vect <- aggregate(rbind(
      round(cells_vect_move_left, digits = tolerance),
      round(cells_vect_move_right, digits = tolerance),
      round(cells_vect_keep, digits = tolerance)),
      by = 'Index')
    
    # To add the south and north bits, add a rectangle ranging to abs(180) and the 
    # Gaussian's max (or min) y coordinate
    for (j in tosouth){
      wrld <- vect(ext(c(ext(cells_vect[j])[1],
                         ext(cells_vect[j])[2],
                         -90, ext(cells_vect[j])[3])), crs = '+proj=lonlat')
      cells <- aggregate(rbind(wrld,cells_vect[j]))
      values(cells) <- values(cells_vect[j])
      
      cells_vect <- rbind(cells_vect[cells_vect$Index != j],
                          cells)
      cells_vect <- cells_vect[order(cells_vect$Index)]
    }
    
    # Same for the northern hemisphere
    for (j in tonorth){
      wrld <- vect(ext(c(ext(cells_vect[j])[1],
                         ext(cells_vect[j])[2],
                         ext(cells_vect[j])[4], 90)), crs = '+proj=lonlat')
      cells <- aggregate(rbind(wrld,cells_vect[j]))
      values(cells) <- values(cells_vect[j])
      
      cells_vect <- rbind(cells_vect[cells_vect$Index != j],
                          cells)
      cells_vect <- cells_vect[order(cells_vect$Index)]
    }
    
    # See whether the selected point is inside the generated shapes
    inside <- relate(cells_vect, p_vect, relation = 'intersects')
    
    for (j in 1:(nrow(cells_vect))){
      if (!inside[j]){
        # If it's not, crop a world-sized rectangle to the shapes...
        cell <- cells_vect[cells_vect$Index == j]
        wrld <- vect(ext(c(min(c(-180, ext(cell)[1])),
                           max(c(180, ext(cell)[2])),
                           -90,
                           90
        )))
        
        # Find those that are rectangles-to-hexagons
        new_cell <- disagg(union(wrld, cell))
        tetra <- as.data.table(geom(new_cell))[,.N, by = 'geom'][N < 10]$geom
        tetra <- unique(tetra,
                        as.data.table(geom(new_cell)
                        )[,length(unique(y)),by='geom'
                        ][V1 < 10]$geom)
        
        # Join that to the original
        cell <- aggregate(rbind(cell,new_cell[tetra]))
        
        # And add it to the output shapes
        cell <- rbind(new_cell)
        cell$Index <- j
        cells_vect <- rbind(cells_vect[cells_vect$Index != j],
                            cell)
      }
    }
    
    # Re-add projection information
    crs(cells_vect) <- '+proj=lonlat'
    
    # Make sure we only keep shapes that encompass the selected point
    ins <- relate(cells_vect, p_vect, relation = 'intersects')
    cells_vect <- cells_vect[ins]
    
    # Re-order them
    cells_vect <- cells_vect[order(cells_vect$Index)]
    
    # The final voroni polygon is going to be the intersection of all polygons.
    # Instantiate with one that encompasses the entire world
    shp <- vect(ext(-180,180,-90,90),crs = '+proj=lonlat')
    for (j in 1:nrow(cells_vect)){
      shp <- intersect(shp,cells_vect[j])
    }
    
    # Some of the identified space may already be taken by lower-ranked observations.
    # Erase that
    if (i != 1) {
      shp <- erase(shp, taken)
    }
    # And add this observation to the taken polygon
    taken <- aggregate(rbind(taken,shp))
    
    
    
    # Add to list
    values(shp) <- data.table(ID = i)
    shps[[i]] <- shp
    
    iter = iter + nrow(obs) - i + 1
    
    if (pb) setTxtProgressBar(pb1, iter)
  }  
  
  
  # The highest-ranked observation has everything that hasn't been taken
  wrld <- vect(ext(-180,180,-90,90), crs = '+proj=lonlat')
  largest <- erase(wrld,taken)
  values(largest) <- data.table(ID = i + 1)
  shp <- shps
  shp[[i + 1]] = largest
  
  # Combine list into final polygon
  shp <- do.call(rbind, shp)
  
  obs_vect <- vect(obs, crs = '+proj=lonlat')
  
  # Crop extent if needed
  if (methods::is(clip, 'logical')){
    if (clip) shp <- crop(shp, ext(obs_vect))
  } else if (stringr::str_detect(class(clip),'^Spat')){
    shp <- terra::crop(shp, clip)
  }
  
  shp$Orig_Order <- orig_order
  shp <- shp[order(shp$Orig_Order)]
  values(shp) <- data.table(Order = 1:nrow(shp))
  
  return(shp)
}