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
#' @title Multiplicatively-weighted Voronoi Polygons
#' @param xy Something coercible to a SpatVect points object and with at least 
#' two observations. The point locations of interest. Units must be degrees longitude for 
#' \code{x} and degrees latitude for \code{y} when \code{topology} is either 
#' \code{'geoid'} or \code{'sphere'}.  In the future lines and polygons will be 
#' supported
#' @param w A vector of equal length to xy. The weights corresponding to each 
#' point
#' @param x A character vector indicating the column with 'x' coordinates. Ignored
#' if \code{xy} is a SpatVect points object.
#' @param y A character vector indicating the column with 'y' coordinates. Ignored
#' if \code{xy} is a SpatVect points object.
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
#' corresponding to the underlying topology. Ignored if \code{topology} is
#' not \code{'geoid'}.
#' @param a Equatorial radius. Default is for WGS84. Ignored if \code{topology} is
#' not \code{'geoid'}.
#' @param f Ellipsoidal flattening. Default is for WGS84
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
                      x = 'lon', y = 'lat',
                      topology = 'geoid', a = 6378137, 
                      f = 1/298.257223563, pb = FALSE){
  # Silence CRAN warnings
  N=Rank=.I=d=..p=lon=lat=ell=theta=r_e=w_ratio=r=p1_lon=p1_lat=x0=y0=..segs=NULL
  ..seg=V1=V2=dLon=LonTest=.N=SegID=Hemi_Init=.GRP=NULL
  ..x=..y=ell1=ell2=x_c=y_c=b=m=finite_m=y_left=y_right=x_bottom=x_top=points=NULL
  valid_x_vert=x_vert=plon=plat=..a=..f=pN=NULL
  spat <- c('SpatVector','SpatialPoints','SpatialPointsDataFrame')
  if (any(class(xy) %in% spat)){
    if (!methods::is(xy,'SpatVector')) xy <- vect(xy)
    proj <- crs(xy)
    xy <- geom(xy)
    x <- colnames(xy)[3]
    y <- colnames(xy)[4]
  } else {
    proj <- ''
  }

  if (topology == 'sphere'){
    f <- 0
    topology <- 'geoid'
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

  if (topology == 'planar'){
    obs_vect <- vect(obs[,.(x = get(..x),y = get(..y))],
                     geom = c('x','y'),
                     crs = proj)
    if (stringr::str_detect(class(clip),'^Spat')){
      wrld <- vect(ext(clip), crs = proj)
    } else {
      wrld <- vect(ext(obs_vect), crs = proj)
    }
    for (i in 1:(nrow(obs)-1)){
      # Select lowest-ranking point. We need to find the point immediately between that
      # point and all higher-ranking ones (q) that lies along the Voronoi border (Point 1),
      # and the radius of the Voronoi circle
      p <- obs[Rank == i, .(x = get(..x), y = get(..y), N)]
      p_vect <- vect(p[,.(x, y)], geom = c('x','y'), crs = proj)
      
      # Make sure that we exclude any lower-ranked points with the same weight, since 
      # we'd get a divide by zero problem. We'll deal with these ties later
      cells <- obs[Rank > i & N != p$N, .(x = get(..x), y = get(..y), N)]  
      
      # Find distance between lowest-ranking point and those above it
      cells[, d := sqrt((..p$x - x)^2 + (..p$y - y)^2)]
      
      # Distance to Point 1 is the distance between p and q times the weights ratio
      cells[, ell1 := d * ..p$N / (N + ..p$N)]
      
      # Distance to Point 2 is the ratio of the other point to the weights difference
      cells[, ell2 := d * N / (N - ..p$N)]
      
      # Radius is half the sum of distances
      cells[, r := (ell1 + ell2) / 2]
      
      # Center is along the same line but opposite bearing from q -> p
      cells[, theta := atan2(y - ..p$y, x - ..p$x) - 3.14159265359]
      
      # Center will be at above bearing, r - ell1 distance away from p
      cells[, `:=`(x_c = ..p$x + (r - ell1) * cos(theta),
                   y_c = ..p$y + (r - ell1) * sin(theta))]
      
      # Generate circle polygons
      cells_vect <- vect(cells[,.(x = x_c, y = y_c)],
                         crs = '',
                         geom = c('x','y'))
      suppressWarnings(cells_vect <- buffer(cells_vect,
                                            width = cells$r,
                                            quadsegs = prec))
      
      # Now deal with ties. 
      ties <- obs[Rank > i & N == p$N, .(x,y,Rank)]
      
      if (nrow(ties) > 0){
        
        # Get the slope m and intercept b of the perpendicular bisector
        ties[, `:=`(x_c = (x + p$x)/2,
                    y_c = (y + p$y)/2,
                    m = (x - p$x)/(p$y - y))
        ][, b := y_c - m * x_c]
        
        # Get extent of the world
        xmin <- ext(wrld)[1]
        xmax <- ext(wrld)[2]
        ymin <- ext(wrld)[3]
        ymax <- ext(wrld)[4]
        
        # Handle infinite slopes separately
        ties[, finite_m := is.finite(m)]
        
        # Compute intersection points for finite slopes
        ties[finite_m == TRUE, `:=`(
          y_left = m * xmin + b,
          y_right = m * xmax + b,
          x_bottom = (ymin - b) / m,
          x_top = (ymax - b) / m
        )]
        
        # Check if the intersection points are within 'wrld' boundaries
        ties[finite_m == TRUE, `:=`(
          valid_y_left = y_left >= ymin & y_left <= ymax,
          valid_y_right = y_right >= ymin & y_right <= ymax,
          valid_x_bottom = x_bottom >= xmin & x_bottom <= xmax,
          valid_x_top = x_top >= xmin & x_top <= xmax
        )]
        
        # Collect intersection points for finite slopes
        ties[finite_m == TRUE, points := Map(function(vl, xl, vy, xr, vb, xb, vt, xt) {
          pts <- list()
          if (vl) pts <- append(pts, list(c(xmin, xl)))
          if (vy) pts <- append(pts, list(c(xmax, xr)))
          if (vb) pts <- append(pts, list(c(xb, ymin)))
          if (vt) pts <- append(pts, list(c(xt, ymax)))
          if (length(pts) >= 2) {
            pts <- unique(pts)
            do.call(rbind, pts)
          } else {
            NULL
          }
        }, ties$valid_y_left, ties$y_left, ties$valid_y_right, ties$y_right,
        ties$valid_x_bottom, ties$x_bottom, ties$valid_x_top, ties$x_top)]
        
        # Handle vertical lines (infinite slopes)
        ties[finite_m == FALSE, `:=`(
          x_vert = x,
          valid_x_vert = x >= xmin & x <= xmax
        )]
        ties[finite_m == FALSE, points := ifelse(valid_x_vert,
                                                 list(rbind(c(x_vert, ymin), 
                                                            c(x_vert, ymax))), 
                                                 list(NULL))]
        
        # Remove entries without valid intersection points
        ties <- ties[!sapply(points, is.null)]
        
        # Create line geometries from intersection points
        ties <- vect(ties$points, type = "lines")
        
        # Split 'wrld' by all lines at once
        ties <- lapply(1:length(ties), function(x) split(wrld, ties[x]))
        ties <- do.call(rbind,ties)
        ins <- relate(ties, p_vect, 'covers')
        ties <- ties[ins]
        
        cells_vect <- rbind(cells_vect, ties)
      }
      # The final voronoi polygon is going to be the intersection of all polygons.
      # Instantiate with one that encompasses the entire world
      shp <- wrld
      
      for (j in 1:nrow(cells_vect)){
        shp <- intersect(shp,cells_vect[j])
      }
      
      # Some of the identified space may already be taken by lower-ranked observations.
      # Erase that
      if (i != 1) {
        shp <- erase(shp, taken)
      }
      
      # And add this observation to the taken polygon
      taken <- aggregate(rbind(shp,taken))
      
      # Add to list
      values(shp) <- data.table(ID = i)
      shps[[i]] <- shp
      
      iter = iter + nrow(obs) - i + 1
      
      if (pb) setTxtProgressBar(pb1, iter)
    }
  } else {

    # Find all unique combinations of points, sorted such that the one with the 
    # lower weight goes first. Ignore ties
    cells <- CJ(p = obs$Rank, q = obs$Rank
                )[p <q
                  ][obs[,.(p = Rank, plon = get(..x), plat = get(..y), pN = N)], 
                    on = 'p'
                    ][obs[,.(q = Rank, lon = get(..x), lat = get(..y), N)], 
                      on = 'q'
                      ][order(q)
                        ][order(p)]
    cells <- stats::na.omit(cells)

    # Calculate the distance between the points
    cells[, d := geosphere::distGeo(data.table(lon = plon, lat = plat),
                                    data.table(lon, lat),
                                    a = ..a,
                                    f = ..f)]
    
    # Distance to Point 1 is the distance between p and q times the weights ratio
    cells[, ell := d * pN / (N + pN)]

    # Get the initial bearings between the points
    cells[,theta := geosphere::bearing(data.table(lon = plon, lat = plat), 
                                       data.table(lon,lat),
                                       a = ..a,
                                       f = ..f)]

    # And now find the location of Point 1
    nms <- c('p1_lon','p1_lat')
    cells[, (nms) := as.data.table(geosphere::destPoint(data.table(lon = plon,
                                                                   lat = plat), 
                                                        b = theta, 
                                                        d = ell,
                                                        a = ..a,
                                                        f = ..f))]

    # When weights are equal, radius of pairwise voronoi is equal to one-fourth
    # earth's circumference. Linear decay to zero when one's weight equals zero
    # Start by finding the radius of the ellipsoid along the axis defined by the 
    # Great Circle connecting the two points
    cells[, r_e := GreatCircleCircum(plon,plat, theta, n = prec,
                                     a = ..a, 
                                     f = ..f) / 4,
          by = c('p','q')]
    cells[, w_ratio := pN /(pN + N)
    ][, r := r_e * w_ratio]

    # Find origin
    cells[, (c('x0','y0')) := as.data.table(
      geosphere::destPoint(data.table(p1_lon, p1_lat),
                           b = geosphere::bearing(data.table(p1_lon, p1_lat),
                                                  data.table(lon = plon,
                                                             lat = plat), 
                                                  a = ..a, 
                                                  f = ..f),
                           d = abs(r),
                           a = ..a,
                           f = ..f))]
    
    # At what bearings from the Voronoi ellipsoid centroid are we going to calculate
    # the location of the border?
    segs <- seq(from = 0, to = 360, length.out = prec)
    segNames <- paste0('Seg',1:(prec))

    cells <- cells[rep(1:nrow(cells), each = prec)]
    cells$SegID <- paste0(cells$p,'-',cells$q,'_Seg',rep(1:prec,nrow(cells)/length(segs)))
    cells[, b := rep(segs, nrow(cells)/length(segs))
    ][, (c('lon','lat')) := as.data.table(geosphere::destPoint(
      data.table(x0,y0),
      b = b,
      d = w_ratio  * mapply(GreatCircleCircum, x0, y0 , b = b, n = prec,
                                       a = ..a, 
                            f = ..f) / 4,
      a = ..a,
      f = ..f
    ))
    ]
    
    # If the longitude changes signs and the difference is more than 90, then it
    # crossed the international date line
    cells[, dLon := shift(lon) - lon, by = c('p','q')]
    
    # Find the shapes that have points meeting the above conditions
    lonchange <- cells[abs(dLon) > 90
    ][, LonTest := unlist(fifelse(.N != 0, list(1:.N), list(integer(0)))),
      by = c('p','q')
    ]
    
    # If it changes signs, to ensure that the shape is drawn as a single part
    # and contiguity is followed add (or subtract, depending on which hemisphere
    # bearing of 0 gets you) 360 to all values after the sign change, unless
    # the sign changes again
    lplus <- lonchange[LonTest == 1]$SegID
    lminus <- lonchange[LonTest == 2]$SegID
    
    # What's the initial hemisphere
    cells[SegID %in% lplus, Hemi_Init := lon/abs(lon)
    ][, Hemi_Init := unique(stats::na.omit(Hemi_Init)),by = c('p','q')
    ][is.na(Hemi_Init), Hemi_Init := first(lon)/abs(first(lon)),by = c('p','q')]
    
    # Switch those signs
    cells[, dLon := fifelse(SegID %in% lplus, 
                            TRUE, 
                            NA)][, dLon := zoo::na.locf(dLon, 
                                                        na.rm = FALSE),
                                 by = c('p','q')
                            ][dLon == TRUE ,
                              lon := lon - 360 * Hemi_Init 
                            ][,dLon := NULL
                            ]
    
    # And switch back those if there were two sign changes
    cells[, dLon := fifelse(SegID %in% lminus, 
                            TRUE, 
                            NA)][, dLon := zoo::na.locf(dLon, 
                                                        na.rm = FALSE),
                                 by = c('p','q')
                            ][dLon == TRUE,
                              lon := lon + 360 * Hemi_Init
                            ][,dLon := NULL
                            ]
    cells_dt <- copy(cells)
    rm(cells)
    for (i in 1:(nrow(obs)-1)){
      # Convert back to SpatVector
      p <- obs[Rank == i, .(lon = get(..x), lat = get(..y), N)]
      cells_vect <- vect(as.matrix(cells_dt[p == i,.(.GRP,.GRP,lon,lat,FALSE), 
                                         by = 'q'][,-1]),
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
      
      # Before we actually add the extra latitudes we need to make sure that
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
    obs_vect <- vect(obs, crs = '+proj=lonlat', geom = c(x,y))
    wrld <- vect(ext(-180,180,-90,90), crs = '+proj=lonlat') 
  }
  if (pb) setTxtProgressBar(pb1, nrow(obs)^2 / 2 - 1)
  # The highest-ranked observation has everything that hasn't been taken
  largest <- erase(wrld,taken)
  values(largest) <- data.table(ID = i + 1)
  shp <- shps
  shp[[i + 1]] = largest
  
  # Combine list into final polygon
  shp <- do.call(rbind, shp)
  
  # Crop maximum extent if needed
  if (methods::is(clip, 'logical')){
    if (clip) shp <- crop(shp, ext(obs_vect))
  } else if (stringr::str_detect(class(clip),'^Spat')){
    shp <- terra::crop(shp, clip)
  }
  
  shp$Orig_Order <- orig_order
  shp <- shp[order(shp$Orig_Order)]
  values(shp) <- data.table(Order = 1:nrow(shp))
  
  crs(shp) <- proj
  return(shp)
}

