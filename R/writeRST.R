#' Read and write rasters using the fst library
#' 
#' @title Fast read and write rasters
#' @name writeRST
#' @aliases importRST
#' @param x An object of class SpatRaster or Raster. It may not contain
#' layers named 'x' or 'y'
#' @param filename Character. Output filename. Do not use extensions. If null,
#' (the defaul), no file is saved
#' @param object Logical. If true, function returns a data.table passable to 
#' \code{\link[lbmech]{importRST}} that will produce a valid raster (i.e. fixes
#' rounding errors when using \code{\link[terra]{init}}. Default is \code{FALSE}.
#' @param layers Character vector containing the names of the layers to import.
#' Default is \code{layers = NULL} which imports all layers. 
#' @param ... Additional parameters to pass on to \code{\link[fst]{read_fst}} or
#' \code{\link[fst]{write_fst}}
#' @return If \code{object = TRUE}, a data.table with an added \code{$crs} attribute
#' containing the projection information. Special values are stored in the first 
#' three rows:
#' (1) The x and y resolution
#' (2) The x and y coordinates of the corner, and
#' (3) The x and y coordinates of the rounding offset
#' @importFrom terra crs
#' @importFrom terra crs<-
#' @examples  
#' n <- 5
#' dem <- expand.grid(list(x = 1:(n * 100),
#'                         y = 1:(n * 100))) / 100
#' dem <- as.data.table(dem)
#' dem[, z := 250 * exp(-(x - n/2)^2) + 
#'       250 * exp(-(y - n/2)^2)]
#' dem <- rast(dem)
#' ext(dem) <- c(10000, 20000, 30000, 40000)
#' crs(dem) <- "+proj=lcc +lat_1=48 +lat_2=33 +lon_0=-100 +datum=WGS84"
#' writeRST(dem, 'DEM.fst')
#' importRST('DEM.fst')
NULL

#' @rdname writeRST
#' @export
writeRST <- function(x, filename = NULL, object = FALSE, ...){
  
  if ('x' %in% names(x) | 'y' %in% names(x)){
    stop("Raster may not contain layers named 'x' or 'y'. Operation aborted")
  }
  xres <- res(x)[1]
  yres <- res(x)[2]
  xoff <- round(as.numeric(xFromCell(x,1))/xres) * 
    xres - as.numeric(xFromCell(x,1))
  yoff <- round(as.numeric(yFromCell(x,1))/yres) * 
    yres - as.numeric(yFromCell(x,1))
  xcorner <- xFromCell(x,1)
  ycorner <- yFromCell(x,1)
  
  out <- rbind(data.table(x = c(xres,xoff,xcorner),
                          y = c(yres,yoff,ycorner)),
               rastToTable(x), fill = TRUE)
  
  if (!is.null(filename)){
    filename <- normalizePath(filename, mustWork = FALSE)
    fst::write_fst(out, 
                   paste0(filename,'.fst'), compress = 0, ...)
    write(crs(x), paste0(filename,'.fstproj'))
  }
  if (object){
    attributes(out)$crs <- crs(x)
    return(out)
  }
}

#' @rdname writeRST
#' @export
importRST <- function(x, layers = NULL,...){
  ..xres=y=..yres=NULL
  
  if (!is.null(layers)){
    layers = c('x','y',layers)
  }
  
  if (!methods::is(x,'data.table')){
    filename <- x
    filename <- normalizePath(filename, mustWork = FALSE)
    proj <-   as.character(noquote(paste(readLines(paste0(filename,'.fstproj')),
                                         collapse='\n')))
    x <- fst::read_fst(normalizePath(paste0(filename,'.fst')),
                       columns = layers,
                       as.data.table = TRUE)
  } else {
    proj <- attributes(x)$crs
  }
  
  xres <- x[1]$x
  yres <- x[1]$y
  xoff <- x[2]$x
  yoff <- x[2]$y
  xcorner <- x[3]$x
  ycorner <- x[3]$y
  x[, `:=`(x = round(x/..xres)*..xres + xoff,
           y = round(y/..yres)*..yres + yoff)]
  x <- tryCatch(rast(x[!c(1,2,3)],
                     crs = proj,
  ), error = function(e) rast(xmin = xcorner, xmax = xcorner + xres,
                              ymin = ycorner, ymax = ycorner + yres,
                              res = c(xres,yres), names = 'z', vals = NA,
                              crs = proj)
  )
  return(x)
}
