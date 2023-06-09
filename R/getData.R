#' Access various datasets used in vignettes and studies associated with the broader \code{lbmech} project
#' 
#' \code{data = 'baleares-currents'} imports ocean current data around Mallorca, Menorca, and Cabrera on June 13, 2022. Data originally downloaded from 
#' \href{https://data.marine.copernicus.eu/}{E.U. Copernicus Marine Service Information}. See Clementi et al. (2021). 
#' An internet connection is \strong{not} needed, and parameters \code{name} and \code{dir} are ignored.
#'
#' \code{data = 'baleares-gps'} downloads GPS tracks for human hikes in the Balearic Islands
#' in GPX format from \url{https://osf.io/77n9t/}. See Lera et al. (2017).
#' The \code{name} parameter will define the folder name in the \code{dir} directory to which the \code{.gpx} files are saved. Default \code{name = 'gpx'}.
#' \code{getData} will \strong{not} import these tracks, for that, use \code{\link[lbmech]{importGPX}}.  An internet connection is needed.
#'  
#' \code{data = 'baleares-places'} imports a SpatVector with twelve locations on Mallorca, Menorca, and Cabrera in the Balearic Islands. 
#'  
#' @title Access specific datasets
#' @param data A character string indicating the name of the dataset to access. See details below.
#' @param name For values of \code{data} that download data (see details), what will be the file/directory name of the
#' downloaded items?
#' @param dir Directory to which downloaded data will be saved. Default is \code{tempdior()}.
#' @param timeout How many seconds before downloads time out? Default is 999. Temporarily overrides value in \code{getOptions("timeout")}.
#' @importFrom terra vect
#' @return Various, depending on \code{data} selection:
#' 
#' \code{data = 'baleares-currents'} returns a SpatRaster with ocean current surface velocities in m/s.
#' 
#' \code{data = 'baleares-gpx'} does not return any object, but creates a sub-directory \code{name} in directory \code{dir} with 15,373 GPX files,
#' of which 15,371 can be successfully imported using importGPX. 
#' 
#' \code{data = 'baleares-places'} imports a SpatVector with twelve points. 
#' @references 
#' Clementi, E., Aydogdu A., Goglio, A. C., Pistoia J., Escudier R., Drudi M., Grandi A., et al. (2021). 
#' Mediterranean Sea Physics Analysis and Forecast (CMEMS MED-Currents, EAS6 System). 
#' \emph{Copernicus Marine Service}. \url{https://doi.org/10.25423/CMCC/MEDSEA_ANALYSISFORECAST_PHY_006_013_EAS7}.
#' 
#' Lera I., Perez T., Guerrero C., Eguiluz V. M., Juiz C. (2017). 
#' Analysing human mobility patterns of hiking activities through complex network theory. 
#' \emph{PLoS ONE} 12(5): e0177712. \url{https://doi.org/10.1371/journal.pone.0177712}
#' @examples 
#' # Import ocean current data for the Balearic Islands
#' 
#' currents <- getData('baleares-currents')
#' @export
getData <- function(data, name = NULL, dir = tempdir(), timeout = 999){
  
  dir <- normalizePath(dir,mustWork=FALSE)
  
  if (data == 'baleares-currents'){
    message("Clementi, Emanuela; Ali Aydogdu; Ana Chiara Goglio; Jenny Pistoia; Romain Escudier; Massimiliano Drudi; Alessandro Grandi; et al.\n(2021) Ocean Surface Currents around Mallorca, Menorca, and Cabrera, 13 June 2022. Mediterranean Sea Physics Analysis and Forecast, CMEMS MED-Currents EAS6 system. European Union Copernicus Marine Service, https://marine.copernicus.eu/\n")

    return(importRST(baleares.currents))
  } else if (data == 'baleares-gps'){
    
    if (is.null(name)) name <- 'gpx'
    
    message("Lera, Isaac; Toni Perez; Carlos Guerrero; Victor M. Eguiluz; and Carlos Juiz;\n(2017) Dataset of Human Hiking GPS Trajectories on Balearic Islands. Supplemental material to \"Analysing human mobility patterns of hiking activities through complex network theory,\" PLoS ONE 12(5): e0177712. https://doi.org/10.1371/journal.pone.0177712\n")
    message("Dataset will be download and uncompressed but will not be imported.\nPlease use 'importGPX()' for that.") 

    timeout.orig <- getOption('timeout')
    options(timeout=timeout)
    gpxDir <- normalizePath(paste0(dir, "/",name),mustWork = FALSE)
    message(paste0("Output Directory: ",gpxDir))
    
    if (!dir.exists(gpxDir)){
      filePath <- normalizePath(paste0(tempdir(),"/GPXtmp8",
                                       paste(letters[round(stats::runif(9,1,26))],collapse=""),'.tar.gz'),
                                mustWork = FALSE)
      utils::download.file("https://osf.io/download/cjv8u/",
                           destfile = filePath,
                           mode = 'wb')
      archive::archive_extract(filePath, dir=gpxDir)
      unlink(filePath)
    } else {
      warning(paste(gpxDir,"already exists. No action performed."))
    }
    options(timeout=timeout.orig)
  } else if (data == 'baleares-places'){
    return(vect(baleares.places, geom=c('long','lat'), crs = '+proj=longlat')
)
  } else {
    stop("Unknown dataset. Run '?lbmech::getData' to see available datasets.")
  }
}
