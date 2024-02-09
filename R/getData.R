#' Access various datasets used in vignettes and studies associated with the broader \code{lbmech} project
#' 
#' \code{data = 'baleares-currents'} imports ocean current data around Mallorca, Menorca, and Cabrera on June 13, 2022. Data originally downloaded from 
#' \href{https://data.marine.copernicus.eu/}{E.U. Copernicus Marine Service Information}. See Clementi et al. (2021). 
#' An internet connection is \strong{not} needed, and parameters \code{name} and \code{dir} are ignored.
#'
#' \code{data = 'baleares-gps'} downloads GPS tracks for human hikes in the Balearic Islands
#' in GPX format from \url{https://osf.io/77n9t/}. See Lera et al. (2017).
#' The \code{name} parameter will define the folder name in the \code{dir} directory to which the \code{.gpx} files are saved. Default \code{name = 'gpx'}.
#' \code{getData} will \strong{not} import these tracks, for that, use \code{\link[lbmech]{importGPX}}.  An internet connection is required
#'  
#' \code{data = 'baleares-places'} imports a SpatVector with twelve locations on Mallorca, Menorca, and Cabrera in the Balearic Islands. See Clementi et al. (2021). 
#'  
#' \code{data = 'gdp'} downloads country GDP estimates since 1970 from the World Bank and United Nations. An internet connection is required. 
#'  
#' \code{data = 'population'} downloads country population estimates since 1950 from the United Nations. An internet connection is required.
#' 
#' \code{data = 'trade'} downloads nation-wise bilateral trade data in goods and services. Trade in goods from the OECD, and trade in services from the WTO and OECD. An internet connection is required.
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
#' 
#' \code{data = 'gdp'} imports a data.table with seven columns for Year, Country Code, Country Name, the UN's GDP estimate, the World Bank's GDP estimate, the average GDP estimate, and the UN's GNI estimate.
#' 
#' \code{data = 'population'} imports a data.table with three columns for the Country Code, Year, and Population estimate. 
#' 
#' \code{data = 'trade'} imports a data.table with seven columns. The country codes and names for the nation from which and to which trade is conducted, and the total value of the services and goods exchanged.
#' 
#' @references 
#' Clementi, E., Aydogdu A., Goglio, A. C., Pistoia J., Escudier R., Drudi M., Grandi A., et al. (2021). 
#' Mediterranean Sea Physics Analysis and Forecast (CMEMS MED-Currents, EAS6 System). 
#' \emph{Copernicus Marine Service}. \doi{10.25423/CMCC/MEDSEA_ANALYSISFORECAST_PHY_006_013_EAS7}.
#' 
#' Lera I., Perez T., Guerrero C., Eguiluz V. M., Juiz C. (2017). 
#' Analysing human mobility patterns of hiking activities through complex network theory. 
#' \emph{PLoS ONE} 12(5): e0177712. \doi{10.1371/journal.pone.0177712}
#' 
#' Organization for Economic Co-Operation and Development (2022). 
#' Bilateral Trade Database by Industry and End-Use (BTDIxE). \url{https://www.oecd.org/sti/ind/bilateraltradeingoodsbyindustryandend-usecategory.htm}
#' 
#' United Nations (2022).
#' UNData. \url{https://data.un.org/}
#' 
#' World Bank (2022).
#' World Bank DataBank. \url{https://databank.worldbank.org/}
#' 
#' World Trade Organization, Organization for Economic Co-Operation and Development. (2023).
#' WTO-OECD Balanced Trade in Services Dataset (BaTiS) — BPM6.
#' \url{https://www.wto.org/english/res_e/statis_e/trade_datasets_e.htm#BaTis6}
#' 
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
  } else if (data %in% c('baleares-gps','baleares-gpx')){
    
    if (is.null(name)) name <- 'gpx'
    message("Lera, Isaac; Toni Perez; Carlos Guerrero; Victor M. Eguiluz; and Carlos Juiz \n(2017) Dataset of Human Hiking GPS Trajectories on Balearic Islands. Supplemental material to \"Analysing human mobility patterns of hiking activities through complex network theory,\" PLoS ONE 12(5): e0177712. https://doi.org/10.1371/journal.pone.0177712\n")
    message("Dataset will be downloaded and uncompressed but will not be imported.\nPlease use 'importGPX()' for that.") 

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
    return(vect(baleares.places, geom=c('long','lat'), crs = '+proj=longlat'))
  } else if (toupper(data) == "GDP") {
    message("World Bank \n(2022) World Bank DataBank. https://databank.worldbank.org/\n \nUnited Nations \n(2022) UNData. https://data.un.org/\n")
    
    if (is.null(name)) name <- "GDP.gz"
    if (stringr::str_sub(name,-3) != ".gz") name <- paste0(name,".gz")
    
    filePath <- normalizePath(paste0(dir,"/",name),
                              mustWork = FALSE)
    
    if (!file.exists(filePath)){
    utils::download.file("https://github.com/andresgmejiar/world_trade_data/blob/main/GDP.gz?raw=TRUE",
                         destfile = filePath,
                         mode = 'wb')
    }
    return(readRDS(filePath))
    
  } else if (tolower(data) == "population") {
    message("United Nations; (2022) UNData. https://data.un.org/\n")
    
    if (is.null(name)) name <- "population.gz"
    if (stringr::str_sub(name,-3) != ".gz") name <- paste0(name,".gz")
    
    filePath <- normalizePath(paste0(dir,"/",name),
                              mustWork = FALSE)
    
    if (!file.exists(filePath)){
      utils::download.file("https://github.com/andresgmejiar/world_trade_data/blob/main/Population.gz?raw=TRUE",
                           destfile = filePath,
                           mode = 'wb')
    }
    return(readRDS(filePath))
    
  } else if (tolower(data) == "trade") {
    message("Organization for Economic Co-Operation and Development\n(2022) Bilateral Trade Database by Industry and End-Use (BTDIxE). \nhttps://www.oecd.org/sti/ind/bilateraltradeingoodsbyindustryandend-usecategory.htm\n\nWorld Trade Organization; Organization for Economic Co-Operation and Development \n(2023) WTO-OECD Balanced Trade in Services Dataset (BaTiS) — BPM6. \nhttps://www.wto.org/english/res_e/statis_e/trade_datasets_e.htm#BaTis6\n")
    
    if (is.null(name)) name <- "trade.gz"
    if (stringr::str_sub(name,-3) != ".gz") name <- paste0(name,".gz")
    
    filePath <- normalizePath(paste0(dir,"/",name),
                              mustWork = FALSE)
    
    if (!file.exists(filePath)){
      utils::download.file("https://github.com/andresgmejiar/world_trade_data/blob/main/Trade.gz?raw=TRUE",
                           destfile = filePath,
                           mode = 'wb')
    }
    return(readRDS(filePath))
    
  } else {
    stop("Unknown dataset. Run '?lbmech::getData' to see available datasets.")
  }
}
