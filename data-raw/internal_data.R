# Internal data prep

### Colors for the LID plots
lid.cols <-  data.table::data.table(
  Class = c("Local Average-Global High",
            "Local High-Global High",
            "Local High-Global Average",
            "Local High-Global Low",
            "Local Average-Global Low",
            "Local Low-Global Low",
            "Local Low-Global Average",
            "Local Low-Global High",
            "Local Average-Global Average"),
  Color = c("#FF00FF", "#FF0000", "#FF9F00", 
                     "#FFFD9C", "#ade567", "#63FFD5", 
                     "#0080FF", "#4B0076", "#FFFFFF"),
                     x = c(0, 1, 1, 1, 0, -1, -1, -1, 0),
  y = c(1, 1, 0, -1, -1, -1, 0, 1, 0)
)


### Stuff associated with the Baleares Vignettes

## Locations on Mallorca, Menorca and Cabrera
baleares.places <- data.table::fread("data-raw/baleares_places.csv")

## MOTU Command to download Barealic Island ocean velocity
## Note that you have to provide your own username and password
# rd <- normalizePath('data-raw/',mustWork=FALSE)
# command <- paste0('python -m motuclient --motu https://nrt.cmems-du.eu/motu-web/Motu --service-id MEDSEA_ANALYSISFORECAST_PHY_006_013-TDS --product-id cmems_mod_med_phy-cur_anfc_4.2km_P1D-m --longitude-min 1.8343864083631984 --longitude-max 4.736546729949184 --latitude-min 38.85730429069023 --latitude-max 40.484272955821766 --date-min "2022-06-13 00:00:00" --date-max "2022-06-14 00:00:00" --depth-min 1.0182366371154785 --depth-max 1.0182366371154785 --variable uo --variable vo --out-dir ', rd ,' --out-name Ocean_Currents.nc --user #YOUR#USERNAME#HERE# --pwd #YOUR#PASSWORD#HERE#>')

# system(command, intern = TRUE)

baleares.currents <- terra::rast("data-raw/Ocean_Currents.nc")
baleares.currents <- writeRST(baleares.currents,object=TRUE)



### Export objects

usethis::use_data(lid.cols, baleares.places, baleares.currents,
                  overwrite = TRUE, internal = TRUE)
