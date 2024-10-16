# Internal data prep

### Colors for the LID plots
lid.cols <-  data.table::data.table(
  Class = c("In-group Average, Out-group High",
            "In-group High, Out-group High",
            "In-group High, Out-group Average",
            "In-group High, Out-group Low",
            "In-group Average, Out-group Low",
            "In-group Low, Out-group Low",
            "In-group Low, Out-group Average",
            "In-group Low, Out-group High",
            "In-group Average, Out-group Average"),
  Color = c("#FF00FF", "#FF0000", "#FF9F00", 
                     "#FFFD9C", "#ade567", "#63FFD5", 
                     "#0080FF", "#4B0076", "#FFFFFF"),
  x = c(0, 1, 1, 1, 0, -1, -1, -1, 0),
  y = c(1, 1, 0, -1, -1, -1, 0, 1, 0),
  Arrow = c(sprintf("\u2191"), sprintf("\u2197"), sprintf("\u2192"),
            sprintf("\u2198"), sprintf("\u2193"), sprintf("\u2199"),
            sprintf("\u2190"), sprintf("\u2196"), sprintf("\u2022"))
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
baleares.currents <- terra::shift(baleares.currents, 
                                  dy = -res(baleares.currents)[2]*1.25)
baleares.currents <- writeRST(baleares.currents,object=TRUE)

### Export objects
usethis::use_data(lid.cols, baleares.places, baleares.currents,
                  overwrite = TRUE, internal = TRUE)
