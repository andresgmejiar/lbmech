## code to prepare `lid.cols` dataset goes here

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

usethis::use_data(lid.cols, overwrite = TRUE, internal = TRUE)
