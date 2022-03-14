
<!-- README.md is generated from README.Rmd. Please edit that file -->

# lbmech

<!-- badges: start -->
<!-- badges: end -->

`lbmech` is a geospatial package for `R` to calculate time- and
energy-based costs-of-travel for humans and animals moving across the
landscape. While providing a similar functionality to the package
`gdistance`, `lbmech` stores data and performs all linear algebra using
the `data.table` package allowing for in-place modification of objects
greatly increasing processing speed. `lbmech` also modularizes important
aspects of the cost-distance workflow allowing for
computationally-intensive and otherwise prohibitively-large operations
to take place. Moreover, unlike similar tools such as `package` for `R`,
`lbmech` allows for the estimation of various types of energetic losses
(due to kinematic locomotion, work against gravity, basal metabolic
processes) instead of simply the total energetic or metabolic
expenditure.

The example provided in this README contains the examples included in
the individual function documentation generating all data from scratch.
For applied/real-world examples, please see the provided vignettes.

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("andresgmejiar/lbmech")
```

## Example Workflow

``` r
library(lbmech)
#> Loading required package: data.table
#> Warning: package 'data.table' was built under R version 4.1.2
#> Loading required package: raster
#> Loading required package: sp
#> 
#> Attaching package: 'raster'
#> The following object is masked from 'package:data.table':
#> 
#>     shift
# Set random seed for reproducibility's sake
set.seed(5574741)
```

### Part 1: Topographic data sources

The first step in a typical `lbmech` workflow is defining the digital
elevation model (DEM) to define as the topographic data source. This may
be provided in one of two ways:

1.  A RasterLayer object representing the digital elevation model for
    the region-of-interest
    -   In order to ensure that there is a stable file path to the
        source DEM throughout the workflow, the raster must have been
        ‘read in’ using the `raster` function without having been
        further modified. If additional modifications are necessary, use
        the `writeRaster` function to save it to the disk first before
        re-reading it in using `raster`. For example:

``` r
# Generate a DEM
n <- 5
dem <- expand.grid(list(x = 1:(n * 100),
                        y = 1:(n * 100))) / 100
dem <- as.data.table(dem)
dem[, z := 250 * exp(-(x - n/2)^2) + 
      250 * exp(-(y - n/2)^2)]
dem <- rasterFromXYZ(dem)
extent(dem) <- c(10000, 20000, 30000, 40000)
crs(dem) <- "+proj=lcc +lat_1=48 +lat_2=33 +lon_0=-100 +datum=WGS84"

# Export it so it doesn't just exist on the memory
dir <- tempdir()
writeRaster(dem, paste0(dir,"/DEM.tif"),format="GTiff",overwrite=TRUE)


# Import raster
dem <- raster(paste0(dir,"/DEM.tif"))

plot(dem, main = "Digital Elevation Model")
```

<img src="man/figures/README-dem-1.png" width="100%" />

2.  A SpatialPolygonsDataFrame object whose individual polygons
    represent sectors with unique DEM sources stored as a file path or
    URL in the data frame object.
    -   As of the most recent version, `lbmech` supports file types
        readible by `rgdal` and `raster`, as well as such files
        compressed in `gz` and `zip` files–although the latter is likely
        to fail in Unix systems.

Even if you have already downloaded or imported a raster to use as a
topographic data source as in case one above, most of the functions will
expect a SpatialPolygonsDataFrame object in the form of case two. You
can make this using the `makeGrid` function:

``` r
grid <- makeGrid(dem = dem, nx = n, ny = n)

plot(dem, main = "Sectors to divide DEM")
plot(grid,col=NULL,add=TRUE)
#> Warning in rep(col, n, n): 'x' is NULL so the result will be NULL
```

<img src="man/figures/README-makeGrid-1.png" width="100%" />

`lbmech` is specifically designed to deal with **large** regions that
would be prohibitive to analyze if the data is stored exclusively within
the memory. To deal with this issue, `lbmech` will crop any input raster
into an `nx` by `ny` grid and save the sector in its own `gz` file for
case 1. In both cases—to save memory and computational time—sectors are
only cropped or downloaded on an as-needed basis. You can use the
`whichTiles` and `getMap` functions to identify which tile(s) might be
needed, and download or crop any such tiles that haven’t been prepared:

``` r
# Generate five random points that fall within the grid
points <- data.table(x = runif(5, extent(dem)[1], extent(dem)[2]),
                     y = runif(5, extent(dem)[3], extent(dem)[4]))
               
                           
# Run whichTiles and getMap to prepare appropriate sector files
tile_list <- whichTiles(region = points, polys = grid) 
#> Loading required namespace: rgeos
print(tile_list)
#> [1] "SECTOR_17" "SECTOR_5"  "SECTOR_1"  "SECTOR_10"

getMap(tiles = tile_list, polys = grid, dir = dir)
#> [1] "Cropping Tile SECTOR_17 (1 of 4)"
#> [1] "Cropping Tile SECTOR_5 (2 of 4)"
#> [1] "Cropping Tile SECTOR_1 (3 of 4)"
#> [1] "Cropping Tile SECTOR_10 (4 of 4)"
print(list.files(dir,recursive=TRUE, pattern = ".gz$"))
#> [1] "Elevations/SECTOR_1.gz"  "Elevations/SECTOR_10.gz"
#> [3] "Elevations/SECTOR_17.gz" "Elevations/SECTOR_5.gz"
```

Most of the functions later on in the workflow automatically call
`whichTiles` and `getMap` as needed, so you shouldn’t have to use these
functions much.

### Part 2: Velocity data sources

$$\\frac{d\\ell}{dt} = v\_{\\mathrm{max}} e^{-k \| \\frac{dz}{d\\ell} - \\alpha \|}$$
