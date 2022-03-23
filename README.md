-   [lbmech](#lbmech)
    -   [Installation](#installation)
    -   [Example Workflow](#example-workflow)
        -   [Part 1: Topographic Data
            Sources](#part-1-topographic-data-sources)
        -   [Part 2: Velocity Data
            Sources](#part-2-velocity-data-sources)
        -   [Part 3: Preparing the World](#part-3-preparing-the-world)
        -   [Part 4: Getting Costs, Paths, and
            Corridors](#part-4-getting-costs-paths-and-corridors)
    -   [Quick Start](#quick-start)
        -   [For Energetic and Temporal Cost
            Calculations](#for-energetic-and-temporal-cost-calculations)
        -   [Custom Cost Functions and Multivariate
            Costs](#custom-cost-functions-and-multivariate-costs)

<!-- README.md is generated from README.Rmd. Please edit that file -->

# lbmech <img src='man/figures/logo.png' align="right" height="139" />

<!-- badges: start -->
<!-- badges: end -->

`lbmech` is a geospatial package for least-cost path analysis in `R`. It
contains additional tools to calculate time- and energy-based
costs-of-travel for humans and animals moving across the landscape. The
general philosophy behind this package is that least-cost path analyses
should always be simple, the most time-taking parts should be done at
most once, and ideally that costs should be rooted in empirical reality.

In general terms, both `lbmech` and the deterministic functions of the
library `gdistance` provide similar capabilities but with notable
differences in computation and ease-of-use.`gdistance` stores movement
costs as their reciprocal (*conductance* ,
i.e.B ![1/\\textrm{resistance}](https://latex.codecogs.com/png.latex?1%2F%5Ctextrm%7Bresistance%7D "1/\textrm{resistance}"))
in a sparse matrix representing every possible transition between two
cells on a raster. The use of conductance allows for most cell-to-cell
transitions to be zero, since impossible transitions due to distance
have a conductance of
![1/\\infty = 0](https://latex.codecogs.com/png.latex?1%2F%5Cinfty%20%3D%200 "1/\infty = 0").
However, the use of conductance makes it cumbersome to employ functions
where f(0) != 0, requiring the use of index masking. This in turn
encounters integer overflow errors with b