% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/normMDS.R
\name{normMDS}
\alias{normMDS}
\title{Normalized Non-Metric Multidimensional Scaling}
\usage{
normMDS(dists, xy = NULL, id = "ID", ax = 1:2, chir = 3, x = "x", y = "y", ...)
}
\arguments{
\item{dists}{A list (named, ideally) of distance matrices with the same points}

\item{xy}{(Optional) The points in geographic space. This can be a 
data.frame, data.table, SpatVect, or SpatialPointsDataFrame object. If provided, 
axes orientations and chirality will be taken with respect to \code{xy}.}

\item{id}{Character. The name of the column in \code{xy} containing the unique
name or id for each observation. Default is \code{xy = 'ID'}.}

\item{ax}{The index (if numeric) or Unique ID (if character) of the two points
defining the reference axis. The first point will be placed at the origin.
If \code{xy = NULL}, the second point will be placed at v = 0. Default is 
\code{ax = 1:2}.}

\item{chir}{The index (if numeric) or Unique ID (if character) of the point
defining the reference handendess/chirality. This must be a point not in \code{ax}.
Default is \code{ax = 3}}

\item{x}{Character. The name of the column in \code{xy} containing the x 
coordinate for each observation. Default is \code{xy = 'x'}. Ignored if 
\code{xy} is a SpatVect or Spatial object.}

\item{y}{Character. The name of the column in \code{xy} containing the y 
coordinate for each observation. Default is \code{xy = 'y'}. Ignored if 
\code{xy} is a SpatVect or Spatial object.}

\item{...}{Additional parameters to pass to \code{\link[MASS]{isoMDS}}}
}
\value{
A named list, with each entry containing a data.table of estimated
pseudo-coordinates in \emph{(u,y)} space. One unit is equal to the distance
between both points in the \code{ax} parameter. If \code{xy} was provided,
the first entry \code{'ref'} will be the geographic data.
}
\description{
A function to compare implied 2D topologies based on different distance
matrices for the same data points
}
\details{
The function relies on \code{\link[MASS]{isoMDS}} from the MASS package to 
generate pseudo-coordinates for distance matrices using a 2D non-metric 
Multi Dimensional Scaling. Since the output orientation can be unpredictable,
normMDS aligns all topologies based on a common user-defined axis 
(the ends of which will appear equidistant and equioriented on an output plot),
and selects a chiral reference point to ensure that plots are oriented in a 
common orientation.
}
\examples{

# Create fake data
xy <- data.table(ID = c('A','B','C','D','E'),
                 x = runif(5, -100, 100),
                 y = runif(5, -100, 100))
                 
# Calculate a distance matrix
dists <- as.matrix(stats::dist(xy[,.(x,y)], upper = TRUE, diag = TRUE) )
rownames(dists) <- xy$ID
colnames(dists) <- xy$ID               
                 
                 
# Make a list of five fake distance matrices from five fake coordinates
dists <- list()
for (i in 1:5){
  fake <- data.table(ID = c('A','B','C','D','E'),
                 x = runif(5, -100, 100),
                 y = runif(5, -100, 100))
                 
  fake <- as.matrix(stats::dist(fake[,.(x,y)], upper = TRUE, diag = TRUE) )
  rownames(fake) <- xy$ID
  colnames(fake) <- xy$ID 
  dists[[i]] <- fake
}

# Get pseudo-coordinates
normMDS(dists, xy = xy)
}
