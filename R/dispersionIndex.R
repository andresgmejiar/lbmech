#' Calculate weighted versions of the Gini and Inoua (2021) indexes as originally
#' defined in the econometric literature, using the half mean relative distance 
#' method.
#' 
#' @title Calculate Gini and Inoua indexes
#' @name dispersionIndex
#' @aliases gini inoua
#' @param x A vector of values
#' @param w A vector of weights with the same length as \code{x}
#' @param index A character string, either \code{'gini'}, or \code{'inoua'}, representing
#' whether distances are calculated in L1 or L2 space, respectively
#' @param weight.mean Logical. Should the mean values be weighted, or does the 
#' global depend exclusively on the observations? Default is \code{TRUE}.
#' @param inverse Logical. Should the value for the inverse weights be
#' calculated as well using binary decomposition? Default is \code{FALSE}.
#' This rarely makes sense if the weights are population-based, 
#' but it does if they're probability-based. 
#' @param max.cross When processing, what is the maximum number of rows that 
#' an internal data.table can have? This is generally not a concern unless
#' the number of observations approaches \code{sqrt(.Machine$integer.max)}--usually
#' about 2^31 for most systems. Lower values result in a greater number of chunks
#' thus allowing larger data.sets to be calculated
#' @param pb Logical. Should a progress bar be displayed? Default is \code{FALSE}, although
#' if a large dataset is processed that requires adjusting \code{max.cross} this can
#' be useful
#' @importFrom data.table CJ
#' @return A numeric of length 1 (if \code{inverse = FALSE}) or 2 (if \code{inverse = TRUE})
#' representing the requested index.
#' @examples 
#' 
#' # Generate dummy observations
#' x <- runif(10, 0, 100)
#' n <- runif(10, 0, 10)
#' 
#' # Calculate Gini index
#' gini(x)
#' 
#' # Calculate weighted Inoua index
#' inoua(x, w = n)
#' @export
dispersionIndex <- function(x, index = 'gini', w = rep(1,length(x)),
                            weight.mean = TRUE, inverse = FALSE,
                            max.cross = .Machine$integer.max,
                            pb = FALSE){
  # This bit to silence CRAN warnings
  Diff=J=Error=Denom=weight=Index=nonweight=NULL
  
  # The group value is the proportion of the weight
  w <- w / sum(w)
  
  # And the nongroup is its inverse
  nw <- (max(w) - w)/sum(max(w) - w)
  
  # Number of subsets needed given computational size limitations
  n_subsets <- ceiling(length(x)/floor(sqrt(max.cross)))
  subset_size <- floor(sqrt(max.cross))
  
  # We'll iterate over as large a number of 'i's as we can, calculating
  # their relationship to all 'j's. To get the total value, we just add the 
  # values to a global inde as we go along. 
  G <- 0
  NG <- 0
  
  # Progress bar is optionsl
  if (pb) pb1 <- txtProgressBar(min=1, max = n_subsets, style = 3)
  
  
  for(i in seq_len(n_subsets)){
    if (pb) setTxtProgressBar(pb1,i)
    # Iterate over every chunk, subsetting the relevant observations and
    # weights
    subset_x <- x[((i-1)*subset_size+1):min(i*subset_size, length(x))]
    subset_w <- w[((i-1)*subset_size+1):min(i*subset_size, length(w))]
    subset_nw <- nw[((i-1)*subset_size+1):min(i*subset_size, length(nw))]
    
    # Index requires a pairwise comparison, so require a cross-join
    dt <- CJ(I = subset_x, J = x, sorted = FALSE)[, Diff := J - I][]
    dt[, `:=`(weight = rep(w, length(subset_x)),
              nonweight = rep(nw, length(subset_x)))]
    
    if (weight.mean == TRUE){
      # If you weight the global mean, apply weight to denominator
      if (index == 'gini'){
        dt[, Error := abs(Diff)
        ][, Denom := 2* sum(weight) * sum(weight * J) , by = 'I'
        ][, Index := weight * Error / Denom]
        
        G <- G + sum(dt$Index * rep(subset_w, each = length(x)))
      } else if (index == 'inoua'){
        dt[, Error := Diff^2
        ][, Denom := 2* sum(weight) * sum(weight * J^2) , by = 'I'
        ][, Index := weight * Error / Denom] 
        
        G <- G + sum(dt$Index  * rep(subset_w, each = length(x)))
      }
    } else if (weight.mean == FALSE){
      # Otherwise pass
      if (index == 'gini'){
        dt[, Error := abs(Diff)
        ][, Denom := 2* sum(weight) * mean(J) 
        ][, Index := weight * Error / Denom]
        
        G <- G + sum(dt$Index *  rep(subset_w, each = length(x)))
      } else if (index == 'inoua'){
        dt[, Error := Diff^2
        ][, Denom := 2* sum(weight) * mean(J^2)
        ][, Index := weight * Error / Denom] 
        
        G <- G + sum(dt$Index * rep(subset_w, each = length(x)))
      }
    }
    
    if (inverse){
      # Calculate the non-group component if we want to see the 
      # opposite of the weights. 
      if (weight.mean == TRUE){
        if (index == 'gini'){
          dt[, Error := abs(Diff)
          ][, Denom := 2* sum(nonweight) * sum(nonweight * J) , by = 'I'
          ][, Index := nonweight * Error / Denom]
          
          NG <- NG + sum(dt$Index  * rep(subset_nw, each = length(x)))
        } else if (index == 'inoua'){
          dt[, Error := Diff^2
          ][, Denom := 2* sum(nonweight) * sum(nonweight * J^2) 
          ][, Index := nonweight * Error / Denom] 
          
          NG <- NG + sum(dt$Index * rep(subset_nw, each = length(x)))
        }
      } else if (weight.mean == FALSE){
        if (index == 'gini'){
          dt[, Error := abs(Diff)
          ][, Denom := 2* sum(nonweight) * mean(J)
          ][, Index := nonweight * Error / Denom]
          
          NG <- NG + sum(dt$Index * rep(subset_nw, each = length(x)))
        } else if (index == 'inoua'){
          dt[, Error := Diff^2
          ][, Denom := 2* sum(nonweight) * mean(J^2)
          ][, Index := nonweight * Error / Denom] 
          
          NG <- NG + sum(dt$Index * rep(subset_nw, each = length(x)))
        }
      }
    }
  }
  
  # Name the numerics by their index type, return to user
  if (!inverse){
    names(G) <- stringr::str_to_title(index)
    return(G)
  } else if (inverse){
    out <- c(G, NG)
    names(out) <- c(paste0(index, "-G"), paste0(index, "-NG"))
    return(out)
  }
}

#' @rdname dispersionIndex
#' @param ... Parameters to pass on to \code{\link[lbmech]{dispersionIndex}}.
#' @export
gini <- function(...){
  dispersionIndex(index = 'gini', ...)
}

#' @rdname dispersionIndex
#' @param ... Parameters to pass on to \code{\link[lbmech]{dispersionIndex}}.
#' @export
inoua <- function(...){
  dispersionIndex(index = 'inoua', ...)
}
