#' Determine the bandwidth that maximizes the non-group component of inequality.
#'
#' @title Incremental Local Indicators of Dispersion
#' @param x A vector of weights with the same length as \code{x}
#' @param dist A matrix or distance object representing pairwise distances. The
#' distances need not be symmetrical.
#' @param bws A vector containing the representing the bandwidth within neighbors 
#' are considered. If \code{mode = 'adaptive'}, \code{bw} is the number of nearest neighbors. 
#' If \code{mode = 'fixed'}, \code{bw} is the radius of the window in the map units.
#' @param n A vector representing population weights. How much of an impact does a given 
#' observation have on any other observation regardless of its influence as provided
#' for in \code{w}. Default is \code{1} for all. 
#' @param ntrials The number of permutations to perform. Default is 50.
#' @param alpha Threshold for significance. Default is \code{alpha = 0.05}.
#' @param standard The standards matrix with dimensions \code{length(x) x length(x)} used 
#' when calculating \code{lid}. Ignored if none had been originally provided, otherwise
#' required. 
#' @param expect The expectations matrix with dimensions \code{length(x) x length(x)} used 
#' when calculating \code{lid}. Ignored if none had been originally provided, otherwise
#' required. 
#' @param mode One of \code{'adaptive'}, which considers a \code{bw} number of nearest 
#' neighbors; or \code{'fixed'}, which considers a fixed bandwidth window of radius \code{bw}.
#' @param weighting One of \code{'membership'}, which considers binary membership 
#' such that neighbors are weighted 1 and non-neighbors 0; \code{'distance'} which 
#' weighs neighbors according to \code{FUN} with the raw distance matrix providing the
#' distance; or \code{'rank'} which uses the rank-distance (i.e. 1 for nearest neighbor,
#' 2 for second nearest...) as the distance variable.
#' @param FUN The distance function. Default is \code{NULL} for \code{'membership'}, and
#' \code{function(x) offset/(offset + x)} otherwise. Ignored if \code{x} is a vector.
#' @param inf.val When singularities arise, (i.e. whenever the value is 1/0), by what value are
#' they replaced? Default is the \code{FUN} of the lowest non-\code{minval} value.
#' Ignored if \code{x} is a vector.
#' @param def.neigh Numeric. At what distance (in the map units) are observations definitely neighbors?
#' All distances are subtracted by this value, and all resulting distances less than zero are reassigned
#' to \code{minval}. 
#' @param offset What value is added to the denominator to prevent singularities from arising
#' (e.g. whenever the value is 1/0)? Larger values imply smaller distance-decay. This should be
#' a numeric of length one or \code{length(def.neigh)}. Alternatively, \code{offset} can be expressed
#' as a function of \code{def.neigh}. Default is \code{offset = function(x) 2 * x}. 
#' Default is \code{offset = 0}. Ignored if \code{x} is a vector.
#' @param minval When distances are raw, what is the minimum allowable distance?
#' Default is \code{0}. Ignored if \code{x}. Use this if you don't want to offset values otherwise. 
#' @param row.stand Logical or \code{'fuzzy'}. If \code{TRUE} (the default), rows are standardized such 
#' that they sum to one. If \code{'fuzzy'}, rows are standardized as a proportion of the 
#' largest value. 
#' @param var.stand Logical. Should the standards be permuted if a matrix was 
#' provided? Default is \code{FALSE}.
#' @param var.exp Logical. Should the expectations be permuted if a matrix was
#' provided? Default is \code{FALSE}.
#' @param ng.invert Does a higher non-group value imply higher between group inequality?
#' Default is \code{TRUE}. This is ignored if matrixes were not originally provided, as 
#' it is automatically performed.
#' @param max.cross When processing, what is the maximum number of rows that 
#' an internal data.table can have? This is generally not a concern unless
#' the number of observations approaches \code{sqrt(.Machine$integer.max)}--usually
#' about 2^31 for most systems. Lower values result in a greater number of chunks
#' thus allowing larger data.sets to be calculated.
#' @param pb Logical. Should a progress bar be displayed? Default is \code{FALSE}, although
#' if a large dataset is processed that requires adjusting \code{max.cross} this can
#' be useful
#' @param ... Additional parameters to pass on to \code{\link[lbmech]{LID}}.
#' @importFrom data.table as.data.table
#' @importFrom data.table setnames
#' @importFrom data.table shift
#' @return A list with three entries:
#' 
#' (1) \code{index} A named character with the code of the index named by its name
#' 
#' (2) \code{$bws} The bandwidths that significantly optimize the non-group inequality. 
#' Generally, a neighborhood is the first significant peak. 
#' 
#' (3) \code{$stats} A data.table with the global group, non-group, and total values 
#' for each bandwidth, as well as a column indicating whether or not it's significant.
#' @examples 
#' 
#' # Generate dummy observations
#' x <- runif(10, 1, 100)
#' 
#' # Get distance matrix
#' dists <- dist(x)
#' 
#' # Bandwidth sizes from 3 to 5
#' bws <- 3:5
#' 
#' inc <- incrementalLID(x, dist = dists, bws = bws, index = 'gini', type = 'local',
#'                       weighting = 'distance', FUN = function(x) 1/x^2, minval = 1)
#' @export
incrementalLID <- function(x, dist, bws = Inf, def.neigh = 0, 
                           offset = function(x) 2 * x, n = rep(1,length(x)),
                           ntrials = 50, alpha = 0.05, standard = NULL, expect = NULL,
                           mode = 'adaptive', weighting = 'membership', FUN = NULL,
                           inf.val = NULL, row.stand = 'fuzzy', minval = 50,
                           var.stand = FALSE, var.exp = FALSE, ng.invert = TRUE,
                           max.cross = .Machine$integer.max, pb = TRUE,
                           ...){
  # This bit to silence CRAN warnings
  delta_J_NG=dt_plus=dt_minus=d2t=NG_Class=`p J_NG`=NULL
  
  if (class(offset) == 'function'){
    offset = offset(def.neigh)
  } else {
    offset = rep(offset,length(def.neigh))
  }
  
  if (mode == 'membership' & sum(is.infinite(bws)) > 0) {
    stop("argument 'bw' is has infinities")
  } 
  
  if (length(bws) != length(def.neigh)){
    if (length(bws) == 1 & length(def.neigh) != 1) {
      bws <- rep(bws, length(def.neigh))
      dominant = 'def.neigh'
    } else if (length(bws) != 1 & length(def.neigh) == 1){
      def.neigh <- rep(def.neigh, length(bws))
      dominant = 'bws'
    } else {
      stop("'bws' and 'def.neigh' must either be of equal length, or at most one may have length 1.")
    }
  } else {
    dominant = 'both'
  }
  
  # Generate table to which values will be appended
  outTable <- data.table()
  if (pb) pb1 <- txtProgressBar(min = 0, max = length(bws),style = 3)
  i = 0
  
  
  
  for (i in 1:length(bws)){
    if (pb) setTxtProgressBar(pb1,i)
    # Make the weights for a given bandwidth
    
    weights <- makeWeights(dist, bw = bws[i], mode = mode, weighting = weighting, 
                           FUN = FUN, offset = offset[i], def.neigh = def.neigh[i], inf.val = inf.val,
                           row.stand = row.stand, minval = minval)
    
    
    
    # Calculate the LID for the given bandwidth
    lid <- LID(x = x, n = n, w = weights, standard = standard, expect = expect,
               max.cross = max.cross, ...)
    if (is.null(standard)) stand <- attributes(lid$local)$standard
    if (is.null(expect)) expct <- attributes(lid$local)$expectation
    
    # Get the local values and calculate the delta-J statistic for the global
    local <- lid$local
    
    vals <- c("delta_J_G" = sum((lid$local$J_Gi - lid$local$J_NGi) * lid$local$n)/sum(lid$local$n),
              "delta_J_NG"= sum((lid$local$J_NGi - lid$global$J) * lid$local$n)/sum(lid$local$n),
              "delta_J"  =  sum((lid$local$J_i - lid$global$J) * lid$local$n)/sum(lid$local$n))
    
    # Perform the inference for the bandwidth 
    invisible(utils::capture.output(
      {inference <- inferLID(lid, w = weights, ntrials = ntrials, alpha = alpha,
                             standard = standard, expect = expect, 
                             var.stand = var.stand, var.exp = var.exp, 
                             ng.invert = ng.invert,
                             max.cross = max.cross)
      }
    ))
    
    # P values come from the global statistic
    suppressWarnings(ps <- as.data.table(inference$global)[2])
    setnames(ps,names(ps)[1:3], c("p J_G","p J_NG","p J"))
    
    # Append p values to table
    outTable <- rbind(outTable,
                      as.data.table(c(bw = bws[i],def.neigh = def.neigh[i],
                                      offset = offset[i], vals,ps)))
    
    i <- i + 1
  }
  
  outTable[, names(outTable) := lapply(.SD,as.numeric)]
  
  # Calculate second derivative (ish)
  outTable[, `:=`(dt_minus = shift(delta_J_NG) -delta_J_NG , dt_plus = shift(delta_J_NG,type='lead') - delta_J_NG)
  ][, `:=`(dt_plus = dt_plus/abs(dt_plus), dt_minus = dt_minus/abs(dt_minus))
  ][, d2t := dt_plus + dt_minus]
  
  # Identify significant peaks
  outTable[, NG_Class := fifelse(`p J_NG` < alpha, "Significant","Not Significant")]
  
  
  if ((stand != 'self' & expct != 'self') | 
      ((stand == 'matrix' | expct == 'matrix') & ng.invert)){
    if (dominant == 'bw') {
      bw <- outTable[d2t == 2 & NG_Class == 'Significant', ]$bw
    } else if (dominant == 'def.neigh'){
      def.neigh <- outTable[d2t == 2 & NG_Class == 'Significant', ]$def.neigh
    }
    
    if (dominant == 'bw' & length(bw) == 0){
      bw <- outTable[which(delta_J_NG == outTable[NG_Class == 'Significant',
                                                  min(delta_J_NG,na.rm = TRUE)]),bw]
    } 
    if (dominant == 'def.neigh' & length(def.neigh) == 0){
      def.neigh <- outTable[which(delta_J_NG == outTable[NG_Class == 'Significant',
                                                         min(delta_J_NG,na.rm = TRUE)]),def.neigh]
    }
  } else if ((stand == 'self' | expct == 'self') | 
             ((stand == 'matrix' | expct == 'matrix') & !ng.invert)){
    if (dominant == 'bw'){
      bw <- outTable[d2t == -2  & NG_Class == 'Significant',]$bw
    } else if (dominant == 'def.neigh'){
      def.neigh <- outTable[d2t == -2 & NG_Class == 'Significant', ]$def.neigh
    }
    
    if (dominant == 'bw' & length(bw) == 0){
      bw <- outTable[which(delta_J_NG == outTable[NG_Class == 'Significant',
                                                  max(delta_J_NG,na.rm = TRUE)]),bw] 
    } 
    
    if (dominant == 'def.neigh' & length(def.neigh) == 0){
      def.neigh <- outTable[which(delta_J_NG == outTable[NG_Class == 'Significant',
                                                         max(delta_J_NG,na.rm = TRUE)]),def.neigh] 
    }
  }
  
  if (dominant == 'bw'){
    out <- list(index = lid$index,
                bw = bw,
                stats = outTable[,c(1:9,13)])
  } else if (dominant == 'def.neigh'){
    out <- list(index = lid$index,
                def.neigh = def.neigh,
                stats = outTable[,c(1:9,13)])
  }
  
  if (pb) setTxtProgressBar(pb1,i)
  return(out)
}
