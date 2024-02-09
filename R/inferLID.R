
#' Infer whether there exists more or less within- and between-group
#' local and global inequality than would be expected versus if for all observations
#' the values of all other observations were permuted. This tests if local values
#' are significantly above or below what is expected given the global dataset, and
#' if global values are significantly above or below what is expected given an 
#' otherwise random distribution.
#' 
#' The output list can be passed to \code{\link[lbmech]{scatterLID}} to plot
#' the group and non-group components of local inequality based on the significance
#' classes.
#' 
#' @title Infer if dispersion is significant
#' @param lid The list output from the \code{\link[lbmech]{LID}} function.
#' @param w The same spatial weights matrix used in calculating the \code{lid} input.
#' @param ntrials The number of permutations to perform. Default is 999.
#' @param alpha Threshold for significance. Default is \code{alpha = 0.05}.
#' @param standard The standards matrix with dimensions \code{length(x) x length(x)} used 
#' when calculating \code{lid}. Ignored if none had been originally provided, otherwise
#' required. 
#' @param expect The expectations matrix with dimensions \code{length(x) x length(x)} used 
#' when calculating \code{lid}. Ignored if none had been originally provided, otherwise
#' required. 
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
#' @param clear.mem Logical. Should \code{\link[base]{gc}} be run in the middle of the 
#' calculation? Default is \code{clear.mem} but set as \code{TRUE} if memory limits are a concern. 
#' @importFrom data.table fifelse
#' @importFrom data.table data.table
#' @importFrom data.table as.data.table
#' @return A list with the following entries:
#' 
#' (1) \code{$local} A data.table with one column, indicating whether an observation is
#' falls in one of nine categories: Out-group High, Average, or Low for between-group
#' inequality, and In-group High, Average, or Low for within-group inequality based on
#' the significance according to the delta-J statistic in the \code{$stats} data.table. 
#' 
#' (2) \code{$global} A list with four entries, \code{$J_G} for the group component of the 
#' global inequality, \code{$J_NG} for the nongroup, \code{$J} for the total,
#' and $Class, containing the significance class for the global dataset. Each 
#' of the first three entries themselves contain three entriesL
#' \code{$delta}, representing the delta-J statistic, \code{$p}, representing its p-value,
#' and $Class, containing the group/non-group class
#' 
#' (3) \code{$stats} A data.table containing the number of permutations a randomly-calculated
#' \code{$J_Gi}, \code{$J_NGi}, or \code{$J_i} was above or below the real value
#' @examples 
#' 
#' # Generate dummy observations
#' x <- runif(10, 1, 100)
#' 
#' # Get distance matrix
#' dists <- dist(x)
#' 
#' # Get fuzzy weights considering 5 nearest neighbors based on 
#' # inverse square distance
#' weights <- makeWeights(dists, bw = 5, 
#'                        mode = 'adaptive', weighting = 'distance',
#'                        FUN = function(x) 1/x^2, minval = 0.1,
#'                        row.stand = 'fuzzy')
#'                        
#' # Obtain the 'local gini' value
#' lid <- LID(x, w = weights, index = 'gini', type = 'local')
#' 
#' # Infer whether values are significant relative to the spatial distribution
#' # of the neighbots
#' inference <- inferLID(lid, w = weights, ntrials = 100)
#' @export
inferLID <- function(lid, w, ntrials = 999, alpha = 0.05,
                     standard = NULL, expect = NULL, 
                     var.stand = FALSE, var.exp = FALSE, ng.invert = TRUE,
                     max.cross =.Machine$integer.max, pb = TRUE,
                     clear.mem = FALSE){
  # This bit to silence CRAN warnings
  Trial=..xrand=id=J_Gi=J_NGi=J_i=n=GroupClass=dGini_MC=dGini=NULL
  NonGroupClass=dNonGroup_MC=dNonGroup=.N=p=N=IndexClass=Class=NULL
  
  # Extract the values from the lid list object
  x <- as.data.table(lid$local)
  x$id <- 1:nrow(x)
  agg <- lid$global 
  index <- attributes(lid$local)$`function`
  mle <- attributes(lid$local)$`mle`
  
  # Make sure the user provided the right matrices
  if (attributes(lid$local)$`standard` == 'matrix'){
    if (is.null(standard)){ stop('Matrix must be provided if original standard was also a matrix')}
  } else {
    standard <- attributes(lid$local)$`standard`
  }
  if (attributes(lid$local)$`expectation` == 'matrix'){
    if (is.null(expect)) {stop('Matrix must be provided if original expectation was also a matrix')}
  } else {
    expect <- attributes(lid$local)$`expectation`
  }
  
  # Create a table to which we'll add the statistics for each permutation
  ptable <- data.table(Trial = rep(1:ntrials,each=nrow(x)))
  if (pb) pb1 <- txtProgressBar(min=0,max=ntrials+1,style=3)
  
  for (i in 1:ntrials){
    if (pb) setTxtProgressBar(pb1,i)
    
    # If you also vary the standard and expectations, do so
    if (var.stand){
      standard <- sample(standard)
    } 
    if (var.exp){
      expect <- sample(expect)
    }
    
    # Shuffling your neighbors but keeping the locations as the only
    # possible places is effectively shuffling the matrix row-wise while
    # keeping the diagonal constant
    w <- row(w) != col(w)
    w[w] <- stats::ave(w[w], row(w)[w], FUN = sample)
    
    # Perform the LID analysis on the randomized weights
    xrand <- data.table(id = x$id,
                        LID(x$var, w = w, n = x$n, index = index, mle = mle,
                            standard = standard,
                            expect = expect,
                            max.cross = max.cross,
                            clear.mem = clear.mem)$local)
    if (clear.mem) gc()
    # Append results to permutation table
    ptable[Trial == i, names(xrand) := ..xrand]
    if (clear.mem){
      rm(xrand)
      gc()
    }
  }
  
  # Significance for within-group inequality is defined based on the 
  # delta-J_Gi statistic, which is simply the group component minus the non-group
  # component
  group <-  merge(ptable[,.(id, dGini_MC = J_Gi - J_NGi)], 
                  x[,.(id, dGini = J_Gi - J_NGi)],
                  allow.cartesian=TRUE)
  
  # Significance for between-group inequality is defined based on the 
  # delta-J_NGi statistic, which is simply the non-group component minus the mean
  # total inequality
  nongroup <- merge(ptable[,.(id, dNonGroup_MC = J_NGi - sum(J_i * n,na.rm=TRUE)/sum(n,na.rm=TRUE)),by='Trial'],
                    x[,.(id, dNonGroup = J_NGi - sum(J_i*n,na.rm=TRUE)/sum(n,na.rm = TRUE))],
                    allow.cartesian=TRUE)
  
  # Combine into one table
  gini <- cbind(nongroup,group[,-1])
  
  # The within-group class is easy to infer, a higher value always means higher
  # inequality
  gini[, GroupClass := fifelse(dGini_MC > dGini, 'In-group Low','In-group High')]
  
  # The between group is a bit more difficult, and generally depends on whether
  # the self is being included in the non-group, or whether it's being excluded
  if ((standard != 'self' & expect != 'self') | 
      ((standard == 'matrix' | expect == 'matrix') & ng.invert)){
    gini[, NonGroupClass := fifelse(dNonGroup_MC > dNonGroup, 'Out-group High','Out-group Low')]
  } else if ((standard == 'self' | expect == 'self') | 
             ((standard == 'matrix' | expect == 'matrix') & !ng.invert)){
    gini[, NonGroupClass := fifelse(dNonGroup_MC > dNonGroup, 'Out-group Low','Out-group High')]  
  }
  
  # Calculate p values
  groupInference <- gini[,.N,by=c('id',"GroupClass")
  ][,p  := 1 - (1+N)/(1+sum(N)),by='id'][]
  
  nongroupInference <- gini[,.N,by=c('id',"NonGroupClass")
  ][,p  := 1 - (1+N)/(1+sum(N)),by='id'][]
  
  # Out-group p values are based on how extreme the permuted global value
  # is relative to the actual one
  globalp <- ptable[,.(J_Gi = average(J_Gi, w = x$n), 
                       J_NGi = average(J_NGi, w = x$n), 
                       J_i = average(J_i, w = x$n)),by='Trial'
  ][, GroupClass := fifelse(agg$J_G > J_Gi,'High','Low')
  ][, NonGroupClass := fifelse(agg$J_NG > J_NGi,'High','Low')
  ][, IndexClass := fifelse(agg$J < J_i,'High','Low')
  ]
  rm(ptable)
  # Calculate p values
  groupp <- globalp[, .N, by = c('GroupClass')
  ][,p  := 1 - (1+N)/(1+sum(N))][]
  
  nongroupp <- globalp[, .N, by = c('NonGroupClass')
  ][,p  := 1 - (1+N)/(1+sum(N))][]
  
  indexp <- globalp[, .N, by = c('IndexClass')
  ][,p  := 1 - (1+N)/(1+sum(N))][]
  
  # Add the p values to a list
  globalp <- list(J_G = list(delta = agg$J_G - mean(globalp$J_Gi),
                             p = min(groupp$p)),
                  J_NG = list(delta = agg$J_NG - mean(globalp$J_NGi),
                              p = min(nongroupp$p)),
                  J = list(delta = agg$J - mean(globalp$J_i),
                           p = min(indexp$p)))

  # Calculate significance classes for the complete dataset
  # Within-group is straight-forward
  if (!is.na(agg$J_G)){
  if (globalp$J_G$p > alpha){
    globalp$J_G$Class <- 'In-group Average'
  } else if (globalp$J_G$delta > 0) {
    globalp$J_G$Class <- 'In-group High'
  } else {
    globalp$J_G$Class <- 'In-group Low'
  }
  } else {
    globalp$J_G$Class <- 'In-group NaN'
  }
  
  # Across group is again a bit complicated, and once again depends on whether
  # the self is included in the comparisons. 
  if (!is.na(agg$J_NG)){
  if (globalp$J_NG$p > alpha){
    globalp$J_NG$Class <- 'Out-group Average'
  } else if ((standard != 'self' & expect != 'self') | 
             ((standard == 'matrix' | expect == 'matrix') & ng.invert)){
    if (globalp$J_NG$delta > 0) {
      globalp$J_NG$Class <- 'Out-group Low'
    } else {
      globalp$J_NG$Class <- 'Out-group High'
    }
  } else if ((standard == 'self' | expect == 'self') | 
             ((standard == 'matrix' | expect == 'matrix') & !ng.invert)){
    if (globalp$J_NG$delta > 0) {
      globalp$J_NG$Class <- 'Out-group High'
    } else {
      globalp$J_NG$Class <- 'Out-group Low'
    }
  }
  } else {
    globalp$J_NG$Class <- 'Out-group NaN'
  }
   globalp$Class <- paste0(globalp$J_G$Class,", ",globalp$J_NG$Class)
   
  # Add only the significant observations to a new data.table;
  # if an observation is not included for a particular group/nongroup group,
  # call it 'average'
  out <- merge(merge(x,
                     groupInference[p < alpha, .(id,GroupClass)],
                     all=TRUE,by='id'),
               nongroupInference[p < alpha, .(id,NonGroupClass)],
               all=TRUE,by='id')[is.na(GroupClass), GroupClass := 'In-group Average'
               ][is.na(NonGroupClass), NonGroupClass := 'Out-group Average'
               ][, Class := paste(GroupClass,NonGroupClass,sep=', ')
               ][]
  
  # Return to the original order since it tends to get shuffled
  out <- out[order(id)]
  out$id <- NULL
  
  # Calculate the statistic, again...
  gini[,`:=`(delta_G = dGini - dGini_MC,
             delta_NG = dNonGroup - dNonGroup_MC)]
  
  # The output table of p values only needs to include the ones above 0.5, since
  # there're only two options
  stats <- merge(groupInference[order(p),.SD[1],by='id'
  ][,.(id,Group_C = N, Group_p = p, GroupClass)],
  nongroupInference[order(p),.SD[1],by='id'
  ][,.(id,NonGroup_C = N, NonGroup_p = p, NonGroupClass)])
  stats$id <- NULL
  
  out <- list(local = out[,.SD,.SDcols = c('Class')],
              global = globalp,
              stats = stats)
  if (pb) setTxtProgressBar(pb1, i+1)
  return(out)
}
