#' Calculate dispersion indexes according to a given set of standards
#' and expectations (Mejia Ramon and Munson 2023), obtaining group, non-group, 
#' and total values for local observations and the global dataset.
#' 
#' @title Local Indicators of Dispersion
#' @param x A vector values
#' @param w A weights matrix of dimensions \code{length(x) x length(x)}
#' representing that a given observation \code{j} (along the columns)
#' is a part of \code{i} (along the rows)'s group. This can be the output of the 
#' \code{\link[lbmech]{makeWeights}} function. 
#' @param index A character string, either 'gini' (the default), or 'inoua', representing
#' whether distances are calculated in L1 or L2 space, respectively. Alternatively,
#' a numeric representing to what value distances and means are raised to when the
#' index is calculated. \code{index = 1} for Gini, and \code{index = 2} for Inoua. 
#' @param expect Either a character string or a matrix with dimensions \code{length(x) x length(x)},
#' representing the expectation value from which errors are calculated for each observation
#' pair between \code{i} and \code{j}. If \code{expect = 'self'}, the expectation is calculated as \code{i};
#' if \code{expect = 'local'}, the expectation is the neighborhood weighted mean; if 
#' if \code{expect = 'global'}, the expectation is the global mean. Expectations that depend on other
#' metrics (including hypothesis-driven that do not depend on the observed dataset) can
#' be provided by using an appropriate matrix. 
#' @param standard Either a character string or a matrix with dimensions \code{length(x) x length(x)},
#' representing the standard by which errors are judged by each observation
#' pair between \code{i} and \code{j}. If \code{standard = 'self'}, the standard is calculated as \code{i}; if 
#' \code{standard = 'other'}, the standard is calculated as \code{j};
#' if \code{standard = 'local'}, the standard is the neighborhood weighted mean; if 
#' if \code{standard = 'global'}, the standard is the global mean. Standards that depend on other
#' metrics (including hypothesis-driven that do not depend on the observed dataset) can
#' be provided by using an appropriate matrix. 
#' @param n A vector representing population weights. How much of an impact does a given 
#' observation have on any other observation regardless of its influence as provided
#' for in \code{w}. Default is \code{1} for all. 
#' @param fun.name If \code{index != c('gini','inoua',1,2)}, how should the function
#' be named? Default is \code{fun.name = paste0(index,'q')}.
#' @param type A character string, either the name or corresponding code of 
#' a particular standard-expectation pair, as defined in #Link to Mejia Ramon and Munson 2023#
#' @param max.cross When processing, what is the maximum number of rows that 
#' an internal data.table can have? This is generally not a concern unless
#' the number of observations approaches \code{sqrt(.Machine$integer.max)}--usually
#' about 2^31 for most systems. Lower values result in a greater number of chunks
#' thus allowing larger data.sets to be calculated.
#' @param canonical Should the canonical Gini or Inoua value also be calculated?
#' Default is \code{FALSE}, and is ignored if \code{index > 2}.
#' @param pb Logical. Should a progress bar be displayed? Default is \code{FALSE}, although
#' if a large dataset is processed that requires adjusting \code{max.cross} this can
#' be useful
#' @param clear.mem Logical. Should \code{\link[base]{gc}} be run in the middle of the 
#' calculation? Default is \code{clear.mem} but set as \code{TRUE} if memory limits are a concern. 
#' @importFrom data.table CJ
#' @importFrom data.table data.table
#' @return A list with the following entries:
#' 
#' (1) \code{$index} A named character string with the code of the index, named with its name
#' 
#' (2) \code{$local} A data.table, with three columns: \code{G_Gi}, the local group dispersion
#' index; \code{G_NGi}, the local non-group dispersion index; and \code{G_i}, the local
#' total dispersion index. Rows are in the same order as the input vector. This data.table
#' also contains the chosen expectations and standards as hidden attributes to be used by
#' \code{\link[lbmech]{inferLID}}.
#'  
#' (3) \code{$global} A list with three entries: \code{$G_G}, the global group dispersion index;
#' \code{$G_NG}, the global nongroup dispersion index; and \code{$G}, the global
#' total dispersion index. 
#' 
#' (4) \code{$canonical} The canonical Gini or Inoua index, if \code{canonical = TRUE} and 
#' \code{index < 3}. 
#' @examples 
#' 
#' # Generate dummy observations
#' x <- runif(10, 0, 100)
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
#' @export      
LID <- function(x, w, index = 'gini', expect = 'self', standard = 'global',
                n = rep(1, length(x)), 
                fun.name = paste0(index,'q') ,type = 'spatial', 
                max.cross = .Machine$integer.max, canonical = FALSE,
                pb = FALSE, clear.mem = gc()){
  # This bit to silence CRAN warnings
  Error=val_j=val_i=n_g=Denom=n_ng=G_Denom=NG_Denom=stand=G_Error=NG_Error=NULL
  G_i=G_Gi=..n=G_NGi=var=..subset_x=..subset_n=NULL
  # In what L-space are we working? If the user declares 'gini' or 'inoua', then
  # override function name and standard transformation entries
  if (index == 'gini' | index == 1){
    fun.name <- 'gini'
    index <- abs
    stand.trans = function(x) x
  } else if ((index == 'inoua' | index == 2)) {
    fun.name <- 'inoua'
    index <- function(x) x^2
    stand.trans = function(x) x^2
  } else {
    stand.trans = function(x) x^index
    index <- function(x) abs(x)^index
    canonical <- FALSE
  }
  fun.code <- toupper(substr(fun.name,1,1))
  
  # Calculate umber of subsets needed based on computational limits
  n_subsets <- ceiling(length(x)/floor(sqrt(max.cross)))
  subset_size <- floor(sqrt(max.cross))
  
  # Global value will be progressively added these values
  G <- 0
  NG <- 0
  
  # What's the name-code-expectation-standard relationship?
  combs <- data.table(Type = c('Spatial', 'Relative', 'Humble', 'Selfish',
                               'Upscaled', 'Local', 'Local-Radical','Local-Critical',
                               'Absolute', 'Downscaled','Global-Radical','Global-Critical'),
                      Code = c('S','R','H','E','U','L','Y','D','A','W','X','C'),
                      Expect = c('self','self','self','self',
                                 'local','local','local','local',
                                 'global','global','global','global'),
                      Standard = c('global','local','other','self',
                                   'global','local','other','self',
                                   'global','local','other','self'))
  
  # If a type is declared, override other decisions
  if (toupper(type) %in% toupper(combs$Type)){
    expect <- combs[toupper(type) == toupper(combs$Type)]$Expect
    standard <- combs[toupper(type) == toupper(combs$Type)]$Standard
  } else if (toupper(type) %in% toupper(combs$Code)){
    expect <- combs[toupper(type) == toupper(combs$Code)]$Expect
    standard <- combs[toupper(type) == toupper(combs$Code)]$Standard
  }
  
  # Non-group weights are inverse of group
  nw <- 1 - w
  
  # Population weights are multiplied by probability
  w <- t(n * t(w))
  nw <- t(n * t(nw))
  
  # Local values are appended to this data table 
  gini_list <- data.table()
  
  # Progress bar is optionsl
  if (pb) pb1 <- txtProgressBar(min=1,max=n_subsets+1)
  
  for (k in seq(n_subsets)){
    # Iterate over every subset, using as large a number of 'i's as we can, calculating
    # their relationship to all 'j's. To get the total value, we just add the 
    # values to a global index and local data.table as we go along. 
    if(pb) setTxtProgressBar(pb1,k)
    
    # Get the needed values for this subset
    subset_x <- x[((k-1)*subset_size+1):min(k*subset_size, length(x))]
    subset_n <- n[((k-1)*subset_size+1):min(k*subset_size, length(n))]
    subset_w <- w[((k-1)*subset_size+1):min(k*subset_size, ncol(w)),]
    subset_nw <- nw[((k-1)*subset_size+1):min(k*subset_size, ncol(nw)),]
    
    # Value is calculated from a cross-join
    dt <- CJ(J = x, I = subset_x, sorted = FALSE)
    dt <- data.table(I=(1:length(subset_x))[row(subset_w)],
                     J=(1:length(x))[col(subset_w)],w=c(subset_w), 
                     nw = c(subset_nw),
                     val_i = dt$I, val_j = dt$J)
    
    # Number of neighbors is the sum of group (and nongroup) weights
    dt[, `:=`(n_g = sum(w), n_ng = sum(nw)),by='I']
    
    # See the above table, but this giant chunk is just the linear algebra
    # needed to get the correct combination of expectation and standard for 
    # a particular choice. Additionally, it allows the user to deal with matrix
    # inputs. There's probably a less-verbose way of handling this, but this was
    # the quickest I could code quickly. 
    #
    # We start by overriding user-selected types and codes if non-matrix expectations and
    # standards were provided. It then calculates the error term based on the index
    # function, then denominator based on the chosen standard, and the total value
    # based on the expectation. Generally, the 'by' slot of the data.table is only
    # needed if a 'local' value is needed
    # 
    # If a matrix is provided for one or both of an expectation, then whatever the user
    # put into the 'type' and 'code' fields becomes dominant. 
    if (expect == 'self' & standard == 'global'){
      type <- 'Spatial'
      code <- 'S'
      dt[, Error := index(val_j - val_i)
      ][, `:=`(Denom = 2 * mean(stand.trans(val_j)))
      ][, `:=`(G_Gi = w/n_g * Error/Denom,
               G_NGi = nw/n_ng * Error/Denom)]
    } else if (expect == 'self' & standard == 'local'){
      type <- 'Relative'
      code <- 'R'
      dt[, Error := index(val_j - val_i)
      ][, `:=`(G_Denom = 2 * sum(w/n_g * stand.trans(val_j)),
               NG_Denom = 2 * sum(nw/n_ng * stand.trans(val_j))), by = 'I'
      ][, `:=`(G_Gi = w/n_g * Error/G_Denom,
               G_NGi = nw/n_ng * Error/NG_Denom)]
    } else if (expect == 'self' & standard == 'other'){
      type <- 'Humble'
      code <- 'Z'
      dt[, Error := index(val_j - val_i)
      ][, `:=`(Denom = 2 * stand.trans(val_j)),
      ][, `:=`(G_Gi = w/n_g * Error/Denom,
               G_NGi = nw/n_ng * Error/Denom)]
    } else if (expect == 'self' & standard == 'self'){
      type <- 'Selfish'
      code <- 'E'
      dt[, Error := index(val_j - val_i)
      ][, `:=`(Denom = 2 * stand.trans(val_i))
      ][, `:=`(G_Gi = w/n_g * Error/Denom,
               G_NGi = nw/n_ng * Error/Denom)]
    } else if (expect == 'self' & methods::is(standard,'matrix')){
      type <- stringr::str_to_sentence(type)
      code <- substr(type,1,1)
      
      standard <- data.table(J=(1:length(x))[row(standard)],
                             I=(1:length(x))[col(standard)],
                             stand=c(standard))
      dt <- merge(dt,standard,by=c('I','J'))
      dt[, Error := index(val_j - val_i)
      ][, `:=`(Denom = 2 * stand.trans(stand))
      ][, `:=`(G_Gi = w/n_g * Error/Denom,
               G_NGi = nw/n_ng * Error/Denom)]
    } else if ((expect == 'local') & (standard == 'global')){
      type <- 'Upscaled'
      code <- 'U'
      dt[, `:=`(G_Error = index(val_j - sum(w/n_g * val_j)),
                NG_Error = index(val_j - sum(nw/n_ng * val_j))), by = 'I'
      ][, `:=`(Denom = 2 * mean(stand.trans(val_j)))
      ][, `:=`(G_Gi = w/n_g * G_Error/Denom,
               G_NGi = nw/n_ng * NG_Error/Denom)]
    } else if (expect == 'local' & standard == 'local'){
      type <- 'Local'
      code <- 'L'
      dt[, `:=`(G_Error = index(val_j - sum(w/n_g * val_j)),
                NG_Error = index(val_j - sum(nw/n_ng * val_j))), by = 'I'
      ][, `:=`(G_Denom = 2 * sum(w/n_g * stand.trans(val_j)),
               NG_Denom = 2 * sum(nw/n_ng * stand.trans(val_j))), by = 'I'
      ][, `:=`(G_Gi = w/n_g * G_Error/G_Denom,
               G_NGi = nw/n_ng * NG_Error/NG_Denom)]
    } else if (expect == 'local' & standard == 'other'){
      type <- 'Local-Radical'
      code <- 'Y'
      dt[, `:=`(G_Error = index(val_j - sum(w/n_g * val_j)),
                NG_Error = index(val_j - sum(nw/n_ng * val_j))), by = 'I'
      ][, `:=`(Denom = 2 * stand.trans(val_j)),
      ][, `:=`(G_Gi = w/n_g * G_Error/Denom,
               G_NGi = nw/n_ng * NG_Error/Denom)]
    } else if (expect == 'local' & standard == 'self'){
      type <- 'Local-Critical'
      code <- 'D'
      dt[, `:=`(G_Error = index(val_j - sum(w/n_g * val_j)),
                NG_Error = index(val_j - sum(nw/n_ng * val_j))), by = 'I'
      ][, `:=`(Denom = 2 * stand.trans(val_i)),
      ][, `:=`(G_Gi = w/n_g * G_Error/Denom,
               G_NGi = nw/n_ng * NG_Error/Denom)]
    } else if (expect == 'local' & methods::is(standard,'matrix')){
      type <- stringr::str_to_sentence(type)
      code <- substr(type,1,1)
      standard <- data.table(J=(1:length(x))[row(standard)],
                             I=(1:length(x))[col(standard)],
                             stand=c(standard))
      dt <- merge(dt,standard,by=c('I','J'))
      dt[, `:=`(G_Error = index(val_j - sum(w/n_g * val_j)),
                NG_Error = index(val_j - sum(nw/n_ng * val_j))), by = 'I'
      ][, `:=`(Denom = 2 * stand)
      ][, `:=`(G_Gi = w/n_g * G_Error/Denom,
               G_NGi = nw/n_ng * NG_Error/Denom)]
    } else if (expect == 'global' & standard == 'global'){
      type <- 'Absolute'
      code <- 'A'
      dt[, Error := index(val_j - mean(val_j))
      ][, `:=`(Denom = 2 * mean(stand.trans(val_j)))
      ][, `:=`(G_Gi = w/n_g * Error/Denom,
               G_NGi = nw/n_ng * Error/Denom)]
    } else if (expect == 'global' & standard == 'local'){
      type <- 'Downscaled'
      code <- 'W'
      dt[, Error := index(val_j - mean(val_j))
      ][, `:=`(G_Denom = 2 * sum(w/n_g * stand.trans(val_j)),
               NG_Denom = 2 * sum(nw/n_ng * stand.trans(val_j))), by = 'I'
      ][, `:=`(G_Gi = w/n_g * Error/G_Denom,
               G_NGi = nw/n_ng * Error/NG_Denom)]
    } else if (expect == 'global' & standard == 'other'){
      type <- 'Global-Radical'
      code <- 'X'
      dt[, Error := index(val_j - mean(val_j))
      ][, `:=`(Denom = 2 * stand.trans(val_j)),
      ][, `:=`(G_Gi = w/n_g * Error/Denom,
               G_NGi = nw/n_ng * Error/Denom)]
    } else if (expect == 'global' & standard == 'self'){
      type <- 'Global-Critical'
      code <- 'C'
      dt[, Error := index(val_j - mean(val_j))
      ][, `:=`(Denom = 2 * stand.trans(val_i)),
      ][, `:=`(G_Gi = w/n_g * Error/Denom,
               G_NGi = nw/n_ng * Error/Denom)]
    } else if (expect == 'global' & methods::is(standard,'matrix')){
      type <- stringr::str_to_sentence(type)
      code <- substr(type,1,1)
      standard <- data.table(J=(1:length(x))[row(standard)],
                             I=(1:length(x))[col(standard)],
                             stand=c(standard))
      dt <- merge(dt,standard,by=c('I','J'))
      dt[, Error := index(val_j - mean(val_j))
      ][, `:=`(Denom = 2 * stand)
      ][, `:=`(G_Gi = w/n_g * Error/Denom,
               G_NGi = nw/n_ng * Error/Denom)]
    } else if (methods::is(expect,'matrix') & standard == 'global'){
      type <- stringr::str_to_sentence(type)
      code <- substr(type,1,1)
      expect <- data.table(J=(1:length(x))[row(expect)],
                           I=(1:length(x))[col(expect)],
                           stand=c(expect))
      dt <- merge(dt,expect,by=c('I','J'))
      dt[, Error := index(val_j - expect)
      ][, `:=`(Denom = 2 * mean(stand.trans(val_j)))
      ][, `:=`(G_Gi = w/n_g * Error/Denom,
               G_NGi = nw/n_ng * Error/Denom)]
    } else if (methods::is(expect,'matrix') & standard == 'local'){
      type <- stringr::str_to_sentence(type)
      code <- substr(type,1,1)
      expect <- data.table(J=(1:length(x))[row(expect)],
                           I=(1:length(x))[col(expect)],
                           stand=c(expect))
      dt <- merge(dt,expect,by=c('I','J'))
      dt[, Error := index(val_j - expect)
      ][, `:=`(G_Denom = 2 * sum(w/n_g * stand.trans(val_j)),
               NG_Denom = 2 * sum(nw/n_ng * stand.trans(val_j))), by = 'I'
      ][, `:=`(G_Gi = w/n_g * Error/G_Denom,
               G_NGi = nw/n_ng * Error/NG_Denom)]
    } else if (methods::is(expect,'matrix') & standard == 'other'){
      type <- stringr::str_to_sentence(type)
      code <- substr(type,1,1)
      expect <- data.table(J=(1:length(x))[row(expect)],
                           I=(1:length(x))[col(expect)],
                           stand=c(expect))
      dt <- merge(dt,expect,by=c('I','J'))
      dt[, Error := index(val_j - expect)
      ][, `:=`(Denom = 2 * stand.trans(val_j)),
      ][, `:=`(G_Gi = w/n_g * Error/Denom,
               G_NGi = nw/n_ng * Error/Denom)]
    } else if (methods::is(expect,'matrix') & standard == 'self'){
      type <- stringr::str_to_sentence(type)
      code <- substr(type,1,1)
      expect <- data.table(J=(1:length(x))[row(expect)],
                           I=(1:length(x))[col(expect)],
                           stand=c(expect))
      dt <- merge(dt,expect,by=c('I','J'))
      dt[, Error := index(val_j - expect)
      ][, `:=`(Denom = 2 * stand.trans(val_i)),
      ][, `:=`(G_Gi = w/n_g * Error/Denom,
               G_NGi = nw/n_ng * Error/Denom)]
    } else if (methods::is(expect,'matrix') & methods::is(standard,'matrix')){
      type <- stringr::str_to_sentence(type)
      code <- substr(type,1,1)
      expect <- data.table(J=(1:length(x))[row(expect)],
                           I=(1:length(x))[col(expect)],
                           stand=c(expect))
      dt <- merge(dt,expect,by=c('I','J'))
      standard <- data.table(J=(1:length(x))[row(standard)],
                             I=(1:length(x))[col(standard)],
                             stand=c(standard))
      dt <- merge(dt,standard,by=c('I','J'))
      dt[, Error := index(val_j - expect)
      ][, `:=`(Denom = 2 * stand)
      ][, `:=`(G_Gi = w/n_g * Error/Denom,
               G_NGi = nw/n_ng * Error/Denom)]
    } else {
      stop("Unknown expectation or standard")
    }
    
    # The local values are the weighted judgements cast by each i 
    dt[,G_i := (G_Gi * n_g/sum(..n) + G_NGi * n_ng/sum(..n))]
    gini <- dt[,lapply(.SD,sum),.SDcols = c("G_Gi","G_NGi","G_i"),by='I'
    ][,.SD,.SDcols = c("G_Gi","G_NGi","G_i")]
    gini[, var := ..subset_x][, n := ..subset_n]
    
    # Append values to the table
    gini_list <- rbind(gini_list, gini)
    if (clear.mem){
      rm(subset_x,subset_n,subset_w,subset_nw,dt)
      gc()
    }
  }
  gini <- gini_list
  rm(gini_list)
  if (clear.mem) gc()
  
  # If a matrix was provided, define the variable that will be added as an 
  # attribute to the output local data.table
  if (methods::is(standard,'matrix')){
    standard <- 'matrix'
  } 
  if (methods::is(expect,'matrix')){
    expect <- 'matrix'
  }
  
  # Index name is the index name plus the type, likewise for the code
  index_name <- paste0(code,fun.code)
  names(index_name) <- paste(type,stringr::str_to_sentence(fun.name))
  
  # Assign the names, expectations, and standards to the output table
  attributes(gini)[c("function","expectation","standard")] <- c(fun.name,expect,standard)
  
  # Output list contains the name, the local table, and the global values
  # (which themselves are just weighted means)
  out <- list(index = index_name,
              local = gini,
              global = list(
                G_G = sum(gini$G_Gi * gini$n/sum(gini$n)),
                G_NG = sum(gini$G_NGi * gini$n/sum(gini$n)),
                G = sum(gini$G_i * gini$n/sum(gini$n))))
  
  # Calculate the canonical value if appropriate, add to list
  if (canonical){
    out$canonical <- as.numeric(dispersionIndex(x, 
                                                w = n,
                                                index= stringr::str_remove(
                                                  tolower(names(index_name)),
                                                  "[a-z\\-]+ "),
                                                max.cross = max.cross))
  }
  if(pb) setTxtProgressBar(pb1,k+1)
  if (clear.mem) gc()
  return(out)
  
}