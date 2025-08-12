




# kmeans function with a set centroid value


#' 
#' kmeans with a set centroid value
#' one dimensional values (do not need rn)
#' and only one of the centroids is fixed
#' also probably bugs, did this by memory...
#' @param values (vector) of numeric values to assign to clusters
#' @param k (numeric) number of clusters
#' @param fixed_center (numeric) single value to fix one center at
#' @param initial_centroid (vector) of numeric values (of length k=#clusters)
#' @param max_iter (numeric) integer of max number of kmeans update iterations
#' @param threshold (numeric) stopping threshold when updating centroid
#'                  calculated by mean((new_centroid - centroid)**2)
#' @returns list with named results
#'     centroid (vector of length k indicating centers for each cluster)
#'     cluster_assignment (vector with elements 1:k of cluster assignments)
#'     centroid_change (change in the centroid update of the last update)
#' @examples     
#' # good example: one cluster is centered at 0 
#' test = myKmeans(values = c(rnorm(80, mean = -1, sd = 1),
#'                            rnorm(100, mean = 0, sd = 1),
#'                            rnorm(150, mean = 3, sd = .8)),
#'                 k = 3,
#'                 fixed_center = 0,
#'                 initial_centroid = NULL,
#'                 max_iter = 100,
#'                 threshold = .0001)
#' test$centroid; test$centroid_change
#' plot(test$cluster_assignment)
#' 
#' 
#' 
#' # questionable example: no cluster is centered at 0 
#' test = myKmeans(values = c(rnorm(80, mean = -1, sd = 1),
#'                            rnorm(100, mean = 1, sd = 1),
#'                            rnorm(150, mean = 3, sd = .8)),
#'                 k = 3,
#'                 fixed_center = 0,
#'                 initial_centroid = NULL,
#'                 max_iter = 100,
#'                 threshold = .0001)
#' test$centroid; test$centroid_change
#' plot(test$cluster_assignment)
myKmeans <- function(values, k, 
                     fixed_center=0,
                     initial_centroid=NULL, 
                     max_iter = 1000,
                     threshold = .001) {
  
  # debugging
  # values = c(rnorm(80, mean = -1, sd = 1),
  #            rnorm(100, mean = 0, sd = 1),
  #            rnorm(150, mean = 3, sd = .8))
  # k = 3
  # fixed_center = 0
  # initial_centroid = NULL
  # max_iter = 10
  # threshold = .0001
  
  # check inputs
  # TODO: assertthat statements
  if(is.null(initial_centroid)) { # should use regular kmeans centroid result
      initial_centroid = seq(from = quantile(values, probs =   1/k), 
                             to   = quantile(values, probs = 1-1/k), 
                             length.out = k)
  }
  
  # initialize variables
  iter = 0
  centroid_change = 999
  centroid = fixOneCentroid(centroid=initial_centroid, fixed_center=fixed_center)
  
  
  # kmeans updates
  while(iter < max_iter && centroid_change > threshold) {
    # print(sprintf('%d: %.2f change, centroid: [%s]', 
    #               iter, centroid_change,
    #               paste0(round(centroid, 2), collapse = ', ')))
    
    
    cluster_assignment = updateClusterAssignment(values=values, centroid=centroid)
    new_centroid       = updateCentroid(values=values, cluster_assignment=cluster_assignment, k=k)
    new_centroid       = fixOneCentroid(centroid=new_centroid, fixed_center=fixed_center)
    
    
    centroid_change    = mean((new_centroid - centroid)**2)
    centroid = new_centroid
    iter = iter + 1
  }
  
  # warnings if results are questionable
  if(iter >= max_iter) {
    warning('kmeans updating exceeded max_iter(ations)')
  }
  
  if(centroid_change > threshold) {
    warning('centroid change did not meet threshold')
  }
  
  return(list(centroid=centroid,
              cluster_assignment=cluster_assignment,
              centroid_change=centroid_change))
  
}

#' Fix one of the centroid values to the specified center
#' by finding the closest centroid value and changing that
#' to the fixed_center
#' @param centroid (vector) of numeric values (of length k=#clusters)
#' @param fixed_center (numeric) single value to fix one center at
fixOneCentroid <- function(centroid, fixed_center) {
  centroid[which.min(abs(centroid - fixed_center))] = fixed_center
  centroid
}


#' update the centroid based on current values' cluster assignment
#' @param values (vector) of numeric values to assign to clusters
#' @param cluster_assignment (vector) of integers (1:k) denoting cluster 
#'                           assignment of each of the values
#' @param k (numeric) integer
updateCentroid <- function(values, cluster_assignment, k) {
  # hmm can we do without loading in dplyr
  centroid = rep(NA, k)
  cluster_assignments = unique(cluster_assignment)
  for(k_ in cluster_assignments) {
    centroid[k_] = mean(values[which(cluster_assignment == k_)])
  }
  
  # if any centers are NA (none were assigned to this cluster)
  # randomly choose value in the middle of values
  if(anyNA(centroid)) {
    centersNA = which(is.na(centroid))
    for(k_ in centersNA) {
      centroid[k_] = runif(n = 1, min = quantile(values, probs = .3), max = quantile(values, probs = .7))
    }
  }
  
  
  return(centroid)
}


#' Update cluster assignment based on the centroid
#' assigns each value to the closest clusters' centroid value
#' @param values (vector) of numeric values to assign to clusters
#' @param centroid (vector) of numeric values (of length k=#clusters)
updateClusterAssignment <- function(values, centroid) {
  distances = sapply(X = centroid,
                     FUN = function(center) {abs(values - center)})
  
  # head(distances);  values[1] - centroid # check
  
  cluster_assignment = apply(X = distances, 
                             FUN = which.min, 
                             MARGIN = 1)
  
  # cluster_assignment[1:3]; head(distances, 3)
  
  return(cluster_assignment)
}






