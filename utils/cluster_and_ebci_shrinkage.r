# util helper functions for clustering and ebci shrinkage

# require(dplyr)  # dataframe manipulation
# require(ggplot2)# plotting
# require(ebci)   # robust shrinkage
# require(mlpack) # kmeans
# require(mclust) # mixture modeling



# ==============================================================================
# ============================ MAIN CALLS ======================================
# ==============================================================================

#' Clusters and performs EBCI with separate datasets for training the clusters
#' and shrinkage
#'
#' Performs the specified clustering on the tstat (unless ebci_formula specifies
#' otherwise) of the train_idx test set. Then, performs EBCI shrinkage on the
#' remaining indices.
#'
#'
#' @param df_cluster (dataframe) with values to cluster on
#'        must have the columns of the input: value_colname
#' @param df_shrink (dataframe) with values to shrink on
#'        must have the columns of the input: value_colname, se_colname,
#' @param value_colname (string) column name of the value
#'                      (to cluster and shrink)
#' @param se_colname (string) column name of the standard error of the value
#' @param cluster_method (string) one of c('kmeans', 'mclust', 'preclustered')
#'                       for the method of clustering
#'        if 'preclustered', must input parameter precluster
#'        'kmeans' from package mlpack
#'        'mclust' from package mclust
#' @param k (integer) number of clusters, NA for some default chosen by BIC
#' @param precluster (vector or matrix) needed if cluster_method='preclustered'
#'
#'           cluster assignment (vector of length #tests ith values 1, ..., k)
#'           or
#'           cluster probability (matrix of #tests x k, each ro must sum to 1)
#' @param precluster_centers (vector) of cluster centers (the point to shrink
#'                           towards for each cluster)
#'
#' @param ebci_formula (string) ignore! actually, do not allow different formulas
#'                              formula for EBCI call,
#'                              leave NA for 'tstat - shrinkage_point ~ 0'
#'     NOTE!! if not NA, no care will be taken for any transformations done.
#'     e.g. 'tstat - shrinkage_point ~ 0' needs shrinkage_point to be added
#'          back to get the original estimate, but this is not coded for custom
#'          formulas
#' @param alpha (numeric) (0,1) significance level of ebci calculations
#'              (used for the length of the CI)
#' @param other_covariates (character vector) of column names of effects_df to
#'        keep the info in the final dataframe (e.g. test: type of test)
#' @output dataframe with columns:
#'             cluster, pr1, pr2, ..., prk,
#'             (test, other covariates, ...),
#'             unshrunk_estimate, unshrunk_estimate_se,
#'             shrunk_estimate, lower_ci, upper_ci
#'             unshrunk_lower_ci, unshrunk_upper_ci
cluster_and_ebci_w_split <- function(df_cluster, df_shrink,
                                     value_colname, se_colname,
                                     cluster_method,
                                     k,
                                     precluster=NA, precluster_centers=NA,
                                     # ebci_formula = NA,
                                     alpha=.1, other_covariates = c('test'), se_threshold=.0001) {
  
  # effects_df = poisson_effects |> head(100);
  # cluster_method = 'kmeans';
  # k=2; ebci_formula = NA;
  # precluster=NA; precluster_centers=NA;
  # alpha = .1;
  # other_covariates = c('test');
  # se_threshold = .0001
  
  assertthat::assert_that(all(c(value_colname) %in% colnames(df_cluster)),
                          msg = "'value_colname' must be columns of df_cluster")
  assertthat::assert_that(all(c(value_colname, se_colname) %in% colnames(df_shrink)),
                          msg = "'value_colname' and 'se_colname' must be columns of df_shrink")
  
  assertthat::assert_that(cluster_method %in% c('myKmeans', 'kmeans', 'mclust', 'preclustered'),
                          msg = "cluster_method must be one of {'myKmeans', 'kmeans', 'mclust', 'preclustered'}")
  assertthat::assert_that( is.na(k) || ((k %% 1 == 0) && k > 0),
                           msg = "if k is not NA, k must be a whole number and positive")
  assertthat::assert_that( cluster_method != 'preclustered' || (!all(is.na(precluster)) && !all(is.na(precluster_centers))),
                           msg = "if cluster_method is 'preclustered', precluster and precluster_centers must be provided")
  assertthat::assert_that(  cluster_method != 'kmeans' || !is.na(k),
                            msg = "when cluster_method='kmeans', k must be specified (cannot be NA)")
  
  
  
  
  
  
  # === Mixture Model/Clustering
  
  
  # cluster + get: shrinkage_point, classification, cl_prob_matrix
  switch (cluster_method,
          'mclust' = {
            # mclust: clustering - specify k number of clusters or use BIC to choose
            if(is.na(k)) {
              BIC = mclust::mclustBIC(  data = df_cluster[,value_colname], verbose = FALSE)
              # plot(BIC); summary(BIC)
              
              clustfit = mclust::Mclust(data = df_cluster[,value_colname], x = BIC, verbose = FALSE)
            } else {
              clustfit = mclust::Mclust(data = df_cluster[,value_colname], G = k  , verbose = FALSE)
            }
            
            
            
            classification_train  = clustfit$classification # classification vector
            cl_prob_matrix_train  = clustfit$z
            # shrink towards weighted average of the k clusters' mean (weights=probs)
            #         shrinkage points = cluster probabilities %*% cluster centers
            #             (#tests x 1) =  (#tests x #clusters) %*% (#clusters x 1)
            shrinkage_point_train = clustfit$z %*% clustfit$parameters$mean
            
            # get assignments on test dataset
            clustfit_test = mclust::predict.Mclust(clustfit, newdata = df_shrink[,value_colname])
            classification_test = clustfit_test$classification
            cl_prob_matrix_test = clustfit_test$z
            shrinkage_point_test= clustfit_test$z %*% clustfit$parameters$mean
            
          },
          'myKmeans' = {
            
            mlpackKmeans =  mlpack::kmeans(cluster = k,
                                        input = df_cluster[,value_colname, drop = FALSE],
                                        kmeans_plus_plus = TRUE, 
                                        max_iterations = 100) # decrease max iter to save time
            
            kmeansfit = myKmeans(values = df_cluster[,value_colname], 
                                 k = k, 
                                 fixed_center = 0, 
                                 initial_centroid = mlpackKmeans$centroid |> as.vector(), 
                                 max_iter = 1000, 
                                 threshold = .0001)
            
            
            classification_train = kmeansfit$cluster_assignment # classification vector
            cl_prob_matrix_train = one_hot_encode(classification_train) # construct a 'prob' matrix, but it is just 1's and 0's
            # shrink towards weighted average of the k clusters' mean (weights=probs)
            #         shrinkage points = cluster probabilities %*% cluster centers
            #             (#tests x 1) =  (#tests x #clusters) %*% (#clusters x 1)
            shrinkage_point_train =  cl_prob_matrix_train %*% kmeansfit$centroid
            
            
            # get assignments on test dataset
            classification_test = get_closest_cluster(values = df_shrink[,value_colname], 
                                                      centroids = kmeansfit$centroid)
            cl_prob_matrix_test = one_hot_encode(classification_test)
            shrinkage_point_test= cl_prob_matrix_test %*% kmeansfit$centroid
            
          },
          'kmeans' = {
            
            kmeansfit =  mlpack::kmeans(cluster = k,
                                        input = df_cluster[,value_colname, drop = FALSE],
                                        kmeans_plus_plus = TRUE)
            
            
            classification_train = kmeansfit$output[,2] + 1 # classification vector
            cl_prob_matrix_train = one_hot_encode(classification_train) # construct a 'prob' matrix, but it is just 1's and 0's
            # shrink towards weighted average of the k clusters' mean (weights=probs)
            #         shrinkage points = cluster probabilities %*% cluster centers
            #             (#tests x 1) =  (#tests x #clusters) %*% (#clusters x 1)
            shrinkage_point_train =  cl_prob_matrix_train %*% kmeansfit$centroid
            
            
            # get assignments on test dataset
            classification_test = get_closest_cluster(values = df_shrink[,value_colname], 
                                                      centroids = kmeansfit$centroid)
            cl_prob_matrix_test = one_hot_encode(classification_test)
            shrinkage_point_test= cl_prob_matrix_test %*% kmeansfit$centroid
            
          },
          'preclustered' = {
            assertthat::assert_that(cluster_method != 'preclustered',
                                    msg = "sorry, cluster_method='preclustered' not implemented")
            # NOT IMPLEMENTED FOR THIS FUNCTION
            # if(length(dim(precluster)) == 2) {
            #   # if precluster is a matrix of cluster probabilities
            #   cl_prob_matrix = precluster
            #   classification = apply(X = cl_prob_matrix, FUN = which.max, MARGIN = 1)
            # 
            # } else {
            #   # if precluster is a vector of classification assignments
            #   classification = precluster
            #   cl_prob_matrix = one_hot_encode(classification)
            # }
            # 
            # 
            # # shrink towards weighted average of the k clusters' mean (weights=probs)
            # #         shrinkage points = cluster probabilities %*% cluster centers
            # #             (#tests x 1) =  (#tests x #clusters) %*% (#clusters x 1)
            # assertthat::assert_that(ncol(cl_prob_matrix) == length(precluster_centers),
            #                         msg = 'precluster_centers must have length equal to the number of clusters')
            # shrinkage_point =  cl_prob_matrix %*% precluster_centers
            # 
            # 
            # # split into train and test
            # classification_train  = classification[train_idx]
            # cl_prob_matrix_train  = cl_prob_matrix[train_idx, ]
            # shrinkage_point_train = shrinkage_point[train_idx]
            # 
            # classification_test  = classification[test_idx]
            # cl_prob_matrix_test  = cl_prob_matrix[test_idx, ]
            # shrinkage_point_test = shrinkage_point[test_idx]
            # 
            # rm(shrinkage_point, classification, cl_prob_matrix); gc(verbose = FALSE)
            
          }
  )
  
  # # add to effects_df the shrinkage point
  # df_cluster$shrinkage_point   = shrinkage_point_train
  # df_shrink$shrinkage_point    = shrinkage_point_test
  # 
  # # add to effects_df the probability of clusters and hard cluster classification
  # df_cluster = cbind(df_cluster,
  #                    data.frame(cl_prob_matrix_train) |> `colnames<-`(paste0('pr', 1:ncol(cl_prob_matrix_train))),
  #                    data.frame(cluster = classification_train)
  #              )
  # 
  # df_shrink = cbind(df_shrink,
  #                    data.frame(cl_prob_matrix_test) |> `colnames<-`(paste0('pr', 1:ncol(cl_prob_matrix_test))),
  #                    data.frame(cluster = classification_test)
  #              )
  # 
  # 
  # rm(shrinkage_point_train, classification_train, cl_prob_matrix_train)
  # rm(shrinkage_point_test,  classification_test,  cl_prob_matrix_test); gc(verbose = FALSE)
  
  # === Shrinkage
  df_shrink = df_shrink[se_colname > se_threshold, ]  # remove tests with too small se
  
  # if(is.na(ebci_formula)) {
  #   ebci_formula = 'tstat - shrinkage_point ~ 0'
  # }
  
  ebci_obj = ebci(formula = sprintf('%s - shrinkage_point ~ 0', value_colname),
                  # formula = 'tstat - shrinkage_point ~ 0', # ebci_formula
                  data    = df_shrink |> 
                    select(all_of(c(value_colname, se_colname))) |> 
                    mutate(shrinkage_point = shrinkage_point_test) |>
                    # add new columns bc the formula is annoying and wants the
                    # se and weights specified not as a string but by the names
                    # of the cols themselves and has issues ... annoying
                    mutate(my_se      = eval(parse(text = se_colname)),      
                           my_weights = eval(parse(text = sprintf('1/%s^2', se_colname)))), 
                  
                  # se = eval(parse(text = se_colname)), 
                  # weights =  eval(parse(text = sprintf('1/%s^2', se_colname))), 
                  se = my_se,
                  weights = my_weights,
                  alpha = alpha)
  
  
  
  ebci_res = cbind(df_shrink |> select(all_of(other_covariates)),
                   data.frame(cl_prob_matrix_test) |> `colnames<-`(paste0('pr', 1:ncol(cl_prob_matrix_test))),
                   data.frame(       cluster = classification_test,
                                     shrinkage_point = shrinkage_point_test,
                                     unshrunk_value = df_shrink[,value_colname],
                                     se = df_shrink[,se_colname],
                                     shrunk_value = ebci_obj$df$th_eb + shrinkage_point_test)
                   
  )
  
  
  
  # # ignore if statement, always add shrinkage point back
  # if(TRUE) {
  #   # if(is.na(ebci_formula)) { # default 'tstat - shrinkage_point ~ 0', add back shrinkage_point
  #   ebci_res = ebci_res |>
  #                   mutate(shrunk_tstat = ebci_obj$df$th_eb + test_df$shrinkage_point)
  # } else { # don't transform the shrunk estimates
  #   ebci_res = ebci_res |>
  #                   mutate(shrunk_tstat = ebci_obj$df$th_eb)
  # }
  
  
  ebci_res =  ebci_res |>
    mutate(lower_ci = shrunk_value - ebci_obj$df$len_eb,
           upper_ci = shrunk_value + ebci_obj$df$len_eb)
  # unshrunk_estimate_upper_ci = unshrunk_estimate + qnorm(alpha/2)*unshrunk_estimate_se,
  # unshrunk_estimate_lower_ci = unshrunk_estimate - qnorm(alpha/2)*unshrunk_estimate_se)
  ebci_res = ebci_res |> mutate(train01 = 0, .before = 1) # was the test in the training idx
  # df_cluster$train01 = 1
  
  
  
  
  # add test_df values
  ebci_res = dplyr::bind_rows(ebci_res,
                              cbind(df_cluster |> select(all_of(other_covariates)),
                                    data.frame(cl_prob_matrix_train) |> `colnames<-`(paste0('pr', 1:ncol(cl_prob_matrix_train))),
                                    data.frame(        cluster = classification_train,
                                                       shrinkage_point = shrinkage_point_train, 
                                                       unshrunk_value = df_cluster[,value_colname],
                                                       se = df_cluster[,se_colname],
                                                       train01 = 1)
                              ) |>
                                select(any_of(colnames(ebci_res))))
  # return(ebci_res)
  return(list(ebci_res=ebci_res,
              ebci_obj=ebci_obj))
}


#' Clusters and performs EBCI
#' 
#' Performs the specified clustering on the tstat (unless ebci_formula specifies
#' otherwise) of the train_idx test set. Then, performs EBCI shrinkage on the
#' remaining indices.
#' 
#'
#' @param effects_df (dataframe) with unshrunk effects (#tests x #other info)
#'        must have columns: estimate, se
#' @param cluster_method (string) one of c('kmeans', 'mclust', 'preclustered')
#'                       for the method of clustering
#'        if 'preclustered', must input parameter precluster
#'        'kmeans' from package mlpack
#'        'mclust' from package mclust
#' @param k (integer) number of clusters, NA for some default chosen by BIC
#' @param precluster (vector or matrix) needed if cluster_method='preclustered'
#'
#'           cluster assignment (vector of length #tests ith values 1, ..., k)
#'           or
#'           cluster probability (matrix of #tests x k, each ro must sum to 1)
#' @param precluster_centers (vector) of cluster centers (the point to shrink
#'                           towards for each cluster)
#' @param train_idx (vector) of indices (for effects_df- \in 1:nrow(effects_df))
#'           for training the clusters on. The remaining indices will be used as 
#'           test samples for ebci shrinkage.
#'           if train_idx=NULL (default), then do not split into train and test
#'           set, and perform clustering and shrinkage on all.
#' 
#' @param ebci_formula (string) ignore! actually, do not allow different formulas
#'                              formula for EBCI call,
#'                              leave NA for 'tstat - shrinkage_point ~ 0'
#'     NOTE!! if not NA, no care will be taken for any transformations done.
#'     e.g. 'tstat - shrinkage_point ~ 0' needs shrinkage_point to be added
#'          back to get the original estimate, but this is not coded for custom
#'          formulas
#' @param alpha (numeric) (0,1) significance level of ebci calculations
#'              (used for the length of the CI)
#' @param other_covariates (character vector) of column names of effects_df to
#'        keep the info in the final dataframe (e.g. test: type of test)
#' @output dataframe with columns:
#'             cluster, pr1, pr2, ..., prk,
#'             (test, other covariates, ...),
#'             unshrunk_estimate, unshrunk_estimate_se,
#'             shrunk_estimate, lower_ci, upper_ci
#'             unshrunk_lower_ci, unshrunk_upper_ci
#' @example
#' train_idx = 1:1000
#' # Testing: kmeans
#' poisson_effects_ebci_kmeans2 =
#'   cluster_and_ebci(effects_df = poisson_effects |> head(100),
#'                    cluster_method = 'kmeans',
#'                    k=2, 
#'                    train_idx = train_idx,
#'                    # ebci_formula = NA,
#'                    precluster=NA, precluster_centers=NA,
#'                    alpha = .1,
#'                    other_covariates = c('test'),
#'                    se_threshold = .0001)
#'
#' poisson_effects_ebci_kmeans3 =
#'   cluster_and_ebci(effects_df = poisson_effects |> head(100),
#'                    cluster_method = 'kmeans',
#'                    k=3,  
#'                    train_idx = train_idx,
#'                    # ebci_formula = NA,
#'                    precluster=NA, precluster_centers=NA,
#'                    alpha = .1,
#'                    other_covariates = c('test'),
#'                    se_threshold = .0001)
#'
#' # Testing: mclust
#' poisson_effects_ebci_mclust2 =
#'   cluster_and_ebci(effects_df = poisson_effects |> head(100),
#'                    cluster_method = 'mclust',
#'                    k=2,  
#'                    train_idx = train_idx,
#'                    # ebci_formula = NA,
#'                    precluster=NA, precluster_centers=NA,
#'                    alpha = .1,
#'                    other_covariates = c('test'),
#'                    se_threshold = .0001)
#'
#' # Testing: mclust + k=NA (detect #clusters with BIC)
#' poisson_effects_ebci_mclustNA =
#'   cluster_and_ebci(effects_df = poisson_effects |> head(100),
#'                    cluster_method = 'mclust',
#'                    k=NA,  
#'                    train_idx = train_idx,
#'                    # ebci_formula = NA,
#'                    precluster=NA, precluster_centers=NA,
#'                    alpha = .1,
#'                    other_covariates = c('test'),
#'                    se_threshold = .0001)
#'
#' # Testing: precluster w/ vector of classifications
#' sceptre_effects_kmeans2 =  mlpack::kmeans(
#'                cluster = 2,
#'                input = sceptre_effects |> head(100) |> select(estimate),
#'                kmeans_plus_plus = TRUE)
#'
#' sceptre_effects_ebci_precl =
#'   cluster_and_ebci(effects_df = sceptre_effects |> head(100),
#'                    cluster_method = 'preclustered',
#'                    k=NA,  
#'                    train_idx = train_idx,
#'                    # ebci_formula = NA,
#'                    precluster=sceptre_effects_kmeans2$output[, 2] + 1,
#'                    precluster_centers=sceptre_effects_kmeans2$centroid,
#'                    alpha = .1,
#'                    other_covariates = c('test'),
#'                    se_threshold = .0001)
#'
#' # Testing: precluster w/ matrix of probabilities
#' sceptre_effects_mclust = mclust::Mclust(data = sceptre_effects  |> head(100)  |> select(tstat), verbose = FALSE)
#'
#' sceptre_effects_ebci_precl =
#'     cluster_and_ebci(effects_df = sceptre_effects |> head(100),
#'                    cluster_method = 'preclustered',
#'                    k=NA,  
#'                    train_idx = train_idx,
#'                    # ebci_formula = NA,
#'                    precluster=sceptre_effects_mclust$z,
#'                    precluster_centers=sceptre_effects_mclust$parameters$mean,
#'                    alpha = .1,
#'                    other_covariates = c('test'),
#'                    se_threshold = .0001)
cluster_and_ebci <- function(effects_df,
                             cluster_method,
                             k,
                             precluster=NA, precluster_centers=NA,
                             train_idx=NULL,
                             # ebci_formula = NA,
                             alpha=.1, other_covariates = c('test'), se_threshold=.0001) {
  
  # effects_df = poisson_effects |> head(100);
  # cluster_method = 'kmeans';
  # k=2; ebci_formula = NA;
  # precluster=NA; precluster_centers=NA;
  # alpha = .1;
  # other_covariates = c('test');
  # se_threshold = .0001
  
  assertthat::assert_that(all(c('estimate', 'se') %in% colnames(effects_df)),
                          msg = "'estimate' and 'se' must be columns of effects_df")
  assertthat::assert_that(cluster_method %in% c('kmeans', 'mclust', 'preclustered'),
                          msg = "cluster_method must be one of {'kmeans', 'mclust', 'preclustered'}")
  assertthat::assert_that( is.na(k) || ((k %% 1 == 0) && k > 0),
                           msg = "if k is not NA, k must be a whole number and positive")
  assertthat::assert_that( cluster_method != 'preclustered' || (!all(is.na(precluster)) && !all(is.na(precluster_centers))),
                           msg = "if cluster_method is 'preclustered', precluster and precluster_centers must be provided")
  assertthat::assert_that(  cluster_method != 'kmeans' || !is.na(k),
                            msg = "when cluster_method='kmeans', k must be specified (cannot be NA)")
  assertthat::assert_that(is.null(train_idx) || (all(train_idx %in% 1:nrow(effects_df))),
                          msg = "train_idx must be NULL or a vector of indices (elements in 1:nrow(effects_df))")
  
  # create tstat (test statistic) column if not present
  if(! 'tstat' %in% colnames(effects_df)) {
    effects_df$tstat = effects_df$estimate / effects_df$se
  }
  
  
  if(is.null(train_idx)) {
    train_idx = 1:nrow(effects_df)
    test_idx  = 1:nrow(effects_idf)
  } else {
    test_idx = setdiff(1:nrow(effects_df), train_idx)
  }
  
  # prob more efficient to keep original df and subset but whatever
  train_df = effects_df[train_idx, ]
  test_df  = effects_df[test_idx, ]
  rm(effects_df); gc() # remove effects_df to save some space
  
  # === Mixture Model/Clustering
  
  
  # cluster + get: shrinkage_point, classification, cl_prob_matrix
  switch (cluster_method,
          'mclust' = {
            # mclust: clustering - specify k number of clusters or use BIC to choose
            if(is.na(k)) {
              BIC = mclust::mclustBIC(data = train_df$tstat, verbose = FALSE)
              # plot(BIC); summary(BIC)
              
              clustfit = mclust::Mclust(data = train_df$tstat, x = BIC, verbose = FALSE)
            } else {
              clustfit = mclust::Mclust(data = train_df$tstat, G = k  , verbose = FALSE)
            }
            
            
            
            classification_train  = clustfit$classification # classification vector
            cl_prob_matrix_train  = clustfit$z
            # shrink towards weighted average of the k clusters' mean (weights=probs)
            #         shrinkage points = cluster probabilities %*% cluster centers
            #             (#tests x 1) =  (#tests x #clusters) %*% (#clusters x 1)
            shrinkage_point_train = clustfit$z %*% clustfit$parameters$mean
            
            # get assignments on test dataset
            clustfit_test = mclust::predict.Mclust(clustfit, newdata = test_df$tstat)
            classification_test = clustfit_test$classification
            cl_prob_matrix_test = clustfit_test$z
            shrinkage_point_test= clustfit_test$z %*% clustfit$parameters$mean
            
          },
          'kmeans' = {
            
            kmeansfit =  mlpack::kmeans(cluster = k,
                                        input = train_df |> select(tstat),
                                        kmeans_plus_plus = TRUE)
            
            
            classification_train = kmeansfit$output[,2] + 1 # classification vector
            cl_prob_matrix_train = one_hot_encode(classification_train) # construct a 'prob' matrix, but it is just 1's and 0's
            # shrink towards weighted average of the k clusters' mean (weights=probs)
            #         shrinkage points = cluster probabilities %*% cluster centers
            #             (#tests x 1) =  (#tests x #clusters) %*% (#clusters x 1)
            shrinkage_point_train =  cl_prob_matrix_train %*% kmeansfit$centroid
            
            
            # get assignments on test dataset
            classification_test = get_closest_cluster(values = test_df$tstat, centroids = kmeansfit$centroid)
            cl_prob_matrix_test = one_hot_encode(classification_test) 
            shrinkage_point_test= cl_prob_matrix_test %*% kmeansfit$centroid
            
          },
          'preclustered' = {
            
            if(length(dim(precluster)) == 2) {
              # if precluster is a matrix of cluster probabilities
              cl_prob_matrix = precluster
              classification = apply(X = cl_prob_matrix, FUN = which.max, MARGIN = 1)
              
            } else {
              # if precluster is a vector of classification assignments
              classification = precluster
              cl_prob_matrix = one_hot_encode(classification)
            }
            
            
            # shrink towards weighted average of the k clusters' mean (weights=probs)
            #         shrinkage points = cluster probabilities %*% cluster centers
            #             (#tests x 1) =  (#tests x #clusters) %*% (#clusters x 1)
            assertthat::assert_that(ncol(cl_prob_matrix) == length(precluster_centers),
                                    msg = 'precluster_centers must have length equal to the number of clusters')
            shrinkage_point =  cl_prob_matrix %*% precluster_centers
            
            
            # split into train and test
            classification_train  = classification[train_idx]
            cl_prob_matrix_train  = cl_prob_matrix[train_idx, ]
            shrinkage_point_train = shrinkage_point[train_idx]
            
            classification_test  = classification[test_idx]
            cl_prob_matrix_test  = cl_prob_matrix[test_idx, ]
            shrinkage_point_test = shrinkage_point[test_idx]
            
            rm(shrinkage_point, classification, cl_prob_matrix); gc(verbose = FALSE)
            
          }
  )
  
  # add to effects_df the shrinkage point
  train_df$shrinkage_point = shrinkage_point_train
  test_df$shrinkage_point  = shrinkage_point_test
  
  # add to effects_df the probability of clusters and hard cluster classification
  train_df = cbind(train_df,
                   data.frame(cl_prob_matrix_train) |> `colnames<-`(paste0('pr', 1:ncol(cl_prob_matrix_train))),
                   data.frame(cluster = classification_train)
  )
  
  test_df = cbind(test_df,
                  data.frame(cl_prob_matrix_test) |> `colnames<-`(paste0('pr', 1:ncol(cl_prob_matrix_test))),
                  data.frame(cluster = classification_test)
  )
  
  
  rm(shrinkage_point_train, classification_train, cl_prob_matrix_train)
  rm(shrinkage_point_test,  classification_test,  cl_prob_matrix_test); gc(verbose = FALSE)
  
  # === Shrinkage
  test_df = test_df |> filter(se > se_threshold)  # remove tests with too small se
  
  # if(is.na(ebci_formula)) {
  #   ebci_formula = 'tstat - shrinkage_point ~ 0'
  # }
  
  ebci_obj = ebci(formula = 'tstat - shrinkage_point ~ 0', # ebci_formula
                  data    = test_df |> mutate(se = 1),  # set se=1 bc tstat
                  se = se, weights = 1/se^2, alpha = .1)
  
  
  
  ebci_res = test_df |> select(shrinkage_point,
                               all_of(c(other_covariates,
                                        grep('pr[0-9]*', colnames(test_df), value = TRUE), # paste0('pr', 1:k),
                                        'cluster')),
                               unshrunk_tstat       = tstat,
                               estimate    = estimate,
                               estimate_se = se)
  
  # ignore if statement, always add shrinkage point back
  if(TRUE) {
    # if(is.na(ebci_formula)) { # default 'tstat - shrinkage_point ~ 0', add back shrinkage_point
    ebci_res = ebci_res |>
      mutate(shrunk_tstat = ebci_obj$df$th_eb + test_df$shrinkage_point)
  } else { # don't transform the shrunk estimates
    ebci_res = ebci_res |>
      mutate(shrunk_tstat = ebci_obj$df$th_eb)
  }
  
  
  ebci_res =  ebci_res |>
    mutate(lower_ci = shrunk_tstat - ebci_obj$df$len_eb,
           upper_ci = shrunk_tstat + ebci_obj$df$len_eb) 
  # unshrunk_estimate_upper_ci = unshrunk_estimate + qnorm(alpha/2)*unshrunk_estimate_se,
  # unshrunk_estimate_lower_ci = unshrunk_estimate - qnorm(alpha/2)*unshrunk_estimate_se)
  ebci_res = ebci_res |> mutate(train01 = 0, .before = 1) # was the test in the training idx
  train_df$train01 = 1
  
  
  # add test_df values
  ebci_res = dplyr::bind_rows(ebci_res,
                              train_df |>
                                rename(unshrunk_tstat = tstat,
                                       estimate_se    = se) |>
                                select(any_of(colnames(ebci_res))))
  return(ebci_res)
}




#' gathers unshrunk effects and calculates ebci shrunk effects and 
#' assembles them together in a dataframe
#' 
#' @param effects_df (dataframe) with unshrunk effects (#tests x #other info)
#'        must have columns: estimate, se
#' @param cluster (vector) of cluster membership (length = #tests)
#' @param alpha (numeric) (0,1) significance level of ebci calculations 
#'              (used for the length of the CI)
#' @param other_covariates (character vector) of column names of effects_df to 
#'        keep the info in the final dataframe (e.g. test: type of test)
#' @output dataframe with columns: 
#'             cluster, (test, other covariates, ...), 
#'             unshrunk_estimate, unshrunk_estimate_se, 
#'             shrunk_estimate, lower_ci, upper_ci
#' @example
#' ebci_by_cluster(effects_df = glm_effects, 
#'                 cluster = glm_effect_kmeans3$output[, 2] + 1, 
#'                 alpha = .1, 
#'                 other_covariates = c('test'), 
#'                 se_threshold = .0001)
ebci_by_cluster <- function(effects_df, cluster, alpha=.1, other_covariates = c('test'), se_threshold=.0001) {
  
  effects_df$cluster = cluster # make new col indicating cluster
  effects_df = effects_df |> arrange(cluster) |> # sort by cluster
    filter(se > se_threshold)  # remove tests with too small se
  
  clusters = unique(effects_df$cluster)
  
  # perform ebci on each cluster
  shrunk_estimate = c()
  lower_ci        = c()
  upper_ci        = c()
  for(cl in clusters) {
    cl_ebci = ebci(formula = estimate ~ 1, 
                   data = effects_df |> filter(cluster == cl), 
                   se = se, weights = 1/se^2, alpha = alpha)
    shrunk_estimate = c(shrunk_estimate, cl_ebci$df$th_eb )
    lower_ci        = c(       lower_ci, cl_ebci$df$th_eb - cl_ebci$df$len_eb)
    upper_ci        = c(       upper_ci, cl_ebci$df$th_eb + cl_ebci$df$len_eb)
  }
  
  
  return(
    effects_df |> select(cluster, all_of(other_covariates),
                         unshrunk_estimate    = estimate,
                         unshrunk_estimate_se = se) |>
      mutate(shrunk_estimate      = shrunk_estimate,
             lower_ci             = lower_ci,
             upper_ci             = upper_ci)
  )
} 




# ==============================================================================
# ============================ PLOTTING FNS ====================================
# ==============================================================================

#' Plot glm hist by cluster group and by test type
#'
#' @param effects_df (dataframe) with unshrunk effects (#tests x #other info)
#'        must have columns: estimate, se
#' @param kmeans_result (kmeans object) result of kmeans clustering
#' @param title (string) title of plot
#' 
#' @example
#' sceptre_effects = read.csv('../saves/sceptre_effects.csv')
#' sceptre_effects_kmeans2 =  
#'        mlpack::kmeans(cluster = 2, 
#'                       input = sceptre_effects |> select(estimate),
#'                       kmeans_plus_plus = TRUE)
#' p_sceptre_kmeans2  = plot_kmeans(effects_df = sceptre_effects, 
#'                                  kmeans_result = sceptre_effects_kmeans2, 
#'                                  title = 'SCEPTRE (k=2)')
plot_kmeans <- function(effects_df, kmeans_result,
                        title    = 'Histogram of Unshrunk Estimates') {
  # kmeans_result = poisson_effect_kmeans2
  # effects_df = poisson_effects
  
  
  plot_glm_df = rbind(effects_df |> select(estimate, category = test) |> mutate(facet = 'test'),
                      data.frame(estimate = effects_df$estimate,
                                 category = kmeans_result[['output']][, 2] + 1,
                                 facet = 'kmeans')) |>
    mutate(category = factor(category, levels = c('negative', 'positive', 'discovery', seq(1:nrow(kmeans_result$centroid)))),
           facet = factor(facet, levels = c('test', 'kmeans')))
  
  
  p_kmeans = ggplot(plot_glm_df) +
    geom_histogram(aes(x = estimate, fill = category, y = after_stat(density)), 
                   position = 'dodge',
                   bins = 25) +
    geom_vline(aes(xintercept = 0), color = 'black', linetype = 'solid', linewidth = .2) +
    # geom_vline(data = data.frame(xint=kmeans_result[['centroid']][1, 1], facet="kmeans") |> 
    #              mutate(facet = factor(facet, levels = c('test', 'kmeans'))), 
    #            aes(xintercept = xint), color = 'gray30', linetype = 'dashed') +
    # geom_vline(data = data.frame(xint=kmeans_result[['centroid']][2, 1], facet="kmeans") |> 
    #              mutate(facet = factor(facet, levels = c('test', 'kmeans'))), 
    #            aes(xintercept = xint), color = 'gray30', linetype = 'dashed') +
    # geom_vline(aes(xintercept = sceptre_effect_kmeans$centroid[1, 1]), color = 'gray30', linetype = 'dashed') +
    # geom_vline(aes(xintercept = sceptre_effect_kmeans$centroid[2, 1]), color = 'gray30', linetype = 'dashed') +
    labs(title    = title,
         subtitle = 'grouped by AY test type and kmeans cluster', 
         x        = 'Unshrunk Estimate',
         fill     = 'Group') +
    scale_fill_discrete(type = c( 'darkgreen', 'firebrick', 'gold3', 
                                  gray.colors(kmeans_result$centroid |> nrow(), start = .1, end = .7)), # 'gray52', 'gray22'),
                        guide = guide_legend(direction = 'vertical')) +
    facet_grid(vars(facet)) +
    
    theme_bw()  +
    theme(legend.position = 'right', 
          panel.grid.minor.y = element_blank()) 
  
  return(p_kmeans)
}

#' Plot EBCI results
plot_ebci <- function(ebci_results, 
                             subsample_idx=c(seq(from =    1, to = 2947, by = floor(2947/80)),
                                             seq(from = 2948, to = 3493, by = floor( 546/110)),
                                             seq(from = 3494, to = 6256, by = floor(5980/110)))) {
  
  p_subsampled = ggplot(ebci_results |>
                          slice(c(seq(from =    1, to = 2947, by = floor(2947/80)),
                                  seq(from = 2948, to = 3493, by = floor( 546/110)),
                                  seq(from = 3494, to = 6256, by = floor(5980/110)))) |>
                          mutate(x = 1:n()), 
                        aes(x = x)) +
    # geom_point(aes(y = lower_ci) )+
    # geom_point(aes(y = upper_ci)) +
    geom_point(aes(y = -1.5, color = factor(test, levels = c('negative', 'positive', 'discovery')))) +
    geom_hline(aes(yintercept = 0)) +
    geom_segment(aes(y = lower_ci, yend = upper_ci), 
                 lineend = 'square', linewidth = 1, alpha = .7, color = 'deepskyblue4') +
    geom_segment(aes(y = unshrunk_estimate, yend = .99 * shrunk_estimate + .01 * unshrunk_estimate), 
                 lineend = 'square', linejoin = 'bevel', arrow = arrow(length = unit(0.2,"cm")), linewidth = .3, alpha = .7, color = 'deepskyblue2') +
    geom_point(aes(y = unshrunk_estimate), color = 'deepskyblue2') +
    geom_point(aes(y = shrunk_estimate), color = 'deepskyblue4') +
    scale_x_continuous(expand = c(0, 0), breaks = seq(0, 1000, by = 10), 
                       # limits = c(1, 90)) +
                       # limits = c(65, 175)) +
                       # limits = c(95, 400)) +
                       # limits = c(1, 400)) +
    ) +
    scale_y_continuous(expand = c(0.025, 0), limits = c(-1.5, .75)) +
    scale_color_discrete(type = c('darkgreen', 'firebrick', 'gold3'),
                         guide = guide_legend(direction = 'horizontal')) +
    labs(x = 'AY Test', y = 'Effect Size', title = 'Estimated Effects Before and After Robust EBCI',
         color = 'test') +
    theme_bw() +
    theme(panel.grid.minor.x = element_blank(),
          legend.position = c(.85, .09),
          legend.background = element_rect(color = 'black'))
  
  p_all = ggplot(ebci_results |>
                   # slice(c(seq(from =    1, to = 2947, by = floor(2947/80)),
                   #           seq(from = 2948, to = 3493, by = floor( 546/110)),
                   #           seq(from = 3494, to = 6256, by = floor(5980/110)))) |>
                   mutate(x = 1:n()), 
                 aes(x = x)) +
    # geom_point(aes(y = lower_ci) )+
    # geom_point(aes(y = upper_ci)) +
    geom_point(aes(y = -1.5, color = factor(test, levels = c('negative', 'positive', 'discovery')))) +
    geom_hline(aes(yintercept = 0)) +
    geom_segment(aes(y = lower_ci, yend = upper_ci), 
                 lineend = 'square', linewidth = 1, alpha = .7, color = 'deepskyblue4') +
    geom_segment(aes(y = unshrunk_estimate, yend = .99 * shrunk_estimate + .01 * unshrunk_estimate), 
                 lineend = 'square', linejoin = 'bevel', arrow = arrow(length = unit(0.2,"cm")), linewidth = .3, alpha = .7, color = 'deepskyblue2') +
    geom_point(aes(y = unshrunk_estimate), color = 'deepskyblue2') +
    geom_point(aes(y = shrunk_estimate), color = 'deepskyblue4') +
    scale_x_continuous(expand = c(0, 0), breaks = seq(0, 10000, by = 500), 
                       # limits = c(1, 90)) +
                       # limits = c(65, 175)) +
                       # limits = c(95, 400)) +
                       # limits = c(1, 400)) +
    ) +
    scale_y_continuous(expand = c(0.025, 0), limits = c(-1.5, .75)) +
    scale_color_discrete(type = c('darkgreen', 'firebrick', 'gold3'),
                         guide = guide_legend(direction = 'horizontal')) +
    labs(x = 'AY Test', y = 'Effect Size', title = 'Estimated Effects Before and After Robust EBCI',
         color = 'test') +
    theme_bw() +
    theme(panel.grid.minor.x = element_blank(),
          legend.position = c(.85, .09),
          legend.background = element_rect(color = 'black'))
  
  return(list(p_subsampled=p_subsampled,
              p_all=p_all))
}





#' Test the formula of ebci function in order to be able to shrink towards a 
#' specific point
#' 
#' myTestFormula = 'estimate ~ 0'
#' myTestFormula = 'estimate ~ 0 + shrinkage_point'
#' myTestFormula = 'estimate ~ 0 + offset(shrinkage_point)'
#' myTestFormula = 'estimate - shrinkage_point ~ 0'
#' myTestFormula = 'estimate_centered ~ 0'
#' 
#' form1 = test_ebci_formula(myTestFormula = 'estimate ~ 0'                          , N=20)
#' form2 = test_ebci_formula(myTestFormula = 'estimate ~ 0 + shrinkage_point'        , N=20)
#' form3 = test_ebci_formula(myTestFormula = 'estimate ~ 0 + offset(shrinkage_point)', N=20)
#' form4 = test_ebci_formula(myTestFormula = 'estimate - shrinkage_point ~ 0'        , N=20)
#' form5 = test_ebci_formula(myTestFormula = 'estimate + shrinkage_point ~ 0'        , N=20)
#' form6 = test_ebci_formula(myTestFormula = 'estimate_centered ~ 0'                 , N=20)
#'
#' if(F){
#' pdf('../plots/testEBCIformula.pdf', height = 8, width = 12)
#' gridExtra::grid.arrange(form1$p_shrunk,
#'                        form2$p_shrunk,
#'                        form3$p_shrunk,
#'                        form4$p_shrunk,
#'                        form5$p_shrunk,
#'                        form6$p_shrunk,
#'                        form1$p_original,
#'                        layout_matrix = rbind(c(7, 1, 2),
#'                                              c(7, 3, 4),
#'                                              c(7, 5, 6)),
#'                        widths = c(1, 1, 1))
#' dev.off()
#' }
test_ebci_formula <- function(myTestFormula, N=100) {
  # no arg... just detect
  # assertthat::assert_that(  (!grepl('estimate_centered', myTestFormula)) || (center), 
  #                           msg='if formula includes estimate_centered, then center arg must be TRUE')
  
  center = grepl('estimate_centered', myTestFormula)
  center_minus = grepl('.*-.*~.*', myTestFormula) # e.g. 'estimate - shrinkage_point ~ 0'
  center_plus  = grepl('.*-.*~.*', myTestFormula)
  
  set.seed(12345)
  fake_effects = data.frame(estimate = rnorm(n=N, mean = 0, sd=1),
                            se       = (1:N),
                            shrinkage_point = 1:N - floor(N/2)) |>
    mutate(estimate_centered = estimate - shrinkage_point)
  
  
  test_ebci = ebci(formula = myTestFormula, 
                   data = fake_effects, 
                   se = se, weights = 1/se^2, alpha = .1)
  
  
  
  test_ebci_results = fake_effects  |> 
    select(shrinkage_point,
           unshrunk_estimate    = estimate,
           unshrunk_estimate_se = se)
  if(center || center_minus) { # centered, so add back shrinkage point
    test_ebci_results = test_ebci_results |>
      mutate(shrunk_estimate = test_ebci$df$th_eb + fake_effects$shrinkage_point)
  } else if(center_plus) { # centered by adding, so subtract
    test_ebci_results = test_ebci_results |>
      mutate(shrunk_estimate = test_ebci$df$th_eb - fake_effects$shrinkage_point)
  } else { # uncentered, so don't re add
    test_ebci_results = test_ebci_results |>
      mutate(shrunk_estimate = test_ebci$df$th_eb)
  }
  
  test_ebci_results =  test_ebci_results |> 
    mutate(lower_ci = shrunk_estimate - test_ebci$df$len_eb,
           upper_ci = shrunk_estimate + test_ebci$df$len_eb,
           unshrunk_upper_ci = unshrunk_estimate + qnorm(.05)*unshrunk_estimate_se,
           unshrunk_lower_ci = unshrunk_estimate - qnorm(.05)*unshrunk_estimate_se)
  p_original = 
    ggplot(test_ebci_results |>
             mutate(x = 1:n()), 
           aes(x = x)) +
    geom_hline(aes(yintercept = mean(unshrunk_estimate)), color = 'gray20') +
    geom_line(aes(y = shrinkage_point)) +
    # geom_segment(aes(y = lower_ci, yend = upper_ci),
    #              lineend = 'square', linewidth = .5, alpha = .7, color = 'deepskyblue4') +
    geom_segment(aes(y = unshrunk_lower_ci, yend = unshrunk_upper_ci),
                 lineend = 'square', linewidth = .5, alpha = .7, color = 'deepskyblue2') +
    # geom_segment(aes(y = unshrunk_estimate, yend = .99 * shrunk_estimate + .01 * unshrunk_estimate),
    #              lineend = 'square', linejoin = 'bevel', arrow = arrow(length = unit(0.2,"cm")), linewidth = .3, alpha = .7, color = 'deepskyblue2') +
    geom_point(aes(y = unshrunk_estimate), color = 'deepskyblue2') +
    # geom_point(aes(y = shrunk_estimate), color = 'deepskyblue4') +
    labs(x = 'AY Test', y = 'Effect Size', title = paste0('Original Estimate +- SE*1.64')) +
    theme_bw() +
    theme(panel.grid.minor.x = element_blank(),
          legend.position = c(.85, .09),
          legend.background = element_rect(color = 'black'))
  p_shrunk = ggplot(test_ebci_results |>
                      mutate(x = 1:n()), 
                    aes(x = x)) +
    geom_hline(aes(yintercept = mean(unshrunk_estimate)), color = 'gray20') +
    geom_line(aes(y = shrinkage_point)) +
    geom_segment(aes(y = unshrunk_estimate, yend = .99 * shrunk_estimate + .01 * unshrunk_estimate),
                 lineend = 'square', linejoin = 'bevel', arrow = arrow(length = unit(0.2,"cm")), linewidth = .3, alpha = .7, color = 'deepskyblue2') +
    geom_point(aes(y = unshrunk_estimate), color = 'deepskyblue2') +
    geom_point(aes(y = shrunk_estimate), color = 'deepskyblue4') +
    labs(x = 'AY Test', y = 'Effect Size', title = paste0('', myTestFormula),
         color = 'test type') +
    theme_bw() +
    theme(panel.grid.minor.x = element_blank(),
          legend.position = c(.85, .09),
          legend.background = element_rect(color = 'black'))
  
  return(list(p_original=p_original,
              p_shrunk=p_shrunk,
              df=fake_effects))
}




#' Plot before and after EBCI shrinkage
#'
#' @param shrunk_df (dataframe) of unshrunk, shrunk, and shrunk CIs
#' @param sample_idx (vector) of subsampled indices for plotting
#' @param title (string) title of the plot
#' @param legend_position (vector) of length 2, for legend.position arg of plot
plot_shrinkage <- function(shrunk_df,
                           sample_idx=NULL,
                           title='Shrinkage',
                           legend_position = c(.85, .09)) {
  
  
  # shrunk_df = ebci_results[[method]][[cl_method]]
  # sample_idx = c(seq(from =    1, to = 2947, by = floor(2947/50)), # /#tests to display
  #                seq(from = 2948, to = 3493, by = floor( 546/50)),
  #                seq(from = 3494, to = 6256, by = floor(5980/50)))
  # title = 'asdf'
  
  if(is.null(sample_idx)) {
    sample_idx = 1:nrow(shrunk_df)
  } else if(sample_idx == 'default subsample') {
    # sample_idx = c(seq(from =    1, to = 2947, by = floor(2947/50)), # /#tests to display
    #                seq(from = 2948, to = 3493, by = floor( 546/50)),
    #                seq(from = 3494, to = 6256, by = floor(5980/50)))
    
    # take a sample of 50 from each test type
    sample_idx = shrunk_df |> mutate(idx = 1:n()) |>
      select(idx, test) |>
      group_by(test) |>
      reframe(sample_idx = sample(idx, size = 50)) |> 
      pull(sample_idx) |> 
      sort()
  }
  
  
  p = ggplot(shrunk_df |>
               # slice(seq(from = 1, to = n(), by = floor(n()/400))) |>
               slice(sample_idx) |>
               mutate(test = factor(test, levels = c('negative', 'positive', 'discovery')),
                      x = 1:n()),
             aes(x = x)) +
    # geom_point(aes(y = lower_ci) )+
    # geom_point(aes(y = upper_ci)) +
    geom_point(aes(y = -1.5, color = test)) +
    geom_hline(aes(yintercept = 0)) +
    geom_segment(aes(y = lower_ci, yend = upper_ci),
                 lineend = 'square', linewidth = 1, alpha = .7, color = 'deepskyblue4') +
    geom_segment(aes(y = unshrunk_tstat, yend = .99 * shrunk_tstat + .01 * unshrunk_tstat),
                 lineend = 'square', linejoin = 'bevel', arrow = arrow(length = unit(0.2,"cm")),
                 linewidth = .3, alpha = .7, color = 'deepskyblue2') +
    geom_point(aes(y = unshrunk_tstat), color = 'deepskyblue2') +
    geom_point(aes(y = shrunk_tstat), color = 'deepskyblue4') +
    scale_x_continuous(expand = c(0, 0), breaks = seq(0, 1000, by = 50),
                       # limits = c(1, 90)) +
                       # limits = c(65, 175)) +
                       # limits = c(95, 400)) +
                       # limits = c(1, 400)) +
    ) +
    # scale_y_continuous(expand = c(0.025, 0), limits = c(-1.5, .75)) +
    scale_color_discrete(type = c('darkgreen', 'firebrick', 'gold3'),
                         guide = guide_legend(direction = 'horizontal')) +
    labs(x = 'AY Test', y = '...',
         title = title,
         subtitle = 'Before and After Robust EBCI',
         color = 'test type') +
    theme_bw() +
    theme(panel.grid.minor.x = element_blank(),
          legend.position = 'inside',
          legend.position.inside = legend_position,
          legend.background = element_rect(color = 'black'))
  
  return(p)
}


#' Plot before and after EBCI shrinkage (v2 has different inputs (more flexible))
#'
#' @param shrunk_df (dataframe) of unshrunk, shrunk, and shrunk CIs
#'    dataframe should have the following specified columns for plotting +
#'    'test', 
#' @param unshrunk_colname (string) column name of the unshrunk value
#' @param shrunk_colname (string) column name of the shrunk value
#' @param CIlower_colname (string) column name of the CI lower value
#' @param CIupper_colname (string) column name of the CI upper value
#' @param sample_idx (vector) of subsampled indices for plotting (optional)
#' @param title (string) title of the plot (optional)
#' @param legend_position (vector) of length 2, for legend.position arg of plot (optional)
plot_shrinkage2 <- function(shrunk_df,
                            unshrunk_colname,
                            shrunk_colname,
                            CIlower_colname,
                            CIupper_colname,
                            sample_idx=NULL,
                            title='Shrinkage',
                            legend_position = c(.85, .09)) {
  
  
  # shrunk_df = ebci_results[[method]][[cl_method]]
  # sample_idx = c(seq(from =    1, to = 2947, by = floor(2947/50)), # /#tests to display
  #                seq(from = 2948, to = 3493, by = floor( 546/50)),
  #                seq(from = 3494, to = 6256, by = floor(5980/50)))
  # title = 'asdf'
  
  if(is.null(sample_idx)) {
    sample_idx = 1:nrow(shrunk_df)
  } else if(sample_idx == 'default subsample') {
    # sample_idx = c(seq(from =    1, to = 2947, by = floor(2947/50)), # /#tests to display
    #                seq(from = 2948, to = 3493, by = floor( 546/50)),
    #                seq(from = 3494, to = 6256, by = floor(5980/50)))
    
    # take a sample of 50 from each test type
    set.seed(12345)
    sample_idx = shrunk_df |> mutate(idx = 1:n()) |>
      select(idx, test) |>
      group_by(test) |>
      reframe(sample_idx = sample(idx, size = 50)) |>
      pull(sample_idx) |>
      sort()
  }
  
  # either parse during plotting (better, but harder)  or  create ne cols ith knon names
  
  
  p = ggplot(shrunk_df |>
               # slice(seq(from = 1, to = n(), by = floor(n()/400))) |>
               slice(sample_idx) |>
               mutate(test = factor(test, levels = c('negative', 'positive', 'discovery')),
                      x = 1:n()),
             aes(x = x)) +
    # geom_point(aes(y = lower_ci) )+
    # geom_point(aes(y = upper_ci)) +
    geom_point(aes(y = -1.5, color = test)) +
    geom_hline(aes(yintercept = 0)) +
    geom_segment(aes(y    = eval(parse(text = CIlower_colname)), 
                     yend = eval(parse(text = CIupper_colname))),
                 lineend = 'square', linewidth = 1, alpha = .7, color = 'deepskyblue4') +
    geom_segment(aes(y = eval(parse(text = unshrunk_colname)), 
                     yend = .99 * eval(parse(text = shrunk_colname)) + 
                       .01 * eval(parse(text = unshrunk_colname))),
                 lineend = 'square', linejoin = 'bevel', arrow = arrow(length = unit(0.2,"cm")),
                 linewidth = .3, alpha = .7, color = 'deepskyblue2') +
    geom_point(aes(y = eval(parse(text = unshrunk_colname))), color = 'deepskyblue2') +
    geom_point(aes(y = eval(parse(text =   shrunk_colname))), color = 'deepskyblue4') +
    scale_x_continuous(expand = c(0, 0), breaks = seq(0, 1000, by = 50),
                       # limits = c(1, 90)) +
                       # limits = c(65, 175)) +
                       # limits = c(95, 400)) +
                       # limits = c(1, 400)) +
    ) +
    # scale_y_continuous(expand = c(0.025, 0), limits = c(-1.5, .75)) +
    scale_color_discrete(type = c('darkgreen', 'firebrick', 'gold3'),
                         guide = guide_legend(direction = 'horizontal')) +
    labs(x = 'AY Test', y = '...',
         title = title,
         subtitle = 'Before and After Robust EBCI',
         color = 'test type') +
    theme_bw() +
    theme(panel.grid.minor.x = element_blank(),
          legend.position = 'inside',
          legend.position.inside = legend_position,
          legend.background = element_rect(color = 'black'))
  
  return(p)
}

#' Plot hist of values by cluster group and by test type
#'
#' @param df (dataframe) with unshrunk effects (#tests x #other info)
#'        must have columns: value, test
#' @param value (string) name of the column to plot for histogram (e.g. estimate or tstat)
#' @param clustering (vector) clustering assignment
#'        (maybe TODO: allow matrix input of cluster probabilities)
#' @param title (string) title of plot
#'
#' @example
#' sceptre_effects = read.csv('../saves/sceptre_effects.csv')
#' set.seed(111)
#' plot_clustering(effects_df = sceptre_effects |> mutate(tstat = estimate / se),
#'                 value = 'tstat',
#'                 clustering =  sample(1:3, size = nrow(sceptre_effects), replace = TRUE))
plot_clustering <- function(df, value, clustering,
                            title    = 'Histogram of Unshrunk Estimates',
                            xlabel   = 'unshrunk values') {
  # kmeans_result = poisson_effect_kmeans2
  # effects_df = poisson_effects
  # clustering = ebci_results[['poisson']][['kmeans2']]$cluster
  # title = 'asdf'
  # value = 'tstat'
  
  
  plot_df = rbind(df[, c(value, 'test')] |> rename(category = test) |> mutate(facet = 'test'),
                  data.frame(df[, value, drop=FALSE],
                             category = clustering,
                             facet = 'cluster')) |>
    mutate(category = factor(category,
                             levels = c('negative', 'positive', 'discovery', seq(1:length(unique(clustering))))),
           facet = factor(facet, levels = c('test', 'cluster')))
  
  
  p = ggplot(plot_df) +
    geom_histogram(aes(x = eval(parse(text=value)), fill = category, y = after_stat(density)),
                   position = 'dodge',
                   bins = 25) +
    geom_vline(aes(xintercept = 0), color = 'black', linetype = 'solid', linewidth = .2) +
    # geom_vline(data = data.frame(xint=kmeans_result[['centroid']][1, 1], facet="kmeans") |>
    #              mutate(facet = factor(facet, levels = c('test', 'kmeans'))),
    #            aes(xintercept = xint), color = 'gray30', linetype = 'dashed') +
    # geom_vline(data = data.frame(xint=kmeans_result[['centroid']][2, 1], facet="kmeans") |>
    #              mutate(facet = factor(facet, levels = c('test', 'kmeans'))),
    #            aes(xintercept = xint), color = 'gray30', linetype = 'dashed') +
    # geom_vline(aes(xintercept = sceptre_effect_kmeans$centroid[1, 1]), color = 'gray30', linetype = 'dashed') +
    # geom_vline(aes(xintercept = sceptre_effect_kmeans$centroid[2, 1]), color = 'gray30', linetype = 'dashed') +
    labs(title    = title,
         subtitle = 'grouped by AY test type and cluster',
         x        = xlabel,
         fill     = 'Group') +
    scale_fill_discrete(type = c( 'darkgreen', 'firebrick', 'gold3',
                                  gray.colors(length(unique(clustering)), start = .1, end = .7)), # 'gray52', 'gray22'),
                        guide = guide_legend(direction = 'vertical')) +
    facet_grid(vars(facet)) +
    
    theme_bw()  +
    theme(legend.position = 'right',
          panel.grid.minor.y = element_blank())
  
  return(p)
}


# ==============================================================================
# ============================ helper fns ======================================
# ==============================================================================




one_hot_encode <- function(vec) {
  # vec = c(1, 1, 2, 2, 3, 3, 1, 1, 2)
  vals = sort(unique(vec))
  mat  = matrix(0, nrow = length(vec), ncol = max(vals)) |> `colnames<-`(vals)
  for(v in 1:length(vals)) {
    mat[which(vec == v), v] = 1
  }
  return(mat)
}



#' get the closest cluster (index of centroid) for each input value
#' 
#' @param values (vector) of numeric values
#' @param centroids (vector) of centers
#' @example get_closest_cluster(values = 1:10, centroids = c(0, 3, 5))
get_closest_cluster <- function(values, centroids) {
  # inefficient but whatever (maybe still faster than flexclust pred method)
  apply(X = sapply(X=centroids, FUN = function(cen) {abs(values - cen)}), 
        MARGIN = 1, 
        FUN = which.min)
  
}





