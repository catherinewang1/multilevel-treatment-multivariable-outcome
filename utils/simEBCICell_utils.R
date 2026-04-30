

# Utility functions specifically for simulations with EBCI 
# uses functions from other util files such as:
# - utils/shrink_matrix from utils/matrix_shrinkage.r
# - also needs some create_blocky_matrix(...) function defined somewhere
# - """""""""" display_matrix(...) """""""
#
#
# Main functions are:
#    sim_EBCI_celllevel
#    make_matrix_ebci_plots
# Helper Functions are:
#  


require(biclust)
require(dplyr)
require(ggplot2)
ggplot2::theme_set(theme_bw() + theme(plot.title = element_text(hjust = .5)))

########################################################
##  Helper FUNCTIONS (e.g. should be private functions)
########################################################


#' Helper function for: sim_EBCI_celllevel
#' uses packages MASS and reshape2
#' @param N (integer) sample size of total treated cells 
#' @param N_control (integer) sample size of non-treated/control cells
#' @param Theta (matrix) of effects P x G
#' # @param grna (list) of 2 matrices ($train and $test) showing grna assignments N x P
#' # @param constant_coef (numeric or vector of length G) constant coeff for glm samples 
#' @param cell_distns (vector) of distributions for the cells, some subset of c('pois', 'nb')
#' @param nb_size (numeric) size param for neg bin samples
#' @output counts (list) $<distribution>$<treatment group>$<split>
#' e.g.
#' counts$pois
#'            $treatment
#'                      $train
#'                      $test
#'            $control
#'                      ...
#'       $nb
#'            ...
#' @examples
h1_sim_cell_data <- function(N, N_control, pi_P, Theta, constant_coef, cell_distns, nb_size) {

  constant_coef = runif(n = G, min = .1, max = 2) # overall mean across genes (mean for untreated/control cells)


  # assign perturbation/grna treatment
  grna = list()
  grna$train =  t(rmultinom(n=N, size=1, prob=pi_P)) # n x p (#cells x #perturbations)
  grna$test  =  t(rmultinom(n=N, size=1, prob=pi_P)) 
  
  # 1. Simulate cell-level data -------------------------------------------------------------------------------------
  # mean effect levels for treatment cells (given constant_coef + effects from Theta)
  cell_mean_effects = list()
  cell_mean_effects$train = matrix(rep(constant_coef, times = N), byrow = TRUE, nrow = N) + grna$train %*% Theta
  cell_mean_effects$test  = matrix(rep(constant_coef, times = N), byrow = TRUE, nrow = N) + grna$test  %*% Theta

  # observed counts for each cell across genes follow pois or nb distribution
  counts           = list()
  for(distn_name in cell_distns) {
    counts[[distn_name]] = list()
    counts[[distn_name]]$treatment = list()
    counts[[distn_name]]$control   = list()

    counts[[distn_name]]$treatment$train = matrix(NA, nrow = N, ncol = G) # treated cells # train
    counts[[distn_name]]$treatment$test  = matrix(NA, nrow = N, ncol = G)                 # test 
    counts[[distn_name]]$control$train   = matrix(NA, nrow = N_control, ncol = G) # control cells # train
    counts[[distn_name]]$control$test    = matrix(NA, nrow = N_control, ncol = G)                 # test
    
    if(distn_name == 'pois') {
      for(j in 1:G) {      
        counts$pois$treatment$train[,j] = rpois(  n=N,                lambda = exp(cell_mean_effects$train[,j])) 
        counts$pois$treatment$test[,j]  = rpois(  n=N,                lambda = exp(cell_mean_effects$test[,j])) 

        counts$pois$control$train[,j]   = rpois(  n=N_control,                lambda = exp(constant_coef[j])) 
        counts$pois$control$test[,j]    = rpois(  n=N_control,                lambda = exp(constant_coef[j])) 
      } 
    } else if(distn_name == 'nb') {
      for(j in 1:G) {
        counts$nb$treatment$train[,j]   = rnbinom(n=N, size = nb_size,    mu = exp(cell_mean_effects$train[,j]))
        counts$nb$treatment$test[,j]    = rnbinom(n=N, size = nb_size,    mu = exp(cell_mean_effects$test[,j]))

        counts$nb$control$train[,j]     = rnbinom(n=N_control, size = nb_size,    mu = exp(constant_coef[j]))
        counts$nb$control$test[,j]      = rnbinom(n=N_control, size = nb_size,    mu = exp(constant_coef[j]))
      }
    } else {
      print("h1_sim_cell_data: bad cell_distns input, must be 'pois' and/or 'nb'")
    }   
  }



  # # 1. Simulate cell-level data -------------------------------------------------------------------------------------
  # # mean effect levels for treatment cells (given constant_coef + effects from Theta)
  # cell_mean_effects = list()
  # cell_mean_effects$train = matrix(rep(constant_coef, times = N), byrow = TRUE, nrow = N) + grna$train %*% Theta
  # cell_mean_effects$test  = matrix(rep(constant_coef, times = N), byrow = TRUE, nrow = N) + grna$test  %*% Theta
  
  
  # # observed counts for each cell across genes follow pois or nb distribution
  # counts           = list()
  # counts$pois      = list()
  # counts$nb        = list()

  # counts$pois$treatment = list()
  # counts$pois$control   = list()
  # counts$nb$treatment   = list()
  # counts$nb$control     = list()
  
  # counts$pois$treatment$train = matrix(NA, nrow = N, ncol = G) # treated cells # train
  # counts$pois$treatment$test  = matrix(NA, nrow = N, ncol = G)                 # test 
  # counts$nb$treatment$train   = matrix(NA, nrow = N, ncol = G)                 # train 
  # counts$nb$treatment$test    = matrix(NA, nrow = N, ncol = G)                 # test  
  
  # counts$pois$control$train = matrix(NA, nrow = N_control, ncol = G) # control cells # train
  # counts$pois$control$test  = matrix(NA, nrow = N_control, ncol = G)                 # test
  # counts$nb$control$train   = matrix(NA, nrow = N_control, ncol = G)                 # train 
  # counts$nb$control$test    = matrix(NA, nrow = N_control, ncol = G)                 # test 
  
  # for(j in 1:G) {
  #   counts$pois$treatment$train[,j] = rpois(  n=N,                lambda = exp(cell_mean_effects$train[,j])) 
  #   counts$pois$treatment$test[,j]  = rpois(  n=N,                lambda = exp(cell_mean_effects$test[,j])) 
  #   counts$nb$treatment$train[,j]   = rnbinom(n=N, size = nb_size,    mu = exp(cell_mean_effects$train[,j]))
  #   counts$nb$treatment$test[,j]    = rnbinom(n=N, size = nb_size,    mu = exp(cell_mean_effects$test[,j]))
    
  #   counts$pois$control$train[,j]   = rpois(  n=N_control,                lambda = exp(constant_coef[j])) 
  #   counts$pois$control$test[,j]    = rpois(  n=N_control,                lambda = exp(constant_coef[j])) 
  #   counts$nb$control$train[,j]     = rnbinom(n=N_control, size = nb_size,    mu = exp(constant_coef[j]))
  #   counts$nb$control$test[,j]      = rnbinom(n=N_control, size = nb_size,    mu = exp(constant_coef[j]))
    
  # }

  return(list(counts=counts, grna=grna))
}



#' Helper function for: sim_EBCI_celllevel
#' @param counts (list) of simulated cell counts (output $counts of h1_sim_cell_counts)
#' @param grna (list) of simulated cell grna/perturbation/treatment assignments (output $grna of h1_sim_cell_counts)
#' @param Theta 
#' @param Theta_rownmames (vector) of characters for rowname assignments (length P=#perturbations)
#' @param Theta_colnmames (vector) of characters for colname assignments (length G=#genes)
#' @output 
#' list(est_effects_matrices=est_effects_matrices, 
#'              est_se_matrices=est_se_matrices,
#'            allcells_results=allcells_results)'
#' @examples
h2_est_effects <- function(counts, grna, Theta, Theta_rownames, Theta_colnames) {
    # am i stupid? we could either do: 
  # (1) gene g ~ grna p for each (g, p) or
  # (2) gene g ~ \vec{grna} for each g (treat perturbations as levels in a categorical variable)
  # these should be the same (as long as there aren't other variables, eg confounders) because
  # the only variables are the treatment
  # (but when there are other variables, (2) is better bc you can then use a large sample size for the
  #  estimation of those coefficients. But, I think SCEPTRE does it one by one (1), so i will do it
  #  this way in the simulation)
  # (i mean, if you think that the confounders effect the outcome differently for each population,
  #  cells w pert 1 vs cells w pert 2, then you should do (1). But that seems maybe unlikely, and it
  #  would probably be helpful to be able to pool information from the larger sample size of cells)
  
  N_control = nrow(counts$pois$control$train)
  
  G = ncol(counts$pois$control$train)
  P = ncol(grna$train)

  # if(is.null(Theta_rownames)) {
  #   Theta_rownames = 1:P
  # } 
  # if(is.null(Theta_colnames)) {
  #   Theta_colnames = 1:G
  # }
  

 

  # Do 'naively' at first: for each (i,j) pair, perform glm fit 
  # t0 = Sys.time()
  est_effects_matrices = list()
  est_se_matrices      = list()
  for(est_method in names(counts)) {
  # for(est_method in c('pois', 'nb')) {
    est_effects_matrices[[est_method]] = list()
    est_se_matrices     [[est_method]] = list()
    for(split in c('train', 'test', 'all')) {
      est_effects_matrices[[est_method]][[split]] = matrix(NA, nrow = P, ncol = G, dimnames = list(Theta_rownames, Theta_colnames))  
      est_se_matrices     [[est_method]][[split]] = matrix(NA, nrow = P, ncol = G, dimnames = list(Theta_rownames, Theta_colnames))  
      
      # estimate for each perturbation,gene
      for(i in 1:P) {
        for(j in 1:G) {
          # construct dataframe 
          if(split == 'all') {
            pert_i_idx_train = which(grna[['train']][,i] == 1)
            cur_dfij_train = data.frame(Y = c(counts[[est_method]][['treatment']][['train']][pert_i_idx_train,j], 
                                              counts[[est_method]][['control']][['train']][,j]),
                                        A = c(rep(1, times = length(pert_i_idx_train)), rep(0, times = N_control)))

            pert_i_idx_test = which(grna[['test']][,i] == 1)
            cur_dfij_test = data.frame(Y = c(counts[[est_method]][['treatment']][['test']][pert_i_idx_test,j], 
                                             counts[[est_method]][['control']][['test']][,j]),
                                       A = c(rep(1, times = length(pert_i_idx_test)), rep(0, times = N_control)))
            cur_dfij = rbind(cur_dfij_train, cur_dfij_test)
            rm(pert_i_idx_train, pert_i_idx_test, cur_dfij_train, cur_dfij_test)
          } else { # train or test
            pert_i_idx = which(grna[[split]][,i] == 1)
            cur_dfij = data.frame(Y = c(counts[[est_method]][['treatment']][[split]][pert_i_idx,j], 
                                             counts[[est_method]][['control']][[split]][,j]),
                                  A = c(rep(1, times = length(pert_i_idx)), rep(0, times = N_control)))
            rm(pert_i_idx)
          }
          
          # glm estimation
          if(est_method == 'pois') {
            glm_fit =       glm(   formula = 'Y ~ 1 + A', data = cur_dfij, family = poisson())
          } else if(est_method == 'nb') {
            glm_fit = MASS::glm.nb(formula = 'Y ~ 1 + A', data = cur_dfij)
          }
          
          est_effects_matrices[[est_method]][[split]][i,j] = as.numeric(glm_fit$coefficients['A'])
          est_se_matrices     [[est_method]][[split]][i,j] = summary(glm_fit)$coefficients['A', 'Std. Error']
          
        }
      }
    }
  }
  # t1 = Sys.time(); print(t1 - t0)
  



  # combine results for 'all' using all cells (standard method)
  allcells_results = list()
  for(est_method in names(counts)) {
  # for(est_method in c('pois', 'nb')) {
    cur_df = 
      merge(reshape2::melt(est_effects_matrices[[est_method]][['all'  ]], value.name = 'estimate', varnames = c('grna', 'gene')),
            merge(reshape2::melt(est_effects_matrices[[est_method]][['all'  ]] - qnorm(1 - ALPHA/2) * est_se_matrices[[est_method]][['all'  ]], 
                     value.name = 'lower_ci', varnames = c('grna', 'gene')), 
                  reshape2::melt(est_effects_matrices[[est_method]][['all'  ]] + qnorm(1 - ALPHA/2) * est_se_matrices[[est_method]][['all'  ]], 
                     value.name = 'upper_ci', varnames = c('grna', 'gene')), 
                  by = c('grna', 'gene')), 
            by = c('grna', 'gene'))
    cur_df = merge(cur_df, 
                   reshape2::melt(est_se_matrices[[est_method]][['all'  ]], 
                                  value.name = 'se', varnames = c('grna', 'gene')),
                    by = c('grna', 'gene'))
    cur_df = merge(cur_df, 
                   reshape2::melt(Theta, 
                                  value.name = 'true_theta', varnames = c('grna', 'gene')),
                    by = c('grna', 'gene'))
    cur_df = merge(merge(cur_df, 
                         data.frame(y_idx = 1:length(Theta_rownames), grna = Theta_rownames),
                         by = 'grna'), 
                   data.frame(x_idx = 1:length(Theta_colnames), gene = Theta_colnames),
                   by = 'gene')
    allcells_results[[est_method]] = cur_df
    rm(cur_df)
  }

  


  return(list(est_effects_matrices=est_effects_matrices, 
              est_se_matrices=est_se_matrices,
              allcells_results=allcells_results))
}




#' Helper function for: sim_EBCI_celllevel
#' Helper function for h3_approximate_matrices
#' requires function: create_matrix_completion_softImpute from utils/matrixPrior_utils.R
#' @param mat (matrix) matrix to take approx of
#' @param method (character) name of the method to take approx of
#'   one of these: 
#'    'matcomp_linearreg', 'matcomp_softImpute', 'matdecomp_svd', 'matdecomp_sparsesvd'
#' @param ranks (vector) of integer ranks
#' @param Theta_rownmames (vector) of characters for rowname assignments (length P=#perturbations)
#' @param Theta_colnmames (vector) of characters for colname assignments (length G=#genes)
#' @output mat_approx_res 
#'      - a matrix if no ranks involved
#'      - a list where each < >[[r]] is a matrix for rank r
h3_1_approximate_matrix <- function(mat, method, ranks, Theta_rownames, Theta_colnames) {
  mat_approx_res = list()
  n = nrow(mat); m = ncol(mat)

  # === matrix completion methods (linear reg, softImpute)
  if(method == 'matcomp_linearreg') {
    # --- Linear Reg
    mat_approx_res = matrix(NA, nrow = n, ncol = m, dimnames = list(Theta_rownames, Theta_colnames))
    for(i in 1:n) {
      for(j in 1:m) {
        mat_approx_res[i,j] = linearreg(mat, i, j)
      }
    }
  } else if(method == 'matcomp_softImpute') {
    # --- softImpute
    mat_approx_res = list()
    # estimate_matapprox_plots[[est_method]][[split]][['matcomp_softImpute']] = list()
    for(r in ranks) {
      mat_approx_res[[r]] = matrix(NA, nrow = n, ncol = m, dimnames = list(Theta_rownames, Theta_colnames))
      matcomp_fun = create_matrix_completion_softImpute(X = mat, rank_max = r, type = 'als') 
      for(i in 1:n) {
        for(j in 1:m) {
          mat_approx_res[[r]][i,j] = matcomp_fun(i, j)
        }
      }
    }
  } 
  # === matrix decomposition methods (svd, sparse svd, etc...)
  else if(method == 'matdecomp_svd') {
    mat_approx_res = list()
    svdres = svd(mat, nu = max(ranks), nv = max(ranks))
    
    for(r in ranks) {
      # --- SVD
      cur_matapprox = svdres$u[, 1:r, drop=FALSE] %*% diag(svdres$d[1:r], nrow=r, ncol=r) %*% t(svdres$v[, 1:r, drop=FALSE]) # take up to r eigenvectors
      row.names(cur_matapprox) = Theta_rownames; colnames(cur_matapprox) = Theta_colnames
      mat_approx_res[[r]] = cur_matapprox
      rm(cur_matapprox)      
    }
  } else if(method == 'matdecomp_sparsesvd') {
    mat_approx_res = list()
    cv.out = PMA::SPC.cv(mat, sumabsvs = seq(1.2, min(5, sqrt(n), sqrt(m)), len = 10), trace=FALSE)
    pmd_res = PMA::SPC(mat, sumabsv=cv.out$bestsumabsv, K = max(ranks), trace=FALSE)
    for(r in ranks) {
      # --- Sparse SVD
      cur_matapprox = pmd_res$u[, 1:r, drop=FALSE] %*% diag(pmd_res$d[1:r], nrow=r, ncol=r) %*% t(pmd_res$v[, 1:r, drop=FALSE]) # take up to r eigenvectors
      row.names(cur_matapprox) = Theta_rownames; colnames(cur_matapprox) = Theta_colnames
      mat_approx_res[[r]] = cur_matapprox
      rm(cur_matapprox)
    }
  } 
  # === Clustering methods
  else if(method == 'spectralbiclust' | method == 'spectralbiclust_threshold') {
    # with the way that this is called, the spectral biclustering will be performed 2x as many times it needs to...
    mat_approx_res = list()
    for(r in ranks) {
      # try, or else return null
      inner_biclustspectral = function() {
        spbicl_res = biclust::spectral(mat,
                   # normalization = "bistochastization", # log, irrc, bistochastization (in example, bistochastization made could not find cl??)
                   # normalization = "log", # they recommend log even w their example matrix of negative values?
                   numberOfEigenvalues = r, # here use ranks? or n_best?
                   minr=2, minc=2, withinVar=2, n_clusters = NULL, 
                   n_best = min(r, 3) # #e.vecs to which the data is projected for the final clustering step, recommended values are 2 or 3. but maybe we should do more?
                   )
        row_cl = data.frame(row_idx = row.names(mat), 
                            row_cl  = spbicl_res@info$row_labels)
        col_cl = data.frame(col_idx = colnames(mat), 
                            col_cl  = spbicl_res@info$column_labels)

        bicluster_assignment = merge(reshape2::melt(mat, varnames = c('row_idx', 'col_idx'), value.name = 'value'), 
                                     tidyr::expand_grid(row_cl, col_cl) |> dplyr::select(row_idx, col_idx, row_cl, col_cl), 
                                     by = c('row_idx', 'col_idx'))
        if(method == 'spectralbiclust') {
          biclust_mean = bicluster_assignment |> dplyr::group_by(row_cl, col_cl) |> dplyr::summarise(bicl_mean = mean(value))
          biclust_mat = reshape2::acast(merge(bicluster_assignment, biclust_mean, by = c('row_cl', 'col_cl')) |> 
                                          dplyr::select(row_idx, col_idx, bicl_mean) |> dplyr::arrange(row_idx, col_idx), 
                                        row_idx ~ col_idx, value.var="bicl_mean")
        
        } else if(method == 'spectralbiclust_threshold') {
          biclust_mean = bicluster_assignment |> dplyr::group_by(row_cl, col_cl) |> dplyr::summarise(bicl_mean = mean(value), bicl_sd = sd(value), count = dplyr::n())
          nbiclusters = nrow(biclust_mean) 
          # zero out if pval is sig dif from 0 w bonferonni correction at alpha=.05 (very rough)
          # from normalization procedure and biclustering, each bicluster should have sd around 1, I think?
          biclust_mean = biclust_mean |> dplyr::mutate(zstat = (bicl_mean - 0)/(bicl_sd/sqrt(count)), pval = pnorm(-abs(zstat)), zeroed = as.integer(pval > .05/nbiclusters), bicl_mean_zeroed = (1 - zeroed) * bicl_mean)        
          
          biclust_mat = reshape2::acast(merge(bicluster_assignment, biclust_mean, by = c('row_cl', 'col_cl')) |> 
                                        dplyr::select(row_idx, col_idx, bicl_mean_zeroed) |> dplyr::arrange(row_idx, col_idx), 
                                        row_idx ~ col_idx, value.var="bicl_mean_zeroed")
        }
        row.names(biclust_mat) = Theta_rownames; colnames(biclust_mat) = Theta_colnames
        return(biclust_mat)
      }

      biclust_mat = tryCatch(expr = {inner_biclustspectral()},
                             error = function(e) {'badfit'}
                    )
      # # Without specifying the  number of row and column clusters
      # spbicl_res = biclust::spectral(mat,
      #              # normalization = "irrc", # log, irrc, bistochastization (in example, bistochastization made could not find cl??)
      #              normalization = "log", # they recommend log even w their example matrix of negative values?
      #              numberOfEigenvalues = r, # here use ranks? or n_best?
      #              minr=2, minc=2, 
      #              withinVar=2, # this should be chosen more carefully
      #              n_clusters = NULL, 
      #              n_best = 3 # #e.vecs to which the data is projected for the final clustering step, recommended values are 2 or 3. but maybe we should do more?
      #              )
      

      mat_approx_res[[r]] = biclust_mat
    }    
  }
  # === Basic Shrinkage Points
  else if(method == 'zeros') { 
    # --- zero 0
    mat_approx_res = matrix(0,         nrow = n, ncol = m, dimnames = list(Theta_rownames, Theta_colnames))
  } else if (method == 'average') {
    # --- overall average 
    mat_approx_res = matrix(mean(mat), nrow = n, ncol = m, dimnames = list(Theta_rownames, Theta_colnames))
  } else { # bad method name input
    print('bad approx matrix method name input in function h3_1_approximate_matrix')
  }
  
  return(mat_approx_res)



  




  # === matrix decomposition methods (svd, sparse svd, etc...)

  # estimate_matapprox_plots[[est_method]][[split]][['matdecomp_svd']]       = list()
  # estimate_matapprox_plots[[est_method]][[split]][['matdecomp_sparsesvd']] = list()


}

#' Helper function for: sim_EBCI_celllevel
#' 
#' @param est_effects_matrices
#' @param ranks
#' @param Theta_rownmames (vector) of characters for rowname assignments (length P=#perturbations)
#' @param Theta_colnmames (vector) of characters for colname assignments (length G=#genes)
#' @param matapprox_methods (vector) of characters e.g c('matcomp_linearreg', 'matcomp_softImpute', 'matdecomp_svd', 'matdecomp_sparsesvd', 'spectralbiclust', 'spectralbiclust_threshold', 'zeros', 'average')
#' @output estimate_matapprox
h3_approximate_matrices <- function(est_effects_matrices, ranks, Theta_rownames, Theta_colnames, matapprox_methods=NULL) {
  estimate_matapprox = list()
  
  for(distn_name in names(est_effects_matrices)) {
  # for(distn_name in c('pois', 'nb')) {
    estimate_matapprox      [[distn_name]] = list()
    
    # for(split in c('train', 'test', 'all')) {
    for(split in c('train', 'all')) { # only approx train and all split
      estimate_matapprox[[distn_name]][[split]] = list()
      cur_mat = est_effects_matrices[[distn_name]][[split]] # matrix to make approximations of
      # n = nrow(cur_mat); m = ncol(cur_mat)
      if(is.null(matapprox_methods)) {
        matapprox_methods = c('matcomp_linearreg', 'matcomp_softImpute', 'matdecomp_svd', 'matdecomp_sparsesvd', 'spectralbiclust', 'spectralbiclust_threshold', 'zeros', 'average')
      }
      for(matapprox_method in matapprox_methods) {
        approximated_matrix = h3_1_approximate_matrix(mat = cur_mat,  method = matapprox_method, ranks = ranks,
                                                      Theta_rownames=Theta_rownames, Theta_colnames=Theta_colnames)
        # if(!is.null(approximated_matrix)) { # only add if the result is not null (e.g. the approx orks)
        #   estimate_matapprox[[distn_name]][[split]][[matapprox_method]] = approximated_matrix
        # }
        estimate_matapprox[[distn_name]][[split]][[matapprox_method]] = approximated_matrix
      }
    }
  }
  
  return(estimate_matapprox)
}



#' Helper function for: sim_EBCI_celllevel
#' uses function shrink_matrix from utils...
#' @param estimate_matapprox
#' @param est_se_matrices
#' @param estimate_matapprox
#' @param ALPHA 
#' @param Theta
#' @param Theta_rownmames (vector) of characters for rowname assignments (length P=#perturbations)
#' @param Theta_colnmames (vector) of characters for colname assignments (length G=#genes)
#' @param ranks
#' @param calc_ebci_pvals (boolean) calc ebci_pvals (takes a hile) or not (can just load in the saved ebci_object to calc pvals later)
h4_perform_ebci <- function(est_effects_matrices, est_se_matrices, estimate_matapprox, 
                            ALPHA, Theta, Theta_rownames, Theta_colnames, ranks, calc_ebci_pvals=FALSE) {
  
  # use shrink_matrix from utils/matrix_shrinkage.r
  ebci_params       = list() # estimated parameters from ebci fit
  shrinkage_results = list() # results from shrinkage
  for(est_method in names(est_effects_matrices)) { # pois or nb
    shrinkage_results[[est_method]] = list()
    ebci_params      [[est_method]] = list()
    # do 2 categories of shrinkage: sample split (train towards test), no sample split (all towards all)
    for(splittype in c('samplesplit', 'nosamplesplit')) {
      shrinkage_results[[est_method]][[splittype]] = list()
      ebci_params      [[est_method]][[splittype]] = list()

      if(splittype == 'samplesplit') {
          cur_estimates_mat = est_effects_matrices[[est_method]][['test']] # matrix of estimates to shrink
          cur_se_mat        =      est_se_matrices[[est_method]][['test']] # matrix of se's of est to shrink
      } else if(splittype == 'nosamplesplit') {
          cur_estimates_mat = est_effects_matrices[[est_method]][['all']] # matrix of estimates to shrink
          cur_se_mat        =      est_se_matrices[[est_method]][['all']] # matrix of se's of est to shrink
      }

    
      for(approx_method in names(estimate_matapprox[[est_method]][['train']])) {
        # # could test by using: is.matrix or is.list
        # is.matrix(estimate_matapprox[[est_method]][['matdecomp_sparsesvd']])
        # is.matrix(estimate_matapprox[[est_method]][['matcomp_linearreg']])
        # is.list(estimate_matapprox[[est_method]][['matdecomp_sparsesvd']])
        # is.list(estimate_matapprox[[est_method]][['matcomp_linearreg']])
        
        # if there are not ranks (eg linearreg)
        if(is.matrix(estimate_matapprox[[est_method]][['train']][[approx_method]])) {
          if(splittype == 'samplesplit') {
              cur_shrinkagepoint_mat = estimate_matapprox[[est_method]][['train']][[approx_method]] # matrix to shrink towards
          } else if(splittype == 'nosamplesplit') {
              cur_shrinkagepoint_mat = estimate_matapprox[[est_method]][['all']][[approx_method]] # matrix to shrink towards
          }

          
          
          # shrink matrix
          shrink_mat_res =  
            shrink_matrix(unshrunk_mat = cur_estimates_mat,
                          shrinkpoint_mat = cur_shrinkagepoint_mat,
                          se_mat = cur_se_mat,
                          ALPHA = ALPHA,
                          return_ebci_obj=TRUE,
                          weight_mat=NULL)
          # add true theta values
          cur_shrink_res = merge(shrink_mat_res$ebci_res, 
                                 reshape2::melt(Theta, varnames = c('grna', 'gene'), value.name = 'true_theta'),
                                 by = c('grna', 'gene'))
          # add idx for plotting later
          cur_shrink_res = merge(merge(cur_shrink_res, 
                                       data.frame(y_idx = 1:length(Theta_rownames), grna = Theta_rownames),
                                       by = 'grna'), 
                                 data.frame(x_idx = 1:length(Theta_colnames), gene = Theta_colnames),
                                 by = 'gene')
          if(calc_ebci_pvals) {
            cur_shrink_res$ebci_pvals = h4_1_ebci_pvals(shrinkage_results_df = cur_shrink_res, ebci_obj = shrink_mat_res$ebci_obj)
          }
          
          
          # save result
          shrinkage_results[[est_method]][[splittype]][[approx_method]] = cur_shrink_res
          ebci_params      [[est_method]][[splittype]][[approx_method]] = shrink_mat_res$ebci_obj[c('mu2', 'kappa', 'delta', 'alpha')]
          
        } else { # if there are ranks (eg softImpute, svd, sparse svd, ..) 
          for(r in ranks) { # if there was an error estimating some subset of ranks, then this might cause issues
            if(splittype == 'samplesplit') {
                cur_shrinkagepoint_mat = estimate_matapprox[[est_method]][['train']][[approx_method]][[r]] # matrix to shrink towards
            } else if(splittype == 'nosamplesplit') {
                cur_shrinkagepoint_mat = estimate_matapprox[[est_method]][['all']][[approx_method]][[r]] # matrix to shrink towards
            }

            if(!is.matrix(cur_shrinkagepoint_mat))  { # if mat approx as bad, skip
              next()
            }

            # shrink matrix
            shrink_mat_res = 
              shrink_matrix(unshrunk_mat = cur_estimates_mat,
                            shrinkpoint_mat = cur_shrinkagepoint_mat,
                            se_mat = cur_se_mat,
                            ALPHA = ALPHA,
                            return_ebci_obj=TRUE,
                            weight_mat=NULL)
            # add true theta values
            cur_shrink_res = merge(shrink_mat_res$ebci_res, 
                                   reshape2::melt(Theta, varnames = c('grna', 'gene'), value.name = 'true_theta'),
                                   by = c('grna', 'gene'))
            # add idx for plotting later
            cur_shrink_res = merge(merge(cur_shrink_res, 
                                         data.frame(y_idx = 1:length(Theta_rownames), grna = Theta_rownames),
                                         by = 'grna'), 
                                   data.frame(x_idx = 1:length(Theta_colnames), gene = Theta_colnames),
                                   by = 'gene')

            if(calc_ebci_pvals) {
              cur_shrink_res$ebci_pvals = h4_1_ebci_pvals(shrinkage_results_df = cur_shrink_res, ebci_obj = shrink_mat_res$ebci_obj)
            }
            

            # save result
            shrinkage_results[[est_method]][[splittype]][[approx_method]][[r]] = cur_shrink_res
            ebci_params      [[est_method]][[splittype]][[approx_method]][[r]] = shrink_mat_res$ebci_obj[c('mu2', 'kappa', 'delta', 'alpha')]
            
          }
          
        }
      }
    } 




  }

  return(list(shrinkage_results=shrinkage_results,
              ebci_params=ebci_params))
}

#' Helper function for: sim_EBCI_celllevel
#' 4 get the average pvals (inverted from ebci) 
#' @param shrinkage_results_df (dataframe) result of  shrink_matrix(...)$ebci_res a dataframe w shrinkage results
#' shouldhave columns 'shrunk_value', 'se', 'w_eb' for calculation of inverted ebci pval
#' @param ebci_obj (ebci_object) that has the parameter estimates (mu2, kappa)
h4_1_ebci_pvals <- function(shrinkage_results_df, ebci_obj) {
  average_pvals = rep(NA, times = nrow(shrinkage_results_df))

  cur_mu2   = tryCatch(expr = {ebci_obj$mu2['estimate']},
                       error = function(e) {NA})
  cur_kappa = tryCatch(expr = {ebci_obj$kappa['estimate']},
                       error = function(e) {NA})
  if(!is.numeric(cur_mu2) | !is.numeric(cur_kappa)) {return(average_pvals)} # if ebci est bad, don't calc
  
  for(i in 1:nrow(shrinkage_results_df)) {

    average_pvals[i] = tryCatch(expr = {get_ebci_pvals(thetashrunk     = shrinkage_results_df[i, 'shrunk_value'], 
                             sigma           = shrinkage_results_df[i, 'se'], 
                             web             = shrinkage_results_df[i, 'w_eb'] , 
                             mu2             = cur_mu2, 
                             kappa           = cur_kappa,
                             alpha_threshold = .001, 
                             maxiter         = 10) },
                             error = function(e) {NA}
                    )
    # average_pvals[i] = 
    #   get_ebci_pvals(thetashrunk     = shrinkage_results_df[i, 'shrunk_value'], 
    #                          sigma           = shrinkage_results_df[i, 'se'], 
    #                          web             = shrinkage_results_df[i, 'w_eb'] , 
    #                          mu2             = cur_mu2, 
    #                          kappa           = cur_kappa,
    #                          alpha_threshold = .001, 
    #                          maxiter         = 10) 
  }
  return(average_pvals)
}

#' Helper function for: sim_EBCI_celllevel
#' 5.1 some plots particular for each method
#' uses function plot_shrink_results from utils/matrix_shrinkage.r
#' 
#' @param shrinkage_results
#' @param save_folder
#' @param ranks
h5_1_plots_summary <- function(shrinkage_results, save_folder, ranks) {
  # 5.1 some plots particular for each method
  # use plot_shrink_results from utils/matrix_shrinkage.r
  
  for(est_method in names(shrinkage_results)) { # pois or nb
    for(approx_method in names(shrinkage_results[[est_method]])) {
      # # could test by using: is.matrix or is.list... should use list so it is consistent across different loops
      # is.matrix(shrinkage_results[[est_method]][['matdecomp_sparsesvd']])
      # is.matrix(shrinkage_results[[est_method]][['matcomp_linearreg']])
      # is.list(shrinkage_results[[est_method]][['matdecomp_sparsesvd']])
      # is.list(shrinkage_results[[est_method]][['matcomp_linearreg']])
      # is.data.frame(shrinkage_results[[est_method]][['matdecomp_sparsesvd']])
      # is.data.frame(shrinkage_results[[est_method]][['matcomp_linearreg']])
      
      # if there are not ranks (eg linearreg)
      if(is.data.frame(shrinkage_results[[est_method]][['samplesplit']][[approx_method]])) {
        cur_shrinkage_plots_folder = sprintf('%s%s/', save_folder, approx_method)
        dir.create(sprintf('%spoints/', cur_shrinkage_plots_folder), showWarnings = FALSE, recursive = TRUE)
        cur_shrink_res = shrinkage_results[[est_method]][['samplesplit']][[approx_method]]
        
        plot_shrink_results(shrink_df      = cur_shrink_res, 
                            plot_folder    = cur_shrinkage_plots_folder, 
                            order_rowscols = FALSE, 
                            grna_index     = NULL, 
                            gene_index     = NULL, 
                            unshrunk_ALPHA = ALPHA)
        
        
        
      } else { # if there are ranks (eg softImpute, svd, sparse svd, ..) 
        for(r in ranks) { # if there was an error estimating some subset of ranks, then this might cause issues
          cur_shrinkage_plots_folder = sprintf('%s/%s_%02.f/', save_folder, approx_method, r)
          dir.create(sprintf('%spoints/', cur_shrinkage_plots_folder), showWarnings = FALSE, recursive = TRUE)
          cur_shrink_res = shrinkage_results[[est_method]][['samplesplit']][[approx_method]][[r]]
          
          plot_shrink_results(shrink_df      = cur_shrink_res, 
                              plot_folder    = cur_shrinkage_plots_folder, 
                              order_rowscols = FALSE, 
                              grna_index     = NULL, 
                              gene_index     = NULL, 
                              unshrunk_ALPHA = ALPHA)
        }
      }
    }
  }
}

#' Helper function for: sim_EBCI_celllevel
#' 5.2 matrix plots
#' uses function make_matrix_ebci_plots
#' 
#' @param shrinkage_results
#' @param est_effects_matrices
#' @param estimate_matapprox
#' @param save_folder
#' @param allcells_results
#' @param Theta
#' @param ranks
h5_2_plots_matrix <- function(shrinkage_results, est_effects_matrices, estimate_matapprox, 
                              save_folder, allcells_results, Theta, ranks) {
  # for(est_method in c('pois', 'nb')) {
  for(est_method in names(shrinkage_results)) {
    for(r in ranks) {
      # for(splittype in c('samplesplit', 'nosamplesplit')) {
        # print(sprintf("Part 5.2 %s %s", est_method, r))
        # dir.create(sprintf('%s/%s/', save_folder, splittype))
        make_matrix_ebci_plots(est_method = est_method,   
                               chosen_rank_to_plot = r, 
                               shrinkage_results = shrinkage_results, 
                               est_effects_matrices = est_effects_matrices, 
                               estimate_matapprox = estimate_matapprox, 
                               allcells_results = allcells_results, 
                               # save_folder = sprintf('%s/%s/', save_folder, splittype),     
                               save_folder=save_folder,                         
                               Theta=Theta,
                               color_limits = c(floor(min(Theta-.75)), ceiling(max(Theta+.75))))
      # }
      
    }
  }
}

#' creates a dataframe to help with plotting, this dataframe will have
#' shrinkage results stacked together with columns
#' sim_distn, method, shrinkage_point, weight, shrunk_value, unshrunk_value, lower_ci, upper_ci, etc...
#' @param shrinkage_results
#' @param allcells_results
#' @param ranks
#' @param ALPHA
h5_0_create_plot_df <- function(shrinkage_results, allcells_results, ranks, ALPHA) {
  plot_df = NULL
    for(est_method in names(shrinkage_results)) {
      for(splittype in c('samplesplit', 'nosamplesplit')) {
        plot_df_ = shrinkage_results[[est_method]][[splittype]][['matcomp_linearreg']] |> # no shrinkage/glm, pick any (would be same)
                        dplyr::mutate(method = 'unshrunk', rank = NA) |>
                        dplyr::mutate(shrinkage_point = unshrunk_value,   # set all values to the original estimates
                                      shrunk_value    = unshrunk_value, 
                                      lower_ci        = unshrunk_value - qnorm(1 - ALPHA/2) * se,
                                      upper_ci        = unshrunk_value + qnorm(1 - ALPHA/2) * se) 
        plot_df_ = rbind(plot_df_, 
                       allcells_results[[est_method]] |> # all cells
                         dplyr::mutate(method = 'unshrunkallcells', rank = NA) |>
                         dplyr::mutate(shrinkage_point = NA,   # set all values to the original estimates
                                       weight = NA, 
                                       shrunk_value    = estimate,
                                       unshrunk_value  = estimate,
                                       w_eb            = NA,
                                       ebci_pvals      = NA)  |> 
                         dplyr::select(all_of(colnames(plot_df_))) )
      for(approx_method in names(shrinkage_results[[est_method]][[splittype]])) {
         # if there are not ranks (eg linearreg)
          if(is.data.frame(shrinkage_results[[est_method]][[splittype]][[approx_method]])) {
            plot_df_ = rbind(plot_df_, 
                            shrinkage_results[[est_method]][[splittype]][[approx_method]] |> 
                            dplyr::mutate(method = approx_method, rank = NA))
          } else { # if there are ranks (eg softImpute, svd, sparse svd, ..) 
            for(r in ranks) { # if there was an error estimating some subset of ranks, then this might cause issues
               plot_df_ = rbind(plot_df_, 
                            shrinkage_results[[est_method]][[splittype]][[approx_method]][[r]] |> 
                            dplyr::mutate(method = approx_method, rank = r))
            }
            
          }
      }
      plot_df = rbind(plot_df, plot_df_ |> dplyr::mutate(sim_distn = est_method, split_type = splittype, .before = 1))
      rm(plot_df_)
      }
    }
    
    plot_df$method = factor(plot_df$method, 
                            levels = c("unshrunkallcells" , "unshrunk", "matcomp_linearreg", "matcomp_softImpute", "matdecomp_svd", "matdecomp_sparsesvd", "spectralbiclust", "spectralbiclust_threshold", "zeros", "average"))
    

    return(plot_df)
}

#' Helper function for: sim_EBCI_celllevel
#' @param plot_df
#' @param ranks
#' @param save_folder
#' @param plot_df_is_summ (boolean)
#' @output saves plot at sprintf('%s/shrinkage_mse.pdf', save_folder)
h5_3_plots_mse <- function(plot_df, ranks, save_folder, height=NULL, width=NULL, save_ggplot=FALSE, plot_df_is_summ=FALSE) {
    if(is.null(height)){height=5} # default height=5, width=8
    if(is.null(width) ){ width=8}
    # mse
    # p_mse = ggplot(plot_df |> 
    #          group_by(sim_distn, method, rank) |> 
    #          summarize(mse = mean((shrunk_value - true_theta)^2), .groups = 'drop') |> 
    #          arrange(sim_distn, method, rank) |> 
    #          mutate(methodrank = factor(paste0(method, rank), 
    #                                 levels = c("unshrunkallcellsNA",
    #                                            "unshrunkNA", 
    #                                            "matcomp_linearregNA", 
    #                                            paste0("matcomp_softImpute", ranks), 
    #                                            paste0("matdecomp_svd", ranks),
    #                                            paste0("matdecomp_sparsesvd", ranks)))), 
    #        aes(x = methodrank, y = mse)) +
    #   geom_col(color = 'black', fill = 'gray') +
    #   facet_grid(rows = vars(sim_distn), scales = 'free_y') +
    #   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))

    # p_mse = ggplot(plot_df |> 
    #    filter(method != 'matcomp_linearreg') |> # exclude this... this performs badly
    #    group_by(sim_distn, split_type, method, rank) |> 
    #    summarize(mse = mean((shrunk_value - true_theta)^2)) |> 
    #    arrange(sim_distn, split_type, method, rank) |> 
    #    mutate(methodrank = factor(paste0(method, rank), 
    #                           levels = c("unshrunkallcellsNA",
    #                                      "unshrunkNA", 
    #                                      "matcomp_linearregNA", 
    #                                      paste0("matcomp_softImpute", ranks), 
    #                                      paste0("matdecomp_svd", ranks),
    #                                      paste0("matdecomp_sparsesvd", ranks), 
    #                                      paste0("spectralbiclust", ranks), 
    #                                      paste0("spectralbiclust_threshold", ranks), 
    #                                      "zerosNA", 
    #                                      "averageNA"))), 
    #   aes(x = methodrank, y = mse)) +
    #   geom_col(color = 'black', fill = 'gray') +
    #   facet_grid(rows = vars(sim_distn), cols = vars(split_type), scales = 'free_y') +
    #   labs(title = 'MSE', x = 'method + rank', y = 'mse') +
    #   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5), panel.grid.major.x = element_blank())
    
    # ggsave(filename = sprintf('%s/shrinkage_mse.pdf', save_folder), plot = p_mse, width = 6, height = 6) 

    # p_mse_nonzero = ggplot(plot_df |> 
    #    filter(method != 'matcomp_linearreg') |> # exclude this... this performs badly 
    #    filter(true_theta != 0) |>
    #    group_by(sim_distn, split_type, method, rank) |> 
    #    summarize(mse = mean((shrunk_value - true_theta)^2)) |> 
    #    arrange(sim_distn, split_type, method, rank) |> 
    #    mutate(methodrank = factor(paste0(method, rank), 
    #                           levels = c("unshrunkallcellsNA",
    #                                      "unshrunkNA", 
    #                                      "matcomp_linearregNA", 
    #                                      paste0("matcomp_softImpute", ranks), 
    #                                      paste0("matdecomp_svd", ranks),
    #                                      paste0("matdecomp_sparsesvd", ranks), 
    #                                      paste0("spectralbiclust", ranks), 
    #                                      paste0("spectralbiclust_threshold", ranks), 
    #                                      "zerosNA", 
    #                                      "averageNA"))), 
    #   aes(x = methodrank, y = mse)) +
    #   geom_col(color = 'black', fill = 'gray') +
    #   facet_grid(rows = vars(sim_distn), cols = vars(split_type), scales = 'free_y') +
    #   labs(title = 'MSE', x = 'method + rank', y = 'mse') +
    #   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5), panel.grid.major.x = element_blank())
    
    # ggsave(filename = sprintf('%s/shrinkage_mse_nonzero.pdf', save_folder), plot = p_mse_nonzero, width = 6, height = 6) 
    # return(p_mse)





    # === prep nice labels and colors ===


    # === colors
    # distinct_colors = paletteer::paletteer_d("colorBlindness::paletteMartin")  # library(paletteer)
    # distinct_colors[c(2, 4, 6, 7, 12, 14)]
    # methodrank_colors =  
    #   c(colorRampPalette(c("white", distinct_colors[2]))(2 + 1)[1:2+1], # (original)
    #     colorRampPalette(c("white", distinct_colors[4]))(length(ranks) + 1)[1:length(ranks)+1], # softImpute
    #     colorRampPalette(c("white", distinct_colors[6]))(length(ranks) + 1)[1:length(ranks)+1], # SVD
    #     colorRampPalette(c("white", distinct_colors[7]))(length(ranks) + 1)[1:length(ranks)+1], # Sparse SVD
    #     # colorRampPalette(c("white", distinct_colors[12]))(length(ranks) + 1)[1:length(ranks)+1], # biclustering
    #     colorRampPalette(c("white", distinct_colors[14]))(2 + 1)[1:2+1] # zeros/avg
    #   )
    
    methodrank_colors = create_color_pallete_nicenames(ranks = ranks)

    
    # === labels
    # method_nicenames = c("unshrunkallcells"="Original (Full Dataset)",
    #                      "unshrunk"="Original", 
    #                      "matcomp_linearreg"="Linear Regression", 
    #                      "matcomp_softImpute"="softImpute", 
    #                      "matdecomp_svd"="SVD",
    #                      "matdecomp_sparsesvd"="Sparse SVD", 
    #                      "spectralbiclust"="Spectral Biclustering", 
    #                      "spectralbiclust_threshold"="Spectral Biclustering w/ Thresholding", 
    #                      "zeros"="Zeros", 
    #                      "average"="Average")

    # rank_nicenames = paste0('(rank=', 1:100, ')')



    # methodrank_nicenames <- function(method_name, rank_) {
    #   # method_name = 'zeros'
    #   # rank_ = 2
    #   # rank_ = NA
    #   # methodrank_nicenames(method_name = plot_df_summ$method[1], rank_ = plot_df_summ$rank[1]) |> unname()
    #   # mapply(FUN = methodrank_nicenames, method_name =  plot_df_summ$method[1:4], rank_ = plot_df_summ$rank[1:4]) |> unname()

      
    #   if(is.na(rank_)) {
    #     return(method_nicenames[method_name] |> unname())
    #   } else {
    #     return(paste0(method_nicenames[method_name], ' ', rank_nicenames[rank_]) |> unname())
    #   }
      
    # }


    # make the ordering for a various number of methodrank levels (even if not used): first method according to method_nicenames, then rank NA, 1, 2, ...
    methodrank_nicenames_order = c()
    for(cur_method in names(method_nicenames)) { # requires declaration of this list/vector (this is defined later in this file, right before the function methodrank_nicenames is defined)
      for(cur_rank in c(NA, ranks)) {
        methodrank_nicenames_order = c(methodrank_nicenames_order, methodrank_nicenames(method_name = cur_method, rank_ = cur_rank) |> unname())
      }
    }

    # sim_distn_nicename = c('pois'= 'Poisson', 'nb'='Negative Binomial') # could use these explicitly, or put labels when making the factor levels/labels
    # split_type_nicename = c('nosamplesplit'='Full Dataset', 'samplesplit'='Sample Split')


    # === MSE overall ===
    if(plot_df_is_summ) {
      plot_df_summ = plot_df |> 
        filter(method != 'matcomp_linearreg') |> # exclude this... this performs badly
        group_by(sim_distn, split_type, method, rank) |> 
        summarize(mse = mean(mse), .groups = 'drop') |> 
        mutate(sim_distn  = factor(sim_distn,  levels = c('pois', 'nb'),                   labels = c('Poisson', 'Negative Binomial')),
               split_type = factor(split_type, levels = c('nosamplesplit', 'samplesplit'), labels = c('Full Dataset', 'Sample Split')))
    } else {
      plot_df_summ = plot_df |> 
          filter(method != 'matcomp_linearreg') |> # exclude this... this performs badly
          group_by(sim_distn, split_type, method, rank) |> 
          summarize(mse = mean((shrunk_value - true_theta)^2)) |> 
          mutate(sim_distn  = factor(sim_distn,  levels = c('pois', 'nb'),                   labels = c('Poisson', 'Negative Binomial')),
                 split_type = factor(split_type, levels = c('nosamplesplit', 'samplesplit'), labels = c('Full Dataset', 'Sample Split')))
    }
   

    plot_df_summ$methodrank = mapply(FUN = methodrank_nicenames, method_name =  plot_df_summ$method, rank_ = plot_df_summ$rank) |> unname()
    plot_df_summ$methodrank = factor(plot_df_summ$methodrank, levels = methodrank_nicenames_order)
    
    # take the mse from unshrunk on Full Dataset
    # original_all_mse = plot_df_summ |> filter(method == 'unshrunkallcells')
    original_all_mse = plot_df_summ |> filter(method == 'unshrunk' & split_type == 'Full Dataset')
    original_all_mse = rbind(original_all_mse, original_all_mse |> mutate(split_type = 'Sample Split')) # repeat but change split type for plotting
    
    p_mse = ggplot() +
      geom_col(data = plot_df_summ, aes(x = methodrank, y = mse, fill = methodrank), color = 'black', alpha = 1) + # MSEs
      geom_hline(data = original_all_mse, aes(yintercept = mse), color = '#DB1A1A', linewidth = .7, alpha = .6) +  # MSEs for 'naive'/original method
      facet_grid(rows = vars(sim_distn), cols = vars(split_type), scales = 'free_y', 
                 # labeller = as_labeller(c(sim_distn_nicename, split_type_nicename))
                 labeller = label_value
                 ) +
      scale_fill_discrete(palette = methodrank_colors[names(methodrank_colors) %in% plot_df_summ$methodrank]) +
      labs(title = 'MSE', x = 'method + rank', y = 'mse', fill = 'Method') +
      scale_y_continuous(expand = expansion(mult = c(0, .05))) +
      theme(# axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5), 
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
            panel.grid.major.x = element_blank(), strip.background = element_rect(fill = NA), 
            axis.ticks.x = element_blank())

    ggsave(filename = sprintf('%s/mse_vert.pdf', save_folder), plot = p_mse, width = width, height = height) 
    if(save_ggplot) {
      saveRDS(p_mse, sprintf('%s/ggplot_mse_vert.rds', save_folder)) # save the ggplot objects (to adjust later as needed)
    }
    

    p_mse = ggplot() +
      geom_col(data = plot_df_summ, aes(y = methodrank, x = mse, fill = methodrank), color = 'black', alpha = 1) + # MSEs
      geom_vline(data = original_all_mse, aes(xintercept = mse), color = '#DB1A1A', linewidth = .7, alpha = .6) +  # MSEs for 'naive'/original method
      facet_grid(cols = vars(sim_distn), rows = vars(split_type), scales = 'free', 
                 # labeller = as_labeller(c(sim_distn_nicename, split_type_nicename))
                 labeller = label_value
                 ) +
      scale_fill_discrete(palette = methodrank_colors[names(methodrank_colors) %in% plot_df_summ$methodrank]) +
      labs(title = 'MSE', x = 'MSE', y = 'Method', fill = 'Method') +
      scale_x_continuous(expand = expansion(mult = c(0, .05))) +
      scale_y_discrete(limits = rev) +
      theme(# axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5), 
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
            panel.grid.minor.y  = element_blank(), 
            panel.grid.major.y  = element_blank(),
            panel.grid.major.x  = element_blank(), 
            strip.background = element_rect(fill = NA), 
            axis.ticks.y = element_blank())

    ggsave(filename = sprintf('%s/mse_hori.pdf', save_folder), plot = p_mse, width = width, height = height) 
    if(save_ggplot) {
      saveRDS(p_mse, sprintf('%s/ggplot_mse_hori.rds', save_folder)) # save the ggplot objects (to adjust later as needed)
    }
    

    # === MSE nonzero theta ===
    if(plot_df_is_summ) {
      plot_df_summ = plot_df |> 
        filter(method != 'matcomp_linearreg') |> # exclude this... this performs badly
        filter(true_theta != 0) |>
        group_by(sim_distn, split_type, method, rank) |> 
        summarize(mse = mean(mse), .groups = 'drop') |> 
        mutate(sim_distn  = factor(sim_distn,  levels = c('pois', 'nb'),                   labels = c('Poisson', 'Negative Binomial')),
               split_type = factor(split_type, levels = c('nosamplesplit', 'samplesplit'), labels = c('Full Dataset', 'Sample Split')))
    } else {
      plot_df_summ = plot_df |> 
        filter(method != 'matcomp_linearreg') |> # exclude this... this performs badly
        filter(true_theta != 0) |>
        group_by(sim_distn, split_type, method, rank) |> 
        summarize(mse = mean((shrunk_value - true_theta)^2)) |> 
        mutate(sim_distn  = factor(sim_distn,  levels = c('pois', 'nb'),                   labels = c('Poisson', 'Negative Binomial')),
               split_type = factor(split_type, levels = c('nosamplesplit', 'samplesplit'), labels = c('Full Dataset', 'Sample Split')))
    }
    

    plot_df_summ$methodrank = mapply(FUN = methodrank_nicenames, method_name =  plot_df_summ$method, rank_ = plot_df_summ$rank) |> unname()
    plot_df_summ$methodrank = factor(plot_df_summ$methodrank, levels = methodrank_nicenames_order)
    
    
    # take the mse from unshrunk on Full Dataset
    # original_all_mse = plot_df_summ |> filter(method == 'unshrunkallcells')
    original_all_mse = plot_df_summ |> filter(method == 'unshrunk' & split_type == 'Full Dataset')
    original_all_mse = rbind(original_all_mse, original_all_mse |> mutate(split_type = 'Sample Split')) # repeat but change split type for plotting
    
    p_mse_nonzero = ggplot() +
      geom_col(data = plot_df_summ, aes(x = methodrank, y = mse, fill = methodrank), color = 'black', alpha = 1) + # MSEs
      geom_hline(data = original_all_mse, aes(yintercept = mse), color = '#DB1A1A', linewidth = .7, alpha = .6) +  # MSEs for 'naive'/original method
      facet_grid(rows = vars(sim_distn), cols = vars(split_type), scales = 'free_y', 
                 # labeller = as_labeller(c(sim_distn_nicename, split_type_nicename))
                 labeller = label_value
                 ) +
      scale_fill_discrete(palette = methodrank_colors[names(methodrank_colors) %in% plot_df_summ$methodrank]) +
      labs(title = 'MSE for nonzero Theta', x = 'method + rank', y = 'mse', fill = 'Method') +
      scale_y_continuous(expand = expansion(mult = c(0, .05))) +
      theme(# axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5), 
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
            panel.grid.major.x = element_blank(), strip.background = element_rect(fill = NA), 
            axis.ticks.x = element_blank())

    # p_mse_nonzero = ggplot() +
    #   geom_col(data = plot_df_summ, aes(x = methodrank, y = mse, fill = methodrank), color = 'black', alpha = 1) + # MSEs
    #   geom_hline(data = original_all_mse, aes(yintercept = mse), color = '#DB1A1A', linewidth = .7, alpha = .6) +  # MSEs for 'naive'/original methods
    #   facet_grid(rows = vars(sim_distn), cols = vars(split_type), scales = 'free_y', 
    #              # labeller = as_labeller(c(sim_distn_nicename, split_type_nicename))
    #              labeller = label_value
    #              ) +
    #   scale_fill_discrete(palette = methodrank_colors[names(methodrank_colors) %in% plot_df_summ$methodrank]) +
    #   labs(title = 'MSE for nonzero Theta', x = 'method + rank', y = 'mse', fill = 'Method') +
    #   scale_y_continuous(expand = expansion(mult = c(0, .05))) +
    #   theme(# axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5), 
    #         axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
    #         panel.grid.major.x = element_blank(), strip.background = element_rect(fill = NA), 
    #         axis.ticks.x = element_blank())


    ggsave(filename = sprintf('%s/mse_nonzero_vert.pdf', save_folder), plot = p_mse_nonzero, width = width, height = height) 
    if(save_ggplot) {
      saveRDS(p_mse_nonzero, sprintf('%s/ggplot_mse_nonzero_vert.rds', save_folder)) # save the ggplot objects (to adjust later as needed)
    }
    

    p_mse_nonzero = ggplot() +
      geom_col(data = plot_df_summ, aes(y = methodrank, x = mse, fill = methodrank), color = 'black', alpha = 1) + # MSEs
      geom_vline(data = original_all_mse, aes(xintercept = mse), color = '#DB1A1A', linewidth = .7, alpha = .6) +  # MSEs for 'naive'/original method
      facet_grid(cols = vars(sim_distn), rows = vars(split_type), scales = 'free', 
                 # labeller = as_labeller(c(sim_distn_nicename, split_type_nicename))
                 labeller = label_value
                 ) +
      scale_fill_discrete(palette = methodrank_colors[names(methodrank_colors) %in% plot_df_summ$methodrank]) +
      labs(title = 'MSE for nonzero Theta', x = 'MSE', y = 'Method', fill = 'Method') +
      scale_x_continuous(expand = expansion(mult = c(0, .05))) +
      scale_y_discrete(limits = rev) +
      theme(# axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5), 
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
            panel.grid.minor.y  = element_blank(), 
            panel.grid.major.y  = element_blank(),
            panel.grid.major.x  = element_blank(), 
            strip.background = element_rect(fill = NA), 
            axis.ticks.y = element_blank())

    ggsave(filename = sprintf('%s/mse_nonzero_hori.pdf', save_folder), plot = p_mse_nonzero, width = width, height = height) 
    if(save_ggplot) {
      saveRDS(p_mse_nonzero, sprintf('%s/ggplot_mse_nonzero_hori.rds', save_folder)) # save the ggplot objects (to adjust later as needed)
    }
    
}


#' 
#' @param plot_df_is_summ (boolean) if the given plot_df is already a summary dataframe (e.g. loaded from a saved file)
#'   FALSE: really tall dataframe
#'   TRUE: probably loaded in from prev run of create_summary_bytest_df 
h_plot_miscoverage <- function(plot_df, ranks, ALPHA, save_folder, height=NULL, width=NULL, save_ggplot=FALSE, plot_df_is_summ=FALSE) {
  if(is.null(height)){height=5} # default: height=5, width=8
  if(is.null(width) ){ width=8}
  # CI Coverage (should be 1-alpha proportion) # Mis-coverage rate to compare with mse (lower is better)

  # coverage rate
  # ggplot(plot_df |> 
  #        filter(method != 'matcomp_linearreg') |>
  #        mutate(isTrueThetaCovered = as.integer( lower_ci <= true_theta & true_theta <= upper_ci)) |>
  #        group_by(sim_distn, split_type, method, rank) |> 
  #        summarize(average_coverage = mean(isTrueThetaCovered), .groups = 'drop') |> 
  #        arrange(sim_distn, split_type, method, rank) |> 
  #        mutate(methodrank = factor(paste0(method, rank), 
  #                               levels = c("unshrunkallcellsNA",
  #                                          "unshrunkNA", 
  #                                          "matcomp_linearregNA", 
  #                                          paste0("matcomp_softImpute", ranks), 
  #                                          paste0("matdecomp_svd", ranks),
  #                                          paste0("matdecomp_sparsesvd", ranks), 
  #                                          "zerosNA", 
  #                                          "averageNA"))), 
  #      aes(x = methodrank, y = average_coverage)) +
  #   geom_col(color = 'black', fill = 'gray') +
  #   geom_hline(aes(yintercept = 1 - sim_results$ALPHA), color = 'orange', alpha = .7) +
  #   facet_grid(rows = vars(sim_distn), cols = vars(split_type), scales = 'free_y') +
  #   labs(title = 'Average EBCI Coverage', x = 'method + rank', y = 'average coverage') +
  #   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5), panel.grid.major.x = element_blank())

  # === prep nice labels and colors ===
  # === colors
  methodrank_colors = create_color_pallete_nicenames(ranks = ranks)
  # === labels
  # make the ordering for a various number of methodrank levels (even if not used): 
  # first method according to method_nicenames, then rank NA, 1, 2, ...
  methodrank_nicenames_order = c()
  for(cur_method in names(method_nicenames)) { # requires declaration of this list/vector (this is defined later in this file, right before the function methodrank_nicenames is defined)
    for(cur_rank in c(NA, ranks)) {
      methodrank_nicenames_order = c(methodrank_nicenames_order, methodrank_nicenames(method_name = cur_method, rank_ = cur_rank) |> unname())
    }
  }
      

  # === Miscoverage overall ===
  if(plot_df_is_summ) {
    plot_df_summ = plot_df |> 
      filter(method != 'matcomp_linearreg') |> # exclude this... this performs badly
      group_by(sim_distn, split_type, method, rank) |> 
      summarize(miscoverage_rate = mean(miscoverage_rate), .groups = 'drop') |> 
      mutate(sim_distn  = factor(sim_distn,  levels = c('pois', 'nb'),                   labels = c('Poisson', 'Negative Binomial')),
             split_type = factor(split_type, levels = c('nosamplesplit', 'samplesplit'), labels = c('Full Dataset', 'Sample Split')))
  } else {
    plot_df_summ = plot_df |> 
      filter(method != 'matcomp_linearreg') |> # exclude this... this performs badly
      mutate(isTrueThetaCovered = as.integer( lower_ci <= true_theta & true_theta <= upper_ci)) |>
      group_by(sim_distn, split_type, method, rank) |> 
      summarize(miscoverage_rate = 1 - mean(isTrueThetaCovered), .groups = 'drop') |> 
      mutate(sim_distn  = factor(sim_distn,  levels = c('pois', 'nb'),                   labels = c('Poisson', 'Negative Binomial')),
             split_type = factor(split_type, levels = c('nosamplesplit', 'samplesplit'), labels = c('Full Dataset', 'Sample Split')))
  }
 

  plot_df_summ$methodrank = mapply(FUN = methodrank_nicenames, method_name =  plot_df_summ$method, rank_ = plot_df_summ$rank) |> unname()
  plot_df_summ$methodrank = factor(plot_df_summ$methodrank, levels = methodrank_nicenames_order)
      
  # vertical
  p_miscoverage_vert = ggplot(plot_df_summ, aes(x = methodrank, y = miscoverage_rate, fill = methodrank)) +
      geom_col(color = 'black') +
      geom_hline(aes(yintercept = ALPHA), color = '#DB1A1A', linewidth = .7, alpha = .6) + 
      scale_fill_discrete(palette = methodrank_colors[names(methodrank_colors) %in% plot_df_summ$methodrank]) +
      facet_grid(rows = vars(sim_distn), cols = vars(split_type), scales = 'free_y') +
      labs(title = 'EBCI Miscoverage Rate', x = 'Method', y = 'Miscoverage Rate') +
      scale_y_continuous(expand = expansion(mult = c(0, .05))) +
      theme(# axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5), 
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
            panel.grid.major.y = element_blank(), 
            panel.grid.major.x = element_blank(), 
            strip.background = element_rect(fill = NA), 
            axis.ticks.x = element_blank())

  ggsave(filename = sprintf('%s/ebci_miscoveragerate_vert.pdf', save_folder), plot = p_miscoverage_vert, width = width, height = height) 
  if(save_ggplot) {
    saveRDS(p_miscoverage_vert, sprintf('%s/ggplot_miscoverage_vert.rds', save_folder)) # save the ggplot objects (to adjust later as needed)
  }
  

  # horizantal
  p_miscoverage_hori = ggplot(plot_df_summ, aes(y = methodrank, x = miscoverage_rate, fill = methodrank)) +
      geom_col(color = 'black') +
      geom_vline(aes(xintercept = ALPHA), color = '#DB1A1A', linewidth = .7, alpha = .6) + 
      scale_fill_discrete(palette = methodrank_colors[names(methodrank_colors) %in% plot_df_summ$methodrank]) +
      facet_grid(rows = vars(sim_distn), cols = vars(split_type), scales = 'fixed') +
      labs(title = 'EBCI Miscoverage Rate', x = 'Method', y = 'Miscoverage Rate') +
      scale_y_continuous(expand = expansion(mult = c(0, .05))) +
      scale_y_discrete(limits = rev) +
      theme(# axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5), 
              axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
              panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(),
              panel.grid.minor.x = element_blank(), 
              strip.background = element_rect(fill = NA), 
              axis.ticks.x = element_blank())

  ggsave(filename = sprintf('%s/ebci_miscoveragerate_hori.pdf', save_folder), plot = p_miscoverage_hori, width = width, height = height)  
  if(save_ggplot) {
    saveRDS(p_miscoverage_hori, sprintf('%s/ggplot_miscoverage_hori.rds', save_folder)) # save the ggplot objects (to adjust later as needed)
  } 
  

  return(NULL)
}



fishers_methodf <- function(x) {-2 * sum(log(x, base = exp(1)))}
fishers_pvalf <- function(s, k) {1 - pchisq(s, df = 2*k)}

#' 
h_plot_ebci_pval <- function(TODO) {



}

h5_2_2_plot_matrix_individual <- function(sim_results, color_limits=NULL) {
    dir.create(sprintf('%s/mat/', save_folder))
    if(is.null(color_limits)) {
      # set limits based on slightly expanded true Theta values
      cur_min = min(sim_results$Theta)
      cur_max = max(sim_results$Theta)
      cur_middle = (cur_min + cur_max) / 2
      scaling = 1.3
      
      color_limits = c(min(-1, ceiling(scaling * (cur_min - cur_middle) + cur_middle)),
                       floor(scaling * (cur_max - cur_middle) + cur_middle))
    }
  
    my_heatmap_fill_scale = scale_fill_gradient2(midpoint = 0, limits = color_limits, low = "#3A3A98" , high = "#832424", 
                                                 breaks = c(color_limits[1], 0, color_limits[2]), 
                                                 labels = c(paste0('<', color_limits[1]), 0, paste0('>', color_limits[2]))) 
    
    #' Function to display matrix in a standard way across all these matrices
    my_display_matrix <- function(X) {
      tryCatch(
        expr = {
           display_matrix(X, color_limits[1], color_limits[2]) + my_heatmap_fill_scale + theme(axis.title = element_blank())
        },
        error = function(e) {
           NULL
        })
     
    }
    
    
    # save Theta
    pl = my_display_matrix(sim_results$Theta) + labs(title = TeX(r'(True Effects $\Theta$)')) 
    ggsave(plot = pl, filename = sprintf('%s/mat/Theta.pdf', save_folder), width = 8, height = 8)
    ggsave(plot = pl, filename = sprintf('%s/mat/Theta.png', save_folder),
           width = 8, height = 8, units = 'in', dpi = 300)
    
    
    
    # save initial estimates
    for(est_method in names(sim_results$est_effects_matrices)) {
      for(splittype in names(sim_results$est_effects_matrices[[est_method]])) {
        pl = my_display_matrix(sim_results$est_effects_matrices[[est_method]][[splittype]]) + labs(title = paste0('Estimated Effects (', splittype, ')' ))
        ggsave(plot = pl, filename = sprintf('%s/mat/esteff_%s_%s.pdf', save_folder, substr(est_method, 1, 2), splittype), width = 8, height = 8)
        ggsave(plot = pl, filename = sprintf('%s/mat/esteff_%s_%s.png', save_folder, substr(est_method, 1, 2), splittype),
           width = 8, height = 8, units = 'in', dpi = 300)
      }
    }
    
    
  # save shrinkage targets/matrix approximations
  for(est_method in names(sim_results$estimate_matapprox)) { # this should really be sim_distn... not est_method
    for(splittype in names(sim_results$estimate_matapprox[[est_method]])) { # train test all
      for(approx_method in names(sim_results$estimate_matapprox[[est_method]][[splittype]])) {
        approx_method_shortname = substr(approx_method, start = nchar(approx_method) - 6, stop = nchar(approx_method))
        if(is.matrix(sim_results$estimate_matapprox[[est_method]][[splittype]][[approx_method]])) {
          pl = my_display_matrix(sim_results$estimate_matapprox[[est_method]][[splittype]][[approx_method]]) + 
               labs(title = paste0('Shrinkage Point/Matrix Approximation (', splittype, ' ', approx_method, ' ', ')' ))
          ggsave(plot = pl, filename = sprintf('%s/mat/approx_%s_%s_%s.pdf', save_folder, substr(est_method, 1, 2), splittype, approx_method_shortname), width = 8, height = 8)
          ggsave(plot = pl, filename = sprintf('%s/mat/approx_%s_%s_%s.png', save_folder, substr(est_method, 1, 2), splittype, approx_method_shortname),
                 width = 8, height = 8, units = 'in', dpi = 300)
        } else {
          for(r in sim_results$ranks) {
            pl = my_display_matrix(sim_results$estimate_matapprox[[est_method]][[splittype]][[approx_method]][[r]]) + 
                 labs(title = paste0('Shrinkage Point/Matrix Approximation (', splittype, ' ', approx_method, ' ', r, ')' ))
            ggsave(plot = pl, filename = sprintf('%s/mat/approx_%s_%s_%s_%s.pdf', save_folder, substr(est_method, 1, 2), splittype, approx_method_shortname, r), width = 8, height = 8)
            ggsave(plot = pl, filename = sprintf('%s/mat/approx_%s_%s_%s_%s.png', save_folder, substr(est_method, 1, 2), splittype, approx_method_shortname, r),
                   width = 8, height = 8, units = 'in', dpi = 300)
          }
        }
      }
    }
  }
  
  # save shrunk estimates
  for(est_method in names(sim_results$shrinkage_results)) {
    for(splittype in names(sim_results$shrinkage_results[[est_method]])) { # samplesplit or not
      if(splittype == 'samplesplit') {
        splittype_short =  'spl'
      } else if(splittype == 'nosamplesplit') {
        splittype_short =  'nospl'
      }
      
      for(approx_method in names(sim_results$shrinkage_results[[est_method]][[splittype]])) {
        approx_method_shortname = substr(approx_method, start = nchar(approx_method) - 6, stop = nchar(approx_method))
        if(is.data.frame(sim_results$shrinkage_results[[est_method]][[splittype]][[approx_method]])) {
          res_df = sim_results$shrinkage_results[[est_method]][[splittype]][[approx_method]]
          res_mat = reshape2::acast(res_df |> dplyr::select(x_idx, y_idx, shrunk_value), formula = y_idx ~ x_idx) # convert from tall df to wide mat
          
          pl = my_display_matrix(res_mat) + 
               labs(title = paste0('Shrunk Estimates (', splittype, ' ', approx_method, ' ', ')' ))
          ggsave(plot = pl, filename = sprintf('%s/mat/shrunk_%s_%s_%s.pdf', save_folder, substr(est_method, 1, 2), splittype_short, approx_method_shortname), width = 8, height = 8)
          ggsave(plot = pl, filename = sprintf('%s/mat/shrunk_%s_%s_%s.png', save_folder, substr(est_method, 1, 2), splittype_short, approx_method_shortname),
                 width = 8, height = 8, units = 'in', dpi = 300)
        } else {
          for(r in sim_results$ranks) {
            res_df = sim_results$shrinkage_results[[est_method]][[splittype]][[approx_method]][[r]]
            res_mat = reshape2::acast(res_df |> dplyr::select(x_idx, y_idx, shrunk_value), formula = y_idx ~ x_idx) # convert from tall df to wide mat
            
            pl = my_display_matrix(res_mat) + 
                 labs(title = paste0('Shrunk Estimates (', splittype, ' ', approx_method, ' ', r, ')' ))
            ggsave(plot = pl, filename = sprintf('%s/mat/shrunk_%s_%s_%s_%s.pdf', save_folder, substr(est_method, 1, 2), splittype_short, approx_method_shortname, r), width = 8, height = 8)
            ggsave(plot = pl, filename = sprintf('%s/mat/shrunk_%s_%s_%s_%s.png', save_folder, substr(est_method, 1, 2), splittype_short, approx_method_shortname, r),
                   width = 8, height = 8, units = 'in', dpi = 300)
          }
        }
      }
    }
  }
  
}

########################################################
##  MAIN FUNCTIONS
########################################################




#' make and save ebci matrices plots 
#' plots of matrices
#'  - (heatmap) thetas
#'    - true thetas
#'    - training estimated thetas
#'    - testing estimated thetas
#'    - 
#'  - (heatmap) matrix approx of test split (shrinkage points)
#'    - <across approx methods>
#'  - (heatmap) ebci estimates (shrunk points)
#'    - <across approx methods>
#'  - (heatmap) ebci significance (ebci ci covers 0 or not)
#'    - <across approx methods>
#'  - (points) pseudo p-values/ rejection rates?
#' @param est_method (character) est_method e.g. 'nb' or 'pois'
#' @param chosen_rank_to_plot (integer) chosen rank for plot, must be a valid input for 
#'        e.g. shrinkage_results[['nb']][['matdecomp_svd']][[chosen_rank_to_plot]]
#' @param shrinkage_results (list) of shrinkage results e.g. shrinkage_results  
#' @param est_effects_matrices (list) of estimated matrix effects e.g. est_effects_matrices
#' @param estimate_matapprox (list) of matrix approximations
#' @param allcells_results (list) of dataframes of the results from using all cells (no split)
#' @param Theta
#' @param save_folder (character) save folder path  e.g. "../../plots/simEBCICell/"
#' @returns nothing, but saves plot at
#' @examples 
#' make_matrix_ebci_plots('nb',   1, shrinkage_results, est_effects_matrices, estimate_matapprox, allcells_results, save_folder = "../../plots/simEBCICell/")
#' make_matrix_ebci_plots('nb',   2, shrinkage_results, est_effects_matrices, estimate_matapprox, allcells_results, save_folder = "../../plots/simEBCICell/")
#' make_matrix_ebci_plots('nb',   3, shrinkage_results, est_effects_matrices, estimate_matapprox, allcells_results, save_folder = "../../plots/simEBCICell/")
#' make_matrix_ebci_plots('pois', 1, shrinkage_results, est_effects_matrices, estimate_matapprox, allcells_results, save_folder = "../../plots/simEBCICell/")
#' make_matrix_ebci_plots('pois', 2, shrinkage_results, est_effects_matrices, estimate_matapprox, allcells_results, save_folder = "../../plots/simEBCICell/")
#' make_matrix_ebci_plots('pois', 3, shrinkage_results, est_effects_matrices, estimate_matapprox, allcells_results, save_folder = "../../plots/simEBCICell/")
make_matrix_ebci_plots <- function(est_method, chosen_rank_to_plot,
                                   shrinkage_results, est_effects_matrices, estimate_matapprox, allcells_results,
                                   save_folder, Theta, color_limits = c(-6, 5)) {
    # print('A make_matrix_ebci_plots')
    # color_limits = c(-6, 5)
    my_heatmap_fill_scale = scale_fill_gradient2(midpoint = 0, limits = color_limits, low = "#3A3A98" , high = "#832424", 
                                                 breaks = c(color_limits[1], 0, color_limits[2]), 
                                                 labels = c(paste0('<', color_limits[1]), 0, paste0('>', color_limits[2]))) 
    
    #' Function to display matrix in a standard way across all these matrices
    my_display_matrix <- function(X) {
      display_matrix(X, color_limits[1], color_limits[2]) + my_heatmap_fill_scale
    }
    
    # print('B ')
    #  - (heatmap) thetas -----------------------------------------------------------------------------------------------------------------
    p_ThetaTrue     = my_display_matrix(Theta)                                         + labs(title = TeX(r'($\Theta$)')) 
    p_ThetaEstTrain = my_display_matrix(est_effects_matrices[[est_method]][['train']]) + labs(title = TeX(r'($\hat{\Theta}$ (Train))'))
    p_ThetaEstTest  = my_display_matrix(est_effects_matrices[[est_method]][['test' ]]) + labs(title = TeX(r'($\hat{\Theta}$ (Test))')) 
    p_ThetaEstAll   = my_display_matrix(est_effects_matrices[[est_method]][['all'  ]]) + labs(title = TeX(r'($\hat{\Theta}$ (All))')) 
    
    # gridExtra::grid.arrange(p_ThetaTrue, p_ThetaEstTrain, p_ThetaEstTest,p_ThetaEstAll,
    #                         layout_matrix = matrix(c(NA,  1,  2,  3, 4), byrow = TRUE, ncol = 5))
    
    # print('C ')
    #  - (heatmap) matrix approx of train split (shrinkage points) --------------------------------------------------------------------------
    # estimate_matapprox[[est_method]][['train']] # <- shrinkage points (matrix approx on train split)
    # names(.) > "matcomp_linearreg"   "matcomp_softImpute"  "matdecomp_svd"       "matdecomp_sparsesvd"
    p_matapprox_matcomp_linearreg  = my_display_matrix(estimate_matapprox[[est_method]][['train']][['matcomp_linearreg']]) + 
                                         labs(title = TeX(r'($\tilde{\Theta}$ (Lin Reg))'))
    p_matapprox_matcomp_softImpute = my_display_matrix(estimate_matapprox[[est_method]][['train']][['matcomp_softImpute']][[chosen_rank_to_plot]]) + 
                                         labs(title = TeX(r'($\tilde{\Theta}$ (Soft Impute))'))
    p_matapprox_matdecomp_svd      = my_display_matrix(estimate_matapprox[[est_method]][['train']][['matdecomp_svd']][[chosen_rank_to_plot]]) + 
                                         labs(title = TeX(r'($\tilde{\Theta}$ (SVD))')) 
    p_matapproxmatdecomp_sparsesvd = my_display_matrix(estimate_matapprox[[est_method]][['train']][['matdecomp_sparsesvd']][[chosen_rank_to_plot]]) + 
                                         labs(title = TeX(r'($\tilde{\Theta}$ (Sparse SVD))'))

    p_matapprox_matcomp_linearreg_all  = my_display_matrix(estimate_matapprox[[est_method]][['all']][['matcomp_linearreg']]) + 
                                         labs(title = TeX(r'($\tilde{\Theta}$ (Lin Reg))'))
    p_matapprox_matcomp_softImpute_all = my_display_matrix(estimate_matapprox[[est_method]][['all']][['matcomp_softImpute']][[chosen_rank_to_plot]]) + 
                                         labs(title = TeX(r'($\tilde{\Theta}$ (Soft Impute))'))
    p_matapprox_matdecomp_svd_all      = my_display_matrix(estimate_matapprox[[est_method]][['all']][['matdecomp_svd']][[chosen_rank_to_plot]]) + 
                                         labs(title = TeX(r'($\tilde{\Theta}$ (SVD))')) 
    p_matapproxmatdecomp_sparsesvd_all = my_display_matrix(estimate_matapprox[[est_method]][['all']][['matdecomp_sparsesvd']][[chosen_rank_to_plot]]) + 
                                         labs(title = TeX(r'($\tilde{\Theta}$ (Sparse SVD))'))
    
    # gridExtra::grid.arrange(p_ThetaTrue, p_ThetaEstTrain, p_ThetaEstTest, 
    #                         p_matapprox_matcomp_linearreg, p_matapprox_matcomp_softImpute, p_matapprox_matdecomp_svd, p_matapproxmatdecomp_sparsesvd,
    #                         layout_matrix = matrix(c(NA,  1,  2,  3, NA, 
    #                                                  NA,  4,  5,  6,  7), byrow = TRUE, ncol = 5))
    
    # print('D ')
    #  - (heatmap) ebci estimates (shrunk points) ---------------------------------------------------------------------------------------------
    p_ebci_estimates_linearreg  = my_display_matrix(reshape2::acast(shrinkage_results[[est_method]][['samplesplit']][['matcomp_linearreg']], # turn back into matrix form
                                                                y_idx ~ x_idx, value.var = 'shrunk_value')) + 
                                     labs(title = TeX(r'(EB $\hat{\Theta}$ (Lin Reg))'))
    p_ebci_estimates_softImpute = my_display_matrix(reshape2::acast(shrinkage_results[[est_method]][['samplesplit']][['matcomp_softImpute']][[chosen_rank_to_plot]], 
                                                                y_idx ~ x_idx, value.var = 'shrunk_value')) + 
                                     labs(title = TeX(r'(EB $\hat{\Theta}$ (Soft Impute))'))
    p_ebci_estimates_svd        = my_display_matrix(reshape2::acast(shrinkage_results[[est_method]][['samplesplit']][['matdecomp_svd']][[chosen_rank_to_plot]], 
                                                                y_idx ~ x_idx, value.var = 'shrunk_value')) + 
                                     labs(title = TeX(r'(EB $\hat{\Theta}$ (SVD))'))
    p_ebci_estimates_sparsesvd  = my_display_matrix(reshape2::acast(shrinkage_results[[est_method]][['samplesplit']][['matdecomp_sparsesvd']][[chosen_rank_to_plot]], 
                                                                y_idx ~ x_idx, value.var = 'shrunk_value')) + 
                                     labs(title = TeX(r'(EB $\hat{\Theta}$ (Sparse SVD))'))

    p_ebci_estimates_linearreg_nosplit  = my_display_matrix(reshape2::acast(shrinkage_results[[est_method]][['nosamplesplit']][['matcomp_linearreg']], # turn back into matrix form
                                                                y_idx ~ x_idx, value.var = 'shrunk_value')) + 
                                     labs(title = TeX(r'(EB $\hat{\Theta}$ (Lin Reg))'))
    p_ebci_estimates_softImpute_nosplit = my_display_matrix(reshape2::acast(shrinkage_results[[est_method]][['nosamplesplit']][['matcomp_softImpute']][[chosen_rank_to_plot]], 
                                                                y_idx ~ x_idx, value.var = 'shrunk_value')) + 
                                     labs(title = TeX(r'(EB $\hat{\Theta}$ (Soft Impute))'))
    p_ebci_estimates_svd_nosplit        = my_display_matrix(reshape2::acast(shrinkage_results[[est_method]][['nosamplesplit']][['matdecomp_svd']][[chosen_rank_to_plot]], 
                                                                y_idx ~ x_idx, value.var = 'shrunk_value')) + 
                                     labs(title = TeX(r'(EB $\hat{\Theta}$ (SVD))'))
    p_ebci_estimates_sparsesvd_nosplit  = my_display_matrix(reshape2::acast(shrinkage_results[[est_method]][['nosamplesplit']][['matdecomp_sparsesvd']][[chosen_rank_to_plot]], 
                                                                y_idx ~ x_idx, value.var = 'shrunk_value')) + 
                                     labs(title = TeX(r'(EB $\hat{\Theta}$ (Sparse SVD))'))
    
    # gridExtra::grid.arrange(p_ThetaTrue, p_ThetaEstTrain, p_ThetaEstTest,  # 1 2 3
    #                         p_matapprox_matcomp_linearreg, p_matapprox_matcomp_softImpute, p_matapprox_matdecomp_svd, p_matapproxmatdecomp_sparsesvd, # 4 5 6 7
    #                         p_ebci_estimates_linearreg, p_ebci_estimates_softImpute, p_ebci_estimates_svd, p_ebci_estimates_sparsesvd,                # 8 9 10 11
    #                         layout_matrix = matrix(c(NA,  1,  2,  3,  NA, 
    #                                                  NA,  4,  5,  6,   7, 
    #                                                  NA,  8,  9,  10, 11), byrow = TRUE, ncol = 5))
    
    # print('E ')
    #  - (heatmap) ebci significance (ebci ci covers 0 or not) -----------------------------------------------------------------------------------------------------
    #' Make sig heatmap plots 
    #' @param df (dataframe) ebci results
    #' @param title (string) title
    #' p_ebci_sig_linearreg  = display_matrix(reshape2::acast(shrinkage_results[[est_method]][['matcomp_linearreg']] |>
    #'                                                               mutate(isSig = as.numeric(! ((lower_ci <= 0) & (0 <= upper_ci)))), # turn back into matrix form
    #'                                                             y_idx ~ x_idx, value.var = 'isSig')) + 
    #'                         scale_fill_continuous(breaks = c(0, 1), palette = c("white", "orangered2")) + 
    #'                                  labs(title = TeX(r'($0 \in$ EBCI (Lin Reg))'))
    #' make_ebci_sig_plot(df = shrinkage_results[[est_method]][['matdecomp_sparsesvd']][[chosen_rank_to_plot]] |> mutate(lower_ci = 1, upper_ci = 2), title = 'test')
    #' make_ebci_sig_plot(df = shrinkage_results[[est_method]][['matdecomp_sparsesvd']][[chosen_rank_to_plot]] |> mutate(lower_ci = -1, upper_ci = 1), title = 'test')
    make_ebci_sig_plot <- function(df, title) {
      display_matrix(reshape2::acast(df |> mutate(isSig = as.numeric(! ((lower_ci <= 0) & (0 <= upper_ci)))), # turn back into matrix form
                                     y_idx ~ x_idx, value.var = 'isSig')) +
      scale_fill_continuous(breaks = c(0, 1), limits = c(0, 1), palette = c("white", "orangered2")) +
      labs(title = title)
    }

    # For sample split procedure
    p_ebci_sig_linearreg  = make_ebci_sig_plot(shrinkage_results[[est_method]][['samplesplit']][['matcomp_linearreg']],                          TeX(r'($0 \in$ EBCI (Lin Reg))'))
    p_ebci_sig_softImpute = make_ebci_sig_plot(shrinkage_results[[est_method]][['samplesplit']][['matcomp_softImpute']][[chosen_rank_to_plot]],  TeX(r'($0 \in$ EBCI (Soft Impute))'))
    p_ebci_sig_svd        = make_ebci_sig_plot(shrinkage_results[[est_method]][['samplesplit']][['matdecomp_svd']][[chosen_rank_to_plot]],       TeX(r'($0 \in$ EBCI (SVD))'))
    p_ebci_sig_sparsesvd  = make_ebci_sig_plot(shrinkage_results[[est_method]][['samplesplit']][['matdecomp_sparsesvd']][[chosen_rank_to_plot]], TeX(r'($0 \in$ EBCI (Sparse SVD))'))
    
    # For no sample split procedure
    p_ebci_sig_linearreg_nosplit  = make_ebci_sig_plot(shrinkage_results[[est_method]][['nosamplesplit']][['matcomp_linearreg']],                          TeX(r'($0 \in$ EBCI (Lin Reg))'))
    p_ebci_sig_softImpute_nosplit = make_ebci_sig_plot(shrinkage_results[[est_method]][['nosamplesplit']][['matcomp_softImpute']][[chosen_rank_to_plot]],  TeX(r'($0 \in$ EBCI (Soft Impute))'))
    p_ebci_sig_svd_nosplit        = make_ebci_sig_plot(shrinkage_results[[est_method]][['nosamplesplit']][['matdecomp_svd']][[chosen_rank_to_plot]],       TeX(r'($0 \in$ EBCI (SVD))'))
    p_ebci_sig_sparsesvd_nosplit  = make_ebci_sig_plot(shrinkage_results[[est_method]][['nosamplesplit']][['matdecomp_sparsesvd']][[chosen_rank_to_plot]], TeX(r'($0 \in$ EBCI (Sparse SVD))'))
    
    
    
    # print('F 5 ')
    p_sig_all             = make_ebci_sig_plot(allcells_results[[est_method]], TeX(r'($0 \in$ CI (all))'))
    
    
    
    # print('G ')
    # gridExtra::grid.arrange(p_ThetaTrue, p_ThetaEstTrain, p_ThetaEstTest,  # 1 2 3
    #                         p_matapprox_matcomp_linearreg, p_matapprox_matcomp_softImpute, p_matapprox_matdecomp_svd, p_matapproxmatdecomp_sparsesvd, # 4   5  6  7
    #                         p_ebci_estimates_linearreg, p_ebci_estimates_softImpute, p_ebci_estimates_svd, p_ebci_estimates_sparsesvd,                # 8   9 10 11
    #                         p_ebci_sig_linearreg, p_ebci_sig_softImpute, p_ebci_sig_svd, p_ebci_sig_sparsesvd,                                        # 12 13 14 15
    #                         layout_matrix = matrix(c(NA,  1,  2,  3,  NA, 
    #                                                  NA,  4,  5,  6,   7, 
    #                                                  NA,  8,  9,  10, 11,
    #                                                  NA,  12, 13, 14, 15), byrow = TRUE, ncol = 5))
    
    # put all together ---------------------------------------------------------------------------------------------------------------------------------------
    # grob <- gridExtra::arrangeGrob(p_ThetaTrue, p_ThetaEstTrain, p_ThetaEstTest,  # 1 2 3
    #                               p_matapprox_matcomp_linearreg, p_matapprox_matcomp_softImpute, p_matapprox_matdecomp_svd, p_matapproxmatdecomp_sparsesvd, # 4   5  6  7
    #                               p_ebci_estimates_linearreg, p_ebci_estimates_softImpute, p_ebci_estimates_svd, p_ebci_estimates_sparsesvd,                # 8   9 10 11
    #                               p_ebci_sig_linearreg, p_ebci_sig_softImpute, p_ebci_sig_svd, p_ebci_sig_sparsesvd,                                        # 12 13 14 15
    #                               p_ThetaEstAll, p_sig_all, # 16 17
    #                               layout_matrix = matrix(c(1,  2,  3,  NA, 16,
    #                                                        4,  5,  6,   7, NA, 
    #                                                        8,  9,  10, 11, NA, 
    #                                                        12, 13, 14, 15, 17), byrow = TRUE, ncol = 5))

    grob <- gridExtra::arrangeGrob(p_ThetaTrue                 , p_ThetaEstTrain               , p_ThetaEstTest,                                            # 1 2 3 
                                  p_matapprox_matcomp_linearreg, p_matapprox_matcomp_softImpute, p_matapprox_matdecomp_svd, p_matapproxmatdecomp_sparsesvd, #  4  5  6  7  
                                  p_ebci_estimates_linearreg   , p_ebci_estimates_softImpute   , p_ebci_estimates_svd     , p_ebci_estimates_sparsesvd,     #  8  9 10 11
                                  p_ebci_sig_linearreg         , p_ebci_sig_softImpute         , p_ebci_sig_svd           , p_ebci_sig_sparsesvd,           # 12 13 14 15
                                  

                                  p_ThetaEstAll, p_sig_all,                                                                                                                  # 16 17
                                  p_matapprox_matcomp_linearreg_all, p_matapprox_matcomp_softImpute_all, p_matapprox_matdecomp_svd_all, p_matapproxmatdecomp_sparsesvd_all,  # 18 19 20 21
                                  p_ebci_estimates_linearreg_nosplit, p_ebci_estimates_softImpute_nosplit, p_ebci_estimates_svd_nosplit, p_ebci_estimates_sparsesvd_nosplit, # 22 23 24 25
                                  p_ebci_sig_linearreg_nosplit, p_ebci_sig_softImpute_nosplit, p_ebci_sig_svd_nosplit, p_ebci_sig_sparsesvd_nosplit,                         # 26 27 28 29
                                  layout_matrix = matrix(c(1,  2,  3,  NA,   16, 17, NA, NA,
                                                           4,  5,  6,   7,   18, 19, 20, 21,
                                                           8,  9,  10, 11,   22, 23, 24, 25,
                                                           12, 13, 14, 15,   26, 27, 28, 29          ), byrow = TRUE, ncol = 8))
    # print('H ')
    # gridExtra::grid.arrange(grob)
    ggsave(sprintf('%s/ebcimatrices_%s_rank=%d.pdf', save_folder, est_method, chosen_rank_to_plot), grob, width = 28, # 18, 
                                                                                                          height = 12)
    ggsave(sprintf('%s/ebcimatrices_%s_rank=%d.png', save_folder, est_method, chosen_rank_to_plot), grob, 
                      width = 28, height = 12, units = 'in', dpi = 300)

}











#' Sim pois and nb cells w sample splitting and perform EBCI shrinkage on matrix approximations
#' 
#' Can input a given treatment matrix Theta OR set Theta=NULL and specify P,G,rank to create a Theta
#' @param P (integer) number of perturbations/grnas/treatments  (used if is.null(Theta))
#' @param G (integer) number of genes/outcomes  (used if is.null(Theta))
#' @param rank (integer) rank of treatment effects matrix (used if is.null(Theta))
#' @param Theta (matrix) of treatment effects (either give a Theta or set=NULL to create a Theta using other parameters 
#'                       P,G,rank)
#' @param cell_distns (vector) of distributions for the cells, some subset of c('pois', 'nb')
#' @param N (integer) sample size of total treated cells 
#' @param N_control (integer) sample size of non-treated/control cells
#' @param pi_P (vector) of probabilities of assignment to each perturbation/treatment
#' @param nb_size (numeric or vector) of size parameter in negative binomial model
#' @param ranks (vector) of integers specifyin matrix approx ranks
#' @param ALPHA (numeric) value in [0,1] specifying alpha level for EBCI coverage
#' @param matapprox_methods (vector) of characters e.g c('matcomp_linearreg', 'matcomp_softImpute', 'matdecomp_svd', 'matdecomp_sparsesvd', 'spectralbiclust', 'zeros', 'average')
#' @param save_folder (character)
#' @param make_plots (boolean)
#' @param write_plot_df (boolean) whether to save the plot_df csv file 
#'                      (set FALSE to save space. all the info is already in the saved object: sim_results.rds,
#'                       and remade with the function h5_0_create_plot_df(...))
#' @param parallel (boolean) whether to run repetitions in parallel (uses future.apply::future_lapply(...). 
#'                 if true, make sure to remember to set up the session before calling this function. 
#'                 e.g. library(future.apply); plan(multisession, workers = 4))
#' @output list of
#'  all_sim_results = list(est_effects_matrices,
#'                         est_se_matrices,
#'                         estimate_matapprox,
#'                         shrinkage_results,
#'                         allcells_results,
#'                         Theta)
#' @examples
#' # === SETTING A === (small, easy, fast setting)
#' set.seed(12345)
#' setting_name = 'A'
#' P=5
#' G=10
#' rank=2
#' Theta = create_blocky_matrix(r = rank, n = P, m = G) # display_matrix(Theta)
#' pi_P = runif(P, min = .3, max = .7); pi_P = pi_P / sum(pi_P) # propensity score for each treatment
#' save_folder = sprintf('%s%s/', overall_save_folder, setting_name); dir.create(save_folder)
#' t0 = Sys.time()
#' sim_A_results = sim_EBCI_celllevel(
#'                    ALPHA=ALPHA,            # alpha testing level (for CIs)
#'                    P=P,                    # number of grnas/perturbations
#'                    G=G,                    # number of genes
#'                    rank=rank,              # true rank of theta matrix
#'                    Theta=Theta,            # true effects ( input NULL if want to simulate)
#'                    N=500,                 # number of treated cells
#'                    N_control=500,          # number of control cells
#'                    pi_P=pi_P,              # propensity score for each treatment
#'                    nb_size=3,              # contributes to the overdispersion of NB (smaller = more overdispersed)
#'                    ranks=c(1, 2, 3),       # ranks for matrix approximations
#'                    save_folder=save_folder # folder to save results and plots
#' )
#' t1 = Sys.time(); print(t1 - t0) # 4 mins
sim_EBCI_celllevel <- function(P, G, rank, Theta, cell_distns, 
                               N, N_control, pi_P, nb_size, 
                               ranks, ALPHA, 
                               matapprox_methods=NULL,
                               save_folder=NULL, make_plots=FALSE, repetitions=1, write_plot_df=FALSE, parallel=FALSE, calc_ebci_pvals=FALSE) {
  
  #' @param repetition (numeric) or NULL: will change where the results are saved 
  #'          (e.g. if NULL, save in save_folder, 
  #'                if integer, save in '<save_folder>/<repetition>')
  inner_function <- function(repetition=NULL, return_simresult=FALSE) {
    # 0. Theta
    # create treatment effects if not given
    if(is.null(Theta)) { # simulate true theta
      Theta = create_blocky_matrix(r = rank, n = P, m = G)
    } else { # use the given theta, overwrite input P and G
      P = nrow(Theta)
      G = ncol(Theta)
    }
    Theta_rownames = sprintf("grna%04.f", 1:P) # create names 
    Theta_colnames = sprintf("gene%04.f", 1:G) # prob should be adaptive to #digits ceiling(nrow(cur_estimates_matrix) / 10)
    row.names(Theta) = Theta_rownames
    colnames( Theta) = Theta_colnames
    
    
    # 1. Simulate cell-level data ----------------------------------------------------------------------------------------  
    # observed counts for each cell across genes follow pois or nb distribution
    print('Part 1: ')
    cell_data = h1_sim_cell_data(N=N, N_control=N_control, pi_P=pi_P, Theta=Theta, constant_coef=constant_coef, cell_distns=cell_distns, nb_size=nb_size)
    counts = cell_data$counts; grna = cell_data$grna; rm(cell_data)
    
    # 2. Estimate effects -------------------------------------------------------------------------------------------
    print('Part 2: ')
    est_eff_res = h2_est_effects(counts=counts, grna=grna, Theta=Theta, Theta_rownames=Theta_rownames, Theta_colnames=Theta_colnames)
    
    # est_effects_matrices = est_eff_res$est_effects_matrices
    # est_se_matrices      = est_eff_res$est_se_matrices
    # allcells_results     = est_eff_res$allcells_results
    # rm(est_eff_res)
    
    # 3. Matrix Decomp/Approximation -------------------------------------------------------------------------------------
    print('Part 3: ')
    estimate_matapprox = h3_approximate_matrices(est_effects_matrices=est_eff_res$est_effects_matrices, ranks=ranks, Theta_rownames=Theta_rownames, Theta_colnames=Theta_colnames, matapprox_methods=matapprox_methods)
    
    # 4. Perform EBCI to matrix approx -------------------------------------------------------------------------------------
    print('Part 4: ')
    h4_res = h4_perform_ebci(est_effects_matrices = est_eff_res$est_effects_matrices, 
                                        est_se_matrices      = est_eff_res$est_se_matrices, 
                                        estimate_matapprox   = estimate_matapprox, 
                                        ALPHA=ALPHA, Theta=Theta, Theta_rownames=Theta_rownames, Theta_colnames=Theta_colnames, ranks=ranks, calc_ebci_pvals=FALSE)
    shrinkage_results = h4_res$shrinkage_results
    ebci_params       = h4_res$ebci_params
    rm(h4_res)
    # 5. Plot and save rds  -----------------------------------------------------------------------------------------
    # if specified to save plots (e.g. valid save_folder and make_plots==TRUE)
    print('Part 5: ')
    all_sim_results = list(est_effects_matrices= est_eff_res$est_effects_matrices,
                           est_se_matrices     = est_eff_res$est_se_matrices,
                           estimate_matapprox  = estimate_matapprox,
                           shrinkage_results   = shrinkage_results,
                           ebci_params         = ebci_params,
                           allcells_results    = est_eff_res$allcells_results,
                           Theta               = Theta, # save some of the parameters used for this sim
                           P=P, G=G, N=N, N_control=N_control, pi_P=pi_P, nb_size=nb_size, ranks=ranks, 
                           ALPHA=ALPHA, save_folder=save_folder, repetition=repetition) 
    
    
    if(!is.null(save_folder) && dir.exists(save_folder)) {      
      
      if(is.null(repetition)) { # save in save_folder
        save_folder_rep = save_folder
      } else { # save in inner save_folder/<repetition>/
        save_folder_rep = sprintf('%s/%s/', save_folder, repetition)
        dir.create(save_folder_rep, showWarnings=FALSE) # outer folder should already exist
      }
      
      # save simulated results
      saveRDS(object = all_sim_results, file = sprintf('%s/sim_results.rds', save_folder_rep))
      
      plot_df = h5_0_create_plot_df(shrinkage_results=shrinkage_results, allcells_results=est_eff_res$allcells_results, ranks=ranks, ALPHA=ALPHA)
      if(write_plot_df) {
        write.csv(x = plot_df, file = sprintf('%s/sim_results_df.csv', save_folder_rep), row.names = FALSE)
      }
      
      

      # additionally, if we want to make plots
      if(make_plots) {
        # 5.1 some plots particular for each method, uses plot_shrink_results from utils/matrix_shrinkage.r
        # print('Part 5.1: ')
        # h5_1_plots_summary(shrinkage_results=shrinkage_results, save_folder=save_folder, ranks=ranks)
        
        # 5.2 plots of matrices (e.g. Theta, estimates, approx, shrunk, etc...)
        print('Part 5.2: ')
        h5_2_plots_matrix(shrinkage_results=shrinkage_results, 
                          est_effects_matrices=est_eff_res$est_effects_matrices, 
                          estimate_matapprox=estimate_matapprox, 
                          save_folder=save_folder_rep, 
                          allcells_results=est_eff_res$allcells_results, 
                          Theta=Theta,
                          ranks=ranks)
        
        # 5.3 mse
        print('Part 5.3: ')
        
        h5_3_plots_mse(plot_df=plot_df, ranks=ranks, save_folder=save_folder_rep)
      }


    }
      
    
    # return -------------------------------------------------------------------------
    if(return_simresult) {return(all_sim_results)}
  }
  

  # run inner_function based on number of repetitions
  if(repetitions == 1) {
    inner_function(repetition=NULL,  return_simresult=TRUE)
  } else if((repetitions  > 0) & (repetitions %% 1 == 0)) { # repetitions should be positive integer
    if(parallel) {
      future_lapply(1:repetitions, inner_function, future.seed = TRUE)
    } else {
      for(repetition in 1:repetitions) {
        inner_function(repetition=repetition, return_simresult=FALSE)
      }
    }
    
  } else {
    print('bad input: repetitions should be positive integer')
  }
}




#' Make plots from previously saved results
#' saved results should be a list
#' @param sim_results_filenames (character or vector of characters) to path and filename(s) of sim_results
#' (list) that is the 'all_sim_results' from the function sim_EBCI_celllevel
#' all_sim_results = list(est_effects_matrices= est_eff_res$est_effects_matrices,
#'                         est_se_matrices     = est_eff_res$est_se_matrices,
#'                         estimate_matapprox  = estimate_matapprox,
#'                         shrinkage_results   = shrinkage_results,
#'                         allcells_results    = est_eff_res$allcells_results,
#'                         Theta               = Theta)
#' @param save_folder (character) 
#' @param which_plots (list) of named booleans, naming which plots to make
#'   matrix, mse, miscoverage
#' @param write_plot_df (boolean) whether or not to save the plot_df
#' #@param create_default_plots (boolean) whether or not to create default plots that could have
#' #       been made during sim_EBCI_celllevel call (e.g. if sim_EBCI_celllevel(... make_plots=FALSE), then
#' #       the default plots were not made. Set this function's create_default_plots=TRUE to make these.)
#' @param plot_specs (list) of plot specifications when saving such as the size
#'          $mse $height=5
#'               $width=8    
#'          $miscoverage 
#'               $height=5
#'               $width=8            
#' @param use_saved_summary_df (boolean) whether to load in saved summary_df  (unimplemented.... maybe shouldn't)
#' @output 
make_plots_from_save <- function(sim_results_filenames, save_folder, which_plots, write_plot_df=FALSE, write_summary_df=FALSE, return_plot_df=FALSE, use_saved_summary_df=FALSE, plot_specs=NULL) {

  if(is.null(save_folder) || !dir.exists(save_folder)) {return('bad save_folder input')}
  

  # create a dataframe for plotting
  if(length(sim_results_filenames) == 1) { # --- 1 repetition
    sim_results = readRDS(sim_results_filenames)
    plot_df = h5_0_create_plot_df(shrinkage_results=sim_results$shrinkage_results, allcells_results=sim_results$allcells_results, ranks=sim_results$ranks, ALPHA=sim_results$ALPHA)
  } else {                                 # --- many repetitions
    plot_df = NULL
    for(fn in sim_results_filenames) {
      sim_results = readRDS(fn) 
      plot_df_ = h5_0_create_plot_df(shrinkage_results=sim_results$shrinkage_results, allcells_results=sim_results$allcells_results, ranks=sim_results$ranks, ALPHA=sim_results$ALPHA)
      plot_df_$filename = fn
      plot_df_$repetition = sim_results$repetition
      plot_df = rbind(plot_df, plot_df_); rm(plot_df_)
    }
  }
  ranks = sim_results$ranks 

  # if specified, save the plot_df (this is very large)
  if(write_plot_df) {
    write.csv(x = plot_df, file = sprintf('%s/sim_results_df.csv', save_folder), row.names = FALSE) # save plotting df
  }
  # if specified, save the summary_df (this is reasonable, try to save this)
  if(write_summary_df) {
    
    summary_df = plot_df |> 
       # filter(method != 'matcomp_linearreg') |>
       mutate(isTrueThetaCovered = as.integer( lower_ci <= true_theta & true_theta <= upper_ci)) |>
       group_by(sim_distn, split_type, method, rank) |> 
       summarize(miscoverage_rate = 1 - mean(isTrueThetaCovered), 
                 mse = mean((shrunk_value - true_theta)^2), 
                 .groups = 'drop') |> 
       arrange(sim_distn, split_type, method, rank) |> 
       mutate(methodrank = factor(paste0(method, rank), 
                              levels = c("unshrunkallcellsNA",
                                         "unshrunkNA", 
                                         "matcomp_linearregNA", 
                                         paste0("matcomp_softImpute", ranks), 
                                         paste0("matdecomp_svd", ranks),
                                         paste0("matdecomp_sparsesvd", ranks), 
                                         paste0("spectralbiclust", ranks), 
                                         paste0("spectralbiclust_threshold", ranks), 
                                         "zerosNA", 
                                         "averageNA")))
    
    # TODO: add fishers pval if there are many repetitions
    # calculate fishers pvals

    write.csv(x = summary_df, file = sprintf('%s/sim_summary_df.csv', save_folder), row.names = FALSE) # save plotting df
  }
  
  # 5.2 matrix plots- if multiple reps, will only plot using the LAST repetition e.g. A/1/sim_results.rds
  if(('matrix' %in% names(which_plots)) && which_plots$matrix) {
    # if we want to make defulat plots
    # 5.1 some plots particular for each method, uses plot_shrink_results from utils/matrix_shrinkage.r
    # print('Part 5.1: ')
    # h5_1_plots_summary(shrinkage_results=sim_results$shrinkage_results, save_folder=save_folder, ranks=sim_results$ranks)
    
    
    print('matrix')
    h5_2_plots_matrix(shrinkage_results=sim_results$shrinkage_results, 
                      est_effects_matrices=sim_results$est_effects_matrices, 
                      estimate_matapprox=sim_results$estimate_matapprox, 
                      save_folder=save_folder, 
                      allcells_results=sim_results$allcells_results, 
                      Theta=sim_results$Theta,
                      ranks=sim_results$ranks)
    
    
    
  }
  # plot individual matrices
  if(('matrix_individual' %in% names(which_plots)) && which_plots$matrix_individual) {
    if(length(sim_results_filenames) == 1) { # --- 1 repetition
      sim_results = readRDS(sim_results_filenames)
    } else {                                 # --- many repetitions
      sim_results = readRDS(sim_results_filenames[1])      
    }
    h5_2_2_plot_matrix_individual(sim_results=sim_results, color_limits=NULL)
    rm(sim_results)
  }

  # 5.3 mse
  if(('mse' %in% names(which_plots)) && which_plots$mse) {
    print('mse')
    p_mse = h5_3_plots_mse(plot_df=plot_df, ranks=sim_results$ranks, save_folder=save_folder, height=plot_specs$mse$height, width=plot_specs$mse$width) 
  }
  
  # miscoverage
  if(('miscoverage' %in% names(which_plots)) && which_plots$miscoverage) {
    print('miscoverage')
    p_miscoverage = h_plot_miscoverage(plot_df=plot_df, ranks=sim_results$ranks, ALPHA = sim_results$ALPHA, save_folder=save_folder, height=plot_specs$mse$height, width=plot_specs$mse$width) 
  }
  
  
  if(return_plot_df) {return(plot_df)}
  
}

########################################################
##  OLD FUNCTIONS
########################################################


# #' Sim pois and nb cells w sample splitting and perform EBCI shrinkage on matrix approximations
# sim_EBCI_celllevel_old <- function(ALPHA, P, G, rank, Theta,  N, N_control, pi_P,  nb_size, ranks, save_folder=NULL) {
#   # constant_coef = 1 # overall mean (mean for untreated/control cells)
#   constant_coef = runif(n = G, min = .1, max = 2) # overall mean across genes (mean for untreated/control cells)
  
#   # treatment effects 
#   if(is.null(Theta)) { # simulate true theta
#     Theta = create_blocky_matrix(r = rank, n = P, m = G)
#   } else { # use the given theta, overwrite input P and G
#     P = nrow(Theta)
#     G = ncol(Theta)
#   }
#   Theta_rownames = sprintf("grna%04.f", 1:P) # create names 
#   Theta_colnames = sprintf("gene%04.f", 1:G) # prob should be adaptive to #digits ceiling(nrow(cur_estimates_matrix) / 10)
#   row.names(Theta) = Theta_rownames
#   colnames(Theta) = Theta_colnames
  
  
#   # assign perturbation/grna treatment
#   grna = list()
#   grna$train =  t(rmultinom(n=N, size=1, prob=pi_P)) # n x p (#cells x #perturbations)
#   grna$test =  t(rmultinom(n=N, size=1, prob=pi_P)) 
#   # grna = t(rmultinom(n=N, size=1, prob=pi_P)) # n x p (#cells x #perturbations)
#   # grna_factor = apply(X=grna, MARGIN=1, FUN = function(v) {which(v==1)}) |> 
#   #               factor() # factor version of grna (which is one-hot encoded)
  
  
  
#   # 1. Simulate cell-level data -------------------------------------------------------------------------------------
#   cell_mean_effects = list()
#   cell_mean_effects$train = matrix(rep(constant_coef, times = N), byrow = TRUE, nrow = N) + grna$train %*% Theta
#   cell_mean_effects$test  = matrix(rep(constant_coef, times = N), byrow = TRUE, nrow = N) + grna$test  %*% Theta
  
  
  
#   # observed counts for each cell across genes follow pois or nb distribution
#   counts           = list()
#   counts[['pois']] = list()
#   counts[['nb']]   = list()
  
#   counts$pois$train = matrix(NA, nrow = N, ncol = G) # treated cells # train
#   counts$pois$test  = matrix(NA, nrow = N, ncol = G)                 # test 
#   counts$nb$train   = matrix(NA, nrow = N, ncol = G)                 # train 
#   counts$nb$test    = matrix(NA, nrow = N, ncol = G)                 # test  
  
#   counts$pois_control$train = matrix(NA, nrow = N_control, ncol = G) # control cells # train
#   counts$pois_control$test  = matrix(NA, nrow = N_control, ncol = G)                 # test
#   counts$nb_control$train   = matrix(NA, nrow = N_control, ncol = G)                 # train 
#   counts$nb_control$test    = matrix(NA, nrow = N_control, ncol = G)                 # test 
  
#   for(j in 1:G) {
#     counts$pois$train[,j] = rpois(  n=N,                lambda = exp(cell_mean_effects$train[,j])) 
#     counts$pois$test[,j]  = rpois(  n=N,                lambda = exp(cell_mean_effects$test[,j])) 
#     counts$nb$train[,j]   = rnbinom(n=N, size = nb_size,    mu = exp(cell_mean_effects$train[,j]))
#     counts$nb$test[,j]    = rnbinom(n=N, size = nb_size,    mu = exp(cell_mean_effects$test[,j]))
    
#     counts$pois_control$train[,j] = rpois(  n=N_control,                lambda = exp(constant_coef[j])) 
#     counts$pois_control$test[,j] = rpois(  n=N_control,                lambda = exp(constant_coef[j])) 
#     counts$nb_control$train[,j]    = rnbinom(n=N_control, size = nb_size,    mu = exp(constant_coef[j]))
#     counts$nb_control$test[,j]    = rnbinom(n=N_control, size = nb_size,    mu = exp(constant_coef[j]))
    
#   }
  
  
  
#   # 2. Estimate effects -------------------------------------------------------------------------------------------
  
#   # am i stupid? we could either do: 
#   # (1) gene g ~ grna p for each (g, p) or
#   # (2) gene g ~ \vec{grna} for each g (treat perturbations as levels in a categorical variable)
#   # these should be the same (as long as there aren't other variables, eg confounders) because
#   # the only variables are the treatment
#   # (but when there are other variables, (2) is better bc you can then use a large sample size for the
#   #  estimation of those coefficients. But, I think SCEPTRE does it one by one (1), so i will do it
#   #  this way in the simulation)
#   # (i mean, if you think that the confounders effect the outcome differently for each population,
#   #  cells w pert 1 vs cells w pert 2, then you should do (1). But that seems maybe unlikely, and it
#   #  would probably be helpful to be able to pool information from the larger sample size of cells)
  
  
#   # Do 'naively' at first: for each (i,j) pair, perform glm fit 
#   t0 = Sys.time()
#   est_effects_matrices = list()
#   est_se_matrices = list()
#   for(est_method in c('pois', 'nb')) {
#     est_effects_matrices[[est_method]] = list()
#     est_se_matrices     [[est_method]] = list()
#     for(split in c('train', 'test', 'all')) {
#       est_effects_matrices[[est_method]][[split]] = matrix(NA, nrow = P, ncol = G, dimnames = list(Theta_rownames, Theta_colnames))  
#       est_se_matrices     [[est_method]][[split]] = matrix(NA, nrow = P, ncol = G, dimnames = list(Theta_rownames, Theta_colnames))  
      
      
#       # estimate
#       for(i in 1:P) {
#         for(j in 1:G) {
#           if(split == 'all') {
#             pert_i_idx_train = which(grna[['train']][,i] == 1)
#             cur_dfij_train = data.frame(Y = c(counts[[est_method]][['train']][pert_i_idx_train,j], counts[[paste0(est_method, '_control')]][['train']][,j]),
#                                   A = c(rep(1, times = length(pert_i_idx_train)), rep(0, times = N_control)))
#             pert_i_idx_test = which(grna[['test']][,i] == 1)
#             cur_dfij_test = data.frame(Y = c(counts[[est_method]][['test']][pert_i_idx_test,j], counts[[paste0(est_method, '_control')]][['test']][,j]),
#                                   A = c(rep(1, times = length(pert_i_idx_test)), rep(0, times = N_control)))
#             cur_dfij = rbind(cur_dfij_train, cur_dfij_test)
#           } else { # train or test
#             pert_i_idx = which(grna[[split]][,i] == 1)
#             cur_dfij = data.frame(Y = c(counts[[est_method]][[split]][pert_i_idx,j], counts[[paste0(est_method, '_control')]][[split]][,j]),
#                                   A = c(rep(1, times = length(pert_i_idx)), rep(0, times = N_control)))
#           }
          
#           if(est_method == 'pois') {
#             glm_fit =       glm(   formula = 'Y ~ 1 + A', data = cur_dfij, family = poisson())
#           } else if(est_method == 'nb') {
#             glm_fit   = MASS::glm.nb(formula = 'Y ~ 1 + A', data = cur_dfij)
#           }
          
#           est_effects_matrices[[est_method]][[split]][i,j] = as.numeric(glm_fit$coefficients['A'])
#           est_se_matrices     [[est_method]][[split]][i,j] = summary(glm_fit)$coefficients['A', 'Std. Error']
          
#         }
#       }
#     }
#   }
#   t1 = Sys.time(); print(t1 - t0)
  
#   # combine results for 'all' using all cells (standard method)
#   allcells_results = list()
#   for(est_method in c('pois', 'nb')) {
#     cur_df = 
#       merge(reshape2::melt(est_effects_matrices[[est_method]][['all'  ]], value.name = 'estimate', varnames = c('grna', 'gene')),
#             merge(reshape2::melt(est_effects_matrices[[est_method]][['all'  ]] - qnorm(1 - ALPHA/2) * est_se_matrices[[est_method]][['all'  ]], 
#                      value.name = 'lower_ci', varnames = c('grna', 'gene')), 
#                   reshape2::melt(est_effects_matrices[[est_method]][['all'  ]] + qnorm(1 - ALPHA/2) * est_se_matrices[[est_method]][['all'  ]], 
#                      value.name = 'upper_ci', varnames = c('grna', 'gene')), 
#                   by = c('grna', 'gene')), 
#             by = c('grna', 'gene'))
#     cur_df = merge(cur_df, 
#                    reshape2::melt(est_se_matrices[[est_method]][['all'  ]], 
#                                   value.name = 'se', varnames = c('grna', 'gene')),
#                     by = c('grna', 'gene'))
#     cur_df = merge(cur_df, 
#                    reshape2::melt(Theta, 
#                                   value.name = 'true_theta', varnames = c('grna', 'gene')),
#                     by = c('grna', 'gene'))
#     cur_df = merge(merge(cur_df, 
#                          data.frame(y_idx = 1:length(Theta_rownames), grna = Theta_rownames),
#                          by = 'grna'), 
#                    data.frame(x_idx = 1:length(Theta_colnames), gene = Theta_colnames),
#                    by = 'gene')
#     allcells_results[[est_method]] = cur_df
#     rm(cur_df)
#   }
  
  
#   # 3. Matrix Decomp/Approximation -------------------------------------------------------------------------------------
#   estimate_matapprox = list()
#   estimate_matapprox_plots = list()
#   for(est_method in c('pois', 'nb')) {
#     estimate_matapprox      [[est_method]] = list()
#     estimate_matapprox_plots[[est_method]] = list()
    
#     for(split in c('train', 'test')) {
#       estimate_matapprox      [[est_method]][[split]] = list()
#       estimate_matapprox_plots[[est_method]][[split]] = list()
#       cur_mat = est_effects_matrices[[est_method]][[split]] # matrix to make approximations of
#       n = nrow(cur_mat); m = ncol(cur_mat)
      
      
#       # === matrix completion methods (linear reg, softImpute)
#       # --- Linear Reg
#       estimate_matapprox[[est_method]][[split]][['matcomp_linearreg']] =  
#                        matrix(NA, nrow = n, ncol = m, dimnames = list(Theta_rownames, Theta_colnames))
#       for(i in 1:n) {
#         for(j in 1:m) {
#           estimate_matapprox[[est_method]][[split]][['matcomp_linearreg']][i,j] = linearreg(cur_mat, i, j)
#         }
#       }
#       estimate_matapprox_plots[[est_method]][[split]][['matcomp_linearreg']] = display_matrix(estimate_matapprox[[est_method]][[split]][['matcomp_linearreg']])
      
      
#       # --- softImpute
#       estimate_matapprox      [[est_method]][[split]][['matcomp_softImpute']] = list()
#       estimate_matapprox_plots[[est_method]][[split]][['matcomp_softImpute']] = list()
#       for(r in ranks) {
#         estimate_matapprox[[est_method]][[split]][['matcomp_softImpute']][[r]] = 
#                        matrix(NA, nrow = n, ncol = m, dimnames = list(Theta_rownames, Theta_colnames))
#         matcomp_fun = create_matrix_completion_softImpute(X = cur_mat, rank_max = r, type = 'als') 
#         for(i in 1:n) {
#           for(j in 1:m) {
#             estimate_matapprox[[est_method]][[split]][['matcomp_softImpute']][[r]][i,j] = matcomp_fun(i, j)
#           }
#         }
#         estimate_matapprox_plots[[est_method]][[split]][['matcomp_softImpute']][[r]] =
#                                 display_matrix(estimate_matapprox[[est_method]][[split]][['matcomp_softImpute']][[r]])
#       }
      
#       # === matrix decomposition methods (svd, sparse svd, etc...)
#       estimate_matapprox      [[est_method]][[split]][['matdecomp_svd']]       = list()
#       estimate_matapprox      [[est_method]][[split]][['matdecomp_sparsesvd']] = list()
#       estimate_matapprox_plots[[est_method]][[split]][['matdecomp_svd']]       = list()
#       estimate_matapprox_plots[[est_method]][[split]][['matdecomp_sparsesvd']] = list()
      
#       svdres = svd(cur_mat, nu = max(ranks), nv = max(ranks))
#       cv.out = PMA::SPC.cv(cur_mat, sumabsvs = seq(1.2, min(5, sqrt(n), sqrt(m)), len = 10))
#       pmd_res = PMA::SPC(cur_mat, sumabsv=cv.out$bestsumabsv, K = max(ranks))
#       for(r in ranks) {
#         # --- SVD
#         cur_matapprox = svdres$u[, 1:r, drop=FALSE] %*% diag(svdres$d[1:r], nrow=r, ncol=r) %*% t(svdres$v[, 1:r, drop=FALSE]) # take up to r eigenvectors
#         row.names(cur_matapprox) = Theta_rownames; colnames(cur_matapprox) = Theta_colnames
#         estimate_matapprox[[est_method]][[split]][['matdecomp_svd']][[r]] = cur_matapprox
#         estimate_matapprox_plots[[est_method]][[split]][['matdecomp_svd']][[r]] = display_matrix(cur_matapprox)
#         rm(cur_matapprox)
        
#         # --- Sparse SVD
#         cur_matapprox = pmd_res$u[, 1:r, drop=FALSE] %*% diag(pmd_res$d[1:r], nrow=r, ncol=r) %*% t(pmd_res$v[, 1:r, drop=FALSE]) # take up to r eigenvectors
#         row.names(cur_matapprox) = Theta_rownames; colnames(cur_matapprox) = Theta_colnames
#         estimate_matapprox[[est_method]][[split]][['matdecomp_sparsesvd']][[r]] = cur_matapprox
#         estimate_matapprox_plots[[est_method]][[split]][['matdecomp_sparsesvd']][[r]] = display_matrix(cur_matapprox)
#         rm(cur_matapprox)
#       }
#     }
#   }
  
  
#   # 4. Perform EBCI -------------------------------------------------------------------------------------
#   # EBCI to arbitrary point
#   # use shrink_matrix from utils/matrix_shrinkage.r
#   shrinkage_results = list()
#   for(est_method in names(estimate_matapprox)) {
#     shrinkage_results[[est_method]] = list()
    
#     cur_estimates_mat = est_effects_matrices[[est_method]][['train']] # matrix of estimates to shrink
#     cur_se_mat        = est_se_matrices[[est_method]][['train']]      # matrix of se's of est to shrink
    
#     for(approx_method in names(estimate_matapprox[[est_method]][['test']])) {
#       # # could test by using: is.matrix or is.list
#       # is.matrix(estimate_matapprox[[est_method]][['matdecomp_sparsesvd']])
#       # is.matrix(estimate_matapprox[[est_method]][['matcomp_linearreg']])
#       # is.list(estimate_matapprox[[est_method]][['matdecomp_sparsesvd']])
#       # is.list(estimate_matapprox[[est_method]][['matcomp_linearreg']])
      
#       # if there are not ranks (eg linearreg)
#       if(is.matrix(estimate_matapprox[[est_method]][['test']][[approx_method]])) {
#         cur_shrinkagepoint_mat = estimate_matapprox[[est_method]][['test']][[approx_method]]
        
#         # shrink matrix
#         cur_shrink_res = 
#           shrink_matrix(unshrunk_mat = cur_estimates_mat,
#                       shrinkpoint_mat = cur_shrinkagepoint_mat,
#                       se_mat = cur_se_mat,
#                       ALPHA = ALPHA,
#                       return_ebci_obj=FALSE,
#                       weight_mat=NULL)
#         # add true theta values
#         cur_shrink_res = merge(cur_shrink_res, 
#                                reshape2::melt(Theta, varnames = c('grna', 'gene'), value.name = 'true_theta'),
#                                by = c('grna', 'gene'))
#         # add idx for plotting later
#         cur_shrink_res = merge(merge(cur_shrink_res, 
#                                      data.frame(y_idx = 1:length(Theta_rownames), grna = Theta_rownames),
#                                      by = 'grna'), 
#                                data.frame(x_idx = 1:length(Theta_colnames), gene = Theta_colnames),
#                                by = 'gene')
#         # save result
#         shrinkage_results[[est_method]][[approx_method]] = cur_shrink_res
  
        
#       } else { # if there are ranks (eg softImpute, svd, sparse svd, ..) 
#         for(r in ranks) { # if there was an error estimating some subset of ranks, then this might cause issues
#           cur_shrinkagepoint_mat = estimate_matapprox[[est_method]][['test']][[approx_method]][[r]]
#           # shrink matrix
#           cur_shrink_res = 
#             shrink_matrix(unshrunk_mat = cur_estimates_mat,
#                         shrinkpoint_mat = cur_shrinkagepoint_mat,
#                         se_mat = cur_se_mat,
#                         ALPHA = ALPHA,
#                         return_ebci_obj=FALSE,
#                         weight_mat=NULL)
#           # add true theta values
#           cur_shrink_res = merge(cur_shrink_res, 
#                                  reshape2::melt(Theta, varnames = c('grna', 'gene'), value.name = 'true_theta'),
#                                  by = c('grna', 'gene'))
#           # add idx for plotting later
#           cur_shrink_res = merge(merge(cur_shrink_res, 
#                                        data.frame(y_idx = 1:length(Theta_rownames), grna = Theta_rownames),
#                                        by = 'grna'), 
#                                  data.frame(x_idx = 1:length(Theta_colnames), gene = Theta_colnames),
#                                  by = 'gene')
#           # save result
#           shrinkage_results[[est_method]][[approx_method]][[r]] = cur_shrink_res
          
#         }
        
#       }
#     }
#   }
  
  
#   # 5. some plots -----------------------------------------------------------------------------------------
#   # if specified to save plots (e.g. !is.null(save_folder))
#   all_sim_results = list(est_effects_matrices=est_effects_matrices,
#                          est_se_matrices=est_se_matrices,
#                          estimate_matapprox=estimate_matapprox,
#                          shrinkage_results=shrinkage_results,
#                          allcells_results=allcells_results,
#                          Theta=Theta)
  
  
#   if(!is.null(save_folder) && dir.exists(save_folder)) {
#     saveRDS(object = all_sim_results, file = sprintf('%s/sim_results.rds', save_folder))
    
    
#     # 5.1 some plots particular for each method
#     # use plot_shrink_results from utils/matrix_shrinkage.r
    
#     for(est_method in names(shrinkage_results)) { # pois or nb
#       for(approx_method in names(shrinkage_results[[est_method]])) {
#         # # could test by using: is.matrix or is.list... should use list so it is consistent across different loops
#         # is.matrix(shrinkage_results[[est_method]][['matdecomp_sparsesvd']])
#         # is.matrix(shrinkage_results[[est_method]][['matcomp_linearreg']])
#         # is.list(shrinkage_results[[est_method]][['matdecomp_sparsesvd']])
#         # is.list(shrinkage_results[[est_method]][['matcomp_linearreg']])
#         # is.data.frame(shrinkage_results[[est_method]][['matdecomp_sparsesvd']])
#         # is.data.frame(shrinkage_results[[est_method]][['matcomp_linearreg']])
        
#         # if there are not ranks (eg linearreg)
#         if(is.data.frame(shrinkage_results[[est_method]][[approx_method]])) {
#           cur_shrinkage_plots_folder = sprintf('%s%s/', save_folder, approx_method)
#           dir.create(sprintf('%spoints/', cur_shrinkage_plots_folder), showWarnings = FALSE, recursive = TRUE)
#           cur_shrink_res = shrinkage_results[[est_method]][[approx_method]]
          
#           plot_shrink_results(shrink_df      = cur_shrink_res, 
#                                 plot_folder    = cur_shrinkage_plots_folder, 
#                                 order_rowscols = FALSE, 
#                                 grna_index     = NULL, 
#                                 gene_index     = NULL, 
#                                 unshrunk_ALPHA = ALPHA)
          
          
          
#         } else { # if there are ranks (eg softImpute, svd, sparse svd, ..) 
#           for(r in ranks) { # if there was an error estimating some subset of ranks, then this might cause issues
#             cur_shrinkage_plots_folder = sprintf('%s/%s_%02.f/', save_folder, approx_method, r)
#             dir.create(sprintf('%spoints/', cur_shrinkage_plots_folder), showWarnings = FALSE, recursive = TRUE)
#             cur_shrink_res = shrinkage_results[[est_method]][[approx_method]][[r]]
            
#             plot_shrink_results(shrink_df      = cur_shrink_res, 
#                                 plot_folder    = cur_shrinkage_plots_folder, 
#                                 order_rowscols = FALSE, 
#                                 grna_index     = NULL, 
#                                 gene_index     = NULL, 
#                                 unshrunk_ALPHA = ALPHA)
#           }
#         }
#       }
#     }
    
    
#     # 5.2 matrix plots
#     for(est_method in c('pois', 'nb')) {
#       for(r in ranks) {
#         make_matrix_ebci_plots(est_method = est_method,   
#                                chosen_rank_to_plot = r, 
#                                shrinkage_results = shrinkage_results, 
#                                est_effects_matrices = est_effects_matrices, 
#                                estimate_matapprox = estimate_matapprox, 
#                                save_folder = save_folder, 
#                                allcells_results = allcells_results, 
#                                color_limits = c(floor(min(Theta-.75)), ceiling(max(Theta+.75))))
#       }
#     }
    
#     # 5.3 shrinkage changes
#     plot_df = NULL
#     for(est_method in names(shrinkage_results)) {
#       plot_df_ = shrinkage_results[[est_method]][['matcomp_linearreg']] |> # pick any (would be same)
#                         dplyr::mutate(method = 'unshrunk', rank = NA) |>
#                         dplyr::mutate(shrinkage_point = unshrunk_value,   # set all values to the original estimates
#                                       shrunk_value    = unshrunk_value, 
#                                       lower_ci        = unshrunk_value - qnorm(1 - ALPHA/2) * se,
#                                       upper_ci        = unshrunk_value + qnorm(1 - ALPHA/2) * se) 
#       plot_df_ = rbind(plot_df_, 
#                        allcells_results[[est_method]] |> # all cells
#                          dplyr::mutate(method = 'unshrunkallcells', rank = NA) |>
#                          dplyr::mutate(shrinkage_point = NA,   # set all values to the original estimates
#                                        weight = NA, 
#                                        shrunk_value    = estimate,
#                                        unshrunk_value  = estimate)  |> 
#                          dplyr::select(all_of(colnames(plot_df_))) )
#       for(approx_method in names(shrinkage_results[[est_method]])) {
#          # if there are not ranks (eg linearreg)
#           if(is.data.frame(shrinkage_results[[est_method]][[approx_method]])) {
#             plot_df_ = rbind(plot_df_, 
#                             shrinkage_results[[est_method]][[approx_method]] |> 
#                             dplyr::mutate(method = approx_method, rank = NA))
#           } else { # if there are ranks (eg softImpute, svd, sparse svd, ..) 
#             for(r in ranks) { # if there was an error estimating some subset of ranks, then this might cause issues
#                plot_df_ = rbind(plot_df_, 
#                             shrinkage_results[[est_method]][[approx_method]][[r]] |> 
#                             dplyr::mutate(method = approx_method, rank = r))
#             }
            
#           }
#       }
#       plot_df = rbind(plot_df, plot_df_ |> dplyr::mutate(sim_distn = est_method, .before = 1))
#       rm(plot_df_)
#     }
    
#     plot_df$method = factor(plot_df$method, 
#                             levels = c("unshrunkallcells" , "unshrunk", "matcomp_linearreg", "matcomp_softImpute", "matdecomp_svd", "matdecomp_sparsesvd"))
    
#     # mse
#     p_mse = ggplot(plot_df |> 
#              group_by(sim_distn, method, rank) |> 
#              summarize(mse = mean((shrunk_value - true_theta)^2)) |> 
#              arrange(sim_distn, method, rank) |> 
#              mutate(methodrank = factor(paste0(method, rank), 
#                                     levels = c("unshrunkallcellsNA",
#                                                "unshrunkNA", 
#                                                "matcomp_linearregNA", 
#                                                paste0("matcomp_softImpute", ranks), 
#                                                paste0("matdecomp_svd", ranks),
#                                                paste0("matdecomp_sparsesvd", ranks)))), 
#            aes(x = methodrank, y = mse)) +
#       geom_col(color = 'black', fill = 'gray') +
#       facet_grid(rows = vars(sim_distn), scales = 'free_y')
    
    
#     ggsave(filename = sprintf('%s/shrinkage_mse.pdf', save_folder), plot = p_mse, width = 6, height = 4)

#   }
  
  
  
  
  
  
#   # return -------------------------------------------------------------------------
#   return(all_sim_results)
# }


# =================================
# === EBCI PVAL FUNCTIONS ===
# =================================

#' description: Create the pseudo p-value: Find the smallest alpha that does not cover 0 OR the largest alpha that still covers 0
#'  The CI is shrunk estimate +- cva * w_eb * sigma
#'  because of the way we call the function, we have our original estimate thetahat and shrinkagepoint
#'  and we shrink (thetahat - shrinkagepoint) towards 0 always to then get thetaeb(this is because we have diff shrinkage points for all the points)
#'  and then we add back the shrinkage point shrunkpoint = thetaeb + shrinkagepoint <-- this is the value to to input into this function, because the CI is centered around this point
#' @param thetashrunk (numeric) shrunkpoint = thetaeb + shrinkagepoint where is the thetaeb is the (thetahat - shrinkagepoint) shrunk towards 0
#' @param sigma (numeric)
#' @param web (numeric) shrinkage factor from ebci object result
#' @param mu2 (numeric) estimated mu2 from ebci object result
#' @param kappa (numeric) estimated kappa from ebci object result
#' @param threshold (numeric) stop condition: perform until alpha's changes are < threshold
#'                            this is the max mis calc error for the pseudo-pval (if MAX_ITER is not reached)
#' @param maxiter (numeric) stop condition: maximum number of iterations to perform
#' 
#' @returns pseudopval (numeric) \in [0,1]
#' 
#' 
get_ebci_pvals <- function(thetashrunk, sigma, web, mu2, kappa,
                                   alpha_threshold, maxiter) {
  # # params for running 
  # MAX_ITER_FOR_P = 20    # limit the number of iterations
  # # distance_threshold  = .001 # distance between (my_theta / (my_web * my_sigma)) and cva_pseudop(m2, kappa)
  # alpha_threshold = .0001 # perform until changes in alpha are small (ie this is the rounding of the p-value, there will be at most log2(1/alpha_threshold) iterations)
  # # estimate from overall ebci fit
  # my_mu2 = ebci_obj$mu2[['estimate']]
  # my_kappa = ebci_obj$kappa[['estimate']]
  # # params for each sample i 
  # my_theta = .2 # shrunk estimate: should be ebci_obj$df$th_eb + shrinkage_point (e.g. not the 'raw' the_eb)
  # my_sigma = 1 # initial estimate standard error
  # my_web = .3   # shrinkage factor

  # helpful calcs to just perform once
  sigma2_mu2 = sigma^2 / mu2  # = m2 in cva function input
  theta_websigma = abs(thetashrunk) / (web * sigma) # |theta| / (w_eb * sigma)   normalized estimate, always > 0



  iter = 1
  cur_min = 0; cur_max = 1 # range to search for pseudo p-value
  cur_alpha = .5 # start in the middle
  track_alphas = c(cur_alpha)
  # distance = 10000
  alpha_change = 1
  # while(distance > distance_threshold & iter <= MAX_ITER_FOR_P) { # by theta_websigma - cur_cva_alpha
  while(alpha_change > alpha_threshold & iter <= maxiter) { # by alpha/p-value change
   cur_cva_alpha = ebci::cva(m2 = sigma2_mu2, kappa = kappa, check=FALSE, alpha = cur_alpha)$cv
   # print(sprintf("%d: %.8f [%.2f, %.2f]: %.2f vs %.2f", 
   #                iter, cur_alpha, cur_min, cur_max, theta_websigma, cur_cva_alpha))
   
   # check if 0 is in CI with cur_alpha
   if(theta_websigma - cur_cva_alpha < 0) { # 0 \in CI: increase alpha to make CI smaller
     cur_min = cur_alpha  
     new_alpha = (cur_alpha + cur_max) / 2
   } else { # 0 \not\in CI: decrease alpha to make CI larger
     cur_max = cur_alpha
     new_alpha = (cur_alpha + cur_min) / 2
   }
   
   
   alpha_change = abs(new_alpha - cur_alpha)
   cur_alpha = new_alpha
   
   track_alphas = c(track_alphas, cur_alpha)
   
   # distance = theta_websigma - cur_cva_alpha
   iter = iter + 1
  }

  return(cur_alpha)
}






#' Using saved values from sim EBCI Cell, calculate the inverted/pseudo ebci pvals
#' Will save a nested list (of the sam structure as ebci_params and shrinkage_results) at the
#' location : <save_folder>/sim_result_ebci_pvals.rds
#' 
#' @param ebci_params (list) nested list of ebci_params (parameters estimated from ebci fit)
#' @param shrinkage_results (list) nested list of shrinkage_results (shrinkage results dataframe)
#' @param ranks (vector) of ranks used for matrix approximations
#' @param save_folder (character) path to folder to save result at
#' @example
#' suppressPackageStartupMessages(library(future.apply))
#' plan(multisession, workers = 8)  # or some other plan
#' source('../../utils/simEBCICell_utils.R') 
#' 
#' # testing/checking get_ebci_pvals
#' setting_name = 'A'
#' save_folder = sprintf('../../plots/simEBCICell/%s', setting_name)
#' sim_results = readRDS(sprintf('%s/1/sim_results.rds', save_folder))
#' 
#' # to hopefully speed up parallelization, clear some objects
#' ebci_params = sim_results$ebci_params
#' shrinkage_results = sim_results$shrinkage_results
#' ranks = sim_results$ranks
#' rm(sim_results); gc()
#' 
#' 
#' create_ebci_pvals_from_save(ebci_params=ebci_params, shrinkage_results=shrinkage_results, ranks=ranks, 
#'                             save_folder=save_folder, 
#'                             maxiter=10, alpha_threshold=.001)
#' 
#' 
#' test = readRDS(sprintf('%s/sim_result_ebci_pvals.rds', save_folder))
create_ebci_pvals_from_save <- function(ebci_params, shrinkage_results, ranks, save_folder, maxiter=10, alpha_threshold=.001) {
  
  if(is.null(save_folder) || !dir.exists(save_folder)) {return('bad save_folder input')}
  
  
  # HELPER FUNCTION
  
  #' calc ebci pvals current cur_ebci_params and cur_shrinkage_results
  #' 
  #' (this way, the code doesn't repeat for cases w ranks vs no ranks)
  calc_ebci_pvals_inner <- function(cur_ebci_params, cur_shrinkage_results) {
    # sim_distn = 'nb'
    # split_type = 'samplesplit'
    # approx_method = 'zeros'
    # 
    # maxiter = 10
    # alpha_threshold = .001
    # 
    # cur_ebci_params = ebci_params[[sim_distn]][[split_type]][[approx_method]]
    # cur_ebci_params
    # 
    # 
    # cur_shrinkage_results = shrinkage_results[[sim_distn]][[split_type]][[approx_method]]
    # cur_shrinkage_results |> arrange(gene, grna) # dataframe for each #pert * #gene
    
    inner_func <- function(thetashrunk, sigma, web) { # This created function will take in (thetashrunk, sigma, web)
      # call using the previously specified stopping conditions alpha_threshold and maxiter
      # and call using the current mu2 and kappa estimates
      
      tryCatch(expr = {get_ebci_pvals(thetashrunk=thetashrunk, sigma=sigma, web=web, 
                                              mu2=cur_ebci_params$mu2['estimate'], kappa=cur_ebci_params$kappa['estimate'], 
                                              alpha_threshold=alpha_threshold, maxiter=maxiter)
                      }, 
               error = function(e) {NA})
    
    }
    
    # # loop sequential
    # ebci_pvals_ = rep(NA, times = nrow(cur_shrinkage_results))
    # for(i in 1:100) {
    #   ebci_pvals[i] = inner_func(thetashrunk=cur_shrinkage_results$shrunk_value[i], sigma=cur_shrinkage_results$se[i], web=cur_shrinkage_results$w_eb[i])
    # }
    # 
    # 
    # 
    # # mapply sequential
    # ebci_pvals_ = mapply(FUN = inner_func, 
    #                      thetashrunk=cur_shrinkage_results$shrunk_value[1:100], 
    #                      sigma=cur_shrinkage_results$se[1:100], 
    #                      web=cur_shrinkage_results$w_eb[1:100])

      
    
    # future_mapply parallel
    ebci_pvals_ = future.apply::future_mapply(FUN = inner_func, 
                                              thetashrunk=cur_shrinkage_results$shrunk_value, 
                                              sigma=cur_shrinkage_results$se, 
                                              web=cur_shrinkage_results$w_eb)

    return(ebci_pvals_)
  }
  
  
  # START CALC EBCI PVAL
  ebci_pvals = list()
  
  # list with the same structure as ebci_params and shrinkage_results (e.g. $nb$samplesplit$...)
  for(sim_distn in names(ebci_params)) {
    ebci_pvals[[sim_distn]] = list()
    for(split_type in names(ebci_params[[sim_distn]])) {
      ebci_pvals[[sim_distn]][[split_type]] = list()
      for(approx_method in names(ebci_params[[sim_distn]][[split_type]])) {
        # if there are ranks
        if(!is.data.frame(shrinkage_results[[sim_distn]][[split_type]][[approx_method]])) { 
          # print(approx_method)
          ebci_pvals[[sim_distn]][[split_type]][[approx_method]] = list()
          for(r in ranks) {
            cur_ebci_params       = ebci_params      [[sim_distn]][[split_type]][[approx_method]][[r]]
            cur_shrinkage_results = shrinkage_results[[sim_distn]][[split_type]][[approx_method]][[r]]
            
            ebci_pvals[[sim_distn]][[split_type]][[approx_method]][[r]] = 
              calc_ebci_pvals_inner(cur_ebci_params=cur_ebci_params, cur_shrinkage_results=cur_shrinkage_results)
          }
        } else { # if there are no ranks
            cur_ebci_params       = ebci_params      [[sim_distn]][[split_type]][[approx_method]]
            cur_shrinkage_results = shrinkage_results[[sim_distn]][[split_type]][[approx_method]]
            
            ebci_pvals[[sim_distn]][[split_type]][[approx_method]] =
              calc_ebci_pvals_inner(cur_ebci_params=cur_ebci_params, cur_shrinkage_results=cur_shrinkage_results)
        }
      }
    }
  }
  
  
  saveRDS(object = ebci_pvals, file = sprintf('%s/sim_result_ebci_pvals.rds', save_folder))
}





# =================================
# === PLOTTING HELPER FUNCTIONS ===
# =================================


method_nicenames = c("unshrunkallcells"="Original (Full Dataset)",
                     "unshrunk"="Original", 
                     "matcomp_linearreg"="Linear Regression", 
                     "matcomp_softImpute"="softImpute", 
                     "matdecomp_svd"="SVD",
                     "matdecomp_sparsesvd"="Sparse SVD", 
                     "spectralbiclust"="Spectral Biclustering", 
                     "spectralbiclust_threshold"="Spectral Biclustering w/ Thresholding", 
                     "zeros"="Zeros", 
                     "average"="Average")

rank_nicenames = paste0('(rank=', 1:100, ')')

methodrank_nicenames <- function(method_name, rank_) {
  # method_name = 'zeros'
  # rank_ = 2
  # rank_ = NA
  # methodrank_nicenames(method_name = plot_df_summ$method[1], rank_ = plot_df_summ$rank[1]) |> unname()
  # mapply(FUN = methodrank_nicenames, method_name =  plot_df_summ$method[1:4], rank_ = plot_df_summ$rank[1:4]) |> unname()
  
  if(is.na(rank_)) {
    return(method_nicenames[method_name] |> unname())
  } else {
    return(paste0(method_nicenames[method_name], ' ', rank_nicenames[rank_]) |> unname())
  }
  
}




#' Create a vector of hex code colors from white to main_color of length number_of_colors
#' @param main_color (character) the hex code of the main color 
#' @param number_of_colors (integer) number of colors to return
create_colors <- function(main_color, number_of_colors) {
  colorRampPalette(c("white", main_color))(number_of_colors + 1)[1:number_of_colors+1]
}


#' create a named color pallette to use for plotting using the given ranks and the given 
#' 
#' distinct_colors = paletteer::paletteer_d("colorBlindness::paletteMartin")  # library(paletteer)
#' distinct_colors[c(2, 4, 6, 7, 12, 13, 14)]
#' @param ranks (vector) of pos integers
create_color_pallete_nicenames <- function(ranks) {
  distinct_colors = paletteer::paletteer_d("colorBlindness::paletteMartin")  # library(paletteer)
  distinct_colors[c(2, 4, 6, 7, 12, 13, 14)]


  if(FALSE) { # given the methods, make the list, but we should just make the full list and keep the colors consistent
    # one by one, create color palette based on the present method and ranks combinations
    # Original Estimates
    methodrank_colors = c("Original (Full Dataset)" = colorRampPalette(c("white", distinct_colors[2]))(2 + 1)[2],
                          "Original"                = colorRampPalette(c("white", distinct_colors[2]))(2 + 1)[3])
    if("matcomp_softImpute" %in% plot_df_summ$method) {
      temp_colors = create_colors(main_color = distinct_colors[4], number_of_colors = length(ranks)) # softImpute
      names(temp_colors) = mapply(FUN = methodrank_nicenames, method_name =  "matcomp_softImpute", rank_ = ranks) |> unname()
      methodrank_colors = c(methodrank_colors, temp_colors) ; rm(temp_colors)
    }
    if("matdecomp_svd" %in% plot_df_summ$method) {
      temp_colors = create_colors(main_color = distinct_colors[6], number_of_colors = length(ranks)) # matdecomp_svd
      names(temp_colors) = mapply(FUN = methodrank_nicenames, method_name =  "matcomp_softImpute", rank_ = ranks) |> unname()
      methodrank_colors = c(methodrank_colors, temp_colors) ; rm(temp_colors)
    }
    if("matdecomp_sparsesvd" %in% plot_df_summ$method) {
      temp_colors = create_colors(main_color = distinct_colors[7], number_of_colors = length(ranks)) # matdecomp_sparsesvd
      names(temp_colors) = mapply(FUN = methodrank_nicenames, method_name =  "matcomp_softImpute", rank_ = ranks) |> unname()
      methodrank_colors = c(methodrank_colors, temp_colors) ; rm(temp_colors)
    }
    if("spectralbiclust" %in% plot_df_summ$method) {
      temp_colors = create_colors(main_color = distinct_colors[12], number_of_colors = length(ranks)) # spectralbiclust
      names(temp_colors) = mapply(FUN = methodrank_nicenames, method_name =  "matcomp_softImpute", rank_ = ranks) |> unname()
      methodrank_colors = c(methodrank_colors, temp_colors) ; rm(temp_colors)
    }
    if("spectralbiclust_threshold" %in% plot_df_summ$method) {
      temp_colors = create_colors(main_color = distinct_colors[13], number_of_colors = length(ranks)) # spectralbiclust_threshold
      names(temp_colors) = mapply(FUN = methodrank_nicenames, method_name =  "matcomp_softImpute", rank_ = ranks) |> unname()
      methodrank_colors = c(methodrank_colors, temp_colors) ; rm(temp_colors)
    }
    
    methodrank_colors = c(methodrank_colors, 
                          "Zeros"   = colorRampPalette(c("white", distinct_colors[14]))(2 + 1)[2], # zeros/avg
                          "Average" = colorRampPalette(c("white", distinct_colors[14]))(2 + 1)[3])
  }
  
  # one by one, create color pallete based on the present method and ranks combinations
  # Original Estimates
  methodrank_colors = c("Original (Full Dataset)" = colorRampPalette(c("white", distinct_colors[2]))(2 + 1)[2],
                        "Original"                = colorRampPalette(c("white", distinct_colors[2]))(2 + 1)[3])
  # softImpute
  temp_colors = create_colors(main_color = distinct_colors[4], number_of_colors = length(ranks)) # softImpute
  names(temp_colors) = mapply(FUN = methodrank_nicenames, method_name =  "matcomp_softImpute", rank_ = ranks) |> unname()
  methodrank_colors = c(methodrank_colors, temp_colors) ; rm(temp_colors)
  
  # matdecomp_svd
  temp_colors = create_colors(main_color = distinct_colors[6], number_of_colors = length(ranks)) # matdecomp_svd
  names(temp_colors) = mapply(FUN = methodrank_nicenames, method_name =  "matdecomp_svd", rank_ = ranks) |> unname()
  methodrank_colors = c(methodrank_colors, temp_colors) ; rm(temp_colors)
  
  # matdecomp_sparsesvd
  temp_colors = create_colors(main_color = distinct_colors[7], number_of_colors = length(ranks)) # matdecomp_sparsesvd
  names(temp_colors) = mapply(FUN = methodrank_nicenames, method_name =  "matdecomp_sparsesvd", rank_ = ranks) |> unname()
  methodrank_colors = c(methodrank_colors, temp_colors) ; rm(temp_colors)
  
  # spectralbiclust
  temp_colors = create_colors(main_color = distinct_colors[12], number_of_colors = length(ranks)) # spectralbiclust
  names(temp_colors) = mapply(FUN = methodrank_nicenames, method_name =  "spectralbiclust", rank_ = ranks) |> unname()
  methodrank_colors = c(methodrank_colors, temp_colors) ; rm(temp_colors)
  
  # spectralbiclust_threshold
  temp_colors = create_colors(main_color = distinct_colors[13], number_of_colors = length(ranks)) # spectralbiclust_threshold
  names(temp_colors) = mapply(FUN = methodrank_nicenames, method_name =  "spectralbiclust_threshold", rank_ = ranks) |> unname()
  methodrank_colors = c(methodrank_colors, temp_colors) ; rm(temp_colors)
  
  
  methodrank_colors = c(methodrank_colors, 
                      "Zeros" = colorRampPalette(c("white", distinct_colors[14]))(2 + 1)[2], # zeros/avg
                    "Average" = colorRampPalette(c("white", distinct_colors[14]))(2 + 1)[3])
  methodrank_colors

}




#' create a tall dataframe with the columns:
#'  setting specifications: sim_distn, split_type, approx_method, rank,
#'  shrinkage relevant values: gene, grna, unshrunk_value, shrinkage_point, se, ...,
#'  AND importantly for this function, add the ebci_pvals saved in a different file
#'  if the pvals were calculated in the initial run, there will be ebci_pvals (initial) and ebci_pvals2 (from saved)
#'  else just ebci_pvals (from saved)
#'  @param sim_results (list)
#'  @param sim_result_ebci_vals (list) 
#'  @examples
#'  sim_results           = readRDS(sprintf('%s/sim_results.rds',           "../../plots/simEBCICell/A/9/"))
#'  sim_result_ebci_pvals = readRDS(sprintf('%s/sim_result_ebci_pvals.rds', "../../plots/simEBCICell/A/9/"))
#'  ebci_shrinkage_df = 
create_ebcipvals_df <- function(sim_results, sim_result_ebci_pvals) {
    # need: sim_results, sim_result_ebci_pvals, ranks
    # dataframe with cols: sim_distn, splittype, approx_method, rank
    df = NULL
    for(sim_distn in names(sim_result_ebci_pvals)) {
        for(splittype in names(sim_result_ebci_pvals[[sim_distn]])) {
            for(approx_method in names(sim_result_ebci_pvals[[sim_distn]][[splittype]])) {
                # # how to tell if there are ranks or not 
                # is.list(sim_result_ebci_pvals[[sim_distn]][[splittype]][['zeros']])         # FALSE
                # is.list(sim_result_ebci_pvals[[sim_distn]][[splittype]][['matdecomp_svd']]) # TRUE
                if(!is.list(sim_result_ebci_pvals[[sim_distn]][[splittype]][[approx_method]])) { # if there are not ranks (eg average, zeros, linearreg, ...)
                    df_             = sim_results$shrinkage_results[[sim_distn]][[splittype]][[approx_method]]  |> 
                                      dplyr::mutate(sim_distn=sim_distn, split_type=splittype, method=approx_method, rank=NA, .before=1)
                    if('ebci_pvals' %in% colnames(df_)) { # if already saved during the initial sim, then add as diff colname
                    df_$ebci_pvals2 = sim_result_ebci_pvals[[sim_distn]][[splittype]][[approx_method]] # will have 2 cols (should be same vals though)
                    } else {
                    df_$ebci_pvals  = sim_result_ebci_pvals[[sim_distn]][[splittype]][[approx_method]] 
                    df_$ebci_vals2  = NA
                    }
                } else {                                                                               # if there are ranks (eg svd...)
                    df_ = NULL
                    for(r in sim_results$ranks) {
                        df_r                 = sim_results$shrinkage_results[[sim_distn]][[splittype]][[approx_method]][[r]] |> 
                                               dplyr::mutate(sim_distn=sim_distn, split_type=splittype, method=approx_method, rank=r, .before=1)
                        if('ebci_pvals' %in% colnames(df_r)) { # if already saved during the initial sim, then add as diff colname
                            df_r$ebci_pvals2 = sim_result_ebci_pvals[[sim_distn]][[splittype]][[approx_method]][[r]] # will have 2 cols (should be same vals though)
                        } else {
                            df_r$ebci_pvals  = sim_result_ebci_pvals[[sim_distn]][[splittype]][[approx_method]][[r]] 
                            df_r$ebci_vals2  = NA
                        }
                        
                        df_ = rbind(df_, df_r); rm(df_r)
                    }
                }
                
                # bind together
                df = rbind(df, df_); rm(df_)
            }
        }
    }
    
    return(df)
}






#' create a summary df over all the repetitions
#'
#' will summarize: mean mse, average coverage rate, fishers pval, pval from 1 repetition
#' result should be (#treatments x #genes) x each setting long (e.g. summarize over the repetitions)
#' (e.g. for overall mse and miscoverage means, you need to group and summarize again)
#' @param df (dataframe) dataframe with shrinkage results, needs to have the columns
#' rep
#' settings: sim_distn  split_type     approx_method rank     y_idx x_idx
#' per test vars: gene grna unshrunk_value shrinkage_point se weight shrunk_value lower_ci upper_ci w_eb true_theta ebci_pvals
#'
#' @param ALPHA (numeric) alpha lvel the CIs from ebci were constructed at
create_summary_bytest_df <- function(df, ALPHA, save_folder=NULL) {
  # === Add average mse and average miscoverage rate
  df_summary = df |>
                 mutate(isTrueThetaCovered = as.integer( lower_ci <= true_theta & true_theta <= upper_ci)) |>
                 group_by(sim_distn, split_type, method, rank, gene, grna, true_theta) |>
                 summarize(mse = mean((shrunk_value - true_theta)^2),
                           miscoverage_rate = 1 - mean(isTrueThetaCovered),
                           nreps = n(),
                           .groups = 'drop')

  # === Add the "Naive" method (unshrunk glm)
  df_unshrunk_glm = df |> filter(method == 'zeros' & is.na(rank)) |> # no shrinkage/glm, pick any (would be same)
        mutate(method = 'unshrunk', rank = NA) |>
        mutate(shrinkage_point = unshrunk_value,   # set all values to the original estimates
               shrunk_value    = unshrunk_value,
               lower_ci        = unshrunk_value - qnorm(1 - ALPHA/2) * se,
               upper_ci        = unshrunk_value + qnorm(1 - ALPHA/2) * se,
               isTrueThetaCovered = as.integer( lower_ci <= true_theta & true_theta <= upper_ci)) |>
        group_by(sim_distn, split_type, method, rank, gene, grna, true_theta) |>
        summarize(mse = mean((shrunk_value - true_theta)^2),
                  miscoverage_rate = 1 - mean(isTrueThetaCovered),
                  nreps = n(),
                  .groups = 'drop')



  df_summary = rbind(df_summary, df_unshrunk_glm)
  
  # # === add back true_theta col... but will take more memory space (No, just group by w true_theta in prev)
  # df_summary = merge(df_summary,
  #                    df |> select(sim_distn,   split_type,   method,   rank,   gene,   grna, true_theta) |> distinct(), # it better be the same across?
  #                    by =      c('sim_distn', 'split_type', 'method', 'rank', 'gene', 'grna'),
  #                    all.x = TRUE) |>
  #              relocate(true_theta, .after = 'grna')

  # === add 1 ebci_pval
  chosen_rep = min(df$rep, na.rm = TRUE) # chosen repetition to save the ebci pvals from 1 rep
  ebci_pval_1 = df |> filter(rep == chosen_rep) |> select(sim_distn, split_type, method, rank, gene, grna, ebci_pvals)
  df_summary = merge(df_summary, ebci_pval_1, by = c('sim_distn', 'split_type', 'method', 'rank', 'gene', 'grna'), all.x = TRUE)
  # === add fishers combined ebci_pval
  fishers_df = df |> #  mutate(isTheta0 = (true_theta == 0)) |>
                         group_by(sim_distn, split_type, method, rank, gene, grna) |>
                         summarize(meanp = mean(ebci_pvals),
                                   fisher_stat = fishers_methodf(ebci_pvals),
                                   count = n(), .groups = 'drop')
  fishers_df$fishers_pval = mapply(FUN = fishers_pvalf, s= fishers_df$fisher_stat, k = fishers_df$count)

  df_summary = merge(df_summary, fishers_df, by = c('sim_distn', 'split_type', 'method', 'rank', 'gene', 'grna'), all.x = TRUE)




  # === save in folder if specified
  if(!is.null(save_folder) && dir.exists(save_folder)) {
    saveRDS(df_summary, sprintf('%s/sim_result_ebci_summary.rds', save_folder))
  }

  return(df_summary)
}

