

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
h1_sim_cell_data <- function(N, N_control, pi_P, Theta, constant_coef, nb_size) {

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
  counts$pois      = list()
  counts$nb        = list()

  counts$pois$treatment = list()
  counts$pois$control   = list()
  counts$nb$treatment   = list()
  counts$nb$control     = list()
  
  counts$pois$treatment$train = matrix(NA, nrow = N, ncol = G) # treated cells # train
  counts$pois$treatment$test  = matrix(NA, nrow = N, ncol = G)                 # test 
  counts$nb$treatment$train   = matrix(NA, nrow = N, ncol = G)                 # train 
  counts$nb$treatment$test    = matrix(NA, nrow = N, ncol = G)                 # test  
  
  counts$pois$control$train = matrix(NA, nrow = N_control, ncol = G) # control cells # train
  counts$pois$control$test  = matrix(NA, nrow = N_control, ncol = G)                 # test
  counts$nb$control$train   = matrix(NA, nrow = N_control, ncol = G)                 # train 
  counts$nb$control$test    = matrix(NA, nrow = N_control, ncol = G)                 # test 
  
  for(j in 1:G) {
    counts$pois$treatment$train[,j] = rpois(  n=N,                lambda = exp(cell_mean_effects$train[,j])) 
    counts$pois$treatment$test[,j]  = rpois(  n=N,                lambda = exp(cell_mean_effects$test[,j])) 
    counts$nb$treatment$train[,j]   = rnbinom(n=N, size = nb_size,    mu = exp(cell_mean_effects$train[,j]))
    counts$nb$treatment$test[,j]    = rnbinom(n=N, size = nb_size,    mu = exp(cell_mean_effects$test[,j]))
    
    counts$pois$control$train[,j]   = rpois(  n=N_control,                lambda = exp(constant_coef[j])) 
    counts$pois$control$test[,j]    = rpois(  n=N_control,                lambda = exp(constant_coef[j])) 
    counts$nb$control$train[,j]     = rnbinom(n=N_control, size = nb_size,    mu = exp(constant_coef[j]))
    counts$nb$control$test[,j]      = rnbinom(n=N_control, size = nb_size,    mu = exp(constant_coef[j]))
    
  }

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
  for(est_method in c('pois', 'nb')) {
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
  for(est_method in c('pois', 'nb')) {
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
#' requires function: create_matrix_completion_softImpute from util...<somwhere? matrix_utils.r?>
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
#' @output estimate_matapprox
h3_approximate_matrices <- function(est_effects_matrices, ranks, Theta_rownames=Theta_rownames, Theta_colnames=Theta_colnames) {
  estimate_matapprox = list()
  for(distn_name in c('pois', 'nb')) {
    estimate_matapprox      [[distn_name]] = list()
    
    # for(split in c('train', 'test', 'all')) {
    for(split in c('train', 'all')) { # only approx train and all split
      estimate_matapprox[[distn_name]][[split]] = list()
      cur_mat = est_effects_matrices[[distn_name]][[split]] # matrix to make approximations of
      # n = nrow(cur_mat); m = ncol(cur_mat)
      
      matapprox_methods = c('matcomp_linearreg', 'matcomp_softImpute', 'matdecomp_svd', 'matdecomp_sparsesvd')
      for(matapprox_method in matapprox_methods) {
        estimate_matapprox[[distn_name]][[split]][[matapprox_method]] = h3_1_approximate_matrix(mat = cur_mat,  method = matapprox_method, ranks = ranks,
                                                                                                Theta_rownames=Theta_rownames, Theta_colnames=Theta_colnames)
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
h4_perform_ebci <- function(est_effects_matrices, est_se_matrices, estimate_matapprox, 
                            ALPHA, Theta, Theta_rownames, Theta_colnames, ranks) {
  
  # use shrink_matrix from utils/matrix_shrinkage.r
  shrinkage_results = list()
  for(est_method in names(est_effects_matrices)) {
    shrinkage_results[[est_method]] = list()
    
    # do 2 categories of shrinkage: sample split (train towards test), no sample split (all towards all)
    for(splittype in c('samplesplit', 'nosamplesplit')) {
      shrinkage_results[[est_method]][[splittype]] = list()

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
          cur_shrink_res = 
            shrink_matrix(unshrunk_mat = cur_estimates_mat,
                          shrinkpoint_mat = cur_shrinkagepoint_mat,
                          se_mat = cur_se_mat,
                          ALPHA = ALPHA,
                          return_ebci_obj=FALSE,
                          weight_mat=NULL)
          # add true theta values
          cur_shrink_res = merge(cur_shrink_res, 
                                 reshape2::melt(Theta, varnames = c('grna', 'gene'), value.name = 'true_theta'),
                                 by = c('grna', 'gene'))
          # add idx for plotting later
          cur_shrink_res = merge(merge(cur_shrink_res, 
                                       data.frame(y_idx = 1:length(Theta_rownames), grna = Theta_rownames),
                                       by = 'grna'), 
                                 data.frame(x_idx = 1:length(Theta_colnames), gene = Theta_colnames),
                                 by = 'gene')
          # save result
          shrinkage_results[[est_method]][[splittype]][[approx_method]] = cur_shrink_res
    
          
        } else { # if there are ranks (eg softImpute, svd, sparse svd, ..) 
          for(r in ranks) { # if there was an error estimating some subset of ranks, then this might cause issues
            cur_shrinkagepoint_mat = estimate_matapprox[[est_method]][['train']][[approx_method]][[r]]
            # shrink matrix
            cur_shrink_res = 
              shrink_matrix(unshrunk_mat = cur_estimates_mat,
                          shrinkpoint_mat = cur_shrinkagepoint_mat,
                          se_mat = cur_se_mat,
                          ALPHA = ALPHA,
                          return_ebci_obj=FALSE,
                          weight_mat=NULL)
            # add true theta values
            cur_shrink_res = merge(cur_shrink_res, 
                                   reshape2::melt(Theta, varnames = c('grna', 'gene'), value.name = 'true_theta'),
                                   by = c('grna', 'gene'))
            # add idx for plotting later
            cur_shrink_res = merge(merge(cur_shrink_res, 
                                         data.frame(y_idx = 1:length(Theta_rownames), grna = Theta_rownames),
                                         by = 'grna'), 
                                   data.frame(x_idx = 1:length(Theta_colnames), gene = Theta_colnames),
                                   by = 'gene')
            # save result
            shrinkage_results[[est_method]][[splittype]][[approx_method]][[r]] = cur_shrink_res
            
          }
          
        }
      }
    } 




  }

  return(shrinkage_results)
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
  for(est_method in c('pois', 'nb')) {
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
h5_0_create_plot_df <- function(shrinkage_results, allcells_results, ranks) {
  plot_df = NULL
    for(est_method in names(shrinkage_results)) {
      for(splittype in c('samplesplit', 'nosamplesplit')) {
        plot_df_ = shrinkage_results[[est_method]][[splittype]][['matcomp_linearreg']] |> # pick any (would be same)
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
                                       unshrunk_value  = estimate)  |> 
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
      plot_df = rbind(plot_df, plot_df_ |> dplyr::mutate(sim_distn = est_method, .before = 1))
      rm(plot_df_)
      }
    }
    
    plot_df$method = factor(plot_df$method, 
                            levels = c("unshrunkallcells" , "unshrunk", "matcomp_linearreg", "matcomp_softImpute", "matdecomp_svd", "matdecomp_sparsesvd"))
    

    return(plot_df)
}

#' Helper function for: sim_EBCI_celllevel
#' @param plot_df
#' @param ranks
#' @param save_folder
#' @output saves plot at sprintf('%s/shrinkage_mse.pdf', save_folder)
h5_3_plots_mse <- function(plot_df, ranks, save_folder) {
    # mse
    p_mse = ggplot(plot_df |> 
             group_by(sim_distn, method, rank) |> 
             summarize(mse = mean((shrunk_value - true_theta)^2)) |> 
             arrange(sim_distn, method, rank) |> 
             mutate(methodrank = factor(paste0(method, rank), 
                                    levels = c("unshrunkallcellsNA",
                                               "unshrunkNA", 
                                               "matcomp_linearregNA", 
                                               paste0("matcomp_softImpute", ranks), 
                                               paste0("matdecomp_svd", ranks),
                                               paste0("matdecomp_sparsesvd", ranks)))), 
           aes(x = methodrank, y = mse)) +
      geom_col(color = 'black', fill = 'gray') +
      facet_grid(rows = vars(sim_distn), scales = 'free_y') +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
    
    ggsave(filename = sprintf('%s/shrinkage_mse.pdf', save_folder), plot = p_mse, width = 6, height = 6) 
    return()
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
    
    gridExtra::grid.arrange(p_ThetaTrue, p_ThetaEstTrain, p_ThetaEstTest,p_ThetaEstAll,
                            layout_matrix = matrix(c(NA,  1,  2,  3, 4), byrow = TRUE, ncol = 5))
    
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
    gridExtra::grid.arrange(grob)
    ggsave(sprintf('%s/ebcimatrices_%s_rank=%d.pdf', save_folder, est_method, chosen_rank_to_plot), grob, width = 28, # 18, 
                                                                                                          height = 12)

}











#' Sim pois and nb cells w sample splitting and perform EBCI shrinkage on matrix approximations
#' 
#' Can input a given treatment matrix Theta OR set Theta=NULL and specify P,G,rank to create a Theta
#' @param P (integer) number of perturbations/grnas/treatments  (used if is.null(Theta))
#' @param G (integer) number of genes/outcomes  (used if is.null(Theta))
#' @param rank (integer) rank of treatment effects matrix (used if is.null(Theta))
#' @param Theta (matrix) of treatment effects (either give a Theta or set=NULL to create a Theta using other parameters 
#'                       P,G,rank)
#' @param N (integer) sample size of total treated cells 
#' @param N_control (integer) sample size of non-treated/control cells
#' @param pi_P (vector) of probabilities of assignment to each perturbation/treatment
#' @param nb_size (numeric or vector) of size parameter in negative binomial model
#' @param ranks (vector) of integers specifyin matrix approx ranks
#' @param ALPHA (numeric) value in [0,1] specifying alpha level for EBCI coverage
#' @param save_folder (character)
#' @param make_plots (boolean)
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
sim_EBCI_celllevel <- function(P, G, rank, Theta,  
                               N, N_control, pi_P, nb_size, 
                               ranks, ALPHA, 
                               save_folder=NULL, make_plots=FALSE) {
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
  cell_data = h1_sim_cell_data(N=N, N_control=N_control, pi_P=pi_P, Theta=Theta, constant_coef=constant_coef, nb_size=nb_size)
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
  estimate_matapprox = h3_approximate_matrices(est_effects_matrices=est_eff_res$est_effects_matrices, ranks=ranks, Theta_rownames=Theta_rownames, Theta_colnames=Theta_colnames)
  
  # 4. Perform EBCI to matrix approx -------------------------------------------------------------------------------------
  print('Part 4: ')
  shrinkage_results = h4_perform_ebci(est_effects_matrices = est_eff_res$est_effects_matrices, 
                                      est_se_matrices      = est_eff_res$est_se_matrices, 
                                      estimate_matapprox   = estimate_matapprox, 
                                      ALPHA=ALPHA, Theta=Theta, Theta_rownames=Theta_rownames, Theta_colnames=Theta_colnames, ranks=ranks)
  
  # 5. Plot and save rds  -----------------------------------------------------------------------------------------
  # if specified to save plots (e.g. valid save_folder and make_plots==TRUE)
  print('Part 5: ')
  all_sim_results = list(est_effects_matrices= est_eff_res$est_effects_matrices,
                         est_se_matrices     = est_eff_res$est_se_matrices,
                         estimate_matapprox  = estimate_matapprox,
                         shrinkage_results   = shrinkage_results,
                         allcells_results    = est_eff_res$allcells_results,
                         Theta               = Theta, # save some of the parameters used for this sim
                         P=P, G=G, N=N, N_control=N_control, pi_P=pi_P, nb_size=nb_size, ranks=ranks, save_folder=save_folder) 
  
  
  if(!is.null(save_folder) && dir.exists(save_folder)) {
    # save simulated results
    saveRDS(object = all_sim_results, file = sprintf('%s/sim_results.rds', save_folder))
    plot_df = h5_0_create_plot_df(shrinkage_results=shrinkage_results, allcells_results=est_eff_res$allcells_results, ranks=ranks)
    write.csv(x = plot_df, file = sprintf('%s/sim_results_df.csv', save_folder), row.names = FALSE)

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
                        save_folder=save_folder, 
                        allcells_results=est_eff_res$allcells_results, 
                        Theta=Theta,
                        ranks=ranks)
      
      # 5.3 mse
      print('Part 5.3: ')
      
      h5_3_plots_mse(plot_df=plot_df, ranks=ranks, save_folder=save_folder)
    }


  }
    
  
  # return -------------------------------------------------------------------------
  return(all_sim_results)
}




#' Make plots from previously saved results
#' saved results should be a list
#' @param sim_results (list) that is the 'all_sim_results' from the function sim_EBCI_celllevel
#' all_sim_results = list(est_effects_matrices= est_eff_res$est_effects_matrices,
#'                         est_se_matrices     = est_eff_res$est_se_matrices,
#'                         estimate_matapprox  = estimate_matapprox,
#'                         shrinkage_results   = shrinkage_results,
#'                         allcells_results    = est_eff_res$allcells_results,
#'                         Theta               = Theta)
#' @param save_folder (character) 
#' @param create_default_plots (boolean) whether or not to create default plots that could have
#'        been made during sim_EBCI_celllevel call (e.g. if sim_EBCI_celllevel(... make_plots=FALSE), then
#'        the default plots were not made. Set this function's create_default_plots=TRUE to make these.)
#' @output 
make_plots_from_save <- function(sim_results, save_folder, create_default_plots=FALSE) {

  if(is.null(save_folder) || !dir.exists(save_folder)) {return('bad save_folder input')}

  

  # create a dataframe for plotting
  plot_df = h5_0_create_plot_df(shrinkage_results=sim_results$shrinkage_results, allcells_results=sim_results$est_eff_res$allcells_results, ranks=sim_results$ranks)
  write.csv(x = plot_df, file = sprintf('%s/sim_results_df.csv', save_folder), row.names = FALSE)


  if(create_default_plots) {
    # if we want to make defulat plots
    # 5.1 some plots particular for each method, uses plot_shrink_results from utils/matrix_shrinkage.r
    # print('Part 5.1: ')
    # h5_1_plots_summary(shrinkage_results=sim_results$shrinkage_results, save_folder=save_folder, ranks=sim_results$ranks)
    
    # 5.2 matrix plots
    print('Part 5.2: ')
    h5_2_plots_matrix(shrinkage_results=sim_results$shrinkage_results, 
                      est_effects_matrices=sim_results$est_eff_res$est_effects_matrices, 
                      estimate_matapprox=sim_results$estimate_matapprox, 
                      save_folder=save_folder, 
                      allcells_results=sim_results$est_eff_res$allcells_results, 
                      Theta=sim_results$Theta,
                      ranks=sim_results$ranks)
    
    # 5.3 mse
    print('Part 5.3: ')
    h5_3_plots_mse(plot_df=plot_df, ranks=sim_results$ranks, save_folder=save_folder) 
  }

  


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
