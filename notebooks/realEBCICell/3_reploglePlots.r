# Use saved replogleShrinkage results to make plots
# 
# 
# Requires previously saved shrinkage results
# Specifically: <plot_path>/sceptre_obj_[train|test|all].rds 
# 
#

# Newly saved plots: <plot_path>/    (moved to separate file)
#  - mat/
#  -     estimates/   Estimates (original estimates)
#  -     matapprox/   Approximating Matrices (matrix Approximations)
#  -     shrinkage/   Shrinkage Matrices (matrix of Shrunk estimates)
#  - shrinkage/  


# ======================================================================================================================================
# --------------------------------------------------------------------------------------------------------------------------------------
#                                                                               ========================================================
#       SET PARAMETERS                                                          ========================================================
#                                                                               ========================================================
#                                                                               ========================================================
# --------------------------------------------------------------------------------------------------------------------------------------
# ======================================================================================================================================
# script params
# sceptre_save_path = '../../saves/sceptre/replogle/' # location where sceptre results are located
replogle_save_path= '../../saves/replogle/'         # location to save replogle approximations/shrinkage/etc. results
plot_path         = '../../plots/replogle/EBCI/'    # location to save plots of EBCI analysis on replogle dataset
# dir.exists(sceptre_save_path); dir.exists(replogle_save_path); dir.exists(plot_path)


# source('../../utils/matrix_shrinkage.r')

source('../../utils/get_ebci_pvals.r')






saved_approxmatrices = readRDS(sprintf('%s/approxmatrices.rds', replogle_save_path))
saved_approxmatrices

saved_EBCI_shrinkage_results = readRDS(sprintf('%s/EBCI_shrinkage_results.rds', replogle_save_path))
saved_EBCI_shrinkage_results


# ======================================================================================================================================
# --------------------------------------------------------------------------------------------------------------------------------------
#                            +   │─┐                                            ========================================================
#       SHRINKAGE PLOTS:     |─┐ │ └┐                                           ========================================================
#                            | ──┘  ────                                        ========================================================
#                            +---------+                                        ========================================================
# --------------------------------------------------------------------------------------------------------------------------------------
# ======================================================================================================================================
print("SHRINKAGE PLOTS: plot some shrinkage results")


# dir.create(sprintf('%s/shrinkage/', plot_path), showWarnings=FALSE) # plots/replogle/EBCI/shrinkage/


for(split_type in c('nosamplesplit', 'samplesplit')) {  
  print(sprintf("[%s]   - shrinkage plots with %s type", format( Sys.time(), format = "%F %r"), split_type))
  for(matapprox_method in matapprox_methods) {

    #' temporary function that plots some summary results of the shrinkage
    temp_plot_shrinkage <- function(cur_plot_folder) {
       
        tryCatch(expr = {
           
           # sparse SVD, rank=3
          # plot_folder = '../plots/replogle/shrink/spSVD03/'
          # shrink_df = read.csv("../saves/replogle/shrinkage/replogle_shrink_sparseSVD03.csv")
          # shrink_df = read.csv(sprintf("%s/ebci_shrinkage_df.csv", cur_plot_folder))
          shrink_res = readRDS(sprintf('%s/ebci_shrinkage_res.rds', cur_plot_folder)) # load in the prev saved ebci shrinkage results
          shrink_df = shrink_res$ebci_res; rm(shrink_res)                             # use the dataframe part of the result
          dir.create(sprintf('%s/points/', cur_plot_folder)); dir.create(sprintf('%s/heatmaps/', cur_plot_folder))
          plot_shrink_results(shrink_df=shrink_df, plot_folder=cur_plot_folder, order_rowscols=T, grna_index=grna_index, gene_index=gene_index, unshrunk_ALPHA=ALPHA)
          }, 
          error = function(e){print(sprintf("Errored at: %s", cur_plot_folder))}
        )
        return(NULL)
    }
    if(matapprox_methods_hasranks[[matapprox_method]]) {
      for(r in ranks) {
        temp_plot_shrinkage(cur_plot_folder = sprintf('%s/shrinkage/%s/%s/rank=%02.f/', plot_path, split_type, matapprox_method, r))
      }
    } else {
        temp_plot_shrinkage(cur_plot_folder = sprintf('%s/shrinkage/%s/%s/', plot_path, split_type, matapprox_method))
    }
  }
}





# ======================================================================================================================================
# --------------------------------------------------------------------------------------------------------------------------------------
#                            +   │─┐                                            ========================================================
#       EBCI P-VALUES:       |─┐ │ └┐                                           ========================================================
#                            | ──┘  ────                                        ========================================================
#                            +---------+                                        ========================================================
# --------------------------------------------------------------------------------------------------------------------------------------
# ======================================================================================================================================
print("EBCI P-values: Invert EBCI to p-values")



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





saved_approxmatrices = readRDS(sprintf('%s/approxmatrices.rds', replogle_save_path))
saved_approxmatrices

saved_EBCI_shrinkage_results = readRDS(sprintf('%s/EBCI_shrinkage_results.rds', replogle_save_path))
saved_EBCI_shrinkage_results





# need: ebci_obj fit (estimated parameters)
#       theta shrunk (shrunk estimate) + web (shrinkage factor)

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
      
      
      
      
      
      