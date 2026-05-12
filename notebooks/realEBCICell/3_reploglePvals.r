# Use saved replogleShrinkage results to make plots
# 
# Run in order:
# 1_...
# 2_replogleShrinkage.r
# 3_reploglePvals.r
# 
# Requires previously saved shrinkage results
# Specifically: <plot_path>/sceptre_obj_[train|test|all].rds 
# 
# 
# 



suppressPackageStartupMessages(library(future.apply))
plan(multisession, workers = 20)  # or some other plan
# plan(sequential)
   

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


CALC_EBCI_PVALS     = TRUE   # calculate pvals or not
RECALC_EBCI_PVALS   = FALSE # If there were previously saved pvals, recalculate all (TRUE). Or only calc missing and save time (FALSE)
ASSEMBLE_EBCI_PVALS = TRUE



ALPHA_THRESHOLD = .001
MAXITER = 10
PARALLEL = TRUE




# source('../../utils/matrix_shrinkage.r')

source('../../utils/get_ebci_pvals.r')

replogleShrinkageParams = readRDS(file = sprintf('%s/replogleShrinkageParams.rds', replogle_save_path))




# saved_approxmatrices = readRDS(sprintf('%s/approxmatrices.rds', replogle_save_path))
# saved_approxmatrices
# 
# saved_EBCI_shrinkage_results = readRDS(sprintf('%s/EBCI_shrinkage_results.rds', replogle_save_path))
# saved_EBCI_shrinkage_results


# Use saved replogleShrinkage results to make plots
# 
# Run in order:
# 1_...
# 2_replogleShrinkage.r
# 3_reploglePvals.r
# 
# Requires previously saved shrinkage results
# Specifically: <plot_path>/sceptre_obj_[train|test|all].rds 
# 
# 
# 
# Newly saved plots: <plot_path>/    (moved to separate file)
#  - mat/
#  -     estimates/   Estimates (original estimates)
#  -     matapprox/   Approximating Matrices (matrix Approximations)
#  -     shrinkage/   Shrinkage Matrices (matrix of Shrunk estimates)
#  - shrinkage/  
#' 
#' 
#' # ======================================================================================================================================
#' # --------------------------------------------------------------------------------------------------------------------------------------
#' #                            +   │─┐                                            ========================================================
#' #       SHRINKAGE PLOTS:     |─┐ │ └┐                                           ========================================================
#' #                            | ──┘  ────                                        ========================================================
#' #                            +---------+                                        ========================================================
#' # --------------------------------------------------------------------------------------------------------------------------------------
#' # ======================================================================================================================================
#' print("SHRINKAGE PLOTS: plot some shrinkage results")
#' 
#' 
#' # dir.create(sprintf('%s/shrinkage/', plot_path), showWarnings=FALSE) # plots/replogle/EBCI/shrinkage/
#' 
#' 
#' for(split_type in c('nosamplesplit', 'samplesplit')) {  
#'   print(sprintf("[%s]   - shrinkage plots with %s type", format( Sys.time(), format = "%F %r"), split_type))
#'   for(matapprox_method in matapprox_methods) {
#' 
#'     #' temporary function that plots some summary results of the shrinkage
#'     temp_plot_shrinkage <- function(cur_plot_folder) {
#'        
#'         tryCatch(expr = {
#'            
#'            # sparse SVD, rank=3
#'           # plot_folder = '../plots/replogle/shrink/spSVD03/'
#'           # shrink_df = read.csv("../saves/replogle/shrinkage/replogle_shrink_sparseSVD03.csv")
#'           # shrink_df = read.csv(sprintf("%s/ebci_shrinkage_df.csv", cur_plot_folder))
#'           shrink_res = readRDS(sprintf('%s/ebci_shrinkage_res.rds', cur_plot_folder)) # load in the prev saved ebci shrinkage results
#'           shrink_df = shrink_res$ebci_res; rm(shrink_res)                             # use the dataframe part of the result
#'           dir.create(sprintf('%s/points/', cur_plot_folder)); dir.create(sprintf('%s/heatmaps/', cur_plot_folder))
#'           plot_shrink_results(shrink_df=shrink_df, plot_folder=cur_plot_folder, order_rowscols=T, grna_index=grna_index, gene_index=gene_index, unshrunk_ALPHA=ALPHA)
#'           }, 
#'           error = function(e){print(sprintf("Errored at: %s", cur_plot_folder))}
#'         )
#'         return(NULL)
#'     }
#'     if(matapprox_methods_hasranks[[matapprox_method]]) {
#'       for(r in ranks) {
#'         temp_plot_shrinkage(cur_plot_folder = sprintf('%s/shrinkage/%s/%s/rank=%02.f/', plot_path, split_type, matapprox_method, r))
#'       }
#'     } else {
#'         temp_plot_shrinkage(cur_plot_folder = sprintf('%s/shrinkage/%s/%s/', plot_path, split_type, matapprox_method))
#'     }
#'   }
#' }




# //////////////////////////////////////////////////////////////////////////////////////////////////
# ================================== EBCI P-VALUES: ================================================
# //////////////////////////////////////////////////////////////////////////////////////////////////
if(CALC_EBCI_PVALS) {
  print("CALCULATE EBCI P-VALUES: Invert EBCI to p-values =================================================")
  
  # need: ebci_obj fit (estimated parameters)
  #       theta shrunk (shrunk estimate) + web (shrinkage factor)
  
  
  # saved_approxmatrices = readRDS(sprintf('%s/approxmatrices.rds', replogle_save_path))
  # saved_approxmatrices
  
  # saved_EBCI_shrinkage_results = readRDS(sprintf('%s/EBCI_shrinkage_results.rds', replogle_save_path))
  # saved_EBCI_shrinkage_results
  # shrink_results = saved_EBCI_shrinkage_results
  
  
  shrink_results = readRDS(sprintf('%s/EBCI_shrinkage_results.rds', replogle_save_path))
  
  
  # try to read in previously saved ebcipvals to skip recalc if already present
  if((!RECALC_EBCI_PVALS) && ('EBCI_shrinkage_ebcipvals.rds' %in% list.files(replogle_save_path))) {
    print('     (loading in prev saved EBCI_shrinkage_ebcipvals.rds)')
    ebci_pvals     = readRDS(sprintf('%s/EBCI_shrinkage_ebcipvals.rds', replogle_save_path))
  } else {
    ebci_pvals = list() # save in same format as shrink_results
    for(split_type in names(shrink_results)) { # samplesplit, nosamplesplit
      ebci_pvals[[split_type]] = list()
    }
  }
  
  
  
  for(split_type in names(shrink_results)) { # samplesplit, nosamplesplit
    print(sprintf("[%s]   - EBCI pvals with %s split type", format( Sys.time(), format = "%F %r"), split_type))
    for(approx_method in names(shrink_results[[split_type]])) {
      print(sprintf("[%s]   - EBCI pvals with %s split type, %s method", format( Sys.time(), format = "%F %r"), split_type, approx_method))
      
      # if(!('ebci_res' %in% names(shrink_results[[split_type]][[approx_method]]))) { # if there are ranks
      #   for(r in length(shrink_results[[split_type]][[approx_method]])) {
      #     
      #   }
      # } else {                                                                      # if there are no ranks
      #   
      # }
      
      #' create a temp function that will calc ebci pvals while handling errors + save results somewhere
      #' 
      temp_pval_function <- function(cur_ebci_params, cur_shrinkage_results, save_folder=NULL) {
        cur_pvals = tryCatch(
          expr = {
            get_ebci_pvals_by_shrinkage_results(cur_ebci_params, cur_shrinkage_results, 
                                                alpha_threshold=ALPHA_THRESHOLD, maxiter=MAXITER, 
                                                parallel=PARALLEL)
          },
          error = function(e) {print(sprintf('      Errored at: %s - %s', split_type, approx_method))}
        )
        
        # No : save the pvals as the entire dataframe so its easier to load (e.g. no reassembling everything)
        # Yes: save the pvals as one vector (can just append/reassemble later then)
        if(!is.null(save_folder) && dir.exists(save_folder)) {
          saveRDS(cur_pvals, sprintf('%s/ebci_pvals.rds', save_folder))
        }
        
        
        return(cur_pvals)
      }
      
      if(replogleShrinkageParams$matapprox_methods_hasranks[[approx_method]]) { # if there are ranks
        # force_calc_all_ranks = calc for all ranks, if false, can still calc missing ranks
        if(approx_method %in% names(ebci_pvals[[split_type]])) { # if already exists, don't add
          force_calc_all_ranks = FALSE
        } else {
          ebci_pvals[[split_type]][[approx_method]] = list()
          force_calc_all_ranks = TRUE
        }
       
        for(r in replogleShrinkageParams$ranks) {
          
          if(force_calc_all_ranks ||
             (!is.vector(ebci_pvals[[split_type]][[approx_method]][[r]])) || 
             (length(ebci_pvals[[split_type]][[approx_method]][[r]]) == 0)) { # if already exists, skip
            
            print(sprintf("[%s]     + now calc EBCI pvals with %s split type, %s method, %s rank", format( Sys.time(), format = "%F %r"), split_type, approx_method, r))
            ebci_pvals[[split_type]][[approx_method]][[r]] = temp_pval_function(
              cur_ebci_params       = shrink_results[[split_type]][[approx_method]][[r]][['ebci_obj']],
              cur_shrinkage_results = shrink_results[[split_type]][[approx_method]][[r]][['ebci_res']], #[1:10, ],
              save_folder=sprintf('%s/shrinkage/%s/%s/rank=%02.f/', plot_path, split_type, approx_method, r)
            )
          }
          
          
        }
      } else { # if there are no ranks
        if((!is.vector(ebci_pvals[[split_type]][[approx_method]])) || (length(ebci_pvals[[split_type]][[approx_method]]) == 0)) { # if already exists, skip
          print(sprintf("[%s]     + now calc EBCI pvals with %s split type, %s method", format( Sys.time(), format = "%F %r"), split_type, approx_method))
          ebci_pvals[[split_type]][[approx_method]] = temp_pval_function(
            cur_ebci_params       = shrink_results[[split_type]][[approx_method]][['ebci_obj']],
            cur_shrinkage_results = shrink_results[[split_type]][[approx_method]][['ebci_res']], #[1:10, ],
            save_folder=sprintf('%s/shrinkage/%s/%s/', plot_path, split_type, approx_method)
          )
        }
      }
    }
  }
  
  
  
  
  saveRDS(object = ebci_pvals, file = sprintf('%s/EBCI_shrinkage_ebcipvals.rds', replogle_save_path))
  
}

# //////////////////////////////////////////////////////////////////////////////////////////////////
# ================================== ASSEMBLE SHRINK RESULTS WITH EBCI PVALS =======================
# //////////////////////////////////////////////////////////////////////////////////////////////////
# using previously saved shrink_results and ebci_pvals that are nested lists, make a tall dataframe 
if(ASSEMBLE_EBCI_PVALS) {
  print('ASSEMBLE SHRINK RESULTS WITH EBCI PVALS  =================================================')
  
  shrink_results = readRDS(sprintf('%s/EBCI_shrinkage_results.rds'  , replogle_save_path))
  ebci_pvals     = readRDS(sprintf('%s/EBCI_shrinkage_ebcipvals.rds', replogle_save_path))
  
  selected_colnames = c('grna', 'gene', 'shrinkage_point', 'shrunk_value', 'lower_ci', 'upper_ci', 'w_eb')
  
  df = NULL
  for(split_type in names(shrink_results)) {
    print(sprintf('  %s', split_type))
    for(approx_method in names(shrink_results[[split_type]])) {
      print(sprintf('    -%s', approx_method))
      if(replogleShrinkageParams$matapprox_methods_hasranks[[approx_method]]) { # if there are ranks
        df_ = NULL
        for(r in replogleShrinkageParams$ranks) {
          df_r = shrink_results[[split_type]][[approx_method]][[r]][['ebci_res']] |> # [1:10, ] |>
                 dplyr::select(dplyr::all_of(selected_colnames)) |>
                 dplyr::mutate(split_type=split_type, approx_method=approx_method, rank=r, .after = 'gene')
          df_r$ebci_pvals = ebci_pvals[[split_type]][[approx_method]][[r]]
          df_ = rbind(df_, df_r); rm(df_r)
        }
      } else { # if there are no ranks
        df_ = shrink_results[[split_type]][[approx_method]][['ebci_res']] |> # [1:10, ] |>
              dplyr::select(dplyr::all_of(selected_colnames)) |>
              dplyr::mutate(split_type=split_type, approx_method=approx_method, rank=NA, .after = 'gene')
        
        df_$ebci_pvals = ebci_pvals[[split_type]][[approx_method]]
      }
      
      df = rbind(df, df_); rm(df_)
    }
  }
  
  df = df |> dplyr::mutate(gene          = as.factor(gene), 
                           grna          = as.factor(grna),
                           split_type    = as.factor(split_type), 
                           approx_method = as.factor(approx_method))
  
  
  df_unshrunk =  shrink_results[['nosamplesplit']][['zeros']][['ebci_res']] |>
                 dplyr::select(grna, gene, unshrunk_value, se, weight)
  df_unshrunk = df_unshrunk |> dplyr::mutate(gene = as.factor(gene), grna = as.factor(grna))
  
  
  # save as RDS (i think, a lot less space)
  saveRDS(object = df,          file = sprintf('%s/EBCI_shrinkage_dataframe.rds', replogle_save_path)) 
  saveRDS(object = df_unshrunk, file = sprintf('%s/EBCI_unshrunk_dataframe.rds', replogle_save_path)) 
  
  # a = readRDS(sprintf('%s/EBCI_shrinkage_dataframe.rds', replogle_save_path))
  # b = readRDS(sprintf('%s/EBCI_unshrunk_dataframe.rds', replogle_save_path))
  
  # # does converting grna, gene (characters) into factors decrease size? yes
  # object.size(df_unshrunk)  # 2066952 bytes
  # df_unshrunk = df_unshrunk |> dplyr::mutate(gene = as.factor(gene), grna = as.factor(grna))
  # object.size(df_unshrunk)  # 1676424 bytes


}










      
      
      