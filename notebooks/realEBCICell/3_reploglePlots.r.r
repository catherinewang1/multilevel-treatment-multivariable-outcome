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


