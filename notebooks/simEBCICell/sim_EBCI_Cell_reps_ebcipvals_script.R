# script for creating and saving ebci pvals using saved objects from the simulation using EBCI at the cell level
# saved objects from running the notebooks/simEBCICell/sim_EBCI_Cell_reps_script.R
# 
# Run in order:
# -        sim_EBCI_Cell_reps_script.R          (actual simulation runs)
# - (THIS) sim_EBCI_Cell_reps_ebcipval_script.R (using saved sim results, find the inverted ebcipvals)
# -        sim_EBCI_Cell_reps_summary_df.R      (from the large saves, make a summary df)
# -        sim_EBCI_Cell_reps_plot_script.R     (making some plots)

suppressPackageStartupMessages(library(future.apply))
plan(multisession, workers = 20)  # or some other plan
source('../../utils/simEBCICell_utils.R') 



overall_save_folder = '../../plots/simEBCICell/'

# ============================== SETTING A & E =========================================================

for(setting_name in c('A', 'E')) {
  
  # list all the repetition directory names (likely just 1, ...., 20)
  save_files_rep = list.files(path = sprintf('%s/%s', overall_save_folder, setting_name) , 
                              pattern = 'sim_results.rds', recursive=TRUE, full.names = TRUE)
  save_folder_rep = sub(pattern = 'sim_results.rds', replacement = '', x = save_files_rep) # just the folder names
  for(cur_save_folder_rep in save_folder_rep) {
    print(sprintf('[%s] %s', Sys.time(), cur_save_folder_rep))
    # create and save ebci pvals for each of the repetitions
    tryCatch(expr = {
      sim_results = readRDS(sprintf('%s/sim_results.rds', cur_save_folder_rep))
      
      # to hopefully speed up parallelization, clear some objects
      ebci_params = sim_results$ebci_params
      shrinkage_results = sim_results$shrinkage_results
      ranks = sim_results$ranks
      rm(sim_results); gc()
      
      
      create_ebci_pvals_from_save(ebci_params=ebci_params, shrinkage_results=shrinkage_results, ranks=ranks, 
                                  save_folder=cur_save_folder_rep, 
                                  maxiter=10, alpha_threshold=.001)
      
      rm(save_files_rep, save_folder_rep, ebci_params, shrinkage_results, ranks); gc(verbose = FALSE)
    },  error = function(e) {
      print(sprintf('    ------ Errored at: %s', cur_save_folder_rep))
    })
  }
}



# 
# 
# # ============================== SETTING E, Rep 1 =========================================================
# 
# 
# 
# 
# 
# setting_name = 'E'
# save_folder = sprintf('../../plots/simEBCICell/%s', setting_name)
# sim_results = readRDS(sprintf('%s/1/sim_results.rds', save_folder))
# 
# # to hopefully speed up parallelization, clear some objects
# ebci_params = sim_results$ebci_params
# shrinkage_results = sim_results$shrinkage_results
# ranks = sim_results$ranks
# rm(sim_results); gc()
# 
# 
# create_ebci_pvals_from_save(ebci_params=ebci_params, shrinkage_results=shrinkage_results, ranks=ranks, 
#                             save_folder=save_folder, 
#                             maxiter=10, alpha_threshold=.001)
# 
# 
# 
