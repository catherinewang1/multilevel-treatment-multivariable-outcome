# script for creating and saving ebci pvals using saved objects from the simulation using EBCI at the cell level
# saved objects from running the notebooks/simEBCICell/sim_EBCI_Cell_reps_script.R


suppressPackageStartupMessages(library(future.apply))
plan(multisession, workers = 20)  # or some other plan
source('../../utils/simEBCICell_utils.R') 



# ============================== SETTING A =========================================================


setting_name = 'A'
save_folder = sprintf('../../plots/simEBCICell/%s', setting_name)
sim_results = readRDS(sprintf('%s/1/sim_results.rds', save_folder))

# to hopefully speed up parallelization, clear some objects
ebci_params = sim_results$ebci_params
shrinkage_results = sim_results$shrinkage_results
ranks = sim_results$ranks
rm(sim_results); gc()


create_ebci_pvals_from_save(ebci_params=ebci_params, shrinkage_results=shrinkage_results, ranks=ranks, 
                            save_folder=save_folder, 
                            maxiter=10, alpha_threshold=.001)

rm(setting_name, save_folder, sim_results, ebci_params, shrinkage_results, ranks)

# ============================== SETTING E =========================================================





setting_name = 'A'
save_folder = sprintf('../../plots/simEBCICell/%s', setting_name)
sim_results = readRDS(sprintf('%s/1/sim_results.rds', save_folder))

# to hopefully speed up parallelization, clear some objects
ebci_params = sim_results$ebci_params
shrinkage_results = sim_results$shrinkage_results
ranks = sim_results$ranks
rm(sim_results); gc()


create_ebci_pvals_from_save(ebci_params=ebci_params, shrinkage_results=shrinkage_results, ranks=ranks, 
                            save_folder=save_folder, 
                            maxiter=10, alpha_threshold=.001)
