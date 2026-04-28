# using previously saved sim_results and the ebci_pvals
# 1. creates dataframes with shrinkage results (individual for each rep + summary using fishers)
# 2. creates ebci_pval plots
#
# Run in order:
# -        sim_EBCI_Cell_reps_script.R          (actual simulation runs)
# -        sim_EBCI_Cell_reps_ebcipval_script.R (using saved sim results, find the inverted ebcipvals)
# - (THIS) sim_EBCI_Cell_reps_summary_df.R      (from the large saves, make a summary df)
# -        sim_EBCI_Cell_reps_plot_script.R     (making some plots)

source('../../utils/simEBCICell_utils.R') 

overall_save_folder = '../../plots/simEBCICell/'

# some script parameters
CREATE_INDIVIDUAL_PVAL_DATAFRAMES = TRUE # create individual pval dataframes or not, set to FALSE if already done and we can just load these
CREATE_SUMMARY_DATAFRAME          = TRUE # create overall summary dataframe or not- avg mse, miscoverage, pvals over the runs (e.g. across 20 repetitions)



setting_names = c('A', 'E')


if(CREATE_INDIVIDUAL_PVAL_DATAFRAMES) {
  
  # format the pvals from nested list of matrices into a dataframe
  for(setting_name in setting_names) {
    # list all the repetition directory names (likely just 1, ...., 20)
    save_files_rep = list.files(path = sprintf('%s%s', overall_save_folder, setting_name) , 
                                pattern = 'sim_result_ebci_pvals.rds', recursive=TRUE, full.names = TRUE)
    save_folders_rep = sub(pattern = 'sim_result_ebci_pvals.rds', replacement = '', x = save_files_rep) # just the folder names
    
    for(cur_save_folder_rep in save_folders_rep) {
      print(sprintf('[%s] %s', Sys.time(), cur_save_folder_rep))
      
      # create and save ebci pvals for each of the repetitions
      tryCatch(expr = {
        sim_results           = readRDS(sprintf('%s/sim_results.rds', cur_save_folder_rep))
        sim_result_ebci_pvals = readRDS(sprintf('%s/sim_result_ebci_pvals.rds', cur_save_folder_rep))
        
        sim_result_ebci_pvals_df = create_ebcipvals_df(sim_results=sim_results, sim_result_ebci_pvals=sim_result_ebci_pvals)
        
        
        
        # writing this df as csv (do not, for disk space reasons)
        # write.csv(x = sim_result_ebci_pvals_df, file = sprintf('%s/sim_result_ebci_pvals_df.csv', cur_save_folder_rep))
        # or
        # writing this df as rds (less space? yes much less, eg for setting A,  570KB vs 92 KB)
        saveRDS(object = sim_result_ebci_pvals_df, file = sprintf('%s/sim_result_ebci_pvals_df.rds', cur_save_folder_rep))
        
        
        rm(sim_results, sim_result_ebci_pvals, sim_result_ebci_pvals_df); gc(verbose = FALSE)
      },  error = function(e) {
        print(sprintf('    ------ Errored at: %s', cur_save_folder_rep))
      })
    }
  }
  
  
  
  
  
} 




if(CREATE_SUMMARY_DATAFRAME) {
  
  # create a summary df
  for(setting_name in c('A', 'E')) {
    # list all the repetition directory names (likely just 1, ...., 20)
    save_files_rep = list.files(path = sprintf('%s%s', overall_save_folder, setting_name) , 
                                pattern = 'sim_result_ebci_pvals.rds', recursive=TRUE, full.names = TRUE)
    save_folders_rep = sub(pattern = 'sim_result_ebci_pvals.rds', replacement = '', x = save_files_rep) # just the folder names
    df = NULL
    # for(cur_save_folder_rep in save_folders_rep[1:5]) {
    for(cur_save_folder_rep in save_folders_rep) {
      print(sprintf('[%s] %s', Sys.time(), cur_save_folder_rep))
      sim_results = readRDS(sprintf('%s/sim_results.rds', cur_save_folder_rep)); ALPHA = sim_results$ALPHA; rm(sim_results) # extract ALPHA 
      
      # create and save ebci pvals for each of the repetitions
      tryCatch(expr = {
        cur_rep = gsub(pattern = "\\D", replacement = "", x = cur_save_folder_rep ) |> as.numeric() # extract the numeric vals from path
        df = rbind(df, readRDS(sprintf('%s/sim_result_ebci_pvals_df.rds', cur_save_folder_rep)) |> mutate(rep = cur_rep)) # |> rename(method = approx_method)
      },  error = function(e) {
        print(sprintf('    ------ Errored at: %s', cur_save_folder_rep))
      })
    }
    
    
    cur_df_summary = create_summary_bytest_df(df = df, ALPHA = ALPHA, save_folder = sprintf('%s%s', overall_save_folder, setting_name))
    
    rm(df, cur_df_summary)
    
    
    # dim(cur_df_summary); head(cur_df_summary)
    # sapply(X = cur_df_summary$grna[1:5], FUN = function(x) {gsub(pattern = "\\D", replacement = "", x = x) |> as.numeric()}) # how to get the x or y idx (instead of the character gene and grna)
    # # check the saved file
    # test_df = readRDS(sprintf('%s%s/sim_result_ebci_summary.rds', overall_save_folder, setting_name))
    # dim(test_df); head(test_df)
    # ggplot(test_df |> group_by(sim_distn, split_type, method, rank) |> summarize(mse = mean(mse)) |> mutate(methodrank=paste0(method, rank)),
    #        aes(x = methodrank, y = mse)) +
    #   geom_col() +
    #   facet_grid(rows = vars(sim_distn), cols = vars(split_type), scales = 'free_y') +
    #   theme(# axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
    #         axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    #         panel.grid.major.x = element_blank(), strip.background = element_rect(fill = NA),
    #         axis.ticks.x = element_blank())
    
  }
}









