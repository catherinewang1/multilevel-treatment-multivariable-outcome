# script for creating plots using saved objects from the simulation using EBCI at the cell level
# saved objects from running the notebooks/simEBCICell/sim_EBCI_Cell_reps_script.R
# 
# Run in order:
# -        sim_EBCI_Cell_reps_script.R          (actual simulation runs)
# -        sim_EBCI_Cell_reps_ebcipval_script.R (using saved sim results, find the inverted ebcipvals)
# -        sim_EBCI_Cell_reps_summary_df.R      (from the large saves, make a summary df)
# - (THIS) sim_EBCI_Cell_reps_plot_script.R     (making some plots)

suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(latex2exp))

source('../../utils/matrix_shrinkage.r')
source('../../utils/matrixPrior_utils.R') # may have a different 'create blocky matrix' version
source('../../utils/simEBCICell_utils.R') 


overall_save_folder = '../../plots/simEBCICell/'
setting_names = c('A', 'E')


CREATE_PLOTS_INDIVIDUAL      = TRUE # create these plots or not
CREATE_PLOTS_OVERALL         = TRUE # create these plots or not

NUMBER_OF_INDIVIDUAL_PLOTS = 3




# Make some plots for some individual repetitions (only do a few)
if(CREATE_PLOTS_INDIVIDUAL) { 
  
  for(setting_name in setting_names) {
    # list all the repetition directory names (likely just 1, ...., 20)
    save_files_rep = list.files(path = sprintf('%s%s', overall_save_folder, setting_name) , 
                                pattern = 'sim_result_ebci_pvals.rds', recursive=TRUE, full.names = TRUE)
    save_folders_rep = sub(pattern = 'sim_result_ebci_pvals.rds', replacement = '', x = save_files_rep) # just the folder names
    
    for(cur_save_folder_rep in save_folders_rep[1:min(NUMBER_OF_INDIVIDUAL_PLOTS, length(save_folders_rep))]) { # just plot some individual runs (don't need all)
      print(sprintf('[%s] %s', Sys.time(), cur_save_folder_rep))
      
      # create and save ebci pvals for each of the repetitions
      tryCatch(expr = {
        
        
        cur_plot_specs = list(mse=list(height=5, width=8), miscoverage=list(height=5, width=8))
        if(setting_name == 'E') {cur_plot_specs$matrix_individual = list(color_limits = c(-.5, 1.5))}
        
        make_plots_from_save(sim_results_filenames=sprintf('%s/sim_results.rds', cur_save_folder_rep),
                            save_folder=cur_save_folder_rep,
                            which_plots = list(matrix=FALSE, mse=TRUE, miscoverage=TRUE, matrix_individual=TRUE),
                            write_plot_df = FALSE, write_summary_df=TRUE,
                            plot_specs =cur_plot_specs)
      
      },  error = function(e) {
        print(sprintf('    ------ Errored at: %s', cur_save_folder_rep))
      })
    }
  }
}



# Make plots using previously saved overall summary stats
if(CREATE_PLOTS_OVERALL) {
  
  
  
  for(setting_name in setting_names) {
    print(sprintf('[%s] %s', Sys.time(), cur_save_folder_rep))
    
    sim_results = readRDS(sprintf('%s%s/1/sim_results.rds', overall_save_folder, setting_name)) # try to make sure rep 1 is present, need ranks, ALPHA
    summary_df  = readRDS(sprintf('%s%s/sim_result_ebci_summary.rds', overall_save_folder, setting_name))
    
    h_plot_miscoverage(plot_df = summary_df, 
                       ranks   = sim_results$ranks, 
                       ALPHA   = sim_results$ALPHA, 
                       save_folder = sprintf('%s%s/', overall_save_folder, setting_name), 
                       height=7, width=12, 
                       save_ggplot=TRUE, plot_df_is_summ=TRUE) 
    
    h5_3_plots_mse(plot_df = summary_df, 
                   ranks   = sim_results$ranks, 
                   save_folder = sprintf('%s%s/', overall_save_folder, setting_name), 
                   height=7, width=12, 
                   save_ggplot=FALSE, plot_df_is_summ=TRUE)
    
    
    
    
    # methodrank_colors = create_color_pallete_nicenames(ranks = sim_results$ranks)
    # # make the ordering for a various number of methodrank levels (even if not used): first method according to method_nicenames, then rank NA, 1, 2, ...
    # methodrank_nicenames_order = c()
    # for(cur_method in names(method_nicenames)) { # requires declaration of this list/vector (this is defined later in this file, right before the function methodrank_nicenames is defined)
    #   for(cur_rank in c(NA, sim_results$ranks)) {
    #     methodrank_nicenames_order = c(methodrank_nicenames_order, methodrank_nicenames(method_name = cur_method, rank_ = cur_rank) |> unname())
    #   }
    # }
    # plot_df_summ = summary_df |> 
    #   filter(method != 'matcomp_linearreg') |> # exclude this... this performs badly
    #   filter(true_theta != 0) |>
    #   group_by(sim_distn, split_type, method, rank) |> 
    #   summarize(mse = mean(mse), .groups = 'drop') |> 
    #   mutate(sim_distn  = factor(sim_distn,  levels = c('pois', 'nb'),                   labels = c('Poisson', 'Negative Binomial')),
    #          split_type = factor(split_type, levels = c('nosamplesplit', 'samplesplit'), labels = c('Full Dataset', 'Sample Split')))
    # 
    # 
    # summary_df |> filter(method == 'unshrunk')
    # 
    # plot_df_summ$methodrank = mapply(FUN = methodrank_nicenames, method_name =  plot_df_summ$method, rank_ = plot_df_summ$rank) |> unname()
    # plot_df_summ$methodrank = factor(plot_df_summ$methodrank, levels = methodrank_nicenames_order)
    # 
    # # take the mse from unshrunk on Full Dataset
    # original_all_mse = plot_df_summ |> filter(method == 'unshrunk' & split_type == 'Full Dataset')
    # original_all_mse = rbind(original_all_mse, original_all_mse |> mutate(split_type = 'Sample Split')) # repeat but change split type for plotting
    # 
    # 
    # 
    # 
    # ggplot() +
    #   geom_col(data = plot_df_summ, aes(x = methodrank, y = mse, fill = methodrank), color = 'black', alpha = 1) + # MSEs
    #   geom_hline(data = original_all_mse, aes(yintercept = mse), color = '#DB1A1A', linewidth = .7, alpha = .6) +  # MSEs for 'naive'/original method
    #   facet_grid(rows = vars(sim_distn), cols = vars(split_type), scales = 'free_y', 
    #              # labeller = as_labeller(c(sim_distn_nicename, split_type_nicename))
    #              labeller = label_value
    #   ) +
    #   scale_fill_discrete(palette = methodrank_colors[names(methodrank_colors) %in% plot_df_summ$methodrank]) +
    #   labs(title = 'MSE', x = 'method + rank', y = 'mse', fill = 'Method') +
    #   scale_y_continuous(expand = expansion(mult = c(0, .05))) +
    #   theme(# axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5), 
    #     axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
    #     panel.grid.major.x = element_blank(), strip.background = element_rect(fill = NA), 
    #     axis.ticks.x = element_blank())
    # 
    # 
    # 
    # 
    # 
    }
}






# === OLD ==========================================================================================
# setting_name = 'A'
# save_folder = sprintf('../../plots/simEBCICell/%s', setting_name)
# # sim_results = readRDS(sprintf('%s/1/sim_results.rds', save_folder))
# 
# sim_results_filenames = grep(pattern = 'sim_results.rds', x = list.files(save_folder, recursive = TRUE, full.names = TRUE), value = TRUE)
# make_plots_from_save(sim_results_filenames=sim_results_filenames,
#                      save_folder=save_folder,
#                      which_plots = list(matrix=FALSE, mse=TRUE, miscoverage=TRUE),
#                      write_plot_df = FALSE, write_summary_df=TRUE, 
#                      plot_specs =list(mse=list(height=5, width=8), miscoverage=list(height=5, width=8)))
# 
# setting_name = 'E'
# save_folder = sprintf('../../plots/simEBCICell/%s', setting_name)
# # sim_results = readRDS(sprintf('%s/1/sim_results.rds', save_folder))
# 
# sim_results_filenames = grep(pattern = 'sim_results.rds', x = list.files(save_folder, recursive = TRUE, full.names = TRUE), value = TRUE)
# make_plots_from_save(sim_results_filenames=sim_results_filenames,
#                      save_folder=save_folder,
#                      which_plots = list(matrix=FALSE, mse=TRUE, miscoverage=TRUE),
#                      write_plot_df = FALSE, write_summary_df=TRUE, 
#                      plot_specs =list(mse=list(height=6, width=12), miscoverage=list(height=6, width=12)))


# sim_results = readRDS(sprintf('%s/sim_results.rds', "../../plots/simEBCICell/A/1/"))
# sim_results$estimate_matapprox[['pois']][['train']][['spectralbiclust']][[1]]
# ggsave(NULL, filename = sprintf('%s/test.pdf', cur_save_folder_rep))