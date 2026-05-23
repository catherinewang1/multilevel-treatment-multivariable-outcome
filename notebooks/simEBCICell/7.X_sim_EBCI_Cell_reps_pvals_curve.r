# script for creating plots using saved objects from the simulation using EBCI at the cell level
# saved objects from running the notebooks/simEBCICell/sim_EBCI_Cell_reps_script.R
# 
# Run in order:
# -        sim_EBCI_Cell_reps_script.R          (actual simulation runs)
# -        sim_EBCI_Cell_reps_ebcipval_script.R (using saved sim results, find the inverted ebcipvals)
#  --- [only required to run up to here] ---
# -        sim_EBCI_Cell_reps_summary_df.R      (from the large saves, make a summary df)
# -        sim_EBCI_Cell_reps_plot_script.R     (making some plots)
# - (THIS) 7.X_sim_EBCI_Cell_reps_pvals_curve.r (combine and make pval plots by averaging curve)
#
#
# 


suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
ggplot2::theme_set(theme_bw() + theme(plot.title = element_text(hjust = .5), 
                                      strip.background = element_rect(fill = 'white')))

# regular qqplot..
# Average them across??
# 
# I think we will have to load in gigantic matrix...



source('../../utils/simEBCICell_utils.R')

overall_save_folder = '../../plots/simEBCICell/'

# some script parameters
CREATE_INDIVIDUAL_PVAL_DATAFRAMES = TRUE # create individual pval dataframes or not, set to FALSE if already done and we can just load these
CREATE_SUMMARY_DATAFRAME          = TRUE # create overall summary dataframe or not- avg mse, miscoverage, pvals over the runs (e.g. across 20 repetitions)



setting_names = c('A', 'E')
# setting_names = c('A')
# setting_names = c('E')

RUN_PVAL_CURVE = TRUE


# plot size
height=6; width=14


if(RUN_PVAL_CURVE) {
  
  # create a summary df
  for(setting_name in setting_names) {
    
    save_files_rep = list.files(path = sprintf('%s%s', overall_save_folder, setting_name) , 
                                pattern = 'sim_result_ebci_pvals.rds', recursive=TRUE, full.names = TRUE)
    save_folders_rep = sub(pattern = 'sim_result_ebci_pvals.rds', replacement = '', x = save_files_rep) # just the folder names
   
    df = NULL
    # for(cur_save_folder_rep in save_folders_rep[1:5]) {
    for(cur_save_folder_rep in save_folders_rep) {
      print(sprintf('[%s] %s', Sys.time(), cur_save_folder_rep))

      # create and save ebci pvals for each of the repetitions
      tryCatch(expr = {
        cur_rep = gsub(pattern = "\\D", replacement = "", x = cur_save_folder_rep ) |> as.numeric() # extract the numeric vals from path
        
        df_one = readRDS(sprintf('%s/sim_result_ebci_pvals_df.rds', cur_save_folder_rep)) |> 
          mutate(rep = cur_rep) |> 
          filter(method != 'matcomp_linearreg') |> # exclude this... this performs badly
          mutate(isTheta0 = (true_theta == 0),
                 isTheta0Named = sapply(FUN = function(b) {if(b){'Null'} else {'Alt'}}, X = isTheta0)) |> 
          group_by(rep, sim_distn, split_type, method, rank, isTheta0Named) |> 
          arrange(ebci_pvals) |>
          mutate(theoretical = (1:n())/n())
        
        
        df = rbind(df, df_one)
      },  error = function(e) {
        print(sprintf('    ------ Errored at: %s', cur_save_folder_rep))
      })
    }
    
    
    # === labels and colors
    # make the ordering for a various number of methodrank levels (even if not used): 
    # first method according to method_nicenames, then rank NA, 1, 2, ...
    sim_results = readRDS(sprintf('%s/sim_results.rds', save_folders_rep[1])); ranks = sim_results$ranks; rm(sim_results)  # get the ranks (should be the same)
    
    methodrank_colors = create_color_pallete_nicenames(ranks = ranks)
    methodrank_nicenames_order = c()
    for(cur_method in names(method_nicenames)) { # requires declaration of this list/vector (this is defined later in this file, right before the function methodrank_nicenames is defined)
      for(cur_rank in c(NA, ranks)) {
        methodrank_nicenames_order = c(methodrank_nicenames_order, methodrank_nicenames(method_name = cur_method, rank_ = cur_rank) |> unname())
      }
    }
    
    
    temp_df =  df |> 
      group_by(sim_distn, split_type, method, rank, isTheta0Named, theoretical) |> 
      summarize(avg_ebci_pval_curve = mean(ebci_pvals), .groups = 'drop') |>
      mutate(sim_distn  = factor(sim_distn,  levels = c('pois', 'nb'),                   labels = c('Poisson', 'Negative Binomial')),
             split_type = factor(split_type, levels = c('nosamplesplit', 'samplesplit'), labels = c('Full Dataset', 'Sample Split'))) 
    
    temp_df$methodrank = mapply(FUN = methodrank_nicenames, method_name =  temp_df$method, rank_ = temp_df$rank) |> unname()
    temp_df$methodrank = factor(temp_df$methodrank, levels = methodrank_nicenames_order)
    
    
    
    
    
    ggplot(temp_df, 
           aes(x = theoretical, y = avg_ebci_pval_curve,
               group = methodrank, color = methodrank, fill = methodrank)) +
      geom_abline(aes(slope = 1, intercept = 0), color = 'black', linewidth = 1) +
      geom_line(alpha = .8, linewidth = .8, key_glyph = 'rect') +
      coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
      labs(title = 'Averaged QQ-plot of Inverted EBCI p-values vs Unif(0,1)') + 
      scale_color_discrete(palette = methodrank_colors[names(methodrank_colors) %in% temp_df$methodrank]) +
      facet_grid(rows = vars(sim_distn), cols = vars(split_type, isTheta0Named), scales = "fixed") +
      theme(panel.grid.major.x = element_blank(), strip.background = element_rect(fill = NA))
    
    
    ggsave(filename = sprintf('%s/ebci_pvals_avgqqcurve.pdf', sprintf('%s%s/', overall_save_folder, setting_name)), width = width, height = height)
    
    # subsampled version
    set.seed(12345)
    ggplot(temp_df |> group_by(sim_distn, split_type, method, rank) |> sample_frac(size = .25), # E: 25000 max (grna-gene pairs) , 
           aes(x = theoretical, y = avg_ebci_pval_curve,
               group = methodrank, color = methodrank, fill = methodrank)) +
      geom_abline(aes(slope = 1, intercept = 0), color = 'black', linewidth = 1) +
      geom_line(alpha = .8, linewidth = .8, key_glyph = 'rect') +
      coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
      labs(title = 'Averaged QQ-plot of Inverted EBCI p-values vs Unif(0,1)') + 
      scale_color_discrete(palette = methodrank_colors[names(methodrank_colors) %in% temp_df$methodrank]) +
      facet_grid(rows = vars(sim_distn), cols = vars(split_type, isTheta0Named), scales = "fixed") +
      theme(panel.grid.major.x = element_blank(), strip.background = element_rect(fill = NA))
    
    
    ggsave(filename = sprintf('%s/ebci_pvals_avgqqcurve_sampled.pdf', sprintf('%s%s/', overall_save_folder, setting_name)), width = width, height = height)
    
    
    
    rm(df, ranks)
  }
}





# ///////////////////////////////////////////////////////////////////////////////////////////////////////
# ============== TRASH ==================================================================================
# ///////////////////////////////////////////////////////////////////////////////////////////////////////
if(F) {
  
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


save_folder = sprintf('%s%s/', overall_save_folder, setting_name)
df_one = readRDS(sprintf('%s/sim_result_ebci_pvals_df.rds', cur_save_folder_rep)) |> mutate(rep = cur_rep, .before = 'sim_distn')


# prep dataframe
# add qqplot theoretical values / 'rejection_rate' /
temp_df = df_one |> 
  filter(method != 'matcomp_linearreg') |> # exclude this... this performs badly
  mutate(isTheta0 = (true_theta == 0),
         isTheta0Named = sapply(FUN = function(b) {if(b){'Null'} else {'Alt'}}, X = isTheta0)) |> 
  group_by(rep, sim_distn, split_type, method, rank, isTheta0Named) |> 
  arrange(ebci_pvals) |>
  mutate(theoretical = (1:n())/n()) |> 
  mutate(sim_distn  = factor(sim_distn,  levels = c('pois', 'nb'),                   labels = c('Poisson', 'Negative Binomial')),
         split_type = factor(split_type, levels = c('nosamplesplit', 'samplesplit'), labels = c('Full Dataset', 'Sample Split'))) 

temp_df$methodrank = mapply(FUN = methodrank_nicenames, method_name =  temp_df$method, rank_ = temp_df$rank) |> unname()
temp_df$methodrank = factor(temp_df$methodrank, levels = methodrank_nicenames_order)


ggplot(temp_df,
       aes(x = theoretical, y = ebci_pvals,
           group = methodrank, color = methodrank, fill = methodrank)) +
  geom_abline(aes(slope = 1, intercept = 0)) +
  geom_line(alpha = .8, linewidth = .8) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(title = 'Averaged QQ-plot of Inverted EBCI p-values vs Unif(0,1)') + 
  scale_color_discrete(palette = methodrank_colors[names(methodrank_colors) %in% temp_df$methodrank]) +
  facet_grid(rows = vars(sim_distn), cols = vars(split_type, isTheta0Named), scales = "fixed") +
  theme(panel.grid.major.x = element_blank(), strip.background = element_rect(fill = NA))

ggsave(filename = sprintf('%s/ebci_pvals_avgqqplotcurve.pdf', save_folder), width = width, height = height) 





# ggplot(data.frame(sample = runif(100, 0, .1)), 
#        aes(sample=sample)) +
#   geom_qq(distribution = stats::qunif)
# 
# 


df = NULL
for(cur_save_folder_rep in save_folders_rep[1:5]) {
  # for(cur_save_folder_rep in save_folders_rep) {
  print(sprintf('[%s] %s', Sys.time(), cur_save_folder_rep))
  # sim_results = readRDS(sprintf('%s/sim_results.rds', cur_save_folder_rep)); ALPHA = sim_results$ALPHA; ranks = sim_results$ranks; rm(sim_results) # extract ALPHA 
  
  # create and save ebci pvals for each of the repetitions
  tryCatch(expr = {
    cur_rep = gsub(pattern = "\\D", replacement = "", x = cur_save_folder_rep ) |> as.numeric() # extract the numeric vals from path
    
    
    df_one = readRDS(sprintf('%s/sim_result_ebci_pvals_df.rds', cur_save_folder_rep)) |> 
      mutate(rep = cur_rep) |> 
      filter(method != 'matcomp_linearreg') |> # exclude this... this performs badly
      mutate(isTheta0 = (true_theta == 0),
             isTheta0Named = sapply(FUN = function(b) {if(b){'Null'} else {'Alt'}}, X = isTheta0)) |> 
      group_by(rep, sim_distn, split_type, method, rank, isTheta0Named) |> 
      arrange(ebci_pvals) |>
      mutate(theoretical = (1:n())/n())
    
    
    df = rbind(df, df_one)
  },  error = function(e) {
    print(sprintf('    ------ Errored at: %s', cur_save_folder_rep))
  })
}



temp_df =  df |> 
  group_by(sim_distn, split_type, method, rank, isTheta0Named, theoretical) |> 
  summarize(avg_ebci_pval_curve = mean(ebci_pvals), .groups = 'drop') |>
  mutate(sim_distn  = factor(sim_distn,  levels = c('pois', 'nb'),                   labels = c('Poisson', 'Negative Binomial')),
         split_type = factor(split_type, levels = c('nosamplesplit', 'samplesplit'), labels = c('Full Dataset', 'Sample Split'))) 

temp_df$methodrank = mapply(FUN = methodrank_nicenames, method_name =  temp_df$method, rank_ = temp_df$rank) |> unname()
temp_df$methodrank = factor(temp_df$methodrank, levels = methodrank_nicenames_order)


ggplot(temp_df, 
       aes(x = theoretical, y = avg_ebci_pval_curve,
           group = methodrank, color = methodrank, fill = methodrank)) +
  geom_abline(aes(slope = 1, intercept = 0), color = 'black', linewidth = 1) +
  geom_line(alpha = .8, linewidth = .8, key_glyph = 'rect') +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(title = 'Averaged QQ-plot of Inverted EBCI p-values vs Unif(0,1)') + 
  scale_color_discrete(palette = methodrank_colors[names(methodrank_colors) %in% temp_df$methodrank]) +
  facet_grid(rows = vars(sim_distn), cols = vars(split_type, isTheta0Named), scales = "fixed") +
  theme(panel.grid.major.x = element_blank(), strip.background = element_rect(fill = NA))


ggsave(filename = sprintf('%s/ebci_pvals_averageqqplotcurve.pdf', save_folder), width = width, height = height) 


}





