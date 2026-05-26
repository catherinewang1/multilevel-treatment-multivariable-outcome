# script for creating plots using saved objects from the simulation using EBCI at the cell level
# saved objects from running the notebooks/simEBCICell/sim_EBCI_Cell_reps_script.R
# 
# Run in order:
# -        1_sim_EBCI_Cell_reps_script.R              (run simulation: sim data, estimate, and shrink)
# -        2_sim_EBCI_Cell_reps_ebcipval_script.R     (using saved sim results, calculate the inverted EBCI pvals)
# -        3_sim_EBCI_Cell_reps_summary_df.R          (from the large saves, make a summary df)
# -        4.1_sim_EBCI_Cell_reps_plot_script.R       (plot MSE, Miscoverage, and pval (one and fishers) plots)
# -        4.2_sim_EBCI_Cell_reps_plot_pvals_curve.r  (plot average pval curve qq plot)
# - (THIS) 4.3_sim_EBCI_Cell_reps_plot_CIlength.r     (plot CI length plot)


# make CI length plots
# comparing CI lengths before and after shrinkage (should decrease)
# 







# source('../../utils/matrix_shrinkage.r')
# source('../../utils/matrixPrior_utils.R') # may have a different 'create blocky matrix' version
source('../../utils/simEBCICell_utils.R') 



suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(latex2exp))
ggplot2::theme_set(theme_bw() + theme(plot.title = element_text(hjust = .5), 
                                      strip.background = element_rect(fill = 'white')))


overall_save_folder = '../../plots/simEBCICell/'
setting_names = c('A', 'E')
# setting_names = c('A')
# setting_names = c('E')

setting_settings = list('A' = list(SAMPLE_N = FALSE,  # sample if too large, sample min(5000, original #) 
                                     xylim  = c(0, 1.1),
                                   xbreaks = c(0, .25, .5, .75, 1)),
                        'E' = list(SAMPLE_N = FALSE, # if true, need to fix so that it samples the edges more often
                                      xylim = c(.1, .85), 
                                   xbreaks  = c(0, .25, .5, .75)))


CREATE_PLOTS_INDIVIDUAL      = TRUE # create these plots or not
CREATE_PLOTS_OVERALL         = TRUE # create these plots or not
NUMBER_OF_INDIVIDUAL_PLOTS = 3



plot_height_sep = 4
plot_width_sep = 6
# ///////////////////////////////////////////////////////////////////////////////////
# === CREATE PLOTTING FUNCTION ======================================================
# ///////////////////////////////////////////////////////////////////////////////////

#' Create and save CI length plot
#' @param plot_df (dataframe) must have cols
#'    unshrunk_CIlength ,  shrunk_CIlength already created
#' @param ranks (vector) of integer ranks, used to create nice labels and colors
plot_CIlength <- function(plot_df, 
                          ranks, 
                          save_name=NULL,
                          return_plot=TRUE, # set to TRUE when debugging to see, else can be FALSE to speed up in IDE
                          SAMPLE_N = FALSE,
                          xylim = c(0, 1.3), 
                          xbreaks =  c(0, .25, .5, .75, 1, 1.25),
                          plot_height = 5,
                          plot_width   = 8) {
  
  methodrank_colors = create_color_pallete_nicenames(ranks = ranks)
  methodrank_nicenames_order = c()
  for(cur_method in names(method_nicenames)) { # requires declaration of this list/vector (this is defined later in this file, right before the function methodrank_nicenames is defined)
    for(cur_rank in c(NA, ranks)) {
      methodrank_nicenames_order = c(methodrank_nicenames_order, methodrank_nicenames(method_name = cur_method, rank_ = cur_rank) |> unname())
    }
  }
  method_colors = methodrank_colors[grepl(sprintf('(rank=%s)|(Zero)|(Average)', max(ranks)), names(methodrank_colors))]
  names(method_colors) = gsub(pattern = sprintf(" \\(rank=%s\\)", max(ranks)), replacement = "", x = names(method_colors))
  
  matapprox_methods_nicenames = c( "matcomp_softImpute"         =  "softImpute", 
                                   "matdecomp_svd"              =  "SVD", 
                                   "matdecomp_sparsesvd"        =  "Sparse SVD", 
                                   "spectralbiclust"            =  "Spectral Biclustering", 
                                   "spectralbiclust_threshold"  =  "Spectral BiClustering (w threshold)", 
                                   "zeros"                      =  "Zero", 
                                   "average"                    =  "Average" )
  if(SAMPLE_N) {
    set.seed(12345)
    chosen_grna_gene = plot_df |> select(grna, gene) |> distinct() 
    chosen_grna_gene = chosen_grna_gene |> slice_sample(n = min(5000, nrow(chosen_grna_gene)))
    
    temp_df = merge(chosen_grna_gene, plot_df, by = c('grna', 'gene'), all.x = TRUE, all.y = FALSE)
  } else {
    temp_df = plot_df 
  }
  
  
  temp_df = temp_df |>
    mutate(sim_distn  = factor(sim_distn,  levels = c('pois', 'nb'),                   labels = c('Poisson', 'Negative Binomial')),
           split_type = factor(split_type, levels = c('nosamplesplit', 'samplesplit'), labels = c('No Sample Split', 'Sample Split'))) 
  temp_df$methodrank = mapply(FUN = methodrank_nicenames, method_name =  temp_df$method, rank_ = temp_df$rank) |> unname()
  temp_df$methodrank = factor(temp_df$methodrank, levels = methodrank_nicenames_order)
  temp_df = temp_df |>
    mutate(method = factor(method, levels = names(matapprox_methods_nicenames), labels = matapprox_methods_nicenames))
  
  p = ggplot(temp_df |> 
           mutate(rank = as.factor(sapply(rank, FUN = function(r){if(is.na(r)){return(max(ranks))}else{return(r)}}))),
         aes(x = unshrunk_CIlength, y = shrunk_CIlength,
             group = methodrank,
             color = method, 
             alpha = rank
         )) + 
    # geom_point(key_glyph = 'rect') + 
    geom_abline(aes(slope = 1, intercept = 0), color = 'black', alpha = .7) +
    geom_line(key_glyph = 'rect') +
    coord_cartesian(xlim = xylim,
                    ylim = xylim,
                    expand = c(0, 0)) +
    # coord_cartesian(expand = c(0, 0)) +
    facet_grid(cols = vars(split_type), 
               rows = vars(sim_distn), scales = 'free') +
    scale_color_discrete(palette = method_colors) +
    scale_x_continuous(breaks = xbreaks, 
                       labels = sapply(xbreaks, FUN = function(x){if(x==0){return('0')} else {return(x)}})) +
    # scale_x_continuous(breaks = c(0, .25, .5, .75, 1, 1.25),
    #                    labels = c('0', .25, .5, .75, 1.00, 1.25)) +
    labs(title = 'CI lengths before and after shrinkage', 
         x = 'Unshrunk CI Length', y = 'Shrunk CI Length', 
         color = 'Matrix\nApproximation\nMethod', 
         alpha = 'rank') +
    guides(color = guide_legend(order = 1), 
           alpha = guide_legend(order = 2)) +
    theme(panel.grid.minor = element_blank(), 
          panel.grid.major = element_blank(),
          strip.text = element_text(size = 12),
          plot.title = element_text(size = 16),
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 8),
          legend.key.size = unit(.4, 'cm'))
  
  if(!is.null(save_name)) {
    # change save_name argument when calling this function
    # save_name = sprintf('%s/CIlength', setting_name) # overall averaged CI lengths
    # save_name = sprintf('%s/%s/CIlength', setting_name, cur_rep) # individual rep
    
    ggsave(plot = p, filename = sprintf('%s/%s.pdf', overall_save_folder, save_name), height = plot_height, width = plot_width)
    ggsave(plot = p, filename = sprintf('%s/%s.png', overall_save_folder, save_name), height = plot_height, width = plot_width, dpi = 300)
  }
  
  if(return_plot) {
    return(p)
  }
}







# ///////////////////////////////////////////////////////////////////////////////////
# === CREATE PLOTS  =================================================================
# ///////////////////////////////////////////////////////////////////////////////////



# Make some plots for some individual repetitions (only do a few)
if(CREATE_PLOTS_INDIVIDUAL) { 
  
  for(setting_name in setting_names) {
    # list all the repetition directory names (likely just 1, ...., 20) with the saved dfs
    save_files_rep = list.files(path = sprintf('%s%s', overall_save_folder, setting_name) , 
                                pattern = 'sim_result_ebci_pvals_df.rds', recursive=TRUE, full.names = TRUE)
    save_folders_rep = sub(pattern = 'sim_result_ebci_pvals_df.rds', replacement = '', x = save_files_rep) # just the folder names
    
    for(cur_save_folder_rep in save_folders_rep[1:min(NUMBER_OF_INDIVIDUAL_PLOTS, length(save_folders_rep))]) { # just plot some individual runs (don't need all)
      print(sprintf('[%s] %s', Sys.time(), cur_save_folder_rep))
      
      
      # load in and plot individual CI length plot
      tryCatch(expr = {
        
        cur_rep = gsub(pattern = "\\D", replacement = "", x = cur_save_folder_rep ) |> as.numeric() # extract the numeric vals from path
        sim_results = readRDS(sprintf('%s/sim_results.rds', cur_save_folder_rep)); ALPHA = sim_results$ALPHA; ranks = sim_results$ranks; rm(sim_results) # extract ALPHA, ranks
        
        
        df_one = readRDS(sprintf('%s/sim_result_ebci_pvals_df.rds', cur_save_folder_rep))
        df_one = df_one |> 
                    filter(method != 'matcomp_linearreg') |>
                    mutate(rep = cur_rep, .before = 1) |> 
                    mutate(unshrunk_CIlength = 2*qnorm(1 - ALPHA/2)*se) |> 
                    mutate(shrunk_CIlength = upper_ci - lower_ci)

         
        plot_CIlength(plot_df = df_one, 
                      ranks = ranks, 
                      save_name = sprintf('%s/%s/CIlength', setting_name, cur_rep),
                      SAMPLE_N = setting_settings[[setting_name]]$SAMPLE_N,
                      xylim    = setting_settings[[setting_name]]$xylim,
                      xbreaks = setting_settings[[setting_name]]$xbreaks)
        
        
        
      },  error = function(e) {
        print(sprintf('    ------ Plot Ind Errored at: %s', cur_save_folder_rep))
      })
      
      rm(df_one, cur_rep, ALPHA, ranks) # should load in new every time
      
    }
  }
}


gc()

# Make plots with averaged CIs
if(CREATE_PLOTS_OVERALL) {
  for(setting_name in setting_names) {
    
    # list all the repetition directory names (likely just 1, ...., 20) with the saved dfs
    save_files_rep = list.files(path = sprintf('%s%s', overall_save_folder, setting_name) , 
                                pattern = 'sim_result_ebci_pvals_df.rds', recursive=TRUE, full.names = TRUE)
    save_folders_rep = sub(pattern = 'sim_result_ebci_pvals_df.rds', replacement = '', x = save_files_rep) # just the folder names
    
    df = NULL
    for(cur_save_folder_rep in save_folders_rep) { 
      print(sprintf('[%s] %s', Sys.time(), cur_save_folder_rep))
      
      
      # load in and plot individual CI length plot
      tryCatch(expr = {
        
        cur_rep = gsub(pattern = "\\D", replacement = "", x = cur_save_folder_rep ) |> as.numeric() # extract the numeric vals from path
        sim_results = readRDS(sprintf('%s/sim_results.rds', cur_save_folder_rep)); ALPHA = sim_results$ALPHA; ranks = sim_results$ranks; rm(sim_results) # extract ALPHA, ranks
        
        
        df_one = readRDS(sprintf('%s/sim_result_ebci_pvals_df.rds', cur_save_folder_rep))
        df_one = df_one |> 
          filter(method != 'matcomp_linearreg') |>
          mutate(rep = cur_rep, .before = 1) |> 
          mutate(unshrunk_CIlength = 2*qnorm(1 - ALPHA/2)*se) |> 
          mutate(shrunk_CIlength = upper_ci - lower_ci)
        
        df = rbind(df, df_one)
      },  error = function(e) {
        print(sprintf('    ------ Plot Ind Errored at: %s', cur_save_folder_rep))
      })
      
      rm(df_one, cur_rep, ALPHA, ranks) # should load in new every time
    }
    
    
    
   
    
    
    # ranks
    ranks = df$rank |> unique()
    ranks = ranks[!is.na(ranks)] |> sort()
    
    # take Average CI by AY test
    df_summ = df |> 
      group_by(sim_distn, split_type, method, rank, gene, grna) |> 
      summarize(unshrunk_CIlength = mean(unshrunk_CIlength), 
                shrunk_CIlength = mean(  shrunk_CIlength), 
                count = n(), 
                .groups = 'drop') 
    gc()
    
    if(setting_name == 'A') { # just customize the inputs here
      plot_CIlength(plot_df = df_summ, 
                    ranks = ranks, 
                    save_name = sprintf('%s/CIlength', setting_name), 
                    SAMPLE_N = FALSE,
                    xylim    = c(.075, .5),
                    xbreaks = c(0, .1, .2, .3, .4))
      
      # just plot each individually... then assemble together seems to be the easiest
      # pois, nosamplesplit
      p1 = plot_CIlength(plot_df = df_summ |> filter(sim_distn == 'pois' & split_type == 'nosamplesplit'), 
                    ranks = ranks, 
                    save_name = sprintf('%s/CIlength_pois_nosamplesplit', setting_name), 
                    SAMPLE_N = FALSE,
                    xylim    = c(.05, .4),
                    xbreaks = c(.1, .2, .3, .4), 
                    plot_height = plot_height_sep, plot_width = plot_width_sep)
      # pois, samplesplit
      p2 = plot_CIlength(plot_df = df_summ |> filter(sim_distn == 'pois' & split_type == 'samplesplit'), 
                    ranks = ranks, 
                    save_name = sprintf('%s/CIlength_pois_samplesplit', setting_name), 
                    SAMPLE_N = FALSE,
                    xylim    = c(.07, .65),
                    xbreaks = c(.2, .4, .6), 
                    plot_height = plot_height_sep, plot_width = plot_width_sep)
      # nb, nosamplesplit
      p3 = plot_CIlength(plot_df = df_summ |> filter(sim_distn == 'nb' & split_type == 'nosamplesplit'), 
                    ranks = ranks, 
                    save_name = sprintf('%s/CIlength_nb_nosamplesplit', setting_name), 
                    SAMPLE_N = FALSE,
                    xylim    = c(.07, .5),
                    xbreaks = c(.1, .2, .3, .4, .5), 
                    plot_height = plot_height_sep, plot_width = plot_width_sep)
      # nb, samplesplit
      p4 = plot_CIlength(plot_df = df_summ |> filter(sim_distn == 'nb' & split_type == 'samplesplit'), 
                    ranks = ranks, 
                    save_name = sprintf('%s/CIlength_nb_samplesplit', setting_name), 
                    SAMPLE_N = FALSE,
                    xylim    = c(.16, .7),
                    xbreaks = c(0, .2, .4, .6), 
                    plot_height = plot_height_sep, plot_width = plot_width_sep)
      
      
      
      
      A_grob = gridExtra::arrangeGrob(
        p1 + facet_grid(NULL) + theme(legend.position = 'none', plot.title = element_blank(), axis.title = element_blank()),
        p2 + facet_grid(NULL) + theme(legend.position = 'none', plot.title = element_blank(), axis.title = element_blank()),
        p3 + facet_grid(NULL) + theme(legend.position = 'none', plot.title = element_blank(), axis.title = element_blank()),
        p4 + facet_grid(NULL) + theme(legend.position = 'none', plot.title = element_blank(), axis.title = element_blank()),
        layout_matrix = matrix(c(1, 2, 
                                 3, 4), byrow=TRUE, nrow=2)
      )
      
      ggsave(plot = A_grob, filename = sprintf('%s/%s/CIlength_assembled.png', overall_save_folder, setting_name), height = 7, width = 8, dpi = 300)
      ggsave(plot = A_grob, filename = sprintf('%s/%s/CIlength_assembled.pdf', overall_save_folder, setting_name), height = 7, width = 8)
      
      rm(p1, p2, p3, p4, A_grob)
      
      
    } else if(setting_name == 'E') {
      
      
      plot_CIlength(plot_df = df_summ, 
                    ranks = ranks, 
                    save_name = sprintf('%s/CIlength', setting_name), 
                    SAMPLE_N = FALSE,
                    xylim    = c(.075, .5),
                    xbreaks = c(0, .2, .4))
      
      
      
      # just plot each individually... then assemble together seems to be the easiest
      # pois, nosamplesplit
      
      E_SAMPLE_N = FALSE
      
      # p1 = plot_CIlength(plot_df = df_summ |> filter(sim_distn == 'pois' & split_type == 'nosamplesplit'), 
      #               ranks = ranks, 
      #               save_name = sprintf('%s/CIlength_pois_nosamplesplit', setting_name), 
      #               SAMPLE_N = E_SAMPLE_N,
      #               xylim    = c(0, .35),
      #               xbreaks = c(.1, .2, .3), 
      #               plot_height = plot_height_sep, plot_width = plot_width_sep)
      
      
      p1 = plot_CIlength(plot_df = df_summ |> filter(sim_distn == 'pois' & split_type == 'nosamplesplit'), 
                    ranks = ranks, 
                    save_name = sprintf('%s/CIlength_pois_nosamplesplit', setting_name), 
                    SAMPLE_N = E_SAMPLE_N,
                    xylim    = c(0, .35),
                    xbreaks = c(.1, .2, .3), 
                    plot_height = plot_height_sep, plot_width = plot_width_sep) +
        coord_cartesian(
                        # xlim = c(0, .35),
                        xlim = c(.1, .35),
                        ylim = c(0, .35), 
                        expand = c(0, 0, 0, 0))
      
      # pois, samplesplit
      p2 = plot_CIlength(plot_df = df_summ |> filter(sim_distn == 'pois' & split_type == 'samplesplit'), 
                    ranks = ranks, 
                    save_name = sprintf('%s/CIlength_pois_samplesplit', setting_name), 
                    SAMPLE_N = E_SAMPLE_N,
                    xylim    = c(.1, .5),
                    xbreaks = c(.1, .2, .3, .4, .5), 
                    plot_height = plot_height_sep, plot_width = plot_width_sep)
      # nb, nosamplesplit
      p3 = plot_CIlength(plot_df = df_summ |> filter(sim_distn == 'nb' & split_type == 'nosamplesplit'), 
                    ranks = ranks, 
                    save_name = sprintf('%s/CIlength_nb_nosamplesplit', setting_name), 
                    SAMPLE_N = E_SAMPLE_N,
                    xylim    = c(0, .47),
                    xbreaks = c(.1, .2, .3, .4, .5), 
                    plot_height = plot_height_sep, plot_width = plot_width_sep) +
            coord_cartesian(
              # xlim = c(0, .47),
              xlim = c(.24, .47),
              ylim = c(0, .47), 
              expand = c(0, 0, 0, 0))
      # nb, samplesplit
      p4 = plot_CIlength(plot_df = df_summ |> filter(sim_distn == 'nb' & split_type == 'samplesplit'), 
                    ranks = ranks, 
                    save_name = sprintf('%s/CIlength_nb_samplesplit', setting_name), 
                    SAMPLE_N = E_SAMPLE_N,
                    xylim    = c(.2, .7),
                    xbreaks = c(0, .2, .4, .6), 
                    plot_height = plot_height_sep, plot_width = plot_width_sep)
      
      # could change to different xylim allowed...
      # p + coord_cartesian(xlim = c(.35, .7), ylim = c(.2, .7))
      
      
      E_grob = gridExtra::arrangeGrob(
                    p1 + facet_grid(NULL) + theme(legend.position = 'none', plot.title = element_blank(), axis.title = element_blank()) + geom_abline(aes(intercept = 0, slope = 1), color = 'black',linetype = 'solid', linewidth = 1.25),
                    p2 + facet_grid(NULL) + theme(legend.position = 'none', plot.title = element_blank(), axis.title = element_blank()) + geom_abline(aes(intercept = 0, slope = 1), color = 'black',linetype = 'solid', linewidth = 1.25),
                    p3 + facet_grid(NULL) + theme(legend.position = 'none', plot.title = element_blank(), axis.title = element_blank()) + geom_abline(aes(intercept = 0, slope = 1), color = 'black',linetype = 'solid', linewidth = 1.25),
                    p4 + facet_grid(NULL) + theme(legend.position = 'none', plot.title = element_blank(), axis.title = element_blank()) + geom_abline(aes(intercept = 0, slope = 1), color = 'black',linetype = 'solid', linewidth = 1.25),
                    layout_matrix = matrix(c(1, 2, 
                                             3, 4), byrow=TRUE, nrow=2)
              )
      
      ggsave(plot = E_grob, filename = sprintf('%s/%s/CIlength_assembled.png', overall_save_folder, setting_name), height = 7, width = 8, dpi = 300)
      ggsave(plot = E_grob, filename = sprintf('%s/%s/CIlength_assembled.pdf', overall_save_folder, setting_name), height = 7, width = 8)
      
      
      rm(p1, p2, p3, p4, E_grob)
      
    } else { # use specified settings
      plot_CIlength(plot_df = df_summ, 
                    ranks = ranks, 
                    save_name = sprintf('%s/CIlength', setting_name), 
                    SAMPLE_N = setting_settings[[setting_name]]$SAMPLE_N,
                    xylim    = setting_settings[[setting_name]]$xylim,
                    xbreaks = setting_settings[[setting_name]]$xbreaks)
    }
  
    
    
    
  }
}






# # arranging plots together
# test_df = penguins |> filter(species %in% c('Adelie', 'Gentoo') & island %in% c('Torgersen', 'Biscoe'))
# 
# head(test_df)
# 
# 
# ggplot(test_df |> filter(), 
#        aes(x = bill_len, y = bill_dep)) +
#   geom_point() + 
#   facet_grid(cols = vars(species), rows = vars(island))
# 
# 
# 
# p1 = ggplot(test_df, 
#             aes(x = bill_len, y = bill_dep)) +
#   geom_point() + 
#   facet_grid(cols = vars(species), rows = vars(island))
# p2 = ggplot(test_df, 
#             aes(x = bill_len, y = bill_dep)) +
#   geom_point() + 
#   facet_grid(cols = vars(species), rows = vars(island))
# 
# p3 = ggplot(test_df, 
#             aes(x = bill_len, y = bill_dep)) +
#   geom_point() + 
#   facet_grid(cols = vars(species), rows = vars(island))
# p4 = ggplot(test_df, 
#             aes(x = bill_len, y = bill_dep)) +
#   geom_point() + 
#   facet_grid(cols = vars(species), rows = vars(island))
# 
# 
# gridExtra::arrange.grid()









# ///////////////////////////////////////////////////////////////////////////////////
# === TRASH =========================================================================
# ///////////////////////////////////////////////////////////////////////////////////


if(F) {
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
        
        cur_rep = gsub(pattern = "\\D", replacement = "", x = cur_save_folder_rep ) |> as.numeric() # extract the numeric vals from path
        sim_results = readRDS(sprintf('%s/sim_results.rds', cur_save_folder_rep)); ALPHA = sim_results$ALPHA; ranks = sim_results$ranks; rm(sim_results) # extract ALPHA, ranks
     
        
        df_one = readRDS(sprintf('%s/sim_result_ebci_pvals_df.rds', cur_save_folder_rep))
        df_one = df_one |> mutate(unshrunk_CIlength = 2*qnorm(1 - ALPHA/2)*se) |> relocate(unshrunk_CIlength, .after = last_col()) # add 'original' CI lengths
        df_one = df_one |> mutate(shrunk_CIlength = upper_ci - lower_ci)
        df_one = df_one |> mutate(methodrank = paste0(method, rank)) |> 
                 mutate(split_type    = factor(   split_type,  levels = c('nosamplesplit', 'samplesplit'), labels = c('No Sample Split', 'Sample Split')))
        
        p_CI = ggplot(df_one |> mutate(rank = as.factor(sapply(rank, FUN = function(r){if(is.na(r)){return(30)}else{return(r)}}))),
                      aes(x = unshrunk_CIlength, y = shrunk_CIlength,
                          group = methodrank,
                          color = method, 
                          alpha = rank
                      )) + 
          # geom_point(key_glyph = 'rect') + 
          geom_line(key_glyph = 'rect') +
          geom_abline(aes(slope = 1, intercept = 0)) +
          coord_cartesian(xlim = c(0, 1.3), 
                          ylim = c(0, 1.3),
                          # ylim = c(.2, 1), 
                          expand = c(0, 0)) +
          facet_grid(cols = vars(split_type)) +
          # scale_color_discrete(palette = method_colors) +
          scale_x_continuous(breaks = c(0, .25, .5, .75, 1, 1.25),
                             labels = c('0', .25, .5, .75, 1.00, 1.25)) + 
          labs(title = 'CI lengths before and after shrinkage', 
               x = 'Unshrunk CI Length', y = 'Shrunk CI Length', color = 'Matrix\nApproximation\nMethod') + #, alpha = 'rank')
          theme(panel.grid.minor  = element_blank(), 
                strip.text = element_text(size = 12),
                plot.title = element_text(size = 16),
                legend.title = element_text(size = 10),
                legend.text = element_text(size = 8),
                legend.key.size = unit(.4, 'cm'))
        p_CI

        
        
        
        
      },  error = function(e) {
        print(sprintf('    ------ Errored at: %s', cur_save_folder_rep))
      })
    }
  }
}



# Make plots with averaged CIs??
if(CREATE_PLOTS_OVERALL) {
  for(setting_name in setting_names) {
      
    # collect dataframe?
    # list all the repetition directory names (likely just 1, ...., 20)
    save_files_rep = list.files(path = sprintf('%s%s', overall_save_folder, setting_name) , 
                                pattern = 'sim_result_ebci_pvals.rds', recursive=TRUE, full.names = TRUE)
    save_folders_rep = sub(pattern = 'sim_result_ebci_pvals.rds', replacement = '', x = save_files_rep) # just the folder names
    
    df = NULL
    for(cur_save_folder_rep in save_folders_rep) { # just plot some individual runs (don't need all)
      print(sprintf('[%s] %s', Sys.time(), cur_save_folder_rep))
      
      # create and save ebci pvals for each of the repetitions
      tryCatch(expr = {
        cur_rep = gsub(pattern = "\\D", replacement = "", x = cur_save_folder_rep ) |> as.numeric() # extract the numeric vals from path
        sim_results = readRDS(sprintf('%s/sim_results.rds', cur_save_folder_rep)); ALPHA = sim_results$ALPHA; rm(sim_results) # extract ALPHA 
        
        df_one = readRDS(sprintf('%s/sim_result_ebci_pvals_df.rds', cur_save_folder_rep))
        df_one = df_one |> 
                   filter(method != 'matcomp_linearreg') |>
                   mutate(rep = cur_rep, .before = 1) |> 
                   mutate(unshrunk_CIlength = 2*qnorm(1 - ALPHA/2)*se) |> 
                   mutate(shrunk_CIlength = upper_ci - lower_ci)
        
        
        
        
        
        
        df = rbind(df, df_one)
        rm(df_one, ALPHA) # should load in new every time
      },  error = function(e) {
        print(sprintf('    ------ Errored collecting dataframe at: %s', cur_save_folder_rep))
      })
      
      
      
      
      
    }
    
    
    
    sim_results = readRDS(sprintf('%s/sim_results.rds', save_folders_rep[1])); ALPHA = sim_results$ALPHA; ranks = sim_results$ranks; rm(sim_results)  # get the ranks (should be the same)
    
    
    
    
    # take Average CI by AY test
    df_summ = df |> 
      group_by(sim_distn, split_type, method, rank, gene, grna) |> 
      summarize(unshrunk_CIlength = mean(unshrunk_CIlength), 
                shrunk_CIlength = mean(  shrunk_CIlength), 
                count = n(), 
                .groups = 'drop') 
    
    temp_df = df_summ |>
      mutate(sim_distn  = factor(sim_distn,  levels = c('pois', 'nb'),                   labels = c('Poisson', 'Negative Binomial')),
             split_type = factor(split_type, levels = c('nosamplesplit', 'samplesplit'), labels = c('Full Dataset', 'Sample Split'))) 
    temp_df$methodrank = mapply(FUN = methodrank_nicenames, method_name =  temp_df$method, rank_ = temp_df$rank) |> unname()
    temp_df$methodrank = factor(temp_df$methodrank, levels = methodrank_nicenames_order)
    
    
    # plot
    ggplot(df_summ |> 
             mutate(methodrank = paste0(method, rank)) |> 
             mutate(rank = as.factor(sapply(rank, FUN = function(r){if(is.na(r)){return(30)}else{return(r)}}))),
           aes(x = unshrunk_CIlength, y = shrunk_CIlength,
               group = methodrank,
               color = method, 
               alpha = rank
           )) + 
      # geom_point(key_glyph = 'rect') + 
      geom_line(key_glyph = 'rect') +
      geom_abline(aes(slope = 1, intercept = 0)) +
      coord_cartesian(xlim = c(0, 1.3), 
                      ylim = c(0, 1.3),
                      # ylim = c(.2, 1), 
                      expand = c(0, 0)) +
      facet_grid(cols = vars(split_type)) +
      # scale_color_discrete(palette = method_colors) +
      scale_x_continuous(breaks = c(0, .25, .5, .75, 1, 1.25),
                         labels = c('0', .25, .5, .75, 1.00, 1.25)) + 
      labs(title = 'CI lengths before and after shrinkage', 
           x = 'Unshrunk CI Length', y = 'Shrunk CI Length', color = 'Matrix\nApproximation\nMethod') + #, alpha = 'rank')
      theme(panel.grid.minor  = element_blank(), 
            strip.text = element_text(size = 12),
            plot.title = element_text(size = 16),
            legend.title = element_text(size = 10),
            legend.text = element_text(size = 8),
            legend.key.size = unit(.4, 'cm'))
    
    
    
    
  }
}



}





