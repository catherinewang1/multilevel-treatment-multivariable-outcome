# Using previously saved objects, make some final shrinkage CI and point visualization plots 
# Run in order:
# 1_...
# 2_...
# 3_...
# ...
# 4.2_reploglePlots_methoddemonstration.R
# 4.3_reploglePlots_showcaseMatrices.R
# 4.4_reploglePlots_shrinkViz.r
# 
# Require previously saved objects at:
# 
# Will save plots at: 
#    <plot_path>/final/method_demonstration/
#                                          /samplesplit
#                                          /nosamplesplit
# Then, I/you should make a figure using these saved plots demonstrating the
# EBCI method. e.g. cells --> SCEPTRE est --> mat approx --> shrink 
#   



# libraries/packages
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
ggplot2::set_theme(theme_bw() + 
                     theme(plot.title = element_text(hjust = .5), 
                           plot.subtitle = element_text(hjust = .5),
                           strip.background = element_rect(fill = 'white')))





# /////////////////////////////////////////////////////////////////////////////////////////////////////
# ================== Setup: Script Params ========================================================================
# /////////////////////////////////////////////////////////////////////////////////////////////////////
# here, set the script params

# sceptre_save_path = '../../saves/sceptre/replogle/' # location where sceptre results are located
replogle_save_path= '../../saves/replogle/'         # location to save replogle approximations/shrinkage/etc. results
plot_path         = '../../plots/replogle/EBCI/'    # location to save plots of EBCI analysis on replogle dataset

# limit the upper and lower values
# MY_COLOR_LIMITS = c(-3, 4)
# MY_COLOR_LIMITS = c(-2.5, 2.5)
MY_COLOR_LIMITS = c(-2, 2) # default color limits if not specified


myRed  = colorRampPalette(RColorBrewer::brewer.pal(n = 7, name = "RdBu"))(7)[1]
myBlue = colorRampPalette(RColorBrewer::brewer.pal(n = 7, name = "RdBu"))(7)[7]


SUBSAMPLE_N = 10 # number of AY (grna-gene) tests to display


# Select some settings to showcase:
selectedSettings = data.frame(matrix(
  c(# settings for sample splitting
    'samplesplit',                         'SVD',  3,
    'samplesplit',                   'sparseSVD', 20,
    'samplesplit',        'sparseSVD_autoparams', 20,
    'samplesplit',                  'softImpute',  5,
    'samplesplit',             'spectralbiclust',  5,
    'samplesplit',   'spectralbiclust_threshold',  5,
    'samplesplit',                     'average', NA,
    'samplesplit',                       'zeros', NA,
    # settings for no sample splitting
    'nosamplesplit',                         'SVD',  3,
    'nosamplesplit',                   'sparseSVD', 20,
    'nosamplesplit',        'sparseSVD_autoparams', 20,
    'nosamplesplit',                  'softImpute',  5,
    'nosamplesplit',             'spectralbiclust',  5,
    'nosamplesplit',   'spectralbiclust_threshold',  5,
    'nosamplesplit',                     'average', NA,
    'nosamplesplit',                       'zeros', NA),
  ncol = 3, byrow = TRUE)) 
colnames(selectedSettings) = c('split_type', 'approx_method', 'rank')
selectedSettings$rank = as.numeric(selectedSettings$rank)
selectedSettings



# /////////////////////////////////////////////////////////////////////////////////////////////////////
# ================== Load =============================================================================
# /////////////////////////////////////////////////////////////////////////////////////////////////////

# previously saved objects
replogleShrinkageParams = readRDS(file = sprintf('%s/replogleShrinkageParams.rds', replogle_save_path))

gene_index = read.csv(sprintf('%s/gene_index.csv', replogle_save_path))
grna_index = read.csv(sprintf('%s/grna_index.csv', replogle_save_path))

df_shrunk   = readRDS(sprintf('%s/EBCI_shrinkage_dataframe.rds', replogle_save_path)) 
df_unshrunk = readRDS(sprintf('%s/EBCI_unshrunk_dataframe.rds', replogle_save_path))  
df_unshrunk = merge(merge(df_unshrunk, gene_index, by = 'gene'), grna_index, by = 'grna') 
df_unshrunk = df_unshrunk |> dplyr::mutate(gene = as.factor(gene), grna = as.factor(grna))


# construct df combining results for plotting
df = merge(df_unshrunk, df_shrunk, by = c('grna', 'gene', 'split_type'), all.y = TRUE)

dim(df); dim(df_shrunk) # should have same number of rows 
head(df); head(df_shrunk) # should have added cols: gene_idx, grna_idx, unshrunk_value, se, weight
gc()




# get the original values from train dataset- (does not appear in unshrunk or shrunk dfs...)

# load in as matrix form:
estse_matrices = readRDS(sprintf('%s/estse_matrices.rds'  , replogle_save_path))
# df_train = reshape2::melt(estse_matrices$est_matrices$train , varnames = c('grna', 'gene'), value.name = 'unshrunk_value')
# df_train_se = reshape2::melt(estse_matrices$se_matrices$train , varnames = c('grna', 'gene'), value.name = 'se')
df_train = merge(reshape2::melt(estse_matrices$est_matrices$train , varnames = c('grna', 'gene'), value.name = 'unshrunk_value'), 
                 reshape2::melt(estse_matrices$se_matrices$train , varnames = c('grna', 'gene'), value.name = 'se'), 
                 by = c('grna', 'gene'))
df_train = merge(merge(df_train, gene_index, by = 'gene'), grna_index, by = 'grna') |>  # add in the previously saved orderings
  dplyr::mutate(gene = as.factor(gene), grna = as.factor(grna))                # make gene and grna into factors
rm(estse_matrices); gc()





split_type_nicename = list('nosamplesplit' = 'No Sample Split', 
                           'samplesplit' = 'Sample Split')
approx_method_nicename = list('SVD' = 'SVD',
                              'sparseSVD' = 'Sparse SVD',
                              'sparseSVD_autoparams' = 'Sparse SVD\n(CV params)',
                              'softImpute' = 'Soft Impute', 
                              'spectralbiclust' = 'Spectral Biclustering', 
                              'spectralbiclust_threshold' = 'Spectral Biclustering\n(threshold)',
                              'average' = 'Average', 
                              'zeros'   = 'Zeros')




set.seed(12345)
# get the same gene/grna tests across all methods and ranks 
sample_tests = df |> select(grna, gene) |> distinct() |> slice_sample(n = SUBSAMPLE_N) |> mutate(x = 1:n())
df_sample = merge(df, sample_tests, all.x = FALSE, all.y = TRUE)
df_sample = merge(df_sample, selectedSettings, by = c('split_type', 'approx_method', 'rank'), all.x = FALSE, all.y = TRUE)
# df_sample = df_sample |> mutate(split_type = factor(split_type, levels = names(split_type_nicename), labels = split_type_nicename),
#                                           approx_method = factor(approx_method, levels = names(approx_method_nicename), labels = approx_method_nicename))
# df_sample |> group_by(split_type, approx_method, rank) |> summarize(count = n()) # should all be nrow(sample_tests)




# /////////////////////////////////////////////////////////////////////////////////////////////////////
# ================== Setup Directory & Fn =============================================================
# /////////////////////////////////////////////////////////////////////////////////////////////////////
dir.create(sprintf('%s/final/', plot_path))
dir.create(sprintf('%s/final/shrinkViz/', plot_path))    # 
dir.create(sprintf('%s/final/shrinkViz/individual', plot_path)) # create subfolder for these, can create keynote fig outside this

plot_heights = 3
plot_widths  = 5
plot_dpi     = 600

# NUM_GRNA_DISPLAY = 50
# NUM_GENE_DISPLAY = 1000










#' Plot the shrinkage visualization for 
#' @param plot_df (dataframe) of shrinkage results 
#'  should have x: an integer 1, 2, ..., indicating position on x axis, each should represent a single grna-gene test
#'  lower_ci and upper_ci: EBCI CIs at 90%
#'  shrunk_value: shrunk value
#'  unshrunk_value and se: original unshrunk results
plot_shrinkViz <- function(plot_df,
                           plot_subtitle=NULL,
                           save_name=NULL, 
                           ylim = c(-1, 1.5), 
                           ybreaks = seq(-20, 20, by = .5),
                           xbreaks = seq(0, 100, by = 5)) {
  delta = .32 # distance between unshrunk and shrunk points and bars on the x axis e.g. .25
  CI_barlengths = .08
  CI_linewidth = .4
  # unshrunk_color = 'deepskyblue2'
  # shrunk_color   = 'deepskyblue4'  
  
  
  # plot_df = plot_df |> # no facets so we don't need to do this?
  #           mutate(   split_type = factor(   split_type, levels = names(split_type_nicename),    labels = split_type_nicename),
  #                  approx_method = factor(approx_method, levels = names(approx_method_nicename), labels = approx_method_nicename))
  
  
  p = ggplot(plot_df |> mutate(CI_unshrunk_lower = unshrunk_value - qnorm(.95)*se, # assuming ALPHA = .1 !! qnorm(1 - (ALPHA/2))
                               CI_unshrunk_upper = unshrunk_value + qnorm(.95)*se),
         aes(x = x)) +
    # geom_point(aes(y = lower_ci) )+
    # geom_point(aes(y = upper_ci)) +
    # geom_point(aes(y = -1.5)) +
    # geom_hline(aes(yintercept = 0)) +
    # ---- Shrunk CIs ---
    geom_segment(aes(y = lower_ci, x = x-CI_barlengths+delta, xend = x+CI_barlengths+delta), linewidth = CI_linewidth, alpha = 1, color = 'deepskyblue4') +
    geom_segment(aes(y = upper_ci, x = x-CI_barlengths+delta, xend = x+CI_barlengths+delta), linewidth = CI_linewidth, alpha = 1, color = 'deepskyblue4') +
    geom_segment(aes(x = x + delta, y = lower_ci, yend = upper_ci),
                 lineend = 'butt', linewidth = 1, alpha = .7, color = 'deepskyblue4') +
    # --- Unshrunk CIs ---
    geom_segment(aes(y = CI_unshrunk_lower, x = x-CI_barlengths, xend = x+CI_barlengths), linewidth = CI_linewidth, alpha = 1, color = 'cyan3') +
    geom_segment(aes(y = CI_unshrunk_upper, x = x-CI_barlengths, xend = x+CI_barlengths), linewidth = CI_linewidth, alpha = 1, color = 'cyan3') +
    geom_segment(aes(x = x, 
                     y    = CI_unshrunk_lower, 
                     yend = CI_unshrunk_upper),
                 lineend = 'butt', linewidth = 1, alpha = .7, color = 'cyan3') + # , arrow = arrow(angle = 90, length = unit(0.03, "npc"))) +
    # geom_segment(aes(y    = eval(parse(text = CIlower_colname)), 
    #                  yend = eval(parse(text = CIupper_colname))),
    #              lineend = 'square', linewidth = 1, alpha = .7, color = 'deepskyblue4') +
    
    # --- Arrow showing shrinkage
    geom_segment(aes(x = x, xend = x+.2,
                     y =    unshrunk_value, 
                     yend = .99 * shrunk_value + .01 * unshrunk_value),
                 lineend = 'square', linejoin = 'bevel', 
                 arrow = arrow(length = unit(0.1,"cm"), angle = 40, type = 'closed'),
                 linewidth = .3, alpha = .7, color = 'cyan3') +
    # --- Unshrunk Point---
    geom_point(aes(x = x        , y = unshrunk_value),           color = 'cyan3', size = 2) +
    # --- Shrunk Point --- 
    geom_point(aes(x = x + delta, y = shrunk_value), shape=18,   color = 'deepskyblue4', size = 2.5) +
    # --- Shrinkage Point---
    geom_point(aes(x = x + delta, y = shrinkage_point), shape=5, color = 'dodgerblue4', size = 1.7, stroke = .7) +
    scale_x_continuous(expand = c(0, 0.07), breaks = xbreaks) +
    scale_y_continuous(expand = c(0, 0), breaks = ybreaks) +
    # --- Labels/
    coord_cartesian(ylim = ylim) +
    labs(x = 'Perturbation-Gene Test', y = 'Estimates',
         title = 'Before and After Robust EBCI',
         subtitle = plot_subtitle) +
    # facet_grid(cols = vars(split_type), 
    #            rows = vars(approx_method)) +
    theme(panel.grid.minor.x = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_blank(),
          legend.background = element_rect(color = 'black'))
  
  
  
  if(!is.null(save_name)) {
    ggsave(plot = p, filename = sprintf('%s/final/shrinkViz/individual/%s.pdf', plot_path, save_name), height = 3, width = 5)
    ggsave(plot = p, filename = sprintf('%s/final/shrinkViz/individual/%s.pdf', plot_path, save_name), height = plot_heights, width = plot_widths)
    ggsave(plot = p, filename = sprintf('%s/final/shrinkViz/individual/%s.png', plot_path, save_name), height = plot_heights, width = plot_widths, dpi = plot_dpi)
  }
  
  
  return(p)
}



# /////////////////////////////////////////////////////////////////////////////////////////////////////
# ================== Shrink Vizualization =============================================================
# /////////////////////////////////////////////////////////////////////////////////////////////////////




for(i in 1:nrow(selectedSettings)) {
  # i = 3
  
  # set up df for plots: matrix approximation + shrinkage
  plot_df = df_sample |> filter(split_type == selectedSettings[i, 'split_type'] & 
                             approx_method == selectedSettings[i, 'approx_method'] & 
                           (is.na(selectedSettings[i, 'rank']) | (rank == selectedSettings[i, 'rank'])))
  
  plot_subtitle = sprintf('%s w/ %s%s', 
                          split_type_nicename[[selectedSettings[i, 'split_type']]],
                          approx_method_nicename[[selectedSettings[i, 'approx_method']]],
                          if(is.na(selectedSettings[i, 'rank'])) {''} else {sprintf(' (rank=%02.f)', selectedSettings[i, 'rank'])})
  plot_subtitle = gsub('\n', ' ', plot_subtitle) # remove \n new lines
  plot_filename = sprintf('%s_matapprox_%s%s', 
                          selectedSettings[i, 'split_type'], 
                          selectedSettings[i, 'approx_method'],
                          if(is.na(selectedSettings[i, 'rank'])) {''} else {sprintf('%02.f', selectedSettings[i, 'rank'])})
  
  
  
  plot_shrinkViz(plot_df = plot_df, 
                 plot_subtitle = plot_subtitle, 
                 save_name = plot_filename, 
                 ylim = c(-2.5, 1.2), 
                 xbreaks = NULL)
                 # xbreaks = seq(from = 0, to = 10, by = 2))
  
  
}















# /////////////////////////////////////////////////////////////////////////////////////////////////////
# ================== One Large Faceted Plot =============================================================
# /////////////////////////////////////////////////////////////////////////////////////////////////////



# === Shrinkage CIs ===

delta = .32 # distance between unshrunk and shrunk points and bars on the x axis e.g. .25
CI_barlengths = .08
CI_linewidth = .4


p_CI = ggplot(df_sample |> 
                mutate(   split_type = factor(   split_type, levels = names(split_type_nicename),    labels = split_type_nicename),
                       approx_method = factor(approx_method, levels = names(approx_method_nicename), labels = approx_method_nicename)) |> 
                mutate(CI_unshrunk_lower = unshrunk_value - qnorm(.95)*se, # assuming ALPHA = .1 !! qnorm(1 - (ALPHA/2))
                       CI_unshrunk_upper = unshrunk_value + qnorm(.95)*se)
              # |> filter(rank == 'rank=15' & approxmethod == 'lowrank') #|> 
              # slice(seq(from = 1, to = n(), by = floor(n()/400))) |>
              # slice(sample_idx) |>
              # mutate(# test = factor(test, levels = c('negative', 'positive', 'discovery')),
              #        x = 1:n()),
              ,
              aes(x = x)) +
  # geom_point(aes(y = lower_ci) )+
  # geom_point(aes(y = upper_ci)) +
  # geom_point(aes(y = -1.5)) +
  # geom_hline(aes(yintercept = 0)) +
  # ---- Shrunk CIs ---
  geom_segment(aes(y = lower_ci, x = x-CI_barlengths+delta, xend = x+CI_barlengths+delta), linewidth = CI_linewidth, alpha = 1, color = 'deepskyblue4') +
  geom_segment(aes(y = upper_ci, x = x-CI_barlengths+delta, xend = x+CI_barlengths+delta), linewidth = CI_linewidth, alpha = 1, color = 'deepskyblue4') +
  geom_segment(aes(x = x + delta, y = lower_ci, yend = upper_ci),
               lineend = 'square', linewidth = 1, alpha = .7, color = 'deepskyblue4') +
  # --- Unshrunk CIs ---
  geom_segment(aes(y = CI_unshrunk_lower, x = x-CI_barlengths, xend = x+CI_barlengths), linewidth = CI_linewidth, alpha = 1, color = 'cyan3') +
  geom_segment(aes(y = CI_unshrunk_upper, x = x-CI_barlengths, xend = x+CI_barlengths), linewidth = CI_linewidth, alpha = 1, color = 'cyan3') +
  geom_segment(aes(x = x, 
                   y    = CI_unshrunk_lower, 
                   yend = CI_unshrunk_upper),
               lineend = 'square', linewidth = 1, alpha = .7, color = 'cyan3') +
  # geom_segment(aes(y    = eval(parse(text = CIlower_colname)), 
  #                  yend = eval(parse(text = CIupper_colname))),
  #              lineend = 'square', linewidth = 1, alpha = .7, color = 'deepskyblue4') +
  
  # --- Arrow showing shrinkage
  geom_segment(aes(x = x, xend = x+.2,
                   y =    unshrunk_value, 
                   yend = .99 * shrunk_value + .01 * unshrunk_value),
               lineend = 'square', linejoin = 'bevel', arrow = arrow(length = unit(0.2,"cm")),
               linewidth = .3, alpha = .7, color = 'cyan3') +
  # geom_segment(aes(y = eval(parse(text = unshrunk_colname)), 
  #                  yend = .99 * eval(parse(text = shrunk_colname)) + 
  #                    .01 * eval(parse(text = unshrunk_colname))),
  #              lineend = 'square', linejoin = 'bevel', arrow = arrow(length = unit(0.2,"cm")),
  #              linewidth = .3, alpha = .7, color = 'deepskyblue2') +
  # --- Unshrunk Point---
  geom_point(aes(x = x      , y = unshrunk_value), color = 'cyan3') +
  # --- Shrunk Point ---
  geom_point(aes(x = x + delta, y = shrunk_value), shape=18, color = 'deepskyblue4') +
  # --- Shrinkage Point---
  geom_point(aes(x = x + delta, y = shrinkage_point), shape=5, color = 'dodgerblue4') +
  # geom_point(aes(y = eval(parse(text = unshrunk_colname))), color = 'deepskyblue2') +
  # geom_point(aes(y = eval(parse(text =   shrunk_colname))), color = 'deepskyblue4') +
  scale_x_continuous(expand = c(0, 0), breaks = seq(0, 1000, by = 1),
                     # limits = c(1, 90)) +
                     # limits = c(65, 175)) +
                     # limits = c(95, 400)) +
                     # limits = c(1, 400)) +
  ) +
  # scale_y_continuous(expand = c(0.025, 0)) +
  coord_cartesian(ylim = c(-1, 1.5)) +
  labs(x = 'Perturbation-Gene Test', y = 'Estimates',
       title = 'Before and After Robust EBCI') +
  # facet_grid(# rows = vars(approxmethod),
  #   cols = vars(rank), 
  #   rows = vars(approx_method)) +
  facet_grid(# rows = vars(approxmethod),
    cols = vars(split_type), 
    rows = vars(approx_method)) +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        # legend.position = 'inside',
        # legend.position.inside = legend_position,
        legend.background = element_rect(color = 'black'),
        strip.background.y = element_rect(fill = NULL), 
        strip.text.x = element_text(           margin = margin(.25,  0, .25,  0, "cm")), 
        strip.text.y = element_text(vjust = 1, margin = margin(.08, .4, .08, .4, "cm")))


p_CI

ggsave(plot = p_CI, filename = sprintf('%s/final/shrinkViz/individual/CI_faceted.pdf', plot_path), height = 11, width = 9)
ggsave(plot = p_CI, filename = sprintf('%s/final/shrinkViz/individual/CI_faceted.png', plot_path), height = 11, width = 9, dpi = 300)










