# Using previously saved objects, make some final matrix plots to showcase the
# final results of the EBCI method (e.g. matrices at various steps)
# Run in order:
# 1_...
# 2_...
# 3_...
# ...
# 4.2_reploglePlots_methoddemonstration.R
# 4.3_reploglePlots_showcaseMatrices.R
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



# Select some settings to showcase:
selectedSettings = data.frame(matrix(
  c(# settings for sample splitting
    'samplesplit',                         'SVD',  1,
    'samplesplit',                         'SVD',  3,
    'samplesplit',                         'SVD',  5,
    'samplesplit',                         'SVD',  10,
    'samplesplit',                         'SVD',  20,
    'samplesplit',                         'SVD',  30,
    'samplesplit',                   'sparseSVD', 20,
    'samplesplit',        'sparseSVD_autoparams',  1,
    'samplesplit',        'sparseSVD_autoparams',  3,
    'samplesplit',        'sparseSVD_autoparams',  5,
    'samplesplit',        'sparseSVD_autoparams',  10,
    'samplesplit',        'sparseSVD_autoparams',  20,
    'samplesplit',        'sparseSVD_autoparams',  30,
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
# rm(estse_matrices); gc()



# /////////////////////////////////////////////////////////////////////////////////////////////////////
# ================== Setup Directory & Fn =============================================================
# /////////////////////////////////////////////////////////////////////////////////////////////////////
dir.create(sprintf('%s/final/', plot_path))
dir.create(sprintf('%s/final/shrinkPointPlots/', plot_path))    # 
dir.create(sprintf('%s/final/shrinkPointPlots/individual/', plot_path)) # create subfolder for these matrices, can create keynote fig outside this




NUM_TESTS_DISPLAY = 1000

plot_heights = 4
plot_widths  = 4
plot_dpi     = 300

final_plot_shrink_points <- function(plot_df, 
                                     plot_subtitle = NULL,
                                     save_name=NULL, 
                                     extra_theme_settings = theme()) {
  
  
  p = ggplot(plot_df , 
             aes(x = unshrunk_value - shrinkage_point, y = shrunk_value - shrinkage_point)) +
    geom_abline(aes(slope = 0, intercept = 0), color = 'gray') +
    geom_abline(aes(slope = 1, intercept = 0)) +
    geom_point(color= 'turquoise4', alpha = .4, size = 1.5, stroke = NA) +
    coord_cartesian(xlim = c(-.5, .5), 
                    ylim = c(-.5, .5)) +
    labs(title = 'Centered Shrunk vs Unshrunk Estimates', 
         subtitle = plot_subtitle,
         x = 'unshrunk estimate - shrinkage point', 
         y = 'shrunk estimate - shrinkage point') +
    theme(panel.grid = element_blank()) + 
    extra_theme_settings
  
  if(!is.null(save_name)) {
    ggsave(plot = p, filename = sprintf('%s/final/shrinkPointPlots/individual/%s.pdf', plot_path, save_name), height = plot_heights, width = plot_widths)
    ggsave(plot = p, filename = sprintf('%s/final/shrinkPointPlots/individual/%s.png', plot_path, save_name), height = plot_heights, width = plot_widths, dpi = plot_dpi)
  }
  
  
}








# /////////////////////////////////////////////////////////////////////////////////////////////////////
# ================== Points: Shrunk vs Original =====================================================
# /////////////////////////////////////////////////////////////////////////////////////////////////////

for(i in 1:nrow(selectedSettings)) {
  # i = 3
  
  # set up df for plots: matrix approximation + shrinkage
  plot_df = df |> filter(split_type == selectedSettings[i, 'split_type'] & 
                           approx_method == selectedSettings[i, 'approx_method'] & 
                           (is.na(selectedSettings[i, 'rank']) | (rank == selectedSettings[i, 'rank'])))
  
  plot_subtitle = sprintf('%s%s', 
                          selectedSettings[i, 'approx_method'],
                          if(is.na(selectedSettings[i, 'rank'])) {''} else {sprintf(' (rank=%02.f)', selectedSettings[i, 'rank'])})
  
  
  # matrix approximations
  plot_filename = sprintf('%s_%s%s', 
                          selectedSettings[i, 'split_type'], 
                          selectedSettings[i, 'approx_method'],
                          if(is.na(selectedSettings[i, 'rank'])) {''} else {sprintf('%02.f', selectedSettings[i, 'rank'])})
  
  
  final_plot_shrink_points(plot_df = plot_df  |> sample_n(size = NUM_TESTS_DISPLAY),
                          plot_subtitle = plot_subtitle, 
                          save_name = plot_filename)
  
}











# /////////////////////////////////////////////////////////////////////////////////////////////////////
# ================== Trash =============================================================================
# /////////////////////////////////////////////////////////////////////////////////////////////////////

# # uncentered
# ggplot(plot_df |> sample_n(size = 400), 
#        aes(x = unshrunk_value , y = shrunk_value )) +
#   geom_abline(aes(slope = 0, intercept = 0), color = 'gray') +
#   geom_abline(aes(slope = 1, intercept = 0)) +
#   geom_point() +
#   coord_cartesian(xlim = c(-1, 1), ylim = c(-.5, .5))


