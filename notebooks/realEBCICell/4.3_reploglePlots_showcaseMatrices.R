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
dir.create(sprintf('%s/final/showcaseMatrices/', plot_path))    # 
dir.create(sprintf('%s/final/showcaseMatrices/mat', plot_path)) # create subfolder for these matrices, can create keynote fig outside this

matrixplot_heights = 6
matrixplot_widths  = 3
matrixplot_dpi     = 300

NUM_GRNA_DISPLAY = 50
NUM_GENE_DISPLAY = 1000


#' @param plot_df (dataframe) plotting dataframe where grna_idx and gene_idx indicate the y and x positions respectively
#' @param plot_title (character) title of the plot 
#' @param value_name (character) name of the column for fill colors in the plot (e.g. unshrunk  or shrinkage point etc...)
#' @param color_limit_inequality (vector) of 2 booleans, whether or not to include the < and > signs in the fill legend
#'         (ie set to c(FALSE, TRUE) for standard error plots)
final_plot_matrix <- function(plot_df, 
                              plot_title, 
                              fill_colname, fill_plotname,
                              plot_subtitle = NULL,
                              save_name=NULL, 
                              color_limits=MY_COLOR_LIMITS,
                              color_limit_inequality = c(TRUE, TRUE),
                              extra_theme_settings = theme()) {
  # color breaks for plotting
  color_breaks = sort(union(color_limits, seq(from = round(color_limits[1]), to = round(color_limits[2]))))
  color_breaks_label = color_breaks
  if(color_limit_inequality[1]) {
    color_breaks_label[which.min(color_breaks)] = sprintf('<%.0f', min(color_breaks))
  } else {
    color_breaks_label[which.min(color_breaks)] = sprintf('%.0f', min(color_breaks))
  }
  if(color_limit_inequality[2]) {
    color_breaks_label[which.max(color_breaks)] = sprintf('>%.0f', max(color_breaks))
  } else {
    color_breaks_label[which.max(color_breaks)] = sprintf('%.0f', max(color_breaks))
  }
  # color_breaks_label[which.min(color_breaks)] = sprintf('<%.0f', min(color_breaks))
  # color_breaks_label[which.max(color_breaks)] = sprintf('>%.0f', max(color_breaks))
  
  # set up df and make plot
  plot_df$temp_fill = plot_df[, fill_colname] # just add a new col here, so that it's easier to call during plotting
  p = ggplot(plot_df) +
    geom_raster(aes(x = gene_idx, y = grna_idx, fill = temp_fill)) +
    coord_cartesian(expand = c(0, 0, 0, 0)) + 
    scale_fill_gradient2(limits = color_limits, # set color limits
                         oob=scales::squish, # if outside lims, set to limits
                         midpoint = 0, 
                         high = myRed, low = myBlue, mid = 'white',
                         breaks = color_breaks,
                         labels = color_breaks_label) +
    labs(title = plot_title, 
         subtitle = plot_subtitle,
         y = "grna", x = "gene", 
         fill = fill_plotname) +
    theme(axis.ticks = element_blank(), 
          axis.text = element_blank(), 
          legend.position = 'bottom',
          legend.key.height = unit(.3, 'cm'),
          legend.key.width  = unit(1, 'cm'),
          legend.text = element_text(size = 7), 
          legend.title = element_text(size = 10, vjust = 1), 
          legend.justification.bottom = 'center') + # 'left')
    extra_theme_settings
  
  if(!is.null(save_name)) {
    ggsave(plot = p, filename = sprintf('%s/final/showcaseMatrices/mat/%s.pdf', plot_path, save_name), height = matrixplot_heights, width = matrixplot_widths)
    ggsave(plot = p, filename = sprintf('%s/final/showcaseMatrices/mat/%s.png', plot_path, save_name), height = matrixplot_heights, width = matrixplot_widths, dpi = matrixplot_dpi)
  }
  
  
  return(p)
}



# /////////////////////////////////////////////////////////////////////////////////////////////////////
# ================== MAT ORIGINAL SCEPTRE ESTIMATES and SE ============================================
# /////////////////////////////////////////////////////////////////////////////////////////////////////

# might have to customize these

# original train SCEPTRE (for sample split method)- Different from the others!   /////////////////////
# estimates
final_plot_matrix(plot_df = df_train,
                  plot_title = "Unshrunk Estimate", 
                  plot_subtitle = "Train Split", 
                  fill_colname = 'unshrunk_value', fill_plotname = "SCEPTRE\nEstimate", 
                  save_name = 'samplesplit_unshrunk_train', 
                  # color_limits = c(-2, 2),
                  extra_theme_settings = theme(legend.title = element_text(size = 9, vjust = 1.7, hjust = .5),
                                                 plot.title = element_text(size = 12, vjust = .5),
                                               legend.key.width  = unit(.95, 'cm')))
# se
final_plot_matrix(plot_df = df_train,
                  plot_title = "SE of Unshrunk Estimate", 
                  plot_subtitle = "Train Split", 
                  fill_colname = 'se', fill_plotname = "standard error\n(se)", 
                  save_name = 'samplesplit_unshrunk_train_se',
                  color_limits = c(0, 2),
                  color_limit_inequality = c(FALSE, TRUE),
                  extra_theme_settings = theme(legend.title = element_text(size = 9, vjust = 1.7, hjust = .5),
                                                 plot.title = element_text(size = 12, vjust = .5),
                                               legend.key.width  = unit(.85, 'cm')))


# original test SCEPTRE (for sample split method) ///////////////////////////////////////////////////
# estimates
final_plot_matrix(plot_df = df_unshrunk |> filter(grna_idx <= NUM_GRNA_DISPLAY & gene_idx <= NUM_GENE_DISPLAY) |> 
                    filter(split_type == 'samplesplit'),
                  plot_title = "Unshrunk Estimate", 
                  plot_subtitle = "Test Split", 
                  fill_colname = 'unshrunk_value', fill_plotname = "SCEPTRE\nEstimate", 
                  save_name = 'samplesplit_unshrunk_test',
                  # color_limits = c(-2, 2),
                  extra_theme_settings = theme(legend.title = element_text(size = 9, vjust = 1.7, hjust = .5),
                                                 plot.title = element_text(size = 12, vjust = .5),
                                               legend.key.width  = unit(.95, 'cm')))
# se
final_plot_matrix(plot_df = df_unshrunk |> filter(grna_idx <= NUM_GRNA_DISPLAY & gene_idx <= NUM_GENE_DISPLAY) |> 
                    filter(split_type == 'samplesplit'),
                  plot_title = "SE of Unshrunk Estimate", 
                  plot_subtitle = "Test Split",
                  fill_colname = 'se', fill_plotname = "standard error\n(se)", 
                  save_name = 'samplesplit_unshrunk_test_se',
                  color_limits = c(0, 2),
                  color_limit_inequality = c(FALSE, TRUE),
                  extra_theme_settings = theme(legend.title = element_text(size = 9, vjust = 1.7, hjust = .5),
                                                 plot.title = element_text(size = 12, vjust = .5),
                                               legend.key.width  = unit(.85, 'cm')))


# original all SCEPTRE (for no sample split method) ////////////////////////////////////////////////
# estimates
final_plot_matrix(plot_df = df_unshrunk |> filter(grna_idx <= NUM_GRNA_DISPLAY & gene_idx <= NUM_GENE_DISPLAY) |>
                    filter(split_type == 'nosamplesplit'),
                  plot_title = "Unshrunk Estimate", 
                  plot_subtitle = "Full Dataset",
                  fill_colname = 'unshrunk_value', fill_plotname = "SCEPTRE\nEstimate", 
                  save_name = 'nosamplesplit_unshrunk_all',
                  # color_limits = c(-2, 2),
                  extra_theme_settings = theme(legend.title = element_text(size = 9, vjust = 1.7, hjust = .5),
                                                 plot.title = element_text(size = 12, vjust = .5),
                                               legend.key.width  = unit(.95, 'cm')))
# se
final_plot_matrix(plot_df = df_unshrunk |> filter(grna_idx <= NUM_GRNA_DISPLAY & gene_idx <= NUM_GENE_DISPLAY) |>
                    filter(split_type == 'nosamplesplit'),
                  plot_title = "SE of Unshrunk Estimate", 
                  plot_subtitle = "Full Dataset",
                  fill_colname = 'unshrunk_value', fill_plotname = "standard error\n(se)", 
                  save_name = 'nosamplesplit_unshrunk_all_se',
                  color_limits = c(0, 2),
                  color_limit_inequality = c(FALSE, TRUE),
                  extra_theme_settings = theme(legend.title = element_text(size = 9, vjust = 1.7, hjust = .5),
                                                 plot.title = element_text(size = 12, vjust = .5),
                                          legend.key.width  = unit(.85, 'cm')))





# /////////////////////////////////////////////////////////////////////////////////////////////////////
# ================== MAT chosen showcase settings =====================================================
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
  plot_filename = sprintf('%s_matapprox_%s%s', 
                          selectedSettings[i, 'split_type'], 
                          selectedSettings[i, 'approx_method'],
                          if(is.na(selectedSettings[i, 'rank'])) {''} else {sprintf('%02.f', selectedSettings[i, 'rank'])})
  
  final_plot_matrix(plot_df = plot_df,
                    plot_title = "Shrinkage Point", 
                    plot_subtitle = plot_subtitle, 
                    fill_colname = 'shrinkage_point', fill_plotname = "value", 
                    save_name = plot_filename, 
                    )
  
  # shrunk estimates
  plot_filename = sprintf('%s_shrunk_%s%s', 
                          selectedSettings[i, 'split_type'], 
                          selectedSettings[i, 'approx_method'],
                          if(is.na(selectedSettings[i, 'rank'])) {''} else {sprintf('%02.f', selectedSettings[i, 'rank'])})
  
  final_plot_matrix(plot_df = plot_df,
                    plot_title = "Shrunk Estimate", 
                    plot_subtitle = plot_subtitle, 
                    fill_colname = 'shrunk_value', fill_plotname = "value", 
                    save_name = plot_filename, 
  )



}






