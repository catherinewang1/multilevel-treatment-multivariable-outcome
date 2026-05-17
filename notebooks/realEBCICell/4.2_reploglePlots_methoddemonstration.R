# Using previously saved objects, make some final matrix plots/figures to demonstrate
# the EBCI method (e.g. matrices at various steps)
# Run in order:
# 1_
# 2_...
# ...
# 4.2_reploglePlots_methoddemonstration.R
# 
# Require previously saved objects at:
# replogleShrinkageParams = readRDS(file = sprintf('%s/replogleShrinkageParams.rds', replogle_save_path))
# 
# gene_index = read.csv(sprintf('%s/gene_index.csv', replogle_save_path))
# grna_index = read.csv(sprintf('%s/grna_index.csv', replogle_save_path))
# 
# df_shrunk   = readRDS(sprintf('%s/EBCI_shrinkage_dataframe.rds', replogle_save_path))
# df_unshrunk = readRDS(sprintf('%s/EBCI_unshrunk_dataframe.rds', replogle_save_path))
# 
# estse_matrices = readRDS(sprintf('%s/estse_matrices.rds'  , replogle_save_path))
# 
# Will save plots at: 
#    <plot_path>/final/method_demonstration/
#                                          /samplesplit
#                                          /nosamplesplit
# Then, I/you should make a figure using these saved plots demonstrating the
# EBCI method. e.g. cells --> SCEPTRE est --> mat approx --> shrink 
# for both samplesplit and nosamplesplit option
#   


# libraries/packages
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
ggplot2::set_theme(theme_bw() + 
                     theme(plot.title = element_text(hjust = .5), 
                           strip.background = element_rect(fill = 'white')))




# /////////////////////////////////////////////////////////////////////////////////////////////////////
# ================== Setup: Script Params ========================================================================
# /////////////////////////////////////////////////////////////////////////////////////////////////////
# here, set the script params

# sceptre_save_path = '../../saves/sceptre/replogle/' # location where sceptre results are located
replogle_save_path= '../../saves/replogle/'         # location to save replogle approximations/shrinkage/etc. results
plot_path         = '../../plots/replogle/EBCI/'    # location to save plots of EBCI analysis on replogle dataset

# limit the upper and lower values
# color_limits = c(-3, 4)
# color_limits = c(-2.5, 2.5)
color_limits = c(-2, 2)


myRed  = colorRampPalette(RColorBrewer::brewer.pal(n = 7, name = "RdBu"))(7)[1]
myBlue = colorRampPalette(RColorBrewer::brewer.pal(n = 7, name = "RdBu"))(7)[7]


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


# color breaks for plotting
color_breaks = sort(union(color_limits, seq(from = round(color_limits[1]), to = round(color_limits[2]))))
color_breaks_label = color_breaks
color_breaks_label[which.min(color_breaks)] = sprintf('<%.0f', min(color_breaks))
color_breaks_label[which.max(color_breaks)] = sprintf('>%.0f', max(color_breaks))



# get the original values from train dataset- (does not appear in unshrunk or shrunk dfs...)

# load in as matrix form:
estse_matrices = readRDS(sprintf('%s/estse_matrices.rds'  , replogle_save_path))
df_train = reshape2::melt(estse_matrices$est_matrices$train , varnames = c('grna', 'gene'), value.name = 'unshrunk_value')
df_train = merge(merge(df_train, gene_index, by = 'gene'), grna_index, by = 'grna') |>  # add in the previously saved orderings
           dplyr::mutate(gene = as.factor(gene), grna = as.factor(grna))                # make gene and grna into factors
rm(estse_matrices); gc()







# /////////////////////////////////////////////////////////////////////////////////////////////////////
# ================== MAT Demonstrating Method ===========================================
# /////////////////////////////////////////////////////////////////////////////////////////////////////
dir.create(sprintf('%s/final/', plot_path))
dir.create(sprintf('%s/final/method_demonstration/', plot_path))
dir.create(sprintf('%s/final/method_demonstration/samplesplit/', plot_path))
dir.create(sprintf('%s/final/method_demonstration/nosamplesplit/', plot_path))

demonstration_heights = 6
demonstration_widths  = 3
demonstration_dpi     = 300

NUM_GRNA_DISPLAY = 50
NUM_GENE_DISPLAY = 1000


#' @param plot_df (dataframe) plotting dataframe where grna_idx and gene_idx indicate the y and x positions respectively
#' @param plot_title (character) title of the plot 
#' @param value_name (character) name of the column for fill colors in the plot (e.g. unshrunk  or shrinkage point etc...)
final_plot_matrix <- function(plot_df, 
                              plot_title, 
                              fill_colname, fill_plotname,
                              save_name=NULL) {
  
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
         y = "grna", x = "gene", 
         fill = fill_plotname) +
    theme(axis.ticks = element_blank(), 
          axis.text = element_blank(), 
          legend.position = 'bottom',
          legend.key.height = unit(.3, 'cm'),
          legend.key.width  = unit(1, 'cm'),
          legend.text = element_text(size = 7), 
          legend.title = element_text(size = 10, vjust = 1), 
          legend.justification.bottom = 'center') # 'left')
  
  if(!is.null(save_name)) {
    ggsave(plot = p, filename = sprintf('%s/final/method_demonstration/%s.pdf', plot_path, save_name), height = demonstration_heights, width = demonstration_widths)
    ggsave(plot = p, filename = sprintf('%s/final/method_demonstration/%s.png', plot_path, save_name), height = demonstration_heights, width = demonstration_widths, dpi = demonstration_dpi)
  }
  
  
  return(p)
}


# --------- > sample split method ------------------------------------------------
# ///////// sample split method ////////////////////////////////////////////////


# original train estimates (for sample split method)- Different from the others!
final_plot_matrix(plot_df = df_train,
                  plot_title = "Unshrunk Estimates on Train Split", 
                  fill_colname = 'unshrunk_value', fill_plotname = "value", 
                  save_name = 'samplesplit/unshrunk_train')

# original test estimates (for sample split method)
final_plot_matrix(plot_df = df_unshrunk |> filter(grna_idx <= NUM_GRNA_DISPLAY & gene_idx <= NUM_GENE_DISPLAY) |> 
                                           filter(split_type == 'samplesplit'),
                  plot_title = "Unshrunk Estimates on Test Split", 
                  fill_colname = 'unshrunk_value', fill_plotname = "value", 
                  save_name = 'samplesplit/unshrunk_test')

# mat approx for SVD 3
final_plot_matrix(plot_df = df |> filter(grna_idx <= NUM_GRNA_DISPLAY & gene_idx <= NUM_GENE_DISPLAY) |>
                                  filter(split_type == 'samplesplit' & approx_method == 'SVD' & rank == 3),
                  plot_title = "SVD on Train Split", 
                  fill_colname = 'shrinkage_point', fill_plotname = "value", 
                  save_name = 'samplesplit/matapprox_SVD03')


# # mat approx for SVD 5
# final_plot_matrix(plot_df = df |> filter(grna_idx <= NUM_GRNA_DISPLAY & gene_idx <= NUM_GENE_DISPLAY) |>
#                                   filter(split_type == 'samplesplit' & approx_method == 'SVD' & rank == 5),
#                   plot_title = "SVD on Train Split", 
#                   fill_colname = 'shrinkage_point', fill_plotname = "value", 
#                   save_name = 'samplesplit/matapprox_SVD05')

# mat approx for sparseSVD 20
final_plot_matrix(plot_df = df |> filter(grna_idx <= NUM_GRNA_DISPLAY & gene_idx <= NUM_GENE_DISPLAY) |>
                    filter(split_type == 'samplesplit' & approx_method == 'sparseSVD' & rank == 20),
                  plot_title = "Sparse SVD on Train Split", 
                  fill_colname = 'shrinkage_point', fill_plotname = "value", 
                  save_name = 'samplesplit/matapprox_SpSVD20')

# # mat approx for sparseSVD 30
# final_plot_matrix(plot_df = df |> filter(grna_idx <= NUM_GRNA_DISPLAY & gene_idx <= NUM_GENE_DISPLAY) |>
#                     filter(split_type == 'samplesplit' & approx_method == 'sparseSVD' & rank == 30),
#                   plot_title = "Sparse SVD on Train Split", 
#                   fill_colname = 'shrinkage_point', fill_plotname = "value", 
#                   save_name = 'samplesplit/matapprox_SpSVD30')

# # mat approx for sparse softImpute (doesn't seem to change based on the rank parameters)
# final_plot_matrix(plot_df = df |> filter(grna_idx <= NUM_GRNA_DISPLAY & gene_idx <= NUM_GENE_DISPLAY) |>
#                     filter(split_type == 'samplesplit' & approx_method == 'softImpute' & rank == 5),
#                   plot_title = "softImpute on Train Split", 
#                   fill_colname = 'shrinkage_point', fill_plotname = "value", 
#                   save_name = 'samplesplit/matapprox_softImpute05')

# # mat approx for sparse spectral Biclust 
# final_plot_matrix(plot_df = df |> filter(grna_idx <= NUM_GRNA_DISPLAY & gene_idx <= NUM_GENE_DISPLAY) |>
#                     filter(split_type == 'samplesplit' & approx_method == 'spectralbiclust' & rank == 5),
#                   plot_title = "Spectral Biclustering on Train Split", 
#                   fill_colname = 'shrinkage_point', fill_plotname = "value", 
#                   save_name = 'samplesplit/matapprox_spectralBiclust05')


# mat approx for sparse spectral Biclust w thresholding
final_plot_matrix(plot_df = df |> filter(grna_idx <= NUM_GRNA_DISPLAY & gene_idx <= NUM_GENE_DISPLAY) |>
                    filter(split_type == 'samplesplit' & approx_method == 'spectralbiclust_threshold' & rank == 5),
                  plot_title = "Spectral Biclustering w/ Thresholding on Train Split", 
                  fill_colname = 'shrinkage_point', fill_plotname = "value", 
                  save_name = 'samplesplit/matapprox_spectralBiclustThresh05')

# mat approx for zeros
final_plot_matrix(plot_df = df |> filter(grna_idx <= NUM_GRNA_DISPLAY & gene_idx <= NUM_GENE_DISPLAY) |>
                    filter(split_type == 'samplesplit' & approx_method == 'zeros'),
                  plot_title = "Zeros on Train Split", 
                  fill_colname = 'shrinkage_point', fill_plotname = "value", 
                  save_name = 'samplesplit/matapprox_zeros')



# shrinkage estimates for SVD 3
final_plot_matrix(plot_df = df |> filter(grna_idx <= NUM_GRNA_DISPLAY & gene_idx <= NUM_GENE_DISPLAY) |>
                    filter(split_type == 'samplesplit' & approx_method == 'SVD' & rank == 3),
                  plot_title = "Shrunk Estimates using SVD", 
                  fill_colname = 'shrunk_value', fill_plotname = "value", 
                  save_name = 'samplesplit/shrunk_SVD03')

# # shrinkage estimates for SVD 5
# final_plot_matrix(plot_df = df |> filter(grna_idx <= NUM_GRNA_DISPLAY & gene_idx <= NUM_GENE_DISPLAY) |>
#                     filter(split_type == 'samplesplit' & approx_method == 'SVD' & rank == 5),
#                   plot_title = "Shrunk Estimates using SVD", 
#                   fill_colname = 'shrunk_value', fill_plotname = "value", 
#                   save_name = 'samplesplit/shrunk_SVD05')

# # shrinkage estimates for sparse SVD 3
# final_plot_matrix(plot_df = df |> filter(grna_idx <= NUM_GRNA_DISPLAY & gene_idx <= NUM_GENE_DISPLAY) |>
#                     filter(split_type == 'samplesplit' & approx_method == 'sparseSVD' & rank == 3),
#                   plot_title = "Shrunk Estimates using Sparse SVD", 
#                   fill_colname = 'shrunk_value', fill_plotname = "value", 
#                   save_name = 'samplesplit/shrunk_SpSVD05')

# shrinkage estimates for sparse SVD 20
final_plot_matrix(plot_df = df |> filter(grna_idx <= NUM_GRNA_DISPLAY & gene_idx <= NUM_GENE_DISPLAY) |>
                    filter(split_type == 'samplesplit' & approx_method == 'sparseSVD' & rank == 20),
                  plot_title = "Shrunk Estimates using Sparse SVD", 
                  fill_colname = 'shrunk_value', fill_plotname = "value", 
                  save_name = 'samplesplit/shrunk_SpSVD20')

# # shrinkage estimates for sparse spectral Biclust 
# final_plot_matrix(plot_df = df |> filter(grna_idx <= NUM_GRNA_DISPLAY & gene_idx <= NUM_GENE_DISPLAY) |>
#                     filter(split_type == 'samplesplit' & approx_method == 'spectralbiclust' & rank == 5),
#                   plot_title = "Shrunk Estimates using Spectral Biclustering", 
#                   fill_colname = 'shrunk_value', fill_plotname = "value", 
#                   save_name = 'samplesplit/shrunk_spectralBiclust05')


# shrinkage estimates for sparse spectral Biclust w thresholding
final_plot_matrix(plot_df = df |> filter(grna_idx <= NUM_GRNA_DISPLAY & gene_idx <= NUM_GENE_DISPLAY) |>
                    filter(split_type == 'samplesplit' & approx_method == 'spectralbiclust_threshold' & rank == 5),
                  plot_title = "Shrunk Estimates using Spectral Biclustering w/ Thresholding", 
                  fill_colname = 'shrunk_value', fill_plotname = "value", 
                  save_name = 'samplesplit/shrunk_spectralBiclustThresh05')



# shrinkage estimates for zeros
final_plot_matrix(plot_df = df |> filter(grna_idx <= NUM_GRNA_DISPLAY & gene_idx <= NUM_GENE_DISPLAY) |>
                    filter(split_type == 'samplesplit' & approx_method == 'zeros'),
                  plot_title = "Shrunk Estimates using Zeros", 
                  fill_colname = 'shrunk_value', fill_plotname = "value", 
                  save_name = 'samplesplit/shrunk_zeros')



# --------- > no sample split method ------------------------------------------------
# ///////// no sample split method ////////////////////////////////////////////////


# original all estimates (for no sample split method)
final_plot_matrix(plot_df = df_unshrunk |> filter(grna_idx <= NUM_GRNA_DISPLAY & gene_idx <= NUM_GENE_DISPLAY) |>
                                           filter(split_type == 'nosamplesplit'),
                  plot_title = "Unshrunk Estimates on Full Data", 
                  fill_colname = 'unshrunk_value', fill_plotname = "value", 
                  save_name = 'nosamplesplit/unshrunk_all')

# mat approx for SVD 3
final_plot_matrix(plot_df = df |> filter(grna_idx <= NUM_GRNA_DISPLAY & gene_idx <= NUM_GENE_DISPLAY) |>
                    filter(split_type == 'nosamplesplit' & approx_method == 'SVD' & rank == 3),
                  plot_title = "SVD on Full Data",
                  fill_colname = 'shrinkage_point', fill_plotname = "value",
                  save_name = 'nosamplesplit/matapprox_SVD03')

# # mat approx for SVD 5
# final_plot_matrix(plot_df = df |> filter(grna_idx <= NUM_GRNA_DISPLAY & gene_idx <= NUM_GENE_DISPLAY) |>
#                     filter(split_type == 'nosamplesplit' & approx_method == 'SVD' & rank == 5),
#                   plot_title = "SVD on Full Data", 
#                   fill_colname = 'shrinkage_point', fill_plotname = "value", 
#                   save_name = 'nosamplesplit/matapprox_SVD05')


# mat approx for sparse SVD 
final_plot_matrix(plot_df = df |> filter(grna_idx <= NUM_GRNA_DISPLAY & gene_idx <= NUM_GENE_DISPLAY) |>
                    filter(split_type == 'nosamplesplit' & approx_method == 'sparseSVD' & rank == 20),
                  plot_title = "Sparse SVD on Full Data",
                  fill_colname = 'shrinkage_point', fill_plotname = "value",
                  save_name = 'nosamplesplit/matapprox_SpSVD20')


# mat approx for spectral Biclust 
final_plot_matrix(plot_df = df |> filter(grna_idx <= NUM_GRNA_DISPLAY & gene_idx <= NUM_GENE_DISPLAY) |>
                    filter(split_type == 'nosamplesplit' & approx_method == 'spectralbiclust' & rank == 5),
                  plot_title = "Spectral Biclustering on Full Data", 
                  fill_colname = 'shrinkage_point', fill_plotname = "value", 
                  save_name = 'nosamplesplit/matapprox_spectralBiclust05')

# # mat approx for sparse spectral Biclust w/ thresholding
# final_plot_matrix(plot_df = df |> filter(grna_idx <= NUM_GRNA_DISPLAY & gene_idx <= NUM_GENE_DISPLAY) |>
#                     filter(split_type == 'nosamplesplit' & approx_method == 'spectralbiclust_threshold' & rank == 5),
#                   plot_title = "Spectral Biclustering on Full Data", 
#                   fill_colname = 'shrinkage_point', fill_plotname = "value", 
#                   save_name = 'nosamplesplit/matapprox_spectralBiclustThresh05')

# mat approx for sparse spectral Biclust 
final_plot_matrix(plot_df = df |> filter(grna_idx <= NUM_GRNA_DISPLAY & gene_idx <= NUM_GENE_DISPLAY) |>
                    filter(split_type == 'nosamplesplit' & approx_method == 'zeros'),
                  plot_title = "Zeros", 
                  fill_colname = 'shrinkage_point', fill_plotname = "value", 
                  save_name = 'nosamplesplit/matapprox_zeros')

# shrinkage estimates for SVD 3
final_plot_matrix(plot_df = df |> filter(grna_idx <= NUM_GRNA_DISPLAY & gene_idx <= NUM_GENE_DISPLAY) |>
                    filter(split_type == 'nosamplesplit' & approx_method == 'SVD' & rank == 3),
                  plot_title = "Shrunk Estimates using SVD", 
                  fill_colname = 'shrunk_value', fill_plotname = "value", 
                  save_name = 'nosamplesplit/shrunk_SVD03')

# shrinkage estimates for sparse SVD 
final_plot_matrix(plot_df = df |> filter(grna_idx <= NUM_GRNA_DISPLAY & gene_idx <= NUM_GENE_DISPLAY) |>
                    filter(split_type == 'nosamplesplit' & approx_method == 'sparseSVD' & rank == 20),
                  plot_title = "Shrunk Estimates using Sparse SVD", 
                  fill_colname = 'shrunk_value', fill_plotname = "value", 
                  save_name = 'nosamplesplit/shrunk_SpSVD20')

# # shrinkage estimates for SVD 5
# final_plot_matrix(plot_df = df |> filter(grna_idx <= NUM_GRNA_DISPLAY & gene_idx <= NUM_GENE_DISPLAY) |>
#                     filter(split_type == 'nosamplesplit' & approx_method == 'SVD' & rank == 5),
#                   plot_title = "Shrunk Estimates using SVD", 
#                   fill_colname = 'shrunk_value', fill_plotname = "value", 
#                   save_name = 'nosamplesplit/shrunk_SVD05')

# shrinkage estimates for sparse spectral Biclust 
final_plot_matrix(plot_df = df |> filter(grna_idx <= NUM_GRNA_DISPLAY & gene_idx <= NUM_GENE_DISPLAY) |>
                    filter(split_type == 'nosamplesplit' & approx_method == 'spectralbiclust' & rank == 5),
                  plot_title = "Shrunk Estimates using Spectral Biclustering", 
                  fill_colname = 'shrunk_value', fill_plotname = "value", 
                  save_name = 'nosamplesplit/shrunk_spectralBiclust05')

# shrinkage estimates for zeros
final_plot_matrix(plot_df = df |> filter(grna_idx <= NUM_GRNA_DISPLAY & gene_idx <= NUM_GENE_DISPLAY) |>
                    filter(split_type == 'nosamplesplit' & approx_method == 'zeros'),
                  plot_title = "Shrunk Estimates using Spectral Biclustering", 
                  fill_colname = 'shrunk_value', fill_plotname = "value", 
                  save_name = 'nosamplesplit/shrunk_zeros')








# /////////////////////////////////////////////////////////////////////////////////////////////////////
# ================== Trash / Old  ===========================================
# /////////////////////////////////////////////////////////////////////////////////////////////////////


if(F) {

  
  
  # # checking
  # names(estse_matrices)
  # names(estse_matrices$est_matrices)
  # estse_matrices$est_matrices$train  |> dim()
  # colnames(estse_matrices$est_matrices$train )
  # row.names(estse_matrices$est_matrices$train )
  # head(df_train); dim(df_train)
  # df_train$grna |> unique() |> length()  # should be 50   different grnas
  # df_train$gene |> unique() |> length()  # should be 1000 different genes
  
  # # wrong loaded in object
  # df_sceptre_train = shrink_results[['samplesplit']][['zeros']][['ebci_res']] |>
  #   dplyr::select(grna, gene, unshrunk_value, se, weight) |> 
  #   dplyr::mutate(split_type == 'samplesplit')
  # 
  # 
  # shrink_results = readRDS(sprintf('%s/EBCI_shrinkage_results.rds'  , replogle_save_path))
  # df_sceptre_train = shrink_results[['samplesplit']][['zeros']][['ebci_res']] |>
  #   dplyr::select(grna, gene, unshrunk_value, se, weight) |> 
  #   dplyr::mutate(split_type == 'samplesplit')
  # approxmatrices = readRDS(sprintf('%s/approxmatrices.rds'  , replogle_save_path))
  # names(approxmatrices)
  # names(approxmatrices$all)
  # names(shrink_results)
  # 
  # rm(shrink_results);
  
# original test estimates (for sample split method)
plot_df = df_unshrunk |> filter(grna_idx <= NUM_GRNA_DISPLAY & gene_idx <= NUM_GENE_DISPLAY)
ggplot(plot_df) +
  geom_raster(aes(x = gene_idx, y = grna_idx, fill = unshrunk_value)) +
  coord_cartesian(expand = c(0, 0, 0, 0)) + 
  scale_fill_gradient2(limits = color_limits, # set color limits
                       oob=scales::squish, # if outside lims, set to limits
                       midpoint = 0, 
                       high = myRed, low = myBlue, mid = 'white',
                       breaks = color_breaks,
                       labels = color_breaks_label) +
  labs(title = "Estimates on Test Split", 
       y = "grna", x = "gene", 
       fill = "estimate") +
  theme(axis.ticks = element_blank(), 
        axis.text = element_blank(), 
        legend.position = 'bottom',
        legend.key.height = unit(.3, 'cm'),
        legend.key.width  = unit(1, 'cm'),
        legend.text = element_text(size = 7))


ggsave(sprintf('%s/final/method_demonstration/unshrunk_test.pdf', plot_path), height = demonstration_heights, width = demonstration_widths)
ggsave(sprintf('%s/final/method_demonstration/unshrunk_test.png', plot_path), height = demonstration_heights, width = demonstration_widths, dpi = 300)


# mat approx for SVD 5
plot_df = df |> 
  filter(grna_idx <= NUM_GRNA_DISPLAY & gene_idx <= NUM_GENE_DISPLAY) |>
  filter(split_type == 'samplesplit' & approx_method == 'SVD' & rank == 3)
ggplot(plot_df) +
  geom_raster(aes(x = gene_idx, y = grna_idx, fill = shrinkage_point)) +
  coord_cartesian(expand = c(0, 0, 0, 0)) + 
  scale_fill_gradient2(limits = color_limits, # set color limits
                       oob=scales::squish, # if outside lims, set to limits
                       midpoint = 0, 
                       high = myRed, low = myBlue, mid = 'white',
                       breaks = color_breaks,
                       labels = color_breaks_label) +
  labs(title = "SVD (rank=3) on Train Split", 
       y = "grna", x = "gene", 
       fill = "value") +
  theme(axis.ticks = element_blank(), 
        axis.text = element_blank(), 
        legend.position = 'bottom',
        legend.key.height = unit(.3, 'cm'),
        legend.key.width  = unit(.7, 'cm'),
        legend.text = element_text(size = 7), 
        legend.title = element_text(size = 10, vjust = 1), legend.justification.bottom = 'left')


ggsave(sprintf('%s/final/method_demonstration/matapprox_svd03.pdf', plot_path), height = demonstration_heights, width = demonstration_widths)
ggsave(sprintf('%s/final/method_demonstration/matapprox_svd03.png', plot_path), height = demonstration_heights, width = demonstration_widths, dpi = demonstration_dpi)



}

