# Use saved replogleShrinkage results to make plots
# 
# Run in order:
# 1_...
# 2_replogleShrinkage.r
# 3_reploglePvals.r
# 4_reploglePlots.r
# 
# Requires previously saved shrinkage results
# Specifically: <plot_path>/sceptre_obj_[train|test|all].rds 
#               <replogle_save_path>/EBCI_shrinkage_dataframe.rds
# 
# 
#
#



# create some final plots to showcase shrinkage results
# while performing shrinkage (scripts 1-3) plots were made, but not finalized
# create specific better looking plots here
# all that is needed should be in the saved dataframe
# sprintf('%s/EBCI_shrinkage_dataframe.rds', replogle_save_path)







# First describe which plots we need
# 
# EDA? not really using original counts
#    - maybe showing data cleaning/picking genes and grnas
# Unshrunk Estimates using SCEPTRE
#    - train, test, all
# Matrix Approximations
#    - 
# Shrunk Estimates 
#    - as matrices
#    - original Estimates to shrinkage point (with CIs)
#    - CI lengths before and after
#    - hist or qqplot (against unif) of inverted pvals
#    - (no) There is no 'MSE', but if we assume most effects are 0, MSE then?
#    - 

# //////////////////////////////////////////////////////////////////////////////////////////////////
# ================================== SCRIPT PARAMETERS =============================================
# //////////////////////////////////////////////////////////////////////////////////////////////////

# sceptre_save_path = '../../saves/sceptre/replogle/' # location where sceptre results are located
replogle_save_path= '../../saves/replogle/'         # location to save replogle approximations/shrinkage/etc. results
plot_path         = '../../plots/replogle/EBCI/'    # location to save plots of EBCI analysis on replogle dataset


assertthat::assert_that(dir.exists(replogle_save_path))
assertthat::assert_that(dir.exists(plot_path))

dir.create(sprintf('%s/final/', plot_path)) # where to put more finalized plots



# //////////////////////////////////////////////////////////////////////////////////////////////////
# ================================== LOAD ==========================================================
# //////////////////////////////////////////////////////////////////////////////////////////////////
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
ggplot2::theme_set(theme_bw() + theme(plot.title = element_text(hjust = .5), strip.background = element_rect(fill = 'white')))

# scripts?


# data/values/objects
replogleShrinkageParams = readRDS(file = sprintf('%s/replogleShrinkageParams.rds', replogle_save_path))

gene_index = read.csv(sprintf('%s/gene_index.csv', replogle_save_path))
grna_index = read.csv(sprintf('%s/grna_index.csv', replogle_save_path))

df_shrunk   = readRDS(sprintf('%s/EBCI_shrinkage_dataframe.rds', replogle_save_path))
df_unshrunk = readRDS(sprintf('%s/EBCI_unshrunk_dataframe.rds', replogle_save_path))
df_unshrunk = merge(merge(df_unshrunk, gene_index, by = 'gene'),                            # add grna and gene indices (the x and y coords for plotting)
                    grna_index, by = 'grna')
df_unshrunk$unshrunk_CIlength = 2*qnorm(1 - replogleShrinkageParams$ALPHA/2)*df_unshrunk$se # add 'original' CI lengths






head(df_shrunk)
head(df_unshrunk)










# //////////////////////////////////////////////////////////////////////////////////////////////////
# ================================== REPLOGLE COLORS =============================================
# //////////////////////////////////////////////////////////////////////////////////////////////////




matapprox_methods_nicenames = c( "softImpute"                 =  "softImpute", 
                                 "SVD"                        =  "SVD", 
                                 "sparseSVD"                  =  "Sparse SVD", 
                                 "sparseSVD_autoparams"       =  "Sparse SVD (CV params)", 
                                 "spectralbiclust"            =  "Spectral Biclustering", 
                                 "spectralbiclust_threshold"  =  "Spectral BiClustering (w threshold)", 
                                 "zeros"                      =  "Zero", 
                                 "average"                    =  "Average" )
rank_nicenames = paste0('(rank=', 1:100, ')')


methodrank_nicenames_replogle <- function(method_name, rank_) {
  # method_name = 'zeros'
  # rank_ = 2
  # rank_ = NA
  # methodrank_nicenames(method_name = plot_df_summ$method[1], rank_ = plot_df_summ$rank[1]) |> unname()
  # mapply(FUN = methodrank_nicenames, method_name =  plot_df_summ$method[1:4], rank_ = plot_df_summ$rank[1:4]) |> unname()
  
  if(is.na(rank_)) {
    return(matapprox_methods_nicenames[method_name] |> unname())
  } else {
    return(paste0(matapprox_methods_nicenames[method_name], ' ', rank_nicenames[rank_]) |> unname())
  }
  
}






#' Create a vector of hex code colors from white to main_color of length number_of_colors
#' @param main_color (character) the hex code of the main color 
#' @param number_of_colors (integer) number of colors to return
create_colors_replogle <- function(main_color, number_of_colors) {
  colorRampPalette(c("white", main_color))(number_of_colors + 1)[1:number_of_colors+1]
}


#'  use same as simEBCI?
create_color_pallete_replogle <- function(ranks) {
  distinct_colors = paletteer::paletteer_d("colorBlindness::paletteMartin")  # library(paletteer)
  distinct_colors[c(2, 4, 6, 7, 12, 13, 14)]
  
  
  
  # methodrank_colors = c("Original (Full Dataset)" = colorRampPalette(c("white", distinct_colors[2]))(2 + 1)[2],
  #                       "Original"                = colorRampPalette(c("white", distinct_colors[2]))(2 + 1)[3])
  
  methodrank_colors = c('Unshrunk' = distinct_colors[2])
  
  # softImpute
  temp_colors = create_colors_replogle(main_color = distinct_colors[4], number_of_colors = length(ranks)) # softImpute
  names(temp_colors) = mapply(FUN = methodrank_nicenames_replogle, method_name =  "softImpute", rank_ = ranks) |> unname()
  methodrank_colors = c(methodrank_colors, temp_colors) ; rm(temp_colors)
  
  # matdecomp_svd
  temp_colors = create_colors_replogle(main_color = distinct_colors[6], number_of_colors = length(ranks)) # matdecomp_svd
  names(temp_colors) = mapply(FUN = methodrank_nicenames_replogle, method_name =  "SVD", rank_ = ranks) |> unname()
  methodrank_colors = c(methodrank_colors, temp_colors) ; rm(temp_colors)
  
  # matdecomp_sparsesvd
  temp_colors = create_colors_replogle(main_color = distinct_colors[7], number_of_colors = length(ranks)) # matdecomp_sparsesvd
  names(temp_colors) = mapply(FUN = methodrank_nicenames_replogle, method_name =  "sparseSVD", rank_ = ranks) |> unname()
  methodrank_colors = c(methodrank_colors, temp_colors) ; rm(temp_colors)
  
  # matdecomp_sparsesvd autoparams
  temp_colors = create_colors_replogle(main_color = distinct_colors[9], number_of_colors = length(ranks)) # matdecomp_sparsesvd
  names(temp_colors) = mapply(FUN = methodrank_nicenames_replogle, method_name =  "sparseSVD_autoparams", rank_ = ranks) |> unname()
  methodrank_colors = c(methodrank_colors, temp_colors) ; rm(temp_colors)
  
  # spectralbiclust
  temp_colors = create_colors_replogle(main_color = distinct_colors[12], number_of_colors = length(ranks)) # spectralbiclust
  names(temp_colors) = mapply(FUN = methodrank_nicenames_replogle, method_name =  "spectralbiclust", rank_ = ranks) |> unname()
  methodrank_colors = c(methodrank_colors, temp_colors) ; rm(temp_colors)
  
  # spectralbiclust_threshold
  temp_colors = create_colors_replogle(main_color = distinct_colors[13], number_of_colors = length(ranks)) # spectralbiclust_threshold
  names(temp_colors) = mapply(FUN = methodrank_nicenames_replogle, method_name =  "spectralbiclust_threshold", rank_ = ranks) |> unname()
  methodrank_colors = c(methodrank_colors, temp_colors) ; rm(temp_colors)
  
  
  methodrank_colors = c(methodrank_colors, 
                        "Zero" = colorRampPalette(c("white", distinct_colors[14]))(2 + 1)[2], # zeros/avg
                        "Average" = colorRampPalette(c("white", distinct_colors[14]))(2 + 1)[3])
  methodrank_colors
  
}



# //////////////////////////////////////////////////////////////////////////////////////////////////
# ================================== SETUP DATAFRAME =============================================
# //////////////////////////////////////////////////////////////////////////////////////////////////







# df for all the matrices 
df = merge(df_unshrunk, df_shrunk, by = c('gene', 'grna')) 
df = df |> relocate(unshrunk_CIlength, .after = last_col()) |> mutate(shrunk_CIlength = upper_ci - lower_ci)


# create colors: for each approx method- then change saturation by rank
ranks = unique(df_shrunk$rank)                # get unique ranks present in df
ranks = ranks[!is.na(ranks)] # rm NA ranks
ranks = sort(ranks)

methodrank_nicenames_order = c()
for(cur_method in names(matapprox_methods_nicenames)) { # requires declaration of this list/vector (this is defined later in this file, right before the function methodrank_nicenames is defined)
  for(cur_rank in c(NA, ranks)) {
    methodrank_nicenames_order = c(methodrank_nicenames_order, methodrank_nicenames_replogle(method_name = cur_method, rank_ = cur_rank) |> unname())
  }
}


# make some cols as factors with nice names?
df$methodrank = mapply(FUN = methodrank_nicenames_replogle, method_name =  as.character(df$approx_method), rank_ = df$rank) |> unname()
df = df |> mutate(methodrank    = factor(   methodrank, levels = methodrank_nicenames_order)) |> relocate(methodrank, .after = 'rank')
df = df |> mutate(approx_method = factor(approx_method, levels = names(matapprox_methods_nicenames), labels = matapprox_methods_nicenames))
df = df |> mutate(split_type    = factor(   split_type,  levels = c('nosamplesplit', 'samplesplit'), labels = c('No Sample Split', 'Sample Split')))

methodrank_colors = create_color_pallete_replogle(ranks = ranks)


# //////////////////////////////////////////////////////////////////////////////////////////////////
# ================================== CI Lengths ====================================================
# //////////////////////////////////////////////////////////////////////////////////////////////////

assertthat::assert_that(dir.exists(sprintf('%s/final/', plot_path)))
dir.create(sprintf('%s/final/CI/', plot_path))


# CI lengths before and after
df_sample = df |>
              # filter(approx_method != 'Sparse SVD (CV params)') |> 
              group_by(split_type, approx_method, rank) |> sample_n(size=10) # if too many points, just do a few for plotting
p_CI = ggplot(df_sample,
       aes(x = unshrunk_CIlength, y = shrunk_CIlength,
           # color = approx_method, 
           # alpha = rank, 
           # group = paste0(approx_method, rank)
           color = methodrank
           )) + 
  # geom_point(key_glyph = 'rect') + 
  geom_line(key_glyph = 'rect', alpha = .85) +
  geom_abline(aes(slope = 1, intercept = 0)) +
  coord_cartesian(xlim = c(0, 1.3), ylim = c(.2, 1), 
                  expand = c(0, 0)) +
  facet_grid(cols = vars(split_type)) +
  scale_color_discrete(palette = methodrank_colors[names(methodrank_colors) %in% df_sample$methodrank]) +
  labs(title = 'CI lengths before and after shrinkage', 
       x = 'Unshrunk CI Length', y = 'Shrunk CI Length', color = 'Matrix Approximation Method') + #, alpha = 'rank')
  theme(panel.grid.minor  = element_blank(), 
        strip.text = element_text(size = 12),
        plot.title = element_text(size = 14),
        legend.text = element_text(size = 6),
        legend.key.size = unit(.4, 'cm')
        )

p_CI
ggsave(filename = sprintf('%s/final/CI/CIlengths.pdf', plot_path), plot = p_CI, width = 10, height = 4)
ggsave(filename = sprintf('%s/final/CI/CIlengths.png', plot_path), plot = p_CI, width = 10, height = 4, dpi = 300)






# //////////////////////////////////////////////////////////////////////////////////////////////////
# ================================== Point Plots =============================================
# //////////////////////////////////////////////////////////////////////////////////////////////////









