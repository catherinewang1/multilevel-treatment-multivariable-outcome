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
  distinct_colors[c(2, 4, 6, 7, 9, 12, 13, 14)]
  
  
  
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
df = merge(df_unshrunk, df_shrunk, by = c('gene', 'grna', 'split_type')) 
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
              group_by(split_type, approx_method, rank) |> sample_n(size=400) # if too many points, just do a few for plotting
# p_CI = ggplot(df_sample,
#        aes(x = unshrunk_CIlength, y = shrunk_CIlength,
#            # color = approx_method, 
#            # alpha = rank, 
#            # group = paste0(approx_method, rank)
#            color = methodrank
#            )) + 
#   # geom_point(key_glyph = 'rect') + 
#   geom_line(key_glyph = 'rect', alpha = .85) +
#   geom_abline(aes(slope = 1, intercept = 0)) +
#   coord_cartesian(xlim = c(0, 1.3), 
#                   ylim = c(0, 1.3),
#                   # ylim = c(.2, 1), 
#                   expand = c(0, 0)) +
#   facet_grid(cols = vars(split_type)) +
#   scale_color_discrete(palette = methodrank_colors[names(methodrank_colors) %in% df_sample$methodrank]) +
#   labs(title = 'CI lengths before and after shrinkage', 
#        x = 'Unshrunk CI Length', y = 'Shrunk CI Length', color = 'Matrix Approximation Method') + #, alpha = 'rank')
#   theme(panel.grid.minor  = element_blank(), 
#         strip.text = element_text(size = 12),
#         plot.title = element_text(size = 14),
#         legend.text = element_text(size = 6),
#         legend.key.size = unit(.4, 'cm')
#         )
# 
# p_CI
# ggsave(filename = sprintf('%s/final/CI/CIlengths.pdf', plot_path), plot = p_CI, width = 10, height = 4)
# ggsave(filename = sprintf('%s/final/CI/CIlengths.png', plot_path), plot = p_CI, width = 10, height = 4, dpi = 300)
# 



method_colors = methodrank_colors[grepl('(rank=30)|(Zero)|(Average)', names(methodrank_colors))]
names(method_colors) = gsub(pattern = " \\(rank=30\\)", replacement = "", x = names(method_colors))
p_CI = ggplot(df_sample |> mutate(rank = as.factor(sapply(rank, FUN = function(r){if(is.na(r)){return(30)}else{return(r)}}))),
       aes(x = unshrunk_CIlength, y = shrunk_CIlength,
           group = methodrank,
           color = approx_method, 
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
  scale_color_discrete(palette = method_colors) +
  scale_x_continuous(breaks = c(0, .25, .5, .75, 1, 1.25),
                     labels = c('0', .25, .5, .75, 1.00, 1.25)) + 
  labs(title = 'CI lengths before and after shrinkage', 
       x = 'Unshrunk CI Length', y = 'Shrunk CI Length', color = 'Matrix\nApproximation\nMethod') + #, alpha = 'rank')
  theme(panel.grid.minor  = element_blank(), 
        strip.text = element_text(size = 12),
        plot.title = element_text(size = 16),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        legend.key.size = unit(.4, 'cm')
  )

p_CI
ggsave(filename = sprintf('%s/final/CI/CIlengths_beforeafter.pdf', plot_path), plot = p_CI, width = 8, height = 4)
ggsave(filename = sprintf('%s/final/CI/CIlengths_beforeafter.png', plot_path), plot = p_CI, width = 8, height = 4, dpi = 300)





# //////////////////////////////////////////////////////////////////////////////////////////////////
# ================================== CI violinplots ================================================
# //////////////////////////////////////////////////////////////////////////////////////////////////

# add unshrunk CI lengths
method_colors_w_unshrunk = c(methodrank_colors['Unshrunk'], method_colors)
df_sample_w_unshrunk = rbind(df_sample, 
                             df_sample |> filter(approx_method == 'Zero') |> mutate(approx_method = 'Unshrunk', 
                                                       methodrank = 'Unshrunk',
                                                       shrunk_CIlength = unshrunk_CIlength, 
                                                       shrunk_value=NA, lower_ci=NA, upper_ci=NA, w_eb=NA, ebci_pvals=NA)
)
df_sample_w_unshrunk = df_sample_w_unshrunk |> mutate(methodrank    = factor(   methodrank, levels = c('Unshrunk', methodrank_nicenames_order))) |>
                                               mutate(approx_method = factor(approx_method, levels = names(method_colors_w_unshrunk)))
ggplot(df_sample_w_unshrunk |> mutate(rank = sapply(rank, FUN = function(r){if(is.na(r)){return('30')}else{return(r)}}) |> 
                                             factor(levels = sort(unique(df$rank)))),
       aes(x = methodrank,
           y = shrunk_CIlength,
           group = methodrank,
           fill = approx_method, 
           alpha = rank
       )) + 
  # geom_point(key_glyph = 'rect') + 
  # geom_histogram(aes(y = after_stat(density)), key_glyph = 'rect', position = 'dodge', binwidth = .2) +
  # geom_density(alpha = .4) +
  geom_violin(bounds = c(0, Inf), width=1.4, position='identity') +
  coord_cartesian(
                  ylim = c(0, 2),
                  # ylim = c(.2, 1),
                  expand = c(0, 0)) +
  # facet_grid(cols = vars(split_type, approx_method), rows = vars(rank)) +
  # facet_grid(cols = vars(split_type), rows = vars(split_type)) +
  facet_wrap(facets = vars(split_type), nrow = 2) +
  scale_fill_discrete(palette = method_colors_w_unshrunk) +
  scale_alpha_discrete(guide = 'none') +
  scale_y_continuous(breaks = seq(from = 0, to = 2, by = .5)) +
  labs(title = 'CI Lengths', 
       x = 'Method', y = 'CI Length', fill = 'Matrix\nApproximation\nMethod') + #, alpha = 'rank')
  theme(panel.grid.minor  = element_blank(), 
        panel.grid.major.x = element_blank(),
        strip.text = element_text(size = 12),
        plot.title = element_text(size = 16),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        legend.key.size = unit(.4, 'cm'), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)
  )


ggsave(filename = sprintf('%s/final/CI/CIlengths_violinplot.pdf', plot_path), width = 10, height = 9)
ggsave(filename = sprintf('%s/final/CI/CIlengths_violinplot.png', plot_path), width = 10, height = 9, dpi = 300)


# only display 3 ranks
ggplot(df_sample_w_unshrunk |> mutate(rank = sapply(rank, FUN = function(r){if(is.na(r)){return('30')}else{return(r)}}) |> 
                                        factor(levels = sort(unique(df$rank)))) |> 
         filter(rank %in% c( '10', '20', '30')),
       aes(x = methodrank,
           y = shrunk_CIlength,
           group = methodrank,
           fill = approx_method, 
           alpha = rank
       )) + 
  # geom_point(key_glyph = 'rect') + 
  # geom_histogram(aes(y = after_stat(density)), key_glyph = 'rect', position = 'dodge', binwidth = .2) +
  # geom_density(alpha = .4) +
  geom_violin(bounds = c(0, Inf), width=1.4, position='identity') +
  coord_cartesian(
    ylim = c(0, 2),
    # ylim = c(.2, 1),
    expand = c(0, 0)) +
  # facet_grid(cols = vars(split_type, approx_method), rows = vars(rank)) +
  # facet_grid(cols = vars(split_type), rows = vars(split_type)) +
  facet_wrap(facets = vars(split_type), nrow = 2) +
  scale_fill_discrete(palette = method_colors_w_unshrunk) +
  scale_alpha_discrete(guide = 'none') +
  scale_y_continuous(breaks = seq(from = 0, to = 2, by = .5)) +
  labs(title = 'CI Lengths', 
       x = 'Method', y = 'CI Length', fill = 'Matrix\nApproximation\nMethod') + #, alpha = 'rank')
  theme(panel.grid.minor  = element_blank(), 
        panel.grid.major.x = element_blank(),
        strip.text = element_text(size = 12),
        plot.title = element_text(size = 16),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        legend.key.size = unit(.4, 'cm'), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)
  )


ggsave(filename = sprintf('%s/final/CI/CIlengths_violinplot2.pdf', plot_path), width = 10, height = 9)
ggsave(filename = sprintf('%s/final/CI/CIlengths_violinplot2.png', plot_path), width = 10, height = 9, dpi = 300)


# //////////////////////////////////////////////////////////////////////////////////////////////////
# ================================== EBCI pvals (maybe put in separate file 4.5...) =============================================
# //////////////////////////////////////////////////////////////////////////////////////////////////
# sample p-value better (e.g. not just random from each group. pick the same AY tests)



method_colors_w_unshrunk = c(methodrank_colors['Unshrunk'], method_colors)


set.seed(12345)
subset_AYtests = df_unshrunk |> filter(split_type == 'samplesplit') |> sample_n(size = 300) |> select(gene, grna)

# construct sample df
df_sample_shrunk = merge(subset_AYtests, df, by=c('grna', 'gene'), all.x=TRUE, all.y=FALSE) |>
                    dplyr::rename(pval = ebci_pvals)

# add original pvals (could go back to SCEPTRE p-value, but this take more effort... some code below for it)
df_sample_unshrunk = merge(subset_AYtests, df_unshrunk, by=c('grna', 'gene'), all.x=TRUE, all.y=FALSE)  |>
                     dplyr::mutate(pval = 2*(1 - pnorm(abs(unshrunk_value / se) ))) |>
                     dplyr::mutate(approx_method = 'Unshrunk', 
                                   rank          = NA,
                                   methodrank    = 'Unshrunk') |> 
                     dplyr::mutate(split_type    = factor(split_type,  levels = c('nosamplesplit', 'samplesplit'), 
                                                                       labels = c('No Sample Split', 'Sample Split')))


df_sample = dplyr::bind_rows(df_sample_shrunk, df_sample_unshrunk)

# colnames(df_sample_shrunk)
# colnames(df_sample_unshrunk)
# colnames(df_sample)
# head(df_sample)
# df_sample |> filter(approx_method == 'Unshrunk')


df_sample = df_sample |> mutate(methodrank    = factor(   methodrank, levels = c('Unshrunk', methodrank_nicenames_order))) |>
                         mutate(approx_method = factor(approx_method, levels = names(method_colors_w_unshrunk)))
ggplot(df_sample |> mutate(rank = sapply(rank, FUN = function(r){if(is.na(r)){return('30')}else{return(r)}}) |> 
                                        factor(levels = sort(unique(df$rank)))),
       aes(x = methodrank,
           # y = -log(pval+.00001),
           y = pval,
           group = methodrank,
           fill = approx_method, 
           alpha = rank
       )) + 
  # geom_point(key_glyph = 'rect') + 
  # geom_histogram(aes(y = after_stat(density)), key_glyph = 'rect', position = 'dodge', binwidth = .2) +
  # geom_density(alpha = .4) +
  geom_violin(#bounds = c(0, Inf), 
              bounds = c(0, 1),
              width=1.8, position='identity') +
  coord_cartesian(
    ylim = c(0, 1),
    # ylim = c(.2, 1),
    expand = c(0, 0)) +
  # facet_grid(cols = vars(split_type, approx_method), rows = vars(rank)) +
  # facet_grid(cols = vars(split_type), rows = vars(split_type)) +
  facet_wrap(facets = vars(split_type), nrow = 2) +
  scale_fill_discrete(palette = method_colors_w_unshrunk) +
  scale_alpha_discrete(guide = 'none') +
  # scale_y_continuous(breaks = seq(from = 0, to = 2, by = .5)) +
  labs(title = 'p-values', 
       x = 'Method', y = 'p-value', fill = 'Matrix\nApproximation\nMethod') + #, alpha = 'rank')
  theme(panel.grid.minor  = element_blank(), 
        panel.grid.major.x = element_blank(),
        strip.text = element_text(size = 12),
        plot.title = element_text(size = 16),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        legend.key.size = unit(.4, 'cm'), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)
  )


ggsave(filename = sprintf('%s/final/CI/pvals_violinplot.pdf', plot_path), width = 14, height = 8)
ggsave(filename = sprintf('%s/final/CI/pvals_violinplot.png', plot_path), width = 14, height = 8, dpi = 300)





ggplot(df_sample |> mutate(rank = sapply(rank, FUN = function(r){if(is.na(r)){return('30')}else{return(r)}}) |> 
                             factor(levels = sort(unique(df$rank)))) |>
                    filter(rank %in% c('5', '10', '30')),
       aes(x = methodrank,
           # y = -log(pval+.00001),
           y = pval,
           group = methodrank,
           fill = approx_method, 
           alpha = rank
       )) + 
  # geom_point(key_glyph = 'rect') + 
  # geom_histogram(aes(y = after_stat(density)), key_glyph = 'rect', position = 'dodge', binwidth = .2) +
  # geom_density(alpha = .4) +
  geom_violin(#bounds = c(0, Inf), 
    bounds = c(0, 1),
    width=1.6, position='identity') +
  coord_cartesian(
    # ylim = c(0, 1),
    ylim = c(0, .3),
    expand = c(0, 0)) +
  # facet_grid(cols = vars(split_type, approx_method), rows = vars(rank)) +
  # facet_grid(cols = vars(split_type), rows = vars(split_type)) +
  facet_wrap(facets = vars(split_type), nrow = 2) +
  scale_fill_discrete(palette = method_colors_w_unshrunk) +
  scale_alpha_discrete(guide = 'none') +
  # scale_y_continuous(breaks = seq(from = 0, to = 2, by = .5)) +
  labs(title = 'p-values', 
       x = 'Method', y = 'p-value', fill = 'Matrix\nApproximation\nMethod') + #, alpha = 'rank')
  theme(panel.grid.minor  = element_blank(), 
        panel.grid.major.x = element_blank(),
        strip.text = element_text(size = 12),
        plot.title = element_text(size = 16),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        legend.key.size = unit(.4, 'cm'), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)
  )

ggsave(filename = sprintf('%s/final/CI/pvals_violinplot2.pdf', plot_path), width = 14, height = 8)
ggsave(filename = sprintf('%s/final/CI/pvals_violinplot2.png', plot_path), width = 14, height = 8, dpi = 300)








# //////////////////////////////////////////////////////////////////////////////////////////////////
# ================================== Trash =========================================================
# //////////////////////////////////////////////////////////////////////////////////////////////////

if(F) {
  # faceted hist of pvals
  ggplot(df_sample |>  mutate(rank = sapply(rank, FUN = function(r){if(is.na(r)){return('30')}else{return(r)}}) |> 
                                factor(levels = sort(unique(df$rank)))) |> 
           filter(rank %in% c('5', '30')),
         aes(x = pval,
             fill = approx_method, 
             alpha = rank)) +
    geom_histogram(aes(y = after_stat(density))) +
    facet_grid(cols = vars(split_type), rows = vars(methodrank)) +
    scale_fill_discrete(palette = method_colors_w_unshrunk) +
    scale_alpha_discrete(guide = 'none') 
}





# there is the SCEPTRE p-value, and then the slighlty regularized (SCEPTRE p-value)-->se-->std norm p-value
# so we could compare to the original (probably the better p-value for analysis)
# or we could compare to the slightly regularized p-value (better to comparison to EBCI)

# //////////////////////////////////////////////////////////////////////////////////////////////////
# ================================== > Load SCEPTRE p-values =========================================
# //////////////////////////////////////////////////////////////////////////////////////////////////
use_SCEPTRE_pvals = FALSE
# have to do some work to load in SCEPTRE p-values
if(use_SCEPTRE_pvals) {
  




sceptre_save_path = '../../saves/sceptre/replogle/' # location where sceptre results are located
replogle_save_path= '../../saves/replogle/'         # location to save replogle approximations/shrinkage/etc. results
plot_path         = '../../plots/replogle/EBCI/'    # location to save plots of EBCI analysis on replogle dataset


gene_index = read.csv(sprintf('%s/gene_index.csv', replogle_save_path))
grna_index = read.csv(sprintf('%s/grna_index.csv', replogle_save_path))


effects = list()
for(split in c('all', 'train', 'test')) {
  sceptre_obj = readRDS(sprintf('%s/sceptre_obj_%s.rds', sceptre_save_path, split))
  effects[[split]] = dplyr::bind_rows(
    sceptre_obj@calibration_result |> mutate(test = 'negative')  |> relocate(test),
    sceptre_obj@power_result       |> mutate(test = 'positive')  |> relocate(test),
    sceptre_obj@discovery_result   |> mutate(test = 'discovery') |> relocate(test)) |>
    mutate(estimate = log_2_fold_change) |> 
    rename(pvalue = p_value,
           grna = grna_target,
           gene = response_id) |>
    mutate(se = mapply(FUN = spline_se, 
                       p   = pvalue, 
                       mu  = estimate) ) |> 
    mutate(tstat = estimate / se) # make a tstat column = estimate / se
  rm(sceptre_obj); gc()
}



sceptre_pvals = list()
for(split_type in c('samplesplit', 'nosamplesplit')) {
  effects_df_split = effects[[split]] |> dplyr::filter((test == 'discovery') & (gene %in% gene_index$gene) & (grna %in% grna_index$grna)) 
  
  # # following code from 2_replogleShrinkage.r... but we don't need to put it into matrix form...
  # # pval 
  # # matrix form: grna x gene
  # matrices_se_split = effects_df_split  |> 
  #   dplyr::select(grna, gene, p) |>
  #   tidyr::pivot_wider(names_from = gene, values_from = p) |>
  #   tibble::column_to_rownames(var='grna') |> as.matrix()
  # # order based on grna_index and gene_index
  # p_matrices[[split]] = matrices_se_split[ grna_index |> arrange(grna_idx) |> pull(grna), 
  #                                           gene_index |> arrange(gene_idx) |> pull(gene)]
  # 
  
  
  
  sceptre_pvals[[split]] = effects_df_split |> dplyr::select(grna, gene, p) |> dplyr::rename(sceptre_pvals = p)
  
}


}










