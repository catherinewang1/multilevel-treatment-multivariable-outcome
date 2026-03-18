




# === Libraries
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(RColorBrewer))
# suppressPackageStartupMessages(library(ebci)) # robust emp shrinkage
suppressPackageStartupMessages(library(assertthat)) 
source('../utils/matrix_shrinkage.r')


plot_folder = '../plots/matrixSampSplit/shrinkMatrix/'
save_folder = '../saves/matrixSampSplit/shrinkMatrix/'

dir.create(sprintf('%s/heatmaps/', plot_folder))


# === Load
shrink_dfs = list()
shrink_dfs[['lowrank']] = list()
shrink_dfs[['sparseSVD']] = list()
for(filename in list.files(save_folder)) {
  
  approxmethod = strsplit(x = filename, split = '_')[[1]][[2]]
  rank = gsub(pattern = '\\D', x = filename, replacement='') |> as.numeric()

  shrink_dfs[[approxmethod]][[rank]] = read.csv(sprintf("%s%s", save_folder, filename))
}


# === Matrix Heatmaps ===

# Plot shrinkage results as heatmap matrices


# get an ordering for genes/grna (marginal clustering)
shrink_df = shrink_dfs[['sparseSVD']][[5]]

shrink_matrix = shrink_df |> # filter(grna_idx %in% 1:400 & gene_idx %in% 1:400) |> 
  dplyr::select(grna, gene, shrunk_value) |>
  tidyr::pivot_wider(names_from = gene, values_from = shrunk_value) |>
  tibble::column_to_rownames(var='grna')

mat = shrink_df |> # filter(grna_idx %in% 1:400 & gene_idx %in% 1:400) |> 
  dplyr::select(grna, gene, unshrunk_value) |>
  tidyr::pivot_wider(names_from = gene, values_from = unshrunk_value) |>
  tibble::column_to_rownames(var='grna')


matscaled = as.matrix(scale(mat))
matscaled[is.nan(matscaled)] = 0
row_order = hclust(dist(matscaled))$order
column_order = hclust(dist(t(matscaled)))$order
grna_index = data.frame(grna = rownames(mat)[row_order],
                        grna_idx  = 1:nrow(mat))
gene_index = data.frame(gene = colnames(mat)[column_order],
                        gene_idx  = 1:ncol(mat))




# color breaks for plotting
color_limits = c(-2, 2)
color_breaks = sort(union(color_limits, seq(from = round(color_limits[1]), to = round(color_limits[2]))))
color_breaks_label = color_breaks
color_breaks_label[which.min(color_breaks)] = sprintf('<%.1f', min(color_breaks))
color_breaks_label[which.max(color_breaks)] = sprintf('>%.1f', max(color_breaks))




# create plot_df
plot_df = NULL
for(approxmethod in c('lowrank', 'sparseSVD')) {
  for(rank in c(5, 10, 15, 30)) {
    plot_df = rbind(plot_df, 
                    shrink_dfs[[approxmethod]][[rank]] |> 
                      dplyr::mutate(approxmethod = approxmethod,
                                    rank = sprintf('rank=%02.f', rank)))
  }
  # add original estimates (all the loaded dfs should have these) (2x for plotting)
  plot_df = rbind(plot_df,
                  shrink_dfs[[approxmethod]][[rank]] |>
                    dplyr::mutate(approxmethod = approxmethod,
                                  rank = 'unshrunk') |>
                    dplyr::mutate(shrinkage_point = unshrunk_value,   # set all values to the original estimates
                                  shrunk_value    = unshrunk_value, 
                                  lower_ci        = unshrunk_value - qnorm(.95) * se,
                                  upper_ci        = unshrunk_value + qnorm(.95) * se) 
  )
  
}

# add idx's- (from gene or grna name --> number on the axis)
plot_df = merge(gene_index, plot_df, by = 'gene')
plot_df = merge(grna_index, plot_df, by = 'grna')

plot_df$rank = factor(plot_df$rank, levels = sort(unique(plot_df$rank), decreasing = TRUE))  





# plot
p = ggplot(plot_df) +
  geom_raster(aes(x = grna_idx, y = gene_idx, fill = shrunk_value)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x='grna', y = 'gene', fill = NULL, title = 'Shrunk Value') +
  scale_fill_gradient2(limits = color_limits, # set color limits
                       oob=scales::squish, # if outside lims, set to limits
                       midpoint = 0,
                       high = myRed, low = myBlue, mid = 'white',
                       # low  = brewer.pal(n = 9, name = "RdBu")[9],
                       # high = brewer.pal(n = 9, name = "RdBu")[1],
                       breaks = color_breaks,
                       labels = color_breaks_label) +
  facet_grid(rows = vars(approxmethod),
             cols = vars(rank)) +
  theme_bw() +
  theme(strip.background = element_rect(fill = 'white'), 
        panel.spacing = unit(.2, 'lines'),
        axis.ticks = element_blank(), 
        axis.text = element_blank(), 
        legend.position = 'bottom', 
        # legend.key.size = unit(.5, 'cm'),
        legend.key.height = unit(.3, 'cm'),
        legend.key.width  = unit(1.75, 'cm'),
        legend.text = element_text(size = 7))


ggsave(plot = p, filename = sprintf('%sheatmaps/shrunkvalue.pdf', plot_folder), height = 7, width = 7)


p = ggplot(plot_df) +
  geom_raster(aes(x = grna_idx, y = gene_idx, fill = shrinkage_point)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x='grna', y = 'gene', fill = NULL, title = 'Shrinkage Point') +
  scale_fill_gradient2(limits = color_limits, # set color limits
                       oob=scales::squish, # if outside lims, set to limits
                       midpoint = 0,
                       high = myRed, low = myBlue, mid = 'white',
                       # low  = brewer.pal(n = 9, name = "RdBu")[9],
                       # high = brewer.pal(n = 9, name = "RdBu")[1],
                       breaks = color_breaks,
                       labels = color_breaks_label) +
  facet_grid(rows = vars(approxmethod),
             cols = vars(rank)) +
  theme_bw() +
  theme(strip.background = element_rect(fill = 'white'), 
        panel.spacing = unit(.2, 'lines'),
        axis.ticks = element_blank(), 
        axis.text = element_blank(), 
        legend.position = 'bottom', 
        # legend.key.size = unit(.5, 'cm'),
        legend.key.height = unit(.3, 'cm'),
        legend.key.width  = unit(1.75, 'cm'),
        legend.text = element_text(size = 7))

ggsave(plot = p, filename = sprintf('%sheatmaps/shrinkagepoint.pdf', plot_folder), height = 7, width = 7)

p = ggplot(plot_df |> mutate(isSig = as.numeric(!(lower_ci <= 0 & 0 <= upper_ci)))) +
  geom_raster(aes(x = grna_idx, y = gene_idx, fill = isSig)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x='grna', y = 'gene', fill = NULL, title = 'Significant at alpha=.1') +
  scale_fill_gradient(
    high = myRed, low = 'white',
    # low  = brewer.pal(n = 9, name = "RdBu")[9],
    # high = brewer.pal(n = 9, name = "RdBu")[1],
    breaks = c(0, 1),
    labels = c('not sig', 'sig')) +
  facet_grid(rows = vars(approxmethod),
             cols = vars(rank)) +
  theme_bw() +
  theme(strip.background = element_rect(fill = 'white'), 
        panel.spacing = unit(.2, 'lines'),
        axis.ticks = element_blank(), 
        axis.text = element_blank(), 
        legend.position = 'bottom', 
        # legend.key.size = unit(.5, 'cm'),
        legend.key.height = unit(.3, 'cm'),
        legend.key.width  = unit(1.75, 'cm'),
        legend.text = element_text(size = 7))

ggsave(plot = p, filename = sprintf('%sheatmaps/significant.pdf', plot_folder), height = 7, width = 7)






# === Shrinkage Changes ===
dir.create(sprintf('%spoints/', plot_folder))
set.seed(12345)
subsample_idx = sample(1:nrow(plot_df), 10000)

# plot shrunk vs unshrunk
p = ggplot(plot_df[subsample_idx, ] |> filter(rank != 'unshrunk')) +
  geom_abline(aes(slope = 1, intercept = 0), color = 'gray', alpha = .6) +
  geom_point(aes(x = unshrunk_value, y = shrunk_value, color = se)) +
  labs(title = 'Shrunk vs Unshrunk Estimates', x = 'unshrunk', y = 'shrunk') +
  coord_cartesian(xlim = c(-2, 2), ylim = c(-2, 2)) +
  facet_grid(rows = vars(approxmethod),
             cols = vars(rank)) +
  theme_classic()

ggsave(plot = p, filename = sprintf('%spoints/shrunkvsunshrunk.pdf', plot_folder), height = 6, width = 7)


# plot shrunk vs unshrunk centered
p = ggplot(plot_df[subsample_idx, ] |> filter(rank != 'unshrunk')) +
  geom_abline(aes(slope = 1, intercept = 0), color = 'gray') +
  geom_point(aes(x = unshrunk_value - shrinkage_point, y = shrunk_value - shrinkage_point, color = se)) +
  labs(title = 'Shrunk vs Unshrunk Estimates Centered', x = 'unshrunk', y = 'shrunk') +
  coord_cartesian(xlim = c(-2, 2), ylim = c(-2, 2)) +
  facet_grid(rows = vars(approxmethod),
             cols = vars(rank)) +
  theme_classic()

ggsave(plot = p, filename = sprintf('%spoints/shrunkvsunshrunkcentered.pdf', plot_folder), height = 6, width = 7)


# plot unshrunk vs shrinkage point
p = ggplot(plot_df[subsample_idx, ] |> filter(rank != 'unshrunk')) +
  geom_abline(aes(slope = 1, intercept = 0), color = 'gray') +
  geom_point(aes(x = shrunk_value, y = shrinkage_point, color = se)) +
  labs(title = 'Shrinkage Point vs Shrunk Estimates', x = 'unshrunk', y = 'shrunk') +
  coord_cartesian(xlim = c(-2, 2), ylim = c(-2, 2)) +
  facet_grid(rows = vars(approxmethod),
             cols = vars(rank)) +
  theme_classic()

ggsave(plot = p, filename =sprintf('%spoints/shrinkagepointvsshrunk.pdf', plot_folder), height = 6, width = 7)




# === Shrinkage CIs ===
set.seed(12345)
# get the same gene/grna tests across all methods and ranks (to be better visually, should do to prev plots too)
sample_tests = plot_df |> filter(rank == 'unshrunk' & approxmethod == 'lowrank') |> slice_sample(n = 15) |> select(grna, gene) |> mutate(x = 1:n())
plot_df_sample = merge(plot_df, sample_tests, all.x = FALSE, all.y = TRUE)

# plot_df_sample |> group_by(approxmethod, rank) |> summarize(count = n()) # should all be nrow(sample_tests)


p = ggplot(plot_df_sample # |> filter(rank == 'rank=15' & approxmethod == 'lowrank') #|> 
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
  geom_segment(aes(x = x + .25, y = lower_ci, yend = upper_ci),
               lineend = 'square', linewidth = 1, alpha = .7, color = 'deepskyblue4') +
  # --- Unshrunk CIs ---
  geom_segment(aes(x = x, 
                   y    = unshrunk_value - qnorm(.95)*se, 
                   yend = unshrunk_value + qnorm(.95)*se),
               lineend = 'square', linewidth = 1, alpha = .7, color = 'deepskyblue2') +
  # geom_segment(aes(y    = eval(parse(text = CIlower_colname)), 
  #                  yend = eval(parse(text = CIupper_colname))),
  #              lineend = 'square', linewidth = 1, alpha = .7, color = 'deepskyblue4') +
  
  # --- Arrow showing shrinkage
  geom_segment(aes(x = x, xend = x+.2,
                   y =    unshrunk_value, 
                   yend = .99 * shrunk_value + .01 * unshrunk_value),
               lineend = 'square', linejoin = 'bevel', arrow = arrow(length = unit(0.2,"cm")),
               linewidth = .3, alpha = .7, color = 'deepskyblue2') +
  # geom_segment(aes(y = eval(parse(text = unshrunk_colname)), 
  #                  yend = .99 * eval(parse(text = shrunk_colname)) + 
  #                    .01 * eval(parse(text = unshrunk_colname))),
  #              lineend = 'square', linejoin = 'bevel', arrow = arrow(length = unit(0.2,"cm")),
  #              linewidth = .3, alpha = .7, color = 'deepskyblue2') +
  # --- Unshrunk Point---
  geom_point(aes(x = x      , y = unshrunk_value), color = 'deepskyblue2') +
  # --- Shrunk Point ---
  geom_point(aes(x = x + .25, y = shrunk_value), shape=18, color = 'deepskyblue4') +
  # geom_point(aes(y = eval(parse(text = unshrunk_colname))), color = 'deepskyblue2') +
  # geom_point(aes(y = eval(parse(text =   shrunk_colname))), color = 'deepskyblue4') +
  scale_x_continuous(expand = c(0, 0), breaks = seq(0, 1000, by = 10),
                     # limits = c(1, 90)) +
                     # limits = c(65, 175)) +
                     # limits = c(95, 400)) +
                     # limits = c(1, 400)) +
  ) +
  # scale_y_continuous(expand = c(0.025, 0)) +
  coord_cartesian(ylim = c(-1, 1)) +
  labs(x = 'AY Test', y = 'Estimates',
       title = 'Before and After Robust EBCI') +
  facet_grid(rows = vars(approxmethod),
             cols = vars(rank)) +
  theme_bw() +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        # legend.position = 'inside',
        # legend.position.inside = legend_position,
        legend.background = element_rect(color = 'black'))


ggsave(plot = p, filename = sprintf('%spoints/CIs.pdf', plot_folder), height = 6, width = 15)


# 90% CI lengths before and after 
set.seed(12345)
p = ggplot(plot_df[sample(x=1:nrow(plot_df), size = 10000), ] |> 
             mutate(unshrunk_ci_length = 2*qnorm(.95)*se,
                    shrunk_ci_length = upper_ci - lower_ci)) +
  geom_abline(aes(intercept = 0, slope = 1)) +
  geom_point(aes(x = unshrunk_ci_length, y = shrunk_ci_length), alpha = .8) +
  labs(x = 'Unshrunk CI Lengths', y = 'Shrunk CI Lengths', title = 'CI Lengths') +
  facet_grid(rows = vars(approxmethod),
             cols = vars(rank)) +
  theme_classic()


ggsave(plot = p, filename =sprintf('%spoints/CIlengths.pdf', plot_folder), height = 6, width = 7)

