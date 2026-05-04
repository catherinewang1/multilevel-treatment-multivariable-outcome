# Use saved SCEPTRE results to perform EBCI shrinkage
# 
# 
# Requires previously saved SCEPTRE analysis results under <sceptre_save_path>/
# Specifically: <sceptre_save_path>/sceptre_obj_[train|test|all].rds 
# 
# Sections:
#   LOAD
#   PREPARE
#   INITIAL SCEPTRE RESULT PLOTS
#   PERFORM MATRIX APPROXIMATIONS
#   PERFORM EBCI SHRINKAGE
#   SHRINKAGE PLOTS
# 
# 
# Require previously saved objects: <sceptre_save_path>/
#  - sceptre_obj_test.rds
#  - sceptre_obj_train.rds
#  - sceptre_obj_all.rds
#
# Newly saved objects: <replogle_save_path>/
#  - gene_index.csv
#  - grna_index.csv
#  - estse_matrices.rds (matrices of estimates (log_2_fold_change) and se's)
#  - approxmatrices.rds (matrices of approximations of estimates)
#  - EBCI_shrinkage.rds (matrix of EBCI shrinkage results)
# Newly saved plots: <plot_path>/
#  - mat/
#  -     estimates/   Estimates (original estimates)
#  -     matapprox/   Approximating Matrices (matrix Approximations)
#  -     shrinkage/   Shrinkage Matrices (matrix of Shrunk estimates)
#  - shrinkage/  


# script params
sceptre_save_path = '../../saves/sceptre/replogle/' # location where sceptre results are located
replogle_save_path= '../../saves/replogle/'         # location to save replogle approximations/shrinkage/etc. results
plot_path         = '../../plots/replogle/EBCI/'    # location to save plots of EBCI analysis on replogle dataset
# dir.exists(sceptre_save_path); dir.exists(replogle_save_path); dir.exists(plot_path)



# Miscoverage level of the EBCI 
ALPHA = .1


# NUM_GENES = 1000; NUM_GRNAS = 50           # number of grna-gene pairs to test
NUM_GENES = 10; NUM_GRNAS = 5           # number of grna-gene pairs to test (Test code is running first)
NA_THRESH_GENE = .05; NA_THRESH_GRNA = .05 # threshold for cleaning grnas/genes

ORDER_METHOD = 'spectralsvd'   # method to create row/col ordering for plots


# should move this to the top of script
matapprox_methods          = c('softImpute', 'SVD', 'sparseSVD', 'sparseSVD_autoparams', 'spectralbiclust', 'spectralbiclust_threshold', 'zeros', 'average')
matapprox_methods_hasranks = c('softImpute'=TRUE, 'SVD'=TRUE, 'sparseSVD'=TRUE, 'spectralbiclust'=TRUE, 'spectralbiclust_threshold'=TRUE, 'zeros'=FALSE, 'average'=FALSE)
# ranks = c(1, 3, 5, 10, 20, 30)
ranks = c(1, 2)
matapprox_methods = matapprox_methods[c('SVD', 'sparseSVD', 'spectralbiclust', 'zeros')] # select a subset for testing




# libraries
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(ebci)) # robust emp shrinkage
suppressPackageStartupMessages(library(assertthat)) 

# our code (perform shrinkage + plots + etc...)
source('../../utils/matrix_approx_fns.r') # should move more functions here for general purpose
source('../../utils/matrix_shrinkage.r')






# ======================================================================================================================================
#                                                                               ========================================================
#       LOAD:                                                                   ========================================================
#                                                                               ========================================================
# ======================================================================================================================================

# --------------------------------------------------------------------------------------------------------------------------
# load saved effects
effects = list()
for(split in c('all', 'train', 'test')) {
  sceptre_obj = readRDS(sprintf('%s/sceptre_obj_%s.rds', sceptre_save_path, split))
  effects[[split]] = bind_rows(
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


# effects[['all']] |> head()
# effects[['all']] |> filter(test == 'discovery') |> group_by(gene) |> summarize(count = n()) |>pull(count) |> range()
# effects[['all']] |> filter(test == 'discovery') |> group_by(grna) |> summarize(count = n()) |>pull(count) |> range()
# 
# effects[['all']] |> filter(test == 'discovery') |> pull(estimate) |> is.na() |> mean()
# 
# effects[['all']] |> filter(test == 'discovery') |> filter(is.na(estimate))



# exclude some genes/grnas that have too many no qc pass/NA estimates



# This is stricter, this filters for genes > thresh & grna > thresh separately. 
# but what we reallt want is just some set of {(grna, gene)} that meet the threshold in this set
# viz

# effects[['all']] |> filter(test == 'discovery') |> group_by(gene) |> summarize(rateNA = mean(is.na(estimate))) |> pull(rateNA) |> hist(breaks=seq(0, 1, by = .02))
# effects[['all']] |> filter(test == 'discovery') |> group_by(grna) |> summarize(rateNA = mean(is.na(estimate))) |> pull(rateNA) |> hist(breaks=seq(0, 1, by = .02))
# 
# 
# exclude_gene = effects[['all']] |> filter(test == 'discovery') |> group_by(gene) |> summarize(NArate = mean(is.na(estimate))) |> filter(NArate > .02) |> pull(gene)
# exclude_grna = effects[['all']] |> filter(test == 'discovery') |> group_by(grna) |> summarize(NArate = mean(is.na(estimate))) |> filter(NArate > .02) |> pull(grna)
# effects_df = effects_df |> filter(!(grna %in% exclude_grna) & !(gene %in% exclude_gene))

# ======================================================================================================================================
#                                                                               ========================================================
#       PREPARE:                                                                ========================================================
#          1. select grnas and genes based on signal and NA                     ========================================================
#          2. order selected grnas and genes (some (bi)/clustering)             ========================================================
#          3. construct effect and se matrices subset and ordered               ========================================================
# ======================================================================================================================================
print(sprintf("[%s]     - prep", Sys.time()))

# TODO: Also filter tests based on too large of se's
# --------------------------------------------------------------------------------------------------------------------------
# 1. important genes and grnas selected by a 'score'
effects_df = effects[['all']] |> filter(test == 'discovery') # | test == 'positive') 
topgenes = effects_df |> group_by(gene) |> summarize(sum_effects = sum(abs(estimate), na.rm = TRUE),
                                                     sd_effects = sd(estimate, na.rm = TRUE),
                                                     sig_effects = sum(significant, na.rm=TRUE),
                                                     sum_tstats  = sum(abs(estimate/se), na.rm = TRUE),
                                                     NArate = mean(is.na(estimate))) |>
  mutate(score = scale(sum_effects) + (1/3)* scale(sd_effects) +  sig_effects) |> 
  arrange(desc(score)) |> 
  arrange(desc(sum_tstats))
topgrnas = effects_df |> group_by(grna) |> summarize(sum_effects = sum(abs(estimate), na.rm = TRUE),
                                                     sd_effects = sd(estimate, na.rm = TRUE),
                                                     sig_effects = sum(significant, na.rm=TRUE),
                                                     sum_tstats  = sum(abs(estimate/se), na.rm = TRUE),
                                                     NArate = mean(is.na(estimate))) |>
  mutate(score = scale(sum_effects) + (0)* scale(sd_effects) +  sig_effects) |> 
  arrange(desc(score)) |> 
  arrange(desc(sum_tstats))


# Do a few simple steps to filter out genes and grnas
# # moved to beginning, define script params at the beginning
# NUM_GENES = 1000; NUM_GRNAS = 50
# NA_THRESH_GENE = .05; NA_THRESH_GRNA = .05

## 1. Take top 2.5*NUM_GENES and 4*NUM_GRNAS (could change mult factor)
cur_genes = topgenes |> arrange(desc(score)) |> head(4*NUM_GENES) |> pull(gene)
cur_grnas = topgrnas |> arrange(desc(score)) |> head(4*NUM_GRNAS) |> pull(grna)
## 2. Remove those with marginal NA missing > threshold
curr_effects = effects_df |> filter(gene %in% cur_genes & grna %in% cur_grnas) 
cur_genes = curr_effects |> group_by(gene) |> summarise(NArate = mean(is.na(estimate)), count = n()) |> filter(NArate <= NA_THRESH_GENE) |> pull(gene)
cur_grnas = curr_effects |> group_by(grna) |> summarise(NArate = mean(is.na(estimate)), count = n()) |> filter(NArate <= NA_THRESH_GRNA) |> pull(grna)
## 3. Take top NUM_GENES and NUM_GRNAS
cur_genes = topgenes |> filter(gene %in% cur_genes) |> arrange(desc(score)) |> head(NUM_GENES) |> pull(gene)
cur_grnas = topgrnas |> filter(grna %in% cur_grnas) |> arrange(desc(score)) |> head(NUM_GRNAS) |> pull(grna)

chosen_genes = cur_genes
chosen_grnas = cur_grnas


if(F) {
  # try out method to get top genes and grnas that is aware of NA 
  num_genes = 1000; num_grnas = 20
  NA_thresh = .02
  cur_genes = topgenes |> arrange(desc(score)) |> head(num_genes) |> pull(gene)
  cur_grnas = topgenes |> arrange(desc(score)) |> head(num_grnas) |> pull(gene)
  
  # slow part...
  
  curr_effects = effects_df |> filter(gene %in% cur_genes & grna %in% cur_grnas) 
  
  curr_effects |> group_by(gene) |> summarise(NArate = mean(is.na(estimate)), count = n()) |> filter(NArate <= NA_thresh)
}

  
  
# --------------------------------------------------------------------------------------------------------------------------
# 2. order the selected grnas and genes 

# prep matrix of effect sizes
# effects_df = effects_df |> filter((gene %in% (topgenes |> arrange(desc(score)) |> head(2000) |> pull(gene)))
#                                   & 
#                                     (grna %in% (topgrnas |> arrange(desc(score)) |> head(50) |> pull(grna))))
effects_df = effects_df |> filter((gene %in% chosen_genes)
                                  & 
                                  (grna %in% chosen_grnas))


estimates_matrix = effects_df  |> # filter(grna_idx %in% 1:400 & gene_idx %in% 1:400) |> 
  dplyr::select(grna, gene, estimate) |>
  tidyr::pivot_wider(names_from = gene, values_from = estimate) |>
  tibble::column_to_rownames(var='grna')
estimates_matrix[is.na(estimates_matrix) ] = 0 # fill NA (from se) to 0 (e.g. not significant)

matscaled = as.matrix(scale(estimates_matrix))
matscaled[is.nan(matscaled)] = 0


# choose some ordering for the rows and columns (mainly for fixed matrix s.t. the blocks exist)
# ORDER_METHOD = 'hierarchichal' # hier cl on margins separately
# ORDER_METHOD = 'spectralsvd'   # CUSTOM SPECTRAL CL ON SVD # move to beginning of script
# ORDER_METHOD = 'biclust'   # something from biclust package (doesn't work)
if(ORDER_METHOD == 'hierarchichal') {
  # HIERARHICHAL CLUSTERING MARGINAL
  # order each margin separately w hclust
  row_order = hclust(dist(matscaled))$order
  column_order = hclust(dist(t(matscaled)))$order
  grna_index = data.frame(grna = rownames(matscaled)[row_order],
                          grna_idx  = 1:nrow(matscaled))
  gene_index = data.frame(gene = colnames(matscaled)[column_order],
                          gene_idx  = 1:ncol(matscaled))
} else if(ORDER_METHOD == 'spectralsvd') {
  # CUSTOM SPECTRAL CL ON SVD
  # I will just try basically spectral decomp and then order by eigenvalues???
  
  # K-Means Clustering on left and right svd eigenvecs
  svdres = svd(matscaled)
  plot(svdres$d)
  
  k = min(5, NUM_GENES-1, NUM_GRNAS-1) # when testing smaller sets, k cannot be too large
  row_clusters = kmeans(svdres$u[, 1:k], centers = k)$cluster
  col_clusters = kmeans(svdres$v[, 1:k], centers = k)$cluster
  grna_index = data.frame(grna = row.names(matscaled),
                          cl   = row_clusters) |>
    arrange(cl) |>
    mutate(grna_idx = 1:n()) # |> select(grna, grna_idx)
  gene_index = data.frame(gene = colnames(matscaled),
                          cl   = col_clusters) |>
    arrange(cl) |>
    mutate(gene_idx = 1:n()) # |> select(gene, gene_idx)
} else if(ORDER_METHOD == 'biclust') {
  # NONE OF THESE WORK (package: biclust)
  # BICLUSTER: spectral doesn't work, clusters not distinct enough
  biclres = biclust::spectral(matscaled,
                              normalization="log", numberOfEigenvalues=6, 
                              minr=2, minc=2, withinVar=1, n_clusters = NULL, n_best = 3)
  # Warning message:
  #   In biclust::spectral(matscaled, normalization = "log", numberOfEigenvalues = 6,  :
  #                          No biclusters found
  # BICLUSTER: Plaid - doesn't work, clusters not distinct enough
  biclres = biclust::biclust(matscaled, method=biclust::BCPlaid(), verbose=FALSE)
  
  # BICLUSTER: Quest - doesn't work
  biclres = biclust::biclust(matscaled, method=biclust::BCQuestmet(), quant=0.25, vari=1, ns=10, nd=10, sd=5, 
                             alpha=0.05, number=100)
  
  row_order = c(); column_order = c()
  for(cl in 1:biclres@Number) { # for each cluster, order grnas and genes
    row_order    = c(row_order, which(biclres@RowxNumber[, cl]))
    column_order = c(column_order, which(biclres@NumberxCol[cl,]))
  }
  row_order = c(row_order, setdiff(1:nrow(matscaled), row_order)) # not all rows/cols were assigned cl, add in rest
  column_order = c(column_order, setdiff(1:ncol(matscaled), column_order))
  
  # BICLUSTER:  CC - doesn't work, only finds 1
  biclres = biclust::biclust(matscaled, method=biclust::BCCC(), delta = 1.0, alpha=1.5, number=5)
  biclres@Number
}




# --------------------------------------------------------------------------------------------------------------------------
# 3. construct estimate matrices for all, train, test for chosen grnas and genes
# need grna_index, gene_index (chosen_genes and chosen_grnas should be in this)
est_matrices = list()
se_matrices = list()
for(split in c('all', 'train', 'test')) {
  # tall dataframe (positive tests are within the discoveries (e.g. so we only filter on discoveries to avoid repeitions))
  effects_df_split =  effects[[split]] |> dplyr::filter((test == 'discovery') & (gene %in% chosen_genes) & (grna %in% chosen_grnas)) 
  # estimates
  # matrix form: grna x gene
  matrices_estimate_split = effects_df_split  |> 
    dplyr::select(grna, gene, estimate) |>
    tidyr::pivot_wider(names_from = gene, values_from = estimate) |>
    tibble::column_to_rownames(var='grna') |> as.matrix()
  # order based on grna_index and gene_index
  est_matrices[[split]] = matrices_estimate_split[grna_index |> arrange(grna_idx) |> pull(grna), 
                                                  gene_index |> arrange(gene_idx) |> pull(gene)]
  # se
  # matrix form: grna x gene
  matrices_se_split = effects_df_split  |> 
    dplyr::select(grna, gene, se) |>
    tidyr::pivot_wider(names_from = gene, values_from = se) |>
    tibble::column_to_rownames(var='grna') |> as.matrix()
  # order based on grna_index and gene_index
  se_matrices[[split]] = matrices_se_split[ grna_index |> arrange(grna_idx) |> pull(grna), 
                                            gene_index |> arrange(gene_idx) |> pull(gene)]
  
  rm(effects_df_split, matrices_estimate_split, matrices_se_split); gc()
}

# save the chosen grna_index and gene_index (ordered collection of grna and gene)
gene_index = gene_index |> select(gene, gene_idx) # only keep these cols (e.g. don't need cluster col if added)
grna_index = grna_index |> select(grna, grna_idx)
write.csv(x = gene_index |> select(gene, gene_idx), file=sprintf('%s/gene_index.csv', replogle_save_path), row.names=FALSE)
write.csv(x = grna_index |> select(grna, grna_idx), file=sprintf('%s/grna_index.csv', replogle_save_path), row.names=FALSE)



saveRDS(list(est_matrices=est_matrices, se_matrices=se_matrices,
             DESCRIPTION = 'Matrices of Estimates (SCEPTRE log_2_fc) and SE (standard errors): using the chosen top and sorted grnas (rows) and genes (cols)'), 
        sprintf('%s/estse_matrices.rds', replogle_save_path))


# ======================================================================================================================================
#                                                                               ========================================================
#       INITIAL SCEPTRE RESULT PLOTS:                                           ========================================================
#             - sceptre estimated effects + approx matrices                     ========================================================
# ======================================================================================================================================

# ------------------------------------------------------------------------------------------------------------------------------------------
# First, some plots showing SCEPTRE Results
print(sprintf("[%s]     - some plots", Sys.time()))
dir.create(plot_path, showWarnings=FALSE)
dir.create(sprintf('%s/initialSCEPTRE/', plot_path), showWarnings=FALSE)








# color breaks for plotting ------------------------------------------------
color_limits = c(-2, 2)
color_breaks = sort(union(color_limits, seq(from = round(color_limits[1]), to = round(color_limits[2]))))
color_breaks_label = color_breaks
color_breaks_label[which.min(color_breaks)] = sprintf('<%.1f', min(color_breaks))
color_breaks_label[which.max(color_breaks)] = sprintf('>%.1f', max(color_breaks))



# for all, train, test split
for(data_split in c('train', 'test', 'all')) {
  # add original estimates ---------------------------------------------------
  # (all the loaded dfs should have these) (2x for plotting)
  # plot_df = effects_df |> select(grna, gene, estimate)
  # plot_df = effects[[data_split]] |> filter(test == 'discovery')
  plot_df = effects[[data_split]] |> filter((test == 'discovery') & (gene %in% chosen_genes) & (grna %in% chosen_grnas)) 
  
  # add idx's- (from gene or grna name --> number on the axis)
  plot_df = merge(gene_index, plot_df, by = 'gene', all.y = TRUE, all.x = FALSE)
  plot_df = merge(grna_index, plot_df, by = 'grna', all.y = TRUE, all.x = FALSE)


  # plot_df$rank = factor(plot_df$rank, levels = c('unshrunk', 'shrunk'))  



  # --------------------------------------------------------------------------------------------------------------------------
  # heatmap of diff values side by side
  plot_df_heatmap =
    rbind(plot_df |> select(grna, gene,    estimate)  |> dplyr::rename(val =     estimate) |> dplyr::mutate(type = 'log_2_fold_change'), 
          plot_df |> select(grna, gene, fold_change)  |> dplyr::rename(val =  fold_change) |> dplyr::mutate(type = 'fold_change'), 
          plot_df |> select(grna, gene, significant)  |> dplyr::rename(val =  significant) |> dplyr::mutate(type = 'significant'),
          plot_df |> select(grna, gene,       tstat)  |> dplyr::rename(val =        tstat) |> dplyr::mutate(type = 'tstat'),
          plot_df |> mutate(isna = is.na(estimate))   |> 
                     select(grna, gene,       isna)   |> dplyr::rename(val =         isna) |> dplyr::mutate(type = 'isna')
          )



  

  # add idx's- (from gene or grna name --> number on the axis)
  plot_df_heatmap = merge(gene_index, plot_df_heatmap, by = 'gene')
  plot_df_heatmap = merge(grna_index, plot_df_heatmap, by = 'grna')


  plot_df_heatmap$type = factor(plot_df_heatmap$type, 
                                levels = c('log_2_fold_change', 'fold_change', 'tstat', 'significant', 'isna'),
                                labels = c('Log2 Fold Change',  'Fold Change', 't-stat', 'significant', 'missing'))  



  # plot
  p_heatmap = ggplot(plot_df_heatmap) +
    geom_raster(aes(x = grna_idx, y = gene_idx, fill = val)) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x='grna', y = 'gene', fill = NULL, title = 'Estimates') +
    scale_fill_gradient2(limits = color_limits, # set color limits
                         oob=scales::squish, # if outside lims, set to limits
                         midpoint = 0,
                         high = myRed, low = myBlue, mid = 'white',
                         # low  = brewer.pal(n = 9, name = "RdBu")[9],
                         # high = brewer.pal(n = 9, name = "RdBu")[1],
                         breaks = color_breaks,
                         labels = color_breaks_label) +
    facet_grid(#rows = vars(approxmethod),
      cols = vars(type)) +
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
  p_heatmap
  
  ggsave(plot = p_heatmap, filename = sprintf('%s/initialSCEPTRE/%s.pdf', plot_path, data_split), height = 7, width = 7)

  rm(plot_df, plot_df_heatmap, p_heatmap)
}


# ======================================================================================================================================
# --------------------------------------------------------------------------------------------------------------------------------------
#       PERFORM MATRIX APPROXIMATIONS:                                             =====================================================
#       these are matrices to shrink towards                                       =====================================================
#       approx the log_2_fold matrices using the following methods ┌┌─┬┐─────────┐ =====================================================
#         - low rank (SVD)                                         │├┬┴┤         │ =====================================================
#         - sparse SVD                                             │└┴─┘╔═╗      │ =====================================================
#         - softImpute                                             │    ╚═╝XXXXXX│ =====================================================
#         - spectralbiclust                                        │==     XXXXXX│ ===================================================== 
#         - spectralbiclust_threshold                              │==     XXXXXX│ ===================================================== 
#         - zeros                                                  └─────────────┘ ===================================================== 
#         - average                                                                =====================================================
# --------------------------------------------------------------------------------------------------------------------------------------
# ======================================================================================================================================


# for each split type (no split vs split)
# for each matrix approximation method
# split_type = 'nosamplesplit' # vs 'samplesplit' 

my_methodParams = list('sparseSVD' = list(type = 'standard', 
                                          sumabs = .35,  # should be between 0-1, 
                                          # sumabsu = 4, sumabsv = 4, # between 1 and sqrt(#col or #rows)
                                          niter = 100,
                                          trace = FALSE))


dir.create(sprintf('%s/matapprox/',      plot_path), showWarnings=FALSE) # replogle/EBCI/matapprox/
dir.create(sprintf('%s/matapprox/all',   plot_path), showWarnings=FALSE) # replogle/EBCI/matapprox/all
dir.create(sprintf('%s/matapprox/train', plot_path), showWarnings=FALSE) # replogle/EBCI/matapprox/train

# plot settings- consistent across approxmatrices plots- based on 'all' estimates
# grna_gene_ordering = plot_set_rowcol_order(mat = all, save_folder=sprintf('%s/mat/all/', plot_path)) # rerun diff for this
# # save the chosen grna_index and gene_index (ordered collection of grna and gene)
# gene_index = read.csv(file=sprintf('%s/gene_index.csv', replogle_save_path), row.names=FALSE)
# grna_index = read.csv(file=sprintf('%s/grna_index.csv', replogle_save_path), row.names=FALSE)
grna_gene_ordering = list(gene_idx=gene_index, grna_index=grna_index)

# color breaks for plotting
# color_limits = c(-2, 2) # already defined previously, don't change here

approxmatrices = list() # list of matrix approximations
for(data_split in c('all', 'train')) {
  print(sprintf("[%s] Approximating on %s split", format( Sys.time(), format = "%F %r"), data_split))
  # dir.create(sprintf('%s/%s', plot_path, data_split),  showWarnings = FALSE) 
  approxmatrices[[data_split]] = list()
  mat_to_approx = est_matrices[[data_split]]   
  mat_to_approx[is.na(mat_to_approx)] = 0 # fill NA as 0

  for(matapprox_method in matapprox_methods) {
    dir.create(sprintf('%s/matapprox/%s/%s/', plot_path, data_split, matapprox_method),  showWarnings = FALSE) 

    # use function approx_matrix from utils/matrix_shrinkage.r 
    approxmatrices[[data_split]][[matapprox_method]] = 
      tryCatch(expr = {

        approx_matrix(mat                        = mat_to_approx,
                      method                     = matapprox_method,
                      ranks                      = ranks, 
                      save_folder                = sprintf('%s/matapprox/%s/%s/', plot_path, data_split, matapprox_method), # e.g. replogle/matapprox/train/SVD/
                      save_individual_rank_plots = TRUE,
                      color_limits               = color_limits,
                      methodParams               = if(matapprox_method %in% names(my_methodParams)) my_methodParams[[matapprox_method]] else NULL,
                      grna_gene_ordering         = grna_gene_ordering
                      )    

      },
        error = function(e) {print(sprintf("Errored at: %s, %s", data_split, matapprox_method)); NULL}
      )
  }
}


saveRDS(approxmatrices, sprintf('%s/approxmatrices.rds', replogle_save_path)) # save the matrix approximations (on log_2_fold_change on train and all splits)





# ======================================================================================================================================
# --------------------------------------------------------------------------------------------------------------------------------------
#                                                                               ========================================================
#       PERFORM EBCI SHRINKAGE:                                                 ========================================================
#                                                                               ========================================================
#                                                                               ========================================================
# --------------------------------------------------------------------------------------------------------------------------------------
# ======================================================================================================================================
dir.create(sprintf('%s/shrinkage/', plot_path), showWarnings=FALSE) # plots/replogle/EBCI/shrinkage/


shrink_results = list()
for(split_type in c('nosamplesplit', 'samplesplit')) {
  dir.create(sprintf('%s/shrinkage/%s/', plot_path, split_type), showWarnings=FALSE)

  shrink_results[[split_type]] = list()

  if(split_type == 'nosamplesplit') { # shrink ALL on ALL
    mat_est_to_shrink =   est_matrices[['all']]
    mat_se_to_shrink  =    se_matrices[['all']]
    mats_to_shrink_to = approxmatrices[['all']] # list of matrices
  } else if(split_type == 'samplesplit') { # shrink TEST on TRAIN
    mat_est_to_shrink =   est_matrices[['test']]
    mat_se_to_shrink  =    se_matrices[['test']] 
    mats_to_shrink_to = approxmatrices[['train']] # list of matrices
  } else {
    print('(replogleShrinkage: shrinkage) bad split_type')
  }
  
  
  for(matapprox_method in matapprox_methods) {
    dir.create(sprintf('%s/shrinkage/%s/%s/', plot_path, split_type, matapprox_method), showWarnings=FALSE) # plots/

    #' temporary function that shrinks the current matrices in environment to shrink to 
    #' varying mats_to_shrink_to
    temp_shrink_matrix <- function(mat_to_shrink_to, save_filename) {
       
        shrinkResult = tryCatch(expr = {
            # shrink_matrix( unshrunk_mat    = mat_est_to_shrink,
            #                shrinkpoint_mat = mat_to_shrink_to,
            #                se_mat          = mat_se_to_shrink,
            #                weight_mat      = (1/mat_se_to_shrink)**4,
            #                ALPHA           = ALPHA, 
            #                return_ebci_obj = TRUE) # ~  mins
          
          shrink_matrix( unshrunk_mat    = mat_est_to_shrink[1:10, 1:50],
                         shrinkpoint_mat = mat_to_shrink_to[1:10, 1:50],
                         se_mat          = mat_se_to_shrink[1:10, 1:50],
                         weight_mat      = (1/mat_se_to_shrink[1:10, 1:50])**4,
                         ALPHA           = ALPHA, 
                         return_ebci_obj = TRUE) # ~  mins

          }, 
          error = function(e){print(sprintf("Errored at: %s, %s", split_type, matapprox_method)); NULL}
        )
          
        # plot(shrinkResult$ebci_res$shrinkage_point, shrinkResult$ebci_res$shrunk_value)
        # write.csv(x = shrinkResult$ebci_res, file = sprintf("../saves/replogle/shrinkage/replogle_shrink_sparseSVD03.csv")) 
        write.csv(x = shrinkResult$ebci_res, file = save_filename) 
        return(shrinkResult)
    }
    if(matapprox_methods_hasranks[[matapprox_method]]) {
      shrink_results[[matapprox_method]] = list()
      for(r in ranks) {
        dir.create(sprintf('%s/shrinkage/%s/%s/rank=%02.f/', plot_path, split_type, matapprox_method, r), showWarnings=TRUE) # rank folder
        shrinkResult = temp_shrink_matrix(
          mat_to_shrink_to = mats_to_shrink_to[[matapprox_method]][[r]]$approxmatrices,
          save_filename = sprintf('%s/shrinkage/%s/%s/rank=%02.f/ebci_shrinkage_df.csv', plot_path, split_type, matapprox_method, r)
        )
        shrink_results[[matapprox_method]][[r]] = shrinkResult$ebci_res

      } 
    } else {
      shrinkResult = temp_shrink_matrix(
        mat_to_shrink_to = mats_to_shrink_to[[matapprox_method]]$approxmatrices,
        # save_filename = sprintf('%s/shrinkage/a.csv', plot_path)
        save_filename = sprintf('%s/shrinkage/%s/%s/ebci_shrinkage_df.csv', plot_path, split_type, matapprox_method)
      )
      shrink_results[[matapprox_method]] = shrinkResult$ebci_res
    }
  }
}
   




saveRDS(shrink_results, sprintf('%s/EBCI_shrinkage_results.rds', replogle_save_path)) # save the shrinkage results 







# ======================================================================================================================================
# --------------------------------------------------------------------------------------------------------------------------------------
#                            +   │─┐                                            ========================================================
#       SHRINKAGE PLOTS:     |─┐ │ └┐                                           ========================================================
#                            | ──┘  ────                                        ========================================================
#                            +---------+                                        ========================================================
# --------------------------------------------------------------------------------------------------------------------------------------
# ======================================================================================================================================


# dir.create(sprintf('%s/shrinkage/', plot_path), showWarnings=FALSE) # plots/replogle/EBCI/shrinkage/


for(split_type in c('nosamplesplit', 'samplesplit')) {  
  for(matapprox_method in matapprox_methods) {

    #' temporary function that plots some summary results of the shrinkage
    temp_plot_shrinkage <- function(cur_plot_folder) {
       
        tryCatch(expr = {
           
           # sparse SVD, rank=3
          # plot_folder = '../plots/replogle/shrink/spSVD03/'
          # shrink_df = read.csv("../saves/replogle/shrinkage/replogle_shrink_sparseSVD03.csv")
          shrink_df = read.csv(sprintf("%s/ebci_shrinkage_df.csv", cur_plot_folder))
          dir.create(sprintf('%s/points/', cur_plot_folder)); dir.create(sprintf('%s/heatmaps/', cur_plot_folder))
          plot_shrink_results(shrink_df=shrink_df, plot_folder=cur_plot_folder, order_rowscols=T, grna_index=grna_index, gene_index=gene_index, unshrunk_ALPHA=ALPHA)
          }, 
          error = function(e){print(sprintf("Errored at: %s", cur_plot_folder))}
        )
        return(NULL)
    }
    if(matapprox_methods_hasranks[[matapprox_method]]) {
      for(r in ranks) {
        temp_plot_shrinkage(cur_plot_folder = sprintf('%s/shrinkage/%s/%s/rank=%02.f/', plot_path, split_type, matapprox_method, r))
      }
    } else {
        temp_plot_shrinkage(cur_plot_folder = sprintf('%s/shrinkage/%s/%s/', plot_path, split_type, matapprox_method))
    }
  }
}









# ======================================================================================================================================
# --------------------------------------------------------------------------------------------------------------------------------------
#                  ==^==                                                             ===================================================
#       OLD:       |[[[|                                                             ===================================================
#                  |[[[|                                                             ===================================================
#                  '---'                                                             ===================================================
# --------------------------------------------------------------------------------------------------------------------------------------
# ====================================================================================================================================== 



if(F) {
  



# --------------------------------------------------------------------------------------------------------------------------------------
# perf sparse SVD on log_2_fc----------------------------------------------------------------------------------------------------------- 
# --------------------------------------------------------------------------------------------------------------------------------------
ranks = c(1, 3, 10, 20, 30)
approxmatrices_sparseSVD = 
  approx_matrix(mat = estimates_matrix, 
                method = 'sparseSVD', 
                ranks = ranks,
                methodParams=list(type = 'standard', 
                                  sumabs = .3,  # should be between 0-1, 
                                  # sumabsu = 4, sumabsv = 4, # between 1 and sqrt(#col or #rows)
                                  niter = 100,
                                  trace = FALSE)) 


approxmatrix_df = estimates_matrix |> as.matrix() |>
  reshape2::melt() |> 
  dplyr::rename(grna = Var1, gene = Var2, val = value) |>
  dplyr::mutate(type = 'original')
for(r in ranks) {
  approxmatrix_df = rbind(approxmatrix_df, 
                          approxmatrices_sparseSVD$approxmatrices[[r]] |> 
                            reshape2::melt() |> 
                            dplyr::rename(grna = Var1, gene = Var2, val = value) |>
                            dplyr::mutate(type = sprintf('rank=%02.f', r)))
}

approxmatrix_df$type = factor(approxmatrix_df$type, levels = c('original', sort(sprintf('rank=%02.f', ranks), decreasing = TRUE)))

# plot
# add idx's- (from gene or grna name --> number on the axis)
plot_df_heatmap = merge(gene_index, approxmatrix_df, by = 'gene')
plot_df_heatmap = merge(grna_index, plot_df_heatmap, by = 'grna')
p_heatmap = ggplot(plot_df_heatmap) +
  geom_raster(aes(x = grna_idx, y = gene_idx, fill = val)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x='grna', y = 'gene', fill = NULL, title = 'Estimates') +
  scale_fill_gradient2(limits = color_limits, # set color limits
                       oob=scales::squish, # if outside lims, set to limits
                       midpoint = 0,
                       high = myRed, low = myBlue, mid = 'white',
                       # low  = brewer.pal(n = 9, name = "RdBu")[9],
                       # high = brewer.pal(n = 9, name = "RdBu")[1],
                       breaks = color_breaks,
                       labels = color_breaks_label) +
  facet_grid(#rows = vars(approxmethod),
    cols = vars(type)) +
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
p_heatmap

ggsave(plot = p_heatmap, filename = sprintf('%s/sceptre_estimates_sparseSVD.pdf', plot_path), height = 7, width = 10)




# --------------------------------------------------------------------------------------------------------------------------
# perf lowrank approx on log_2_fc 
ranks = c(1, 3, 10, 20, 30)
estimates_matrix[is.na(estimates_matrix)] = 0 # fill with 0?
approxmatrices_sparseSVD = 
  approx_matrix(mat = estimates_matrix, 
                method = 'lowrank', 
                ranks = ranks) 

approxmatrix_df = estimates_matrix |> as.matrix() |>
  reshape2::melt() |> 
  dplyr::rename(grna = Var1, gene = Var2, val = value) |>
  dplyr::mutate(type = 'original')
for(r in ranks) {
  approxmatrix_df = rbind(approxmatrix_df, 
                          approxmatrices_sparseSVD$approxmatrices[[r]] |> 
                            reshape2::melt() |> 
                            dplyr::rename(grna = Var1, gene = Var2, val = value) |>
                            dplyr::mutate(type = sprintf('rank=%02.f', r)))
}

approxmatrix_df$type = factor(approxmatrix_df$type, levels = c('original', sort(sprintf('rank=%02.f', ranks), decreasing = TRUE)))

# plot
# add idx's- (from gene or grna name --> number on the axis)
plot_df_heatmap = merge(gene_index, approxmatrix_df, by = 'gene')
plot_df_heatmap = merge(grna_index, plot_df_heatmap, by = 'grna')
p_heatmap = ggplot(plot_df_heatmap) +
  geom_raster(aes(x = grna_idx, y = gene_idx, fill = val)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x='grna', y = 'gene', fill = NULL, title = 'Estimates') +
  scale_fill_gradient2(limits = color_limits, # set color limits
                       oob=scales::squish, # if outside lims, set to limits
                       midpoint = 0,
                       high = myRed, low = myBlue, mid = 'white',
                       # low  = brewer.pal(n = 9, name = "RdBu")[9],
                       # high = brewer.pal(n = 9, name = "RdBu")[1],
                       breaks = color_breaks,
                       labels = color_breaks_label) +
  facet_grid(#rows = vars(approxmethod),
    cols = vars(type)) +
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
p_heatmap

ggsave(plot = p_heatmap, filename = sprintf('%s/sceptre_estimates_lowrank.pdf', plot_path), height = 7, width = 10)




# ======================================================================================================================================
#                                                                               ========================================================
#       SHRINKAGE:                                                              ========================================================
#                                                                               ========================================================
# ======================================================================================================================================


dir.create('../saves/replogle/shrinkage/')

# --------------------------------------------------------------------------------------------------------------------------
# Sparse SVD
approxmatrices_sparseSVD = 
  approx_matrix(mat = est_matrices[['train']], 
                method = 'sparseSVD', 
                ranks = c(3, 10),
                methodParams=list(type = 'standard', 
                                  sumabs = .35,  # should be between 0-1, 
                                  # sumabsu = 4, sumabsv = 4, # between 1 and sqrt(#col or #rows)
                                  niter = 100,
                                  trace = FALSE)) 

shrink2Grnas = 
  shrink_matrix(   unshrunk_mat = est_matrices[['test']][c(1,2), ],
                   shrinkpoint_mat = approxmatrices_sparseSVD$approxmatrices[[3]][c(1,2), ],
                   se_mat = se_matrices[['test']][c(1,2), ],
                   weight_mat = (1/se_matrices[['test']][c(1,2), ])**4,
                   ALPHA = ALPHA, 
                   return_ebci_obj = TRUE) # ~  mins

# plot(shrink2Grnas$ebci_res$shrinkage_point, shrink2Grnas$ebci_res$shrunk_value)


write.csv(x = shrink2Grnas$ebci_res, file = sprintf("../saves/replogle/shrinkage/replogle_shrink_2grna.csv"))


shrinkSparse = 
  shrink_matrix(   unshrunk_mat = est_matrices[['test']],
                   shrinkpoint_mat = approxmatrices_sparseSVD$approxmatrices[[3]],
                   se_mat = se_matrices[['test']],
                   weight_mat = (1/se_matrices[['test']])**4,
                   ALPHA = ALPHA, 
                   return_ebci_obj = TRUE) # ~  mins
# plot(shrinkSparse$ebci_res$shrinkage_point, shrinkSparse$ebci_res$shrunk_value)


write.csv(x = shrinkSparse$ebci_res, file = sprintf("../saves/replogle/shrinkage/replogle_shrink_sparseSVD03.csv"))


shrinkSparse = 
  shrink_matrix(   unshrunk_mat = est_matrices[['test']],
                   shrinkpoint_mat = approxmatrices_sparseSVD$approxmatrices[[10]],
                   se_mat = se_matrices[['test']],
                   weight_mat = (1/se_matrices[['test']])**4,
                   ALPHA = ALPHA, 
                   return_ebci_obj = TRUE) # ~  mins
# plot(shrinkSparse$ebci_res$shrinkage_point, shrinkSparse$ebci_res$shrunk_value)


write.csv(x = shrinkSparse$ebci_res, file = sprintf("../saves/replogle/shrinkage/replogle_shrink_sparseSVD10.csv"))



# --------------------------------------------------------------------------------------------------------------------------
# Low Rank
est_mat_train_fillna = est_matrices[['train']]
est_mat_train_fillna[is.na(est_mat_train_fillna)] = 0
approxmatrices_lowrank = 
  approx_matrix(mat = est_mat_train_fillna, 
                method = 'lowrank', 
                ranks = c(3, 10)) 
rm(est_mat_train_fillna); gc()


shrinkLowrank = 
  shrink_matrix(   unshrunk_mat = est_matrices[['test']],
                   shrinkpoint_mat = approxmatrices_lowrank$approxmatrices[[3]],
                   se_mat = se_matrices[['test']],
                   weight_mat = (1/se_matrices[['test']])**4,
                   ALPHA = ALPHA, 
                   return_ebci_obj = TRUE) # ~  mins
# plot(shrinkLowrank$ebci_res$shrinkage_point, shrinkLowrank$ebci_res$shrunk_value)


write.csv(x = shrinkLowrank$ebci_res, file = sprintf("../saves/replogle/shrinkage/replogle_shrink_lowrank03.csv"))




shrinkLowrank = 
  shrink_matrix(   unshrunk_mat = est_matrices[['test']],
                   shrinkpoint_mat = approxmatrices_lowrank$approxmatrices[[10]],
                   se_mat = se_matrices[['test']],
                   weight_mat = (1/se_matrices[['test']])**4,
                   ALPHA = ALPHA, 
                   return_ebci_obj = TRUE) # ~  mins
# plot(shrinkLowrank$ebci_res$shrinkage_point, shrinkLowrank$ebci_res$shrunk_value)


write.csv(x = shrinkLowrank$ebci_res, file = sprintf("../saves/replogle/shrinkage/replogle_shrink_lowrank10.csv"))

# ======================================================================================================================================
#                                                                               ========================================================
#       PLOTS:                                                                  ========================================================
#                                                                               ========================================================
# ======================================================================================================================================





# # load gene_index and grna_index df's to save time
# gene_index = read.csv(sprintf('%s/gene_index.csv', replogle_save_path))
# grna_index = read.csv(sprintf('%s/grna_index.csv', replogle_save_path))


# TODO: fix... somehow the  ordered gene and grna sets are not the same as the ones w estimates???
# test
shrink_df = read.csv("../saves/replogle/shrinkage/replogle_shrink_sparseSVD03.csv")
all(grna_index$grna %in% shrink_df$grna) # this should be true... but somehow it is not??
# reorder for those grna/gene present....

gene_index = gene_index |> filter(gene %in% shrink_df$gene) |> arrange(gene_idx) |> mutate(gene_idx = 1:n())
grna_index = grna_index |> filter(grna %in% shrink_df$grna) |> arrange(grna_idx) |> mutate(grna_idx = 1:n())





# start


# sparse SVD, rank=3
plot_folder = '../plots/replogle/shrink/spSVD03/'
shrink_df = read.csv("../saves/replogle/shrinkage/replogle_shrink_sparseSVD03.csv")

dir.create(sprintf('%s/points/', plot_folder), recursive = T); dir.create(sprintf('%s/heatmaps/', plot_folder))
plot_shrink_results(shrink_df=shrink_df, plot_folder=plot_folder, order_rowscols=T, grna_index=grna_index, gene_index=gene_index, unshrunk_ALPHA=ALPHA)

# sparse SVD, rank=10
plot_folder = '../plots/replogle/shrink/spSVD10/'
shrink_df = read.csv("../saves/replogle/shrinkage/replogle_shrink_sparseSVD10.csv")

dir.create(sprintf('%s/points/', plot_folder), recursive = T); dir.create(sprintf('%s/heatmaps/', plot_folder))
plot_shrink_results(shrink_df=shrink_df, plot_folder=plot_folder, order_rowscols=T, grna_index=grna_index, gene_index=gene_index, unshrunk_ALPHA=ALPHA)

# lowrank, rank=3
plot_folder = '../plots/replogle/shrink/lowrank03/'
shrink_df = read.csv("../saves/replogle/shrinkage/replogle_shrink_lowrank03.csv")

dir.create(sprintf('%s/points/', plot_folder), recursive = T); dir.create(sprintf('%s/heatmaps/', plot_folder))
plot_shrink_results(shrink_df=shrink_df, plot_folder=plot_folder, order_rowscols=T, grna_index=grna_index, gene_index=gene_index, unshrunk_ALPHA=ALPHA)

# lowrank, rank=10
plot_folder = '../plots/replogle/shrink/lowrank10/'
shrink_df = read.csv("../saves/replogle/shrinkage/replogle_shrink_lowrank10.csv")

dir.create(sprintf('%s/points/', plot_folder), recursive = T); dir.create(sprintf('%s/heatmaps/', plot_folder))
plot_shrink_results(shrink_df=shrink_df, plot_folder=plot_folder, order_rowscols=T, grna_index=grna_index, gene_index=gene_index, unshrunk_ALPHA=ALPHA)







}



