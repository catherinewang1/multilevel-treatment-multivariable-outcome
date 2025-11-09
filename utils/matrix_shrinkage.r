# Functions for performing matrix approximation + ebci shrinkage
# Analysis steps:
#   0. Load and prep data
#   1. Get a lower structure approximation of effects matrix
#   2. Shrink towards lower structure
#   3. Make some plots


myRed  = colorRampPalette(brewer.pal(n = 7, name = "RdBu"))(7)[1]
myBlue = colorRampPalette(brewer.pal(n = 7, name = "RdBu"))(7)[7]


#' From an effects dataframe where rows = test, cols = effects, other vals...
#' make matrices that are #perturbations x #genes shaped 
#' e.g. effects_matrix can be used in image(), heatmap(), ...
#' @param effects_df (dataframe)
#' @param save_plots_filepath (string) or NULL for no plots
make_matrices <- function(effects_df, save_plots_filepath=NULL) {
  # save_plots_filepath = '../plots/matrix/'
  # Some tests missing bc failed quality check (from sceptre)
  # fill estimate = 0, and se = 999 (for now)
  foldchange_matrix = effects_df |> # filter(grna_idx %in% 1:400 & gene_idx %in% 1:400) |> 
    dplyr::select(grna, gene, fold_change) |>
    tidyr::pivot_wider(names_from = gene, values_from = fold_change) |>
    tibble::column_to_rownames(var='grna')
  foldchange_matrix[is.na(foldchange_matrix)] = 0
  
  # fill estimate = 0, and se = 999 (for now)
  estimates_matrix = effects_df |> # filter(grna_idx %in% 1:400 & gene_idx %in% 1:400) |> 
    dplyr::select(grna, gene, estimate) |>
    tidyr::pivot_wider(names_from = gene, values_from = estimate) |>
    tibble::column_to_rownames(var='grna')
  estimates_matrix[is.na(estimates_matrix)] = 0
  
  se_matrix = effects_df |> # filter(grna_idx %in% 1:400 & gene_idx %in% 1:400) |> 
    dplyr::select(grna, gene, se) |>
    tidyr::pivot_wider(names_from = gene, values_from = se) |>
    tibble::column_to_rownames(var='grna')
  se_matrix[is.na(se_matrix)] = 999
  
  tstat_matrix = effects_df |> # filter(grna_idx %in% 1:400 & gene_idx %in% 1:400) |> 
    dplyr::select(grna, gene, tstat) |>
    tidyr::pivot_wider(names_from = gene, values_from = tstat) |>
    tibble::column_to_rownames(var='grna')
  tstat_matrix[is.na(tstat_matrix) ] = 0 # fill NA (from se) to 0 (e.g. not significant)
  
  significant_matrix = effects_df |>
    dplyr::select(grna, gene, significant) |>
    tidyr::pivot_wider(names_from = gene, values_from = significant) |>
    tibble::column_to_rownames(var='grna')
  significant_matrix[is.na(significant_matrix) ] = 0 # fill NA (from se) to 0 (e.g. not significant)
  
  
  # vis
  if(!is.null(save_plots_filepath)) {
    require(RColorBrewer)
    # Fold Change
    p = pheatmap::pheatmap(foldchange_matrix |> as.matrix(), treeheight_row = 0, treeheight_col = 0, 
                           color = colorRampPalette((brewer.pal(n = 7, name = "YlOrRd")))(100),
                           show_rownames = FALSE, show_colnames = FALSE)
    ggsave(plot = p, filename = sprintf('%s/foldchange.pdf', save_plots_filepath))
    
    # Estimates- SOME COLORS ARE CUT OFF TO MIN/MAX! range(estimates_matrix) hist(estimates_matrix |> as.matrix() |> as.vector())
    p = pheatmap::pheatmap(estimates_matrix |> as.matrix(), treeheight_row = 0, treeheight_col = 0, 
                           color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100), 
                           breaks = seq(from = -2, to = 2, length.out = 101),
                           show_rownames = FALSE, show_colnames = FALSE)
    ggsave(plot = p, filename = sprintf('%s/estimates.pdf', save_plots_filepath))
    
    
    # SE range(tstat_matrix)
    p = pheatmap::pheatmap(se_matrix |> as.matrix(), treeheight_row = 0, treeheight_col = 0, 
                           color = colorRampPalette((brewer.pal(n = 7, name = "YlOrRd")))(100), 
                           show_rownames = FALSE, show_colnames = FALSE)
    ggsave(plot = p, filename = sprintf('%s/se.pdf', save_plots_filepath))
    
    
    # T statistic- SOME COLORS ARE CUT OFF TO MIN/MAX! range(tstat_matrix) hist(tstat_matrix |> as.matrix() |> as.vector())
    p = pheatmap::pheatmap(tstat_matrix |> as.matrix(), treeheight_row = 0, treeheight_col = 0, 
                           color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
                           breaks = seq(from = -6, to = 6, length.out = 101),
                           show_rownames = FALSE, show_colnames = FALSE)
    ggsave(plot = p, filename = sprintf('%s/tstat.pdf', save_plots_filepath))
    
    # Significance (0/1)
    p = pheatmap::pheatmap(significant_matrix |> as.matrix(), treeheight_row = 0, treeheight_col = 0, 
                           color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100), 
                           show_rownames = FALSE, show_colnames = FALSE)
    ggsave(plot = p, filename = sprintf('%s/significance.pdf', save_plots_filepath))
  }
  
  
  
  
  matrices = list(foldchange = foldchange_matrix,
                  estimates = estimates_matrix,
                  se = se_matrix,
                  tstat = tstat_matrix,
                  significance = significant_matrix)
  return(matrices)
}



#' Get an approximation of the given matrix (mat) with the given ranks (ranks)
#' Can be a low rank matrix (SVD) or a sparse low rank matrix (sparse SVD (PMA::PMD))
#' @param mat (matrix) matrix to approximate
#' @param method (character) method for approximating the matrix
#'         either 'lowrank' or 'sparseSVD' 
#' @param ranks (vector) of positive integers indicating the ranks
#' @param save_plots_filepath (string) or NULL for no plots
#' @param save_individual_rank_plots (boolean) if save_plots_filepath is specified, 
#'     should e also save plots for every rank individually? 
#'     if TRUE, saves plots like 'rank=10.pdf'
#'     if FALSE, does not save plots for every rank, and only saves singular
#'     values and all ranks together
#' @param color_limits (vector) of length 2, for the limits of the colors 
#'     default = c(-2, 2), try to use integers... it looks ugly ow
#' @param methodParams (list) of other parameters used in the approximation method
#' @example
#' lowrankmatrices = approx_lowrank(
#'        mat = matrices$estimates,
#'        ranks = c(1, 3, 5, 10, 15),
#'        save_plots_filepath = '../plots/matrix/lowrank/estimate/',
#'        color_limits = c(-2, 2)
#' )
approx_matrix <- function(mat, method, ranks, save_plots_filepath=NULL, save_individual_rank_plots=FALSE, color_limits = c(-2, 2), methodParams=NULL) {
  
  
  # mat = matrices$estimates|> as.matrix() 
  # mat |> as.matrix() |> as.vector() |> hist() # <- set limits to -2, 2
  # method = 'sparseSVD'
  # ranks = c(1, 3, 5, 10, 20)
  # save_plots_filepath = "../plots/matrix/sparseSVD/"
  # color_limits = c(-2, 2)
  # methodParams = list(type = 'standard', 
  #                     sumabs = .2, 
  #                     # sumabsu = 4, sumabsv = 4,
  #                     niter = 50,
  #                     trace = FALSE)
  
  
  
  # create some objects used for plotting later
  # - get a row and col order using clustering (marginally), for plotting later
  # - tall version of mat matrix
  if(!is.null(save_plots_filepath)) {
    matscaled = as.matrix(scale(mat))
    matscaled[is.nan(matscaled)] = 0
    row_order = hclust(dist(matscaled))$order
    column_order = hclust(dist(t(matscaled)))$order
    
    
    grna_index = data.frame(grna = rownames(mat)[row_order],
                            grna_idx  = 1:nrow(mat))
    gene_index = data.frame(gene = colnames(mat)[column_order],
                            gene_idx  = 1:ncol(mat))
    
    rm(matscaled, row_order, column_order)
    
    
    mat_tall = mat |> as.data.frame() |> tibble::rownames_to_column('grna') |> 
      tidyr::pivot_longer(cols = 2:(ncol(mat)+1), 
                          names_to = 'gene', values_to = 'value') |>
      mutate(group = 'original')
    
    
    mat_tall = merge(gene_index, mat_tall, by = 'gene')
    mat_tall = merge(grna_index, mat_tall, by = 'grna')
    
    plot_df_all = mat_tall
    
    # color breaks for plotting
    color_breaks = sort(union(color_limits, seq(from = round(color_limits[1]), to = round(color_limits[2]))))
    color_breaks_label = color_breaks
    color_breaks_label[which.min(color_breaks)] = sprintf('<%.1f', min(color_breaks))
    color_breaks_label[which.max(color_breaks)] = sprintf('>%.1f', max(color_breaks))
  }
  
  
  if(method == 'lowrank') {
    # low rank approximation ---------------------------
    svd_res = svd(x = mat, nu = max(ranks), nv = max(ranks))
    u = svd_res$u
    d = svd_res$d
    v = svd_res$v
    
    rm(svd_res)
  } else if(method == 'sparseSVD') {
    pmd_res = do.call(PMA::PMD,c(list(x=as.matrix(mat), K =  max(ranks)), methodParams)) 
    
    u = pmd_res$u
    d = pmd_res$d
    v = pmd_res$v
    
    rm(pmd_res)
  } else {
    return('Bad Method Input')
  }
  
  
  if(!is.null(save_plots_filepath)) { # plot the singular values
    # svd_res$d |> plot()
    ggplot(NULL, aes(x = 1:length(d),
                     y = d)) +
      geom_point() +
      labs(x = 'singular value #', y = 'singular value', 
           title = 'Singular Values') +
      theme_bw() +
      theme(panel.grid.minor = element_blank())
    ggsave(sprintf('%s/singvals.pdf', save_plots_filepath), height = 4, width = 6)
    
  }
  
  
  
  
  # rank = r <= max(ranks)
  approxmatrices = list()
  for(r in ranks) { # there is a faster ay by adding the ne eigenvecs to prev...but it is ok
    # approximation should be U D V^T
    if(r == 1) { # r drops the vector/matrix form...
      diag_singvals = matrix(d[1], nrow=1, ncol=1)   
    } else {
      diag_singvals = diag(d[1:r]) 
    }
    
    approx = u[, 1:r, drop=FALSE] %*% diag_singvals %*% t(v[, 1:r, drop=FALSE])
    colnames(approx) = colnames(mat)
    rownames(approx) = rownames(mat)
    
    approxmatrices[[r]] = approx
    
    # plot original vs approx (side by side)? or just the approx? i think just 1 at a time
    if(!is.null(save_plots_filepath)) {
      # create df for plotting 
      plot_df = approx |> as.data.frame() |> tibble::rownames_to_column('grna') |> 
        tidyr::pivot_longer(cols = 2:(ncol(approx)+1), 
                            names_to = 'gene', values_to = 'value') |>
        mutate(group = sprintf('rank=%02d', r))
      plot_df = merge(gene_index, plot_df, by = 'gene')
      plot_df = merge(grna_index, plot_df, by = 'grna')
      
      # one plot at a time (individual plots for eah rank)
      if(save_individual_rank_plots) {
        ggplot(plot_df) +
          geom_raster(aes(x = grna_idx, y = gene_idx, fill = value)) +
          scale_x_continuous(expand = c(0, 0)) +
          scale_y_continuous(expand = c(0, 0)) +
          labs(x='grna', y = 'gene', fill = NULL) +
          scale_fill_gradient2(limits = color_limits, # set color limits
                               oob=scales::squish, # if outside lims, set to limits
                               midpoint = 0, 
                               high = myRed, low = myBlue, mid = 'white',
                               breaks = color_breaks,
                               labels = color_breaks_label) +
          facet_wrap(vars(group), nrow=1) +
          theme_bw() +
          theme(strip.background = element_rect(fill = 'white'), 
                axis.ticks = element_blank(), 
                axis.text = element_blank(), 
                legend.position = 'bottom',
                legend.key.height = unit(.3, 'cm'),
                legend.key.width  = unit(1, 'cm'),
                legend.text = element_text(size = 7))
        
        ggsave(sprintf('%s/rank=%02d.pdf', save_plots_filepath, r), height = 6, width = 3)
      }
      
      
      plot_df_all = rbind(plot_df_all, plot_df)
      
      
    }
  }
  
  
  # plot original and all lowrank approximations side by side
  if(!is.null(save_plots_filepath)) {
    
    plot_df_all$group = factor(
      plot_df_all$group, 
      levels = c('original', sprintf('rank=%02d', sort(ranks, decreasing=TRUE)))
    )
    
    
    
    ggplot(plot_df_all) +
      geom_raster(aes(x = grna_idx, y = gene_idx, fill = value)) +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0)) +
      labs(x='grna', y = 'gene', fill = NULL) +
      scale_fill_gradient2(limits = color_limits, # set color limits
                           oob=scales::squish, # if outside lims, set to limits
                           midpoint = 0,
                           high = myRed, low = myBlue, mid = 'white',
                           # low  = brewer.pal(n = 9, name = "RdBu")[9],
                           # high = brewer.pal(n = 9, name = "RdBu")[1],
                           breaks = color_breaks,
                           labels = color_breaks_label) +
      facet_wrap(vars(group), nrow=1) +
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
    
    ggsave(sprintf('%s/rankall.pdf', save_plots_filepath), height = 6, width = 1.5*(length(ranks) + 1))
    
    
    # return the plot_df_all if it is created... (e.g. if the save plot file path is not NULL)
    # (want to return bc it is useful for shrinkage fn... should prob just make a function for
    #  matrix --> dataframe...)
    return(list(approxmatrices =approxmatrices,
                approxdf=plot_df_all))
  }
  
  
  
  
  return(list(approxmatrices=approxmatrices))
  
}



#' @param mat (matrix) to plot
#' @param row_order (vector) of integers 
#'    c(4, 2, 8) indicates 1st row is rowidx=4
#' @param column_order (vector) of integers
#' @param color_limits (vector) of length 2, for the limits of the colors 
plot_matrix <- function(mat, row_order=NULL, column_order=NULL, color_limits = c(-2, 2)) {
  # plot using the given row and col ordering
  if(!is.null(row_order) & !is.null(column_order)) {
    grna_index = data.frame(grna = rownames(mat)[row_order],
                            grna_idx  = 1:nrow(mat))
    gene_index = data.frame(gene = colnames(mat)[column_order],
                            gene_idx  = 1:ncol(mat))
  } else { # plot matrix as is
    grna_index = data.frame(grna = rownames(mat),
                            grna_idx  = 1:nrow(mat))
    gene_index = data.frame(gene = colnames(mat),
                            gene_idx  = 1:ncol(mat))
  }
  
  # color breaks for plotting
  color_breaks = sort(union(color_limits, seq(from = round(color_limits[1]), to = round(color_limits[2]))))
  color_breaks_label = color_breaks
  color_breaks_label[which.min(color_breaks)] = sprintf('<%.1f', min(color_breaks))
  color_breaks_label[which.max(color_breaks)] = sprintf('>%.1f', max(color_breaks))
  
  # create df for plotting 
  plot_df = mat |> as.data.frame() |> tibble::rownames_to_column('grna') |> 
    tidyr::pivot_longer(cols = 2:(ncol(mat)+1), 
                        names_to = 'gene', values_to = 'value')
  plot_df = merge(gene_index, plot_df, by = 'gene')
  plot_df = merge(grna_index, plot_df, by = 'grna')
  
  # plot matrix
  ggplot(plot_df) +
    geom_raster(aes(x = grna_idx, y = gene_idx, fill = value)) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x='grna', y = 'gene', fill = NULL) +
    scale_fill_gradient2(limits = color_limits, # set color limits
                         oob=scales::squish, # if outside lims, set to limits
                         midpoint = 0, 
                         high = myRed, low = myBlue, mid = 'white',
                         breaks = color_breaks,
                         labels = color_breaks_label) +
    theme_bw() +
    theme(strip.background = element_rect(fill = 'white'), 
          axis.ticks = element_blank(), 
          axis.text = element_blank(), 
          legend.position = 'bottom',
          legend.key.height = unit(.3, 'cm'),
          legend.key.width  = unit(1, 'cm'),
          legend.text = element_text(size = 7))
}





#' Shrink unshrunk matrix to shrinkage point matrix using ebci method 
#' se mat should have the unshrunk matrix estimates' se
#' @param unshrunk_mat (matrix) of unshrunk estimates
#' @param shrinkpoint_mat (matrix) of points to shrink to (e.g. lower dim mat)
#' @param se_mat (matrix) of standard errors for unshrunk matrices' estimates
#' @param ALPHA (numeric) in [0,1] test level for constructing ebci CI's
#' @param return_ebci_obj (boolean) whether to return the ebci object  
#' @param weight_mat (matrix) optional weight matrix for ebci shrinkage, else 1/n()
#' @output dataframe of shrinkage results with at least the following columns
#' "grna"            "gene"            "unshrunk_value"  "shrinkage_point" "se"
#' "shrunk_value"    "lower_ci"        "upper_ci" 
#' @example 
#'  shrink_matrix(   unshrunk_mat = matrices$estimates,
#'              shrinkpoint_mat = sparseSVD_estimate$approxmatrices[[15]],
#'                       se_mat = se_mat,
#'              ALPHA = .1) # ~15 mins for 50x2000 
shrink_matrix <- function(unshrunk_mat,
                          shrinkpoint_mat,
                          se_mat,
                          ALPHA,
                          return_ebci_obj=FALSE,
                          weight_mat=NULL) {
  
  
  unshrunk       = unshrunk_mat|> as.data.frame() |> 
    tibble::rownames_to_column('grna') |>
    tidyr::pivot_longer(cols = 2:(ncol(unshrunk_mat)+1),
                        names_to = 'gene', values_to = 'unshrunk_value')
  
  
  shrinkagepoint = shrinkpoint_mat |> as.data.frame() |>
    tibble::rownames_to_column('grna') |>
    tidyr::pivot_longer(cols = 2:(ncol(shrinkpoint_mat)+1),
                        names_to = 'gene', values_to = 'shrinkage_point')
  
  se             = se_mat |> as.data.frame() |> 
    tibble::rownames_to_column('grna') |> 
    tidyr::pivot_longer(cols = 2:(ncol(se_mat)+1), 
                        names_to = 'gene', values_to = 'se')
  

  ebci_data = merge(merge(unshrunk, shrinkagepoint, 
                          by = c('grna', 'gene')),
                    se, by = c('grna', 'gene'))
  
  # add weights 
  if(!is.null(weight_mat)) {
    weights_df = weight_mat |> as.data.frame() |> 
      tibble::rownames_to_column('grna') |> 
      tidyr::pivot_longer(cols = 2:(ncol(se_mat)+1), 
                          names_to = 'gene', values_to = 'weight')
    ebci_data = merge(ebci_data, weights_df,  by = c('grna', 'gene'))
  } else {
    ebci_data$weight = 1/nrow(ebci_data)
  }
  
  # remove tests with missing values (using unshrunk_value column - from unshrunk_mat input)
  ebci_data = ebci_data |> filter(!is.na(unshrunk_value))
  
  # t0 = Sys.time()
  ebci_obj = ebci::ebci(formula = sprintf('%s - shrinkage_point ~ 0', 'unshrunk_value'),
                        # formula = 'tstat - shrinkage_point ~ 0', # ebci_formula
                        data    = ebci_data, # |> dplyr::mutate(my_weights = 1/n()), 
                        # se = eval(parse(text = se_colname)), 
                        # weights =  eval(parse(text = sprintf('1/%s^2', se_colname))), 
                        # TODO: set weights, prob 1/n
                        se = se,
                        weights = weight,
                        # weights = 1/nrow(ebci_data),
                        alpha = ALPHA)
  # t1 = Sys.time(); print(t1 - t0)
  
  ebci_res = ebci_data |> 
    # head(nrow(ebci_obj$df)) |>
    dplyr::mutate(shrunk_value = ebci_obj$df$th_eb + shrinkage_point,
                  lower_ci = shrunk_value - ebci_obj$df$len_eb,
                  upper_ci = shrunk_value + ebci_obj$df$len_eb)
  if(return_ebci_obj) {
    return(list(ebci_res = ebci_res,
                ebci_obj = ebci_obj))
  } else {
    return(ebci_res)
  }
  
}







#' old
#' @param df (dataframe) of ebci results, should have the following columns:
#' "grna"            "gene"            "unshrunk_value"  "shrinkage_point" "se"
#' "shrunk_value"    "lower_ci"        "upper_ci" 
#' @output list of plots. you can plot them all together with 
#'         gridExtra::grid.arrange(grobs = somePlots)
shrink_matrix_plots0 <- function(df) {
  
  plots = list()
  
  plots$se_hist = ggplot(df, aes(x = se)) + geom_histogram() + theme_classic()
  
  plots$shrinkagepoint_hist = ggplot(df, aes(x = shrinkage_point)) + geom_histogram() + theme_classic()
  
  
  # plot shrunk vs unshrunk
  plots$shrunk_unshrunk_point = 
    ggplot(df) +
    geom_abline(aes(slope = 1, intercept = 0), color = 'gray') +
    geom_point(aes(x = unshrunk_value, y = shrunk_value, color = se)) +
    theme_classic()
  
  # plot shrunk vs unshrunk centered
  plots$shrunk_unshrunk_centered_point = 
    ggplot(df) +
    geom_abline(aes(slope = 1, intercept = 0), color = 'gray') +
    geom_point(aes(x = unshrunk_value - shrinkage_point, y = shrunk_value - shrinkage_point, color = se)) +
    theme_classic()
  
  # plot unshrunk vs shrinkage point
  plots$unshrunk_shrinkagepoint_point =
    ggplot(df) +
    geom_abline(aes(slope = 1, intercept = 0), color = 'gray') +
    geom_point(aes(x = unshrunk_value, y = shrinkage_point, color = se)) +
    theme_classic()
  
  return(plots)
  
}



#' Plot shrinkage results using saved ebci shrinkage res df
#' @param shrink_df (dataframe) of ebci results, should have the following columns:
#' "grna"            "gene"            "unshrunk_value"  "shrinkage_point" "se"
#' "shrunk_value"    "lower_ci"        "upper_ci" 
#' @param plot_folder (character) of plot_folder to save to
#' @param order_rowscols (boolean) whether to order rows and cols or not
#'     if TRUE, either use grna_index and gene_index to order
#'              or if those are unspecified, perform marginal hierarchical clustering
#' @param grna_index (dataframe) specifying grna (row) order. should have 2 columns
#'     'grna' for grna name, and 'grna_idx' for that grna's index
#' @param gene_index (dataframe) specifying gene (col) order. should have 2 columns
#'     'gene' for gene name, and 'gene_idx' for that gene's index
#' @output list of plots. you can plot them all together with 
#'         gridExtra::grid.arrange(grobs = somePlots)
plot_shrink_results <- function(shrink_df, plot_folder, order_rowscols=TRUE, grna_index=NULL, gene_index=NULL) {
  
  # reorder columns and rows
  if(order_rowscols) {
    if(is.null(grna_index) | is.null(gene_index)) { # use given
      # get an ordering for genes/grna (marginal clustering) ---------------------
      # shrink_matrix = shrink_df |> # filter(grna_idx %in% 1:400 & gene_idx %in% 1:400) |> 
      #   dplyr::select(grna, gene, shrunk_value) |>
      #   tidyr::pivot_wider(names_from = gene, values_from = shrunk_value) |>
      #   tibble::column_to_rownames(var='grna')
      
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
    }  # else use pre-specified ordering 
    
  } else {
    grna_index = data.frame(grna = shrink_df$grna|>unique()|>sort()) |> mutate(grna_idx = 1:n())
    gene_index = data.frame(gene = shrink_df$gene|>unique()|>sort()) |> mutate(gene_idx = 1:n())
  }
  
  
  
  
  
  # color breaks for plotting ------------------------------------------------
  color_limits = c(-2, 2)
  color_breaks = sort(union(color_limits, seq(from = round(color_limits[1]), to = round(color_limits[2]))))
  color_breaks_label = color_breaks
  color_breaks_label[which.min(color_breaks)] = sprintf('<%.1f', min(color_breaks))
  color_breaks_label[which.max(color_breaks)] = sprintf('>%.1f', max(color_breaks))
  
  
  
  
  # add original estimates ---------------------------------------------------
  # (all the loaded dfs should have these) (2x for plotting)
  plot_df = rbind(shrink_df |> 
                    dplyr::mutate(rank = 'shrunk'),
                  shrink_df |>
                    dplyr::mutate(rank = 'unshrunk') |>
                    dplyr::mutate(shrinkage_point = unshrunk_value,   # set all values to the original estimates
                                  shrunk_value    = unshrunk_value, 
                                  lower_ci        = unshrunk_value - qnorm(.95) * se,
                                  upper_ci        = unshrunk_value + qnorm(.95) * se) 
  )
  
  # add idx's- (from gene or grna name --> number on the axis)
  plot_df = merge(gene_index, plot_df, by = 'gene')
  plot_df = merge(grna_index, plot_df, by = 'grna')
  
  
  plot_df$rank = factor(plot_df$rank, levels = c('unshrunk', 'shrunk'))  
  
  
  
  # === Heatmap ===
  
  # heatmap of diff values side by side
  plot_df_heatmap =
    rbind(shrink_df |> dplyr::select(grna, gene, unshrunk_value)  |> dplyr::rename(val =  unshrunk_value) |> dplyr::mutate(type = 'unshrunk'), 
          shrink_df |> dplyr::select(grna, gene, shrinkage_point) |> dplyr::rename(val = shrinkage_point) |> dplyr::mutate(type = 'shrinkagepoint'), 
          shrink_df |> dplyr::select(grna, gene, shrunk_value)    |> dplyr::rename(val =    shrunk_value) |> dplyr::mutate(type = 'shrunk'),
          shrink_df |> dplyr::mutate(significant = (upper_ci <= 0) | (0 <= lower_ci)) |> 
            dplyr::select(grna, gene, significant)    |> dplyr::rename(val =    significant) |> dplyr::mutate(type = 'significant'))
  
  # add true effect if known (e.g. simulations)
  if('true_effect' %in% colnames(shrink_df)) {
    plot_df_heatmap =
      rbind(plot_df_heatmap, 
            shrink_df |> dplyr::select(grna, gene, true_effect)  |> dplyr::rename(val =  true_effect) |> dplyr::mutate(type = 'trueeffect'))
    # also add in coverage of true effect
    plot_df_heatmap =
      rbind(plot_df_heatmap, 
            shrink_df |> 
              dplyr::mutate(true_effect_covered = (lower_ci <= true_effect) & (true_effect <= upper_ci)) |> 
              dplyr::select(grna, gene, true_effect_covered) |> dplyr::rename(val =  true_effect_covered) |> dplyr::mutate(type = 'trueeffectcovered'))
    
    
  }
  
  
  
  # add idx's- (from gene or grna name --> number on the axis)
  plot_df_heatmap = merge(gene_index, plot_df_heatmap, by = 'gene')
  plot_df_heatmap = merge(grna_index, plot_df_heatmap, by = 'grna')
  
  
  plot_df_heatmap$type = factor(plot_df_heatmap$type, levels = c('trueeffect', 'unshrunk', 'shrinkagepoint', 'shrunk', 'significant', 'trueeffectcovered'))  
  
  
  
  # plot
  p_heatmap = ggplot(plot_df_heatmap) +
    geom_raster(aes(x = grna_idx, y = gene_idx, fill = val)) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x='grna', y = 'gene', fill = NULL, title = 'Estimates After Shrinking') +
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
  
  
  ggsave(plot = p_heatmap, filename = sprintf('%sheatmap.pdf', plot_folder), height = 7, width = 7)
  
  
  
  
  
  
  
  
  
  # === Shrinkage Changes ===
  
  set.seed(12345)
  subsample_idx = sample(1:nrow(plot_df), min(10000, nrow(plot_df)))
  
  # plot shrunk vs unshrunk
  p1 = ggplot(plot_df[subsample_idx, ] |> filter(rank != 'unshrunk')) +
    geom_abline(aes(slope = 1, intercept = 0), color = 'gray', alpha = .6) +
    geom_point(aes(x = unshrunk_value, y = shrunk_value, color = se)) +
    labs(title = 'Shrunk vs Unshrunk Estimates', x = 'unshrunk', y = 'shrunk') +
    coord_cartesian(xlim = c(-2, 2), ylim = c(-2, 2)) +
    facet_grid(cols = vars(rank)) +
    theme_classic()
  
  # ggsave(plot = p1, filename = sprintf('%spoints/shrunkvsunshrunk.pdf', plot_folder), height = 6, width = 7)
  
  
  # plot shrunk vs unshrunk centered
  p2 = ggplot(plot_df[subsample_idx, ] |> filter(rank != 'unshrunk')) +
    geom_abline(aes(slope = 1, intercept = 0), color = 'gray') +
    geom_point(aes(x = unshrunk_value - shrinkage_point, y = shrunk_value - shrinkage_point, color = se)) +
    labs(title = 'Shrunk vs Unshrunk Estimates Centered', x = 'unshrunk - shrinkage point', y = 'shrunk - shrinkage point') +
    coord_cartesian(xlim = c(-2, 2), ylim = c(-2, 2)) +
    facet_grid(# rows = vars(approxmethod),
      cols = vars(rank)) +
    theme_classic()
  
  # ggsave(plot = p2, filename = sprintf('%spoints/shrunkvsunshrunkcent.pdf', plot_folder), height = 6, width = 7)
  
  
  # plot unshrunk vs shrinkage point
  p3 = ggplot(plot_df[subsample_idx, ] |> filter(rank != 'unshrunk')) +
    geom_abline(aes(slope = 1, intercept = 0), color = 'gray') +
    geom_point(aes(x = shrunk_value, y = shrinkage_point, color = se)) +
    labs(title = 'Shrinkage Point vs Shrunk Estimates', x = 'shrunk', y = 'shrinkage point') +
    coord_cartesian(xlim = c(-2, 2), ylim = c(-2, 2)) +
    facet_grid(# rows = vars(approxmethod),
      cols = vars(rank)) +
    theme_classic()
  
  # ggsave(plot = p3, filename =sprintf('%spoints/shrinkptvsshrunk.pdf', plot_folder), height = 6, width = 7)
  
  
  
  
  
  
  
  # === Shrinkage CIs ===
  set.seed(12345)
  # get the same gene/grna tests across all methods and ranks (to be better visually, should do to prev plots too)
  sample_tests = plot_df |> filter(rank == 'unshrunk') |> slice_sample(n = 15) |> select(grna, gene) |> mutate(x = 1:n())
  plot_df_sample = merge(plot_df, sample_tests, all.x = FALSE, all.y = TRUE)
  
  # plot_df_sample |> group_by(approxmethod, rank) |> summarize(count = n()) # should all be nrow(sample_tests)
  
  
  p_CI = ggplot(plot_df_sample # |> filter(rank == 'rank=15' & approxmethod == 'lowrank') #|> 
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
    # --- Shrinkage Point---
    geom_point(aes(x = x + .25, y = shrinkage_point), shape=5, color = 'deepskyblue4') +
    # geom_point(aes(y = eval(parse(text = unshrunk_colname))), color = 'deepskyblue2') +
    # geom_point(aes(y = eval(parse(text =   shrunk_colname))), color = 'deepskyblue4') +
    scale_x_continuous(expand = c(0, 0), breaks = seq(0, 1000, by = 10),
                       # limits = c(1, 90)) +
                       # limits = c(65, 175)) +
                       # limits = c(95, 400)) +
                       # limits = c(1, 400)) +
    ) +
    # scale_y_continuous(expand = c(0.025, 0)) +
    coord_cartesian(ylim = c(-2, 2)) +
    labs(x = 'AY Test', y = 'Estimates',
         title = 'Before and After Robust EBCI') +
    facet_grid(# rows = vars(approxmethod),
      cols = vars(rank)) +
    theme_bw() +
    theme(panel.grid.minor.x = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_blank(),
          # legend.position = 'inside',
          # legend.position.inside = legend_position,
          legend.background = element_rect(color = 'black'))
  
  
  ggsave(plot = p_CI, filename = sprintf('%spoints/CIs.pdf', plot_folder), height = 6, width = 15)
  
  
  # 90% CI lengths before and after 
  set.seed(12345)
  p_CI_len = ggplot(plot_df[sample(x=1:nrow(plot_df), size = min(10000, nrow(plot_df))), ] |> 
                      mutate(unshrunk_ci_length = 2*qnorm(.95)*se,
                             shrunk_ci_length = upper_ci - lower_ci)) +
    geom_abline(aes(intercept = 0, slope = 1)) +
    geom_point(aes(x = unshrunk_ci_length, y = shrunk_ci_length), alpha = .8) +
    labs(x = 'Unshrunk CI Lengths', y = 'Shrunk CI Lengths', title = 'CI Lengths') +
    facet_grid(#rows = vars(approxmethod),
      cols = vars(rank)) +
    theme_classic()
  
  
  ggsave(plot = p_CI_len, filename =sprintf('%spoints/CIlen.pdf', plot_folder), height = 6, width = 7)
  
  
  # plot all together
  pdf(sprintf('%sall.pdf', plot_folder), height = 12, width = 12)
  gridExtra::grid.arrange(p_heatmap + theme(legend.position = 'right', legend.key.height = unit(1.75, 'cm'), legend.key.width  = unit(.3, 'cm')), 
                          p1, p2, p3, 
                          p_CI, p_CI_len, layout_matrix = matrix(c(1, 1, 1,
                                                                   1, 1, 1,
                                                                   2, 3, 4,
                                                                   5, 5, 6), byrow = T, nrow=4))
  dev.off()
  
}







#' reorder matrix rows and columns to a specified order 
#' (created for reordering test matrices to train matrices order or just
#'  to the same order)
#' @param mat (matrix) matrix to reorder
#' @param rorder (vector) of strings of row names 
#' @param corder (vector) of strings of column order
#' test = mapply(FUN = my_reorder_rc, mat = matrices, 
#'               MoreArgs = list(rorder = chosengrnas, corder = chosengenes), 
#'               SIMPLIFY=FALSE)
#' 
#' dim(test)
#' length(test)
#' test$estimates |> row.names()
#' test$estimates |> colnames() |> head(100)
#' test$estimates[1:10, 1:10]
#' test$foldchange |> colnames() |> head(100)
#' # make sure order of ros and cols the same 
#' for now, just do alphabetical order (chosengrnas and chosengenes)
#' # make sure the numbers are still the same
#' matrices_r      =  mapply(FUN = my_reorder_rc, mat = matrices, 
#'                         MoreArgs = list(rorder = chosengrnas, corder = chosengenes), 
#'                         SIMPLIFY=FALSE)
#' matrices_test_r =  mapply(FUN = my_reorder_rc, mat = matrices_test, 
#'                          MoreArgs = list(rorder = chosengrnas, corder = chosengenes), 
#'                          SIMPLIFY=FALSE)
#' 
#' mean(unlist(as.vector(matrices$estimates))); mean(unlist(as.vector(matrices_r$estimates)))
#' mean(unlist(as.vector(matrices_test$estimates))); mean(unlist(as.vector(matrices_test_r$estimates)))
#' mean(unlist(as.vector(matrices$foldchange))); mean(unlist(as.vector(matrices_r$foldchange)))
#' mean(unlist(as.vector(matrices_test$foldchange))); mean(unlist(as.vector(matrices_test_r$foldchange)))
my_reorder_rc <- function(mat, rorder, corder) {
  # first check that the set of rows and cols are the same in mat and the inputs
  if(!(all(row.names(mat) %in% rorder) & all(rorder %in% row.names(mat)))) {
    print('Bad mat row names or rorder input')
  }
  if(!(all(colnames(mat) %in% corder) & all(corder %in% colnames(mat)))) {
    print('Bad mat col names or corder input')
  }
  
  # reorder
  return(mat[rorder, corder])
}


#' Spline that does this:
#' 
#'       |                . 
#'       |              .  <- linear, y=x
#'       |            .   
#'epsilon+  -  -  - + <- match 1st deriv
#'       |        . |
#'       |      .     <- quadratic
#' delta +.  .      |
#'       |
#'       -----------+
#'       0          epsilon
#'        
#' To ensure the function does not decrease on [0, epsilon],
#' the parameters must satisfy delta <= epsilon / 2
#' @param delta (numeric) function parameter
#' @param epsilon (numeric) function parameter
#' @example 
#' x = seq(0, 1, length.out = 200) 
#' my_spline = make_spline(delta = .3, epsilon = .6) 
#' plot(x, sapply(X = x, FUN = my_spline), type = 'l')
#' 
#' my_spline = make_spline(delta = .3, epsilon = .4) <- bad params
make_spline <- function(delta, epsilon) {
  a = (delta / epsilon**2) 
  b = (1 - 2 * delta / epsilon)
  fn <- function(x) {
    if(is.na(x)) {
      return(x)
    } else if(x < epsilon) {
      return(a * x**2 + b * x + delta)
    } else {
      return(x)
    }
  }
  return(fn)
} 




my_spline = make_spline(delta = .01, epsilon = .03)

# create new se values (hopefully more robust)
spline_se <- function(p, mu) {
  
  # x = seq(0, 1, length.out = 200) ; plot(x, sapply(X = x, FUN = my_spline), type = 'l')
  my_spline(abs(mu)) / my_spline(abs(qnorm(p = (1/2)*p, mean = 0, sd = 1)))
}





