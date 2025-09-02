# functions for performing SCEPTRE analysis,
# specifically for the cleary dataset 

# Example code =================================================================
# Run like this:
# t0 = Sys.time()
# 
# # label each cell as part of training or testing dataset now
# set.seed(12345)
# cell_names = clearyrds[[]] |> rownames()
# N_cells = length(cell_names)
# cell_train = sample(cell_names, size = floor(N_cells / 2), replace = FALSE) |> sort()
# cell_test  = setdiff(cell_names, cell_train) |> sort()
# 
# # create sceptre objects
# sceptre_obj_train = prepare_sceptre(seurat_obj = clearyrds, cell_subset = cell_train, seed = 12345)
# sceptre_obj_test  = prepare_sceptre(seurat_obj = clearyrds, cell_subset = cell_test , seed = 12345)
# 
# # grna_target_data_frame and response_names combined from train and test
# my_grna_target_data_frame = rbind(sceptre_obj_train$grna_target_data_frame,
#                                   sceptre_obj_test$grna_target_data_frame) |> 
#   distinct()
# my_response_names = unique(c(sceptre_obj_train$response_names,
#                              sceptre_obj_test$response_names))
# 
# my_discovery_pairs = create_discovery_pairs(grna_target_data_frame = my_grna_target_data_frame, 
#                                             response_names = my_response_names,
#                                             NUM_GENES_PER_GRNA=15, NUM_GRNA=NULL,
#                                             seed = 12345)
# 
# my_positive_control_pairs = construct_positive_control_pairs(sceptre_obj_train$sceptre_object)
# 
# # perform sceptre 
# sceptre_obj_train$sceptre_object = 
#   perform_sceptre(sceptre_object = sceptre_obj_train$sceptre_object, 
#                   discovery_pairs = my_discovery_pairs,
#                   positive_control_pairs = my_positive_control_pairs,
#                   save_dir_name = '../saves/sceptre/trainsplit',
#                   save_obj_name = '../saves/sceptre/sceptre_object_train.rds' )
# 
# 
# sceptre_obj_test$sceptre_object = 
#   perform_sceptre(sceptre_object = sceptre_obj_test$sceptre_object, 
#                   discovery_pairs = my_discovery_pairs,
#                   positive_control_pairs = my_positive_control_pairs,
#                   save_dir_name = '../saves/sceptre/testsplit',
#                   save_obj_name = '../saves/sceptre/sceptre_object_test.rds' )
# 
# 
# # get effect estimates and se (inferred from pval)
# sceptre_effects_train = 
#   construct_sceptre_effects_df(sceptre_object = sceptre_obj_train$sceptre_object, 
#                                save_file = '../saves/sceptre/sceptre_effects_train.csv')
# 
# sceptre_effects_test = 
#   construct_sceptre_effects_df(sceptre_object = sceptre_obj_test$sceptre_object, 
#                                save_file = '../saves/sceptre/sceptre_effects_test.csv')
# 
# t1 = Sys.time(); print(t1 - t0)
# ==============================================================================








#' prepare a sceptre object using the specified 
#' seurat object (assumes same structure as clearyrds object from:
#' '../../genData/cleary/GSM6858447_KO_conventional.rds'))
#' @param seurat_obj (object of class seurat)
#' @param cell_subset (vector) of cell names to subset for if wanted
#'                    NULL means no further subsetting
#' @param seed (integer) set the seed for rng 
#' @examples
#' clearyrds = readRDS('../../genData/cleary/GSM6858447_KO_conventional.rds')
#' 
#' set.seed(12345)
#' cell_names = clearyrds[[]] |> rownames()
#' N_cells = length(cell_names)
#' cell_train = sample(cell_names, size = floor(N_cells / 2), replace = FALSE) |> sort()
#' cell_test  = setdiff(cell_names, train_idx) |> sort()
#' 
#' sceptre_obj_train = prepare_sceptre(seurat_obj = clearyrds, cell_subset = cell_train, seed = 12345)
#' sceptre_obj_test  = prepare_sceptre(seurat_obj = clearyrds, cell_subset = cell_test , seed = 12345)
prepare_sceptre <- function(seurat_obj, cell_subset=NULL, seed=12345) {
  set.seed(seed)
  
  # ============================= Prepare Inputs =================================
  # cellid's with 1 grna & (NOT safe-targeting perturbations)
  cells1grna = seurat_obj[[]] |> 
    filter(Total_number_of_guides == 1) |>   # cells w/ 1 perturbation
    filter(Guides_collapsed_by_gene != 'safe-targeting') |> # targeting or non-targeting
    rownames()
  
  if(is.null(cell_subset)) {cell_subset = cells1grna} else {cell_subset = intersect(cell_subset, cells1grna)}
  
  # cell confounders =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  covariate_names = c('Total_RNA_count', '10X_channel', 'S_score', 'G2M_score', 'Cell_cycle_phase') # hmm library size here or later?
  # cell_covariates = cleary$obs[cells1grna, ] |> select(covariate_names) # using h5
  extra_covariates = seurat_obj[[]][cell_subset, ]|> select(all_of(covariate_names)) # using rds #(75043 or some subset) x 5
  
  
  # response matrix =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  response_matrix = seurat_obj[['RNA']]$counts[, cell_subset] # genes x cells, 16952 x (75043 or some subset)
  
  # grna matrix =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  grna_matrix = seurat_obj[['perturbations']]$counts[, cell_subset] # grnas x cells, 600 x (75043 or some subset)
  
  
  # grna_target_data_frame =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  
  # comments: cleary grna counts' rownames already collapse by targeted gene (600) (e.g. instead of individual perturbation)
  # clearyrds[[]] |> filter(Total_number_of_guides == 1) |> pull(Guides_collapsed_by_gene) |> unique() |> length()
  # clearyrds[[]] |> filter(Total_number_of_guides == 1) |> distinct(Guides, Guides_collapsed_by_gene) |> group_by(Guides_collapsed_by_gene) |> summarize(count = n()) |> arrange(desc(count))
  
  grna_target_data_frame = seurat_obj[[]] |>
    filter(Total_number_of_guides == 1) |>
    mutate(grna_id = Guides_collapsed_by_gene, grna_target = Guides_collapsed_by_gene) |>
    select(grna_id, grna_target) |>
    distinct() # 600  x  2
  
  # grna_target_data_frame = clearyrds[[]] |> 
  #                          filter(Total_number_of_guides == 1) |>
  #                          mutate(grna_id = Guides, grna_target = Guides_collapsed_by_gene) |>
  #                          select(grna_id, grna_target) |>
  #                          distinct() # 600  x  2
  
  
  # Fix naming issue?? =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  # caused by: grna_matrix rownames are already collapsed by target gene. So 'non-targeting' is a grna_id, 
  # which is not liked by sceptre https://github.com/Katsevich-Lab/sceptre/blob/main/R/check_functions.R#L1
  # I think the preferred format is: NT1, NT2, ... and then grna_target_data_frame is NT1, NT2, ...\\ non-targeting, non-targeting, non-targeting...
  
  # c('non-targeting' %in% grna_target_data_frame$grna_id, 'non-targeting' %in% grna_target_data_frame$grna_target) # TRUE, TRUE
  # 'non-targeting' %in% rownames(grna_matrix) # there is an issue when 'non-targeting' is a rowname of grna_matrix (indicating grna_id) 
  
  rownames(grna_matrix)[which(rownames(grna_matrix) == 'non-targeting')] = 'non-targeting1'
  grna_target_data_frame[which(grna_target_data_frame$grna_id == 'non-targeting'), 'grna_id'] = 'non-targeting1' 
  
  # SCEPTRE (`check_calibration_check_inputs`) wants >=2 non-targeting grnas... 
  # we only have 1 type. For now, just randomly separate them, but must ignore the results
  
  
  #     manipulate grna_matrix
  nt2_idx = (grna_matrix['non-targeting1', ] == 1) |> which()
  nt2_idx = sample(x = nt2_idx, size = floor(length(nt2_idx)/2), replace = FALSE)
  
  nt2_grna = matrix(0, nrow = 1, ncol = ncol(grna_matrix)) |> as("sparseMatrix")
  nt2_grna[1, nt2_idx] = 1
  rownames(nt2_grna) = c('non-targeting2')
  
  grna_matrix2 = rbind(grna_matrix, nt2_grna) 
  grna_matrix2['non-targeting1', nt2_idx] = 0
  
  grna_matrix = grna_matrix2
  
  #     manipulate grna_target_data_frame
  grna_target_data_frame = rbind(grna_target_data_frame,
                                 data.frame(grna_id = 'non-targeting2',
                                            grna_target = 'non-targeting'))
  
  
  
  
  
  # rownames(grna_matrix)[which(! rownames(grna_matrix) %in% grna_target_data_frame$grna_id)]
  
  # response_names =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  response_names = rownames(response_matrix)
  
  
  
  
  
  # ============================= SCEPTRE ========================================
  
  
  sceptre_object <- sceptre::import_data(
    response_matrix        = response_matrix,
    grna_matrix            = grna_matrix,
    grna_target_data_frame = grna_target_data_frame,
    moi                    = 'low',
    extra_covariates       = extra_covariates,
    response_names         = response_names
  )
  
  
  return(list(        sceptre_object = sceptre_object,
                      grna_target_data_frame = grna_target_data_frame,
                      response_names = response_names))
}

#' create dataframe of discovery pairs (perturbation - genes) for testing
#' @param grna_target_data_frame (dataframe) of grna and their target info
#'                               (should be the input to sceptre's fn)
#' (potential bug if input is from training subset and is missing or has extra
#'  for testing subset... when inputting the arg, maybe take union of 
#'  test and train)
#' @param response_names (vector) of string of gene response names
#' @param NUM_GENES_PER_GRNA (integer) for each grna (perturbation), how many 
#'                           gene responses do we consider
#' @param NUM_GRNA (integer) how many grna (perturbation)s do we consider
#'                 default NULL will use all grna's
#' @param seed (integer) seed for rng
#' @param all_combos (boolean) if true, choose randomly NUM_GRNA grnas and 
#'  NUM_GENES_PER_GRNA genes. Then choose all combinations to fill a full
#'  NUM_GRNA x NUM_GENES_PER_GRNA matrix 
#' @param must_include_grna (vector) if all_combos=TRUE
#' @param must_include_gene (vector) if all_combos=TRUE
#' @param prioritize_genes (vector) if all_combos=TRUE, then fill remaining 
#'  genes with top genes of these genes (input is given in order of importance)
#'  potential bug if this does not include enough genes (e.g. result will not 
#'  give the specified number of NUM_GENES_PER_GRNA and/or bug out, index out
#'  of bounds)
create_discovery_pairs <- function(grna_target_data_frame, 
                                   response_names,
                                   NUM_GENES_PER_GRNA=10, NUM_GRNA=NULL,
                                   seed = 12345, 
                                   all_combos = FALSE,
                                   must_include_grna = NULL,
                                   must_include_gene = NULL,
                                   prioritize_genes  = NULL) {
  
  assertthat::assert_that(is.null(NUM_GRNA) || (NUM_GRNA > 0 && NUM_GRNA%%1 == 0),
                          msg = "NUM_GRNA must be NULL or a positive integer")
  assertthat::assert_that(NUM_GENES_PER_GRNA > 0 && NUM_GENES_PER_GRNA%%1 == 0,
                          msg = "NUM_GENES_PER_GRNA must be a positive integer")
  assertthat::assert_that(seed > 0 && seed%%1 == 0,
                          msg = "seed must be a positive integer")
  
  set.seed(seed)
  
  # just take some grnas and some genes
  possible_grnas = grna_target_data_frame |> 
    filter(grna_target != 'non-targeting' & grna_target != 'safe-targeting') |> 
    pull(grna_target)
  
  if(is.null(NUM_GRNA)) {NUM_GRNA = length(possible_grnas) } # ~600
  
  if(all_combos) {
    # choose grnas
    if(is.null(must_include_grna)) {
      chosen_grna  = sample(x = possible_grnas, size = NUM_GRNA, replace = FALSE)
    } else {
      chosen_grna = c(must_include_grna,
                      sample(x = setdiff(possible_grnas, must_include_grna), 
                             size = NUM_GRNA - length(must_include_grna), replace = FALSE))
    }
    
    # choose genes
    chosen_genes = c()
    if(!is.null(must_include_gene)) {
      chosen_genes = c(chosen_genes, must_include_gene)
    } 
    
    if(is.null(prioritize_genes)) {
      chosen_genes = c(chosen_genes, 
                       sample(x = setdiff(response_names, chosen_genes), 
                              size = NUM_GENES_PER_GRNA, 
                              replace = FALSE))
    } else {
      chosen_genes = c(chosen_genes, 
                       (setdiff(prioritize_genes, chosen_genes))[1:(NUM_GENES_PER_GRNA - length(chosen_genes))])
    }
    
    # combine arrange in dataframe (row=(grna, gene))
    discovery_pairs = data.frame(
      grna_target = rep(chosen_grna,  each = length(chosen_genes)),
      response_id = rep(chosen_genes, times = length(chosen_grna))
    )
    
  } else {
    discovery_pairs = NULL
    for(i in 1:NUM_GRNA) {
      discovery_pairs = rbind(discovery_pairs,
                              data.frame(grna_target = possible_grnas[i], 
                                         response_id = sample(x    = setdiff(response_names, possible_grnas[i]), 
                                                              size = NUM_GENES_PER_GRNA, replace = FALSE)))
    }
    # discovery_pairs[which(! discovery_pairs$grna_target %in% grna_target_data_frame$grna_target), ]
  }

  return(discovery_pairs)
}



#' perform sceptre
#' 
#' @param sceptre_object (sceptre object)
#' @param discovery_pairs (dataframe) of grna_targets to gene names
#' @param positive_control_pairs (dataframe) of positive controls
#' @param save_dir_name (string) path to save sceptre info
#' @param save_obj_name (string) filename to save sceptre object
#' @param parallel (boolean) run sceptre command's parallel arg
perform_sceptre <- function(sceptre_object, discovery_pairs, positive_control_pairs,
                            save_dir_name = '../saves/sceptre/default',
                            save_obj_name = '../saves/sceptre/default/sceptre_object_test.rds',
                            parallel=TRUE) {
  
  
  # # construct positive grna-gene pairs to analyze
  # positive_control_pairs <- construct_positive_control_pairs(sceptre_object)
  
  
  # construct discovery grna-gene pairs to analyze
  # DOES NOT WORK -- names of genes are in a diff format (ENSG00000243485 vs TXNDC17)
  # could convert but will just make own pairs
  # discovery_pairs <- construct_cis_pairs(
  #   sceptre_object, 
  #   positive_control_pairs = positive_control_pairs
  # )
  
  # rownames(response_matrix)[1:10]
  # grna_target_data_frame$grna_id[1:10]
  # response_names[1:10]
  # 
  # head(gene_position_data_frame_grch38)
  # gene_position_data_frame_grch38
  # grna_target_data_frame
  # 
  # (rownames(response_matrix) %in% gene_position_data_frame_grch38$response_id) |> table()
  # (rownames(response_matrix) %in% gene_position_data_frame_grch38$response_name) |> table()
  # data(highmoi_example_data)
  # data(grna_target_data_frame_highmoi)
  # # import data
  # sceptre_object <- import_data(
  #   response_matrix = highmoi_example_data$response_matrix,
  #   grna_matrix = highmoi_example_data$grna_matrix,
  #   grna_target_data_frame = grna_target_data_frame_highmoi,
  #   moi = "high",
  #   extra_covariates = highmoi_example_data$extra_covariates,
  #   response_names = highmoi_example_data$gene_names
  # )
  # positive_control_pairs <- construct_positive_control_pairs(sceptre_object)
  # discovery_pairs <- construct_cis_pairs(sceptre_object,
  #   positive_control_pairs = positive_control_pairs,
  #   distance_threshold = 5e6
  # )
  
  
  
  
  
  # apply the pipeline functions to the sceptre_object in order
  sceptre_object <- sceptre_object |> # |> is R's base pipe, similar to %>%
    set_analysis_parameters(discovery_pairs = discovery_pairs, 
                            positive_control_pairs = positive_control_pairs,
                            side = "both") |>
    assign_grnas(parallel = parallel,
                 method = 'maximum', 
                 min_grna_n_umis_threshold = 1, 
                 umi_fraction_threshold = .8) |>
    run_qc() 
  
  
  
  sceptre_object = sceptre_object |>
    run_calibration_check(parallel=parallel)
  
  sceptre_object = sceptre_object |>
    run_power_check(parallel=parallel)
  
  sceptre_object = sceptre_object  |>
    run_discovery_analysis(parallel=parallel)
  
  
  
  
  # ==============================================================================
  # 8. Write outputs to directory
  # https://timothy-barry.github.io/sceptre-book/sceptre.html#sec-sceptre_write_outputs_to_directory
  # ==============================================================================
  dir.create(save_dir_name, recursive = TRUE)
  # list.files('../../saves/sceptre')
  write_outputs_to_directory(
    sceptre_object = sceptre_object, 
    directory = save_dir_name
  )
  
  
  
  
  # print(sceptre_object) # output suppressed for brevity
  # plot(sceptre_object) 
  
  saveRDS(sceptre_object, file = save_obj_name)
  
  
  return(sceptre_object)
}






#' From a p-value and mean, get the standard error,
#' assuming that the p-value is calculated by:
#' 2*pnorm(-(mean - 0)/se, mean = 0, sd = 1)
#' 
#' @param p (numeric) the p-value
#' @param mu (numeric) the mean
#' @examples
#' estimated_p = .025
#' estimated_mu = .2
#' z = qnorm(p = 2*estimated_p, mean = 0, sd = 1)
#' 
#' estimated_se = abs(estimated_mu / z)
#' 
#' # double check
#' pnorm(-abs(estimated_mu) / estimated_se, mean = 0, sd = 1)
#' 
#' get_estimated_se(estimated_p, estimated_mu)
#' pnorm(-abs(estimated_mu) / estimated_se)
#' 
#' mapply(FUN = get_se, p = c(.02, .03, .6), mu = c(1, 2, .1))
get_se <- function(p, mu) {
  abs(mu / qnorm(p = (1/2)*p, mean = 0, sd = 1))
} 



#' get effect estimates and se (inferred from pval)
#' and save if wanted
#' 
#' @param sceptre_object (sceptre object) with the resulting effect estimates 
#'                       and p values
#' @param save_file (string) path to save the effects dataframe at 
#'                  NULL = do not save
construct_sceptre_effects_df <- function(sceptre_object,
                                         save_file = NULL) {
  # get effect estimates and se (inferred from pval)
  sceptre_effects = rbind(sceptre_object@calibration_result |> mutate(test = 'negative') |> relocate(test),
                          sceptre_object@power_result |> mutate(test = 'positive') |> relocate(test),
                          sceptre_object@discovery_result |> mutate(test = 'discovery') |> relocate(test), 
                          fill = TRUE) |>
    rename(estimate = log_2_fold_change) # |>
    # filter(pass_qc)
  
  sceptre_effects$se = 
    mapply(FUN = get_se, 
           p   = sceptre_effects$p_value, 
           mu  = sceptre_effects$estimate) 
  
  # some p=1 --> se=infty. Just replace with a large number for now
  sceptre_effects$se = pmin(sceptre_effects$se, 999)
  
  # rename columns 
  sceptre_effects = sceptre_effects |> rename(pvalue = p_value,
                                              grna = grna_target,
                                              gene = response_id)
  
  if(!is.null(save_file)) {
    write.csv(sceptre_effects, file = save_file)
  }
  return(sceptre_effects)
}




#' Get GLM estimates on specified tests using information from sceptre objects result
#'# non-targeting idx
#' nt_idx = which(sceptre_object@grna_matrix[[1]]['non-targeting1', ] == 1  |
#'                sceptre_object@grna_matrix[[1]]['non-targeting2', ] == 1  )
#' 
#' # negative
#' which(sceptre_object@calibration_result$pass_qc)
#' get_glm_estimate(test_idx = 1,
#'                  test_info = sceptre_object@calibration_result,
#'                  test_type = 'negative',
#'                  grna = sceptre_object@grna_matrix,
#'                  gene = sceptre_object@response_matrix,
#'                  covariates = sceptre_object@covariate_data_frame,
#'                  nt_idx= nt_idx,
#'                  which_glm = which_glm
#'                  )
#' 
#' # positive
#' which(sceptre_object@power_result$pass_qc)
#' get_glm_estimate(test_idx = 1,
#'                  test_info = sceptre_object@power_result,
#'                  test_type = 'positive',
#'                  grna = sceptre_object@grna_matrix,
#'                  gene = sceptre_object@response_matrix,
#'                  covariates = sceptre_object@covariate_data_frame,
#'                  nt_idx= nt_idx,
#'                  which_glm = which_glm
#'                  )
#' # discovery
#' which(sceptre_object@calibration_result$pass_qc)
#' get_glm_estimate(test_idx = 1,
#'                  test_info = sceptre_object@discovery_result,
#'                  test_type = 'discovery',
#'                  grna = sceptre_object@grna_matrix,
#'                  gene = sceptre_object@response_matrix,
#'                  covariates = sceptre_object@covariate_data_frame,
#'                  nt_idx= nt_idx,
#'                  which_glm = which_glm
#'                  )
get_glm_estimate <- function(test_idx, test_info, test_type, grna, gene, covariates, nt_idx, which_glm) {
  grna_name = test_info[[test_idx, 'grna_target']] # grna_name
  gene_name = test_info[[test_idx, 'response_id']] # gene_name
  
  # targeting idx
  grna_idx = which(grna[[1]][grna_name, ] == 1)
  
  # for calibration/negative tests, nt are used for controls AND perturbations,
  # take the set difference for controls (should not change anything if perturbation is not non-targeting)
  control_idx = setdiff(nt_idx, grna_idx)
  
  # gene response
  # gene_response = sceptre_object@response_matrix[[1]][gene_name, c(grna_idx, control_idx)] 
  gene_response = gene[[1]][gene_name, c(grna_idx, control_idx)] 
  
  # construct df for analysis
  df = cbind(data.frame(perturb = c(rep(1, length(grna_idx)), rep(0, length(control_idx))),
                        count   = gene_response),
             covariates[c(grna_idx, control_idx), ])
  
  # colnames(df)
  # which covariates do we want to include??
  # "perturb"            "count"              "response_n_nonzero" "response_n_umis"    "response_p_mito"    "grna_n_nonzero"    
  # "grna_n_umis"        "Total_RNA_count"    "10X_channel"        "S_score"            "G2M_score"          "Cell_cycle_phase"  
  
  # negative binomial or poisson fit depending on which_glm value
  tryCatch(
    {switch(which_glm,
            'negative_binomial'={fit = MASS::glm.nb(count ~ perturb + log(Total_RNA_count) + `10X_channel` + S_score + G2M_score + Cell_cycle_phase, 
                                                    data = df)},
            'poisson'          ={fit =          glm(count ~ perturb + log(Total_RNA_count) + `10X_channel` + S_score + G2M_score + Cell_cycle_phase, 
                                                    data = df, 
                                                    family = poisson())},
            stop("'which_glm' should have value 'negative_binomial' or 'poisson'"))
      
      newrow = data.frame(summary(fit)$coefficients[2, c(1, 2, 4), drop=FALSE],
                          warning=NA)
    },
    error=function(e) { # erred
      print(e)
      assign('newrow',
             data.frame('Estimate'=NA, 'Std..Error'=NA, 'Pr...z..'=NA,
                        warning=NA),
             envir = globalenv())
    }# ,
    # warning=function(w) { # warning, but still can get results
    #     # print(w)
    #     newrow = data.frame(summary(fit)$coefficients[2, c(1, 2, 4), drop=FALSE],
    #                         warning=TRUE)
    # }
  )
  
  
  res = data.frame(test = test_type,
                   gene = gene_name, 
                   grna = grna_name,
                   newrow)
  rm(fit, newrow)
  
  return(res)
}


#' get glm estimates all at once
#' @examples
#' sceptre_obj_train = readRDS('../saves/sceptre/sceptre_object_train.rds')
#' get_glm_estimate_df(sceptre_object = sceptre_obj_train$sceptre_object,
#'                     which_glm = 'negative_binomial',
#'                     save_file = '../saves/negativebinomial_effects_train.csv')
get_glm_estimate_df <- function(sceptre_object, which_glm, save_file=NULL) {
  # non-targeting idx
  nt_idx = which(sceptre_object@grna_matrix[[1]]['non-targeting1', ] == 1  |
                   sceptre_object@grna_matrix[[1]]['non-targeting2', ] == 1  )
  
  # get glm effects from the function
  df_effects = dplyr::bind_rows( 
    # negative
    lapply(X = which(sceptre_object@calibration_result$pass_qc), 
           FUN = get_glm_estimate, 
           test_info = sceptre_object@calibration_result,
           test_type = 'negative',
           grna = sceptre_object@grna_matrix,
           gene = sceptre_object@response_matrix,
           covariates = sceptre_object@covariate_data_frame,
           nt_idx= nt_idx,
           which_glm = which_glm)
    ,
    # positive
    lapply(X = which(sceptre_object@power_result$pass_qc), 
           FUN = get_glm_estimate, 
           test_info = sceptre_object@power_result,
           test_type = 'positive',
           grna = sceptre_object@grna_matrix,
           gene = sceptre_object@response_matrix,
           covariates = sceptre_object@covariate_data_frame,
           nt_idx= nt_idx,
           which_glm = which_glm)
    ,
    # discovery
    lapply(X = which(sceptre_object@discovery_result$pass_qc), 
           FUN = get_glm_estimate, 
           test_info = sceptre_object@discovery_result,
           test_type = 'discovery',
           grna = sceptre_object@grna_matrix,
           gene = sceptre_object@response_matrix,
           covariates = sceptre_object@covariate_data_frame,
           nt_idx= nt_idx,
           which_glm = which_glm)) |>
    `colnames<-` (c('test', 'gene', 'grna', 'estimate', 'se', 'pvalue', 'warning')) |>
    `rownames<-` (NULL)
  if(!is.null(save_file)) {write.csv(df_effects, file = save_file)}
  return(df_effects)
}



