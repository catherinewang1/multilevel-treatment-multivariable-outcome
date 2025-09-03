# code from shrinkEstimates.qmd but in .r file to run it more easily

# =================== Start ========================================================================
print(sprintf("[%s] START: SCEPTRE to get effect sizes", Sys.time()))


# =================== load ====================================================
print(sprintf("[%s]    - load libraries and data", Sys.time()))
suppressPackageStartupMessages(library(dplyr))
# devtools::install_github("katsevich-lab/sceptre")
suppressPackageStartupMessages(library(sceptre))

source('../utils/perform_sceptre_cleary.r')

gene_dev_df = read.csv('../saves/gene_deviance.csv') # df of gene importance (cols: idx, deviance, gene_name)
clearyrds = readRDS('../../../genData/cleary/GSM6858447_KO_conventional.rds')




# =================== SCEPTRE ====================================================

print(sprintf("[%s]    - start SCEPTRE", Sys.time()))

t0 = Sys.time() 

# label each cell as part of training or testing dataset now
set.seed(12345)
cell_names = clearyrds[[]] |> rownames()
N_cells = length(cell_names)
cell_train = sample(cell_names, size = floor(N_cells / 2), replace = FALSE) |> sort()
cell_test  = setdiff(cell_names, cell_train) |> sort()

# create sceptre objects
sceptre_obj_train = prepare_sceptre(seurat_obj = clearyrds, cell_subset = cell_train, seed = 12345)
sceptre_obj_test  = prepare_sceptre(seurat_obj = clearyrds, cell_subset = cell_test , seed = 12345)
sceptre_obj_all   = prepare_sceptre(seurat_obj = clearyrds, cell_subset = NULL      , seed = 12345)

# grna_target_data_frame and response_names combined from train and test
my_grna_target_data_frame = rbind(sceptre_obj_train$grna_target_data_frame,
                                  sceptre_obj_test$grna_target_data_frame,
                                  sceptre_obj_all$grna_target_data_frame) |> 
  distinct()
my_response_names = unique(c(sceptre_obj_train$response_names,
                             sceptre_obj_test$response_names,
                             sceptre_obj_all$response_names))

# positive tests w/ perturbs > THRESHOLD_CELLS_PER_PERTURB = 150 cells per perturb
my_positive_control_pairs = construct_positive_control_pairs(sceptre_obj_all$sceptre_object)
THRESHOLD_CELLS_PER_PERTURB = 100
grna_counts = clearyrds[['perturbations']]$counts |> apply(MARGIN=1, FUN=sum)
my_positive_control_pairs = my_positive_control_pairs |> filter(grna_target %in% names(grna_counts[grna_counts > THRESHOLD_CELLS_PER_PERTURB])) # 286 when #cells>150, 450 when >100


# # random genes per grna, ~ 9000 ... ~ 21 mins
# my_discovery_pairs = create_discovery_pairs(grna_target_data_frame = my_grna_target_data_frame, 
#                                            response_names = my_response_names,
#                                            NUM_GENES_PER_GRNA=15, NUM_GRNA=NULL,
#                                            seed = 12345)

# vs random genes and grna, take all combos: 200 x 100 = 20000 tests, 34 mins
# 500 x 300 = 150000 tests, 3.3 hours
# 286 x 2000 = 572000 tests...
# much faster in parallel, can do more: 450 grnas x 3000 genes
my_discovery_pairs = create_discovery_pairs(grna_target_data_frame = my_grna_target_data_frame, 
                                            response_names = my_response_names,
                                            NUM_GENES_PER_GRNA=10000, NUM_GRNA=nrow(my_positive_control_pairs),
                                            seed = 12345, 
                                            all_combos = TRUE,
                                            must_include_grna = my_positive_control_pairs$grna_target,
                                            must_include_gene = my_positive_control_pairs$response_id,
                                            prioritize_genes = gene_dev_df |> arrange(desc(deviance)) |> pull(gene_name))

write.csv(my_positive_control_pairs, '../saves/sceptre/positive_control_pairs.csv', row.names = FALSE)
write.csv(my_discovery_pairs       , '../saves/sceptre/discovery_pairs.csv'       , row.names = FALSE)


# dim(my_discovery_pairs)
# dim(my_discovery_pairs |> distinct())
# my_discovery_pairs |> group_by(grna_target) |> summarize(count = n()) |> pull(count) |> unique()
# my_discovery_pairs |> group_by(response_id) |> summarize(count = n()) |> pull(count) |> unique()
# gene_dev_df$gene_name[1:100] %in% my_positive_control_pairs[1:300, ]$response_id
# gene_dev_df$gene_name[1:100] %in% my_discovery_pairs$response_id


# perform sceptre 
sceptre_obj_train$sceptre_object = 
  perform_sceptre(sceptre_object = sceptre_obj_train$sceptre_object, 
                  discovery_pairs = my_discovery_pairs,
                  positive_control_pairs = my_positive_control_pairs,
                  save_dir_name = '../saves/sceptre/trainsplit',
                  save_obj_name = '../saves/sceptre/sceptre_object_train.rds',
                  parallel=TRUE)


sceptre_obj_test$sceptre_object = 
  perform_sceptre(sceptre_object = sceptre_obj_test$sceptre_object, 
                  discovery_pairs = my_discovery_pairs,
                  positive_control_pairs = my_positive_control_pairs,
                  save_dir_name = '../saves/sceptre/testsplit',
                  save_obj_name = '../saves/sceptre/sceptre_object_test.rds',
                  parallel=TRUE)

sceptre_obj_all$sceptre_object = 
  perform_sceptre(sceptre_object = sceptre_obj_all$sceptre_object, 
                  discovery_pairs = my_discovery_pairs,
                  positive_control_pairs = my_positive_control_pairs,
                  save_dir_name = '../saves/sceptre/nosplit',
                  save_obj_name = '../saves/sceptre/sceptre_object_all.rds',
                  parallel=TRUE)


# get effect estimates and se (inferred from pval)
sceptre_effects_train = 
  construct_sceptre_effects_df(sceptre_object = sceptre_obj_train$sceptre_object, 
                               save_file = '../saves/sceptre/sceptre_effects_train.csv')

sceptre_effects_test = 
  construct_sceptre_effects_df(sceptre_object = sceptre_obj_test$sceptre_object, 
                               save_file = '../saves/sceptre/sceptre_effects_test.csv')

sceptre_effects_all = 
  construct_sceptre_effects_df(sceptre_object = sceptre_obj_all$sceptre_object, 
                               save_file = '../saves/sceptre/sceptre_effects_all.csv')

t1 = Sys.time(); print(t1 - t0)



# =================== END ==========================================================================
print(sprintf("[%s] END", Sys.time()))


