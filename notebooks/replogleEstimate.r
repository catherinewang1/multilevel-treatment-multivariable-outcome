# SCEPTRE on Replogle rpe1 dataset 
# https://plus.figshare.com/articles/dataset/_Mapping_information-rich_genotype-phenotype_landscapes_with_genome-scale_Perturb-seq_Replogle_et_al_2022_processed_Perturb-seq_datasets/20029387

# Params
# Load data (anndata format h5ad)
data_path = '../../../genData/replogle'
# data_file = 'K562_essential_raw_singlecell_01.h5ad'
data_file = 'rpe1_raw_singlecell_01.h5ad'

# Quality Control Params
MIN_GRNA_COUNT = 100            # perturbations must be taken by at least X cells
MIN_N_NONZERO_PER_GENE = 175000 # genes must have >= X cells w nonzero counts

# SCEPTRE Params
PARALLEL=TRUE # sceptre parellel-set true on server

# ==============================================================================
#                    Quality Control- grnas and genes
# ==============================================================================

print(sprintf("[%s] START: Quality Control", Sys.time()))


data = anndata::read_h5ad(sprintf('%s/%s', data_path, data_file))
grna_count = data$obs$gene |> table()
counts = data$X
# colnames(counts) = data$var$gene_name
colnames(counts) = row.names(data$var)


print(sprintf("[%s]   - filter grnas and genes", Sys.time()))
# prep
## for each perturbation, how many cells received this perturbation
grna_count = data$obs$gene_id |> table() # grna_count[grna_count < 400] |> hist(xlim = c(0, 400), breaks = seq(from = 0, to = 400, by = 10)) |> abline(v=MIN_GRNA_COUNT, col='red')
## for each gene, how many cells have nonzero counts
n_nonzero_per_gene = apply(counts, MARGIN = 2, FUN = function(x) {sum(x == 0)}) # hist(n_nonzero_per_gene)  |> abline(v=MIN_N_NONZERO_PER_GENE, col='red')

# choose genes and perturbation
chosen_grnas = grna_count[MIN_GRNA_COUNT <= grna_count] |> names()        
chosen_grnas = chosen_grnas[chosen_grnas != 'nan'] # remove 'nan'               # 782 out of 2392 grnas
# chosen_grnas = c(chosen_grnas[1:2], "non-targeting" )

chosen_genes_mask= which(MIN_N_NONZERO_PER_GENE <= n_nonzero_per_gene)
chosen_genes_mask = union(chosen_genes_mask, which(colnames(counts) %in% chosen_grnas)) |> sort() # add genes w/ a grna target even if gene doesn't meet n_nonzero thresh
# chosen_grnas[which(!chosen_grnas %in% colnames(counts))]
# counts[1:4, "ENSG00000033327"]
# chosen_genes_mask = chosen_genes_mask[1:2]
# chosen_genes     = data$var[chosen_genes_mask, 'gene_name'] |> as.character() # 4627 out of 8749 genes
chosen_genes     = data$var[chosen_genes_mask, ] |> row.names() |> as.character() # 4627 out of 8749 genes




# save data in format for SCEPTRE:
# - response_mat.csv
# - grna_mat.csv
chosen_cells = (data$obs$gene_id %in% chosen_grnas) #  224719 out of 247914 cells

print(sprintf("[%s]   - construct response_mat and grna_mat", Sys.time()))
# SCEPTRE: response_matrix: genes as ros x cells as cols, roname should be gene names
response_mat = counts[chosen_cells, chosen_genes] |> t() # filter for genes and cells receiving chosen pert 6803 x 224719
# SCEPTRE: grna_matrix: mat of gRNA UMI counts, gRNAs as ros x cells in cols
grna_assignment = data$obs$gene_id[chosen_cells] |> as.character()
grna_mat = model.matrix(~0+grna_assignment)
row.names(grna_mat) = row.names(data$obs)[chosen_cells]
colnames(grna_mat) = gsub(pattern = "grna_assignment", replacement = "", x = colnames(grna_mat))
grna_mat =t(grna_mat) # 1610 x 224719

print(sprintf("[%s]   - save response_mat and grna_mat", Sys.time()))
# write.csv(response_mat[1:5, 1:5], sprintf('%s/sceptre_format/response_mat.csv', data_path))
# write.csv(    grna_mat, sprintf('%s/sceptre_format/grna_mat.csv',     data_path))
saveRDS(response_mat, sprintf('%s/sceptre/response_mat.rds', data_path))
saveRDS(    grna_mat, sprintf('%s/sceptre/grna_mat.rds',     data_path))


# ==============================================================================
#                              SCEPTRE
# ==============================================================================

print(sprintf("[%s] START: SCEPTRE", Sys.time()))

# START: SCEPTRE on Replogle rpe1 experiment
library(sceptre)
library(ondisc)
# library(Matrix)


print(sprintf("[%s]   - load and prep data", Sys.time()))
data = anndata::read_h5ad(sprintf('%s/%s', data_path, data_file))
response_mat = readRDS(sprintf('%s/sceptre/response_mat.rds', data_path))
grna_mat     = readRDS(sprintf('%s/sceptre/grna_mat.rds',     data_path))



# grna_target_data_frame: give no information about the grna... (data$var doesn't have all grna's info... data$obs doesn't have th desired info: chr, start, ...) could pull from another dataset but...
grna_target_data_frame = data.frame(grna_id     = row.names(grna_mat),
                                    grna_target = row.names(grna_mat))



# randomly distribute non-targeting to non-targeting1 - 5... this is what sceptre wants...

## modify grna_target_data_frame
grna_target_data_frame = rbind(grna_target_data_frame |> dplyr::filter(grna_id != 'non-targeting'), 
                               data.frame(grna_id = paste0('non-targeting', 1:5), 
                                          grna_target = 'non-targeting'))
## modify grna_mat- remove original 'non-targeting' and split nt cells into 5 nt's
set.seed(1235667)
nt_idx       = which(grna_mat['non-targeting', ] == 1) |> sample(replace=FALSE) # shuffle order
ncell_per_nt = length(nt_idx)/5

nt_mat = matrix(0, nrow = 5, ncol = ncol(grna_mat))
row.names(nt_mat) = paste0('non-targeting', 1:5)
nt_mat[1, nt_idx[                   1:(1*ncell_per_nt)]] = 1
nt_mat[2, nt_idx[(1*ncell_per_nt + 1):(2*ncell_per_nt)]] = 1
nt_mat[3, nt_idx[(2*ncell_per_nt + 1):(3*ncell_per_nt)]] = 1
nt_mat[4, nt_idx[(3*ncell_per_nt + 1):(4*ncell_per_nt)]] = 1
nt_mat[5, nt_idx[(4*ncell_per_nt + 1):length(nt_idx)]]   = 1

# dim(grna_mat)
# grna_mat[-which(row.names(grna_mat) == 'non-targeting'), ] |> dim()
grna_mat2 = rbind(grna_mat[-which(row.names(grna_mat) == 'non-targeting'), ] , 
                  nt_mat)  |> as.matrix()   
# row.names(grna_mat2)[grepl('non-targeting', row.names(grna_mat2))]


# extra_covariates: covariates in nb fit, data$obs has potential covariates: mitopercent, UMI_count, z_gemgroup_UMI, core_scale_factor, core_adjusted_UMI_count
assertthat::assert_that(all(colnames(response_mat) %in% row.names(data$obs))) # all chosen cell names are in data$obs df
extra_covariates = data$obs[colnames(response_mat), ] |> dplyr::select(mitopercent, core_adjusted_UMI_count) 



assertthat::assert_that(all(grna_target_data_frame$grna_id %in% row.names(grna_mat2)))



# check inputs

# dim(response_mat)
# row.names(response_mat)
# response_mat[1:2, 1:20]
# rowSums(response_mat)
# dim(grna_mat2)
# row.names(grna_mat2)
# rowSums(grna_mat2)
# head(grna_target_data_frame)
# head(extra_covariates)

print(sprintf("[%s]   - create SCEPTRE obj", Sys.time()))
# sceptre
sceptre_object <- sceptre::import_data(
  response_matrix        = response_mat,
  grna_matrix            = grna_mat2,
  grna_target_data_frame = grna_target_data_frame,
  moi                    = 'low',
  extra_covariates       = extra_covariates,
  response_names         = row.names(response_mat)
)

print(sprintf("[%s]   - construct positive and discovery pairs", Sys.time()))
# positive_control_pairs <- construct_positive_control_pairs(sceptre_object) # <- does not work for some reason
positive_control_pairs = data.frame(grna_target = row.names(grna_mat),
                                    response_id = row.names(grna_mat)) |> 
  dplyr::filter(grna_target != 'non-targeting') |>
  dplyr::filter(grna_target %in% row.names(response_mat))

# discovery_pairs <- construct_cis_pairs(
#   sceptre_object, 
#   positive_control_pairs = positive_control_pairs
# )

discovery_pairs<-construct_trans_pairs( # use trans instead of cis!
  sceptre_object,
  positive_control_pairs = positive_control_pairs
)


print(sprintf("[%s]       + subset pos and disc pairs for testing script", Sys.time()))
positive_control_pairs = positive_control_pairs |> dplyr::slice_sample(n=300)  # test
discovery_pairs        = discovery_pairs        |> dplyr::slice_sample(n=1000) 


print(sprintf("[%s]   - set analysis params, assign grnas, and run qc", Sys.time()))
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


print(sprintf("[%s]   - run analysis calibration (negatives)", Sys.time()))
sceptre_object = sceptre_object |>
  run_calibration_check(parallel=PARALLEL)

print(sprintf("[%s]   - run analysis power (positives)", Sys.time()))
sceptre_object = sceptre_object |>
  run_power_check(parallel=PARALLEL)

print(sprintf("[%s]   - run analysis discovery (discovery)", Sys.time()))
sceptre_object = sceptre_object  |>
  run_discovery_analysis(parallel=PARALLEL)





print(sprintf("[%s]   - write outputs", Sys.time()))
print(sprintf("[%s]       + sceptre auto outputs", Sys.time()))
# Write outputs to directory
dir.create('../saves/sceptre/replogle/', recursive=TRUE)
list.files('../saves/sceptre/replogle/outputs/')
write_outputs_to_directory(
  sceptre_object = sceptre_object, 
  directory = '../saves/sceptre/replogle/outputs/'
)

print(sprintf("[%s]       + positive control and discovery pairs", Sys.time()))
write.csv(positive_control_pairs, '../saves/sceptre/replogle/positive_control_pairs.csv', row.names = FALSE)
write.csv(discovery_pairs       , '../saves/sceptre/replogle/discovery_pairs.csv'       , row.names = FALSE)

print(sprintf("[%s]       + sceptre object", Sys.time()))
saveRDS(sceptre_object, file = '../saves/sceptre/replogle/sceptre_object.rds')


print(sprintf("[%s] END", Sys.time()))

# ==============================================================================




# ==============================================================================
#                              TRASH
# ==============================================================================

# 
# 
# # Trash
# # numbers diff bc use rpe1 instead of K562
# remotes::install_github("mojaveazure/seurat-disk")
# data_path = '../../../genData/replogle'
# data_file = 'K562_essential_raw_singlecell_01.h5ad'
# data_file = 'rpe1_raw_singlecell_01.h5ad'
# 
# 
# data = anndata::read_h5ad(sprintf('%s/%s', data_path, data_file))
# # data = SeuratDisk::LoadH5Seurat(sprintf('%s/%s', data_path, data_file))
# 
# 
# data$n_obs
# 
# # guessing: cells
# data$obs |> dim() # 310385     11
# data$obs |> head()
# 
# 
# grna_count = data$obs$gene |> table()
# sort(grna_count, decreasing = TRUE)
# grna_count |> hist()
# grna_count[grna_count < 400] |> hist(xlim = c(0, 400), breaks = seq(from = 0, to = 400, by = 10)) # 
# grna_count[75 < grna_count & grna_count < 400] |> hist(xlim = c(0, 400), breaks = seq(from = 0, to = 400, by = 10)) # 
# 
# # guessing: genes
# data$var |> dim() # 8563   12
# data$var |> head()
# 
# 
# 
# # counts: cells x genes
# data$X |> dim() # 310385   8563
# 
# counts = data$X 
# 
# # for each gene, how many cells have nonzero counts
# n_nonzero_per_gene = apply(counts, MARGIN = 2, FUN = function(x) {sum(x == 0)})
# # for each gene, avg gene expression over cells 
# avg_count_per_gene = apply(counts, MARGIN = 2, FUN = function(x) {mean(x)})
# 
# 
# hist(n_nonzero_per_gene)
# 
# 
# 
# 
# # Trash
# 
# 
# 
# # sceptre- Ergan- somehow this doesn't work...
# 
# print('Create sceptre Dataset...')
# sceptre_object <- import_data(
#   response_matrix = response_mat,
#   grna_matrix = grna_mat2,
#   grna_target_data_frame = grna_target_data_frame,
#   moi = 'low',
#   extra_covariates = extra_covariates,
#   use_ondisc = TRUE,
#   directory_to_write = disc_path
# )
# print('Object Created!')
# # positive_control_pairs <- construct_positive_control_pairs(sceptre_object) # <- does not work for some reason
# positive_control_pairs = data.frame(grna_target = row.names(grna_mat),
#                                     response_id = row.names(grna_mat)) |> 
#   dplyr::filter(grna_target != 'non-targeting') |>
#   dplyr::filter(grna_target %in% row.names(response_mat))
# print('Begin pairs!')
# # discovery_pairs <- construct_cis_pairs(
# #   sceptre_object, 
# #   positive_control_pairs = positive_control_pairs
# # )
# 
# discovery_pairs<-construct_trans_pairs( # use trans instead of cis!
#   sceptre_object,
#   positive_control_pairs = positive_control_pairs
# )
# 
# # apply the pipeline functions to the sceptre_object in order
# print('Apply the Pipeline...')
# sceptre_object <- sceptre_object |>
#   set_analysis_parameters(discovery_pairs, positive_control_pairs)
# 
# 
# 
# 
# sceptre_object |> print()
# 
# # options(error = recover)
# sceptre_object@cells_in_use
# browser()
# 
# sceptre_object |>
#   run_calibration_check()
# 
# 
# nt_idxs_new <- stats::setNames(lapply(grna_assignments_raw$indiv_nt_grna_idxs, 
#                                       function(v) {
#                                         update_idxs(v, cells_in_use, n_cells)
#                                       }), names(grna_assignments_raw$indiv_nt_grna_idxs))
# 
# nt_idxs_new[vapply(nt_idxs_new, length, FUN.VALUE = integer(1)) != 
#               0L]
# vapply(nt_idxs_new, length, FUN.VALUE = integer(1)) != 
#   0L
# 
# sceptre_object = sceptre_object
# run_power_check() |>
#   run_discovery_analysis()
# print('Begin Saving...')
# 
# write_outputs_to_directory(sceptre_object, disc_path)
# 
# 
#  
# grna_id = row.names(grna_mat)
# 
# (!(grna_id %in% row.names(data$var))) |> sum()
# (!(grna_id %in% data$obs$gene_id)) |> sum() # NOT all grna's have a corresponding gene that is measured????
# sum(! (data$obs |> dplyr::pull(gene) |> unique()) %in% data$var$gene_name)
# sum(! (data$obs |> dplyr::filter(gene_id %in% grna_id) |> dplyr::pull(gene) |> unique()) %in% data$var$gene_name)
# 
# data$obs |> dplyr::filter(gene_id %in% grna_id) |> dplyr::pull(gene)
# 
# 
# 
# 
# grna_target_data_frame = data$var[grna_id[1609:1610], ] |> dplyr::mutate(grna_id = grna_id, gene_id = grna_id, .before = 1)
# 
# 
# # change from gene names to gen id... gene names are not unique?????
# grna_id = row.names(grna_mat)
# data$obs |> head()
# data$var |> head()
# 
# gene_meta = data$var
# row.names(gene_meta) = gene_meta$gene_name
# data$var |> dim()
# data$var |> dplyr::group_by(gene_name) |> dplyr::summarize(count = dplyr::n()) |> dplyr::arrange(desc(count))
# # duplicate gene?? HSPA14
# data$var |> dplyr::filter(gene_name == 'HSPA14')
# 
# 'HSPA14' %in% grna_id
# grna_id[grna_id == 'HSPA14']
# 
# # from Ergan's SCEPTRE pipeline- other data is saved in different format
# 
# 
# 
# # first column as the colnames: row.names=1
# response_matrix<-Matrix(as.matrix(read.csv(paste0(data_path, 'response_mat.csv'), row.names = 1, check.names = FALSE), sparse=True))
# 
# 
# 
# # rownames like ENSG000....
# print(rownames(response_matrix)[1:10])
# print(colnames(response_matrix)[1:10])
# 
# grna_matrix<-Matrix(as.matrix(read.csv(paste0(data_path, 'grna_mat.csv'), row.names=1, check.names=FALSE), sparse=True))
# 
# # rownames like NOCL2, non-targeting_03751|non-targeting_00758....
# print(rownames(grna_matrix)[1:10])
# print(colnames(grna_matrix)[1:10])
# grna_id<-rownames(grna_matrix)
# grna_target<-ifelse(grepl("^non-targeting", grna_id),
#                     "non-targeting",
#                     gene_corresp$gene_id[match(grna_id, gene_corresp$gene_name)]
# )
# grna_chr <- ifelse(
#   grepl("^non-targeting", grna_id),
#   NA,
#   gene_corresp$chr[match(grna_id, gene_corresp$gene_name)]
# )
# grna_start <- ifelse(
#   grepl("^non-targeting", grna_id),
#   NA,
#   gene_corresp$start[match(grna_id, gene_corresp$gene_name)]
# )
# grna_end <- ifelse(
#   grepl("^non-targeting", grna_id),
#   NA,
#   gene_corresp$end[match(grna_id, gene_corresp$gene_name)]
# )
# 
# 
# meta<-data.frame(grna_id=grna_id, grna_target=grna_target, chr=grna_chr, start=grna_start, end=grna_end)
# print(meta[c(1:3, 2216:2219), ])
# rownames(meta)<-NULL
# 
# print('Create sceptre Dataset...')
# sceptre_object <- import_data(
#   response_matrix = response_matrix,
#   grna_matrix = grna_matrix,
#   grna_target_data_frame = meta,
#   moi = 'low',
#   use_ondisc = TRUE,
#   directory_to_write = disc_path
# )
# print('Object Created!')
# positive_control_pairs <- construct_positive_control_pairs(sceptre_object)
# print('Begin pairs!')
# # discovery_pairs <- construct_cis_pairs(
# #   sceptre_object, 
# #   positive_control_pairs = positive_control_pairs
# # )
# 
# discovery_pairs<-construct_trans_pairs( # use trans instead of cis!
#   sceptre_object,
#   positive_control_pairs = positive_control_pairs
# )
# 
# # apply the pipeline functions to the sceptre_object in order
# print('Apply the Pipeline...')
# sceptre_object <- sceptre_object |> # |> is R's base pipe, similar to %>%
#   set_analysis_parameters(discovery_pairs, positive_control_pairs) |>
#   run_calibration_check() |>
#   run_power_check() |>
#   run_discovery_analysis()
# print('Begin Saving...')
# save_folder<-'/home/eshang/diffusion_and_protein/repogle/data/Tim_res'
# write_outputs_to_directory(sceptre_object, save_folder)
# 
# 
# 
# 
# 
# 
# 
