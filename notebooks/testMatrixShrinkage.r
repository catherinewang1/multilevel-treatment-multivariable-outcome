# Test ebci shrinkage to random values
# (util functions should be in utils/matrix_shrinkage.r)
# (saved files should be in saves/)


print(sprintf("[%s] START: Test ebci shrinkage to random values", Sys.time()))




# === Settings/Parameters ===
ALPHA = .1



# === Load and Prep Data ===

print(sprintf("[%s] Load and Prep Data", Sys.time()))


## load
print(sprintf("[%s]     - load", Sys.time()))

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(ebci)) # robust emp shrinkage
suppressPackageStartupMessages(library(assertthat)) 


source('../utils/matrix_shrinkage.r')





dir.create('../plots/matrixSampSplit/shrinkMatrixTest/') # test shrinkage




# load estimated effects
effects_list = list()
for(est_method in c('sceptre')) {
  effects_list[[est_method]] = list()
  for(testsplit in c('train', 'test')) { 
    
    if(est_method == 'sceptre') {path = 'sceptre/sceptre'} 
    else if(est_method == 'nb') {path = 'negativebinomial'}
    else {path = est_method}
    
    effects_list[[est_method]][[testsplit]] = 
      read.csv(sprintf('../saves/%s_effects_%s.csv',
                       path, testsplit)) |> 
      # filter(!is.na(estimate)) |>   # remove tests with NA estimates
      mutate(tstat = estimate / se) |> # make a tstat column = estimate / se
      mutate(se_old = se, 
             se = mapply(FUN = spline_se, 
                         p   = pvalue, 
                         mu  = estimate) )

    
  }
}

# compare se=est/z to smooth se
ggplot(effects_list$sceptre$train |> slice_sample(n=1000), 
       aes(x = se_old, y = se)) + 
       geom_abline(aes(slope = 1, intercept = 0)) +
       scale_x_log10() +
       scale_y_log10() +
       geom_point()


## prep
print(sprintf("[%s]     - prep", Sys.time()))
### important genes and grnas selected by a 'score'
effects_df = effects_list[['sceptre']][['train']] |> filter(test == 'discovery') # | test == 'positive') 
topgenes = effects_df |> group_by(gene) |> summarize(sum_effects = sum(abs(estimate), na.rm = TRUE),
                                                     sd_effects = sd(estimate, na.rm = TRUE),
                                                     sig_effects = sum(significant, na.rm=TRUE),
                                                     sum_tstats  = sum(abs(estimate/se), na.rm = TRUE)) |>
  mutate(score = scale(sum_effects) + (1/3)* scale(sd_effects) +  sig_effects) |> 
  arrange(desc(score)) |> 
  arrange(desc(sum_tstats))
topgrnas = effects_df |> group_by(grna) |> summarize(sum_effects = sum(abs(estimate), na.rm = TRUE),
                                                     sd_effects = sd(estimate, na.rm = TRUE),
                                                     sig_effects = sum(significant, na.rm=TRUE),
                                                     sum_tstats  = sum(abs(estimate/se), na.rm = TRUE)) |>
  mutate(score = scale(sum_effects) + (0)* scale(sd_effects) +  sig_effects) |> 
  arrange(desc(score)) |> 
  arrange(desc(sum_tstats))


### Construct matrices of estimate and se (and others, but not used, plots some in the process)
# construct matrices (explore toptstat, but use topscore)

# score
effects_df = effects_df |> filter((gene %in% (topgenes |> arrange(desc(score)) |> head(2000) |> pull(gene)))
                                  & 
                                    (grna %in% (topgrnas |> arrange(desc(score)) |> head(50) |> pull(grna))))


matrices = make_matrices(effects_df = effects_df) #, save_plots_filepath = '../plots/matrixSampSplit/topGeneGrna/topscoreTrain/')

# test (using train's top scores)
chosengenes = effects_df |> pull(gene) |> unique() |> sort()
chosengrnas = effects_df |> pull(grna) |> unique() |> sort()

effects_df_test = effects_list[['sceptre']][['test']] |> 
  dplyr::filter(test == 'discovery' & (gene %in% chosengenes) & (grna %in% chosengrnas))

matrices_test = make_matrices(effects_df = effects_df_test) #, save_plots_filepath = '../plots/matrixSampSplit/topGeneGrna/topscoreTest/')

# make sure order of ros and cols the same (it should be ok w/o, but to be safe)
# for now, just do alphabetical order (chosengrnas and chosengenes)
matrices      =  mapply(FUN = my_reorder_rc, mat = matrices,
                        MoreArgs = list(rorder = chosengrnas, corder = chosengenes),
                        SIMPLIFY=FALSE)
matrices_test =  mapply(FUN = my_reorder_rc, mat = matrices_test,
                        MoreArgs = list(rorder = chosengrnas, corder = chosengenes),
                        SIMPLIFY=FALSE)



# stop('Stop after making matrices and plots comparing toptstat vs topscore')

# max se set to more reasonable value (not 999)... idk why this takes long...
se_mat_test = apply(X = matrices_test$se, MARGIN=2,
                    FUN = function(x){pmin(x, quantile(matrices_test$se |> unlist(), .9))})


assertthat::assert_that(quantile(matrices$se |> unlist(), .9) < 999, 
                        msg = "Fix: set limit for se to be < 999 (e.g. set quantile smaller)")



# choose some ordering for the rows and columns (mainly for fixed matrix s.t. the blocks exist)
matscaled = as.matrix(scale(matrices$estimates))
matscaled[is.nan(matscaled)] = 0
row_order = hclust(dist(matscaled))$order
column_order = hclust(dist(t(matscaled)))$order
grna_index = data.frame(grna = rownames(matscaled)[row_order],
                        grna_idx  = 1:nrow(matscaled))
gene_index = data.frame(gene = colnames(matscaled)[column_order],
                        gene_idx  = 1:ncol(matscaled))



# ==============================================================================
# Start: EBCI shrinkage 
# - sim data w random theta
# - sim data w fixed theta
# - sim data w fixed theta + incorrect se's (1.25x, 1.5x, 2x)
# - real data to fixed matrix
# - real data to dim reduced matrix (low rank)
# - real data to dim reduced matrix (sparse)
# ==============================================================================
## perform ebci shrinkage to fixed matrix --------------------------------------
nrow = dim(matrices_test$estimates)[1]
ncol = dim(matrices_test$estimates)[2]


fixed_matrix = c(rep(1, 20 ), rep(0, nrow-20 )) %*% t(c(rep(1, 750), rep(0, ncol-750))) * .5 +
               c(rep(0, 20 ), rep(1, 10), rep(0, nrow-30 )) %*% t(c(rep(0, 750), rep(1,400), rep(0, ncol-1150))) * (-.25)
colnames(fixed_matrix) = colnames(se_mat_test[row_order, column_order])
rownames(fixed_matrix) = rownames(se_mat_test[row_order, column_order])

dim(fixed_matrix); image(fixed_matrix)

shrinkFixed = 
  shrink_matrix(   unshrunk_mat = matrices_test$estimates[row_order, column_order],
                   shrinkpoint_mat = fixed_matrix,
                   se_mat = se_mat_test[row_order, column_order],
                   weight_mat = 1/se_mat_test[row_order, column_order]**4,
                   ALPHA = ALPHA, 
                   return_ebci_obj = TRUE) # ~15 mins


write.csv(x = shrinkFixed$ebci_res,
          file = sprintf("../saves/matrixSampSplit/shrinkMatrixTest/shrinkMatrix_fixed.csv"))

## perform ebci shrinkage of sim est w/ random theta to fixed matrix --------------------------------------
# random theta (parametric case), should be ideal case

# 30 x 1000
set.seed(12345)

true_theta    = matrix(rnorm(30*1000, mean = 0, sd = 1), nrow = 30, ncol = 1000)
true_theta = as.vector(true_theta)

shrinkage_points = matrix(0, nrow = 30, ncol = 1000)
shrinkage_points[1:10, 1:200] = .7
shrinkage_points[11:15, 201:500] = -.1
shrinkage_points[1:3, 800:815] = .3
shrinkage_points[20:30, 801:1000] = .9
shrinkage_points[20:30, 601:800] = .1
shrinkage_points = as.vector(shrinkage_points)

se_vec = se_mat_test |> as.vector(); se_vec = se_vec[se_vec<1.5] 
sim_se = sample(se_vec, size = length(as.vector(true_theta)))

sim_estimates = mapply(FUN = rnorm, n = 1, mean = true_theta, sd = sim_se)

sim_rownames = sprintf('grna%02.f', 1:30)
sim_colnames = sprintf('gene%03.f', 1:1000)
sim_estimates    = matrix(sim_estimates,    nrow = 30, dimnames = list(sim_rownames, sim_colnames)) 
true_theta       = matrix(true_theta,       nrow = 30, dimnames = list(sim_rownames, sim_colnames)) 
shrinkage_points = matrix(shrinkage_points, nrow = 30, dimnames = list(sim_rownames, sim_colnames)) 
sim_se           = matrix(sim_se,           nrow = 30, dimnames = list(sim_rownames, sim_colnames)) 



true_theta_df = reshape2::melt(true_theta, varnames = c("grna", "gene"), value.name = 'true_effect')

image(true_theta); image(shrinkage_points); image(sim_se); image(sim_estimates)


# use actual se's (100%)
shrinkSimEst = 
  shrink_matrix(   unshrunk_mat    = sim_estimates,
                   shrinkpoint_mat = shrinkage_points,
                   se_mat          = sim_se,
                   weight_mat      = 1/sim_se**4,
                   ALPHA = ALPHA, 
                   return_ebci_obj = TRUE) # ~15 mins

write.csv(x = merge(shrinkSimEst$ebci_res, true_theta_df, by = c('grna', 'gene')), # add true effect
          file = sprintf("../saves/matrixSampSplit/shrinkMatrixTest/shrinkMatrix_sim_randtheta.csv"))



# random theta, mixture of 0's and other norm(0,1)
# basically a mixture of: zero-inflated distribution + N(0,1) + N(-.5, 1) + N(.6, 1)
set.seed(12345)
true_theta    = matrix(0, nrow = 30, ncol = 1000) # mostly 0's
true_theta[ 1:10,    1:400] = matrix(rnorm(10*400, mean = 0, sd = 1), nrow = 10, ncol = 400)
true_theta[11:15,  401:600] = matrix(rnorm(5*200, mean = -.5, sd = 1), nrow = 5, ncol = 200)
true_theta[16:25,  601:900] = matrix(rnorm(10*300, mean = .6, sd = 1), nrow = 10, ncol = 300)
# image(true_theta)
true_theta = as.vector(true_theta)

shrinkage_points = matrix(0, nrow = 30, ncol = 1000)
shrinkage_points[1:10, 1:200] = .7
shrinkage_points[11:15, 201:500] = -.1
shrinkage_points[1:3, 800:815] = .3
shrinkage_points[20:30, 801:1000] = .9
shrinkage_points[20:30, 601:800] = .1
shrinkage_points = as.vector(shrinkage_points)

se_vec = se_mat_test |> as.vector(); se_vec = se_vec[se_vec<1.5] 
sim_se = sample(se_vec, size = length(as.vector(true_theta)))

sim_estimates = mapply(FUN = rnorm, n = 1, mean = true_theta, sd = sim_se)

sim_rownames = sprintf('grna%02.f', 1:30)
sim_colnames = sprintf('gene%03.f', 1:1000)
sim_estimates    = matrix(sim_estimates,    nrow = 30, dimnames = list(sim_rownames, sim_colnames)) 
true_theta       = matrix(true_theta,       nrow = 30, dimnames = list(sim_rownames, sim_colnames)) 
shrinkage_points = matrix(shrinkage_points, nrow = 30, dimnames = list(sim_rownames, sim_colnames)) 
sim_se           = matrix(sim_se,           nrow = 30, dimnames = list(sim_rownames, sim_colnames)) 



true_theta_df = reshape2::melt(true_theta, varnames = c("grna", "gene"), value.name = 'true_effect')

image(true_theta); image(shrinkage_points); image(sim_se); image(sim_estimates)


# use actual se's (100%)
shrinkSimEst = 
  shrink_matrix(   unshrunk_mat    = sim_estimates,
                   shrinkpoint_mat = shrinkage_points,
                   se_mat          = sim_se,
                   weight_mat      = 1/sim_se**4,
                   ALPHA = ALPHA, 
                   return_ebci_obj = TRUE) # ~15 mins

write.csv(x = merge(shrinkSimEst$ebci_res, true_theta_df, by = c('grna', 'gene')), # add true effect
          file = sprintf("../saves/matrixSampSplit/shrinkMatrixTest/shrinkMatrix_sim_randtheta0.csv"))



## perform ebci shrinkage of sim est to fixed matrix --------------------------------------


# simulate data that is similar to the actual estimated se's


# 30 x 1000
set.seed(12345)

true_theta    = c(rep(1, 10 ), rep(0, 20 )) %*% t(c(rep(1, 200), rep(0, 800))) * .5 +
                c(rep(0, 10 ), rep(1, 10), rep(0, 10 )) %*% t(c(rep(0, 200), rep(1,500), rep(0, 300))) * (-.25)

true_theta[1:3, 800:815] = .75
true_theta[21:25, 700:800] = 1
true_theta = as.vector(true_theta)

shrinkage_points = matrix(0, nrow = 30, ncol = 1000)
shrinkage_points[1:10, 1:200] = .7
shrinkage_points[11:15, 201:500] = -.1
shrinkage_points[1:3, 800:815] = .3
shrinkage_points[20:30, 801:1000] = .9
shrinkage_points[20:30, 601:800] = .1
shrinkage_points = as.vector(shrinkage_points)

se_vec = se_mat_test |> as.vector(); se_vec = se_vec[se_vec<1.5] 
sim_se = sample(se_vec, size = length(as.vector(true_theta)))

sim_estimates = mapply(FUN = rnorm, n = 1, mean = true_theta, sd = sim_se)

sim_rownames = sprintf('grna%02.f', 1:30)
sim_colnames = sprintf('gene%03.f', 1:1000)
sim_estimates    = matrix(sim_estimates,    nrow = 30, dimnames = list(sim_rownames, sim_colnames)) 
true_theta       = matrix(true_theta,       nrow = 30, dimnames = list(sim_rownames, sim_colnames)) 
shrinkage_points = matrix(shrinkage_points, nrow = 30, dimnames = list(sim_rownames, sim_colnames)) 
sim_se           = matrix(sim_se,           nrow = 30, dimnames = list(sim_rownames, sim_colnames)) 



true_theta_df = reshape2::melt(true_theta, varnames = c("grna", "gene"), value.name = 'true_effect')

image(true_theta); image(shrinkage_points); image(sim_se); image(sim_estimates)


# use actual se's (100%)
shrinkSimEst = 
  shrink_matrix(   unshrunk_mat    = sim_estimates,
                   shrinkpoint_mat = shrinkage_points,
                   se_mat          = sim_se,
                   weight_mat      = 1/sim_se**4,
                   ALPHA = ALPHA, 
                   return_ebci_obj = TRUE) # ~15 mins

write.csv(x = merge(shrinkSimEst$ebci_res, true_theta_df, by = c('grna', 'gene')), # add true effect
          file = sprintf("../saves/matrixSampSplit/shrinkMatrixTest/shrinkMatrix_simestimates.csv"))


# increase the se estimates by 1.25x... 125% 
shrinkSimEst = 
  shrink_matrix(   unshrunk_mat    = sim_estimates,
                   shrinkpoint_mat = shrinkage_points,
                   se_mat          = (sim_se*1.25),
                   weight_mat      = 1/(sim_se*1.25)**4,
                   ALPHA = ALPHA, 
                   return_ebci_obj = TRUE) # ~15 mins

write.csv(x = merge(shrinkSimEst$ebci_res, true_theta_df, by = c('grna', 'gene')), # add true effect
          file = sprintf("../saves/matrixSampSplit/shrinkMatrixTest/shrinkMatrix_simest_125.csv"))

# increase the se estimates by 1.50x... 150% 
shrinkSimEst = 
  shrink_matrix(   unshrunk_mat    = sim_estimates,
                   shrinkpoint_mat = shrinkage_points,
                   se_mat          = (sim_se*1.50),
                   weight_mat      = 1/(sim_se*1.50)**4,
                   ALPHA = ALPHA, 
                   return_ebci_obj = TRUE) # ~15 mins

write.csv(x = merge(shrinkSimEst$ebci_res, true_theta_df, by = c('grna', 'gene')), # add true effect
          file = sprintf("../saves/matrixSampSplit/shrinkMatrixTest/shrinkMatrix_simest_150.csv"))

# increase the se estimates by 2x... 200% 
shrinkSimEst = 
  shrink_matrix(   unshrunk_mat    = sim_estimates,
                   shrinkpoint_mat = shrinkage_points,
                   se_mat          = (sim_se*2),
                   weight_mat      = 1/(sim_se*2)**4,
                   ALPHA = ALPHA, 
                   return_ebci_obj = TRUE) # ~15 mins

write.csv(x = merge(shrinkSimEst$ebci_res, true_theta_df, by = c('grna', 'gene')), # add true effect
          file = sprintf("../saves/matrixSampSplit/shrinkMatrixTest/shrinkMatrix_simest_200.csv"))

# decrease the se estimates by 0.75x... 75% 
shrinkSimEst = 
  shrink_matrix(   unshrunk_mat    = sim_estimates,
                   shrinkpoint_mat = shrinkage_points,
                   se_mat          = (sim_se*.75),
                   weight_mat      = 1/(sim_se*.75)**4,
                   ALPHA = ALPHA, 
                   return_ebci_obj = TRUE) # ~15 mins

write.csv(x = merge(shrinkSimEst$ebci_res, true_theta_df, by = c('grna', 'gene')), # add true effect
          file = sprintf("../saves/matrixSampSplit/shrinkMatrixTest/shrinkMatrix_simest_075.csv"))

# decrease the se estimates by 0.5x... 50% 
shrinkSimEst = 
  shrink_matrix(   unshrunk_mat    = sim_estimates,
                   shrinkpoint_mat = shrinkage_points,
                   se_mat          = (sim_se*.50),
                   weight_mat      = 1/(sim_se*.50)**4,
                   ALPHA = ALPHA, 
                   return_ebci_obj = TRUE) # ~15 mins

write.csv(x = merge(shrinkSimEst$ebci_res, true_theta_df, by = c('grna', 'gene')), # add true effect
          file = sprintf("../saves/matrixSampSplit/shrinkMatrixTest/shrinkMatrix_simest_050.csv"))


## perform ebci shrinkage to random matrix point -------------------------------
print(sprintf("[%s]     - perform ebci shrinkage to random", Sys.time()))

nrow = dim(matrices_test$estimates)[1]
ncol = dim(matrices_test$estimates)[2]
set.seed(123456)
random_matrix = matrix(rnorm(n = nrow*ncol), nrow = nrow, ncol = ncol)
colnames(random_matrix) = colnames(se_mat_test)
rownames(random_matrix) = rownames(se_mat_test)

shrinkRandom = 
  shrink_matrix(   unshrunk_mat = matrices_test$estimates,
                   shrinkpoint_mat = random_matrix,
                   se_mat = se_mat_test,
                   weight_mat = 1/se_mat_test[c(1,2), ]**4,
                   ALPHA = ALPHA, 
                   return_ebci_obj = TRUE) # ~15 mins



# shrinkDFPlots = shrink_matrix_plots(df = shrinkRandom$ebci_res)
# 
# pdf(sprintf("../plots/matrixSampSplit/shrinkMatrix/shrinkMatrix__r=%02.f.pdf",approxMethod, rank))
# gridExtra::grid.arrange(grobs = shrinkDFPlots)
# dev.off()

write.csv(x = shrinkRandom$ebci_res,
          file = sprintf("../saves/matrixSampSplit/shrinkMatrixTest/shrinkMatrix_random.csv"))





## perform ebci shrinkage on subset (2, 10, 25 perturbations, many genes) --------------
## - how does shrinkage change with increasing number of tests
## - All whould be 50 perturbations 

## perform ebci shrinkage on subset (2 perturbations, many genes) --------------
print(sprintf("[%s]     - perform ebci shrinkage to 2 perturbations", Sys.time()))



approxmatrices_sparseSVD = 
  approx_matrix(mat = matrices$estimates, 
                method = 'sparseSVD', 
                ranks = c(10),
                methodParams=list(type = 'standard', 
                                  sumabs = .35,  # should be between 0-1, 
                                  # sumabsu = 4, sumabsv = 4, # between 1 and sqrt(#col or #rows)
                                  niter = 100,
                                  trace = FALSE)) 

shrink2Grnas = 
  shrink_matrix(   unshrunk_mat = matrices_test$estimates[c(1,2), ],
                   shrinkpoint_mat = approxmatrices_sparseSVD$approxmatrices[[10]][c(1,2), ],
                   se_mat = se_mat_test[c(1,2), ],
                   weight_mat = (1/se_mat_test[c(1,2), ])**4,
                   ALPHA = ALPHA, 
                   return_ebci_obj = TRUE) # ~  mins



write.csv(x = shrink2Grnas$ebci_res, file = sprintf("../saves/matrixSampSplit/shrinkMatrixTest/shrinkMatrix_2grna.csv"))




## perform ebci shrinkage on subset (10 perturbations, many genes) --------------
print(sprintf("[%s]     - perform ebci shrinkage to 10 perturbations", Sys.time()))



approxmatrices_sparseSVD = 
  approx_matrix(mat = matrices$estimates, 
                method = 'sparseSVD', 
                ranks = c(10),
                methodParams=list(type = 'standard', 
                                  sumabs = .35,  # should be between 0-1, 
                                  # sumabsu = 4, sumabsv = 4, # between 1 and sqrt(#col or #rows)
                                  niter = 100,
                                  trace = FALSE)) 

shrink2Grnas = 
  shrink_matrix(   unshrunk_mat = matrices_test$estimates[1:10, ],
                   shrinkpoint_mat = approxmatrices_sparseSVD$approxmatrices[[10]][1:10, ],
                   se_mat = se_mat_test[1:10, ],
                   weight_mat = (1/se_mat_test[1:10, ])**4,
                   ALPHA = ALPHA, 
                   return_ebci_obj = TRUE) # ~  mins



write.csv(x = shrink2Grnas$ebci_res, file = sprintf("../saves/matrixSampSplit/shrinkMatrixTest/shrinkMatrix_10grna.csv"))


## perform ebci shrinkage on subset (25 perturbations, many genes) --------------
print(sprintf("[%s]     - perform ebci shrinkage to 25 perturbations", Sys.time()))



approxmatrices_sparseSVD = 
  approx_matrix(mat = matrices$estimates, 
                method = 'sparseSVD', 
                ranks = c(10),
                methodParams=list(type = 'standard', 
                                  sumabs = .35,  # should be between 0-1, 
                                  # sumabsu = 4, sumabsv = 4, # between 1 and sqrt(#col or #rows)
                                  niter = 100,
                                  trace = FALSE)) 

shrink2Grnas = 
  shrink_matrix(   unshrunk_mat = matrices_test$estimates[1:25, ],
                   shrinkpoint_mat = approxmatrices_sparseSVD$approxmatrices[[10]][1:25, ],
                   se_mat = se_mat_test[1:25, ],
                   weight_mat = (1/se_mat_test[1:25, ])**4,
                   ALPHA = ALPHA, 
                   return_ebci_obj = TRUE) # ~  mins



write.csv(x = shrink2Grnas$ebci_res, file = sprintf("../saves/matrixSampSplit/shrinkMatrixTest/shrinkMatrix_25grna.csv"))


## perform ebci shrinkage on subset (50 perturbations, many genes) --------------
print(sprintf("[%s]     - perform ebci shrinkage to 50 perturbations", Sys.time()))



approxmatrices_sparseSVD = 
  approx_matrix(mat = matrices$estimates, 
                method = 'sparseSVD', 
                ranks = c(10),
                methodParams=list(type = 'standard', 
                                  sumabs = .35,  # should be between 0-1, 
                                  # sumabsu = 4, sumabsv = 4, # between 1 and sqrt(#col or #rows)
                                  niter = 100,
                                  trace = FALSE)) 

shrink2Grnas = 
  shrink_matrix(   unshrunk_mat = matrices_test$estimates[1:50, ],
                   shrinkpoint_mat = approxmatrices_sparseSVD$approxmatrices[[10]][1:50, ],
                   se_mat = se_mat_test[1:50, ],
                   weight_mat = (1/se_mat_test[1:50, ])**4,
                   ALPHA = ALPHA, 
                   return_ebci_obj = TRUE) # ~  mins



write.csv(x = shrink2Grnas$ebci_res, file = sprintf("../saves/matrixSampSplit/shrinkMatrixTest/shrinkMatrix_50grna.csv"))


# LOW RANK
## perform ebci shrinkage on subset (50 perturbations, many genes) --------------
print(sprintf("[%s]     - perform ebci shrinkage to 50 perturbations", Sys.time()))



approxmatrices_lowrank = 
  approx_matrix(mat = matrices$estimates, 
                method = 'lowrank', 
                ranks = c(10)) 

shrink2Grnas = 
  shrink_matrix(   unshrunk_mat = matrices_test$estimates[1:50, ],
                   shrinkpoint_mat = approxmatrices_sparseSVD$approxmatrices[[10]][1:50, ],
                   se_mat = se_mat_test[1:50, ],
                   weight_mat = (1/se_mat_test[1:50, ])**4,
                   ALPHA = ALPHA, 
                   return_ebci_obj = TRUE) # ~  mins



write.csv(x = shrink2Grnas$ebci_res, file = sprintf("../saves/matrixSampSplit/shrinkMatrixTest/shrinkMatrix_50grna_lowrank.csv"))


print(sprintf("[%s] END", Sys.time()))

# === === === === === === === === === === === === === === === === === === === === === === === === === === === === === === === === === === === === ===
# === === === === === === === === === === === === === === === === === === === === === === === === === === === === === === === === === === === === ===
# === === === === === === === === === === === === === === === === === === === === === === === === === === === === === === === === === === === === ===


shrink2Grnas$ebci_obj$mu2
shrink2Grnas$ebci_obj$kappa
shrink2Grnas$ebci_obj$df |> colnames()




shrinkRandom$ebci_obj$mu2
shrink2Grnas$ebci_obj$mu2

shrinkRandom$ebci_obj$kappa
shrink2Grnas$ebci_obj$kappa

shrinkRandom$ebci_obj$delta
shrink2Grnas$ebci_obj$delta

shrinkRandom$ebci_obj$delta
shrink2Grnas$ebci_obj$delta



# check the point estimate against their formula

# w_eb_i = mu_2 /(mu_2 + sigma_i^2)

# === random 
i = 103
mu2 = shrinkRandom$ebci_obj$mu2[1]
sigma_i = shrinkRandom$ebci_obj$df[i, 'se']
th_us = shrinkRandom$ebci_obj$df[i, 'th_us']
th_eb = shrinkRandom$ebci_obj$df[i, 'th_eb']
w_eb  = shrinkRandom$ebci_obj$df[i, 'w_eb']


print(paste0(w_eb, '  ', mu2 / (mu2 + sigma_i^2))) # should be =
print(paste0(th_eb, '  ', w_eb * th_us)) # should be equal =

ggplot(shrinkRandom$ebci_obj$df,
       aes(x = w_eb)) + geom_histogram()

ggplot(shrinkRandom$ebci_obj$df,
       aes(x = se)) + geom_histogram()

# === 2 grnas
i = 103
mu2 = shrink2Grnas$ebci_obj$mu2[1]
sigma_i = shrink2Grnas$ebci_obj$df[i, 'se']
th_us = shrink2Grnas$ebci_obj$df[i, 'th_us']
th_eb = shrink2Grnas$ebci_obj$df[i, 'th_eb']
w_eb  = shrink2Grnas$ebci_obj$df[i, 'w_eb']


print(paste0(w_eb, '  ', mu2 / (mu2 + sigma_i^2))) # should be =
print(paste0(th_eb, '  ', w_eb * th_us)) # should be equal =



ggplot(shrink2Grnas$ebci_obj$df,
       aes(x = th_us, y = th_eb)) + geom_point()
ggplot(shrink2Grnas$ebci_obj$df,
       aes(x = w_eb)) + geom_histogram()

ggplot(shrink2Grnas$ebci_obj$df,
       aes(x = se)) + geom_histogram()


# the se's should be the same... so the only change when calculating 
# w_eb should be mu2... somehow the mu2 is so much smaller in 2grnas...
# .39 vs .0023 (even smaller with spline_se 0.0007741335...)
# in the estimate for mu2= max(...\eps ..., ...) 
# the second term doesn't change with shrinkage points
# the first term changes with \hat\eps = Y - shrinkage point
# 

shrinkRandom$ebci_obj$mu2

# calculate the second term
n = shrinkRandom$ebci_obj$df |> nrow()
2 * sum(shrinkRandom$ebci_obj$df[, 'se']**4) / sum(shrinkRandom$ebci_obj$df[, 'se']**2)
2 * sum(shrinkRandom$ebci_obj$df[, 'se']**2) / sum(shrinkRandom$ebci_obj$df[, 'se']**1)
2 * sum(shrinkRandom$ebci_obj$df[, 'se']**0) / sum(shrinkRandom$ebci_obj$df[, 'se']**0)
2 * sum((1/n)**2 * shrinkRandom$ebci_obj$df[, 'se']**4) / sum((1/n) * shrinkRandom$ebci_obj$df[, 'se']**2)

?ebci::cva


ggplot(shrinkRandom$ebci_obj$df,
        aes(x = th_us - 0)) + geom_histogram()

ggplot(shrink2Grnas$ebci_obj$df,
       aes(x = th_us**2, y = se**2)) + geom_point()


mean(shrink2Grnas$ebci_obj$df$th_us**2)

shrink2Grnas$ebci_obj$mu2
shrink2Grnas$ebci_obj$kappa

#' #' estimate mu2 by hand
#' #' we dont have any X's
#' #' @param Y (vector) of numeric values, the Y_i's 
#' #'     in our case, this is unshrunk value - shrinkage point
#' #'     bc formula = unshrunk val - shrinkage point - 0
#' #'     which may be incorrect!
#' #' @param se (vector) of numeric values, the standard errors_i's
#' #'     this is \hat{\sigma}_i  in the formulas
#' #' @param weights (vector) of numeric values, the weights to put on each sample
#' #' @param X (matrix) of numeric values, the covariates
#' #' @param delta (vector) of coefs for covariates  
#' #' @returns estimated mu2 (numeric value)
#' #' @example
#' #' # this actually gives the same... 0.0007741335
#' #' my_est_mu2(Y = unlist(matrices_test$estimates[c(1,2), ] - 
#' #'                         approxmatrices_sparseSVD$approxmatrices[[10]][c(1,2), ]),
#' #'            se = as.vector(se_mat_test[c(1,2), ]),
#' #'            weights = rep(1/length(unlist(se_mat_test[c(1,2), ])), 
#' #'                          times = length(unlist(se_mat_test[c(1,2), ]))))
#' #' 
#' my_est_mu2 <- function(Y, se, weights, X=NULL, delta=NULL) {
#'   # Y = matrices_test$estimates[1:50, 1] - approxmatrices_sparseSVD$approxmatrices[[10]][1:50, 1]
#'   # se = se_mat_test[1:50, 1]
#'   # weights = rep(1/length(se), times = length(se))
#'   # X=NULL; delta=NULL
#'   
#'   if(is.null(X) | is.null(delta)) {
#'     eps = Y
#'   } else {
#'     eps = Y - X %*% delta
#'   }
#'   
#'   
#'   first = (sum(weights * (eps*eps - se*se))
#'              /
#'            sum(weights))
#'     
#'   
#'   second = 2 * (sum(weights * weights * (se^4))
#'                   /
#'                (sum(weights) * sum(weights * se * se)))
#'   
#'   # print(sprintf('first = %.6f, second = %.6f', first, second))
#'   return(max(first, second))
#' }
#' 
#' 
#' #' Estimate kappa by hand
#' #' @param Y (vector) of numeric values, the Y_i's 
#' #'     in our case, this is unshrunk value - shrinkage point
#' #'     bc formula = unshrunk val - shrinkage point - 0
#' #'     which may be incorrect!
#' #' @param se (vector) of numeric values, the standard errors_i's
#' #'     this is \hat{\sigma}_i  in the formulas
#' #' @param weights (vector) of numeric values, the weights to put on each sample
#' #' @param mu2 (numeric) estimate of mu2 (from a prev step)
#' #' @param X (matrix) of numeric values, the covariates
#' #' @param delta (vector) of coefs for covariates  
#' #' @example
#' #' my_est_kappa(Y = matrices_test$estimates[1:50, 1] - approxmatrices_sparseSVD$approxmatrices[[10]][1:50, 1],
#' #'              se = se_mat_test[1:50, 1],
#' #'              mu2 = shrink2Grnas$ebci_obj$mu2['estimate'],
#' #'              weights = rep(1/length(se), times = length(se)))
#' #' 
#' #' # this actually gives the same...2848206
#' #' my_est_kappa(Y = unlist(matrices_test$estimates[c(1,2), ] - 
#' #'                           approxmatrices_sparseSVD$approxmatrices[[10]][c(1,2), ]),
#' #'              se = as.vector(se_mat_test[c(1,2), ]),
#' #'              # mu2 = 5,
#' #'              mu2 = shrink2Grnas$ebci_obj$mu2['estimate'],
#' #'              # mu2 = shrink2Grnas$ebci_obj$mu2['uncorrected_estimate'],
#' #'              weights = rep(1/length(unlist(se_mat_test[c(1,2), ])), 
#' #'                            times = length(unlist(se_mat_test[c(1,2), ]))))
#' my_est_kappa <- function(Y, se, weights, mu2, X=NULL, delta=NULL) {
#'   # Y = matrices_test$estimates[1:50, 1] - approxmatrices_sparseSVD$approxmatrices[[10]][1:50, 1]
#'   # se = se_mat_test[1:50, 1]
#'   # weights = rep(1/length(se), times = length(se))
#'   # X=NULL; delta=NULL
#'   # mu2 = .000774
#'   
#'   if(is.null(X) | is.null(delta)) {
#'     eps = Y
#'   } else {
#'     eps = Y - X %*% delta
#'   }
#'   
#'   first = (sum( weights * (eps^4 - 6 * se * se * eps * eps + 3 * se^4))
#'             /
#'             (mu2*mu2 * sum(weights))
#'   )
#'   
#'   
#'   second = 1 + 
#'     (32 * sum(weights * weights * se^8)
#'       /
#'      (mu2*mu2 * sum(weights) * sum(weights * se^4))
#'     )
#'   
#'   
#'   
#'   return(max(first, second))
#' }
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' #' Calc the eb shrunk value
#' #' 
#' #' @param Y (vector) of numeric values, the Y_i's 
#' #'     in our case, this is unshrunk value - shrinkage point
#' #'     bc formula = unshrunk val - shrinkage point - 0
#' #'     which may be incorrect!
#' #' @param se (vector) of numeric values, the standard errors_i's
#' #'     this is \hat{\sigma}_i  in the formulas
#' #' @param weights (vector) of numeric values, the weights to put on each sample
#' #' @param mu2 (numeric) estimate of mu2 (from a prev step)
#' #' @param kappa (numeric) estimate of kappa (from a prev step)
#' #' @param X (matrix) of numeric values, the covariates
#' #' @param delta (vector) of coefs for covariates 
#' #' 
#' #' @example
#' #' theb = my_est_theb(
#' #'        Y = unlist(matrices_test$estimates[c(1,2), ] - 
#' #'        approxmatrices_sparseSVD$approxmatrices[[10]][c(1,2), ]),
#' #'        se = as.vector(se_mat_test[c(1,2), ]),
#' #'        # mu2 = 5,
#' #'        mu2 = shrink2Grnas$ebci_obj$mu2['estimate'],
#' #'        # mu2 = shrink2Grnas$ebci_obj$mu2['uncorrected_estimate'],
#' #'        kappa = shrink2Grnas$ebci_obj$kappa['estimate'],
#' #'        weights = rep(1/length(unlist(se_mat_test[c(1,2), ])), 
#' #'                      times = length(unlist(se_mat_test[c(1,2), ]))))
#' #'        theb$w_eb |> hist()
#' #'        theb$th_eb |> hist()
#' #'                
#' my_est_theb <- function(Y, se, weights, mu2, kappa, X=NULL, delta=NULL) {
#'   # Y = matrices_test$estimates[1:50, 1] - approxmatrices_sparseSVD$approxmatrices[[10]][1:50, 1]
#'   # se = se_mat_test[1:50, 1]
#'   # weights = rep(1/length(se), times = length(se))
#'   # X=NULL; delta=NULL
#'   # mu2 = .000774
#'   
#'   if(is.null(X) | is.null(delta)) {
#'     c = rep(0, times = length(Y))
#'   } else {
#'     c = X %*% delta
#'   }
#'   
#'   w_eb = mu2 / (mu2 + se*se) # vector
#'   
#'   th_eb = c + w_eb * (Y - c)
#'   
#'   
#'   # return(list(w_eb = w_eb, th_eb = th_eb))
#'   return(data.frame(w_eb = w_eb, th_eb = th_eb))
#' }
#' 
#' 
#' 
#' 
#' 
#' #' Calc the eb shrunk value's confidence interval
#' #' This step takes a while from the ebci::cva call 
#' #' 
#' #' @param th_eb (vector) of eb shrunk estimate (from a prev step)
#' #' @param w_eb (vector) of eb shrinkage factors (from a prev step)
#' #' @param se (vector) of numeric values, the standard errors_i's
#' #'     this is \hat{\sigma}_i  in the formulas
#' #' @param weights (vector) of numeric values, the weights to put on each sample
#' #' @param mu2 (numeric) estimate of mu2 (from a prev step)
#' #' @param kappa (numeric) estimate of kappa (from a prev step)
#' #'
#' #' @example
#' #' th_ci = my_est_ci(th_eb = theb$th_eb,
#' #'                   w_eb  = theb$w_eb,
#' #'                   se = as.vector(se_mat_test[c(1,2), ]),
#' #'                   mu2 = shrink2Grnas$ebci_obj$mu2['estimate'],
#' #'                   kappa = shrink2Grnas$ebci_obj$kappa['estimate'],
#' #'                   alpha = .1)
#' #' 
#' my_est_ci <- function(th_eb, w_eb, se, mu2, kappa, alpha)  {
#'   cv = rep(NA, times = length(se)) 
#'   m2_ = se*se / mu2
#'   for(i in 1:length(cv)) {
#'     # if(i %% 100 == 0){print(i)}
#'      
#'     # !!! This is kind of slow!!! (this is where the time will be spent)
#'     cv[i] = ebci::cva(m2 = m2_[i], kappa=kappa, alpha=alpha)$cv
#'   }
#'   
#'   len_eb = cv * w_eb * se
#'   
#'   return(data.frame(ci_lower = th_eb - len_eb,
#'                     ci_upper = th_eb + len_eb))
#' }
#' 
#' 
#' #' calc ebci by hand to compare
#' #'
#' #' @param Y (vector) of numeric values, the Y_i's 
#' #'     in our case, this is unshrunk value - shrinkage point
#' #'     bc formula = unshrunk val - shrinkage point - 0
#' #'     which may be incorrect!
#' #' @param se (vector) of numeric values, the standard errors_i's
#' #'     this is \hat{\sigma}_i  in the formulas
#' #' @param weights (vector) of numeric values, the weights to put on each sample
#' #' @param mu2 (numeric) estimate of mu2 (from a prev step)
#' #' @param kappa (numeric) estimate of kappa (from a prev step)
#' #' @param X (matrix) of numeric values, the covariates
#' #' @param delta (vector) of coefs for covariates 
#' #' @returns ebci results as a list:
#' #' $mu2: estimate of mu2
#' #' $kappa: estimate of kappa
#' #' $df: data.frame with eb value columns
#' #'     $th_eb: theta estimate
#' #'     $w_eb: shrinkage factor
#' #'     $ci_lower: eb CI lower value
#' #'     $ci_upper: eb CI upper value
#' #' 
#' #' @example
#' #' my_ebci_obj = 
#' #'   my_ebci(Y = unlist(matrices_test$estimates[c(1,2), ] - 
#' #'                        approxmatrices_sparseSVD$approxmatrices[[10]][c(1,2), ]),
#' #'           se = as.vector(se_mat_test[c(1,2), ]),
#' #'           weights = rep(1/length(as.vector(se_mat_test[c(1,2), ])), 
#' #'                         times = length(as.vector(se_mat_test[c(1,2), ]))),
#' #'           alpha = .1)
#' #' 
#' #' 
#' #' my_ebci_obj$mu2
#' #' my_ebci_obj$kappa
#' #' 
#' #' head(my_ebci_obj$df) 
#' #' 
#' my_ebci <- function(Y, se, weights, alpha, X=NULL, delta=NULL) {
#'   
#'   # estimate mu2 = E[(\theta - X^T delta)^2 |X]  and kappa ...
#'   mu2   = my_est_mu2(  Y=Y, se=se, weights=weights,          X=X, delta=delta)
#'   kappa = my_est_kappa(Y=Y, se=se, weights=weights, mu2=mu2, X=X, delta=delta)
#'   # estimate shrunk values
#'   th_df = my_est_theb( Y=Y, se=se, weights=weights, mu2=mu2, kappa=kappa, X=X, delta=delta)
#'   ci_df = my_est_ci(th_eb=th_df$th_eb, w_eb=th_df$w_eb, se=se, mu2=mu2, kappa=kappa, alpha=alpha)
#'   
#'  
#'   
#'   
#'   return(list(mu2 = mu2,
#'               kappa = kappa,
#'               df =  cbind(th_df, ci_df,
#'                           data.frame(Y=Y, se=se, weights=weights))))
#'   
#' }




source('../utils/myEBCI.r')



ALPHA = .1


ebci_data = data.frame(unshrunk_value = unlist(matrices_test$estimates[c(1,2), ]),
                       shrinkage_point = as.vector(approxmatrices_sparseSVD$approxmatrices[[10]][c(1,2), ]),
                       se = as.vector(se_mat_test[c(1,2), ]),
                       se4 = as.vector(se_mat_test[c(1,2), ])^4) |>
            dplyr::mutate(# weights = 1/n(), 
                          # weights = 1/se^2, 
                          weights = (1/se)^4,
                          Y = unshrunk_value - shrinkage_point)

ebci_data |> head()

ebci_obj = ebci::ebci(formula = 'unshrunk_value - shrinkage_point ~ 0',
                      data    = ebci_data, 
                      se = se,
                      weights = weights,
                      alpha = ALPHA)
my_ebci_obj = my_ebci(Y = ebci_data$Y,
                      se = ebci_data$se,
                      weights = ebci_data$weights,
                      alpha = ALPHA)

# Now, try estimating mu2 in a different way:
# (replace Y with the original unshrunk value in the estimation of mu2)
# unfort, the mu2 est is the same?? bc the first term is <0, so it just choose the second which depends only on se
my_est_mu2(  Y=ebci_data$Y, se=ebci_data$se, weights=ebci_data$weights, print=T)
new_mu2   = my_est_mu2(  Y=ebci_data$unshrunk_value, se=ebci_data$se, weights=ebci_data$weights, print=T)
new_kappa = my_est_kappa(Y=ebci_data$unshrunk_value, se=ebci_data$se, weights=ebci_data$weights, mu2=new_mu2, X=NULL, delta=NULL, print=T)

th_df = my_est_theb( Y=ebci_data$Y, se=se, weights=weights, mu2=new_mu2, X=X, delta=delta)
plot(th_df$w_eb, my_ebci_obj_newmu2$df$w_eb)
ci_df = my_est_ci(th_eb=th_df$th_eb, w_eb=th_df$w_eb, se=se, mu2=new_mu2, kappa=kappa, alpha=alpha)
plot(ci_df$len_eb, my_ebci_obj_newmu2$df$len_eb)
hist(ci_df$m2)
hist(ci_df$cv)
sum(ebci_data$weights^2 * ebci_data$se^4)
4000 / sum( ebci_data$se^2)
plot(ebci_data$unshrunk_value, ebci_data$Y); abline(a=0,b=1)


my_ebci_obj_newmu2 = my_ebci(Y = ebci_data$Y,
                             se = ebci_data$se,
                             weights = ebci_data$weights,
                             alpha = ALPHA,
                             mu2 = new_mu2,
                             kappa = new_kappa)





shrink2Grnas$ebci_obj$mu2;ebci_obj$mu2; my_ebci_obj$mu2; my_ebci_obj_newmu2$mu2
shrink2Grnas$ebci_obj$kappa; my_ebci_obj$kappa; my_ebci_obj_newmu2$kappa


shrink2Grnas$ebci_obj$df |> head()
my_ebci_obj$df |> head()

# compare (should be a y=x line)
ggplot(data.frame(ebci    =    ebci_obj$df$th_eb,   # theta estimate
                  my_ebci = my_ebci_obj$df$th_eb),
       aes(x=ebci, y = my_ebci)) + geom_point()

ggplot(data.frame(ebci    =    ebci_obj$df$w_eb,    # shrinkage factor
                  my_ebci = my_ebci_obj$df$w_eb),
       aes(x=ebci, y = my_ebci)) + geom_point()

ggplot(data.frame(ebci    =    ebci_obj$df$len_eb,  # CI length
                  my_ebci = my_ebci_obj$df |> dplyr::mutate(len_eb = (ci_upper - ci_lower) / 2) |> pull(len_eb)),
       aes(x=ebci, y = my_ebci)) + geom_point()


# compare to using different way to calc mu2 (want this to be diff)
ggplot(data.frame(use_Y        =    ebci_obj$df$th_eb,   # theta estimate
                  use_unshrunk = my_ebci_obj_newmu2$df$th_eb),
       aes(x=use_Y, y = use_unshrunk)) + geom_point()

ggplot(data.frame(use_Y        =    ebci_obj$df$w_eb,   # shrinkage factor
                  use_unshrunk = my_ebci_obj_newmu2$df$w_eb),
       aes(x=use_Y, y = use_unshrunk)) + geom_point()

ggplot(data.frame(use_Y        =    ebci_obj$df$len_eb,   # CI length
                  use_unshrunk = my_ebci_obj_newmu2$df$len_eb),
       aes(x=use_Y, y = use_unshrunk)) + geom_point()






# try constructing se's in a different way.
# Given the estimate mu and the p-value p, the direct method of getting the se is:
# abs(mu / qnorm(p = (1/2)*p, mean = 0, sd = 1))


# denominator is small when p is large?
qnorm(p = (1/2)*.99, mean = 0, sd = 1)

mus = seq(-1, 1, length.out = 100)
ps = seq(0, 1, length.out = 100)

mu_ = c()
p_ = c()
se_ = c()
z_ = c()
for(mu in mus) {
  for(p in ps) {
    mu_ = c(mu_, mu)
    p_  = c(p_, p)
    z = qnorm(p = (1/2)*p, mean = 0, sd = 1)
    z_ = c(z_, z)
    se_ = c(se_, abs(mu / z))
  }
}

se_calc_df = data.frame(mu = mu_, p = p_, se = se_, z=z_)


dim(se_calc_df)


ggplot(se_calc_df |> filter(0 <= p & p <= .1), 
       aes(x = mu, y = p, fill = se)) + 
  geom_raster()

ggplot(se_calc_df) +
  geom_point(aes(x = p, y = z))


se_calc_df |> filter(p == 0)


delta = .01
se_calc_df = se_calc_df |> dplyr::mutate(se_mod = (abs(mu) + delta)/ (z + delta))




ggplot(se_calc_df |> filter(0 <= p & p <= .1), 
       aes(x = se, y = se_mod)) + 
  geom_abline(aes(slope = 1, intercept = 0)) +
  geom_point()

 ggplot(se_calc_df |> filter(0 <= p & p <= .1), 
       aes(x = mu, y = p, fill = se_mod)) + 
  geom_raster()

 
 
 
# as n increases, mu hat increases??
second = c()
second_seweight = c()
for(n_ in seq(from = 10, to = 200, by = 10)) {
  sigma = rep(.01, times = n_)
  weights = rep(1/n_, times = n_)
  weights_se = 1/sigma^2
  second = c(2 * sum(weights^2 * sigma^4) / (sum(weights) * sum(weights * sigma^2)),
             second)
  second_seweight = c(2 * sum(weights_se^2 * sigma^4) / (sum(weights_se) * sum(weights_se * sigma^2)),
                      second_seweight)
}
plot(second)
plot(second_seweight)
plot(second, second_seweight)

# as se increases, mu hat decreases
v = c(.001, .005, .01, .05, .1)
second = c()
for(s in v) {
  sigma = rep(s, times = 100)
  weights = rep(1/100, times = 100)
  second = c(2 * sum(weights^2 * sigma^4) / (sum(weights) * sum(weights * sigma^2)),
             second)
}
plot(v, second)


#' #' Spline that does this:
#' #' 
#' #'       |                . 
#' #'       |              .  <- linear, y=x
#' #'       |            .   
#' #'epsilon+  -  -  - + <- match 1st deriv
#' #'       |        . |
#' #'       |      .     <- quadratic
#' #' delta +.  .      |
#' #'       |
#' #'       -----------+
#' #'       0          epsilon
#' #'        
#' #' To ensure the function does not decrease on [0, epsilon],
#' #' the parameters must satisfy delta <= epsilon / 2
#' #' @param delta (numeric) function parameter
#' #' @param epsilon (numeric) function parameter
#' #' @example 
#' #' x = seq(0, 1, length.out = 200) 
#' #' my_spline = make_spline(delta = .3, epsilon = .6) 
#' #' plot(x, sapply(X = x, FUN = my_spline), type = 'l')
#' #' 
#' #' my_spline = make_spline(delta = .3, epsilon = .4) <- bad params
#' make_spline <- function(delta, epsilon) {
#'   fn <- function(x) {
#'     if(x < epsilon) {
#'       return(  (delta / epsilon**2)      * x**2 + 
#'                  (1 - 2 * delta / epsilon) * x + 
#'                  delta                             )
#'     } else {
#'       return(x)
#'     }
#'   }
#'   return(fn)
#' } 


x = seq(0, 1, length.out = 200) 
my_spline = make_spline(delta = .01, epsilon = .03) 
plot(x, sapply(X = x, FUN = my_spline), type = 'l')
 





# === === === === === === === === === === === === === === === === === === === === === === === === === === === === === === === === === === === === ===
# === === === === === === === === === === === === === === === === === === === === === === === === === === === === === === === === === === === === ===
# === === === === === === === === === === === === === === === === === === === === === === === === === === === === === === === === === === === === ===





# === === === === === === === === === === === === === === === === === === === === === 
#                                   PLOT
# === === === === === === === === === === === === === === === === === === === === === 




print(sprintf("[%s] START: Plotting", Sys.time()))


plot_folder = '../plots/matrixSampSplit/shrinkMatrixTest/'
save_folder = '../saves/matrixSampSplit/shrinkMatrixTest/'

dir.create(plot_folder)
dir.create(save_folder)
# dir.create(sprintf('%s/heatmaps/', plot_folder))


# =========================== shrink to fixed matrix

# real data, shrink to arbitrary fixed matrix 
plot_folder = '../plots/matrixSampSplit/shrinkMatrixTest/fixedMatrix/'

dir.create(plot_folder)
dir.create(sprintf('%s/points/', plot_folder))
dir.create(sprintf('%s/heatmaps/', plot_folder))

shrink_df = read.csv("../saves/matrixSampSplit/shrinkMatrixTest/shrinkMatrix_fixed.csv")

plot_shrink_results(shrink_df, plot_folder)



# simulated w random theta ~ N(0,1)
plot_folder = '../plots/matrixSampSplit/shrinkMatrixTest/simMatRTheta/'
shrink_df = read.csv("../saves/matrixSampSplit/shrinkMatrixTest/shrinkMatrix_sim_randtheta.csv")

dir.create(plot_folder)
dir.create(sprintf('%s/points/', plot_folder))
dir.create(sprintf('%s/heatmaps/', plot_folder))


plot_shrink_results(shrink_df=shrink_df, plot_folder=plot_folder, order_rowscols=F)


# simulated w random theta ~ 0 + some N(.., ..)'s
plot_folder = '../plots/matrixSampSplit/shrinkMatrixTest/simMatRTheta0/'
shrink_df = read.csv("../saves/matrixSampSplit/shrinkMatrixTest/shrinkMatrix_sim_randtheta0.csv")

dir.create(plot_folder)
dir.create(sprintf('%s/points/', plot_folder))
dir.create(sprintf('%s/heatmaps/', plot_folder))


plot_shrink_results(shrink_df=shrink_df, plot_folder=plot_folder, order_rowscols=F)


# simulated  w correct se 100%
plot_folder = '../plots/matrixSampSplit/shrinkMatrixTest/simMatrix/'
shrink_df = read.csv("../saves/matrixSampSplit/shrinkMatrixTest/shrinkMatrix_simestimates.csv")

dir.create(plot_folder)
dir.create(sprintf('%s/points/', plot_folder))
dir.create(sprintf('%s/heatmaps/', plot_folder))


plot_shrink_results(shrink_df=shrink_df, plot_folder=plot_folder, order_rowscols=F)

# simulated w incorrect se 125%
plot_folder = '../plots/matrixSampSplit/shrinkMatrixTest/simMatrix125/'
shrink_df = read.csv("../saves/matrixSampSplit/shrinkMatrixTest/shrinkMatrix_simest_125.csv")

dir.create(plot_folder)
dir.create(sprintf('%s/points/', plot_folder))
dir.create(sprintf('%s/heatmaps/', plot_folder))


plot_shrink_results(shrink_df=shrink_df, plot_folder=plot_folder, order_rowscols=F)

# simulated w incorrect se 150%
plot_folder = '../plots/matrixSampSplit/shrinkMatrixTest/simMatrix150/'
shrink_df = read.csv("../saves/matrixSampSplit/shrinkMatrixTest/shrinkMatrix_simest_150.csv")

dir.create(plot_folder)
dir.create(sprintf('%s/points/', plot_folder))
dir.create(sprintf('%s/heatmaps/', plot_folder))


plot_shrink_results(shrink_df=shrink_df, plot_folder=plot_folder, order_rowscols=F)

# simulated w incorrect se 200%
plot_folder = '../plots/matrixSampSplit/shrinkMatrixTest/simMatrix200/'
shrink_df = read.csv("../saves/matrixSampSplit/shrinkMatrixTest/shrinkMatrix_simest_200.csv")

dir.create(plot_folder)
dir.create(sprintf('%s/points/', plot_folder))
dir.create(sprintf('%s/heatmaps/', plot_folder))


plot_shrink_results(shrink_df=shrink_df, plot_folder=plot_folder, order_rowscols=F)

# simulated w incorrect se 75%
plot_folder = '../plots/matrixSampSplit/shrinkMatrixTest/simMatrix075/'
shrink_df = read.csv("../saves/matrixSampSplit/shrinkMatrixTest/shrinkMatrix_simest_075.csv")

dir.create(plot_folder)
dir.create(sprintf('%s/points/', plot_folder))
dir.create(sprintf('%s/heatmaps/', plot_folder))


plot_shrink_results(shrink_df=shrink_df, plot_folder=plot_folder, order_rowscols=F)

# simulated w incorrect se 50%
plot_folder = '../plots/matrixSampSplit/shrinkMatrixTest/simMatrix050/'
shrink_df = read.csv("../saves/matrixSampSplit/shrinkMatrixTest/shrinkMatrix_simest_050.csv")

dir.create(plot_folder)
dir.create(sprintf('%s/points/', plot_folder))
dir.create(sprintf('%s/heatmaps/', plot_folder))


plot_shrink_results(shrink_df=shrink_df, plot_folder=plot_folder, order_rowscols=F)

# real data- sparsePCA - 2 grnas
plot_folder = '../plots/matrixSampSplit/shrinkMatrixTest/2grna/'
shrink_df = read.csv("../saves/matrixSampSplit/shrinkMatrixTest/shrinkMatrix_2grna.csv")

dir.create(plot_folder)
dir.create(sprintf('%s/points/', plot_folder))
dir.create(sprintf('%s/heatmaps/', plot_folder))


plot_shrink_results(shrink_df=shrink_df, plot_folder=plot_folder, order_rowscols=T)


# real data- sparsePCA - 10 grnas
plot_folder = '../plots/matrixSampSplit/shrinkMatrixTest/10grna/'
shrink_df = read.csv("../saves/matrixSampSplit/shrinkMatrixTest/shrinkMatrix_10grna.csv")

dir.create(plot_folder)
dir.create(sprintf('%s/points/', plot_folder))
dir.create(sprintf('%s/heatmaps/', plot_folder))


plot_shrink_results(shrink_df=shrink_df, plot_folder=plot_folder, order_rowscols=T)

# real data- sparsePCA - 25 grnas
plot_folder = '../plots/matrixSampSplit/shrinkMatrixTest/25grna/'
shrink_df = read.csv("../saves/matrixSampSplit/shrinkMatrixTest/shrinkMatrix_25grna.csv")

dir.create(plot_folder)
dir.create(sprintf('%s/points/', plot_folder))
dir.create(sprintf('%s/heatmaps/', plot_folder))


plot_shrink_results(shrink_df=shrink_df, plot_folder=plot_folder, order_rowscols=T)

# real data- sparsePCA - 50 grnas
plot_folder = '../plots/matrixSampSplit/shrinkMatrixTest/50grna/'
shrink_df = read.csv("../saves/matrixSampSplit/shrinkMatrixTest/shrinkMatrix_50grna.csv")

dir.create(plot_folder)
dir.create(sprintf('%s/points/', plot_folder))
dir.create(sprintf('%s/heatmaps/', plot_folder))


plot_shrink_results(shrink_df=shrink_df, plot_folder=plot_folder, order_rowscols=T)


# real data- lowrank - 50 grnas
plot_folder = '../plots/matrixSampSplit/shrinkMatrixTest/lowrank50grna/'
shrink_df = read.csv("../saves/matrixSampSplit/shrinkMatrixTest/shrinkMatrix_50grna_lowrank.csv")

dir.create(plot_folder)
dir.create(sprintf('%s/points/', plot_folder))
dir.create(sprintf('%s/heatmaps/', plot_folder))


plot_shrink_results(shrink_df=shrink_df, plot_folder=plot_folder, order_rowscols=T)



if(F) {
  
  
  
  # put into a function 
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
  
  
  
}



# =========================== shrink to random matrix


plot_folder = '../plots/matrixSampSplit/shrinkMatrixTest/'
save_folder = '../saves/matrixSampSplit/shrinkMatrixTest/'

dir.create(plot_folder)
dir.create(save_folder)
dir.create(sprintf('%s/heatmaps/', plot_folder))

shrink_df = read.csv("../saves/matrixSampSplit/shrinkMatrixTest/shrinkMatrix_random.csv")


# === Matrix Heatmaps ===

# Plot shrinkage results as heatmap matrices


# get an ordering for genes/grna (marginal clustering)
# shrink_df = shrinkDF

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









# add original estimates (all the loaded dfs should have these) (2x for plotting)
plot_df = rbind(shrink_df |> 
                  dplyr::mutate(rank = 'random'),
                shrink_df |>
                  dplyr::mutate(rank = 'unshrunk') |>
                  dplyr::mutate(shrinkage_point = unshrunk_value,   # set all values to the original estimates
                                shrunk_value    = unshrunk_value, 
                                lower_ci        = unshrunk_value - qnorm(.95) * se,
                                upper_ci        = unshrunk_value + qnorm(.95) * se) 
)


plot_df = rbind(shrink_df |> dplyr::mutate(rank = 'random'),
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


plot_df$rank = factor(plot_df$rank, levels = c('unshrunk', 'random'))  


# heatmap of diff values side by side
plot_df_heatmap =
  rbind(shrink_df |> dplyr::select(grna, gene, unshrunk_value)  |> dplyr::rename(val =  unshrunk_value) |> dplyr::mutate(type = 'unshrunk'), 
        shrink_df |> dplyr::select(grna, gene, shrinkage_point) |> dplyr::rename(val = shrinkage_point) |> dplyr::mutate(type = 'shrinkagepoint'), 
        shrink_df |> dplyr::select(grna, gene, shrunk_value)    |> dplyr::rename(val =    shrunk_value) |> dplyr::mutate(type = 'shrunk'),
        shrink_df |> dplyr::mutate(significant = (upper_ci <= 0) | (0 <= lower_ci)) |> 
                     dplyr::select(grna, gene, significant)    |> dplyr::rename(val =    significant) |> dplyr::mutate(type = 'significant'))





# add idx's- (from gene or grna name --> number on the axis)
plot_df_heatmap = merge(gene_index, plot_df_heatmap, by = 'gene')
plot_df_heatmap = merge(grna_index, plot_df_heatmap, by = 'grna')


plot_df_heatmap$type = factor(plot_df_heatmap$type, levels = c('unshrunk', 'shrinkagepoint', 'shrunk', 'significant'))  



# plot
p = ggplot(plot_df_heatmap) +
  geom_raster(aes(x = grna_idx, y = gene_idx, fill = val)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x='grna', y = 'gene', fill = NULL, title = 'Estimates After Shrinking to Random Points ~N(0,1)') +
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


ggsave(plot = p, filename = sprintf('%sheatmaps/heatmap.pdf', plot_folder), height = 7, width = 7)









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
  facet_grid(cols = vars(rank)) +
  theme_classic()

ggsave(plot = p, filename = sprintf('%spoints/shrunkvsunshrunk.pdf', plot_folder), height = 6, width = 7)


# plot shrunk vs unshrunk centered
p = ggplot(plot_df[subsample_idx, ] |> filter(rank != 'unshrunk')) +
  geom_abline(aes(slope = 1, intercept = 0), color = 'gray') +
  geom_point(aes(x = unshrunk_value - shrinkage_point, y = shrunk_value - shrinkage_point, color = se)) +
  labs(title = 'Shrunk vs Unshrunk Estimates Centered', x = 'unshrunk - shrinkage point', y = 'shrunk - shrinkage point') +
  coord_cartesian(xlim = c(-2, 2), ylim = c(-2, 2)) +
  facet_grid(# rows = vars(approxmethod),
             cols = vars(rank)) +
  theme_classic()

ggsave(plot = p, filename = sprintf('%spoints/shrunkvsunshrunkcent.pdf', plot_folder), height = 6, width = 7)


# plot unshrunk vs shrinkage point
p = ggplot(plot_df[subsample_idx, ] |> filter(rank != 'unshrunk')) +
  geom_abline(aes(slope = 1, intercept = 0), color = 'gray') +
  geom_point(aes(x = shrunk_value, y = shrinkage_point, color = se)) +
  labs(title = 'Shrinkage Point vs Shrunk Estimates', x = 'shrunk', y = 'shrinkage point') +
  coord_cartesian(xlim = c(-2, 2), ylim = c(-2, 2)) +
  facet_grid(# rows = vars(approxmethod),
             cols = vars(rank)) +
  theme_classic()

ggsave(plot = p, filename =sprintf('%spoints/shrinkptvsshrunk.pdf', plot_folder), height = 6, width = 7)







# === Shrinkage CIs ===
set.seed(12345)
# get the same gene/grna tests across all methods and ranks (to be better visually, should do to prev plots too)
sample_tests = plot_df |> filter(rank == 'unshrunk') |> slice_sample(n = 15) |> select(grna, gene) |> mutate(x = 1:n())
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


ggsave(plot = p, filename = sprintf('%spoints/CIs.pdf', plot_folder), height = 6, width = 15)


# 90% CI lengths before and after 
set.seed(12345)
p = ggplot(plot_df[sample(x=1:nrow(plot_df), size = 10000), ] |> 
             mutate(unshrunk_ci_length = 2*qnorm(.95)*se,
                    shrunk_ci_length = upper_ci - lower_ci)) +
  geom_abline(aes(intercept = 0, slope = 1)) +
  geom_point(aes(x = unshrunk_ci_length, y = shrunk_ci_length), alpha = .8) +
  labs(x = 'Unshrunk CI Lengths', y = 'Shrunk CI Lengths', title = 'CI Lengths') +
  facet_grid(#rows = vars(approxmethod),
             cols = vars(rank)) +
  theme_classic()


ggsave(plot = p, filename =sprintf('%spoints/CIlengths.pdf', plot_folder), height = 6, width = 7)



# =========================== 2 GRNAs


print(sprintf("[%s] START: Plotting", Sys.time()))


plot_folder = '../plots/matrixSampSplit/shrinkMatrixTest/'
save_folder = '../saves/matrixSampSplit/shrinkMatrixTest/'

dir.create(plot_folder)
dir.create(save_folder)
dir.create(sprintf('%s/heatmaps/', plot_folder))

shrink_df = read.csv("../saves/matrixSampSplit/shrinkMatrixTest/shrinkMatrix_2grna.csv")


# === Matrix Heatmaps ===

# Plot shrinkage results as heatmap matrices


# get an ordering for genes/grna (marginal clustering)
# shrink_df = shrinkDF

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









# add original estimates (all the loaded dfs should have these) (2x for plotting)
plot_df = rbind(shrink_df |> 
                  dplyr::mutate(rank = 'random'),
                shrink_df |>
                  dplyr::mutate(rank = 'unshrunk') |>
                  dplyr::mutate(shrinkage_point = unshrunk_value,   # set all values to the original estimates
                                shrunk_value    = unshrunk_value, 
                                lower_ci        = unshrunk_value - qnorm(.95) * se,
                                upper_ci        = unshrunk_value + qnorm(.95) * se) 
)


plot_df = rbind(shrink_df |> dplyr::mutate(rank = 'random'),
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


plot_df$rank = factor(plot_df$rank, levels = c('unshrunk', 'random'))  


# heatmap of diff values side by side
plot_df_heatmap =
  rbind(shrink_df |> dplyr::select(grna, gene, unshrunk_value)  |> dplyr::rename(val =  unshrunk_value) |> dplyr::mutate(type = 'unshrunk'), 
        shrink_df |> dplyr::select(grna, gene, shrinkage_point) |> dplyr::rename(val = shrinkage_point) |> dplyr::mutate(type = 'shrinkagepoint'), 
        shrink_df |> dplyr::select(grna, gene, shrunk_value)    |> dplyr::rename(val =    shrunk_value) |> dplyr::mutate(type = 'shrunk'),
        shrink_df |> dplyr::mutate(significant = (upper_ci <= 0) | (0 <= lower_ci)) |> 
          dplyr::select(grna, gene, significant)    |> dplyr::rename(val =    significant) |> dplyr::mutate(type = 'significant'))





# add idx's- (from gene or grna name --> number on the axis)
plot_df_heatmap = merge(gene_index, plot_df_heatmap, by = 'gene')
plot_df_heatmap = merge(grna_index, plot_df_heatmap, by = 'grna')


plot_df_heatmap$type = factor(plot_df_heatmap$type, levels = c('unshrunk', 'shrinkagepoint', 'shrunk', 'significant'))  






# plot
p = ggplot(plot_df_heatmap) +
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


ggsave(plot = p, filename = sprintf('%sheatmaps/heatmap_2grna.pdf', plot_folder), height = 7, width = 7)









# === Shrinkage Changes ===
dir.create(sprintf('%spoints/', plot_folder))
set.seed(12345)
subsample_idx = sample(1:nrow(plot_df), 8000)

# plot shrunk vs unshrunk
p = ggplot(plot_df[subsample_idx, ] |> filter(rank != 'unshrunk')) +
  geom_abline(aes(slope = 1, intercept = 0), color = 'gray', alpha = .6) +
  geom_point(aes(x = unshrunk_value, y = shrunk_value, color = se)) +
  labs(title = 'Shrunk vs Unshrunk Estimates', x = 'unshrunk', y = 'shrunk') +
  coord_cartesian(xlim = c(-2, 2), ylim = c(-2, 2)) +
  facet_grid(cols = vars(rank)) +
  theme_classic()

ggsave(plot = p, filename = sprintf('%spoints/shrk_unshrk_2grna.pdf', plot_folder), height = 6, width = 7)


# plot shrunk vs unshrunk centered
p = ggplot(plot_df[subsample_idx, ] |> filter(rank != 'unshrunk')) +
  geom_abline(aes(slope = 1, intercept = 0), color = 'gray') +
  geom_point(aes(x = unshrunk_value - shrinkage_point, y = shrunk_value - shrinkage_point, color = se)) +
  labs(title = 'Shrunk vs Unshrunk Estimates Centered', x = 'unshrunk - shrinkage point', y = 'shrunk - shrinkage point') +
  coord_cartesian(xlim = c(-2, 2), ylim = c(-2, 2)) +
  facet_grid(# rows = vars(approxmethod),
    cols = vars(rank)) +
  theme_classic()

ggsave(plot = p, filename = sprintf('%spoints/shrk_unshrkcen_2grna.pdf', plot_folder), height = 6, width = 7)


# plot unshrunk vs shrinkage point
p = ggplot(plot_df[subsample_idx, ] |> filter(rank != 'unshrunk')) +
  geom_abline(aes(slope = 1, intercept = 0), color = 'gray') +
  geom_point(aes(x = shrunk_value, y = shrinkage_point, color = se)) +
  labs(title = 'Shrinkage Point vs Shrunk Estimates', x = 'shrunk', y = 'shrinkage point') +
  coord_cartesian(xlim = c(-2, 2), ylim = c(-2, 2)) +
  facet_grid(# rows = vars(approxmethod),
    cols = vars(rank)) +
  theme_classic()

ggsave(plot = p, filename =sprintf('%spoints/shrkpt_shrk_2grna.pdf', plot_folder), height = 6, width = 7)







# === Shrinkage CIs ===
set.seed(12345)
# get the same gene/grna tests across all methods and ranks (to be better visually, should do to prev plots too)
sample_tests = plot_df |> filter(rank == 'unshrunk') |> slice_sample(n = 15) |> select(grna, gene) |> mutate(x = 1:n())
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


ggsave(plot = p, filename = sprintf('%spoints/CIs_2grna.pdf', plot_folder), height = 6, width = 15)


# 90% CI lengths before and after 
set.seed(12345)
p = ggplot(plot_df[sample(x=1:nrow(plot_df), size = 8000), ] |> 
             mutate(unshrunk_ci_length = 2*qnorm(.95)*se,
                    shrunk_ci_length = upper_ci - lower_ci)) +
  geom_abline(aes(intercept = 0, slope = 1)) +
  geom_point(aes(x = unshrunk_ci_length, y = shrunk_ci_length), alpha = .8) +
  labs(x = 'Unshrunk CI Lengths', y = 'Shrunk CI Lengths', title = 'CI Lengths') +
  facet_grid(#rows = vars(approxmethod),
    cols = vars(rank)) +
  theme_classic()


ggsave(plot = p, filename =sprintf('%spoints/CIlengths_2grna.pdf', plot_folder), height = 6, width = 7)




