# Run matrix approximation + ebci shrinkage (w/ sample splitting)
# (util functions should be in utils/matrix_shrinkage.r)
# (saved files should be in saves/)
# Analysis steps:
#   0. Load and prep data
#   1. Get a lower structure approximation of effects matrix on train
#   2. Shrink towards lower structure on test
#   3. Make some plots


print(sprintf("[%s] START: Run matrix approximation (w/ sample splitting) + ebci shrinkage", Sys.time()))


# === Settings/Parameters ===
ALPHA = .1

# approximations and shrinkage choices
PLOT_RANKS = c(1, 3, 5, 10, 15, 30)
APPROX_SHRINK_CHOICES = list('1' = list(approx='lowrank', rank= 5),
                             '2' = list(approx='lowrank', rank=10),
                             '3' = list(approx='lowrank', rank=15),
                             '4' = list(approx='lowrank', rank=30),
                             '5' = list(approx='sparseSVD', rank= 5),
                             '6' = list(approx='sparseSVD', rank=10),
                             '7' = list(approx='sparseSVD', rank=15),
                             '8' = list(approx='sparseSVD', rank=30))

APPROX_SHRINK_CHOICES = list('1' = list(approx='lowrank', rank= 5),
                             '7' = list(approx='sparseSVD', rank=15),
                             '8' = list(approx='sparseSVD', rank=30))




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


# source('../utils/myKmeans.R')
# source('../utils/cluster_and_ebci_shrinkage.r')
source('../utils/matrix_shrinkage.r')

# plot_folder = '../saves/matrixSampSplit/'

dir.create('../saves/matrixSampSplit/')
dir.create('../saves/matrixSampSplit/shrinkMatrix/')

dir.create('../plots/matrixSampSplit/')
dir.create('../plots/matrixSampSplit/topGeneGrna/') # compare diff ways to rank top genes/grna
dir.create('../plots/matrixSampSplit/topGeneGrna/toptstat/')
dir.create('../plots/matrixSampSplit/topGeneGrna/topscoreTrain/')
dir.create('../plots/matrixSampSplit/topGeneGrna/topscoreTest/')

dir.create('../plots/matrixSampSplit/approxMatrix/') # diff ways to approximate matrix
dir.create('../plots/matrixSampSplit/approxMatrix/lowrank/')
dir.create('../plots/matrixSampSplit/approxMatrix/lowrank/estimate/')
dir.create('../plots/matrixSampSplit/approxMatrix/lowrank/tstat/')
dir.create('../plots/matrixSampSplit/approxMatrix/lowrank/significance/')
dir.create('../plots/matrixSampSplit/approxMatrix/lowrank/foldchange/')
dir.create('../plots/matrixSampSplit/approxMatrix/lowrank/se/')
dir.create('../plots/matrixSampSplit/approxMatrix/sparseSVD/')
dir.create('../plots/matrixSampSplit/approxMatrix/sparseSVD/estimate/')
dir.create('../plots/matrixSampSplit/approxMatrix/sparseSVD/tstat/')
dir.create('../plots/matrixSampSplit/approxMatrix/sparseSVD/significance/')
dir.create('../plots/matrixSampSplit/approxMatrix/sparseSVD/foldchange/')
dir.create('../plots/matrixSampSplit/approxMatrix/sparseSVD/se/')

dir.create('../plots/matrixSampSplit/shrinkMatrix/') # shrinkage

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
      mutate(tstat = estimate / se) # make a tstat column = estimate / se
  }
}


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
# tstat
matrices = make_matrices(effects_df = effects_df  |> 
                           filter((gene %in% (topgenes |> arrange(desc(sum_tstats)) |> head(2000) |> pull(gene)))
                                  &
                                    (grna %in% (topgrnas |> arrange(desc(sum_tstats)) |> head(50) |> pull(grna)))), 
                         save_plots_filepath = '../plots/matrixSampSplit/topGeneGrna/toptstat/')

# score
effects_df = effects_df |> filter((gene %in% (topgenes |> arrange(desc(score)) |> head(2000) |> pull(gene)))
                                  & 
                                    (grna %in% (topgrnas |> arrange(desc(score)) |> head(50) |> pull(grna))))


matrices = make_matrices(effects_df = effects_df, save_plots_filepath = '../plots/matrixSampSplit/topGeneGrna/topscoreTrain/')

# test (using train's top scores)
chosengenes = effects_df |> pull(gene) |> unique() |> sort()
chosengrnas = effects_df |> pull(grna) |> unique() |> sort()

effects_df_test = effects_list[['sceptre']][['test']] |> 
  dplyr::filter(test == 'discovery' & (gene %in% chosengenes) & (grna %in% chosengrnas))

matrices_test = make_matrices(effects_df = effects_df_test, save_plots_filepath = '../plots/matrixSampSplit/topGeneGrna/topscoreTest/')

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
                    FUN = function(x){pmin(x, quantile(matrices_test$se |> unlist(), .91))})


assertthat::assert_that(quantile(matrices$se |> unlist(), .91) < 999, 
                        msg = "Fix: set limit for se to be < 999 (e.g. set quantile smaller)")



print(sprintf("[%s] Approximate Matrices and Shrink", Sys.time()))

## approximate matrices
print(sprintf("[%s]     - approximate matrices", Sys.time()))
# gather all ranks requested
allranks = c()
for(name in names(APPROX_SHRINK_CHOICES)) {
  allranks = c(allranks, APPROX_SHRINK_CHOICES[[name]]$rank)
}

### just for creating and saving plots: done on train split
if(T) { # Also do for other values like tstat, se, etc...
  
  lowrankres = approx_matrix(mat = matrices$estimates,
                             method = 'lowrank', 
                             ranks = PLOT_RANKS,
                             save_plots_filepath = '../plots/matrixSampSplit/approxMatrix/lowrank/estimate/',
                             color_limits = c(-2, 2)) # matrices$estimates |> as.matrix() |> as.vector() |> hist()
  
  
  # testggplot = ggplot(data = data.frame(x = 1:10, y = 1:10), aes(x=x,y=y)) + geom_point()
  # ggsave(plot = testggplot, filename =  '../plots/matrixSampSplit/approxMatrix/lowrank/estimate/aaaaaaaaaaaa.pdf')
  # ggsave(plot = testggplot, filename =  '../plots/matrixSampSplit/approxMatrix/lowrank/estimate/aaaaaaaaaaaaa.pdf')
  
  
  lowrankres = approx_matrix(mat = matrices$tstat,
                             method = 'lowrank', 
                             ranks = PLOT_RANKS,
                             save_plots_filepath = '../plots/matrixSampSplit/approxMatrix/lowrank/tstat/',
                             color_limits = c(-4, 4)) # matrices$tstat |> as.matrix() |> as.vector() |> hist()
  
  lowrankres = approx_matrix(mat = matrices$significance,
                             method = 'lowrank', 
                             ranks = PLOT_RANKS,
                             save_plots_filepath = '../plots/matrixSampSplit/approxMatrix/lowrank/significance/',
                             color_limits = c(-1, 1))
  
  lowrankres = approx_matrix(mat = matrices$foldchange,
                             method = 'lowrank', 
                             ranks = PLOT_RANKS,
                             save_plots_filepath = '../plots/matrixSampSplit/approxMatrix/lowrank/foldchange/',
                             color_limits = c(0, 2)) # matrices$foldchange |> as.matrix() |> as.vector() |> hist()
  
  lowrankres = approx_matrix(mat = matrices$se,
                             method = 'lowrank', 
                             ranks = PLOT_RANKS,
                             save_plots_filepath = '../plots/matrixSampSplit/approxMatrix/lowrank/se/',
                             color_limits = c(0, 3)) # matrices$se |> as.matrix() |> as.vector() |> hist()
  rm(lowrankres)
  
  sparsesvdres = approx_matrix(mat = matrices$estimates, 
                               method = 'sparseSVD', 
                               ranks = PLOT_RANKS,
                               save_plots_filepath="../plots/matrixSampSplit/approxMatrix/sparseSVD/estimate/", 
                               color_limits = c(-2, 2), 
                               methodParams=list(type = 'standard', 
                                                 sumabs = .35,  # should be between 0-1, 
                                                 # sumabsu = 4, sumabsv = 4, # between 1 and sqrt(#col or #rows)
                                                 niter = 100,
                                                 trace = FALSE)) 
  
  sparsesvdres = approx_matrix(mat = matrices$tstat, 
                               method = 'sparseSVD', 
                               ranks = PLOT_RANKS,
                               save_plots_filepath="../plots/matrixSampSplit/approxMatrix/sparseSVD/tstat/", 
                               color_limits = c(-4, 4), 
                               methodParams=list(type = 'standard', 
                                                 sumabs = .35,  # should be between 0-1, 
                                                 # sumabsu = 4, sumabsv = 4, # between 1 and sqrt(#col or #rows)
                                                 niter = 100,
                                                 trace = FALSE))
  
  sparsesvdres = approx_matrix(mat = matrices$foldchange, 
                               method = 'sparseSVD', 
                               ranks = PLOT_RANKS,
                               save_plots_filepath="../plots/matrixSampSplit/approxMatrix/sparseSVD/foldchange/", 
                               color_limits = c(-4, 4), 
                               methodParams=list(type = 'standard', 
                                                 sumabs = .35,  # should be between 0-1, 
                                                 # sumabsu = 4, sumabsv = 4, # between 1 and sqrt(#col or #rows)
                                                 niter = 100,
                                                 trace = FALSE))
  rm(sparsesvdres)
  
}


### approximate matrices for shrinkage (on train split)
approxmatrices = list()

approxmatrices$lowrank = 
  approx_matrix(mat = matrices$estimates, 
                method = 'lowrank', 
                ranks = allranks) 

approxmatrices$sparseSVD = 
  approx_matrix(mat = matrices$estimates, 
                method = 'sparseSVD', 
                ranks = allranks,
                methodParams=list(type = 'standard', 
                                  sumabs = .35,  # should be between 0-1, 
                                  # sumabsu = 4, sumabsv = 4, # between 1 and sqrt(#col or #rows)
                                  niter = 100,
                                  trace = FALSE)) 




## perform ebci shrinkage to approx matrix point
print(sprintf("[%s]     - perform ebci shrinkage", Sys.time()))
for(name in names((APPROX_SHRINK_CHOICES))) {
  approxMethod = APPROX_SHRINK_CHOICES[[name]]$approx
  rank         = APPROX_SHRINK_CHOICES[[name]]$rank
  print(sprintf("[%s]         + method=%s  rank=%02.f", Sys.time(), approxMethod, rank))
  
  
  shrinkDF = 
    shrink_matrix(   unshrunk_mat = matrices_test$estimates,
                     shrinkpoint_mat = approxmatrices[[approxMethod]]$approxmatrices[[rank]],
                     se_mat = se_mat_test,
                     ALPHA = ALPHA) # ~15 mins
  
  
  
  shrinkDFPlots = shrink_matrix_plots(df = shrinkDF)
  
  pdf(sprintf("../plots/matrixSampSplit/shrinkMatrix/shrinkMatrix_%s_r=%02.f.pdf",approxMethod, rank))
  gridExtra::grid.arrange(grobs = shrinkDFPlots)
  dev.off()
  
  write.csv(x = shrinkDF,
            file = sprintf("../saves/matrixSampSplit/shrinkMatrix/shrinkMatrix_%s_r=%02.f.csv",approxMethod, rank))
  
}




print(sprintf("[%s] END", Sys.time()))



