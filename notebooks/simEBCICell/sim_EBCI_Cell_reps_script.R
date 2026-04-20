# libraries
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(ggplot2))
ggplot2::theme_set(theme_bw())
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(ebci)) # robust emp shrinkage
suppressPackageStartupMessages(library(assertthat)) 
suppressPackageStartupMessages(library(softImpute))
suppressPackageStartupMessages(library(latex2exp))
suppressPackageStartupMessages(library(biclust))
suppressPackageStartupMessages(library(future.apply))
plan(multisession, workers = 10)
# our code
source('../../utils/matrix_shrinkage.r')
source('../matrixPrior/matrixPrior_utils.R') # may have a different 'create blocky matrix' version

source('../../utils/simEBCICell_utils.R') 



make_plots_forall = FALSE
repetitions_forall = 20
# repetitions_forall = 3
matapprox_methods_forall = c('matcomp_linearreg', 'matcomp_softImpute', 'matdecomp_svd', 'matdecomp_sparsesvd', 'spectralbiclust', 'spectralbiclust_threshold', 'zeros', 'average')
cell_distns_forall = c('pois', 'nb')
parallel_forall = TRUE
# simulations


ALPHA = .1 # PARAMETER
overall_save_folder = '../../plots/simEBCICell/' # dir.create(shrinkage_plots_folder)


# === SETTING A === (small, easy, fast setting)
# set.seed(12345)
# setting_name = 'A'
# P=5
# G=10
# rank=2
# ThetaA = create_blocky_matrix(r = rank, n = P, m = G) # display_matrix(Theta)
# pi_P = runif(P, min = .3, max = .7); pi_P = pi_P / sum(pi_P) # propensity score for each treatment
# save_folder = sprintf('%s%s/', overall_save_folder, setting_name); dir.create(save_folder)
# t0 = Sys.time()
# sim_A_results = sim_EBCI_celllevel(
#   ALPHA=ALPHA,            # alpha testing level (for CIs)
#   P=P,                    # number of grnas/perturbations
#   G=G,                    # number of genes
#   rank=rank,              # true rank of theta matrix
#   Theta=ThetaA,            # true effects ( input NULL if want to simulate)
#   cell_distns = cell_distns_forall, # cell_distns (pois/nb)
#   N=500,                 # number of treated cells
#   N_control=500,          # number of control cells
#   pi_P=pi_P,              # propensity score for each treatment
#   nb_size=3,              # contributes to the overdispersion of NB (smaller = more overdispersed)
#   ranks=c(1, 2, 3),       # ranks for matrix approximations
#   matapprox_methods = matapprox_methods_forall, # matrix approx methods
#   save_folder=save_folder,# folder to save results and plots
#   make_plots = make_plots_forall, # whether to create plots (set to FALSE bc time)
#   repetitions=repetitions_forall, # number of repetitions
#   write_plot_df=FALSE,            # whether to save plot_df (set to FALSE bc this large file)
#   parallel=parallel_forall        # whether to run in parallel
# )
# 
# # N=500; N_control=500; nb_size=3; ranks=c(1,2,3)
# t1 = Sys.time(); print(t1 - t0) # 8 sec w/o plots, 1.7 mins with plots??

# 
# 
# 
# # === SETTING B === (slightly larger, ) 
# set.seed(12345)
# setting_name = 'B'
# P=10
# G=20
# rank=3
# Theta = create_blocky_matrix(r = rank, n = P, m = G) # display_matrix(Theta)
# pi_P = runif(P, min = .3, max = .7); pi_P = pi_P / sum(pi_P) # propensity score for each treatment
# save_folder = sprintf('%s%s/', overall_save_folder, setting_name); dir.create(save_folder)
# t0 = Sys.time()
# sim_A_results = sim_EBCI_celllevel(
#   ALPHA=ALPHA,            # alpha testing level (for CIs)
#   P=P,                    # number of grnas/perturbations
#   G=G,                    # number of genes
#   rank=rank,              # true rank of theta matrix
#   Theta=Theta,            # true effects ( input NULL if want to simulate)
#   N=1000,                 # number of treated cells
#   N_control=500,          # number of control cells
#   pi_P=pi_P,              # propensity score for each treatment
#   nb_size=3,              # contributes to the overdispersion of NB (smaller = more overdispersed)
#   ranks=c(1, 2, 3),       # ranks for matrix approximations
#   save_folder=save_folder,# folder to save results and plots
#   make_plots = make_plots_forall
# )
# t1 = Sys.time(); print(t1 - t0) # 1.4  mins without plots, 3.x with plots
# 
# 
# # === SETTING C === (more noise)
# set.seed(12345)
# setting_name = 'C'
# P=10
# G=20
# rank=3
# Theta = create_blocky_matrix(r = rank, n = P, m = G) # display_matrix(Theta)
# pi_P = runif(P, min = .3, max = .7); pi_P = pi_P / sum(pi_P) # propensity score for each treatment
# save_folder = sprintf('%s%s/', overall_save_folder, setting_name); dir.create(save_folder)
# t0 = Sys.time()
# sim_A_results = sim_EBCI_celllevel(
#   ALPHA=ALPHA,            # alpha testing level (for CIs)
#   P=P,                    # number of grnas/perturbations
#   G=G,                    # number of genes
#   rank=rank,              # true rank of theta matrix
#   Theta=Theta,            # true effects ( input NULL if want to simulate)
#   N= 500,                 # number of treated cells
#   N_control=250,          # number of control cells
#   pi_P=pi_P,              # propensity score for each treatment
#   nb_size=3,              # contributes to the overdispersion of NB (smaller = more overdispersed)
#   ranks=c(1, 2, 3),       # ranks for matrix approximations
#   save_folder=save_folder,# folder to save results and plots
#   make_plots = make_plots_forall
# )
# t1 = Sys.time(); print(t1 - t0) # 
# 
# 
# 
# 
# # === SETTING D === (larger matrix)
# set.seed(12345)
# setting_name = 'D'
# P=50
# G=500
# rank=3
# Theta = create_blocky_matrix(r = rank, n = P, m = G) # display_matrix(Theta)
# pi_P = runif(P, min = .3, max = .7); pi_P = pi_P / sum(pi_P) # propensity score for each treatment
# save_folder = sprintf('%s%s/', overall_save_folder, setting_name); dir.create(save_folder)
# t0 = Sys.time()
# sim_A_results = sim_EBCI_celllevel(
#   ALPHA=ALPHA,            # alpha testing level (for CIs)
#   P=P,                    # number of grnas/perturbations
#   G=G,                    # number of genes
#   rank=rank,              # true rank of theta matrix
#   Theta=Theta,            # true effects ( input NULL if want to simulate)
#   N=5000,                 # number of treated cells
#   N_control=500,          # number of control cells
#   pi_P=pi_P,              # propensity score for each treatment
#   nb_size=3,              # contributes to the overdispersion of NB (smaller = more overdispersed)
#   ranks=c(1, 2, 3),       # ranks for matrix approximations
#   save_folder=save_folder,# folder to save results and plots
#   make_plots = make_plots_forall
# )
# t1 = Sys.time(); print(t1 - t0) # 
# 


# === SETTING E === (larger matrix, larger noise by smaller sample size)
set.seed(12345)
setting_name = 'E'
P=50
G=500
rank=3
Theta = create_blocky_matrix(r = rank, n = P, m = G) # display_matrix(Theta)
pi_P = runif(P, min = .3, max = .7); pi_P = pi_P / sum(pi_P) # propensity score for each treatment
save_folder = sprintf('%s%s/', overall_save_folder, setting_name); dir.create(save_folder)
t0 = Sys.time()
sim_A_results = sim_EBCI_celllevel(
  ALPHA=ALPHA,            # alpha testing level (for CIs)
  P=P,                    # number of grnas/perturbations
  G=G,                    # number of genes
  rank=rank,              # true rank of theta matrix
  Theta=Theta,            # true effects ( input NULL if want to simulate)
  cell_distns = cell_distns_forall, # cell_distns (pois/nb)
  N=2000,                 # number of treated cells
  N_control=200,          # number of control cells
  pi_P=pi_P,              # propensity score for each treatment
  nb_size=3,              # contributes to the overdispersion of NB (smaller = more overdispersed)
  ranks=c(1, 2, 3),       # ranks for matrix approximations
  matapprox_methods = matapprox_methods_forall, # matrix approx methods
  save_folder=save_folder,# folder to save results and plots
  make_plots = make_plots_forall, # whether to create plots (set to FALSE bc time)
  repetitions=repetitions_forall, # number of repetitions
  write_plot_df=FALSE,            # whether to save plot_df (set to FALSE bc this large file)
  parallel=parallel_forall        # whether to run in parallel
)
t1 = Sys.time(); print(t1 - t0) #
# 
# 
# # === SETTING F === (larger matrix, largerer noise by smallerer sample size)
# set.seed(12345)
# setting_name = 'F'
# P=50
# G=500
# rank=3
# Theta = create_blocky_matrix(r = rank, n = P, m = G) # display_matrix(Theta)
# pi_P = runif(P, min = .3, max = .7); pi_P = pi_P / sum(pi_P) # propensity score for each treatment
# save_folder = sprintf('%s%s/', overall_save_folder, setting_name); dir.create(save_folder)
# t0 = Sys.time()
# sim_A_results = sim_EBCI_celllevel(
#   ALPHA=ALPHA,            # alpha testing level (for CIs)
#   P=P,                    # number of grnas/perturbations
#   G=G,                    # number of genes
#   rank=rank,              # true rank of theta matrix
#   Theta=Theta,            # true effects ( input NULL if want to simulate)
#   N=1000,                 # number of treated cells
#   N_control=100,          # number of control cells
#   pi_P=pi_P,              # propensity score for each treatment
#   nb_size=3,              # contributes to the overdispersion of NB (smaller = more overdispersed)
#   ranks=c(1, 2, 3),       # ranks for matrix approximations
#   save_folder=save_folder,# folder to save results and plots
#   make_plots = make_plots_forall
# )
# t1 = Sys.time(); print(t1 - t0) #

# === SETTING G === (larger noise?, rank increase)
set.seed(12345)
setting_name = 'G'
P=50
G=500
rank=5
Theta = create_blocky_matrix(r = rank, n = P, m = G) # display_matrix(Theta)
pi_P = runif(P, min = .3, max = .7); pi_P = pi_P / sum(pi_P) # propensity score for each treatment
save_folder = sprintf('%s%s/', overall_save_folder, setting_name); dir.create(save_folder)
t0 = Sys.time()
sim_A_results = sim_EBCI_celllevel(
  ALPHA=ALPHA,            # alpha testing level (for CIs)
  P=P,                    # number of grnas/perturbations
  G=G,                    # number of genes
  rank=rank,              # true rank of theta matrix
  Theta=Theta,            # true effects ( input NULL if want to simulate)
  cell_distns = cell_distns_forall, # cell_distns (pois/nb)
  N=1000,                 # number of treated cells
  N_control=500,          # number of control cells
  pi_P=pi_P,              # propensity score for each treatment
  nb_size=3,              # contributes to the overdispersion of NB (smaller = more overdispersed)
  ranks=c(1, 3, 5),       # ranks for matrix approximations
  matapprox_methods = matapprox_methods_forall, # matrix approx methods
  save_folder=save_folder,# folder to save results and plots
  make_plots = make_plots_forall, # whether to create plots (set to FALSE bc time)
  repetitions=repetitions_forall, # number of repetitions
  write_plot_df=FALSE,            # whether to save plot_df (set to FALSE bc this large file)
  parallel=parallel_forall        # whether to run in parallel
)
t1 = Sys.time(); print(t1 - t0) # 


