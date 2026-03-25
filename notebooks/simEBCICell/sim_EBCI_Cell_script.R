






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

# our code
source('../../utils/matrix_shrinkage.r')
source('../matrixPrior/matrixPrior_utils.R') # may have a different 'create blocky matrix' version

source('../../utils/simEBCICell_utils.R') 



# simulations


ALPHA = .1 # PARAMETER
overall_save_folder = '../../plots/simEBCICell/' # dir.create(shrinkage_plots_folder)


# === SETTING A === (small, easy, fast setting)
set.seed(12345)
setting_name = 'A'
P=5
G=10
rank=2
ThetaA = create_blocky_matrix(r = rank, n = P, m = G) # display_matrix(Theta)
pi_P = runif(P, min = .3, max = .7); pi_P = pi_P / sum(pi_P) # propensity score for each treatment
save_folder = sprintf('%s%s/', overall_save_folder, setting_name); dir.create(save_folder)
t0 = Sys.time()
sim_A_results = sim_EBCI_celllevel(
  ALPHA=ALPHA,            # alpha testing level (for CIs)
  P=P,                    # number of grnas/perturbations
  G=G,                    # number of genes
  rank=rank,              # true rank of theta matrix
  Theta=ThetaA,            # true effects ( input NULL if want to simulate)
  N=500,                 # number of treated cells
  N_control=500,          # number of control cells
  pi_P=pi_P,              # propensity score for each treatment
  nb_size=3,              # contributes to the overdispersion of NB (smaller = more overdispersed)
  ranks=c(1, 2, 3),       # ranks for matrix approximations
  save_folder=save_folder,# folder to save results and plots
  make_plots = TRUE
)

# N=500; N_control=500; nb_size=3; ranks=c(1,2,3)
t1 = Sys.time(); print(t1 - t0) # 8 sec w/o plots, 1.7 mins with plots??




# === SETTING B === (slightly larger, ) 
set.seed(12345)
setting_name = 'B'
P=10
G=20
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
  N=1000,                 # number of treated cells
  N_control=500,          # number of control cells
  pi_P=pi_P,              # propensity score for each treatment
  nb_size=3,              # contributes to the overdispersion of NB (smaller = more overdispersed)
  ranks=c(1, 2, 3),       # ranks for matrix approximations
  save_folder=save_folder,# folder to save results and plots
  make_plots = TRUE
)
t1 = Sys.time(); print(t1 - t0) # 1.4  mins without plots, 3.x with plots


# === SETTING C === (even larger, )
set.seed(12345)
setting_name = 'C'
P=25
G=200
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
  N=2500,                 # number of treated cells
  N_control=500,          # number of control cells
  pi_P=pi_P,              # propensity score for each treatment
  nb_size=3,              # contributes to the overdispersion of NB (smaller = more overdispersed)
  ranks=c(1, 2, 3),       # ranks for matrix approximations
  save_folder=save_folder,# folder to save results and plots
  make_plots = TRUE
)
t1 = Sys.time(); print(t1 - t0) # 




# === SETTING D === (even larger? )
set.seed(12345)
setting_name = 'D'
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
  N=5000,                 # number of treated cells
  N_control=500,          # number of control cells
  pi_P=pi_P,              # propensity score for each treatment
  nb_size=3,              # contributes to the overdispersion of NB (smaller = more overdispersed)
  ranks=c(1, 2, 3),       # ranks for matrix approximations
  save_folder=save_folder,# folder to save results and plots
  make_plots = TRUE
)
t1 = Sys.time(); print(t1 - t0) # 



# === SETTING E === (smaller sample size (larger noise))
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
  N=2000,                 # number of treated cells
  N_control=200,          # number of control cells
  pi_P=pi_P,              # propensity score for each treatment
  nb_size=3,              # contributes to the overdispersion of NB (smaller = more overdispersed)
  ranks=c(1, 2, 3),       # ranks for matrix approximations
  save_folder=save_folder,# folder to save results and plots
  make_plots = TRUE
)
t1 = Sys.time(); print(t1 - t0) #




# === SETTING F === (rank increase)
set.seed(12345)
setting_name = 'F'
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
  N=5000,                 # number of treated cells
  N_control=500,          # number of control cells
  pi_P=pi_P,              # propensity score for each treatment
  nb_size=3,              # contributes to the overdispersion of NB (smaller = more overdispersed)
  ranks=c(1, 3, 5),       # ranks for matrix approximations
  save_folder=save_folder,# folder to save results and plots
  make_plots = TRUE
)
t1 = Sys.time(); print(t1 - t0) # 


