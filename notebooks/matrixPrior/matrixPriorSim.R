################################## Run Simulations ################################################# 

library(ggplot2)
library(latex2exp)
ggplot2::set_theme(theme_bw())
library(PMA)
library(softImpute)

source('./matrixPrior_utils.R')


save_folder = "C:/Users/Cathe/Documents/School/Genetic Hypothesis Testing using Negative Controls/multilevel-treatment-multivariable-outcome-github/multilevel-treatment-multivariable-outcome/plots/matrixPrior/sim/"




# M = create_lowrank_matrix(r=r, n=n, m=m)




# A: Easy case- large signal, small noise
r = 5
n = 15
m = 10
set.seed(12345)
M        = create_blocky_matrix(r=r, n=n, m=m)     # the signal
noise_sd = runif(n*m, min=.1, max=.5)              # sd of the noise


grobA = sim_and_plot_samplesplit(M=M, noise_sd=noise_sd, svd_maxrank=r,
                                 save_folder=sprintf('%s/A/', save_folder))

grobA_nosamplesplit = sim_and_plot_nosamplesplit(M=M, noise_sd=noise_sd/sqrt(2), svd_maxrank=r,
                                                 save_folder=sprintf('%s/A/', save_folder))


# B: Medium case- ----- signal, medium noise
r = 5
n = 15
m = 10
set.seed(12345)
M        = create_blocky_matrix(r=r, n=n, m=m)    # the signal
noise_sd = runif(n*m, min=.3, max=4)              # sd of the noise


grobB = sim_and_plot_samplesplit(M=M, noise_sd=noise_sd, svd_maxrank=r,
                     save_folder=sprintf('%s/B/', save_folder))

grobB_nosamplesplit = sim_and_plot_nosamplesplit(M=M, noise_sd=noise_sd/sqrt(2), svd_maxrank=r,
                                                 save_folder=sprintf('%s/B/', save_folder))

# C: Hard case- ----- signal, large noise
r = 5
n = 15
m = 10
set.seed(12345)
M        = create_blocky_matrix(r=r, n=n, m=m)  # the signal
noise_sd = runif(n*m, min=3, max=8)            # sd of the noise


grobC = sim_and_plot_samplesplit(M=M, noise_sd=noise_sd, svd_maxrank=r,
                     save_folder=sprintf('%s/C/', save_folder))

grobC_nosamplesplit = sim_and_plot_nosamplesplit(M=M, noise_sd=noise_sd/sqrt(2), svd_maxrank=r,
                                                 save_folder=sprintf('%s/C/', save_folder))




# D: Hard case- strong signal, large noise, large matrix
r = 5
n = 100
m = 50
set.seed(12345)
M        = create_blocky_matrix(r=r, n=n, m=m)*5  # the signal
noise_sd = runif(n*m, min=5, max=15)            # sd of the noise


grobD = sim_and_plot_samplesplit(M=M, noise_sd=noise_sd, svd_maxrank=r,
                     save_folder=sprintf('%s/D/', save_folder))

grobD_nosamplesplit = sim_and_plot_nosamplesplit(M=M, noise_sd=noise_sd/sqrt(2), svd_maxrank=r,
                                                 save_folder=sprintf('%s/D/', save_folder))

# E: Hard case- strong signal, larger noise, large matrix
r = 5
n = 100
m = 50
set.seed(12345)
M        = create_blocky_matrix(r=r, n=n, m=m)*5  # the signal
noise_sd = runif(n*m, min=30, max=50)             # sd of the noise


grobE = sim_and_plot_samplesplit(M=M, noise_sd=noise_sd, svd_maxrank=r,
                     save_folder=sprintf('%s/E/', save_folder))

grobE_nosamplesplit = sim_and_plot_nosamplesplit(M=M, noise_sd=noise_sd/sqrt(2), svd_maxrank=r,
                                                 save_folder=sprintf('%s/E/', save_folder))



















################################## Test Matrix Completion ################################################# 
# test our implementation, hich is simplified bc e only have 1 element missing, 
# so e use some linear regression procedure


r = 5
n = 15
m = 10
set.seed(12345)
# M = create_lowrank_matrix(r=r, n=n, m=m)
M = create_blocky_matrix(r=r, n=n, m=m) 
svd_maxrank = r
noise_sd_s = seq(from = 0, to = .5, length.out = 5)


M_obs_plots = list()
noise_label_plots = list()
X_mcomp_list = list()
X_mcomp_plots = list()
X_svd_plots = list()
X_ssvd_plots = list()

for(noise_sd_s_ in noise_sd_s) {
  print(noise_sd_s_)
  noise_sd = rep(noise_sd_s_, times = n*m) 
  M_obs    = M + matrix(rnorm(n = length(noise_sd), mean = 0, sd = noise_sd),
                     nrow = n, ncol = m, byrow = TRUE) # display_matrix(M_obs) 
  
  X_mcomp_list[[as.character(noise_sd_s_)]] = matrix(NA, nrow = n, ncol = m)
  noise_label_plots[[as.character(noise_sd_s_)]] <- local({
    noise_sd_s_ <- noise_sd_s_
    ggplot() + geom_text(aes(x=0, y=0, label = sprintf('noise sd:\n%.2f', noise_sd_s_)), size = 15) + theme_void()
  })
  # noise_sd_s_text = sprintf('noise sd:\n%.2f', noise_sd_s_)
  # noise_label_plots[[as.character(noise_sd_s_)]] = ggplot() + geom_text(aes(x=0, y=0, label = noise_sd_s_text), size = 15) + theme_void()
  for(i in 1:n) {
    for(j in 1:m) {
      # LR: X is linearreg
      X_mcomp_list[[as.character(noise_sd_s_)]][i,j] = linearreg(M_obs, i, j)
    }
  }

  M_obs_plots[[as.character(noise_sd_s_)]] = display_matrix(M_obs)
  X_mcomp_plots[[as.character(noise_sd_s_)]] = display_matrix(X_mcomp_list[[as.character(noise_sd_s_)]])
  
  # === Construct X matrices from decomp (svd, sparse svd, etc...)
  svdres = svd(M_obs, nu = svd_maxrank, nv = svd_maxrank)
  X_svd  = svdres$u %*% diag(svdres$d[1:svd_maxrank]) %*% t(svdres$v)
  X_svd_plots[[as.character(noise_sd_s_)]] = display_matrix(X_svd)
  
  
  cv.out = PMA::SPC.cv(M_obs, sumabsvs = seq(1.2, min(5, sqrt(n), sqrt(m)), len = 10))
  pmd_res = PMA::SPC(M_obs, sumabsv=cv.out$bestsumabsv, K = svd_maxrank)
  X_sparse_svd = pmd_res$u %*% diag(pmd_res$d) %*% t(pmd_res$v)
  X_ssvd_plots[[as.character(noise_sd_s_)]] = display_matrix(X_sparse_svd)
  
  rm(noise_sd_s_)
}



# plot labels
p_label_Mobs  = ggplot() + geom_text(aes(x=0, y=0, label = 'Observed\nMatrix'), size = 15) + theme_void()
p_label_Mcomp = ggplot() + geom_text(aes(x=0, y=0, label = 'Matrix\nCompletion'), size = 15) + theme_void()
p_label_SVD = ggplot() + geom_text(aes(x=0, y=0, label = 'SVD'), size = 15) + theme_void()
p_label_SSVD = ggplot() + geom_text(aes(x=0, y=0, label = 'Sparse SVD'), size = 15) + theme_void()


layout_matrix = cbind( matrix(c(NA, 5*length(noise_sd_s)+(1:4)),
                              nrow = 5),
                       matrix(1:(5*length(noise_sd_s)), 
                              nrow = 5, byrow = TRUE))

grob <- gridExtra::arrangeGrob(grobs = c(noise_label_plots, 
                                         M_obs_plots, 
                                         X_mcomp_plots, 
                                         X_svd_plots,
                                         X_ssvd_plots,
                                         p_label_Mobs, p_label_Mcomp, p_label_SVD, p_label_SSVD),
                               layout_matrix = layout_matrix)



gridExtra::grid.arrange(grob)

ggsave(sprintf('%s/test_matcomp_linearreg.pdf', save_folder), grob, width = 24, height = 16)



























################################## Trash ###########################################################


if(F) {
  # Testing out softImput package
  library(softImpute)
  
  
  # example code from  https://cran.r-project.org/web/packages/softImpute/refman/softImpute.html#softImpute-package
  set.seed(101)
  n=200
  p=100
  J=50
  np=n*p
  missfrac=0.3
  x=matrix(rnorm(n*J),n,J)%*%matrix(rnorm(J*p),J,p)+matrix(rnorm(np),n,p)/5
  ix=seq(np)
  imiss=sample(ix,np*missfrac,replace=FALSE)
  xna=x
  xna[imiss]=NA
  xnaC=as(xna,"Incomplete")
  ### here we do it a different way to demonstrate Incomplete
  ### In practise the observed values are stored in this market-matrix format.
  i = row(xna)[-imiss]
  j = col(xna)[-imiss]
  xnaC=Incomplete(i,j,x=x[-imiss])
  # why is this such a weird input format
  
  
  
  
  # vignette
  # https://cran.r-project.org/web/packages/softImpute/vignettes/softImpute.html
  set.seed(1011)
  x=matrix(rnorm(30),6,5)
  x[sample(1:30,10,replace=FALSE)]=NA
  x
  
  fits=softImpute(x,trace=TRUE,type="svd")
  fits
  fits2=softImpute(x,rank.max=3,lambda=1.9,trace=TRUE,type="svd")
  complete(x,fits2)
  
  
  xc=biScale(x,col.scale=FALSE,row.scale=FALSE,trace=TRUE)
  fits3=softImpute(xc,rank.max=3,lambda=1,type="svd")
  fits3$d
  complete(x,fits3,unscale=TRUE)
  
  
  
  
  # choosing lambda parameter
  xs=as(x,"Incomplete")
  lam0=lambda0(xs)
  lam0
  
  fit0=softImpute(xs,lambda=lam0+.2)
  fit0$d
  lamseq=exp(seq(from=log(lam0),to=log(1),length=10))
  lamseq
  
  fits=as.list(lamseq)
  ranks=as.integer(lamseq)
  rank.max=2
  warm=NULL
  for( i in seq(along=lamseq)){
    fiti=softImpute(xs,lambda=lamseq[i],rank=rank.max,warm=warm)
    ranks[i]=sum(round(fiti$d,4)>0)
    rank.max=min(ranks[i]+2,4)
    warm=fiti
    fits[[i]]=fiti
    cat(i,"lambda=",lamseq[i],"rank.max",rank.max,"rank",ranks[i],"\n")
  }
  
  
  
  # min code for running
  library(softImpute)
  set.seed(1011)
  x=matrix(rnorm(30),6,5)
  x[sample(1:30,10,replace=FALSE)]=NA
  x
  
  # choosing lambda parameter
  xs=as(x,"Incomplete")
  lam0=lambda0(xs)
  lam0

  lamseq=exp(seq(from=log(lam0),to=log(1),length=10))
  lamseq
  
  fits=as.list(lamseq)
  ranks=as.integer(lamseq)
  rank.max=2
  warm=NULL
  for( i in seq(along=lamseq)){
    fiti=softImpute::softImpute(xs,lambda=lamseq[i],rank=rank.max,warm=warm)
    ranks[i]=sum(round(fiti$d,4)>0)
    rank.max=min(ranks[i]+2,4)
    warm=fiti
    fits[[i]]=fiti
    cat(i,"lambda=",lamseq[i],"rank.max",rank.max,"rank",ranks[i],"\n")
  }
  
  xc = softImpute::complete(x,fiti,unscale=TRUE)
  
  
  # to speed up computation, we can try to reuse the same parameters rather than retry every time
  
  #' Perform Matrix Completion for leave one out X
  #' @param X (matrix)
  #' @param type (character) 'svd' or 'als' for the type of softImpute method used
  create_matrix_completion_softImpute <- function(X, rank_max, type = 'svd') {
    # debug
    # X=matrix(rnorm(30),6,5)
    # X[sample(1:30,10,replace=FALSE)]=NA
    
    
    
    
    
    # Part 1: Set up hyperparameters/initial start to reuse ?? 
    # choosing lambda parameter ??? idk, just do (1/2) * lam0 for now??
    X_dgCMatrix=as(X,"Incomplete")
    lam0=lambda0(X_dgCMatrix)
    # 
    # lamseq=exp(seq(from=log(lam0),to=log(1),length=10))
    # 
    # fits=as.list(lamseq)
    # ranks=as.integer(lamseq)
    # rank.max=2
    # warm=NULL
    # nuclear_norms = as.list(lamseq)
    # matrix_completed = as.list(lamseq)
    # for( i in seq(along=lamseq)){
    #   fiti=softImpute::softImpute(X_dgCMatrix,lambda=lamseq[i],rank=rank.max,warm=warm)
    #   ranks[i]=sum(round(fiti$d,4)>0)
    #   rank.max=min(ranks[i]+2,4)
    #   warm=fiti
    #   fits[[i]]=fiti
    #   # cat(i,"lambda=",lamseq[i],"rank.max",rank.max,"rank",ranks[i],"\n")
    #   matrix_completed[i] = softImpute::complete(X,fiti,unscale=TRUE)
    #   nuclear_norms[i] = sum(fiti$d)
    # }
    # 
    # # complete matrix
    # X_completed = softImpute::complete(X,fiti,unscale=TRUE)
    # X_completed = matrix(NA, nrow=nrow(X), ncol=ncol(X)) 
    # for(i in 1:nrow(X)) {
    #   for(j in 1:ncol(X)) {
    #     
    #   }
    # }
    # 
    
    
    
    
    matrix_completion_softImpute <- function(i, j) {
      X_ij = X
      X_ij[i,j] = NA
      # use the same fit for i,j element complete (used for sample splitting version)
      fit = softImpute::softImpute(X_ij, lambda = (lam0 + 1)/2, rank=rank_max, warm=NULL, type = type)
      
      
      X_completed = softImpute::complete(X_ij, fit, unscale=TRUE)
      
      return(X_completed[i,j])
    }
    
    return(matrix_completion_softImpute)
    
   
  }
  
  
  # naive 
  r = 5
  n = 15
  m = 10
  set.seed(12345)
  # M = create_lowrank_matrix(r=r, n=n, m=m)
  M = create_blocky_matrix(r=r, n=n, m=m)
  noise_sd = runif(n*m, min=.1, max=(max(M) - min(M)) / 12) * (1/sqrt(2))
  M_obs = M + matrix(rnorm(n = length(noise_sd), mean = 0, sd = noise_sd),
                     nrow = n, ncol = m, byrow = TRUE) # display_matrix(M_obs) 
  svd_maxrank = r
  
  matcomp_fun_M    = create_matrix_completion_softImpute(X = M,     rank_max = svd_maxrank, type = 'als')
  matcomp_fun_Mobs = create_matrix_completion_softImpute(X = M_obs, rank_max = svd_maxrank, type = 'als')
  
  M_comp     = matrix(NA, nrow = n, ncol = m)
  M_obs_comp = matrix(NA, nrow = n, ncol = m)
  
  for(i in 1:nrow(M)) {
    for(j in 1:ncol(M)) {
      M_comp[i,j]     = matcomp_fun_M(i,j)
      M_obs_comp[i,j] = matcomp_fun_Mobs(i,j)
    }
  }
  
  M_comp
  
  M_comp - M
  display_matrix( M)
  display_matrix((M_comp - M)/pmax(M, 1))
  
  display_matrix( M_obs)
  display_matrix(M_obs_comp)
  display_matrix((M_obs_comp - M)/pmax(M, 1))
  
  
  

}


if(F) {
  set.seed(12345)
  
  # construct low rank matrix 
  r = 5
  n = 20
  m = 10
  U = matrix(rnorm(n = r*n, mean = 0, sd = 1), 
             nrow = n)
  V = matrix(rnorm(n = r*m, mean = 0, sd = 1), 
             nrow = m)
  Sigma = diag(runif(n = r, min = .1, max = 5))
  M = U %*% Sigma %*% t(V)
}





# old- one run of testing, easy case
if(F) {
  library(ggplot2)
  library(latex2exp)
  ggplot2::set_theme(theme_bw())
  r = 5
  n = 15
  m = 10
  set.seed(12345)
  # M = create_lowrank_matrix(r=r, n=n, m=m)
  M = create_blocky_matrix(r=r, n=n, m=m) * 3
  
  noise_sd = runif(n*m, min=.1, max=.5)
  noise_sd_matrix = matrix(noise_sd, nrow = n, ncol = m, byrow = TRUE)
  
  M_split1 = M + matrix(rnorm(n = length(noise_sd), mean = 0, sd = noise_sd),
                        nrow = n, ncol = m, byrow = TRUE)
  M_split2 = M + matrix(rnorm(n = length(noise_sd), mean = 0, sd = noise_sd),
                        nrow = n, ncol = m, byrow = TRUE)
  
  library(PMA)
  # should tune parameters using cv
  # pmd_res = PMA::PMD.cv(M_split1, type = 'standard')
  pmd_res = PMA::PMD(x=M_split1, type = 'standard')
  
  X_sparse_svd = pmd_res$u%*%t(pmd_res$v)
  
  
  
  svdres = svd(M_split1, nu = maxrank, nv = maxrank)
  X_svd  = svdres$u %*% diag(svdres$d[1:maxrank]) %*% t(svdres$v)
  
  
  X      = matrix(NA, nrow = n, ncol = m) 
  tstats = matrix(NA, nrow = n, ncol = m)
  pvals  = matrix(NA, nrow = n, ncol = m)
  pvals_standard = matrix(NA, nrow = n, ncol = m)
  pvals_svd   = matrix(NA, nrow = n, ncol = m)
  for(i in 1:n) {
    for(j in 1:m) {
      
      # X is custom fill fn
      # X[i,j] = average(M_split1, i, j)
      X[i,j] = linearreg(M_split1, i, j)
      tstat_and_pval = calc_tstat_and_pval(
        Yij      = M_split2[i,j], 
        sij      = noise_sd_matrix[i,j], 
        thetaij  = X[i,j], 
        thetaij0 = 0
      )
      tstats[i,j] = tstat_and_pval$tstat
      pvals[i,j]  = tstat_and_pval$pval
      
      
      # X is SVD
      tstat_and_pval = calc_tstat_and_pval(
        Yij      = M_split2[i,j], 
        sij      = noise_sd_matrix[i,j], 
        thetaij  = X_svd[i,j], 
        thetaij0 = 0
      )
      pvals_svd[i,j] = tstat_and_pval$pval
      
      
      
      # Standard: 
      pvals_standard[i,j] = 2*(1 - pnorm(abs(M_split2[i,j] - 0) /noise_sd_matrix[i,j]))
    }
  }
  
  # X_df = data.frame(i = rep())
  # mapply(FUN = average, i=1:5, j = 1:3, MoreArgs = list(M=M))
  
  
  
  p_M = display_matrix(M) + labs(title = 'M')
  p_M1 = display_matrix(M_split1) + labs(title = 'M Split 1')
  p_M2 = display_matrix(M_split2) + labs(title = 'M Split 2')
  
  
  p_Xsvd = display_matrix(X_svd) + labs(title = TeX(r'($\tilde{\Theta}$ (SVD))'))
  
  # display_matrix(X)
  # display_matrix(tstats)
  # display_matrix(pvals) + scale_fill_gradientn(colours = c('red', 'orange', 'white'), values = c(0, .05, 1))
  
  p_pval_Xsvd = display_matrix(pvals_svd) + scale_fill_gradientn(colours = c('red', 'orange', 'white'), values = c(0, .05, 1)) + labs(title = TeX(r'($\tilde{\Theta}$ (SVD) p-vals)'))
  p_pval_standard = display_matrix(pvals_standard) + scale_fill_gradientn(colours = c('red', 'orange', 'white'), values = c(0, .05, 1)) + labs(title = 'Standard p-vals')
  
  
  
  
  grob <- gridExtra::arrangeGrob(p_M, p_M1, p_M2,
                                 p_Xsvd, 
                                 p_pval_Xsvd, p_pval_standard,
                                 layout_matrix = matrix(c(1, 2, 3,
                                                          4, 5, 6), byrow = TRUE, nrow = 2))
  # ggsave(sprintf('%s/all.pdf', save_folder), grob, width = 24, height = 18)
  gridExtra::grid.arrange(grob)
}
