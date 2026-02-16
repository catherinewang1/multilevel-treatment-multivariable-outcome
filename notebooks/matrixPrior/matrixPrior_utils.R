






#' creates a nxm matrix of low rank r
#' @param r (integer) rank
#' @param n (integer) number of rows
#' @param m (integer) number of cols
create_lowrank_matrix <- function(r, n, m) {
  # this is NOT the svd, but it helps create a low rank matrix
  U = matrix(rnorm(n = r*n, mean = 0, sd = 1), 
             nrow = n)
  V = matrix(rnorm(n = r*m, mean = 0, sd = 1), 
             nrow = m)
  Sigma = diag(runif(n = r, min = .1, max = 5))
  M = U %*% Sigma %*% t(V)
  
  return(M)
}


#' creates a nxm matrix of low rank r that is blocky on the diagonal
#' 
#' @param r (integer) rank
#' @param n (integer) number of rows
#' @param m (integer) number of cols
create_blocky_matrix <- function(r, n, m) {

  
  U = matrix(0, nrow = n, ncol = r)
  V = matrix(0, nrow = m, ncol = r)
  for(r_ in 1:r) {
    
    u_idx = ((n/r)*(r_-1) + 1):((n/r)*r_) # some rounding
    v_idx = ((m/r)*(r_-1) + 1):((m/r)*r_)
    # make blocks of the same sign
    sign_ = 2 * rbinom(n = 1, size = 1, prob = .5) - 1
    U[u_idx, r_] = sign_ * abs(rnorm(n = length(u_idx), mean = 2/sqrt(n/r), sd = 1/sqrt(n/r)))
    V[v_idx, r_] = abs(rnorm(n = length(v_idx), mean = 2/sqrt(m/r), sd = 1/sqrt(n/r)))
  }
  
  Sigma = diag(runif(n = r, min = 1, max = 5))
  
  M = U %*% Sigma %*% t(V)
  
  return(M)
}







average <- function(M, i, j) {
  # take the average of the rows
  row_avg = mean(M[i, -j])
  # take the average of the cols
  col_avg = mean(M[i, -j])
  
  return(mean(row_avg, col_avg))
}

linearreg <- function(M, i, j) {
  M_tilde = M[-i, -j] # M without the ith row and jth col
  
  y = M[-i, j]
  
  lmfit = lm(y ~ M_tilde + 0)
  
  summary(lmfit)
  
  coefs = lmfit$coefficients 
  coefs[is.na(coefs)] = 0 
  Xij = sum(coefs * M[i, -j])
  
  # if exact, Xij == M[i,j]. But in practice, not exact
  
  return(Xij)
}

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
  
  
  
  #' @param i (integer) row index
  #' @param j (integer) col index
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



calc_tstat_and_pval <- function(Yij, sij, thetaij, thetaij0) {
  sign_thetaij_minus_thetaij0 = (2*as.integer(thetaij >= thetaij0) - 1) # 1 if thetaij - thetaij0 >= 0, -1 otherwise
                                                                       # replaces sign(thetaij - thetaij0) which will give 0 if equal
  # calculate test statistic
  tstat = sign_thetaij_minus_thetaij0 * (Yij - (1/2) * (thetaij + thetaij0))
  # calculate the p-value (tstat - theoretical mean under null) / theoretical sd under null (then z stat N(0,1))
  zstat = sign_thetaij_minus_thetaij0 * (Yij - thetaij0) / sij
  pval = 1 - pnorm(zstat) # 1-sided, where large positive val --> reject
  return(list(tstat=tstat,
              pval=pval))
}


display_matrix <- function(X) {
  ggplot(data=reshape2::melt(X, c("x", "y"), value.name = "z"),aes(x=x,y=y,fill=z)) + 
    scale_x_continuous(expand=c(0, 0)) +
    scale_y_continuous(expand=c(0, 0), transform = 'reverse') +
    scale_fill_gradient2(midpoint=0) + 
    geom_tile()
}



#' @param P matrix of pvals to plot
#' @param M true params (need to know which are actually 0) or matrix where 0s are in the right pos
plot_pvals <- function(P, M) {
  # P = pvals_svd
  # P = pvals_standard
  dfP = reshape2::melt(P, c("x", "y"), value.name = "pval")
  dfM = reshape2::melt(M, c("x", "y"), value.name = "theta")
  df = merge(dfP, dfM, by = c('x', 'y'))
  
  # null (theta = 0)
  p_null = ggplot(df |> dplyr::filter(theta == 0) |> 
           dplyr::arrange(pval) |> 
           dplyr::mutate(rej_rate = (1:dplyr::n())/dplyr::n()),
         aes(x = pval, y = rej_rate)) + 
    geom_abline(aes(slope = 1, intercept = 0)) +
    geom_point(alpha = .9, size = .6) +
    # scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    labs(title = 'Rejection Rate by p-value (qqplot)',
         subtitle = '(when Theta=0)',
         x = 'p-value', y = 'Rejection Rate')
  
  # alt (theta != 0)
  p_alt = ggplot(df |> dplyr::filter(theta != 0) |> 
           dplyr::arrange(pval) |> 
           dplyr::mutate(rej_rate = (1:dplyr::n())/dplyr::n()),
         aes(x = pval, y = rej_rate)) + 
    geom_abline(aes(slope = 1, intercept = 0)) +
    geom_point(alpha = .9, size = .6) +
    # scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    labs(title = 'Rejection Rate by p-value (qqplot)',
         subtitle = '(when Theta!=0)',
         x = 'p-value', y = 'Rejection Rate')
  
  
  return(list(p_null=p_null, p_alt=p_alt))
}


#' @param P list of matrices of pvals to plot, named
#' @param M true params (need to know which are actually 0) or matrix where 0s are in the right pos
plot_pvals_many <- function(Ps, M, method_colors=NULL) {
  # P = pvals_svd
  # P = pvals_standard
  
  df = NULL
  for(p_name in names(Ps)) {
    P = Ps[[p_name]]
    dfP = reshape2::melt(P, c("x", "y"), value.name = "pval")
    dfM = reshape2::melt(M, c("x", "y"), value.name = "theta")
    df_ = merge(dfP, dfM, by = c('x', 'y'))
    df_$method = p_name
    df = rbind(df, df_)
    
  }
  
  df$method = factor(df$method, levels = names(Ps))

  
  # null (theta = 0)
  p_null = ggplot(df |> dplyr::filter(theta == 0) |> 
                    dplyr::group_by(method) |>
                    dplyr::arrange(pval) |> 
                    dplyr::mutate(rej_rate = (1:dplyr::n())/dplyr::n()),
                  aes(x = pval, y = rej_rate, group = method, color = method)) + 
    geom_abline(aes(slope = 1, intercept = 0)) +
    geom_point(alpha = .8, size = .6, key_glyph = "rect") +
    geom_line() +
    # scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_color_discrete(palette = method_colors) + # NULL uses default colors
    labs(title = 'Rejection Rate by p-value (qqplot)',
         subtitle = '(when Theta=0)',
         x = 'p-value', y = 'Rejection Rate')
  
  # alt (theta != 0)
  p_alt = ggplot(df |> dplyr::filter(theta != 0) |> 
                   dplyr::group_by(method) |>
                   dplyr::arrange(pval) |> 
                   dplyr::mutate(rej_rate = (1:dplyr::n())/dplyr::n()),
                 aes(x = pval, y = rej_rate, group = method, color = method)) + 
    geom_abline(aes(slope = 1, intercept = 0)) +
    geom_point(alpha = .8, size = .6, key_glyph = "rect") +
    geom_line() +
    # scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_color_discrete(palette = method_colors) + # NULL uses default colors
    labs(title = 'Rejection Rate by p-value (qqplot)',
         subtitle = '(when Theta!=0)',
         x = 'p-value', y = 'Rejection Rate')
  
  
  return(list(p_null=p_null, p_alt=p_alt))
}


#' Simulate values and perform tests
#' to compare w sample split version, make sure to change noise_sd matrix by decreasing by 1/\sqrt{2}
#'
#' Saves the plots in  '<save_folder>/all_nosamplesplit.pdf' and returns the plot grob
#' @param M (matrix)
#' @param noise_sd (vector) of sd's for noise of M
#' @param svd_maxrank (integer) positive integer indicating rank for svd decompositions
#' @param save_folder (character) path to folder to save plots at '<save_folder>/all_nosamplesplit.pdf'
#' @return plots as grob (use gridExtra::grid.arrange(grob) to view)
sim_and_plot_nosamplesplit <- function(M, noise_sd, svd_maxrank, save_folder) {
  # # debug
  # r = 5
  # n = 15
  # m = 10
  # set.seed(12345)
  # # M = create_lowrank_matrix(r=r, n=n, m=m)
  # M = create_blocky_matrix(r=r, n=n, m=m)
  # # noise_sd = runif(n*m, min=.1, max=.5) * (1/sqrt(2))
  # noise_sd = runif(n*m, min=.1, max=(max(M) - min(M)) / 12) * (1/sqrt(2))
  # svd_maxrank = 5

  dir.create(save_folder, recursive = FALSE, showWarnings = FALSE)
  n = nrow(M)
  m = ncol(M)
  
  # === Create Observed Matrix M_obs
  noise_sd_matrix = matrix(noise_sd, nrow = n, ncol = m, byrow = TRUE)
  
  M_obs = M + matrix(rnorm(n = length(noise_sd), mean = 0, sd = noise_sd),
                        nrow = n, ncol = m, byrow = TRUE) # display_matrix(M_obs) 
  
  # === Construct X matrices (svd, sparse svd, etc...)

  # svd
  svdres = svd(M_obs, nu = svd_maxrank, nv = svd_maxrank)
  X_svd  = svdres$u %*% diag(svdres$d[1:svd_maxrank]) %*% t(svdres$v)
  
  
  # sparse svd
  cv.out = PMA::SPC.cv(M_obs, sumabsvs = seq(1.2, min(5, sqrt(n), sqrt(m)), len = 10))
  pmd_res = PMA::SPC(M_obs, sumabsv=cv.out$bestsumabsv, K = svd_maxrank)
  X_sparse_svd = pmd_res$u %*% diag(pmd_res$d) %*% t(pmd_res$v)

  #  single element matrix completion methods

  X_mcomp_linearreg  = matrix(NA, nrow = n, ncol = m) 
  X_mcomp_softImpute = matrix(NA, nrow = n, ncol = m) 
  matcomp_fun = create_matrix_completion_softImpute(X = M_obs, rank_max = svd_maxrank, type = 'als')
  for(i in 1:n) {
    for(j in 1:m) {
      X_mcomp_linearreg[i,j]  = linearreg(M_obs, i, j)
      X_mcomp_softImpute[i,j] = matcomp_fun(i, j)
    }
  }
  
  # === Calculate p-values
  pvals_standard            = matrix(NA, nrow = n, ncol = m) # standard method
  pvals_true                = matrix(NA, nrow = n, ncol = m) # LR Method: use TRUE matrix (unrealistic setting, a sanity check)
  pvals_svd                 = matrix(NA, nrow = n, ncol = m) # LR Method: use matrix from svd
  pvals_ssvd                = matrix(NA, nrow = n, ncol = m) # LR Method: use matrix from sparse svd
  pvals_mcomp_linearreg     = matrix(NA, nrow = n, ncol = m) # LR Method: use matrix from single element matrix completion linearreg
  pvals_mcomp_softImpute    = matrix(NA, nrow = n, ncol = m) # LR Method: use matrix from single element matrix completion softImpute
  
  for(i in 1:n) {
    for(j in 1:m) {
      
      # Standard: 
      pvals_standard[i,j] = 2 * (1 - pnorm(abs(M_obs[i,j]  - 0) /noise_sd_matrix[i,j]  ))
      
      # LR: X is true M
      tstat_and_pval = calc_tstat_and_pval(
        Yij      = M_obs[i,j], 
        sij      = noise_sd_matrix[i,j], 
        thetaij  = M[i,j], 
        thetaij0 = 0
      )
      # pvals_true[i,j] = runif(n=1, min = .001, max = .999) # test code
      pvals_true[i,j] = tstat_and_pval$pval ; rm(tstat_and_pval)
      
      # LR: X is SVD
      tstat_and_pval = calc_tstat_and_pval(
        Yij      = M_obs[i,j], 
        sij      = noise_sd_matrix[i,j], 
        thetaij  = X_svd[i,j], 
        thetaij0 = 0
      )
      pvals_svd[i,j] = tstat_and_pval$pval ; rm(tstat_and_pval)
      
      # LR: X is Sparse SVD
      tstat_and_pval = calc_tstat_and_pval(
        Yij      = M_obs[i,j], 
        sij      = noise_sd_matrix[i,j], 
        thetaij  = X_sparse_svd[i,j], 
        thetaij0 = 0
      )
      pvals_ssvd[i,j] = tstat_and_pval$pval ; rm(tstat_and_pval)
      
      # LR: X is linearreg
      tstat_and_pval = calc_tstat_and_pval(
        Yij      = M_obs[i,j],
        sij      = noise_sd_matrix[i,j],
        thetaij  = X_mcomp_linearreg[i,j],
        thetaij0 = 0
      )
      pvals_mcomp_linearreg[i,j]  = tstat_and_pval$pval ; rm(tstat_and_pval)
      
      # LR: X is softImpute
      tstat_and_pval = calc_tstat_and_pval(
        Yij      = M_obs[i,j],
        sij      = noise_sd_matrix[i,j],
        thetaij  = X_mcomp_softImpute[i,j],
        thetaij0 = 0
      )
      pvals_mcomp_softImpute[i,j]  = tstat_and_pval$pval ; rm(tstat_and_pval)
    }
  }
  
  
  
  p_M    = display_matrix(M) + labs(title = 'M')
  p_Mobs = display_matrix(M_obs) + labs(title = 'M Observed')
  
  
  p_Xsvd  = display_matrix(X_svd) + labs(title = TeX(r'($\tilde{\Theta}$ (SVD))'))
  p_Xssvd = display_matrix(X_sparse_svd) + labs(title = TeX(r'($\tilde{\Theta}$ (Sparse SVD))'))
  p_Xmcomp_linearreg  = display_matrix(X_mcomp_linearreg ) + labs(title = TeX(r'($\tilde{\Theta}$ (Mat. Comp.- Linear Reg))'))
  p_Xmcomp_softImpute = display_matrix(X_mcomp_softImpute) + labs(title = TeX(r'($\tilde{\Theta}$ (Mat. Comp.- Soft Impute))'))
  # display_matrix(X)
  # display_matrix(tstats)
  # display_matrix(pvals) + scale_fill_gradientn(colours = c('red', 'orange', 'white'), values = c(0, .05, 1))
  
  p_pval_standard = display_matrix(pvals_standard) + scale_fill_gradientn(colours = c('red', 'orange', 'white'), values = c(0, .05, 1)) + labs(title = 'Standard p-vals')
  p_pval_Xsvd     = display_matrix(pvals_svd)  + scale_fill_gradientn(colours = c('red', 'orange', 'white'), values = c(0, .05, 1)) + labs(title = TeX(r'($\tilde{\Theta}$ (SVD) p-vals)'))
  p_pval_Xssvd    = display_matrix(pvals_ssvd) + scale_fill_gradientn(colours = c('red', 'orange', 'white'), values = c(0, .05, 1)) + labs(title = TeX(r'($\tilde{\Theta}$ (Sparse SVD) p-vals)'))
  p_pval_Xmcomp_linearreg  = display_matrix(pvals_mcomp_linearreg ) + scale_fill_gradientn(colours = c('red', 'orange', 'white'), values = c(0, .05, 1)) + labs(title = TeX(r'($\tilde{\Theta}$ (Mat. Comp.- Linear Reg) p-vals)'))
  p_pval_Xmcomp_softImpute = display_matrix(pvals_mcomp_softImpute) + scale_fill_gradientn(colours = c('red', 'orange', 'white'), values = c(0, .05, 1)) + labs(title = TeX(r'($\tilde{\Theta}$ (Mat. Comp.- Soft Impute) p-vals)'))

  
  p_pvals_svd_qqplots       = plot_pvals(P = pvals_svd      , M = M)
  p_pvals_ssvd_qqplots      = plot_pvals(P = pvals_ssvd     , M = M)
  p_pvals_standard_qqplots  = plot_pvals(P = pvals_standard , M = M)
  
  
  # p_pvals_all_qqplots = plot_pvals_many(P = list('SVD'        = pvals_svd,
  #                                                'Sparse SVD' = pvals_ssvd,
  #                                                'Standard'   = pvals_standard) ,
  #                                       M = M, 
  #                                       method_colors = c('dodgerblue3', 'palegreen3', 'orangered2'))
  
  p_pvals_all_qqplots = plot_pvals_many(P = list('LR: True'        = pvals_true,
                                                 'LR: SVD'         = pvals_svd,
                                                 'LR: Sparse SVD'  = pvals_ssvd,
                                                 'LR: Mat Comp [Linear Reg]'     = pvals_mcomp_linearreg,
                                                 'LR: Mat Comp [Soft Impute]'    = pvals_mcomp_softImpute,
                                                 'Standard'        = pvals_standard) ,
                                        M = M,  
                                        method_colors = c('mediumpurple3', 'dodgerblue2', 'dodgerblue4', 'palegreen2', 'palegreen4', 'orangered2'))
  
  
  # grob <- gridExtra::arrangeGrob(p_M, p_M1, p_M2,
  #                                p_Xsvd, 
  #                                p_pval_Xsvd, p_pval_standard,
  #                                p_pvals_svd_qqplots$p_null, p_pvals_svd_qqplots$p_alt,
  #                                p_pvals_standard_qqplots$p_null, p_pvals_standard_qqplots$p_alt,
  #                                layout_matrix = matrix(c(1, 2, 3,
  #                                                         4, 5, 6), byrow = TRUE, nrow = 2))
  # grob <- gridExtra::arrangeGrob(p_M, p_M1, p_M2, # 1,2,3
  #                                p_Xsvd,           # 4
  #                                p_pval_Xsvd, p_pval_standard, # 5 , 6
  #                                p_pvals_svd_qqplots$p_null, p_pvals_svd_qqplots$p_alt, # 7, 8
  #                                p_pvals_standard_qqplots$p_null, p_pvals_standard_qqplots$p_alt, # 9, 10
  #                                layout_matrix = matrix(c(NA, 1,  2, 3, 
  #                                                         4,  5, 7, 8,
  #                                                         NA, 6, 9, 10), byrow = TRUE, nrow = 3))
  # grob <- gridExtra::arrangeGrob(p_M, p_M1, p_M2, # 1,2,3
  #                                p_Xsvd, p_Xssvd,          # 4 5
  #                                p_pval_Xsvd, p_pval_Xssvd, p_pval_standard, # 6 7 8
  #                                p_pvals_svd_qqplots$p_null, p_pvals_svd_qqplots$p_alt, # 9 10
  #                                p_pvals_ssvd_qqplots$p_null, p_pvals_ssvd_qqplots$p_alt, # 11 12
  #                                p_pvals_standard_qqplots$p_null, p_pvals_standard_qqplots$p_alt, # 13, 14
  #                                layout_matrix = matrix(c(NA, 1,  2, 3, 
  #                                                         4,  6,  9, 10,
  #                                                         5,  7,  11, 12,
  #                                                         NA, 8,  13, 14), byrow = TRUE, nrow = 4))
  # grob <- gridExtra::arrangeGrob(p_M, p_Mobs,                # 1 , 2
  #                                p_Xsvd, p_Xssvd, p_Xmcomp,          # 3 4 5 
  #                                p_pval_Xsvd, p_pval_Xssvd, p_pval_Xmcomp, p_pval_standard, # 6 7 8 9 
  #                                p_pvals_all_qqplots$p_null, p_pvals_all_qqplots$p_alt,     # 10 11
  #                                # p_pvals_svd_qqplots$p_null, p_pvals_svd_qqplots$p_alt, # 
  #                                # p_pvals_ssvd_qqplots$p_null, p_pvals_ssvd_qqplots$p_alt, # 
  #                                # p_pvals_standard_qqplots$p_null, p_pvals_standard_qqplots$p_alt, # 
  #                                layout_matrix = matrix(c(NA, 1,  2, NA, 
  #                                                         3, 4, 5, NA,
  #                                                         6, 7, 8, 9,
  #                                                         10, 10, 11, 11,
  #                                                         10, 10, 11, 11), byrow = TRUE, nrow = 5))
  
  grob <- gridExtra::arrangeGrob(p_M, p_Mobs,                # 1 , 2
                                 p_Xsvd, p_Xssvd, p_Xmcomp_linearreg, p_Xmcomp_softImpute,         # 3 4 5 6
                                 p_pval_Xsvd, p_pval_Xssvd, p_pval_Xmcomp_linearreg, p_pval_Xmcomp_softImpute, p_pval_standard, # 7 8 9 10 11 
                                 p_pvals_all_qqplots$p_null, p_pvals_all_qqplots$p_alt,     #  12 13
                                 # p_pvals_svd_qqplots$p_null, p_pvals_svd_qqplots$p_alt, # 
                                 # p_pvals_ssvd_qqplots$p_null, p_pvals_ssvd_qqplots$p_alt, # 
                                 # p_pvals_standard_qqplots$p_null, p_pvals_standard_qqplots$p_alt, # 
                                 layout_matrix = matrix(c(NA, 1,   2, NA, NA,
                                                           3, 4,   5,  6, NA,
                                                           7, 8,   9, 10, 11,
                                                          12, 12, 13, 13, NA,
                                                          12, 12, 13, 13, NA), byrow = TRUE, nrow = 5))
   
  
  
  gridExtra::grid.arrange(grob)
  ggsave(sprintf('%s/all_nosamplesplit.pdf', save_folder), grob, width = 18, height = 12)
  return(grob)
  
}



#' @param M (matrix)
#' @param noise_sd (vector) of sd's for noise of M
#' @param svd_maxrank (integer) positive integer indicating rank for svd decompositions
#' @param save_folder (character) path to folder to save plots at '<save_folder>/all.pdf'
sim_and_plot_samplesplit <- function(M, noise_sd, svd_maxrank, save_folder) {
  # # debug
  # r = 5
  # n = 15
  # m = 10
  # set.seed(12345)
  # # M = create_lowrank_matrix(r=r, n=n, m=m)
  # M = create_blocky_matrix(r=r, n=n, m=m) * 3
  # 
  # noise_sd = runif(n*m, min=.1, max=.5)
  # svd_maxrank = 5
  
  dir.create(save_folder, recursive = FALSE, showWarnings = FALSE)
  n = nrow(M)
  m = ncol(M)
  
  
  # === Create Observed Matrix M_obs
  noise_sd_matrix = matrix(noise_sd, nrow = n, ncol = m, byrow = TRUE)
  
  M_split1 = M + matrix(rnorm(n = length(noise_sd), mean = 0, sd = noise_sd),
                        nrow = n, ncol = m, byrow = TRUE)
  M_split2 = M + matrix(rnorm(n = length(noise_sd), mean = 0, sd = noise_sd),
                        nrow = n, ncol = m, byrow = TRUE)
  
  
  
  
  # === Construct X matrices (svd, sparse svd, etc...) 
  
  # svd
  svdres = svd(M_split1, nu = svd_maxrank, nv = svd_maxrank)
  X_svd  = svdres$u %*% diag(svdres$d[1:svd_maxrank]) %*% t(svdres$v)
  
  
  # sparse svd
  cv.out = PMA::SPC.cv(M_split1, sumabsvs = seq(1.2, min(5, sqrt(n), sqrt(m)), len = 10))
  pmd_res = PMA::SPC(M_split1, sumabsv=cv.out$bestsumabsv, K = svd_maxrank)
  X_sparse_svd = pmd_res$u %*% diag(pmd_res$d) %*% t(pmd_res$v)
  
  #  single element matrix completion methods- Just use split2??
  X_mcomp_linearreg  = matrix(NA, nrow = n, ncol = m) 
  X_mcomp_softImpute = matrix(NA, nrow = n, ncol = m) 
  matcomp_fun = create_matrix_completion_softImpute(X = M_split2, rank_max = svd_maxrank, type = 'als')
  for(i in 1:n) {
    for(j in 1:m) {
      X_mcomp_linearreg[i,j]  = linearreg(M_split2, i, j)
      X_mcomp_softImpute[i,j] = matcomp_fun(i, j)
    }
  }
  
  # === Calculate p-values TODO: change to split2
  pvals_standard            = matrix(NA, nrow = n, ncol = m) # standard method
  pvals_true                = matrix(NA, nrow = n, ncol = m) # LR Method: use TRUE matrix (unrealistic setting, a sanity check)
  pvals_svd                 = matrix(NA, nrow = n, ncol = m) # LR Method: use matrix from svd
  pvals_ssvd                = matrix(NA, nrow = n, ncol = m) # LR Method: use matrix from sparse svd
  pvals_mcomp_linearreg     = matrix(NA, nrow = n, ncol = m) # LR Method: use matrix from single element matrix completion linearreg
  pvals_mcomp_softImpute    = matrix(NA, nrow = n, ncol = m) # LR Method: use matrix from single element matrix completion softImpute
  
  for(i in 1:n) {
    for(j in 1:m) {
      
      
      # Standard: AVG Y AND DECREASE SD! (use the sd from no sample splitting) sd / sqrt(2) (assuming even split into 2 parts)
      pvals_standard[i,j] = 2*(1 - pnorm(abs((M_split1[i,j] + M_split2[i,j])/2  - 0) /(noise_sd_matrix[i,j]/sqrt(2))  ))
      
      
      # LR: X is true M
      tstat_and_pval = calc_tstat_and_pval(
        Yij      = M_split2[i,j], 
        sij      = noise_sd_matrix[i,j], 
        thetaij  = M[i,j], 
        thetaij0 = 0
      )
      # pvals_true[i,j] = runif(n=1, min = .001, max = .999) # test code
      pvals_true[i,j] = tstat_and_pval$pval ; rm(tstat_and_pval)
      
      # LR: X is SVD
      tstat_and_pval = calc_tstat_and_pval(
        Yij      = M_split2[i,j], 
        sij      = noise_sd_matrix[i,j], 
        thetaij  = X_svd[i,j], 
        thetaij0 = 0
      )
      pvals_svd[i,j] = tstat_and_pval$pval ; rm(tstat_and_pval)
      
      # LR: X is Sparse SVD
      tstat_and_pval = calc_tstat_and_pval(
        Yij      = M_split2[i,j], 
        sij      = noise_sd_matrix[i,j], 
        thetaij  = X_sparse_svd[i,j], 
        thetaij0 = 0
      )
      pvals_ssvd[i,j] = tstat_and_pval$pval ; rm(tstat_and_pval)
      
      # LR: X is linearreg
      tstat_and_pval = calc_tstat_and_pval(
        Yij      = M_split2[i,j],
        sij      = noise_sd_matrix[i,j],
        thetaij  = X_mcomp_linearreg[i,j],
        thetaij0 = 0
      )
      pvals_mcomp_linearreg[i,j]  = tstat_and_pval$pval ; rm(tstat_and_pval)
      
      # LR: X is softImpute
      tstat_and_pval = calc_tstat_and_pval(
        Yij      = M_split2[i,j],
        sij      = noise_sd_matrix[i,j],
        thetaij  = X_mcomp_softImpute[i,j],
        thetaij0 = 0
      )
      pvals_mcomp_softImpute[i,j]  = tstat_and_pval$pval ; rm(tstat_and_pval)
    }
  }
  
  
  
  p_M    = display_matrix(M) + labs(title = 'M')
  p_Msplit1 = display_matrix(M_split1) + labs(title = 'M Split1')
  p_Msplit2 = display_matrix(M_split2) + labs(title = 'M Split2')
  
  
  p_Xsvd  = display_matrix(X_svd) + labs(title = TeX(r'($\tilde{\Theta}$ (SVD))'))
  p_Xssvd = display_matrix(X_sparse_svd) + labs(title = TeX(r'($\tilde{\Theta}$ (Sparse SVD))'))
  p_Xmcomp_linearreg  = display_matrix(X_mcomp_linearreg ) + labs(title = TeX(r'($\tilde{\Theta}$ (Mat. Comp.- Linear Reg))'))
  p_Xmcomp_softImpute = display_matrix(X_mcomp_softImpute) + labs(title = TeX(r'($\tilde{\Theta}$ (Mat. Comp.- Soft Impute))'))
  # display_matrix(X)
  # display_matrix(tstats)
  # display_matrix(pvals) + scale_fill_gradientn(colours = c('red', 'orange', 'white'), values = c(0, .05, 1))
  
  p_pval_standard = display_matrix(pvals_standard) + scale_fill_gradientn(colours = c('red', 'orange', 'white'), values = c(0, .05, 1)) + labs(title = 'Standard p-vals')
  p_pval_Xsvd     = display_matrix(pvals_svd)  + scale_fill_gradientn(colours = c('red', 'orange', 'white'), values = c(0, .05, 1)) + labs(title = TeX(r'($\tilde{\Theta}$ (SVD) p-vals)'))
  p_pval_Xssvd    = display_matrix(pvals_ssvd) + scale_fill_gradientn(colours = c('red', 'orange', 'white'), values = c(0, .05, 1)) + labs(title = TeX(r'($\tilde{\Theta}$ (Sparse SVD) p-vals)'))
  p_pval_Xmcomp_linearreg  = display_matrix(pvals_mcomp_linearreg ) + scale_fill_gradientn(colours = c('red', 'orange', 'white'), values = c(0, .05, 1)) + labs(title = TeX(r'($\tilde{\Theta}$ (Mat. Comp.- Linear Reg) p-vals)'))
  p_pval_Xmcomp_softImpute = display_matrix(pvals_mcomp_softImpute) + scale_fill_gradientn(colours = c('red', 'orange', 'white'), values = c(0, .05, 1)) + labs(title = TeX(r'($\tilde{\Theta}$ (Mat. Comp.- Soft Impute) p-vals)'))
  
  
  p_pvals_svd_qqplots       = plot_pvals(P = pvals_svd      , M = M)
  p_pvals_ssvd_qqplots      = plot_pvals(P = pvals_ssvd     , M = M)
  p_pvals_standard_qqplots  = plot_pvals(P = pvals_standard , M = M)
  
  
  # p_pvals_all_qqplots = plot_pvals_many(P = list('SVD'        = pvals_svd,
  #                                                'Sparse SVD' = pvals_ssvd,
  #                                                'Standard'   = pvals_standard) ,
  #                                       M = M, 
  #                                       method_colors = c('dodgerblue3', 'palegreen3', 'orangered2'))
  
  p_pvals_all_qqplots = plot_pvals_many(P = list('LR: True'        = pvals_true,
                                                 'LR: SVD'         = pvals_svd,
                                                 'LR: Sparse SVD'  = pvals_ssvd,
                                                 'LR: Mat Comp [Linear Reg]'     = pvals_mcomp_linearreg,
                                                 'LR: Mat Comp [Soft Impute]'    = pvals_mcomp_softImpute,
                                                 'Standard (No Samp Split)'        = pvals_standard) ,
                                        M = M,  
                                        method_colors = c('mediumpurple3', 'dodgerblue2', 'dodgerblue4', 'palegreen2', 'palegreen4', 'orangered2'))
  
  
  
  grob <- gridExtra::arrangeGrob(p_M, p_Msplit1, p_Msplit2,               # 1 , 2 3
                                 p_Xsvd, p_Xssvd, p_Xmcomp_linearreg, p_Xmcomp_softImpute,         #  4 5 6 7
                                 p_pval_Xsvd, p_pval_Xssvd, p_pval_Xmcomp_linearreg, p_pval_Xmcomp_softImpute, p_pval_standard, #  8 9 10 11 12
                                 p_pvals_all_qqplots$p_null, p_pvals_all_qqplots$p_alt,     #  12 13
                                 # p_pvals_svd_qqplots$p_null, p_pvals_svd_qqplots$p_alt, # 
                                 # p_pvals_ssvd_qqplots$p_null, p_pvals_ssvd_qqplots$p_alt, # 
                                 # p_pvals_standard_qqplots$p_null, p_pvals_standard_qqplots$p_alt, # 
                                 layout_matrix = matrix(c(NA,  1,  2,  3, NA,
                                                           4,  5,  6,  7, NA,
                                                           8,  9, 10, 11, 12,
                                                          13, 13, 14, 14, NA,
                                                          13, 13, 14, 14, NA), byrow = TRUE, nrow = 5))
  gridExtra::grid.arrange(grob)
  ggsave(sprintf('%s/all_samplesplit.pdf', save_folder), grob, width = 15, height = 12)
  return(grob)
}




#' @param M (matrix)
#' @param noise_sd (vector) of sd's for noise of M
#' @param svd_maxrank (integer) positive integer indicating rank for svd decompositions
#' @param save_folder (character) path to folder to save plots at '<save_folder>/all.pdf'
sim_and_plot_old <- function(M, noise_sd, svd_maxrank, save_folder) {
  # # debug
  # r = 5
  # n = 15
  # m = 10
  # set.seed(12345)
  # # M = create_lowrank_matrix(r=r, n=n, m=m)
  # M = create_blocky_matrix(r=r, n=n, m=m) * 3
  # 
  # noise_sd = runif(n*m, min=.1, max=.5)
  # svd_maxrank = 5
  
  dir.create(save_folder, recursive = FALSE, showWarnings = FALSE)
  n = nrow(M)
  m = ncol(M)
  
  noise_sd_matrix = matrix(noise_sd, nrow = n, ncol = m, byrow = TRUE)
  
  M_split1 = M + matrix(rnorm(n = length(noise_sd), mean = 0, sd = noise_sd),
                        nrow = n, ncol = m, byrow = TRUE)
  M_split2 = M + matrix(rnorm(n = length(noise_sd), mean = 0, sd = noise_sd),
                        nrow = n, ncol = m, byrow = TRUE)
  
  
  svdres = svd(M_split1, nu = svd_maxrank, nv = svd_maxrank)
  X_svd  = svdres$u %*% diag(svdres$d[1:svd_maxrank]) %*% t(svdres$v)
  
  # should use cv to tune hyperparams
  # cv.out = PMA::PMD.cv(M_split1, type="standard", sumabss=seq(0.1, 0.6, len=20))
  # # pmd_res = PMA::PMD(x=M_split1, type = 'standard', K = svd_maxrank)
  # pmd_res = PMA::PMD(x=M_split1, type = 'standard', K = svd_maxrank, sumabs=cv.out$bestsumabs, v=cv.out$v.init)
  # X_sparse_svd = pmd_res$u %*% diag(pmd_res$d) %*% t(pmd_res$v)
  
  cv.out = PMA::SPC.cv(M_split1, sumabsvs = seq(1.2, min(5, sqrt(n), sqrt(m)), len = 10))
  pmd_res = PMA::SPC(M_split1, sumabsv=cv.out$bestsumabsv, K = svd_maxrank)
  X_sparse_svd = pmd_res$u %*% diag(pmd_res$d) %*% t(pmd_res$v)
  
  
  # X      = matrix(NA, nrow = n, ncol = m) 
  # tstats = matrix(NA, nrow = n, ncol = m)
  # pvals  = matrix(NA, nrow = n, ncol = m)
  pvals_standard = matrix(NA, nrow = n, ncol = m)
  pvals_svd      = matrix(NA, nrow = n, ncol = m)
  pvals_ssvd     = matrix(NA, nrow = n, ncol = m)   
  for(i in 1:n) {
    for(j in 1:m) {
      
      # X is custom fill fn
      # X[i,j] = average(M_split1, i, j)
      # X[i,j] = linearreg(M_split1, i, j)
      # tstat_and_pval = calc_tstat_and_pval(
      #   Yij      = M_split2[i,j], 
      #   sij      = noise_sd_matrix[i,j], 
      #   thetaij  = X[i,j], 
      #   thetaij0 = 0
      # )
      # tstats[i,j] = tstat_and_pval$tstat
      # pvals[i,j]  = tstat_and_pval$pval
      
      
      # X is SVD
      tstat_and_pval = calc_tstat_and_pval(
        Yij      = M_split2[i,j], 
        sij      = noise_sd_matrix[i,j], 
        thetaij  = X_svd[i,j], 
        thetaij0 = 0
      )
      pvals_svd[i,j] = tstat_and_pval$pval
      
      # X is Sparse SVD
      tstat_and_pval = calc_tstat_and_pval(
        Yij      = M_split2[i,j], 
        sij      = noise_sd_matrix[i,j], 
        thetaij  = X_sparse_svd[i,j], 
        thetaij0 = 0
      )
      pvals_ssvd[i,j] = tstat_and_pval$pval
      
      # Standard: AVG Y AND DECREASE SD! (use the sd from no sample splitting) sd / sqrt(2) (assuming even split into 2 parts)
      pvals_standard[i,j] = 2*(1 - pnorm(abs((M_split1[i,j] + M_split2[i,j])/2  - 0) /(noise_sd_matrix[i,j]/sqrt(2))  ))
    }
  }
  
  
  
  
  
  p_M = display_matrix(M) + labs(title = 'M')
  p_M1 = display_matrix(M_split1) + labs(title = 'M Split 1')
  p_M2 = display_matrix(M_split2) + labs(title = 'M Split 2')
  
  
  p_Xsvd  = display_matrix(X_svd) + labs(title = TeX(r'($\tilde{\Theta}$ (SVD))'))
  p_Xssvd = display_matrix(X_sparse_svd) + labs(title = TeX(r'($\tilde{\Theta}$ (Sparse SVD))'))
  # display_matrix(X)
  # display_matrix(tstats)
  # display_matrix(pvals) + scale_fill_gradientn(colours = c('red', 'orange', 'white'), values = c(0, .05, 1))
  
  p_pval_Xsvd  = display_matrix(pvals_svd) + scale_fill_gradientn(colours = c('red', 'orange', 'white'), values = c(0, .05, 1)) + labs(title = TeX(r'($\tilde{\Theta}$ (SVD) p-vals)'))
  p_pval_Xssvd = display_matrix(pvals_ssvd) + scale_fill_gradientn(colours = c('red', 'orange', 'white'), values = c(0, .05, 1)) + labs(title = TeX(r'($\tilde{\Theta}$ (Sparse SVD) p-vals)'))
  p_pval_standard = display_matrix(pvals_standard) + scale_fill_gradientn(colours = c('red', 'orange', 'white'), values = c(0, .05, 1)) + labs(title = 'Standard p-vals')
  
  
  
  p_pvals_svd_qqplots       = plot_pvals(P = pvals_svd      , M = M)
  p_pvals_ssvd_qqplots      = plot_pvals(P = pvals_ssvd     , M = M)
  p_pvals_standard_qqplots  = plot_pvals(P = pvals_standard , M = M)
  
  
  p_pvals_all_qqplots = plot_pvals_many(P = list('SVD'        = pvals_svd,
                                                 'Sparse SVD' = pvals_ssvd,
                                                 'Standard'   = pvals_standard) ,
                                        M = M, 
                                        method_colors = c('dodgerblue3', 'palegreen3', 'orangered2'))
  
  
  
  
  # grob <- gridExtra::arrangeGrob(p_M, p_M1, p_M2,
  #                                p_Xsvd, 
  #                                p_pval_Xsvd, p_pval_standard,
  #                                p_pvals_svd_qqplots$p_null, p_pvals_svd_qqplots$p_alt,
  #                                p_pvals_standard_qqplots$p_null, p_pvals_standard_qqplots$p_alt,
  #                                layout_matrix = matrix(c(1, 2, 3,
  #                                                         4, 5, 6), byrow = TRUE, nrow = 2))
  # grob <- gridExtra::arrangeGrob(p_M, p_M1, p_M2, # 1,2,3
  #                                p_Xsvd,           # 4
  #                                p_pval_Xsvd, p_pval_standard, # 5 , 6
  #                                p_pvals_svd_qqplots$p_null, p_pvals_svd_qqplots$p_alt, # 7, 8
  #                                p_pvals_standard_qqplots$p_null, p_pvals_standard_qqplots$p_alt, # 9, 10
  #                                layout_matrix = matrix(c(NA, 1,  2, 3, 
  #                                                         4,  5, 7, 8,
  #                                                         NA, 6, 9, 10), byrow = TRUE, nrow = 3))
  # grob <- gridExtra::arrangeGrob(p_M, p_M1, p_M2, # 1,2,3
  #                                p_Xsvd, p_Xssvd,          # 4 5
  #                                p_pval_Xsvd, p_pval_Xssvd, p_pval_standard, # 6 7 8
  #                                p_pvals_svd_qqplots$p_null, p_pvals_svd_qqplots$p_alt, # 9 10
  #                                p_pvals_ssvd_qqplots$p_null, p_pvals_ssvd_qqplots$p_alt, # 11 12
  #                                p_pvals_standard_qqplots$p_null, p_pvals_standard_qqplots$p_alt, # 13, 14
  #                                layout_matrix = matrix(c(NA, 1,  2, 3, 
  #                                                         4,  6,  9, 10,
  #                                                         5,  7,  11, 12,
  #                                                         NA, 8,  13, 14), byrow = TRUE, nrow = 4))
  grob <- gridExtra::arrangeGrob(p_M, p_M1,  p_M2, # 1,2,3
                                 p_Xsvd, p_Xssvd,          # 4 5
                                 p_pval_Xsvd, p_pval_Xssvd, p_pval_standard, # 6 7 8
                                 p_pvals_svd_qqplots$p_null, p_pvals_svd_qqplots$p_alt, # 9 10
                                 p_pvals_ssvd_qqplots$p_null, p_pvals_ssvd_qqplots$p_alt, # 11 12
                                 p_pvals_standard_qqplots$p_null, p_pvals_standard_qqplots$p_alt, # 13, 14
                                 p_pvals_all_qqplots$p_null, p_pvals_all_qqplots$p_alt, # 15, 16
                                 layout_matrix = matrix(c(NA, 1,  2, 3, 
                                                          4,  6,  9, 10,
                                                          5,  7,  11, 12,
                                                          NA, 8,  13, 14,
                                                          NA, NA, 15, 16), byrow = TRUE, nrow = 5))
  gridExtra::grid.arrange(grob)
  ggsave(sprintf('%s/all.pdf', save_folder), grob, width = 15, height = 12)
  return(grob)
}



