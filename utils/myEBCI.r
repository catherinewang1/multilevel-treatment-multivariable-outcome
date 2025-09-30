# Functions to perform EBCI shrinkage 




#' estimate mu2 by hand
#' we dont have any X's
#' @param Y (vector) of numeric values, the Y_i's 
#'     in our case, this is unshrunk value - shrinkage point
#'     bc formula = unshrunk val - shrinkage point - 0
#'     which may be incorrect!
#' @param se (vector) of numeric values, the standard errors_i's
#'     this is \hat{\sigma}_i  in the formulas
#' @param weights (vector) of numeric values, the weights to put on each sample
#' @param X (matrix) of numeric values, the covariates
#' @param delta (vector) of coefs for covariates  
#' @returns estimated mu2 (numeric value)
#' @example
#' # this actually gives the same... 0.0007741335
#' my_est_mu2(Y = unlist(matrices_test$estimates[c(1,2), ] - 
#'                         approxmatrices_sparseSVD$approxmatrices[[10]][c(1,2), ]),
#'            se = as.vector(se_mat_test[c(1,2), ]),
#'            weights = rep(1/length(unlist(se_mat_test[c(1,2), ])), 
#'                          times = length(unlist(se_mat_test[c(1,2), ]))))
#' 
my_est_mu2 <- function(Y, se, weights, X=NULL, delta=NULL, print=FALSE) {
  # Y = matrices_test$estimates[1:50, 1] - approxmatrices_sparseSVD$approxmatrices[[10]][1:50, 1]
  # se = se_mat_test[1:50, 1]
  # weights = rep(1/length(se), times = length(se))
  # X=NULL; delta=NULL
  
  if(is.null(X) | is.null(delta)) {
    eps = Y
  } else {
    eps = Y - X %*% delta
  }
  
  
  first = (sum(weights * (eps*eps - se*se))
           /
             sum(weights))
  
  
  second = 2 * (sum(weights * weights * (se^4))
                /
                  (sum(weights) * sum(weights * se * se)))
  
  if(print) {print(sprintf('first = %.6f, second = %.6f', first, second))}
  
  return(max(first, second))
}


#' Estimate kappa by hand
#' @param Y (vector) of numeric values, the Y_i's 
#'     in our case, this is unshrunk value - shrinkage point
#'     bc formula = unshrunk val - shrinkage point - 0
#'     which may be incorrect!
#' @param se (vector) of numeric values, the standard errors_i's
#'     this is \hat{\sigma}_i  in the formulas
#' @param weights (vector) of numeric values, the weights to put on each sample
#' @param mu2 (numeric) estimate of mu2 (from a prev step)
#' @param X (matrix) of numeric values, the covariates
#' @param delta (vector) of coefs for covariates  
#' @example
#' my_est_kappa(Y = matrices_test$estimates[1:50, 1] - approxmatrices_sparseSVD$approxmatrices[[10]][1:50, 1],
#'              se = se_mat_test[1:50, 1],
#'              mu2 = shrink2Grnas$ebci_obj$mu2['estimate'],
#'              weights = rep(1/length(se), times = length(se)))
#' 
#' # this actually gives the same...2848206
#' my_est_kappa(Y = unlist(matrices_test$estimates[c(1,2), ] - 
#'                           approxmatrices_sparseSVD$approxmatrices[[10]][c(1,2), ]),
#'              se = as.vector(se_mat_test[c(1,2), ]),
#'              # mu2 = 5,
#'              mu2 = shrink2Grnas$ebci_obj$mu2['estimate'],
#'              # mu2 = shrink2Grnas$ebci_obj$mu2['uncorrected_estimate'],
#'              weights = rep(1/length(unlist(se_mat_test[c(1,2), ])), 
#'                            times = length(unlist(se_mat_test[c(1,2), ]))))
my_est_kappa <- function(Y, se, weights, mu2, X=NULL, delta=NULL, print=FALSE) {
  # Y = matrices_test$estimates[1:50, 1] - approxmatrices_sparseSVD$approxmatrices[[10]][1:50, 1]
  # se = se_mat_test[1:50, 1]
  # weights = rep(1/length(se), times = length(se))
  # X=NULL; delta=NULL
  # mu2 = .000774
  
  if(is.null(X) | is.null(delta)) {
    eps = Y
  } else {
    eps = Y - X %*% delta
  }
  
  first = (sum( weights * (eps^4 - 6 * se * se * eps * eps + 3 * se^4))
           /
             (mu2*mu2 * sum(weights))
  )
  
  
  second = 1 + 
    (32 * sum(weights * weights * se^8)
     /
       (mu2*mu2 * sum(weights) * sum(weights * se^4))
    )
  
  
  if(print) {print(sprintf('first = %.6f, second = %.6f', first, second))}
  return(max(first, second))
}







#' Calc the eb shrunk value
#' 
#' @param Y (vector) of numeric values, the Y_i's 
#'     in our case, this is unshrunk value - shrinkage point
#'     bc formula = unshrunk val - shrinkage point - 0
#'     which may be incorrect!
#' @param se (vector) of numeric values, the standard errors_i's
#'     this is \hat{\sigma}_i  in the formulas
#' @param weights (vector) of numeric values, the weights to put on each sample
#' @param mu2 (numeric) estimate of mu2 (from a prev step)
#' @param kappa (numeric) estimate of kappa (from a prev step)
#' @param X (matrix) of numeric values, the covariates
#' @param delta (vector) of coefs for covariates 
#' 
#' @example
#' theb = my_est_theb(
#'        Y = unlist(matrices_test$estimates[c(1,2), ] - 
#'        approxmatrices_sparseSVD$approxmatrices[[10]][c(1,2), ]),
#'        se = as.vector(se_mat_test[c(1,2), ]),
#'        # mu2 = 5,
#'        mu2 = shrink2Grnas$ebci_obj$mu2['estimate'],
#'        # mu2 = shrink2Grnas$ebci_obj$mu2['uncorrected_estimate'],
#'        kappa = shrink2Grnas$ebci_obj$kappa['estimate'],
#'        weights = rep(1/length(unlist(se_mat_test[c(1,2), ])), 
#'                      times = length(unlist(se_mat_test[c(1,2), ]))))
#'        theb$w_eb |> hist()
#'        theb$th_eb |> hist()
#'                
my_est_theb <- function(Y, se, weights, mu2, X=NULL, delta=NULL) {
  # Y = matrices_test$estimates[1:50, 1] - approxmatrices_sparseSVD$approxmatrices[[10]][1:50, 1]
  # se = se_mat_test[1:50, 1]
  # weights = rep(1/length(se), times = length(se))
  # X=NULL; delta=NULL
  # mu2 = .000774
  
  if(is.null(X) | is.null(delta)) {
    c = rep(0, times = length(Y))
  } else {
    c = X %*% delta
  }
  
  w_eb = mu2 / (mu2 + se*se) # vector
  
  th_eb = c + w_eb * (Y - c)
  
  
  # return(list(w_eb = w_eb, th_eb = th_eb))
  return(data.frame(w_eb = w_eb, th_eb = th_eb))
}





#' Calc the eb shrunk value's confidence interval
#' This step takes a while from the ebci::cva call 
#' 
#' @param th_eb (vector) of eb shrunk estimate (from a prev step)
#' @param w_eb (vector) of eb shrinkage factors (from a prev step)
#' @param se (vector) of numeric values, the standard errors_i's
#'     this is \hat{\sigma}_i  in the formulas
#' @param weights (vector) of numeric values, the weights to put on each sample
#' @param mu2 (numeric) estimate of mu2 (from a prev step)
#' @param kappa (numeric) estimate of kappa (from a prev step)
#'
#' @example
#' th_ci = my_est_ci(th_eb = theb$th_eb,
#'                   w_eb  = theb$w_eb,
#'                   se = as.vector(se_mat_test[c(1,2), ]),
#'                   mu2 = shrink2Grnas$ebci_obj$mu2['estimate'],
#'                   kappa = shrink2Grnas$ebci_obj$kappa['estimate'],
#'                   alpha = .1)
#' 
my_est_ci <- function(th_eb, w_eb, se, mu2, kappa, alpha)  {
  cv = rep(NA, times = length(se)) 
  m2_ = se*se / mu2
  for(i in 1:length(cv)) {
    # if(i %% 100 == 0){print(i)}
    
    # !!! This is kind of slow!!! (this is where the time will be spent)
    cv[i] = ebci::cva(m2 = m2_[i], kappa=kappa, alpha=alpha)$cv
  }
  
  len_eb = cv * w_eb * se
  
  return(data.frame(m2 = m2_, cv=cv,
                    len_eb = len_eb,
                    ci_lower = th_eb - len_eb,
                    ci_upper = th_eb + len_eb))
}


#' calc ebci by hand to compare
#'
#' @param Y (vector) of numeric values, the Y_i's 
#'     in our case, this is unshrunk value - shrinkage point
#'     bc formula = unshrunk val - shrinkage point - 0
#'     which may be incorrect!
#' @param se (vector) of numeric values, the standard errors_i's
#'     this is \hat{\sigma}_i  in the formulas
#' @param weights (vector) of numeric values, the weights to put on each sample
#' @param mu2 (numeric) estimate of mu2 (from a prev step)
#' @param kappa (numeric) estimate of kappa (from a prev step)
#' @param X (matrix) of numeric values, the covariates
#' @param delta (vector) of coefs for covariates 
#' @param mu2 (numeric) overwrite estimation of mu2 with this value
#'  (e.g. estimate mu2 in a different way, but do eb shrinkage w/ this estimate)
#' @param kappa (numeric) overwrite estimation of kappa with this value
#'  (e.g. estimate kappa in a different way, but do eb shrinkage w/ this estimate)
#' @returns ebci results as a list:
#' $mu2: estimate of mu2
#' $kappa: estimate of kappa
#' $df: data.frame with eb value columns
#'     $th_eb: theta estimate
#'     $w_eb: shrinkage factor
#'     $ci_lower: eb CI lower value
#'     $ci_upper: eb CI upper value
#' 
#' @example
#' my_ebci_obj = 
#'   my_ebci(Y = unlist(matrices_test$estimates[c(1,2), ] - 
#'                        approxmatrices_sparseSVD$approxmatrices[[10]][c(1,2), ]),
#'           se = as.vector(se_mat_test[c(1,2), ]),
#'           weights = rep(1/length(as.vector(se_mat_test[c(1,2), ])), 
#'                         times = length(as.vector(se_mat_test[c(1,2), ]))),
#'           alpha = .1)
#' 
#' 
#' my_ebci_obj$mu2
#' my_ebci_obj$kappa
#' 
#' head(my_ebci_obj$df) 
#' 
my_ebci <- function(Y, se, weights, alpha, X=NULL, delta=NULL, 
                    mu2=NULL, kappa=NULL) {
  
  # estimate mu2 = E[(\theta - X^T delta)^2 |X]  and kappa ...
  # or use input if specified
  if(is.null(mu2)) {
    mu2   = my_est_mu2(  Y=Y, se=se, weights=weights,          X=X, delta=delta)
  } 
  if(is.null(kappa)) {
    kappa = my_est_kappa(Y=Y, se=se, weights=weights, mu2=mu2, X=X, delta=delta)
  }
  
  # estimate shrunk values
  th_df = my_est_theb( Y=Y, se=se, weights=weights, mu2=mu2, X=X, delta=delta)
  ci_df = my_est_ci(th_eb=th_df$th_eb, w_eb=th_df$w_eb, se=se, mu2=mu2, kappa=kappa, alpha=alpha)
  
  
  
  
  return(list(mu2 = mu2,
              kappa = kappa,
              df =  cbind(th_df, ci_df,
                          data.frame(Y=Y, se=se, weights=weights))))
  
}













# example of how to call 
if(F){
  
  # call fns individually
  # uses obj in testMatrixShrinkage
  #  - matrices_test
  #  - approxmatrices_sparseSVD
  #  - se_mat_test
  #  
  #  - shrink2Grnas <- compare to this result from ebci
  
  
  my_est_mu2(Y = matrices_test$estimates[1:50, 1] - approxmatrices_sparseSVD$approxmatrices[[10]][1:50, 1],
             se = se_mat_test[1:50, 1],
             weights = rep(1/length(se_mat_test[1:50, 1]), times = length(se_mat_test[1:50, 1])))
  
  # this actually gives the same... 0.0007741335
  my_est_mu2(Y = unlist(matrices_test$estimates[c(1,2), ] - 
                          approxmatrices_sparseSVD$approxmatrices[[10]][c(1,2), ]),
             se = as.vector(se_mat_test[c(1,2), ]),
             weights = rep(1/length(unlist(se_mat_test[c(1,2), ])), 
                           times = length(unlist(se_mat_test[c(1,2), ]))))
  
  
  # as sample size increases, the estimate decreases???
  for(j in seq(from = 50, to= 1500, by = 50)) {
    my_est_mu2(Y = unlist(matrices_test$estimates[c(1,2), 1:j] - 
                            approxmatrices_sparseSVD$approxmatrices[[10]][c(1,2), 1:j]),
               se = as.vector(se_mat_test[c(1,2), 1:j]),
               weights = rep(1/length(as.vector(se_mat_test[c(1,2), 1:j])), 
                             times = length(as.vector(se_mat_test[c(1,2), 1:j]))))
  }
  
  
  
  my_est_kappa(Y = matrices_test$estimates[1:50, 1] - approxmatrices_sparseSVD$approxmatrices[[10]][1:50, 1],
               se = se_mat_test[1:50, 1],
               mu2 = shrink2Grnas$ebci_obj$mu2['estimate'],
               weights = rep(1/length(se), times = length(se)))
  
  # this actually gives the same...2848206
  my_est_kappa(Y = unlist(matrices_test$estimates[c(1,2), ] - 
                            approxmatrices_sparseSVD$approxmatrices[[10]][c(1,2), ]),
               se = as.vector(se_mat_test[c(1,2), ]),
               # mu2 = 5,
               mu2 = shrink2Grnas$ebci_obj$mu2['estimate'],
               # mu2 = shrink2Grnas$ebci_obj$mu2['uncorrected_estimate'],
               weights = rep(1/length(unlist(se_mat_test[c(1,2), ])), 
                             times = length(unlist(se_mat_test[c(1,2), ]))))
  
  theb = my_est_theb(Y = unlist(matrices_test$estimates[c(1,2), ] - 
                                  approxmatrices_sparseSVD$approxmatrices[[10]][c(1,2), ]),
                     se = as.vector(se_mat_test[c(1,2), ]),
                     # mu2 = 5,
                     mu2 = shrink2Grnas$ebci_obj$mu2['estimate'],
                     # mu2 = shrink2Grnas$ebci_obj$mu2['uncorrected_estimate'],
                     kappa = shrink2Grnas$ebci_obj$kappa['estimate'],
                     weights = rep(1/length(unlist(se_mat_test[c(1,2), ])), 
                                   times = length(unlist(se_mat_test[c(1,2), ]))))
  theb$w_eb |> hist()
  theb$th_eb |> hist()
  
  
  th_ci = my_est_ci(th_eb = theb$th_eb,
                    w_eb  = theb$w_eb,
                    se = as.vector(se_mat_test[c(1,2), ]),
                    mu2 = shrink2Grnas$ebci_obj$mu2['estimate'],
                    kappa = shrink2Grnas$ebci_obj$kappa['estimate'],
                    alpha = .1)
  
  
  
  my_ebci_obj = 
    my_ebci(Y = unlist(matrices_test$estimates[c(1,2), ] - 
                         approxmatrices_sparseSVD$approxmatrices[[10]][c(1,2), ]),
            se = as.vector(se_mat_test[c(1,2), ]),
            weights = rep(1/length(as.vector(se_mat_test[c(1,2), ])), 
                          times = length(as.vector(se_mat_test[c(1,2), ]))),
            alpha = .1)
  
  
  my_ebci_obj$mu2
  my_ebci_obj$kappa
  
  head(my_ebci_obj$df) 
  
  
  
}