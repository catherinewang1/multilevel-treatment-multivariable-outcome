
#' description: Create the pseudo p-value: Find the smallest alpha that does not cover 0 OR the largest alpha that still covers 0
#'  The CI is shrunk estimate +- cva * w_eb * sigma
#'  because of the way we call the function, we have our original estimate thetahat and shrinkagepoint
#'  and we shrink (thetahat - shrinkagepoint) towards 0 always to then get thetaeb(this is because we have diff shrinkage points for all the points)
#'  and then we add back the shrinkage point shrunkpoint = thetaeb + shrinkagepoint <-- this is the value to to input into this function, because the CI is centered around this point
#' @param thetashrunk (numeric) shrunkpoint = thetaeb + shrinkagepoint where is the thetaeb is the (thetahat - shrinkagepoint) shrunk towards 0
#' @param sigma (numeric)
#' @param web (numeric) shrinkage factor from ebci object result
#' @param mu2 (numeric) estimated mu2 from ebci object result
#' @param kappa (numeric) estimated kappa from ebci object result
#' @param alpha_threshold (numeric) stop condition: perform until alpha's changes are < threshold
#'                            this is the max mis calc error for the pseudo-pval (if MAX_ITER is not reached)
#' @param maxiter (numeric) stop condition: maximum number of iterations to perform
#' 
#' @returns pseudopval (numeric) \in [0,1]
#' 
#' 
get_ebci_pvals <- function(thetashrunk, sigma, web, mu2, kappa,
                           alpha_threshold, maxiter) {
  # # params for running 
  # MAX_ITER_FOR_P = 20    # limit the number of iterations
  # # distance_threshold  = .001 # distance between (my_theta / (my_web * my_sigma)) and cva_pseudop(m2, kappa)
  # alpha_threshold = .0001 # perform until changes in alpha are small (ie this is the rounding of the p-value, there will be at most log2(1/alpha_threshold) iterations)
  # # estimate from overall ebci fit
  # my_mu2 = ebci_obj$mu2[['estimate']]
  # my_kappa = ebci_obj$kappa[['estimate']]
  # # params for each sample i 
  # my_theta = .2 # shrunk estimate: should be ebci_obj$df$th_eb + shrinkage_point (e.g. not the 'raw' the_eb)
  # my_sigma = 1 # initial estimate standard error
  # my_web = .3   # shrinkage factor

  # helpful calcs to just perform once
  sigma2_mu2 = sigma^2 / mu2  # = m2 in cva function input
  theta_websigma = abs(thetashrunk) / (web * sigma) # |theta| / (w_eb * sigma)   normalized estimate, always > 0



  iter = 1
  cur_min = 0; cur_max = 1 # range to search for pseudo p-value
  cur_alpha = .5 # start in the middle
  track_alphas = c(cur_alpha)
  # distance = 10000
  alpha_change = 1
  # while(distance > distance_threshold & iter <= MAX_ITER_FOR_P) { # by theta_websigma - cur_cva_alpha
  while(alpha_change > alpha_threshold & iter <= maxiter) { # by alpha/p-value change
   cur_cva_alpha = ebci::cva(m2 = sigma2_mu2, kappa = kappa, check=FALSE, alpha = cur_alpha)$cv
   # print(sprintf("%d: %.8f [%.2f, %.2f]: %.2f vs %.2f", 
   #                iter, cur_alpha, cur_min, cur_max, theta_websigma, cur_cva_alpha))
   
   # check if 0 is in CI with cur_alpha
   if(theta_websigma - cur_cva_alpha < 0) { # 0 \in CI: increase alpha to make CI smaller
     cur_min = cur_alpha  
     new_alpha = (cur_alpha + cur_max) / 2
   } else { # 0 \not\in CI: decrease alpha to make CI larger
     cur_max = cur_alpha
     new_alpha = (cur_alpha + cur_min) / 2
   }
   
   
   alpha_change = abs(new_alpha - cur_alpha)
   cur_alpha = new_alpha
   
   track_alphas = c(track_alphas, cur_alpha)
   
   # distance = theta_websigma - cur_cva_alpha
   iter = iter + 1
  }

  return(cur_alpha)
}



#' Use previously saved shrinkage results (object and dataframe) to 
#' get the ebci inverted pvalues. 
#' params named kind of badly...
#' @param cur_ebci_params (list)
#' @param cur_shrinkage_results (dataframe) 
#' @param parallel (boolean) run in parallel or not, default=FALSE
#'        requires parallel setup using future.apply package:
#'            suppressPackageStartupMessages(library(future.apply))
#'            plan(multisession, workers = 20)  # or some other plan
#' @param alpha_threshold (numeric) stop condition: perform until alpha's changes are < threshold
#'                            this is the max mis calc error for the pseudo-pval (if MAX_ITER is not reached)
#'        (pass into get_ebci_pvals call)
#' @param maxiter (numeric) stop condition: maximum number of iterations to perform
#'        (pass into get_ebci_pvals call)
#' 
get_ebci_pvals_by_shrinkage_results <- function(cur_ebci_params, cur_shrinkage_results, alpha_threshold=.001, maxiter=10, parallel=FALSE) {
  # sim_distn = 'nb'
  # split_type = 'samplesplit'
  # approx_method = 'zeros'
  # 
  # maxiter = 10
  # alpha_threshold = .001
  # 
  # cur_ebci_params = ebci_params[[sim_distn]][[split_type]][[approx_method]]
  # cur_ebci_params
  # 
  # 
  # cur_shrinkage_results = shrinkage_results[[sim_distn]][[split_type]][[approx_method]]
  # cur_shrinkage_results |> arrange(gene, grna) # dataframe for each #pert * #gene
  
  inner_func <- function(thetashrunk, sigma, web) { # This created function will take in (thetashrunk, sigma, web)
    # call using the previously specified stopping conditions alpha_threshold and maxiter
    # and call using the current mu2 and kappa estimates
    
    tryCatch(expr = {get_ebci_pvals(thetashrunk=thetashrunk, sigma=sigma, web=web, 
                                    mu2=cur_ebci_params$mu2['estimate'], kappa=cur_ebci_params$kappa['estimate'], 
                                    alpha_threshold=alpha_threshold, maxiter=maxiter)
    }, 
    error = function(e) {NA})
    
  }
  
  # # loop sequential
  # ebci_pvals_ = rep(NA, times = nrow(cur_shrinkage_results))
  # for(i in 1:100) {
  #   ebci_pvals[i] = inner_func(thetashrunk=cur_shrinkage_results$shrunk_value[i], sigma=cur_shrinkage_results$se[i], web=cur_shrinkage_results$w_eb[i])
  # }
  # 
  # 
  # 
  # # mapply sequential
  # ebci_pvals_ = mapply(FUN = inner_func, 
  #                      thetashrunk=cur_shrinkage_results$shrunk_value[1:100], 
  #                      sigma=cur_shrinkage_results$se[1:100], 
  #                      web=cur_shrinkage_results$w_eb[1:100])
  
  
  
  # future_mapply parallel
  ebci_pvals_ = future.apply::future_mapply(FUN = inner_func, 
                                            thetashrunk=cur_shrinkage_results$shrunk_value, 
                                            sigma=cur_shrinkage_results$se, 
                                            web=cur_shrinkage_results$w_eb)
  
  return(ebci_pvals_)
}

