# simulation for Zero-Inflated FAB p-values




# \theta_i  
#            = 0              wp     \pi_0 
#         \sim N(\mu, \tau^2) wp 1 - \pi_0


#' simulate N groups each of size G with variance s (singular val or vec)
#' where each Y_ij = \hat \theta_ij \sim N(\theta_i, s_i^2)
#' and \theta_i  
#'            = 0              wp     \pi_0 
#'         \sim N(\mu, \tau^2) wp 1 - \pi_0
#' @param N (integer) the number of groups
#' @param G (integer) the number of observations per group
#' @param s (numeric) >= 0, the sd of Y_ij | \theta_i
#' @param mu (numeric) the mean of \theta_i
#' @param tau (numeric) >= 0, the sd of \theta_i
#' @param pi0 (numeric) in [0,1], the probability of \theta_i = 0
#' @return dataframe with N x G rows and the following columns:
#'    i, j, Y, s, theta
sim_values <- function(N, G, s, mu, tau, pi0) {
  # # testing
  # N = 50
  # G = 2
  # s = 1
  # s = rbeta(N, 4, 2) + .1
  # mu = 1
  # tau = 2
  # pi0 = .3
  
  # sim thetas
  zero_mask = rbinom(N, size = 1, prob = 1-pi0) # 1 wp (1-pi0)
  theta = zero_mask * rnorm(N, mean = mu, sd = tau)
  
  
  # sim Y
  Y = rep(NA, N * G)
  if(length(s) == 1) {
    for(i in 1:N) {
      Y[(G*(i-1) + 1):(G*i)] = rnorm(G, mean = theta[i], s = s)
    }
  } else if(length(s) == N) {
    for(i in 1:N) {
      Y[(G*(i-1) + 1):(G*i)] = rnorm(G, mean = theta[i], s = s[i])
    }
  } else {
    return("Bad s input (must be single numeric val or vector of length G)")
  }

  
  
  return(data.frame(i = rep(1:N, each  = G), 
                    j = rep(1:G, times = N), 
                    Y = Y, 
                    s = if(length(s) == 1) rep(s, times = N*G) else rep(s, each = G), 
                    theta = theta))
}






#' estimate the prior G's parameters using an EM algorithm
#' Consider each i=1,...,N belonging to a Z=0 or Z=1 group
#' when Z=0, theta=0, and when Z=1, theta\sim N(\mu, \tau^2)
#' 
#' 
#' @param Y (vector) of observed Y values
#' @param S (vector) of Y standard deviations
#' @param num_steps (integer) number of iterations for EM '
#'    (TODO: later make it adaptive,
#'           args --> max_steps, 
#'                    threshold (keep iterating until all params |est - prev est| < threshold))
#' @param initial_mu (numeric)  initial value for param mu
#' @param initial_tau (numeric) initial value for param tau
#' @param initial_pi0 (numeric) initial value for param pi0
est_G_params <- function(Y, S, 
                         num_steps,
                         initial_mu = 0,
                         initial_tau = 1,
                         initial_pi0 = .5) {
  # # test
  # set.seed(12345)
  # sim_df = sim_values(N=500, G=1, s=1, mu=1, tau=2, pi0=.3)
  # Y = sim_df$Y
  # S = sim_df$s
  # num_steps = 100
  # initial_mu = -10
  # initial_tau = 20
  # initial_pi0 = .6
  
  
  
  assertthat::assert_that(length(Y) == length(S), msg = 'Y and S args must be vectors of the same length')
  assertthat::assert_that((num_steps %% 1 == 0) && num_steps > 0, msg = 'K arg must be positive integer')
  
  params = data.frame(step = rep(NA, num_steps),
                      mu   = rep(NA, num_steps),
                      tau  = rep(NA, num_steps),
                      pi0  = rep(NA, num_steps))
  # initial params
  # cur_mu = 0 # does setting cur_mu = 0 make it indistinguishable? I think its ok bc of the sd diff, 
  # cur_tau = 1
  # cur_pi0 = .5 # .5
  # params[1, ] = c(1, cur_mu, cur_tau, cur_pi0)
  
  
  params[1, ] = c(1, initial_mu, initial_tau, initial_pi0)
  
  for(k in 2:num_steps) {
    
    # E step
    P =  mapply(FUN = est_G_params_Estep,
                y=Y, s=S,
                MoreArgs = list(mu =params[k-1,  'mu'], 
                                tau=params[k-1, 'tau'], 
                                pi0=params[k-1, 'pi0']))
    # hist(P)
    
    
    # M step
    M_step_result = est_G_params_Mstep(Y=Y, S=S, P=P)
    
    
    # update params
    # cur_mu  = M_step_result[1]
    # cur_tau = M_step_result[2]
    # cur_pi0 = M_step_result[3]
    
    params[k, ] = c(k, M_step_result[1], 
                       M_step_result[2], 
                       M_step_result[3])
  }
  
  
  # # display parameter estimates over iterations
  # p_mu  = ggplot(params, aes(x = step, y =  mu)) + geom_line() + labs(title = 'mu')
  # p_tau = ggplot(params, aes(x = step, y = tau)) + geom_line() + labs(title = 'tau')
  # p_pi0 = ggplot(params, aes(x = step, y = pi0)) + geom_line() + labs(title = 'pi0')
  # gridExtra::grid.arrange(p_mu, p_tau, p_pi0)
  # 
  # 
  # # final estimates
  # params[nrow(params), ]
  
  
  return(params)
  
}


#' EM algorithm's E Step
#' using the currently parameters, calculate P(Z=0|Y=y)
#' (later, to speed up a tiny bit, should parameterize as tau2=tau^2)
est_G_params_Estep <- function(y, s, mu, tau, pi0) {
  # calculate the P(Z=0|Y=y)
  # numerator: P(Y=y|Z=0)pi0 (cancel out a 1/sqrt(2pi) factor)
  
  numerator = pi0 * (1/s) * exp(-((y^2)/(2*s^2)))
  
  
  # denominator: P(Y=y|Z=0)pi0 + P(Y=y|Z=1)(1-pi0)
  den1 = 1/(sqrt(tau^2 + s^2))
  den2 = exp(- (
              (tau^2*y^2 + s^2*mu^2 - (tau^2*y+s^2*mu)^2/(tau^2+s^2)) 
                    / 
              (2*s^2*tau^2)
            ))
  denominator = numerator + (1-pi0) * den1 * den2
  
  return(numerator / denominator)
}


#' EM algorithm's M Step
#' using the currently estimated P(Z=0|Y=y), estimate parameters
#' @param Y (vector) of Y values                 (same length as P & S)
#' @param S (vector) of Y's standard deviations  (same length as Y & P)
#' @param P (vector) of probabilities P(Z=0|Y=y) (same length as Y & S)
est_G_params_Mstep <- function(Y, S, P) {
  normalization = sum(1-P)
  
  mu = sum((1 - P)*Y) / normalization
  
  tau = (sum((1-P) * ((Y-mu)^2 - S*S))/normalization) |> sqrt()
  
  pi0 = mean(P)
  
  return(c(mu, tau, pi0))
}





#' Calculate the test statistic (when G=1, 1 obs per group):
#' T_j = \int L(\theta_j) dP_G(\theta_j)
#'       -------------------------------
#'             L(\theta_j0)
#' (likelihood for G>1 will look slightly diff but essentially the same 
#'  e.g. prob use average and scale the s value)
#' @param Yj (numeric) Yj value (observed value/estimated theta)
#' @param sj (numeric) sj value (sd of Y|theta)
#' @param thetaj0 (numeric) H_0 theta value testing against (typically 0)
#' @param mu (numeric) estimated mu parameter
#' @param tau (numeric) estimated tau parameter
#' @param pi0 (numeric) estimated pi0 parameter
calc_T_stat_one <- function(Yj, sj, thetaj0, mu, tau, pi0) {
  # # test
  # set.seed(12345)
  # sim_df = sim_values(N=50, G=1, s=1, mu=1, tau=2, pi0=.3)
  # 
  # Yj = sim_df[1, 'Y']
  # sj = sim_df[1, 's']
  # thetaj0 = 0
  # # test, use true params
  # mu = 1
  # tau = 2
  # pi0 = .3
  
  # Calculate Denominator
  denominator = (1/sj) * exp(-((Yj-thetaj0)^2)/(2*sj^2))
  
  # Calculate Numerator
  num1 = pi0 * (1/sj) * exp(-(Yj^2)/(2*sj^2))
  
  num2a = -(tau^2 * Yj^2 + sj^2 * mu^2 - (tau^2 * Yj + sj^2 * mu)^2 / (tau^2 + sj^2)) / (2 * sj^2 * tau^2)
  num2b = (1/sqrt(tau^2 + sj^2))
  num2 = (1-pi0) * exp(num2a) * num2b
  
  numerator = num1 + num2
  
  # T statistic
  return(numerator / denominator)
}



#' Used to estimate the p-value by drawing Y_j samples from the Null
#' For observed T statistic T_j, the p-value is 
#'      p(\theta_{j0}) = P(T_j^* > T_j)
#'       where T_j^* has distribution according to the null N(\theta_{j0}, \sigma_j^2) 
#' @param B (integer) number of repetitions 
#'         (TODO: could make it adaptive (keep going until change < threshold))
#' @param sj (numeric) sj value (sd of Y|theta)
#' @param thetaj0 (numeric) H_0 theta value testing against (typically 0)
#' @param mu (numeric) estimated mu parameter
#' @param tau (numeric) estimated tau parameter
#' @param pi0 (numeric) estimated pi0 parameter
sample_tstat_null <- function(B, sj, thetaj0, mu, tau, pi0) {
  # # debug
  # set.seed(12345)
  # sim_df = sim_values(N=50, G=1, s=1, mu=1, tau=2, pi0=.3)
  # 
  # Yj = sim_df[1, 'Y']
  # sj = sim_df[1, 's']
  # thetaj0 = 0
  # # test, use true params
  # mu = 1
  # tau = 2
  # pi0 = .3
  # 
  # B = 10000
  
  
  
  # draw Ys
  Y = rnorm(n = B, mean = thetaj0, sd = sj)
  # calculate tstats when H0:thetaj=thetaj0
  tstat_null = sapply(X = Y, FUN = calc_T_stat_one,
                      sj=sj, thetaj0=thetaj0, mu=mu, tau=tau, pi0=pi0)
  return(tstat_null)
}




#' Estimate the p-value by drawing Y_j samples using the estimated prior \hat{G} 
#' For observed T statistic T_j, the p-value is 
#'      p(\theta_{j0}) = P(T_j^* > T_j)
#'       where T_j^* has distribution according to the null N(\theta_{j0}, \sigma_j^2) 
#' @param tj (numeric) observed T_j statistic value
#' @param B (integer) number of repetitions 
#'         (TODO: could make it adaptive (keep going until change < threshold))
#' @param sj (numeric) sj value (sd of Y|theta)
#' @param thetaj0 (numeric) H_0 theta value testing against (typically 0)
#' @param mu (numeric) estimated mu parameter
#' @param tau (numeric) estimated tau parameter
#' @param pi0 (numeric) estimated pi0 parameter
est_pval <- function(tj, B, sj, thetaj0, mu, tau, pi0) {
  # # debug
  # set.seed(12345)
  # sim_df = sim_values(N=50, G=1, s=1, mu=1, tau=2, pi0=.3)
  # 
  # Yj = sim_df[1, 'Y']
  # sj = sim_df[1, 's']
  # thetaj0 = 0
  # # test, use true params
  # mu = 1
  # tau = 2
  # pi0 = .3
  # 
  # t = calc_T_stat_one(Yj=Yj, sj = sj, thetaj0 = thetaj0,
  #                     # test, use true params
  #                     mu = mu, tau = tau, pi0 = pi0)
  # B = 10000
  
  
  # do in another fn to reuse values
  # # draw Ys
  # Y = rnorm(n = B, mean = thetaj0, sd = sj)
  # # calculate tstats when H0:thetaj=thetaj0
  # tstat_null = sapply(X = Y, FUN = calc_T_stat_one,
  #                     sj=sj, thetaj0=thetaj0, mu=mu, tau=tau, pi0=pi0)
  
  
  tstat_null = sample_tstat_null(B, sj, thetaj0, mu, tau, pi0)
  # estimate p-value (empirical distribution)
  phat = mean(tstat_null > tj)
  
  # TODO: could make this adaptive to B ... not for now, just choose a large enough B
  return(phat)
}



#' perform simulation and plot results under specified sim settings
#' 
#' @param save_folder (character)
#' @param N (integer)
#' @param G (integer)
#' @param s (numeric or numeric  vector)
#' @param mu (numeric)
#' @param tau (numeric)
#' @param pi0 (numeric)
#' @param thetaj0 (numeric)
#' @param B (integer)
sim_and_plot <- function(save_folder,
                         N, G, s, mu, tau, pi0,
                         thetaj0, B) {
  dir.create(save_folder, showWarnings = TRUE)
  
  sim_df = sim_values(N=N, G=G, s=s, mu=mu, tau=tau, pi0=pi0)
  # (sim_df$theta == 0) |> mean()
  
  
  # Standard t-statistics + p-values (compare Y to N(0, s^2))
  Tstat_standard =              sim_df[ ,'Y'] / sim_df[ ,'s']
  p_standard     = 2*pnorm(-abs(sim_df[ ,'Y'] / sim_df[ ,'s']))
  
  
  
  # True parameters
  # calculate the T statistic using true parameters
  Tstat_trueparams = mapply(FUN = calc_T_stat_one,
                            Yj = sim_df[ ,'Y'],
                            sj = sim_df[ , 's'], 
                            MoreArgs = list(thetaj0 = thetaj0,
                                            mu      = mu, 
                                            tau     = tau,
                                            pi0     = pi0))
  
  # estimate the p-values (using simulation) using true parameters
  phat_trueparams = mapply(FUN = est_pval,
                           tj = Tstat_trueparams,
                           sj = sim_df[ , 's'], 
                           MoreArgs = list(B       = B,
                                           thetaj0 = thetaj0,
                                           mu      = mu, 
                                           tau     = tau,
                                           pi0     = pi0))
  
  # Estimated (using EM alg) parameters
  # estimate the prior, G, parameters
  estparams = est_G_params(Y = sim_df[ , 'Y'], 
                           S = sim_df[ , 's'], 
                           num_steps=100,
                           initial_mu = 0,
                           initial_tau = 1,
                           initial_pi0 = .5)
  
  estmu  = estparams[nrow(estparams),  'mu']
  esttau = estparams[nrow(estparams), 'tau']
  estpi0 = estparams[nrow(estparams), 'pi0']
  
  
  # calculate the T statistic using true parameters
  Tstat_estparams = mapply(FUN = calc_T_stat_one,
                           Yj = sim_df[ ,'Y'],
                           sj = sim_df[ , 's'], 
                           MoreArgs = list(thetaj0 = thetaj0,
                                           mu      = estmu, 
                                           tau     = esttau,
                                           pi0     = estpi0))
  
  # estimate the p-values (using simulation) using true parameters
  phat_estparams = mapply(FUN = est_pval,
                          tj = Tstat_trueparams,
                          sj = sim_df[ , 's'], 
                          MoreArgs = list(B       = B,
                                          thetaj0 = thetaj0,
                                          mu      = estmu, 
                                          tau     = esttau,
                                          pi0     = estpi0))
  
  
  
  # combine to dataframe
  simres_df = rbind(cbind(sim_df,
                          data.frame(method = 'standard',
                                     Tstat  = Tstat_standard,
                                     pval   = p_standard)),
                    cbind(sim_df,
                          data.frame(method = 'ZIFAB_trueG',
                                     Tstat  = Tstat_trueparams,
                                     pval   = phat_trueparams)),
                    cbind(sim_df,
                          data.frame(method = 'ZIFAB_EMG',
                                     Tstat  = Tstat_estparams,
                                     pval   = phat_estparams))
  )
  
  simres_df$method = factor(simres_df$method,
                            levels = c('standard', 'ZIFAB_trueG', 'ZIFAB_EMG'),
                            labels = c('standard', 'ZIFAB (true G)', 'ZIFAB (EM G)'))
  method_colors = c('blue', 'orange', 'brown')
  method_colors = RColorBrewer::brewer.pal(3, name = 'Set2')
  
  
  # test EM estimation performance (same parameters, change sample size)
  
  # Estimated (using EM alg) parameters
  sim_df = sim_values(N=5000, G=G, s=s, mu=mu, tau=tau, pi0=pi0)
  Ns = c(25, 50, 100, 250, 500, 1000, 1500, 2000, 2500, 3000, 4000, 5000)
  EM_est = NULL
  for(n in Ns) {
    t0 = Sys.time()
    EMparams = est_G_params(Y = sim_df[1:n, 'Y'], 
                            S = sim_df[1:n, 's'], 
                            num_steps=100,
                            initial_mu = 0,
                            initial_tau = 1,
                            initial_pi0 = .5)
    t1 = Sys.time()
    EM_est = rbind(EM_est, 
                   data.frame(n   = n,
                              mu  = EMparams[nrow(EMparams),  'mu'],
                              tau = EMparams[nrow(EMparams), 'tau'],
                              pi0 = EMparams[nrow(EMparams), 'pi0'], 
                              time = difftime(t1, t0, units = 'sec')))
  }
  
  
  # display EM parameter estimates over iterations
  p_mu_samplesize  = ggplot() + geom_line(data=EM_est, aes(x = n, y =  mu)) + labs(title = 'mu')  + geom_hline(aes(yintercept =  mu), color = 'orange')
  p_tau_samplesize = ggplot() + geom_line(data=EM_est, aes(x = n, y = tau)) + labs(title = 'tau') + geom_hline(aes(yintercept = tau), color = 'orange')
  p_pi0_samplesize = ggplot() + geom_line(data=EM_est, aes(x = n, y = pi0)) + labs(title = 'pi0') + geom_hline(aes(yintercept = pi0), color = 'orange')
  grob <- gridExtra::grid.arrange(p_mu_samplesize, p_tau_samplesize, p_pi0_samplesize)
  ggsave(sprintf('%s/estparam_EM_bysamplesize.pdf', save_folder), grob, width = 4, height = 6)
  
  
  
  # Visualizations
  
  # display EM parameter estimates over iterations
  p_mu  = ggplot() + geom_line(data=estparams, aes(x = step, y =  mu)) + labs(title = 'mu')  + geom_hline(aes(yintercept =  mu), color = 'orange')
  p_tau = ggplot() + geom_line(data=estparams, aes(x = step, y = tau)) + labs(title = 'tau') + geom_hline(aes(yintercept = tau), color = 'orange')
  p_pi0 = ggplot() + geom_line(data=estparams, aes(x = step, y = pi0)) + labs(title = 'pi0') + geom_hline(aes(yintercept = pi0), color = 'orange')
  grob <- gridExtra::grid.arrange(p_mu, p_tau, p_pi0)
  ggsave(sprintf('%s/estparam_EM_byiterations.pdf', save_folder), grob, width = 4, height = 6)
  
  
  
  
  
  
  
  
  # EDA
  
  # Prior: Histogram of theta
  p1 = ggplot(sim_df, aes(x = theta)) +
    geom_histogram(binwidth = .5) +
    scale_x_continuous(limits = c(-5, 6)) +
    labs(title = 'Histogram of Theta')
  p1 
  ggsave(filename = sprintf('%s/theta_hist.pdf', save_folder), height = 4, width = 4)
  
  
  # Posterior: Histogram of Y
  p2 = ggplot(sim_df, aes(x = Y)) +
    geom_histogram(binwidth = .5) +
    scale_x_continuous(limits = c(-5, 6)) +
    labs(title = 'Histogram of Y')
  p2
  ggsave(filename = sprintf('%s/Y_hist.pdf', save_folder), height = 4, width = 4)
  
  
  # Results
  
  
  # Points: Test Statistic vs true Theta
  p3 = ggplot(simres_df |> dplyr::filter(method != 'standard'), aes(x = theta, y = Tstat, group = method, color = method)) +
    geom_point(alpha = .7, size = .6) +
    scale_y_log10() +
    # scale_x_continuous(limits = c(-4, 4)) +
    scale_color_manual(values = method_colors[-1]) +
    labs(title = 'Test Statistic vs Theta',
         x ='Theta', y = 'Test Statistic (log scale)')
  p3
  ggsave(filename = sprintf('%s/tstat_theta_pts.pdf', save_folder), height = 4, width = 6)
  
  
  
  
  # Points: Test Statistic vs Y
  p4 = ggplot(simres_df |> dplyr::filter(method != 'standard'), 
              aes(x = Y, y = Tstat, group = method, color = method)) +
    geom_point(alpha = .7, size = 1) +
    scale_y_log10() +
    # scale_x_continuous(limits = c(-4, 4)) +
    scale_color_manual(values = method_colors[-1]) +
    labs(title = 'Test Statistic vs Y',
         x ='Y', y = 'Test Statistic (log scale)')
  p4
  ggsave(filename = sprintf('%s/tstat_Y_pts.pdf', save_folder), height = 4, width = 6)
  
  
  
  # Points: Estimated p-value vs true Theta
  p5 = ggplot(simres_df, aes(x = theta, y = pval, group = method, color = method)) +
    geom_point(alpha = .7, size = 1) +
    scale_color_manual(values = method_colors) +
    scale_x_continuous(limits = c(-4, 4)) +
    labs(title = 'p-value vs Theta',
         x ='Theta', y = 'p-value')
  p5
  ggsave(filename = sprintf('%s/pval_theta_pts.pdf', save_folder), height = 4, width = 6)
  
  
  # Points: Estimated p-value vs Y
  p6 = ggplot(simres_df, 
              aes(x = Y, y = pval, group = method, color = method)) +
    geom_point(alpha = .7, size = 1) +
    scale_color_manual(values = method_colors) +
    scale_x_continuous(limits = c(-4, 4)) +
    labs(title = 'p-value vs Y',
         x ='Y', y = 'p-value')
  p6
  ggsave(filename = sprintf('%s/pval_Y_pts.pdf', save_folder), height = 4, width = 6)
  
  
  
  # For theta=0, histogram of pvals
  p7 = ggplot(simres_df |> dplyr::filter(theta == 0),
              aes(x = pval, group = method, fill = method)) +
    geom_histogram(breaks = seq(from = -.1, to = 1.1, by = .05), position = 'dodge') +
    scale_fill_manual(values = method_colors) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    labs(title = 'Histogram of p-values',
         subtitle = '(when Theta=0)',
         x = 'p-value')
  p7
  ggsave(filename = sprintf('%s/pval_hist_0.pdf', save_folder), height = 4, width = 6)
  
  # For theta!=0, histogram of pvals
  p8 = ggplot(simres_df |> dplyr::filter(theta != 0),
              aes(x = pval, group = method, fill = method)) +
    geom_histogram(breaks = seq(from = -.1, to = 1.1, by = .05), position = 'dodge') +
    scale_fill_manual(values = method_colors) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    labs(title = 'Histogram of p-values',
         subtitle = '(when Theta!=0)',
         x = 'p-value')
  p8
  ggsave(filename = sprintf('%s/pval_hist_not0.pdf', save_folder), height = 4, width = 6)
  
  
  
  # For theta=0, qqplot of phat
  p9 = ggplot(simres_df |> dplyr::filter(theta == 0) |> 
                dplyr::group_by(method) |>
                dplyr::arrange(pval) |> 
                dplyr::mutate(rej_rate = (1:dplyr::n())/dplyr::n()),
              aes(x = pval, y = rej_rate, group = method, color = method)) +
    geom_abline(aes(slope = 1, intercept = 0)) +
    geom_point(alpha = .9, size = .6) +
    # scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_color_manual(values = method_colors) +
    labs(title = 'Rejection Rate by p-value (qqplot)',
         subtitle = '(when Theta=0)',
         x = 'p-value', y = 'Rejection Rate')
  p9
  ggsave(filename = sprintf('%s/pval_qqplot_0.pdf', save_folder), height = 4, width = 6)
  
  
  # For theta!=0, qqplot of phat
  p10 = ggplot(simres_df |> dplyr::filter(theta != 0) |> 
                 dplyr::group_by(method) |>
                 dplyr::arrange(pval) |> 
                 dplyr::mutate(rej_rate = (1:dplyr::n())/dplyr::n()),
               aes(x = pval, y = rej_rate, group = method, color = method)) +
    geom_abline(aes(slope = 1, intercept = 0)) +
    geom_point(alpha = .9, size = .6) +
    # scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_color_manual(values = method_colors) +
    labs(title = 'Rejection Rate by p-value (qqplot)',
         subtitle = '(when Theta!=0)',
         x = 'p-value', y = 'Rejection Rate')
  p10
  ggsave(filename = sprintf('%s/pval_qqplot_not0.pdf', save_folder), height = 4, width = 6)
  
  
  
  
  grob <- gridExtra::arrangeGrob(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10,
                                 p_mu, p_tau, p_pi0,
                                 p_mu_samplesize, p_tau_samplesize, p_pi0_samplesize,
                                 # ncol = 2, 
                                 layout_matrix = matrix(c(14, 15, 16, NA,
                                                          11, 12, 13,NA,
                                                          1, 2, NA, NA,
                                                          3, 4, 5, 6,
                                                          7, 8, 9, 10), byrow = TRUE, nrow = 5))
  ggsave(sprintf('%s/all.pdf', save_folder), grob, width = 24, height = 18)
  
  
  

  
  
  
  
  
}




library(ggplot2)
library(cowplot)

ggplot2::set_theme(theme_cowplot() +
                     theme(title = element_text(hjust = .5)))


save_folder_overall = "C:/Users/Cathe/Documents/School/Genetic Hypothesis Testing using Negative Controls/multilevel-treatment-multivariable-outcome-github/multilevel-treatment-multivariable-outcome/plots/ZIFAB/sim/"


# simulation with not so well separated groups (leads to poor estimation of G --> worse results)
set.seed(12345)
sim_and_plot(save_folder = paste0(save_folder_overall, "A/"), 
             N=1000, 
             G=1,
             s=1,
             mu=1,
             tau=2,
             pi0=.3,
             thetaj0=0,
             B=10000)


# simulation with  well separated groups (leads to good estimation of G --> good results. but standard is still just ideal, maybe we can get more power??)
set.seed(12345)
sim_and_plot(save_folder = paste0(save_folder_overall, "B/"), 
             N=1000, 
             G=1,
             s=.5,
             mu=3,
             tau=.5,
             pi0=.3,
             thetaj0=0,
             B=10000)


















# old simulation
if(F) {
  

set.seed(12345)
N = 1000
G=1
s=1
mu=1
tau=2
pi0=.3

thetaj0 = 0
B = 10000 # #samples used for estimating p-value

sim_df = sim_values(N=N, G=G, s=s, mu=mu, tau=tau, pi0=pi0)
# (sim_df$theta == 0) |> mean()


# Standard t-statistics + p-values (compare Y to N(0, s^2))
Tstat_standard =              sim_df[ ,'Y'] / sim_df[ ,'s']
p_standard     = 2*pnorm(-abs(sim_df[ ,'Y'] / sim_df[ ,'s']))



# True parameters
# calculate the T statistic using true parameters
Tstat_trueparams = mapply(FUN = calc_T_stat_one,
                          Yj = sim_df[ ,'Y'],
                          sj = sim_df[ , 's'], 
                          MoreArgs = list(thetaj0 = thetaj0,
                                          mu      = mu, 
                                          tau     = tau,
                                          pi0     = pi0))

# estimate the p-values (using simulation) using true parameters
phat_trueparams = mapply(FUN = est_pval,
                         tj = Tstat_trueparams,
                         sj = sim_df[ , 's'], 
                         MoreArgs = list(B       = B,
                                         thetaj0 = thetaj0,
                                         mu      = mu, 
                                         tau     = tau,
                                         pi0     = pi0))

# Estimated (using EM alg) parameters
# estimate the prior, G, parameters
estparams = est_G_params(Y = sim_df[ , 'Y'], 
                         S = sim_df[ , 's'], 
                         num_steps=100,
                         initial_mu = 0,
                         initial_tau = 1,
                         initial_pi0 = .5)

estmu  = estparams[nrow(estparams),  'mu']
esttau = estparams[nrow(estparams), 'tau']
estpi0 = estparams[nrow(estparams), 'pi0']


# calculate the T statistic using true parameters
Tstat_estparams = mapply(FUN = calc_T_stat_one,
                          Yj = sim_df[ ,'Y'],
                          sj = sim_df[ , 's'], 
                          MoreArgs = list(thetaj0 = thetaj0,
                                          mu      = estmu, 
                                          tau     = esttau,
                                          pi0     = estpi0))

# estimate the p-values (using simulation) using true parameters
phat_estparams = mapply(FUN = est_pval,
                         tj = Tstat_trueparams,
                         sj = sim_df[ , 's'], 
                         MoreArgs = list(B       = B,
                                         thetaj0 = thetaj0,
                                         mu      = estmu, 
                                         tau     = esttau,
                                         pi0     = estpi0))



# combine to dataframe
simres_df = rbind(cbind(sim_df,
                     data.frame(method = 'standard',
                                Tstat  = Tstat_standard,
                                pval   = p_standard)),
               cbind(sim_df,
                     data.frame(method = 'ZIFAB_trueG',
                                Tstat  = Tstat_trueparams,
                                pval   = phat_trueparams)),
               cbind(sim_df,
                     data.frame(method = 'ZIFAB_EMG',
                                Tstat  = Tstat_estparams,
                                pval   = phat_estparams))
               )

simres_df$method = factor(simres_df$method,
                          levels = c('standard', 'ZIFAB_trueG', 'ZIFAB_EMG'),
                          labels = c('standard', 'ZIFAB (true G)', 'ZIFAB (EM G)'))
method_colors = c('blue', 'orange', 'brown')
method_colors = RColorBrewer::brewer.pal(3, name = 'Set2')

# Visualizations


library(ggplot2)
library(cowplot)

ggplot2::set_theme(theme_cowplot() +
                     theme(title = element_text(hjust = .5)))




save_folder = "C:/Users/Cathe/Documents/School/Genetic Hypothesis Testing using Negative Controls/multilevel-treatment-multivariable-outcome-github/multilevel-treatment-multivariable-outcome/plots/ZIFAB/sim/"




# display EM parameter estimates over iterations
p_mu  = ggplot() + geom_line(data=estparams, aes(x = step, y =  mu)) + labs(title = 'mu')  + geom_hline(aes(yintercept =  mu), color = 'orange')
p_tau = ggplot() + geom_line(data=estparams, aes(x = step, y = tau)) + labs(title = 'tau') + geom_hline(aes(yintercept = tau), color = 'orange')
p_pi0 = ggplot() + geom_line(data=estparams, aes(x = step, y = pi0)) + labs(title = 'pi0') + geom_hline(aes(yintercept = pi0), color = 'orange')
grob <- gridExtra::grid.arrange(p_mu, p_tau, p_pi0)
ggsave(sprintf('%s/estparam_EM_byiterations.pdf', save_folder), grob, width = 4, height = 6)








# EDA

# Prior: Histogram of theta
p1 = ggplot(sim_df, aes(x = theta)) +
  geom_histogram(binwidth = .5) +
  scale_x_continuous(limits = c(-5, 6)) +
  labs(title = 'Histogram of Theta')
p1 
ggsave(filename = sprintf('%s/theta_hist.pdf', save_folder), height = 4, width = 4)


# Posterior: Histogram of Y
p2 = ggplot(sim_df, aes(x = Y)) +
  geom_histogram(binwidth = .5) +
  scale_x_continuous(limits = c(-5, 6)) +
  labs(title = 'Histogram of Y')
p2
ggsave(filename = sprintf('%s/Y_hist.pdf', save_folder), height = 4, width = 4)


# Results


# Points: Test Statistic vs true Theta
p3 = ggplot(simres_df |> dplyr::filter(method != 'standard'), aes(x = theta, y = Tstat, group = method, color = method)) +
  geom_point(alpha = .7, size = .6) +
  scale_y_log10() +
  # scale_x_continuous(limits = c(-4, 4)) +
  scale_color_manual(values = method_colors[-1]) +
  labs(title = 'Test Statistic vs Theta',
       x ='Theta', y = 'Test Statistic (log scale)')
p3
ggsave(filename = sprintf('%s/tstat_theta_pts.pdf', save_folder), height = 4, width = 6)




# Points: Test Statistic vs Y
p4 = ggplot(simres_df |> dplyr::filter(method != 'standard'), 
            aes(x = Y, y = Tstat, group = method, color = method)) +
  geom_point(alpha = .7, size = 1) +
  scale_y_log10() +
  # scale_x_continuous(limits = c(-4, 4)) +
  scale_color_manual(values = method_colors[-1]) +
  labs(title = 'Test Statistic vs Y',
       x ='Y', y = 'Test Statistic (log scale)')
p4
ggsave(filename = sprintf('%s/tstat_Y_pts.pdf', save_folder), height = 4, width = 6)



# Points: Estimated p-value vs true Theta
p5 = ggplot(simres_df, aes(x = theta, y = pval, group = method, color = method)) +
  geom_point(alpha = .7, size = 1) +
  scale_color_manual(values = method_colors) +
  scale_x_continuous(limits = c(-4, 4)) +
  labs(title = 'p-value vs Theta',
       x ='Theta', y = 'p-value')
p5
ggsave(filename = sprintf('%s/pval_theta_pts.pdf', save_folder), height = 4, width = 6)


# Points: Estimated p-value vs Y
p6 = ggplot(simres_df, 
            aes(x = Y, y = pval, group = method, color = method)) +
  geom_point(alpha = .7, size = 1) +
  scale_color_manual(values = method_colors) +
  scale_x_continuous(limits = c(-4, 4)) +
  labs(title = 'p-value vs Y',
       x ='Y', y = 'p-value')
p6
ggsave(filename = sprintf('%s/pval_Y_pts.pdf', save_folder), height = 4, width = 6)



# For theta=0, histogram of pvals
p7 = ggplot(simres_df |> dplyr::filter(theta == 0),
            aes(x = pval, group = method, fill = method)) +
  geom_histogram(breaks = seq(from = -.1, to = 1.1, by = .05), position = 'dodge') +
  scale_fill_manual(values = method_colors) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  labs(title = 'Histogram of p-values',
       subtitle = '(when Theta=0)',
       x = 'p-value')
p7
ggsave(filename = sprintf('%s/pval_hist_0.pdf', save_folder), height = 4, width = 6)

# For theta!=0, histogram of pvals
p8 = ggplot(simres_df |> dplyr::filter(theta != 0),
            aes(x = pval, group = method, fill = method)) +
  geom_histogram(breaks = seq(from = -.1, to = 1.1, by = .05), position = 'dodge') +
  scale_fill_manual(values = method_colors) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  labs(title = 'Histogram of p-values',
       subtitle = '(when Theta!=0)',
       x = 'p-value')
p8
ggsave(filename = sprintf('%s/pval_hist_not0.pdf', save_folder), height = 4, width = 6)



# For theta=0, qqplot of phat
p9 = ggplot(simres_df |> dplyr::filter(theta == 0) |> 
              dplyr::group_by(method) |>
              dplyr::arrange(pval) |> 
              dplyr::mutate(rej_rate = (1:dplyr::n())/dplyr::n()),
            aes(x = pval, y = rej_rate, group = method, color = method)) +
  geom_abline(aes(slope = 1, intercept = 0)) +
  geom_point(alpha = .9, size = .6) +
  # scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_color_manual(values = method_colors) +
  labs(title = 'Rejection Rate by p-value (qqplot)',
       subtitle = '(when Theta=0)',
       x = 'p-value', y = 'Rejection Rate')
p9
ggsave(filename = sprintf('%s/pval_qqplot_0.pdf', save_folder), height = 4, width = 6)


# For theta!=0, qqplot of phat
p10 = ggplot(simres_df |> dplyr::filter(theta != 0) |> 
              dplyr::group_by(method) |>
              dplyr::arrange(pval) |> 
              dplyr::mutate(rej_rate = (1:dplyr::n())/dplyr::n()),
            aes(x = pval, y = rej_rate, group = method, color = method)) +
  geom_abline(aes(slope = 1, intercept = 0)) +
  geom_point(alpha = .9, size = .6) +
  # scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_color_manual(values = method_colors) +
  labs(title = 'Rejection Rate by p-value (qqplot)',
       subtitle = '(when Theta!=0)',
       x = 'p-value', y = 'Rejection Rate')
p10
ggsave(filename = sprintf('%s/pval_qqplot_not0.pdf', save_folder), height = 4, width = 6)




grob <- gridExtra::arrangeGrob(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10,
                               p_mu, p_tau, p_pi0,
                               # ncol = 2, 
                               layout_matrix = matrix(c(11, 12, 13,NA,
                                                        1, 2, NA, NA,
                                                        3, 4, 5, 6,
                                                        7, 8, 9, 10), byrow = TRUE, nrow = 4))
ggsave(sprintf('%s/all.pdf', save_folder), grob, width = 24, height = 12)











# test EM estimation performance (same parameters, change sample size)

# Estimated (using EM alg) parameters
set.seed(12345)
sim_df = sim_values(N=5000, G=G, s=s, mu=mu, tau=tau, pi0=pi0)
Ns = c(25, 50, 100, 250, 500, 1000, 1500, 2000, 2500, 3000, 4000, 5000)
EM_est = NULL
for(n in Ns) {
  t0 = Sys.time()
  EMparams = est_G_params(Y = sim_df[1:n, 'Y'], 
                           S = sim_df[1:n, 's'], 
                           num_steps=100,
                           initial_mu = 0,
                           initial_tau = 1,
                           initial_pi0 = .5)
  t1 = Sys.time()
  EM_est = rbind(EM_est, 
                 data.frame(n   = n,
                            mu  = EMparams[nrow(EMparams),  'mu'],
                            tau = EMparams[nrow(EMparams), 'tau'],
                            pi0 = EMparams[nrow(EMparams), 'pi0'], 
                            time = difftime(t1, t0, units = 'sec')))
}


# display EM parameter estimates over iterations
p_mu_samplesize  = ggplot() + geom_line(data=EM_est, aes(x = n, y =  mu)) + labs(title = 'mu')  + geom_hline(aes(yintercept =  mu), color = 'orange')
p_tau_samplesize = ggplot() + geom_line(data=EM_est, aes(x = n, y = tau)) + labs(title = 'tau') + geom_hline(aes(yintercept = tau), color = 'orange')
p_pi0_samplesize = ggplot() + geom_line(data=EM_est, aes(x = n, y = pi0)) + labs(title = 'pi0') + geom_hline(aes(yintercept = pi0), color = 'orange')
grob <- gridExtra::grid.arrange(p_mu_samplesize, p_tau_samplesize, p_pi0_samplesize)
ggsave(sprintf('%s/estparam_EM_bysamplesize.pdf', save_folder), grob, width = 4, height = 6)







}











# Old sim
if(F) {
  


set.seed(12345)
N = 500
G=1
s=1
mu=1
tau=2
pi0=.3

thetaj0 = 0
B = 10000 # #samples used for estimating p-value

sim_df = sim_values(N=N, G=G, s=s, mu=mu, tau=tau, pi0=pi0)
# (sim_df$theta == 0) |> mean()


# if same \theta_j0 AND same sj for all j, we can just reuse the estimated tstat null distn... but sj will NOT be the same... so don't do this
tstat_null = sample_tstat_null(B       = B, 
                               sj      = s,
                               thetaj0 = thetaj0,
                               mu      = mu, # test, use true params
                               tau     = tau,
                               pi0     = pi0)

Tstat_trueparams = rep(NA, nrow(sim_df))
Tstat_EMparams   = rep(NA, nrow(sim_df))
phat_trueparams  = rep(NA, nrow(sim_df))
phat_EMparams    = rep(NA, nrow(sim_df))
for(j in 1:nrow(sim_df)) {
  Tstat_trueparams[j] = calc_T_stat_one(
                             Yj = sim_df[j, 'Y'],
                             sj = sim_df[j, 's'],
                             thetaj0 = thetaj0,
                             mu      = mu, # test, use true params
                             tau     = tau,
                             pi0     = pi0)
  # same thetaj0 for all j
  phat[j] = mean(tstat_null > Tstat[j])
  # # code for potentially diff thetaj0 for every j
  # phat[j] = est_pval(t = Tstat[j],
  #                    B = B,
  #                    sj = sim_df[j, 's'],
  #                    thetaj0 = thetaj0,
  #                    mu      = mu, # test, use true params
  #                    tau     = tau,
  #                    pi0     = pi0)
  
}
hist(Tstat)
plot(sim_df$theta, Tstat |> log())
hist(phat)


sim_df$Tstat = Tstat
sim_df$phat = phat






library(ggplot2)
library(cowplot)

ggplot2::set_theme(theme_cowplot() +
                   theme(title = element_text(hjust = .5)))




save_folder = "C:/Users/Cathe/Documents/School/Genetic Hypothesis Testing using Negative Controls/multilevel-treatment-multivariable-outcome-github/multilevel-treatment-multivariable-outcome/plots/ZIFAB/sim/"






# EDA

# Prior: Histogram of theta
p1 = ggplot(sim_df, aes(x = theta)) +
  geom_histogram(binwidth = .5) +
  scale_x_continuous(limits = c(-5, 6)) +
  labs(title = 'Histogram of Theta')
p1 
ggsave(filename = sprintf('%s/theta_hist.pdf', save_folder), height = 4, width = 4)


# Posterior: Histogram of Y
p2 = ggplot(sim_df, aes(x = Y)) +
  geom_histogram(binwidth = .5) +
  scale_x_continuous(limits = c(-5, 6)) +
  labs(title = 'Histogram of Y')
p2
ggsave(filename = sprintf('%s/Y_hist.pdf', save_folder), height = 4, width = 4)


# Results


# Points: Test Statistic vs true Theta
p3 = ggplot(sim_df, aes(x = theta, y = Tstat)) +
  geom_point(alpha = .4, size = 1) +
  scale_y_log10() +
  labs(title = 'Test Statistic vs Theta',
       x ='Theta', y = 'Test Statistic (log scale)')
p3
ggsave(filename = sprintf('%s/tstat_theta_pts.pdf', save_folder), height = 4, width = 4)




# Points: Test Statistic vs Y
p4 = ggplot(sim_df, aes(x = Y, y = Tstat)) +
  geom_point(alpha = .4, size = 1) +
  scale_y_log10() +
  labs(title = 'Test Statistic vs Y',
       x ='Y', y = 'Test Statistic (log scale)')
p4
ggsave(filename = sprintf('%s/tstat_Y_pts.pdf', save_folder), height = 4, width = 4)



# Points: Estimated p-value vs true Theta
p5 = ggplot(sim_df, aes(x = theta, y = phat)) +
  geom_point(alpha = .4, size = 1) +
  labs(title = 'p-value vs Theta',
       x ='Theta', y = 'p-value')
p5
ggsave(filename = sprintf('%s/pval_theta_pts.pdf', save_folder), height = 4, width = 4)


# Points: Estimated p-value vs Y
p6 = ggplot(sim_df, aes(x = Y, y = phat)) +
  geom_point(alpha = .4, size = 1) +
  labs(title = 'p-value vs Y',
       x ='Y', y = 'p-value')
p6
ggsave(filename = sprintf('%s/pval_Y_pts.pdf', save_folder), height = 4, width = 4)



# For theta=0, histogram of phat
p7 = ggplot(sim_df |> dplyr::filter(theta == 0),
       aes(x = phat)) +
  geom_histogram(breaks = seq(from = -.1, to = 1.1, by = .025)) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  labs(title = 'Histogram of p-values',
       subtitle = '(when Theta=0)',
       x = 'p-value')
p7
ggsave(filename = sprintf('%s/pval_hist_0.pdf', save_folder), height = 4, width = 4)

# For theta!=0, histogram of phat
p8 = ggplot(sim_df |> dplyr::filter(theta != 0),
       aes(x = phat)) +
  geom_histogram(breaks = seq(from = -.1, to = 1.1, by = .025)) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  labs(title = 'Histogram of p-values',
       subtitle = '(when Theta!=0)',
       x = 'p-value')
p8
ggsave(filename = sprintf('%s/pval_hist_not0.pdf', save_folder), height = 4, width = 4)



# For theta=0, qqplot of phat
p9 = ggplot(sim_df |> dplyr::filter(theta == 0) |> 
         dplyr::arrange(phat) |> 
         dplyr::mutate(rej_rate = (1:dplyr::n())/dplyr::n()),
       aes(x = phat, y = rej_rate)) +
  geom_abline(aes(slope = 1, intercept = 0)) +
  geom_point() +
  # scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  labs(title = 'Rejection Rate by p-value (qqplot)',
       subtitle = '(when Theta=0)',
       x = 'p-value', y = 'Rejection Rate')
p9
ggsave(filename = sprintf('%s/pval_qqplot_0.pdf', save_folder), height = 4, width = 4)


# For theta!=0, qqplot of phat
p10 = ggplot(sim_df |> dplyr::filter(theta != 0) |> 
         dplyr::arrange(phat) |> 
         dplyr::mutate(rej_rate = (1:dplyr::n())/dplyr::n()),
       aes(x = phat, y = rej_rate)) +
  geom_abline(aes(slope = 1, intercept = 0)) +
  geom_point() +
  # scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  labs(title = 'Rejection Rate by p-value (qqplot)',
       subtitle = '(when Theta!=0)',
       x = 'p-value', y = 'Rejection Rate')
p10
ggsave(filename = sprintf('%s/pval_qqplot_not0.pdf', save_folder), height = 4, width = 4)




grob <- gridExtra::arrangeGrob(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10,
                    # ncol = 2, 
                    layout_matrix = matrix(c(1, 1, 2, 2,
                                              3, 4, 5, 6,
                                              7, 8, 9, 10), byrow = TRUE, nrow = 3))
ggsave(sprintf('%s/all.pdf', save_folder), grob, width = 16, height = 12)



}




# check speed of i^2 vs i**2 vs i*i
# basically same (i**2 slightly faster than i^2) except i*i really bad for large i 
# t0 = Sys.time()
# for(i in 1:99999) {
#   i^2
# }
# print(Sys.time() - t0)
# 
# t0 = Sys.time()
# for(i in 1:99999) {
#   i**2
# }
# print(Sys.time() - t0)
# 
# 
# t0 = Sys.time()
# for(i in 1:99999) {
#   i*i
# }
# print(Sys.time() - t0)
# 
# temp = c(1, 2, 3, 4)
# t0 = Sys.time()
# for(i in 1:100000) {
#   temp * temp
# }
# Sys.time() - t0
# 
# 
# temp = c(1, 2, 3, 4)
# t0 = Sys.time()
# for(i in 1:100000) {
#   temp^2
# }
# Sys.time() - t0
# 
# 




