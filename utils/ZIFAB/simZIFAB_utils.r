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
#' @param EM_iterations (integer) number of iterations for EM '
#'    (TODO: later make it adaptive,
#'           args --> max_steps, 
#'                    threshold (keep iterating until all params |est - prev est| < threshold))
#' @param NR_iterations (integer) number of iterations for Neton-Raphson (updating tau during M-step)
#' @param initial_mu (numeric)  initial value for param mu
#' @param initial_tau (numeric) initial value for param tau
#' @param initial_pi0 (numeric) initial value for param pi0
est_G_params <- function(Y, S, 
                         EM_iterations,
                         NR_iterations,
                         initial_mu = 0,
                         initial_tau = 1,
                         initial_pi0 = .5) {
  # # test
  # set.seed(12345)
  # # sim_df = sim_values(N=500, G=1, s=1, mu=1, tau=2, pi0=.3)
  # sim_df = sim_values(N=1000, G=1, s=5, mu=3, tau=1, pi0=.3)
  # Y = sim_df$Y
  # S = sim_df$s
  # EM_iterations = 100
  # initial_mu = -10
  # initial_tau = 20
  # initial_pi0 = .6
  
  
  
  assertthat::assert_that(length(Y) == length(S), msg = 'Y and S args must be vectors of the same length')
  assertthat::assert_that((EM_iterations %% 1 == 0) && EM_iterations > 0, msg = 'K arg must be positive integer')
  
  params = data.frame(step = rep(NA, EM_iterations),
                      mu   = rep(NA, EM_iterations),
                      tau  = rep(NA, EM_iterations),
                      pi0  = rep(NA, EM_iterations))
  # initial params
  # cur_mu = 0 # does setting cur_mu = 0 make it indistinguishable? I think its ok bc of the sd diff, 
  # cur_tau = 1
  # cur_pi0 = .5 # .5
  # params[1, ] = c(1, cur_mu, cur_tau, cur_pi0)
  
  
  params[1, ] = c(1, initial_mu, initial_tau, initial_pi0)
  
  for(k in 2:EM_iterations) {
    
    # E step
    P =  mapply(FUN = est_G_params_Estep,
                y=Y, s=S,
                MoreArgs = list(mu =params[k-1,  'mu'], 
                                tau=params[k-1, 'tau'], 
                                pi0=params[k-1, 'pi0']))
    # hist(P)
    
    
    # M step
    M_step_result = est_G_params_Mstep(Y=Y, S=S, P=P, NR_iterations = NR_iterations)
    
    
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
#' @param NR_iterations (integer) of number of times to iterate for fixed point
est_G_params_Mstep <- function(Y, S, P, NR_iterations=NA) {
  normalization = sum(1-P)
  w = (1-P)/normalization
  # print(normalization)
  
  # update pi0
  pi0 = mean(P)
  # update mu
  # mu = sum((1 - P)*Y) / normalization
  mu = sum(w * Y)
  
  # update tau2
  # # tau2 = (sum((1-P) * ((Y-mu)^2 - S*S))/normalization)
  # tau2 = sum(w * ((Y - mu)^2 - S*S))
  # if((!is.na(iterations)) & iterations > 0) {
  #   # update tau2 using fixed point iteration method 
  #   for(iter in 1:iterations) {
  #     tau2 = (sum(w * (tau2 + S)^(-2)))     / 
  #       (sum(w * (tau2 + S)^(-2) * ((Y - mu)^2 - S*S)))
  #     
  #     if(tau2 < 0) { # sometimes tau^2 ends up negative, put arb threshold to ensure some values
  #       tau = 1/length(Y) #  min(1/n, sqrt(|tau2|))
  #     }
  #   }
  #   
  # }
  
  tau2 = est_G_params_Mstep_tau2(mu=mu, Y=Y, S=S, P=P, NR_iterations=NR_iterations)
  if(tau2 < 0) { # sometimes tau^2 ends up negative, put arb threshold to ensure some values
    tau = 1/length(Y) #  min(1/n, sqrt(|tau2|))
  } else {
    tau = sqrt(tau2)
  }
  
  
  
  
  
  return(c(mu, tau, pi0))
}





#' Estimating G's M step- updating tau2
#' @param mu (numeric) current estimate of mu
#' @param Y (vector) of Y values                 (same length as P & S)
#' @param S (vector) of Y's standard deviations  (same length as Y & P)
#' @param P (vector) of probabilities P(Z=0|Y=y) (same length as Y & S)
#' @param NR_iterations (integer) of number of times to iterate for fixed point
est_G_params_Mstep_tau2 <- function(mu, Y, S, P, NR_iterations) {
  # could speed up by changing P inputs to normalized w ...
  
  
  # initial guess
  normalization = sum(1-P)
  w = (1-P)/normalization
  tau2_0 = sum(w * ((Y - mu)^2 - S*S))
  
  # Newton-Raphson
  my_dQtau2  = create_dQtau2( mu=mu, Y=Y, S=S, P=P)
  my_d2Qtau2 = create_d2Qtau2(mu=mu, Y=Y, S=S, P=P)
  NR_df = newtonraphson(x = tau2_0, f = my_dQtau2, fprime = my_d2Qtau2, NR_iterations = NR_iterations)
  
  return(NR_df[nrow(NR_df), 'x'])
}



#' Create fn that evaluates the first derivative of Q w.r.t. tau^2 
#' (find zero of this)
#' @param mu (numeric) the mean of G_1 
#' @param Y (vector) of Y values                 (same length as P & S)
#' @param S (vector) of Y's standard deviations  (same length as Y & P)
#' @param P (vector) of probabilities P(Z=0|Y=y) (same length as Y & S)
create_dQtau2 <- function(mu, Y, S, P) {
  w = (1-P) / sum((1-P)) # normalized eights 1-p_i
  dQtau2 <- function(tau2) {
    first   = tau2 * (sum(w * (tau2 + S)^(-2)))     
    second  = (sum(w * (tau2 + S)^(-2) * ((Y - mu)^2 - S*S)))
    return(first - second)
  }
  return(dQtau2)
}


#' Create fn that evaluates the second derivative of Q w.r.t. tau^2 
#' f(x) := first deriv (ant to find root of f) (f(x) = 0)
#' f'(x):= second deriv (use to numerically find root of f) 
#' 
#' x_{t+1} = x_t - (f(x_t)) / (f'(x_t))
#' 
#' @param mu (numeric) the mean of G_1 
#' @param Y (vector) of Y values                 (same length as P & S)
#' @param S (vector) of Y's standard deviations  (same length as Y & P)
#' @param P (vector) of probabilities P(Z=0|Y=y) (same length as Y & S)
create_d2Qtau2 <- function(mu, Y, S, P) {
  w = (1-P) / sum((1-P)) # normalized eights 1-p_i
  d2Qtau2 <- function(tau2) {
    first   = - 1 * (sum(w * (tau2 + S)^(-3) * (tau2 - S)))     
    second  = 2 * (sum(w * (tau2 + S)^(-3) * ((Y - mu)^2 - S*S)))
    return(first + second)
  }
  return(d2Qtau2)
}

#' Newton-Raphson (used for estimating tau2)
#' 
#' x_t+1 = x_t - f(x_t)/fprime(x_t)
#' @param x0 (numeric) initial starting point 
#'   (e.g. tau2_0 = sum(w * ((Y - mu)^2 - S*S)))
#' @param f (function) that takes arg x, and returns f(x)
#'   (e.g. a dQtau2)
#' @param fprime (function) that takes arg x, and returns f'(x)
#'   (e.g. a d2Qtau2)
#' @param NR_iterations (integer) number of iterations (can change to threshold)
newtonraphson <- function(x0, f, fprime, NR_iterations) {
  x = x0
  
  if(is.na(NR_iterations) | NR_iterations == 0) {
    return(data.frame(iter   = 0, 
                      x      = x,
                      fx     = NA,
                      fprimex= NA))
  }
  
  df = NULL
  for(iter in 1:NR_iterations) {
    fx = f(x)
    fprimex = fprime(x)
    df = rbind(df,
               data.frame(iter   = iter, 
                          x      = x,
                          fx     = fx,
                          fprimex= fprimex))
    x = x - fx / fprimex
  }
  return(df)
  
}

#' Calculate the test statistic:
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
calc_T_stat_G <- function(Yj, sj, thetaj0, mu, tau, pi0) {
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






#' Calculate the test statistic (when G=1, 1 obs per group) with ONLY G_1 on numerator:
#' T_j = \int L(\theta_j) dP_G_1(\theta_j)
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
calc_T_stat_G1 <- function(Yj, sj, thetaj0, mu, tau, pi0) {
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
  # num1 = pi0 * (1/sj) * exp(-(Yj^2)/(2*sj^2))
  
  num2a = -(tau^2 * Yj^2 + sj^2 * mu^2 - (tau^2 * Yj + sj^2 * mu)^2 / (tau^2 + sj^2)) / (2 * sj^2 * tau^2)
  num2b = (1/sqrt(tau^2 + sj^2))
  # num2 = (1-pi0) * exp(num2a) * num2b
  num2 = exp(num2a) * num2b
  
  numerator = num2
  
  # T statistic
  return(numerator / denominator)
}





#' Calculate the test statistic (when G=1, 1 obs per group):
#' T_j = log[\int L(\theta_j) dP_G_1(\theta_j)  ]
#'          [-------------------------------    ]
#'          []      L(\theta_j0)                ]
#' (likelihood for G>1 will look slightly diff but essentially the same 
#'  e.g. prob use average and scale the s value)
#' @param Yj (numeric) Yj value (observed value/estimated theta)
#' @param sj (numeric) sj value (sd of Y|theta)
#' @param thetaj0 (numeric) H_0 theta value testing against (typically 0)
#' @param mu (numeric) estimated mu parameter
#' @param tau (numeric) estimated tau parameter
#' @param pi0 (numeric) estimated pi0 parameter
calc_T_stat_logG1 <- function(Yj, sj, thetaj0, mu, tau, pi0) {
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
  
  # T statistic
  return( (Yj/sj + sj*mu/(tau^2))^2)
}

#' Calculate the T statistic as the posterior mean
#' T_j = E(\theta_j | Y_j) = ...
#' 
#' 
#' 
#' @param Yj (numeric) Yj value (observed value/estimated theta)
#' @param sj (numeric) sj value (sd of Y|theta)
#' @param thetaj0 (numeric) H_0 theta value testing against (typically 0)
#' @param mu (numeric) estimated mu parameter
#' @param tau (numeric) estimated tau parameter
#' @param pi0 (numeric) estimated pi0 parameter
calc_T_stat_posteriormean <- function(Yj, sj, thetaj0, mu, tau, pi0) {
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
  
  # TODO: calculate this val
  return(NA)
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
#' @param T_stat_name (character) name of the T statistic: G, G1, posteriormean
sample_tstat_null <- function(B, sj, thetaj0, mu, tau, pi0, T_stat_name) {
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
  if(T_stat_name == 'G') {
    tstat_null = sapply(X = Y, FUN = calc_T_stat_G,
                        sj=sj, thetaj0=thetaj0, mu=mu, tau=tau, pi0=pi0)
  } else if(T_stat_name == 'G1') {
    tstat_null = sapply(X = Y, FUN = calc_T_stat_G1,
                        sj=sj, thetaj0=thetaj0, mu=mu, tau=tau, pi0=pi0)
  } else if(T_stat_name == 'logG1') {
    tstat_null = sapply(X = Y, FUN = calc_T_stat_logG1,
                        sj=sj, thetaj0=thetaj0, mu=mu, tau=tau, pi0=pi0)
  } else if(T_stat_name == 'posteriormean') {
    tstat_null = sapply(X = Y, FUN = calc_T_stat_posteriormean,
                        sj=sj, thetaj0=thetaj0, mu=mu, tau=tau, pi0=pi0)
  } else {
    return('bad T_stat_name argument (one of G, G1, posteriormean)')
  }
  
  
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
#' @param T_stat_name (character) name of the T statistic: G, G1, logG1, posteriormean
calc_pval <- function(tj, B, sj, thetaj0, mu, tau, pi0, T_stat_name) {
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
  # t = calc_T_stat_G(Yj=Yj, sj = sj, thetaj0 = thetaj0,
  #                     # test, use true params
  #                     mu = mu, tau = tau, pi0 = pi0)
  # B = 10000
  
  
  # do in another fn to reuse values
  # # draw Ys
  # Y = rnorm(n = B, mean = thetaj0, sd = sj)
  # # calculate tstats when H0:thetaj=thetaj0
  # tstat_null = sapply(X = Y, FUN = calc_T_stat_G,
  #                     sj=sj, thetaj0=thetaj0, mu=mu, tau=tau, pi0=pi0)
  
  
  if(T_stat_name == 'G' | T_stat_name == 'G1') {
    # use simulated nulls
    
    # TODO: could make this adaptive to B ... not for now, just choose a large enough B
    tstat_null = sample_tstat_null(B, sj, thetaj0, mu, tau, pi0, T_stat_name)
    
    # estimate p-value (empirical distribution)
    pval = mean(tstat_null > tj)
    
  } else if(T_stat_name == 'logG1') {
    # something is not right, sample for no to check
    
    # use exact theoretical distribution
    # under null, T stat is a non-central chi-sq(1df, sj mu /tau^2)
    pval = 1 - pchisq(tj, df = 1, ncp = (sj*mu/(tau^2))^2 )
    
    
    
    tstat_null = sample_tstat_null(B, sj, thetaj0, mu, tau, pi0, T_stat_name)
    pval = mean(tstat_null > tj)
    
    
  } else if(T_stat_name == 'posteriormean') {
    # use exact theoretical distribution
    # under null, T stat is a normal?
    return('calc pval for posteriormean not impl yet')
    
    
  } else {
    return('bad T_stat_name argument (G, G1, logG1, posteriormean)')
  }
  
  
  
  return(pval)
}


#' check EM algorithm converges as sample size increases
#' 
#' @param save_folder (character)
#' @param G (integer)
#' @param s (numeric or numeric  vector)
#' @param mu (numeric)
#' @param tau (numeric)
#' @param pi0 (numeric)
#' @param EM_iterations (integer) number of EM iteration steps
#' @param samplesize_for_numsteps (integer) the sample size used for the run using EM_iterations
#' @param samplesizes (vector) of integers, sample sizes to test out
#' @param nreps (integer) number of times to repeat
check_EM <- function(save_folder,
                     G, s, mu, tau, pi0,
                     EM_iterations, 
                     NR_iterations,
                     samplesize_for_numsteps,
                     samplesizes,
                     nreps) {
  # # debug
  # set.seed(123456)
  # G=1
  # s=5
  # mu=3
  # tau=1
  # pi0=.3
  # EM_iterations=250
  # samplesize_for_numsteps = 1000
  # samplesizes = c(25, 50, 100, 250, 500, 1000, 1500, 2000, 2500, 3000, 4000, 5000)
  # samplesizes = c(25, 50, 100, 250, 500)
  # nreps = 3
  
  
  
  EM_est = NULL; estparams = NULL
  for(rep in 1:nreps) {
    
    # EM over different sample sizes
    sim_df = sim_values(N=max(samplesizes, samplesize_for_numsteps), G=G, s=s, mu=mu, tau=tau, pi0=pi0)
    for(n in samplesizes) {
      t0 = Sys.time()
      EMparams = est_G_params(Y = sim_df[1:n, 'Y'], 
                              S = sim_df[1:n, 's'], 
                              EM_iterations=EM_iterations,
                              NR_iterations=NR_iterations,
                              initial_mu = 0,
                              initial_tau = 1,
                              initial_pi0 = .5)
      t1 = Sys.time()
      EM_est = rbind(EM_est, 
                     data.frame(rep = rep,
                                n   = n,
                                mu  = EMparams[nrow(EMparams),  'mu'],
                                tau = EMparams[nrow(EMparams), 'tau'],
                                pi0 = EMparams[nrow(EMparams), 'pi0'], 
                                time = difftime(t1, t0, units = 'sec')))
    }
    
    # EM over number of steps (using samplesize_for_numsteps)
    estparams = rbind(estparams,
                      est_G_params(Y = sim_df[1:samplesize_for_numsteps, 'Y'], 
                                   S = sim_df[1:samplesize_for_numsteps, 's'], 
                                   EM_iterations=EM_iterations,
                                   NR_iterations=NR_iterations,
                                   initial_mu = 0,
                                   initial_tau = 1,
                                   initial_pi0 = .5) |> dplyr::mutate(rep = rep))
  }
  
  
  
  # display EM parameter estimates over iterations
  p_mu_samplesize  = ggplot() + geom_line(data=EM_est, aes(x = n, y =  mu, group = rep)) + labs(title = 'mu')  + geom_hline(aes(yintercept =  mu), color = 'orange')
  p_tau_samplesize = ggplot() + geom_line(data=EM_est, aes(x = n, y = tau, group = rep)) + labs(title = 'tau') + geom_hline(aes(yintercept = tau), color = 'orange')
  p_pi0_samplesize = ggplot() + geom_line(data=EM_est, aes(x = n, y = pi0, group = rep)) + labs(title = 'pi0') + geom_hline(aes(yintercept = pi0), color = 'orange')
  # grob <- gridExtra::grid.arrange(p_mu_samplesize, p_tau_samplesize, p_pi0_samplesize)
  # ggsave(sprintf('%s/check_EM_bysamplesize.pdf', save_folder), grob, width = 4, height = 6)
  
  
  
  # Visualizations
  
  # display EM parameter estimates over iterations
  p_mu  = ggplot() + geom_line(data=estparams, aes(x = step, y =  mu, group = rep)) + labs(title = 'mu')  + geom_hline(aes(yintercept =  mu), color = 'orange')
  p_tau = ggplot() + geom_line(data=estparams, aes(x = step, y = tau, group = rep)) + labs(title = 'tau') + geom_hline(aes(yintercept = tau), color = 'orange')
  p_pi0 = ggplot() + geom_line(data=estparams, aes(x = step, y = pi0, group = rep)) + labs(title = 'pi0') + geom_hline(aes(yintercept = pi0), color = 'orange')
  # grob <- gridExtra::grid.arrange(p_mu, p_tau, p_pi0)
  # ggsave(sprintf('%s/estparam_EM_byiterations.pdf', save_folder), grob, width = 4, height = 6)
  
  
  
  
  p_EM_label_samplesize = ggplot() + geom_text(aes(x=0,y=0, label = 'EM Estimates \n by sample size'       ), size = 8) + theme_nothing()
  p_EM_label_iteration  = ggplot() + geom_text(aes(x=0,y=0, label = 'EM Estimates \n by EM iteration steps'), size = 8) + theme_nothing()
  
  
  
  allplots = list(p_mu_samplesize  = p_mu_samplesize, 
                  p_tau_samplesize = p_tau_samplesize, 
                  p_pi0_samplesize = p_pi0_samplesize, 
                  p_mu             = p_mu, 
                  p_tau            = p_tau, 
                  p_pi0            = p_pi0, 
                  p_EM_label_samplesize = p_EM_label_samplesize, 
                  p_EM_label_iteration  = p_EM_label_iteration)
  
  
  
  
  grob <- gridExtra::arrangeGrob(p_mu_samplesize, p_tau_samplesize, p_pi0_samplesize, # 1, 2, 3
                                 p_mu, p_tau, p_pi0,   # 4, 5, 6
                                 p_EM_label_samplesize, p_EM_label_iteration, # 7, 8
                                 # ncol = 2, 
                                 layout_matrix = matrix(c(7, 1, 2, 3,
                                                          8, 4, 5, 6), byrow = TRUE, nrow = 2))
  # gridExtra::grid.arrange(grob)
  ggsave(sprintf('%s/checkEM.pdf', save_folder), grob, width = 14, height = 6)
  
  return(EM_est)
  
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
#' @param EM_iterations (integer) number of EM iteration steps
sim_and_plot <- function(save_folder,
                         N, G, s, mu, tau, pi0,
                         thetaj0, B,
                         EM_iterations=100, 
                         NR_iterations=20) {
  # # debug
  # N=1000
  # G=1
  # s=5
  # mu=3
  # tau=1
  # pi0=.3
  # thetaj0=0
  # B=1000
  # EM_iterations=1000
  
  
  dir.create(save_folder, showWarnings = TRUE)
  
  sim_df = sim_values(N=N, G=G, s=s, mu=mu, tau=tau, pi0=pi0)
  # (sim_df$theta == 0) |> mean()
  
  
  
  
  
  
  
  # estimate the prior, G, parameters
  estparams = est_G_params(Y = sim_df[ , 'Y'], 
                           S = sim_df[ , 's'], 
                           EM_iterations=EM_iterations,
                           NR_iterations=NR_iterations,
                           initial_mu = 0,
                           initial_tau = 1,
                           initial_pi0 = .5)
  
  estmu  = estparams[nrow(estparams),  'mu']
  esttau = estparams[nrow(estparams), 'tau']
  estpi0 = estparams[nrow(estparams), 'pi0']
  
  # setup list of diff options
  paramList = list('true' = list(thetaj0 = thetaj0,
                                 mu      = mu, 
                                 tau     = tau,
                                 pi0     = pi0),
                   'EM'   = list(thetaj0 = thetaj0,
                                 mu      = estmu, 
                                 tau     = esttau,
                                 pi0     = estpi0))
  
  
  tstatList = list('G'     = calc_T_stat_G, 
                   'G1'    = calc_T_stat_G1,
                   'logG1' = calc_T_stat_logG1) 
  
  
  # start calc t-stats + pvals
  
  # Standard t-statistics + p-values (compare Y to N(0, s^2))
  Tstat_standard =              sim_df[ ,'Y'] / sim_df[ ,'s']
  p_standard     = 2*pnorm(-abs(sim_df[ ,'Y'] / sim_df[ ,'s']))
  
  simres_df = cbind(sim_df,
                    data.frame(method = 'standard',
                               tstat  = 'standard', 
                               param  = 'NA',
                               Tstat  = Tstat_standard,
                               pval   = p_standard))
  # other methods
  for(T_stat_name in c('G', 'G1', 'logG1')) { # for each t stat type
    for(paramtype in c('true', 'EM')) {     # using either true or EM estimated parameters
      
      cur_Tstats =  mapply(FUN = tstatList[[T_stat_name]],
                           Yj = sim_df[ ,'Y'],
                           sj = sim_df[ , 's'], 
                           MoreArgs = paramList[[paramtype]])
      
      cur_pvals  = mapply(FUN = calc_pval,
                          tj = cur_Tstats,
                          sj = sim_df[ , 's'], 
                          MoreArgs = c(list(B          = B,
                                            T_stat_name=T_stat_name),
                                       paramList[[paramtype]]))
      
      simres_df = rbind(simres_df,
                        cbind(sim_df,
                              data.frame(method = sprintf('ZIFAB_%sG_T%s', paramtype, T_stat_name),
                                         tstat  = T_stat_name, 
                                         param  = paramtype,
                                         Tstat  = cur_Tstats,
                                         pval   = cur_pvals)))
    }
  }
  
  
  
  simres_df$method = factor(simres_df$method,
                            levels = c('standard',     'ZIFAB_trueG_TG',    'ZIFAB_EMG_TG',    'ZIFAB_trueG_TG1',      'ZIFAB_EMG_TG1',     'ZIFAB_trueG_TlogG1',     'ZIFAB_EMG_TlogG1'),
                            labels = c('standard', 'ZIFAB (true G) [G]','ZIFAB (EM G) [G]','ZIFAB (true G) [G1]',  'ZIFAB (EM G) [G1]', 'ZIFAB (true G) [logG1]', 'ZIFAB (EM G) [logG1]'))
  
  # # True parameters
  # # calculate the T statistic using true parameters
  # Tstat_trueparams = mapply(FUN = calc_T_stat_G,
  #                           Yj = sim_df[ ,'Y'],
  #                           sj = sim_df[ , 's'], 
  #                           MoreArgs = list(thetaj0 = thetaj0,
  #                                           mu      = mu, 
  #                                           tau     = tau,
  #                                           pi0     = pi0))
  # 
  # # estimate the p-values (using simulation) using true parameters
  # phat_trueparams = mapply(FUN = calc_pval,
  #                          tj = Tstat_trueparams,
  #                          sj = sim_df[ , 's'], 
  #                          MoreArgs = list(B       = B,
  #                                          thetaj0 = thetaj0,
  #                                          mu      = mu, 
  #                                          tau     = tau,
  #                                          pi0     = pi0,
  #                                          T_stat_name='G'))
  # 
  # # calculate the T statistic w only G1 in the numerator using true parameters
  # TstatG1_trueparams = mapply(FUN = calc_T_stat_G1,
  #                             Yj = sim_df[ ,'Y'],
  #                             sj = sim_df[ , 's'], 
  #                             MoreArgs = list(thetaj0 = thetaj0,
  #                                             mu      = mu, 
  #                                             tau     = tau,
  #                                             pi0     = pi0))
  # 
  # # estimate the p-values w only G1 in the numerator (using simulation) using true parameters
  # phatG1_trueparams = mapply(FUN = calc_pval,
  #                            tj = TstatG1_trueparams,
  #                            sj = sim_df[ , 's'], 
  #                            MoreArgs = list(B       = B,
  #                                            thetaj0 = thetaj0,
  #                                            mu      = mu, 
  #                                            tau     = tau,
  #                                            pi0     = pi0,
  #                                            T_stat_name='G1'))
  #  
  # # Estimated (using EM alg) parameters
  # # estimate the prior, G, parameters
  # estparams = est_G_params(Y = sim_df[ , 'Y'], 
  #                          S = sim_df[ , 's'], 
  #                          EM_iterations=EM_iterations,
  #                          initial_mu = 0,
  #                          initial_tau = 1,
  #                          initial_pi0 = .5)
  # 
  # estmu  = estparams[nrow(estparams),  'mu']
  # esttau = estparams[nrow(estparams), 'tau']
  # estpi0 = estparams[nrow(estparams), 'pi0']
  # 
  # 
  # # calculate the T statistic using true parameters
  # Tstat_estparams = mapply(FUN = calc_T_stat_G,
  #                          Yj = sim_df[ ,'Y'],
  #                          sj = sim_df[ , 's'], 
  #                          MoreArgs = list(thetaj0 = thetaj0,
  #                                          mu      = estmu, 
  #                                          tau     = esttau,
  #                                          pi0     = estpi0))
  # 
  # # estimate the p-values (using simulation) using true parameters
  # phat_estparams = mapply(FUN = calc_pval,
  #                         tj = Tstat_estparams,
  #                         sj = sim_df[ , 's'], 
  #                         MoreArgs = list(B       = B,
  #                                         thetaj0 = thetaj0,
  #                                         mu      = estmu, 
  #                                         tau     = esttau,
  #                                         pi0     = estpi0,
  #                                         T_stat_name = 'G'))
  # 
  # # calculate the T statistic w only G1 in the numerator using est parameters
  # TstatG1_estparams = mapply(FUN = calc_T_stat_G1,
  #                             Yj = sim_df[ ,'Y'],
  #                             sj = sim_df[ , 's'], 
  #                             MoreArgs = list(thetaj0 = thetaj0,
  #                                             mu      = estmu, 
  #                                             tau     = esttau,
  #                                             pi0     = estpi0))
  # 
  # # estimate the p-values w only G1 in the numerator (using simulation) using est parameters
  # phatG1_estparams = mapply(FUN = calc_pval,
  #                            tj = TstatG1_estparams,
  #                            sj = sim_df[ , 's'], 
  #                            MoreArgs = list(B       = B,
  #                                            thetaj0 = thetaj0,
  #                                            mu      = estmu, 
  #                                            tau     = esttau,
  #                                            pi0     = estpi0,
  #                                            T_stat_name='G1'))
  # 
  # 
  # 
  # # combine to dataframe
  # simres_df = rbind(cbind(sim_df,
  #                         data.frame(method = 'standard',
  #                                    Tstat  = Tstat_standard,
  #                                    pval   = p_standard)),
  #                   # ZIFAB true G
  #                   cbind(sim_df,
  #                         data.frame(method = 'ZIFAB_trueG_TG',
  #                                    Tstat  = Tstat_trueparams,
  #                                    pval   = phat_trueparams)),
  #                   cbind(sim_df,
  #                         data.frame(method = 'ZIFAB_trueG_TG1',
  #                                    Tstat  = TstatG1_trueparams,
  #                                    pval   = phatG1_trueparams)),
  #                   # ZIFAB estimated G
  #                   cbind(sim_df,
  #                         data.frame(method = 'ZIFAB_EMG_TG',
  #                                    Tstat  = Tstat_estparams,
  #                                    pval   = phat_estparams)),
  #                   cbind(sim_df,
  #                         data.frame(method = 'ZIFAB_EMG_TG1',
  #                                    Tstat  = TstatG1_estparams,
  #                                    pval   = phatG1_estparams))
  # )
  # 
  # 
  # simres_df$method = factor(simres_df$method,
  #                           levels = c('standard',     'ZIFAB_trueG_TG',    'ZIFAB_EMG_TG',    'ZIFAB_trueG_TG1',      'ZIFAB_EMG_TG1'),
  #                           labels = c('standard', 'ZIFAB (true G) [G]','ZIFAB (EM G) [G]','ZIFAB (true G) [G1]',  'ZIFAB (EM G) [G1]'))
  
  
  
  # method_colors = c('blue', 'orange', 'brown')
  set2_colors = RColorBrewer::brewer.pal(7, name = 'Set2')
  
  barplot( 1:length(set2_colors), col = set2_colors)
  
  
  method_colors = c(set2_colors[2], 
                    colorRampPalette(c('white', set2_colors[1]))(9)[c(4, 9)], 
                    colorRampPalette(c('white', set2_colors[3]))(9)[c(4, 9)],
                    colorRampPalette(c('white', set2_colors[4]))(9)[c(4, 9)])
  barplot( 1:length(method_colors), col = method_colors)
  
  # test EM estimation performance (same parameters, change sample size)
  
  # Estimated (using EM alg) parameters
  sim_df = sim_values(N=5000, G=G, s=sample(s, size = 5000, replace = TRUE), mu=mu, tau=tau, pi0=pi0)
  Ns = c(25, 50, 100, 250, 500, 1000, 1500, 2000, 2500, 3000, 4000, 5000)
  EM_est = NULL
  for(n in Ns) {
    t0 = Sys.time()
    EMparams = est_G_params(Y = sim_df[1:n, 'Y'], 
                            S = sim_df[1:n, 's'], 
                            EM_iterations=EM_iterations,
                            NR_iterations=NR_iterations,
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
  
  
  
  
  p_EM_label_samplesize = ggplot() + geom_text(aes(x=0,y=0, label = 'EM Estimates \n by sample size'       ), size = 8) + theme_nothing()
  p_EM_label_iteration  = ggplot() + geom_text(aes(x=0,y=0, label = 'EM Estimates \n by EM iteration steps'), size = 8) + theme_nothing()
  
  
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
                                 p_mu_samplesize, p_tau_samplesize, p_pi0_samplesize, # 14, 15, 16
                                 p_EM_label_samplesize, p_EM_label_iteration, # 17, 18
                                 # ncol = 2, 
                                 layout_matrix = matrix(c(17, 14, 15, 16,
                                                          18, 11, 12, 13,
                                                          1, 2, NA, NA,
                                                          3, 4, 5, 6,
                                                          7, 8, 9, 10), byrow = TRUE, nrow = 5))
  ggsave(sprintf('%s/all.pdf', save_folder), grob, width = 24, height = 18)
  
  
  
  
  
  
  
  
  
}

