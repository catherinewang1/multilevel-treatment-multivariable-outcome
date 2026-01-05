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









#' estimate the prior G's parameters 
est_G_params <- function(df) {
  
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
#' @param t (numeric) observed T_j statistic value
#' @param B (integer) number of repetitions 
#'         (TODO: could make it adaptive (keep going until change < threshold))
#' @param sj (numeric) sj value (sd of Y|theta)
#' @param thetaj0 (numeric) H_0 theta value testing against (typically 0)
#' @param mu (numeric) estimated mu parameter
#' @param tau (numeric) estimated tau parameter
#' @param pi0 (numeric) estimated pi0 parameter
est_pval <- function(t, B, sj, thetaj0, mu, tau, pi0) {
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
  phat = mean(tstat_null > t)
  
  # TODO: could make this adaptive to B ... not for now, just choose a large enough B
  return(phat)
}




set.seed(12345)
N = 5000
G=1
s=1
mu=1
tau=2
pi0=.3

thetaj0 = 0
B = 10000 # #samples used for estimating p-value

sim_df = sim_values(N=N, G=G, s=s, mu=mu, tau=tau, pi0=pi0)
# (sim_df$theta == 0) |> mean()

# if same \theta_j0 for all j, we can just resuse the estimated tstat null distn
tstat_null = sample_tstat_null(B       = B, 
                               sj      = sim_df[j, 's'],
                               thetaj0 = thetaj0,
                               mu      = mu, # test, use true params
                               tau     = tau,
                               pi0     = pi0)


Tstat = rep(NA, nrow(sim_df))
phat  = rep(NA, nrow(sim_df))
for(j in 1:nrow(sim_df)) {
  Tstat[j] = calc_T_stat_one(Yj = sim_df[j, 'Y'],
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
  scale_x_continuous(limits = c(-6, 7)) +
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
       x ='Theta', y = 'Y (log scale)')
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











