# Run simulation for Zero-Inflated FAB p-values
# functions are defined in utils/ZIFAB/simZIFAB_utils.r
# 
# main functions to call are:
#  - sim_and_plot
#  - check_EM



library(ggplot2)
library(cowplot)

ggplot2::set_theme(theme_cowplot() +
                     theme(title = element_text(hjust = .5)))


source('../../utils/ZIFAB/simZIFAB_utils.r')



# save_folder_overall = "C:/Users/Cathe/Documents/School/Genetic Hypothesis Testing using Negative Controls/multilevel-treatment-multivariable-outcome-github/multilevel-treatment-multivariable-outcome/plots/ZIFAB/sim/"
save_folder_general = '../../plots/'
save_folder_overall = sprintf('%s/ZIFAB/sim/', save_folder_general)
assertthat::assert_that(dir.exists(save_folder_general), msg = 'save_folder_general must exist')
dir.create(save_folder_overall, recursive = TRUE)

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
             B=1000,
             nsplits = 10)


# simulation with  well separated groups 
set.seed(12345)
sim_and_plot(save_folder = paste0(save_folder_overall, "B/"), 
             N=1000, 
             G=1,
             s=.5,
             mu=3,
             tau=.5,
             pi0=.3,
             thetaj0=0,
             B=10000,
             nsplits = 10)

# simulation with varying si 
set.seed(12345)
sim_and_plot(save_folder = paste0(save_folder_overall, "C/"), 
             N=1000, 
             G=1,
             s=1.8*rbeta(1000, 2, 1) + .2, # range from .2 to 2
             mu=1,
             tau=2,
             pi0=.3,
             thetaj0=0,
             B=10000,
             nsplits = 10)


set.seed(12345)
EM_est_C = 
  check_EM(save_folder = paste0(save_folder_overall, "C/"), 
           G=1,
           s=1.8*rbeta(1500, 2, 1) + .2, # range from .2 to 2
           mu=1,
           tau=2,
           pi0=.3,
           EM_iterations=100, 
           NR_iterations=5,
           samplesize_for_numsteps = 1000,
           samplesizes = c(25, 50, 100,  250, 500, 1000, 1500),
           nreps = 6)




# simulation with larger si (s.t. standard has low power)
set.seed(12345)
sim_and_plot(save_folder = paste0(save_folder_overall, "D/"), 
             N=1000, 
             G=1,
             s=5,
             mu=3,
             tau=1,
             pi0=.3,
             thetaj0=0,
             B=10000,
             EM_iterations=1000, 
             NR_iterations=20,
             nsplits = 10)

set.seed(12345)
EM_est_D = 
  check_EM(save_folder = paste0(save_folder_overall, "D/"), 
         G=1,
         s=5,
         mu=3,
         tau=1,
         pi0=.3,
         EM_iterations=1000, 
         NR_iterations=20,
         samplesize_for_numsteps = 1000,
         samplesizes = c(25, 50, 100, 250, 500, 1000, 1500, 2000, 2500, 3000, 4000),
         nreps = 15)




# # NaN's produced in EM
# set.seed(12345)
# sim_df = sim_values(N=5000, G=1, s=5, mu=3, tau=1, pi0=.3)
# 
# estparams = est_G_params(Y = sim_df[ , 'Y'], 
#                          S = sim_df[ , 's'], 
#                          num_steps=1000,
#                          initial_mu = 0,
#                          initial_tau = 1,
#                          initial_pi0 = .5)
# p_mu  = ggplot(estparams, aes(x = step, y =  mu)) + geom_line() + labs(title = 'mu')
# p_tau = ggplot(estparams, aes(x = step, y = tau)) + geom_line() + labs(title = 'tau')
# p_pi0 = ggplot(estparams, aes(x = step, y = pi0)) + geom_line() + labs(title = 'pi0')
# gridExtra::grid.arrange(p_mu, p_tau, p_pi0)
# estparams
# estparams[nrow(estparams), ]


# simulation with larger si (s.t. standard has low power), and less distinction between theta_i groups
set.seed(12345)
sim_and_plot(save_folder = paste0(save_folder_overall, "E/"), 
             N=1000, 
             G=1,
             s=5,
             mu=2,
             tau=2,
             pi0=.3,
             thetaj0=0,
             B=10000,
             num_steps=1000)

# simulation with larger si (s.t. standard has low power), and more distinction between theta_i groups
set.seed(12345)
sim_and_plot(save_folder = paste0(save_folder_overall, "F/"), 
             N=1000, 
             G=1,
             s=5,
             mu=5,
             tau=1,
             pi0=.3,
             thetaj0=0,
             B=10000,
             num_steps=1000)





if(F) {
  # check diff between standard method for p-values vs some linear transformation a(Y-b)
  
  a = 1 
  b = 1
  s = 2
  
  Y_sample =  rnorm(n = 1000, mean = 0, sd = s)
  
  
  pval_standard = 2*pnorm(-abs(Y_sample / s))
  
  
  Tstats_lineartrans = a * (Y_sample - b)
  
  pval_lineartrans  = pnorm(Tstats_lineartrans, mean = -a*b, sd = a^2 * s^2)# why are these not the same
  pval_lineartrans  = 2*pmin(pval_lineartrans, 1- pval_lineartrans)
  pval_lineartrans2 = 2*pnorm(-(abs(Tstats_lineartrans + a*b)/a*s), mean = 0, sd = 1)
  
  
  
  plot(pval_standard, pval_lineartrans)
  plot(pval_lineartrans, pval_lineartrans2)
  plot(pval_standard, pval_lineartrans2)
  
  
  
  
}


if(F) {
  # Checking the distribution of logG1, which we thought to be noncenral chisq, but is not working...
  
  
  set.seed(12345)
  N=1000
  G=1
  s=5
  mu=3
  tau=1
  pi0=.3
  thetaj0=0
  B=1000
  num_steps=1000
  
  sim_df = sim_values(N=N, G=G, s=s, mu=mu, tau=tau, pi0=pi0)
  # (sim_df$theta == 0) |> mean()
  
  j = 3
  tstat_j = calc_T_stat_logG1(Yj = sim_df[j, 'Y'], 
                              sj = s, 
                              thetaj0 = 0, 
                              mu = mu, 
                              tau = tau, 
                              pi0 = pi0)
  
  
  # under the null, the tstat should be noncentral chisq
  Y_sample = rnorm(n = 1000, mean = 0, sd = s)
  
  tstats_null = mapply(FUN = calc_T_stat_logG1,
         Yj = Y_sample,
         sj = s, 
         MoreArgs = list(thetaj0 = 0, 
                         mu = mu, 
                         tau = tau, 
                         pi0 = pi0))
  
  hist(tstats_null)
  
  x = seq(from = .1, to = 300, length.out = 100)
  y = dchisq(x, df = 1, ncp = (s*mu/(tau^2))^2)
  plot(x, y)
  
  hist(rchisq(n = 1000, df = 1, ncp =(s*mu/(tau^2))^2))
  
  (Yj/sj + sj*mu/(tau^2))^2
  
  
  hist(Y_sample / s)
  hist(Y_sample / s + s * mu / (tau^2))
  hist((Y_sample / s + s * mu / (tau^2))^2)
  
  qqplot(tstats_null, rchisq(n = 1000, df = 1, ncp =(s*mu/(tau^2))^2))
  
  
  # suppose Y sim N(0, 5)
  # c = 15
  
  Y_sample = rnorm(1000, mean = 0, sd = 5)
  hist(Y_sample + 5)
  hist(Y_sample + 5)
  
  
  
}



# 
# 
# theme_path = 'C:/Users/Cathe/Downloads/Baby, The Code Shines Bright.tmTheme'
# theme_path = 'C:/Users/Cathe/Downloads/base16-3024-modified.rstheme'
# rstudioapi::addTheme(themePath=theme_path, apply=TRUE, force = TRUE)
# 
# 
# rstudioapi::convertTheme(
#   themePath=theme_path, 
#   add=FALSE, 
#   # outputLocation, 
#   apply=TRUE, 
#   # force, globally
# )


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
Tstat_trueparams = mapply(FUN = calc_T_stat_G,
                          Yj = sim_df[ ,'Y'],
                          sj = sim_df[ , 's'], 
                          MoreArgs = list(thetaj0 = thetaj0,
                                          mu      = mu, 
                                          tau     = tau,
                                          pi0     = pi0))

# estimate the p-values (using simulation) using true parameters
phat_trueparams = mapply(FUN = calc_pval,
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
Tstat_estparams = mapply(FUN = calc_T_stat_G,
                          Yj = sim_df[ ,'Y'],
                          sj = sim_df[ , 's'], 
                          MoreArgs = list(thetaj0 = thetaj0,
                                          mu      = estmu, 
                                          tau     = esttau,
                                          pi0     = estpi0))

# estimate the p-values (using simulation) using true parameters
phat_estparams = mapply(FUN = calc_pval,
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
  Tstat_trueparams[j] = calc_T_stat_G(
                             Yj = sim_df[j, 'Y'],
                             sj = sim_df[j, 's'],
                             thetaj0 = thetaj0,
                             mu      = mu, # test, use true params
                             tau     = tau,
                             pi0     = pi0)
  # same thetaj0 for all j
  phat[j] = mean(tstat_null > Tstat[j])
  # # code for potentially diff thetaj0 for every j
  # phat[j] = calc_pval(t = Tstat[j],
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






