



# Test EM




library(dplyr)
library(ggplot2)
library(cowplot)

ggplot2::set_theme(theme_cowplot() +
                     theme(title = element_text(hjust = .5)))


source('./simZIFAB_utils.r')



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
fixedpoint_iter_steps = 100

sim_df = sim_values(N=N, G=G, s=s, mu=mu, tau=tau, pi0=pi0)




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



# E step
P =  mapply(FUN = est_G_params_Estep,
            y=sim_df$Y, s=sim_df$s,
            MoreArgs = list(mu =0, 
                            tau=1, 
                            pi0=.5))


normalization = sum(1-P)
w = (1-P)/normalization
# print(normalization)


Y= sim_df$Y; S = sim_df$s

# update pi0
pi0 = mean(P)
# update mu
# mu = sum((1 - P)*Y) / normalization
mu = sum(w * Y)


iter_times = NA
iter_times = 10
# update tau2
# tau2 = (sum((1-P) * ((Y-mu)^2 - S*S))/normalization)
tau2_0 = sum(w * ((Y - mu)^2 - S*S))
tau2 = tau2_0
tau2s = c(tau2)
nums = c()
denoms = c()

if((!is.na(iter_times)) & iter_times > 0) {
  # update tau2 using fixed point iteration method 
  for(iter in 1:iter_times) {
    
    my_dQtau2 = create_dQtau2(mu=mu, Y=Y, S=S, P=P)
    gridtau2 = seq(from = .1, to = 10, length.out = 100)
    f_x = sapply(X=gridtau2, FUN = my_dQtau2)
    
    plot(gridtau2, f_x); abline(h = 0); abline(v = tau2_0)
    my_dQtau2(tau2_0)
    
    num   = (sum(w * (tau2 + S)^(-2)))     
    denom = (sum(w * (tau2 + S)^(-2) * ((Y - mu)^2 - S*S)))
    

    
    print(paste0(tau2, ' ', num, ' ', denom))
    
    
    tau2 = num / denom
    
    # if(tau2 < 0) { # sometimes tau^2 ends up negative, put arb threshold to ensure some values
    #   tau = 1/length(Y) #  min(1/n, sqrt(|tau2|))
    # }
    tau2s = c(tau2s, tau2)
    nums   = c(nums, num)
    denoms = c(denoms, denom)
  }
  
}

plot(1:length(tau2s), tau2s)
plot(1:length(nums), nums)
plot(1:length(denoms), denoms)






# Question: For a function like f(x) = (deriv of Q rt tau2), is the proposed solution a root? 
# Ans: No, it only orked for s_i = s (equal variance of Y_i|\theta_i)
# f(tau2) = tau2 * (sum(w * (tau2 + S)^(-2)))     -   (sum(w * (tau2 + S)^(-2) * ((Y - mu)^2 - S*S)))
# tau2_0 = sum(w * ((Y - mu)^2 - S*S))
# To make the problem more general, consider any a_i = ((Y - mu)^2 - S*S)_i
# f(x) = x (sum(w * (x + S)^(-2))) - (sum(w * (x + S)^(-2) * a))
# x_0 = sum(w*a)

f_ <- function(a, S) {
  f_inner <- function(x) {
    x * (sum(w * (x + S)^(-2))) - (sum(w * (x + S)^(-2) * a))
  }
  return(f_inner)
}

# S = sim_df$s
S = sample(c(1, 3, 5), size = nrow(sim_df), replace = TRUE)
normalization = sum(1-P)
w = (1-P)/normalization
x = seq(from = -.5, to = 5, length.out = 1000) # tau2 should not be neg, but can happen numerically

# try out different sequences to see if true/approx for any a
set.seed(13456)
a_s = list(rnorm(n = length(S), mean = .5, sd = 3),
           rnorm(n = length(S), mean = 0, sd = 7),
           runif(n = length(S), min = -3, max = 3),
           runif(n = length(S), min = -8, max = 8),
           runif(n = length(S), min = 0, max = 8))


 

temp_df = NULL# see f_x over a range of x values
a_df = NULL   # see a values ((y-mu)^ - s^)
x_0_df = NULL # see ho close x_0 is a to a root
for(idx in 1:length(a_s)) {
  my_f_ = f_(a = a_s[[idx]], S = S)
  f_x = sapply(X=x, FUN=my_f_)
  x_0 = sum(w*a_s[[idx]])
  # plot(x, f_x, type = 'l'); abline(h = 0); abline(v = x_0, col='red')
  temp_df = rbind(temp_df,
                  data.frame(x   = x, 
                             f_x = f_x,
                             x_0 = x_0, 
                             group= idx))
  a_df = rbind(a_df,
               data.frame(group = idx,
                          a = a_s[[idx]]))
  x_0_df = rbind(x_0_df,
                 data.frame(group = idx,
                            x_0 = x_0,
                            f_x_0 = my_f_(x_0)))
}


x_0_df # 


library(ggplot2)
library(cowplot)
ggplot2::theme_set(theme_cowplot())



ggplot(temp_df, aes(x = x, y = f_x, group=group, linetype = factor(group))) +
  geom_line() +
  geom_hline(aes(yintercept = 0), color = 'gray') +
  geom_vline(aes(xintercept = x_0, linetype = factor(group)), color = 'orange') +
  facet_grid(vars(group))


ggplot(a_df, aes(x = a, group=group)) +
  geom_histogram(binwidth = .5) +
  facet_grid(vars(group))






# Newton-Raphson for estimating tau2

set.seed(12345)
N=1000
G=1
s= sample(c(1, 3, 5), size = N, replace = TRUE)
mu=3
tau=1
pi0=.3


sim_df = sim_values(N=N, G=G, s=s, mu=mu, tau=tau, pi0=pi0)


Y = sim_df$Y
S = sim_df$s

# E step
P =  mapply(FUN = est_G_params_Estep,
            y=Y, s=s,
            MoreArgs = list(mu =0, 
                            tau=1, 
                            pi0=.5))
# check M step
my_dQtau2  = create_dQtau2( mu=mu, Y=Y, S=S, P=P)
my_d2Qtau2 = create_d2Qtau2(mu=mu, Y=Y, S=S, P=P)
gridx = seq(from = .1, to = 5, length.out = 100)
f_x      = sapply(X=gridx, FUN = my_dQtau2)
fprime_x = sapply(X=gridx, FUN = my_d2Qtau2)


plot(gridx, f_x)
plot(gridx, fprime_x)

# try NR for tau
NR_df = newtonraphson(x = 4, f = my_dQtau2, fprime = my_d2Qtau2, iterations = 40)





# try all together

est_G_params_Mstep(Y, S, P, iterations=30)

test_EM = est_G_params(Y, S, 
                         EM_iterations = 10,
                         NR_iterations = 1,
                         initial_mu = 0,
                         initial_tau = 1,
                         initial_pi0 = .5)

test_EM = est_G_params(Y, S, 
             EM_iterations = 100,
             NR_iterations = 30,
             initial_mu = 0,
             initial_tau = 1,
             initial_pi0 = .5)
tail(test_EM)



test_EM = est_G_params(Y, S, 
                       EM_iterations = 15,
                       NR_iterations = NA,
                       initial_mu = 0,
                       initial_tau = 1,
                       initial_pi0 = .5)
tail(test_EM)



taus = c()
for(NR_iter in 0:30) {
  test_EM = est_G_params(Y, S, 
                         EM_iterations = 15,
                         NR_iterations = NR_iter,
                         initial_mu = 0,
                         initial_tau = 1,
                         initial_pi0 = .5)
  taus = c(taus, test_EM[nrow(test_EM), 'tau'])
}

plot(taus)








