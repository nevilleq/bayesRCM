library(gt)
library(MASS)
library(Rcpp)
library(RcppArmadillo)
library(glasso)
library(GIGrvg)
library(tidyverse)
library(bayesRCM)

#Sim data
sim_res <- sim_data()

#Test constant pdf
K <- 1
omega_0 <- sim_res$true_params$omega_0
omega_k <- sim_res$true_params$omega_k[[K]]
tau_k   <- sim_res$true_params$tau_k[K]
alpha_tau <- length(sim_res$true_params$omega_k)
lambda_2  <- alpha_tau / 50


log_tau_post_test <- function(tau_k, omega_k, sigma_0, alpha_tau, lambda_2, m_iter) {
  #Parameters
  #p   <- nrow(omega_k)
  b   <- tau_k + 2
  D   <- sigma_0 * tau_k
  adj <- abs(omega_k) > 0.001 
  tri_adj <- adj
  tri_adj[lower.tri(tri_adj, diag = T)] <- 0
  
  #Log-Gwish (un-normalized)
  log_unnorm_pdf <- (b - 2) / 2 * log(matdet(omega_k)) - mattr(matprod(D, omega_k)) / 2 + #Gwish posterior
    (alpha_tau - 1) * log(tau_k) - lambda_2 * tau_k #Tau prior
  log_gwish_norm <- -BDgraph::gnorm(tri_adj, b = b, D = D, iter = m_iter)
  
  #Combined pdf 
  pdf <- log_unnorm_pdf + log_gwish_norm
  
  #Posterior pdf combo of GWish and inv gamma prior on tau^-1
  # pdf <- -tau_k * p/2 * log(2) - log_multi_gamma(p, v/2) + v/2 * p * log(tau_k) + #\pi(Omega_k | G_k, tau_k, Omega_0)
  #         tau_k * (log(matdet(omega_k)) - log(matdet(omega_0))) - (mattr(matsolve(omega_0, omega_k))/2) * tau_k + #\pi(Omega_k | G_k, tau_k, Omega_0)
  #         (alpha_tau + 1) * log(tau_k) - lambda_2 * tau_k #Gamma prior 
  
  
  #Posterior pdf combo of GWish and exponential prior on tau_k
  # pdf <- -tau_k * p/2 * log(2) - log_multi_gamma(p, v/2) + v/2 * p * log(tau_k) +
  #   tau_k * (log(matdet(omega_k)) - log(matdet(omega_0))) - (mattr(matsolve(omega_0, omega_k))/2 + lambda_2) * tau_k
  
  #Return pdf
  #return(c(pdf = pdf, unnorm = log_unnorm_pdf, norm = log_gwish_norm))
  return(pdf)
}


tibble(
  tau_grid = seq(1, 100, length.out = 100),
  tau_res  = map_dbl(.x = tau_grid, ~log_tau_post_test(.x, omega_k, omega_0, alpha_tau  = subjects, lambda_2 = lambda_2, m_iter = 1000))
) %>%
  ggplot(aes(x = tau_grid, y = tau_res)) +
  geom_point() +
  geom_line() +
  theme_minimal()


window   <- 0.1

quick_sim <- function(window, K = 1) {
tau_k    <- rgamma(100, alpha_tau, lambda_2)
tau_prop <- tau_k * exp(rnorm(1, sd = window))
res_orig <- res_prop <- matrix(NA, nrow = length(tau_k), ncol = 3)

for (k in 1:length(tau_k)) {
res_orig[k, ] <- log_tau_post_test(tau_k[k], omega_k[[K]], omega_0, alpha_tau, lambda_2)
res_prop[k, ] <- log_tau_post_test(tau_prop[k], omega_k[[K]], omega_0, alpha_tau, lambda_2)
}

cbind(res_orig, res_prop) %>%
  as.data.frame() %>%
  mutate(window = window, pdf_dif = V1 - V4, unnorm_dif = V2 - V5, norm_dif = V3 - V6,
         tau_k = tau_k, tau_prop = tau_prop) %>%
  dplyr::select(window:tau_prop) %>%
  return()
}

map_df(.x = c(0.01, 0.1, 0.2, 0.5), ~quick_sim(.x, K = 2)) %>%
  group_by(window) %>%
  summarise(across(where(is.numeric), mean))
