#' Function to update tau_k given an Exp(lambda_2) prior
#'
#' @param tau_k numeric > 0, current value of tau for Kth subject
#' @param omega_k matrix, current subject specific precision matrix
#' @param sigma_0 matrix, current group-level precision matrix
#' @param lambda_2 #numeric > 0, parameter to the exponential prior Exp(lambda_2)
#' @param window #numeric > 0, parameter to set step size for tau update proposals from log-normal
#'
#' @return
#' @export
#'
#' @examples
tau_update <- function(tau_k, omega_k, sigma_0, alpha_tau, lambda_2, window, trunc = 100) {

  #Propose new tau
  #tau_prop  <- tau_k * exp(rnorm(1, sd = window))
  tau_prop   <- tau_k * exp(truncnorm::rtruncnorm(1, a = -trunc, b = trunc, sd = window))
  
  if (tau_prop >= trunc) {
    accept <- FALSE
  } else {
    #Compute MH transition probability
    log_diff  <- log_tau_posterior(tau_prop, omega_k, sigma_0, alpha_tau, lambda_2) - log_tau_posterior(tau_k, omega_k, sigma_0, alpha_tau, lambda_2)
    #prop_diff <- dlnorm(tau_k, log(tau_prop), window, log = TRUE) - dlnorm(tau_prop, log(tau_k), window, log = TRUE)
    prop_diff <- log(EnvStats::dlnormTrunc(tau_k, log(tau_prop), sdlog = window, min = 0, max = trunc)) - 
                 log(EnvStats::dlnormTrunc(tau_prop, log(tau_k), sdlog = window, min = 0, max = trunc))
    if(is.nan(prop_diff)) {
      prop_diff <- -Inf
    }
    #MH Step
    if (log(runif(1)) < (log_diff + prop_diff)) {
      #Accept new proposal
      tau_k  <- tau_prop
      accept <- TRUE
    } else {
      accept <- FALSE
    }
  }

  #Return list of tau_k and acceptance
  return(list(tau_k = tau_k, accept = accept))

}

log_tau_posterior <- function(tau_k, omega_k, sigma_0, alpha_tau, lambda_2, trunc = 100, m_iter = 100) {
  #Parameters
  #p   <- nrow(omega_k)
  b   <- tau_k + 2
  D   <- sigma_0 * tau_k
  nu  <- min(200 - tau_k, 0.0001)
  adj <- abs(omega_k) > 0.001 
  tri_adj <- adj
  tri_adj[lower.tri(tri_adj, diag = T)] <- 0

  #Log-Gwish (un-normalized)
  log_unnorm_pdf <- (b - 2) / 2 * log(matdet(omega_k)) - mattr(matprod(D, omega_k)) / 2 + #Gwish posterior
                    (alpha_tau - 1) * log(nu) - lambda_2 * nu #Tau prior
  log_gwish_norm <- BDgraph::gnorm(tri_adj, b = b, D = D, iter = m_iter)
  
  #If is -Inf then don't use in posterior (will be inf for both proposal & current)
  if (is.infinite(log_gwish_norm)) {
    log_gwish_norm <- 0 #Has no affect on posterior since both constants will be -Inf
  }
  
  #Combined pdf 
  pdf <- log_unnorm_pdf - log_gwish_norm
  
  #Posterior pdf combo of GWish and inv gamma prior on tau^-1
  # pdf <- -tau_k * p/2 * log(2) - log_multi_gamma(p, v/2) + v/2 * p * log(tau_k) + #\pi(Omega_k | G_k, tau_k, Omega_0)
  #         tau_k * (log(matdet(omega_k)) - log(matdet(omega_0))) - (mattr(matsolve(omega_0, omega_k))/2) * tau_k + #\pi(Omega_k | G_k, tau_k, Omega_0)
  #         (alpha_tau + 1) * log(tau_k) - lambda_2 * tau_k #Gamma prior 
    

  #Posterior pdf combo of GWish and exponential prior on tau_k
  # pdf <- -tau_k * p/2 * log(2) - log_multi_gamma(p, v/2) + v/2 * p * log(tau_k) +
  #   tau_k * (log(matdet(omega_k)) - log(matdet(omega_0))) - (mattr(matsolve(omega_0, omega_k))/2 + lambda_2) * tau_k
  
  #Return pdf
  return(pdf)
}

