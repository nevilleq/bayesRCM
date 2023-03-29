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
tau_update <- function(tau_k, omega_k, sigma_0, alpha_tau, lambda_2, window, trunc = c(0, 100)) {

  #Propose new tau
  #tau_prop  <- tau_k * exp(rnorm(1, sd = window))
  tau_prop   <- truncnorm::rtruncnorm(1, mean = tau_k, sd = window, a = trunc[1], b = trunc[2])
  
  if (tau_prop < trunc[1] | tau_prop > trunc[2]) {
    accept <- FALSE
  } else {
    #Compute MH transition probability
    log_diff  <- log_tau_posterior(tau_prop, omega_k, sigma_0, alpha_tau, lambda_2) - log_tau_posterior(tau_k, omega_k, sigma_0, alpha_tau, lambda_2)
    prop_diff <- dtrunc_norm(tau_k, mean = tau_prop, sd = window, a = trunc[1], b = trunc[2], log = TRUE) - 
                  dtrunc_norm(tau_prop, mean = tau_k, sd = window, a = trunc[1], b = trunc[2], log = TRUE)
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

log_tau_posterior <- function(tau_k, omega_k, sigma_0, alpha_tau, lambda_2, trunc = c(0, 100), m_iter = 100) {
  #Parameters
  #p   <- nrow(omega_k)
  b   <- tau_k + 2
  D   <- sigma_0 * tau_k
  nu  <- pmax(trunc[2] - tau_k, trunc[1] + 0.000000001)
  adj <- abs(omega_k) > 0.001 
  tri_adj <- adj
  tri_adj[lower.tri(tri_adj, diag = T)] <- 0

  #Log-Gwish (un-normalized)
  log_unnorm_pdf <- (b - 2) / 2 * log(matdet(omega_k)) - mattr(matprod(D, omega_k)) / 2 + #Gwish prior
                      dtrunc_gamma(nu, shape = alpha_tau, rate = lambda_2, log = TRUE) #Tau prior
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

dtrunc_gamma <- function(x, shape, rate, a = 0, b = Inf, log = FALSE) {
  #If x \not\in (a, b) return 0
  if(x < a | x > b) {
    return(0) #outside domain/support
  }
  
  #PDF gamma scaled by CDF at truncation
  pdf <- dgamma(x, shape, rate) / 
    (pgamma(b, shape, rate) - pgamma(a, shape, rate))
  
  #If log return log(pdf) else return pdf
  if (log) {
    return(log(pdf))
  } else {
    return(pdf)
  }
}

dtrunc_norm <- function(x, mean, sd, a, b, log = FALSE) {
  #If x \not\in (a, b) return 0
  if(x < a | x > b) {
    return(0) #outside domain/support
  }
  
  #PDF normal scaled by CDF at truncation
  pdf <- dnorm(x, mean, sd) / (pnorm(b, mean, sd) - pnorm(a, mean, sd))
  
  #If log return log(pdf) else return pdf
  if (log) {
    return(log(pdf))
  } else {
    return(pdf)
  }
}