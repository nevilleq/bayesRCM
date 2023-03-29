#' Function to update tau_k given an Exp(lambda_2) prior
#'
#' @param alpha_tau numeric > 0, current value of alpha_tau shape parameter in tau_k prior
#' @param tau_k numeric > 0, current value of tau for Kth subject
#' @param lambda_2 numeric > 0, parameter to the exponential prior Exp(lambda_2)
#' @param window numeric > 0, parameter to set step size for tau update proposals from log-normal
#' @param trunc  numeric > 0, truncation for regularization parameter tau_k
#'
#' @return
#' @export
#'
#' @examples
alpha_update <- function(alpha_tau, tau_vec, lambda_2, mu_tau, sigma_tau, window, trunc = c(0, 100)) {
  
  alpha_tau = 50; tau_vec = runif(20, 0, 100); lambda_2 = 1; mu_tau = 30; sigma_tau = 3; trunc = c(0, 100);
  window = 2;

  #Propose new alpha_tau (shape parameter in tau_k reg. gamma dist)
  alpha_prop   <- truncnorm::rtruncnorm(1, mean = alpha_tau, sd = window, a = 0, b = Inf)
  alpha_prop   <- alpha_tau + rnorm(1, sd = window)
  
  if (alpha_prop < 0) {
    accept <- FALSE
  } else {
    #Compute MH transition probability
    log_diff  <- log_alpha_posterior(alpha_prop, tau_vec, lambda_2, mu_tau, sigma_tau, trunc) -
                  log_alpha_posterior(alpha_tau, tau_vec, lambda_2, mu_tau, sigma_tau, trunc)

    prop_diff <- dtrunc_norm(alpha_tau, mean = alpha_prop, sd = window, a = 0, b = Inf, log = TRUE) - 
                  dtrunc_norm(alpha_prop, mean = alpha_tau, sd = window, a = 0, b = Inf, log = TRUE)

    #MH Step
    if (log(runif(1)) < (log_diff + prop_diff)) {
      #Accept new proposal
      alpha_tau  <- alpha_prop
      accept <- TRUE
    } else {
      accept <- FALSE
    }
  }

  #Return list of tau_k and acceptance
  return(list(alpha_tau = alpha_tau, accept = accept))

}

log_alpha_posterior <- function(alpha_tau, tau_vec, lambda_2, mu_tau, sigma_tau, trunc = c(0, 100), type = "mean") {
  #Parameters
  nu    <- pmax(trunc[2] - tau_vec, trunc[1] + 0.000000001)
  k     <- length(nu)
  sigma <- sigma_tau * lambda_2
  
  #Tau prior truncated gamma with
  log_tau <- map_dbl(.x = nu, 
                     ~dtrunc_gamma(.x, shape = alpha_tau, rate = lambda_2,
                                   a = trunc[1], b = trunc[2], log = TRUE))
  
  #Alpha prior log pdf 
  if (type %in% "mean") {
    log_alpha <- dtrunc_norm(alpha_tau, mu_tau, sigma, a = 0, b = Inf, log = TRUE)
  } else if (type %in% "mode") {
    log_alpha <- dtrunc_norm(alpha_tau, mu_tau + 1, sigma, a = 0, b = Inf, log = TRUE)
  } else {
    stop("type must be one of 'mean' or 'mode'")
  }

  #Sum for log pdf
  log_pdf <- sum(log_tau) + k * log_alpha
  
  #Return log pdf
  return(log_pdf)
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