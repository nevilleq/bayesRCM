#' Function to update tau_k given an Exp(lambda_2) prior
#'
#' @param tau_k numeric > 0, current value of tau for Kth subject
#' @param omega_k matrix, current subject specific precision matrix
#' @param omega_0 matrix, current group-level precision matrix
#' @param lambda_2 #numeric > 0, parameter to the exponential prior Exp(lambda_2)
#' @param window #numeric > 0, parameter to set step size for tau update proposals from log-normal
#'
#' @return
#' @export
#'
#' @examples
tau_update <- function(tau_k, omega_k, omega_0, lambda_2, window) {

  #Propose new tau
  tau_prop  <- tau_k * exp(rnorm(1, sd = window))

  #Compute MH transition probability
  log_diff  <- log_tau_posterior(tau_prop, omega_k, omega_0, lambda_2) - log_tau_posterior(tau_k, omega_k, omega_0, lambda_2)
  prop_diff <- dlnorm(tau_k, log(tau_prop), window, log = TRUE) - dlnorm(tau_prop, log(tau_k), window, log = TRUE)

  #MH Step
  if (log(runif(1)) < (log_diff + prop_diff)) {
    #Accept new proposal
    tau_k  <- tau_prop
    accept <- TRUE
  } else {
    accept <- FALSE
  }

  #Return list of tau_k and acceptance
  return(list(tau_k = tau_k, accept = accept))

}

log_tau_posterior <- function(tau_k, omega_k, omega_0, lambda_2) {
  #Parameters
  p <- nrow(omega_k)
  v <- (tau_k + p + 1)

  #Posterior pdf combo of GWish and exp
  #pdf <- log(lambda_2) - lambda_2 * tau_k + log(dGWish) ?
  pdf <- -tau_k * p/2 * log(2) - log_multi_gamma(p, v/2) + v/2 * p * log(tau_k) +
        tau_k * (log(matdet(omega_k)) - log(matdet(omega_0))) - (mattr(matsolve(omega_0, omega_k))/2 + lambda_2) * tau_k

  #Return pdf
  return(pdf)
}

