#' Function to update lambda_2 with Gamma prior
#'
#' @description Function to update 
#'
#' @param alpha_tau numeric > 0, current value of alpha_tau shape parameter in tau_k prior
#' @param tau_k numeric > 0, current value of tau for Kth subject
#' @param lambda_2 numeric > 0, parameter to the exponential prior Exp(lambda_2)
#' @param window numeric > 0, parameter to set step size for tau update proposals from log-normal
#' @param trunc  numeric > 0, truncation for regularization parameter tau_k
#'
#' @return Returns list of MH-step updated (or not) lambda_2 scalar value and a logicial T/F whether the proposal was accepted
#' @export
#'
#' @examples
#' lambda2_update(lambda_2, alpha_tau, tau_vec, mu_tau, sigma_tau, window, trunc = c(0, 100))
lambda2_update <- function(lambda_2, alpha_tau, tau_vec, mu_tau, sigma_tau, window, trunc = c(0, 100), ebayes = FALSE, power = 1) {
  
   # alpha_tau = 50; tau_vec = runif(20, 0, 100); lambda_2 = 1; mu_tau = 30; sigma_tau = 3; trunc = c(0, 100);
   # window = 0.1;

  #Propose new lambda2 (rate parameter in tau_k reg. gamma dist)
  lambda_prop   <- truncdist::rtrunc(spec = "norm", a = 0, b = Inf, n = 1, mean = lambda_2, sd = window)
  #lambda_prop   <- lambda_2 + rnorm(1, sd = window)
  
  if (lambda_prop <= 0) {
    accept <- FALSE
  } else {
    #Compute MH transition probability
    log_diff  <- log_lambda_posterior(lambda_prop, alpha_tau, tau_vec, mu_tau, sigma_tau, trunc) -
                  log_lambda_posterior(lambda_2, alpha_tau, tau_vec, mu_tau, sigma_tau, trunc)

    # prop_diff <- dtrunc_norm(lambda_2, mean = lambda_prop, sd = window, a = 0, b = Inf, log = TRUE) - 
    #               dtrunc_norm(lambda_prop, mean = lambda_2, sd = window, a = 0, b = Inf, log = TRUE)
    prop_diff <- log(truncdist::dtrunc(spec = "norm", a = 0, b = Inf, x = lambda_2, mean = lambda_prop, sd = window)) -
                  log(truncdist::dtrunc(spec = "norm", a = 0, b = Inf, x = lambda_prop, mean = lambda_2, sd = window))

    #MH Step
    if (log(runif(1)) < (log_diff + prop_diff)) {
      #Accept new proposal
      lambda_2 <- lambda_prop
      accept <- TRUE
    } else {
      accept <- FALSE
    }
  }

  #Return list of tau_k and acceptance
  return(list(lambda_2 = lambda_2, accept = accept))

}

log_lambda_posterior <- function(lambda_2, alpha_tau, tau_vec, mu_tau, sigma_tau, trunc = c(0, 100), type = "mode", ebayes = FALSE, power = 1) {
  #Parameters
  nu    <- pmax(trunc[2] - tau_vec, trunc[1] + 0.000000001)
  k     <- length(nu)
  sigma <- sigma_tau * lambda_2
  mu    <- mu_tau * lambda_2
  
  #Tau prior truncated gamma with
  log_tau <- map_dbl(
              .x = nu, 
              ~truncdist::dtrunc(spec = "gamma", x = .x, a = trunc[1], b = trunc[2],
                                 shape = alpha_tau, rate = lambda_2, log = TRUE) #update power transform
             )

  #Alpha prior log pdf 
  if (type %in% "mean") {
    log_alpha <- truncdist::dtrunc(spec = "norm", a = trunc[1], b = (trunc[2] * lambda_2),
                                   x = alpha_tau, mean = mu, sd = sigma) |> log()
  } else if (type %in% "mode") {
    log_alpha <- truncdist::dtrunc(spec = "norm", a = trunc[1], b = (trunc[2] * lambda_2 + 1),
                                   x = alpha_tau, mean = mu + 1, sd = sigma) |> log()
  } else {
    stop("type must be one of 'mean' or 'mode'")
  }
  
  
  #Lambda prior log pdf
  log_lambda <- dgamma(lambda_2, shape = 1, rate = 1/10, log = TRUE) # 1/power or power?
  
  #Sum/*k for log posterior pdf
  if(ebayes) {
    log_pdf <- sum(log_tau) + log_lambda
  } else {
    log_pdf <- sum(log_tau) + log_alpha + log_lambda #no K on 
  }
  
  #Return log pdf
  return(log_pdf)
}
