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
  
  # alpha_tau = 50; tau_vec = runif(20, 0, 100); lambda_2 = 1; mu_tau = 30; sigma_tau = 3; trunc = c(0, 100);
  # window = 2;

  #Propose new alpha_tau (shape parameter in tau_k reg. gamma dist)
  alpha_prop   <- truncdist::rtrunc(spec = "norm", a = trunc[1], b = trunc[2], n = 1, mean = alpha_tau, sd = window)
  #alpha_prop   <- alpha_tau + rnorm(1, sd = window)
  
  if (alpha_prop < 0) {
    accept <- FALSE
  } else {
    #Compute MH transition probability
    log_diff  <- log_alpha_posterior(alpha_prop, tau_vec, lambda_2, mu_tau, sigma_tau, trunc) -
                  log_alpha_posterior(alpha_tau, tau_vec, lambda_2, mu_tau, sigma_tau, trunc)

    # prop_diff <- dtrunc_norm(alpha_tau, mean = alpha_prop, sd = window, a = trunc[1], b = trunc[2], log = TRUE) - 
    #               dtrunc_norm(alpha_prop, mean = alpha_tau, sd = window, a = trunc[1], b = trunc[2], log = TRUE)
    prop_diff <- log(truncdist::dtrunc(spec = "norm", a = trunc[1], b = trunc[2],
                                       x = alpha_tau, mean = alpha_prop, sd = window)) -
                  log(truncdist::dtrunc(spec = "norm", a = trunc[1], b = trunc[2],
                                        x = alpha_prop, mean = alpha_tau, sd = window))
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

log_alpha_posterior <- function(alpha_tau, tau_vec, lambda_2, mu_tau, sigma_tau, trunc = c(0, 100), type = "mode") {
  #Parameters
  nu    <- pmax(trunc[2] - tau_vec, trunc[1] + 0.000000001)
  k     <- length(nu)
  sigma <- sigma_tau * lambda_2
  
  #Tau prior truncated gamma with
  log_tau <- map_dbl(
              .x = nu, 
              ~truncdist::dtrunc(spec = "gamma", x = .x, a = trunc[1], b = trunc[2],
                                 shape = alpha_tau, rate = lambda_2, log = FALSE) |>
               log()
            )
  
  #Alpha prior log pdf 
  if (type %in% "mean") {
    log_alpha <- truncdist::dtrunc(spec = "norm", a = trunc[1], b = trunc[2],
                                   x = alpha_tau, mean = mu_tau, sd = sigma) |> log()
  } else if (type %in% "mode") {
    log_alpha <- truncdist::dtrunc(spec = "norm", a = trunc[1], b = trunc[2],
                                   x = alpha_tau, mean = mu_tau + 1, sd = sigma) |> log()
  } else {
    stop("type must be one of 'mean' or 'mode'")
  }

  #Sum for log pdf
  log_pdf <- sum(log_tau) + k * log_alpha
  
  #Return log pdf
  return(log_pdf)
}
