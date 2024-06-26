#' Independent glasso method for starting values
#'
#' @param x
#' @param lambda
#'
#' @return
#' @export
#'
#' @examples
ind_graphs <- function(x,lambda) {

    # Inputs:
    K = length(x)
    p = dim(x[[1]])[2]
    S = sapply(x, cov, simplify = FALSE)

    # lambda = 0.1

    # create array to store results
    Omegas = rep(list(NA),K)

    # estimate independently using glasso
    for (k in 1:K) {
        Omegas[[k]] = glasso(S[[k]], lambda, penalize.diagonal = TRUE)$wi
	    Omegas[[k]] = (Omegas[[k]] + t(Omegas[[k]]))/2   # for computational stability
    }

    return(Omegas)

}
#' #BIC Calculations for tuning the lambda penalty in ind_graphs
#'
#' @param x
#' @param Omegas
#' @param Gk_est
#'
#' @return
#' @export
#'
#' @examples
bic_cal <- function(x, Omegas, Gk_est=NULL) {


    nk = sapply(x, nrow)
    S  = sapply(x, cov, simplify = FALSE)

    p  = dim(x[[1]])[2]
    K  = length(x)

    if(is.null(Gk_est)) Gk_est = sapply(Omegas,function(x) abs(x) > 0.0001, simplify=FALSE)
    nedges = (sapply(Gk_est,sum)+p)/2

    bic <- mapply(function(x1,x2,x3,x4) (x1-1)*sum(diag(x2%*%x3)) - x1*log(det(x3)) + x4*log(x1),
                  nk, S,Omegas,nedges)

    return(sum(bic))

}

#' Modified BIC for RCM
#'
#' This function calculates the modified BIC for the random covariance model (RCM)
#' @param x List of \eqn{K} data matrices each of dimension \eqn{n_k} x \eqn{p}.
#' @param Omega0 \eqn{p} x \eqn{p} group-level precision matrix estimate.
#' @param Omegas \eqn{p} x \eqn{p} x \eqn{K} array of \eqn{K} number of estimated subject-level precision matrices.
#' @param lambda2 Non-negative scalar. Induces similarity between subject-level matrices and group-level matrix.
#' @param G0_est \eqn{p} x \eqn{p} group-level network estimate.
#' @param Gk_est \eqn{p} x \eqn{p} x \eqn{K} array of \eqn{K} number of estimated subject-level networks.
#' @return Numeric Modified BIC value
#'
#' @export
mbic_cal <- function(x, Omega0, Omegas, lambda2, G0_est=NULL, Gk_est=NULL) {
  # modified BIC
  nk = sapply(x, nrow)
  S = sapply(x, cov, simplify="array")
  
  p = dim(Omegas)[1]
  K = length(x)
  
  if(is.null(Gk_est)) Gk_est = (round(Omegas, 3) != 0) - array(diag(p), c(p, p, K))
  if(is.null(G0_est)) G0_est = (round(Omega0, 3) != 0) - diag(p)
  nedges = apply(Gk_est, 3, sum) / 2
  
  # tmp = lambda2 / nk
  df.r <- (nedges + p) / (1 + lambda2)
  df.f <- (sum(G0_est) / 2 + p) * lambda2 / (1 + lambda2)
  
  mbic <- mapply(function(x1, x2, x3) (x1 - 1) * sum(diag(x2 %*% x3)) - x1 * log(det(x3)), nk, plyr::alply(S, 3), plyr::alply(Omegas, 3))
  
  return(sum(mbic) + (sum(df.r) + df.f) * log(sum(nk)))
}

