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

