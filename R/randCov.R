#' Random Covariance Model
#'
#' This function implements the Random Covariance Model (RCM) for joint estimation of
#' multiple sparse precision matrices. Optimization is conducted using block
#' coordinate descent.
#' @param x List of \eqn{K} data matrices each of dimension \eqn{n_k} x \eqn{p}.
#' @param lambda1 Non-negative scalar. Induces sparsity in subject-level matrices.
#' @param lambda2 Non-negative scalar. Induces similarity between subject-level matrices and group-level matrix.
#' @param lambda3 Non-negative scalar. Induces sparsity in group-level matrix.
#' @return A list of length 2 containing:
#' \enumerate{
#' \item Group-level precision matrix estimate (Omega0).
#' \item \eqn{p} x \eqn{p} x \eqn{K} array of \eqn{K} subject-level precision matrix estimates (Omegas).
#' }
#' @author
#' Lin Zhang
#'
#' @export
#'
#' @references
#' Zhang, Lin, Andrew DiLernia, Karina Quevedo, Jazmin Camchong, Kelvin Lim, and Wei Pan.
#' "A Random Covariance Model for Bi-level Graphical Modeling with Application to Resting-state FMRI Data." 2019. https://arxiv.org/pdf/1910.00103.pdf
randCov <- function(x, lambda1, lambda2, lambda3 = 0) {

  K = length(x)
  p = dim(x[[1]])[2]
  S = sapply(x,cov,simplify="array")


  # criteria for convergence
  delta = 0.0001

  # inital values
  Omega0 = solve(apply(S,1:2,mean) + diag(1,p)*0.01)

  Omegas = sapply(plyr::alply(S,3), function(x1) solve(x1+diag(1,p)*0.01), simplify = 'array' )


  Omega0.old = matrix(0,p,p)
  Omegas.old = array(0,c(p,p,K))
  count = 0

  rho = lambda1/(1+lambda2/K)

  # start BCD algorithm
  while(max(abs(Omega0-Omega0.old))>delta | max(abs(Omegas-Omegas.old))>delta) {

    # record current Omega0 & Omegas
    Omega0.old = Omega0
    Omegas.old = Omegas

    # 1st step:

    sk = sapply(plyr::alply(S,3), function(x1) (solve(Omega0)*lambda2/K + x1)/(1+lambda2/K), simplify=FALSE)
    for(k in 1:K) {
      Omegas[,,k] = glasso::glasso(sk[[k]],rho)$wi
    }

    # 2nd step:
    if(lambda3 == 0) Omega0 = apply(Omegas,1:2,mean) else {

      s0 = apply(Omegas,1:2,mean)
      log <- capture.output({
        Omega0 = spcov_bcd(s0,lambda3)
      })
      # Omega0 = spcov(Omegas.old,s0,lambda,step.size=100)$Sigma

    }


    # record BCD iterations
    count = count+1
    # cat('iteration',count,'done \n')

    if(count>100) {
      cat('Omegas fail to converge for lam1 =',rho,'lam2 =',lambda2/K,'lam3 =',lambda3,'\n')
      break
    }
  }


  res = list(Omega0,Omegas)
  names(res)=c("Omega0","Omegas")
  return(res)
}


#' Covariance Graphical Lasso
#'
#' This function implements the Random Covariance Model (RCM) for joint estimation of
#' multiple sparse precision matrices. Optimization is conducted using block
#' coordinate descent.
#' @param samp_cov \eqn{p} x \eqn{p} sample covariance matrix.
#' @param rho Non-negative scalar. Induces sparsity in covariance matrix.
#' @param initial Initial value for covariance matrix.
#' @return  \eqn{p} x \eqn{p} sparse covariance matrix estimate.
#'
#' @references
#' Wang, Hao.
#' "Two New Algorithms for Solving Covariance Graphical Lasso Based on Coordinate Descent and ECM." 2012. https://arxiv.org/pdf/1205.4120.pdf
spcov_bcd <- function(samp_cov, rho, initial = NULL) {

  p = dim(samp_cov)[1]

  if(is.null(initial)) Sigma=samp_cov+0.01*diag(p) else Sigma = initial

  delta <- 0.0001

  Sigma.old = matrix(0,p,p)
  count2 = 0
  while(max(abs(Sigma-Sigma.old)) > delta) {  # loop 1: Sigma convergence

    Sigma.old <- Sigma
    for(i in 1:p) {  # loop 2: update each row/column of Sigma

      Omega11 <- solve(Sigma[-i,-i])
      beta <- Sigma[-i,i]

      S11 <- samp_cov[-i,-i]
      s12 <- samp_cov[-i,i]
      s22 <- samp_cov[i,i]

      a <- t(beta)%*%Omega11%*%S11%*%Omega11%*%beta - 2*t(s12)%*%Omega11%*%beta + s22

      if(rho == 0) gamma <- a else
        if(c(a) < 10^-10) gamma <- a else
          gamma <- (-1/(2*rho)+(1/(4*rho^2)+c(a)/rho)^0.5)

      V <- Omega11%*%S11%*%Omega11/gamma + rho*Omega11
      u <- t(s12)%*%Omega11/gamma

      beta.old <- 0
      while(max(abs(beta-beta.old)) > delta) {  # loop 3: off-diagonals convergence
        beta.old <- beta
        for(j in 1:(p-1)) {  # loop 4: each element
          temp = u[j] - V[j,-j]%*%beta[-j]
          beta[j] <- sign(temp)*max(0,abs(temp)-rho)/V[j,j]
        }  # loop 4
      }  # loop 3

      Sigma[i,-i] <- t(beta)
      Sigma[-i,i] <- beta
      Sigma[i,i] <- gamma + t(beta)%*%Omega11%*%beta

    }  # loop 2

    # record spcov iterations
    count2 = count2+1

    if(count2>100) {
      cat('Omega0 fails to converge for lam1 =',rho)
      break
    }

  }  # loop 1

  return(Sigma)

}
