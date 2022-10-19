#' Update the regularization/penalty terms lambda 1, 2, 3 on G_k, tau_k, and Omega_0
#'
#' @param lambda_1 numeric > 0, penalty on the cardinality of G_k
#' @param lambda_2 numeric > 0, exponential parameter in tau_k prior, regularization
#' @param lambda_3 numeric > 0, L1 regularization on ||Omega_0||_1
#' @param tau_k numeric > 0, current value of tau for subject k
#' @param G_k matrix, current graph or adjacency matrix for subject k
#' @param omega_0 #matrix, current overall precision matrix for all subjects
#'
#' @return
#' @export
#'
#' @examples
lambda_update <- function(lambda_1, lambda_2, lambda_3,
                          tau_k, adj_k, omega_0,
                          alpha, beta, p) {
  #Lambda 1
  cardinality_k <- (sapply(adj_k, sum) - p)/2
  lam1 <- rgamma(1, alpha[1] + K, rate = beta[1] + sum(cardinality_k))

        lam2 <- rgamma(1,a.lam+K,b.lam+sum(tau_k))

        q0 <- (sum(abs(Omega0)>0.0001)+p)/2
        lam3 <- rgamma(1,a.lam+q0, b.lam+sum(abs(Omega0))/2)

}
