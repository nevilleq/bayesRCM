#' Function to update precision matrix Omega
#'
#' @param df Integer, degrees of freedom of the Gwishart distribution
#' @param D Matrix, location parameter of the Gwishart distribution
#' @param adj Matrix, current adjacency matrix/graph constraining precision matrix Omega
#' @param omega Matrix, current precision matrix Omega
#'
#' @return NULL
#' @export
#'
#' @examples
gwish_update <- function(df, D, adj, omega) {
  # Parameters
  p     <- nrow(omega)
  omega <- omega*adj
  sigma <- matinv(omega)
  iso_node_id <- which(rowSums(adj) == 1) #isolated nodes


  # Sample isolated nodes
  for (cliq_id in iso_node_id) {
      #c = rWishart(1,df=b,matinv(D[cliq_id,cliq_id]))
      c = rchisq(1, df = b) / D[cliq_id,cliq_id]
      omega[cliq_id,cliq_id] = c
      sigma[cliq_id,cliq_id] = 1/c
  }

  # Edge-wise update
  for (i in 1:(p - 1)) { #Loop through edges
    for (j in (i + 1):p) {
      if (adj[i, j] == 1) { #If edge

        cliq_id <- c(i, j) #Clique
        cliq_n  <- 2; #N = 2, recommended in Wang and Li (2012)
        A       <- #A matrix from Wang and Li (2012)
          matrix(
           data = rWishart(
                    n  = 1,
                    df = b + cliq_n - 1 ,
                    matinv(D[cliq_id, cliq_id])),
           nrow = clique_n
          )

        #Elements w & without clique
        o12 <- omega[cliq_id, -cliq_id]
        s12 <- sigma[cliq_id, -cliq_id]
        s22 <- sigma[-cliq_id, -cliq_id]

        #Inverse elements
        s11_inv <- matinv(sigma[cliq_id, cliq_id])
        o22_inv <- Sig22 - matABA(t(Sig12),Sig11.inv)

        #C matrix from Wang and Li (2012)
        C = A + matABA(o12, o22_inv)

        #Inverse of difference between current and new Omega (2x2 block)
        delta = matinv(omega[cliq_id, cliq_id] - C)

        #Update omega (2x2 block)
        omega[cliq_id, cliq_id] = C

        # ## Rank 2 update Sigma
        # sbb = sigma[cliq_id, cliq_id]
        # aa    = matinv(delta - sbb)
        # sigma = sigma + matABA(sigma[ , cliq_id], aa)
      }
    }
  }
  #Return updated omega
  return(omega)
}
