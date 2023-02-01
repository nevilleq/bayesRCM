#' Function to update graph G, constraining the precision matrix Omega
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
graph_update <- function(row_col, df, D, v, S, adj, omega, lambda_1) {
    # Network size
    p <- nrow(omega)

    # reorder the row_col-th row and column to be first in the update and ensure uppertrianguler adj/omega
    # need to work with upper triangular omega and adj, and this permutation is necessary!!!
    reorder   <- c(row_col, setdiff(1:p, row_col))
    backorder <- match(1:p, reorder)

    D     <- D[reorder,reorder]
    S     <- S[reorder,reorder]
    omega <- omega[reorder,reorder]
    adj   <- adj[reorder,reorder]

    # Updated posterior params from data
    b_post <- df + v;
    D_post <- D + S;

    #Adjacency threshold / graph
    omega  <- omega * adj
    accept <- rep(FALSE, (p - 1))
    
    # Sample off-diagonal elements
    i <- 1; # current 1st row/col is old row_col-th row/col (i.e. the one we want to update)
    for (j in 2:p) {
        #1. Propose G'
        #a. calculate the logit of no edge vs an edge (p is prob of having no edge)
        w <- log_H(b_post, D_post, omega, i, j) + lambda_1

        #Obtain probability of edge
        p <- 1 / (1 + exp(w)) #expit --> probability of edge
        ij_cur  <- adj[i, j]
        ij_prop <- runif(1) <= p
        #print(paste0("ij_cur:", ij_cur))
        #print(paste0("ij_prop:", ij_cur))
        
        #b. If it's accepted (yay!)
        if (ij_prop != ij_cur) { #If accepted
            #Record acceptance & make proposal
            accept[j - 1]  <- TRUE
            adj_prop       <- adj #Start with old adjacency/graph
            adj_prop[i, j] <- adj_prop[j,i] <- ij_prop #replace w proposal from above

            tri_adj_prop <- adj_prop
            tri_adj_prop[lower.tri(tri_adj_prop, diag = T)] <- 0
            omega_prop <- BDgraph::rgwish(1, tri_adj_prop, df, D) #Propose via prior
            omega_prop <- (omega_prop + t(omega_prop))/2  # for computational stability


            # step3: update G using MH, accept G' with mh_prob ratio r2
            mh_prob <- log_GWish_NOij_pdf(df, D, omega_prop, i, j, ij_cur) -
                       log_GWish_NOij_pdf(df, D, omega_prop, i, j, ij_prop);

            if (log(runif(1)) < mh_prob) { #If accepted
                adj[i,j] <- adj[j,i] <- ij_cur <- ij_prop


                # step 4: update \omega_ij if G' is accepted via paper
                omega <- gwish_ij_update(b_post, D_post, omega, i, j, ij_cur)

                # If omega update is not symmetric, force symmetry by upper triangle via Matrix package
                if (!Matrix::isSymmetric(omega)) {
                  #print('Asymetric proposition, forcing symmetric with Matrix::forceSymmetric()')
                  omega <- Matrix::forceSymmetric(omega, uplo = "U") #Determined by upper triangle
                }
            }
        }

    } # end of j-loop

    # Shuffle row/col in matrices back to original order (put row_col ii back where it belongs)
    omega <- omega[backorder,backorder]
    adj   <- adj[backorder,backorder]

    #Return precision matrices and graphs per subject
    return(list(omega = omega, adj = adj, accept = mean(accept)))
}




#' Function to update an individual ij entry in G_k via method outline in Wang and Li (2011) "Efficient Gaussian Graphical Model Determination without Approximating Normalizing Constant of the G-Wishart distribution"
#'
#' @param b_post
#' @param D_post
#' @param omega
#' @param i
#' @param j
#' @param ij_cur
#'
#' @return
#' @export
#'
#' @examples
gwish_ij_update <- function(b, D, omega, i, j, ij_cur) {

    #Graph dimension of network/graph
    p <- nrow(omega)

    ## If no edge between (i,j), omega_ij=0
    if (ij_cur == 0) {

        omega[i,j] <- omega[j,i] <- 0;

        o12 <- matrix(omega[j,-j], nrow = 1)
        o22 <- omega[-j,-j]

        C <- matABinvA(o12, o22)

        #omega[j,j] = rWishart(1,b,1/D[j,j]) + c;
        omega[j,j] <- rchisq(1, b)/D[j,j] + C;

    } else { # If e_ij \in E, i.e. edge exists

        #Reorder & take Cholesky decomp
        reorder <- c(setdiff(1:p, c(i,j)), i, j)
        o_pt    <- omega[reorder, reorder]
        
        #Check to make sure symmetric, if not force by upper diag? (or lower?)
        if (!isSymmetric(o_pt)) {
          print("Not symmetric")
          o_pt <- as.matrix(Matrix::forceSymmetric(o_pt, uplo = "U"))
        }
        
        #Check to make sure PD
        if (any(eigen(o_pt)$values < 0)) {
          o_pt <- o_pt |>
            (\(x) {as.matrix(Matrix::nearPD(x)$mat)})() 
        }
        
        #Cholesky decomp (must be Sym, PD)
        #print(o_pt)
        R <- matchol(o_pt)

        #Posterior params
        m_post      <- -R[p - 1, p - 1] * D[i,j] / D[j,j]
        sig_post    <- 1 / sqrt(D[j, j])
        R[p - 1, p] <- rnorm(1) * sig_post + m_post;
        R[p, p]     <- sqrt(rgamma(1, b/2,rate = D[j, j]/2))

        #Update
        omega_update <- t(R[,(p - 1):p]) %*% R[ ,p];
        omega[i,j]   <- omega[j,i] <- omega_update[1]
        omega[j,j]   <- omega_update[2]
    }
    #Return
    return(omega)
}
