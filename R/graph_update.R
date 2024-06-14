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
  
  # k = 8;
  # row_col  = row_col;
  # df       = tau_vec[k] + 2;
  # D        = sigma_0 * tau_vec[k];
  # v        = vk[k];
  # S        = Sk[[k]];
  # adj      = adj_k[[k]];
  # omega    = omega_k_test[[k]];
  # lambda_1 = lambda_1;

  
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
    omegak1 <- omega
    omega  <- omega * adj
    accept <- rep(FALSE, (p - 1))
    
    # Sample off-diagonal elements
    i <- 1; # current 1st row/col is old row_col-th row/col (i.e. the one we want to update)
    for (j in 2:p) {
      print(j)
     # j = 8;
        #1. Propose G'
        #a. calculate the logit of no edge vs an edge (p is prob of having no edge)
        w <- log_H(b_post, D_post, omega, i, j) + lambda_1

        #Obtain probability of edge
        p <- 1 / (1 + exp(w)) #expit --> probability of edge
        print(paste0("p(Edge):", as.numeric(p)))
        ij_cur  <- adj[i, j]
        ij_prop <- runif(1) <= p
        print(paste0("ij_cur:", as.numeric(ij_cur)))
        print(paste0("ij_prop:", as.numeric(ij_prop)))
        
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
            #print("NOij")
            #log prob of current
            tryCatch({
              log_cur <- log_GWish_NOij_pdf(df, D, omega_prop, i, j, ij_cur)
            }, error = function(e) {
              conditionMessage(e)
            })
            #log prob of proprosal 
            tryCatch({
              log_prop <- log_GWish_NOij_pdf(df, D, omega_prop, i, j, ij_prop)
            }, error = function(e) {
              conditionMessage(e)
            })
            
            #Handle cases where null 
            if (is.null(log_cur) & (!is.null(log_prop))) {
              print("current null")
              mh_prob <- 0 #accept proposal
            } else if ((!is.null(log_cur)) & is.null(log_prop)) {
              print("prop null")
              mh_prob <- -Inf #reject
            } else if ((is.null(log_cur)) & is.null(log_prop)) {
              print("both null")
              mh_prob <- -Inf #reject
            } else {
              mh_prob <- log_cur - log_prop; #accept w prob mh_prob
            }

            if (log(runif(1)) < mh_prob) { #If accepted
                adj[i,j] <- adj[j,i] <- ij_cur <- ij_prop


                # step 4: update \omega_ij if G' is accepted via paper
                #print("Update")
                omega_temp <- tryCatch({
                  gwish_ij_update(b_post, D_post, omega, i, j, ij_cur)},
                  error = function(e) {
                  conditionMessage(e)
                })
                if(!is.null(omega_temp)) {
                  omega <- omega_temp
                }

                # If omega update is not symmetric, force symmetry by upper triangle via Matrix package
                # if (!Matrix::isSymmetric(omega)) {
                #   #print('Asymetric proposition, forcing symmetric with Matrix::forceSymmetric()')
                #   omega <- Matrix::forceSymmetric(omega, uplo = "U") #Determined by upper triangle
                # }
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

    #Test
    #i = 1; j=8; ij_cur = TRUE; b = b_post; D=D_post; omega = omega_prop;
  
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
        print(paste0("Det(o_pt): ", det(o_pt)))
        
        #Cholesky decomp (must be Sym, PD)
        if(any(eigen(o_pt)$values < 0)) {
          print("non-PD o_pt, trying again")
          o_ptpd <- o_pt |>
            (\(x) {as.matrix(Matrix::nearPD(x)$mat)})() #Find closest PD matrix
          R <- matchol(o_ptpd)  #cholesky decomp
        } else {
          R <- matchol(o_pt)
        }
        
        # #Cholesky decomp (must be Sym, PD)
        # R <- tryCatch({matchol(o_pt)}, error = function(e) {conditionMessage(e)})
        # if(is.character(R)) {
        #   print("No PD, trying again...")
        #   o_ptpd <- o_pt |>
        #     (\(x) {as.matrix(Matrix::nearPD(x)$mat)})()
        #   R <- matchol(o_ptpd) 
        # }

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

#Source scripts for testing
#Rcpp::sourceCpp("./src/bayesRCM_helper_funcs.cpp")

# log_GWish_NOij_pdf <- function(b, D, Omega, i, j, edgeij) {
#   # Compute log p(Omega\omega(i,j)) up to the normalizing constant of G-Wishart
# 
#   #b = b_post; D = D_post; i = 1; j = 2; edgeij = FALSE; Omega = omega_prop;
# 
#   if (edgeij == 0) {
#     Omega[i, j] <- 0
#     Omega[j, i] <- 0
# 
#     Ome12 <- Omega[j, -j] %>% matrix()
#     Ome22 <- Omega[-j, -j]
# 
#     c <- t(Ome12) %*% solve(Ome22) %*% Ome12
#     Omega_new <- Omega
#     Omega_new[j, j] <- c[1, 1]
# 
#     D1 <- as.matrix(D[j, j])
#     f <- -log_iwishart_InvA_const(b, D1) + (b - 2) / 2 * log(det(Ome22)) - mattr(D %*% Omega_new) / 2
#   } else {
#     Ome12 <- Omega
#     temp_row <- Ome12[i+1, ]
#     Ome12[i+1, ] <- Ome12[j, ]
#     Ome12[j, ]   <- temp_row
#     Ome12 <- Ome12[i:(i+1), -c(i, j)]
# 
#     Ome22 <- Omega[-c(i, j), -c(i, j)]
#     Ome11 <- matrix(c(Omega[i, i], Omega[i, j], Omega[i, j], Omega[j, j]), nrow = 2, byrow = TRUE)
# 
#     A <- Ome11 - Ome12 %*% solve(Ome22) %*% t(Ome12)
#     log_Joint <- (b - 2) / 2 * log(det(Omega)) - sum(diag((D %*% Omega))) / 2
# 
#     D_ij <- matrix(c(D[i, i], D[i, j], D[i, j], D[j, j]), nrow = 2, byrow = TRUE)
#     logK2by2 <- log_dWish(A, b, D_ij)
# 
#     V <- matinv(D_ij)
#     Dii <- as.matrix(1 / V[2, 2])
#     Aii <- as.matrix(A[1, 1])
#     logKii <- log_dWish(Aii, b + 1, Dii)
# 
#     f <- log_Joint + logKii - logK2by2
#   }
#   return(f)
# }
# 
# 
