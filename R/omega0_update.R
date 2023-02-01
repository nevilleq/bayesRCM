#' Function to update group-level sparse covariance matrix Omega_0 / precision matrix Omega_0
#'
#' @param omega_0 Matrix, current precision matrix Omega
#' @param D Matrix, sum_k Omega_k / tau_k
#' @param nu Numeric, sum_k tau_k
#' @param lambda_3 Numeric, penalty on L1 Norm of omega_0
#'
#' @return NULL
#' @export
#'
#' @examples
omega0_update <- function(omega_0, D, nu, lambda_3) {
  #Grab parameter (p dim Omega / ROI network)
  p <- nrow(omega_0) #network/graph size

  #A. Sample Diagonal Elements
  for (j in 1:p) {
      reorder <- c(setdiff(1:p, j), j)
      o_pt    <- omega_0[reorder, reorder]

      o12 <- matrix(o_pt[p,-p], nrow = 1)
      o2i <- matinv(o_pt[-p,-p])

      c <- matABA(o12, o2i)             

      tmp <- matrix(c(matprod(o12, o2i), -1), nrow = 1)
      d   <- matABA(tmp, D[reorder, reorder])   # tmp %*% D[reorder,reorder] %*% t(tmp)
      #D_r <- D[reorder, reorder]

      omega_0[j,j] <- GIGrvg::rgig(1, 1 - nu/2, d, lambda_3) + c #verify rgig documentation
      #Checked with documentation vs. paper, appears to be parameterized the same
  }

  #B. Sample off-Diagonal Elements
  #pct accepted via mcmc step-proposal
  accept <- vector(mode = "logical", length = 0L)

  #For upper off diagnonal elements (reordered)
  for (j in 1:(p - 1)) {
    for (l in (j + 1):p) {
        reorder <- c(setdiff(1:p,c(j,l)),j,l)
        o_pt    <- omega_0[reorder,reorder]

        o12 <- o_pt[(p - 1):p, 1:(p - 2)]
        o2i <- matinv(o_pt[1:(p - 2), 1:(p - 2)])

        C    <- matABA(o12, o2i)  # A %*% B %*% t(A)
        tmp  <- cbind(matprod(o12, o2i), -diag(1, 2))
        D_c  <- matABA(tmp, D[reorder,reorder])

        #Obtain proposal shift and bound
        a     <- omega_0[j,j] - C[1, 1] #L
        b     <- omega_0[l,l] - C[2, 2] #R
        c     <- -C[1, 2]  #proposed change
        bound <- sqrt(a*b) #boundary for proposal(s)

        #Approximate the conditional distribution with step function
        #Proposal from step-proposal distribution, dependent on bound
        #1. Set up stepwise grid(s) via bound & scale down
        step_size     <- 2 * bound / 100
        kappa_grid    <- seq((step_size / 2) - bound, bound, by = step_size) #Between boundary

        #Compute the log_prob over grid, scale down, and obtain normalizing constant relative to stepsize
        log_prob_grid <- g0_log_density(kappa_grid, nu, a, b, D_c) - lambda_3 * abs(kappa_grid - c)
        scale_max     <- max(log_prob_grid) #Scale down by max
        log_prob_grid <- log_prob_grid - scale_max #Scaled grid
        norm_const    <- sum(exp(log_prob_grid) * step_size) #Discrete normalizing constant via stepsize

        #2. If 0 lies between the boundary -- (p_new = density, q_new = proposal, step = update)
        if (abs(c) < bound) {
          #Grab probability for MH-step
          tmp  <- g0_log_density(c, nu, a, b, D_c) - log(lambda_3 / 2) - scale_max - log(norm_const) #logit
          prob <- 1 / (1 + exp(-tmp)) #expit

          #If less than prob, shrink to zero (see 4. below, s-c=0)
          if (runif(1) < prob) {
            #Update sample, proposal, density
            prop  <- c #Set to c, which will shrink update to zero below
            q_new <- 1 #w/probability 1
            p_new <- 1 #w/probability 1
          } else {
            #Else generate update from step proposal
            grid_prop <- sample(kappa_grid, 1, prob = exp(log_prob_grid)) #sample along step function

            #Update sample, proposal, density
            prop  <- grid_prop + (runif(1) - 0.5) * step_size #Discrete MH step proposal
            q_new <- g0_log_density(grid_prop, nu, a, b, D_c) - lambda_3 * abs(grid_prop - c) #step-function sample
            p_new <- g0_log_density(prop, nu, a, b, D_c) - lambda_3 * abs(prop - c) #proposal
          }
        } else { #3. If interval does not include c, do same as else{} above & generate proposal
            #Generate from proposal grid
            grid_prop <- sample(kappa_grid, 1, prob = exp(log_prob_grid)) #sample along step function

            #Update sample, proposal, density
            prop  <- grid_prop + (runif(1) - 0.5) * step_size #Discrete MH step proposal
            q_new <- g0_log_density(grid_prop, nu, a, b, D_c) - lambda_3 * abs(grid_prop - c) #probability of step-function sample
            p_new <- g0_log_density(prop, nu, a, b, D_c) - lambda_3 * abs(prop - c) #probability of proposal
        }

        #4. Calculate the posterior density p_dens, proposal density q_dens of current value
        current <- omega_0[j,l] + c

        if (abs(omega_0[j,l]) < 0.001) {
          #Accept with probability 1
          q_cur <- 1
          p_cur <- 1
        } else {
          #Find closest grid index to current value
          grid_ind <- which.min(abs(kappa_grid - current))
          q_cur    <- g0_log_density(kappa_grid[grid_ind], nu, a, b, D_c) - lambda_3 * abs(kappa_grid[grid_ind] - c)
          p_cur    <- g0_log_density(current, nu, a, b, D_c) - lambda_3 * abs(current - c)
        }

        #5. MH-step
        #Acceptance & Log-ratio
        log_ratio <- p_new - q_new + q_cur - p_cur

        #If accepted
        if(runif(1) < exp(log_ratio)) {
          #Update omega_0 entry
          omega_0[j, l] <- omega_0[l, j] <- prop - c
          accept    <- c(accept, TRUE)
        } else {
          accept    <- c(accept, FALSE)
        }
    } #end l loop
  } #end j loop

  #Return omega_0 and the acceptance rate for this update iteration (of off diagonal elements)
  return(list(omega = omega_0, pct_accept = mean(accept)))
}



#' Log density of omega_{0, (i, j)} posterior given bounds a, b
#'
#' @param kappa
#' @param nu
#' @param a
#' @param b
#' @param D
#'
#' @return
#' @export
#'
#' @examples
g0_log_density <- function(kappa, nu, a, b, D) {
  #calculate log of g function(step, bounds (a,b))
  delta <- kappa^2/(a*b)
  res   <- -nu/2*log(1-delta) - (D[1,1]/a+D[2,2]/b-2*D[1,2]*kappa/(a*b))/(1-delta)/2
  return(res)
}

