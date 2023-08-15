#' Bayesian robust Random Covariance Model (RCM) for heirarchical graphical modeling
#'
#' @param data_list List of $K$ subject's residual BOLD time series array (array of V volumes of matrices in time)
#' @param priors List of prior distribution specifications, defaults to -
#' @param n_samples Integer, number of posterior samples to obtain
#' @param n_burn Integer, number of posterior samples to burn
#' @param n_cores Integer, number of cores to be run in parallel - if empty/NULL will run sequentially
#'
#' @return NULL
#' @export
#'
#' @examples
#' rcm::example_data_list
#' rcm(example_data_list)
rcm <- function(y = data_list, priors = NULL, n_samples = 100, n_burn = 10, n_cores = 4, n_updates = 2) {
  # sim_res <- sim_data(subjects = 20, volumes = 200, rois = 10, alpha_tau = 25, lambda_2 = 0.5,
  #                     prop_true_con = 1/5, n_flip = 1, write = FALSE)
  # y = sim_res$data_list; priors = NULL; n_samples = 100; n_burn = 0; n_cores = 4; n_updates = 2;
  # library(Rcpp)
  # library(RcppArmadillo)
  # library(GIGrvg)
  # library(glasso)
  # library(abind)
  # library(BDgraph)
  # library(bayesRCM)
  # library(doParallel)
  #Set up parallel clusters (later may want to have multiple nodes (chains) to process multiple subjects (cores))
  cl <- parallel::makeCluster(n_cores)
  on.exit(parallel::stopCluster(cl)) #When done stop cluster
  doParallel::registerDoParallel(cl) #Initialize clusters
  `%dopar%` <- foreach::`%dopar%` #Need to pass this function later
    
  #Grab no. of subjects, rois, volumes
  K  <- length(y)
  p  <- ncol(y[[1]])
  vk <- map_dbl(y, nrow)
  Sk <- map(.x = y, ~t(.x) %*% .x)

  #Set lambda 1-3 penalty gamma a, b  hyperparams (uninformative / flat)
  alpha <- c(1, 1, 1) #1 - G_k, 2 - tau_k, 3 - Omega_0 (glasso)
  beta  <- c(1/10, 1/10, 1/10)
  
  #Starting value for alpha_tau & lambda_2 (based off of subjects & a mean of 50)
  alpha_tau <- K
  lambda_2  <- alpha_tau / 50
  mu_tau    <- 50 #Midpoint of truncation, prior mean/mode
  sigma_tau <- 20 #Prior sd, can iteratively update based on tau est

  #MCMC step size/window
  step_tau     <- rep(10, K)
  step_alpha   <- sqrt(2)
  step_lambda2 <- 0.05
  
  #How many updates during burn to adapt window/step-size
  n_updates <- floor(n_samples / n_updates)

  #Initialize estimates for Omega_k, Omega_0
  #Omega_k - tune lambda by bic and then grab G and Omega_0
  lambda_grid <- 10^seq(-1, 0.5, length.out = 20)
  bic         <- vector(mode = "numeric", length = 0)

  #Tune lambda for independent glasso
  for (lam in lambda_grid) {
    omega_k <- ind_graphs(y, lam)
    bic   <- c(bic, bic_cal(y, omega_k))
  }
  
  #Compute initial est. for omega_k, G_k/adj_k via best mBIC, and Omega_0
  omega_k <- ind_graphs(y, lambda_grid[which.min(bic)])
  adj_k   <- map(.x = omega_k, ~abs(.x) > 0.001)
  omega_0 <- Reduce("+", omega_k) / K
  sigma_0 <- solve(omega_0)

  #Initialize estimates for tau_k
  #Tau vector of subject specific regularization param on Omega_0
  tau_vec <- vector(mode = "numeric", length = K)

  #Iterate over each subject, find optimal tau_k based on posterior in 1D
  for (k in 1:K) {
    #print(k)
    #Tau posterior for fixed omega_k, omega_0, and lambda_2 = 0
    f_opt <- function(tau) {
      -1 * log_tau_posterior(tau, omega_k[[k]], sigma_0, alpha_tau, lambda_2, m_iter = 100)
    }
    #Optimize in 1D
    tau_vec[k] <- c(optimize(f_opt, interval = c(0, 100), tol = 0.01)$min)
  }

  #Initialize alpha_tau from tau_vec
  #Alpha_tau posterior for fixed tau_k, lambda_2 = 0, mu_tau, sigma_tau
  f_opt <- function(alpha_tau) {
    -1 * log_alpha_posterior(alpha_tau, tau_vec, lambda_2 = K / 50, mu_tau = mu_tau, sigma_tau = sigma_tau, trunc = c(0, 100), type = "mode")
  }
  #Optimize in 1D
  alpha_tau <- c(optimize(f_opt, lower = 1, upper = 100, tol = 0.01)$min)
  
  #Initialize lambda_2 from tau_vec
  #Alpha_tau posterior for fixed tau_k, lambda_2 = 0, mu_tau, sigma_tau
  f_opt <- function(lambda_2) {
    -1 * log_lambda_posterior(lambda_2, alpha_tau, tau_vec, mu_tau, sigma_tau, trunc = c(0, 100), type = "mode")
  }
  #Optimize in 1D
  lambda_2 <- c(optimize(f_opt, lower = 0, upper = 10, tol = 0.01)$min)
  
  
  #Set up storage for results
  #Omegas
  omegas_res <- array(NA, c(p * (p + 1) / 2, K, n_samples))
  accept_k   <- vector(mode = "numeric", length = K)
  omega0_res <- array(NA, c(p * (p + 1) / 2, n_samples))
  pct_omega_acc <- vector(mode = "integer", length = n_samples)
  pct_k_acc  <- matrix(NA, nrow = n_samples, ncol = K)

  #Taus, Alpha, Lambda_2
  accept_tau     <- matrix(NA, nrow = 0, ncol = K)
  accept_alpha   <- matrix(NA, nrow = 0, ncol = 1)
  accept_lambda2 <- matrix(NA, nrow = 0, ncol = 1)
  step_tau_mat   <- step_tau #Adaptive window for MH tau
  step_alpha_mat <- step_alpha 
  step_lam2_mat  <- step_lambda2
  tau_res        <- array(NA, c(K, n_samples))
  alpha_res      <- vector(mode = "numeric", length = n_samples)
  
  #Lambdas
  lambda_res <- array(NA, c(3, n_samples), dimnames = list(str_c("lambda_", 1:3)))

  #Set timer
  timer   <- 0
  t_start <- proc.time()
  n_iter  <- (n_burn + n_samples)
  #n_iter  <- 3

  #Loop through sampling algorithm n_samples + n_burn # times
  for (t in 1:n_iter) {
    #Print iteration for early testing
    print(paste0("Iteration: ", t))
    
    #Update Lambdas via direct sampling
    #Lambda 1 sparsity-inducing penalty on G_k
    card_k   <- (sapply(adj_k, sum) - p)/2 #Cardinality of G_k / # edges
    lambda_1 <- rgamma(1, alpha[1] + K, rate = beta[1] + sum(card_k))

    #Lambda 2 Gamma rate parameter for df/shrinkage tau_k prior
    #lambda_2 <- rgamma(1, alpha[2] + K + 1, beta[2] + sum(tau_vec))
    lambda_next    <- lambda2_update(lambda_2, alpha_tau, tau_vec,
                                     mu_tau, sigma_tau, window = step_lambda2)
    lambda_2       <- lambda_next$lambda_2
    accept_lambda2 <- rbind(accept_lambda2, lambda_next$accept)

    #Lambda 3 Sparse L-1 penalty on group precision omega_0 prior
    card_0 <- (sum(abs(omega_0) > 0.001) + p) / 2 #Cardinality omega_0 / # Edges or non-zero elements
    lambda_3 <- rgamma(1, alpha[3] + card_0, beta[3] + sum(abs(omega_0))/2)

    #Invert for Covariance & randomly select row_col pair
    sigma_0 <- matinv(omega_0)
    row_col <- sample(1:p, 1)
    
    #Set up foreach:: combine into 2 list, multicombine = TRUE
    # my_combine <- function(x, ...) {
    #   lapply(seq_along(x),
    #          function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
    # }
    
    #Test list for non-parallel loop to debug graph_update
    #omega_k_update <- list()
    #Update G_k, Omega_k via modified BIPS proposal & update scheme - Wang and Li (2012)
    omega_Gk_update <- foreach::foreach(
      k         = 1:K, #Subject index
      #         .combine  = "c", #List output
      #         .multicombine = TRUE,
      #         .init     = list(omega_k = list(), accept = list()),
      .packages = c("bayesRCM", "BDgraph", "Matrix"), #Packages
      .noexport = c("graph_update","gwish_ij_update", "rgwish") #Functions necessary to export
    ) %dopar% {
      #for(k in 1:k) {
      print(paste0("Subject -- ", k))
      #Propose/update G_k, Omega_k via
      update <-
        graph_update(
          row_col  = row_col,
          df       = tau_vec[k] + 2,
          D        = sigma_0 * tau_vec[k],
          v        = vk[k],
          S        = Sk[[k]],
          adj      = adj_k[[k]],
          omega    = omega_k[[k]],
          lambda_1 = lambda_1
        )
      
      #Record acceptance rate
      accept_k[k] <- update$accept
      
      #Upper triangular and averaged transpose for computational stability
      tri_adj <- update$adj
      tri_adj[lower.tri(tri_adj, diag = T)] <- 0
      
      #Update omega_k
      omega_k_update <- BDgraph::rgwish(1, tri_adj, vk[k] + tau_vec[k] + 2, Sk[[k]] + sigma_0 * tau_vec[k])
      omega_k_update <- (omega_k_update + t(omega_k_update)) / 2 # for computational stability
      #omega_k_update[[k]] <- BDgraph::rgwish(1, tri_adj, vk[k] + tau_vec[k] + 2, Sk[[k]] + sigma_0 * tau_vec[k])
      #omega_k_update[[k]] <- (omega_k_update[[k]] + t(omega_k_update[[k]])) / 2 # for computational stability
      
      #Return updated omega_k 
      return(omega_k_update)
    }
    
    #Pass back to omega_k
    omega_k  <- omega_Gk_update
    #omega_k <- omega_k_update
    

    #Tau_k update
    tau_k <-
      map(
        .x = 1:K, #Iterate from index 1 to Ktau_k, omega_k, sigma_0, alpha_tau, lambda_2, window
        ~tau_update(tau_vec[.x], omega_k[[.x]], sigma_0, alpha_tau, lambda_2, step_tau[.x])
      ) #Return list object
    # for (k in 1:K) {
    #   print(paste0("sub: ", k))
    #   tau_update(tau_vec[k], omega_k[[k]], sigma_0, alpha_tau, lambda_2, step_tau[k])
    # }
    accept_tau <- rbind(accept_tau, tau_k %>% map_lgl("accept")) #Pull out acceptance
    tau_vec    <- tau_k %>% map_dbl("tau_k") #Pull out the numeric tau_k list object

    #Adaptive tau window/stepsize ~ variance/sigma in log normal
     if (t %% n_updates == 0 & t <= n_burn) {
        #Compute acceptance rate (colwise mean)
        accept_rate <- apply(accept_tau, 2, mean)
        #For each subject, adjust tau_k proposal (lognormal) step size
         for (k in 1:K) {
           if (accept_rate[k] > 0.75) { #If accepting to many, inc variance of proposal
            step_tau[k] <- min(20, step_tau[k] + 1)
           } else if (accept_rate[k] < 0.5) { #If not accepting enough, dec variance of proposal
             step_tau[k] <- max(1, step_tau[k] - 1)
            }
         }
         step_tau_mat  <- rbind(step_tau_mat, step_tau) #Record adaptive step sizes
         accept_tau    <- matrix(NA, nrow = 0, ncol = K) #Restart acceptance rate tracking
     }
    
    #Update alpha_tau
    alpha_next <- alpha_update(alpha_tau, tau_vec, lambda_2, mu_tau, sigma_tau, window = step_alpha)
    alpha_tau    <- alpha_next$alpha_tau
    accept_alpha <- rbind(accept_alpha, alpha_next$accept)
    
    

    #Update Omega_0 via Wang and Li (2012) + step-proposal distribution
    D       <- apply(mapply('*', omega_k, tau_vec, SIMPLIFY = 'array'), 1:2, sum)
    omega_0 <- omega0_update(omega_0, D, sum(tau_vec), lambda_3)
    pct_accept <- omega_0$pct_accept #Off-diagonal acceptance%
    omega_0 <- omega_0$omega #Precision matrix itself

    #Save those results after burn-in
    if(t > n_burn) {
      t_burn <- t - n_burn
      omegas_res[, , t_burn] <- sapply(omega_k, function(x) x[upper.tri(x, diag = TRUE)])
      omega0_res[, t_burn]   <- omega_0[upper.tri(omega_0, diag = TRUE)]
      tau_res[, t_burn]      <- tau_vec
      alpha_res[t_burn]      <- alpha_tau
      lambda_res[, t_burn]   <- c(lambda_1, lambda_2, lambda_3)
      pct_omega_acc[t_burn]  <- pct_accept
      pct_k_acc[t_burn, ]    <- accept_k
    }

    #Track temporal progress (every 20% progress update)
    if (t %% floor(0.2 * (n_samples + n_burn))) {
      t_now   <- proc.time()
      timer   <- c(timer, (t_now - t_start)[3])
      t_start <- t_now
    }
  }

  #Result list of results
  result <-
    list(
      omega_0    = omega0_res,
      omega_k    = omegas_res,
      omega_acc  = pct_omega_acc,
      tau_k       = tau_res,
      tau_acc     = accept_tau,
      tau_step    = step_tau_mat,
      alpha_tau   = alpha_res,
      alpha_acc   = accept_alpha,
      alpha_step  = step_alpha_mat,
      lambdas     = lambda_res, 
      lambda_acc  = accept_lambda2,
      lambda_step = step_lam2_mat,
      timer       = timer
    )
  
  #Return result
  return(result)
}

