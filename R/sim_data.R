#' Simulate bayesRCM data with underlying G-Wishart data generating mechanism
#'
#' @param subjects K number of subjects
#' @param volumes V number of samples/fMRI observations (length of time-series)
#' @param rois P network dimension
#' @param prop_true_con Proportion of true connections 
#' @param seed Reproducible seed for the simulation
#' @param save Logical, if TRUE write save
#'
#' @return data_list
#' @export 
#'
#' @examples
sim_data <- function(subjects = 20, volumes = 200, rois = 10, prop_true_con = 1/5, trunc = c(0, 100),
                     alpha_tau = 50, lambda_2 = 2/5, n_flip = 1, seed = 4, write = FALSE) {
  
  #subjects = 20; volumes = 100; rois = 10; prop_true_con = 1/2; seed = 1; n_flip = 1;
  #A. G_0 --> Omega_0
  ##1. Start with diagonal p x p matrix
  G_0 <- diag(1, rois)
  
  ##2. To induce relative sparsity, make every 3rd off diagonal an edge
  #G_0[upper.tri(G_0)][seq(1, (subjects * (subjects - 1) / 2), by = 3)] <- TRUE 
  #G_0[upper.tri(G_0)] <- (rbinom(length(G_0[upper.tri(G_0)]), 1, prob = true_con / length(G_0[upper.tri(G_0)])) == 1)
  upper_G0 <- upper.tri(G_0, diag = FALSE)
  set.seed(seed)
  #G_0[upper_G0][seq(1, length(upper_G0), by = floor(1 / prop_true_con))] <- TRUE
  G_0[upper_G0] <- rbinom(sum(upper_G0), 1, prob = prop_true_con)
  
  ##3. Enforce symmetry by the upper trial
  G_0         <- G_0 + t(G_0) - diag(1, rois) #Subtract extra 1 from diagonal
  G_0_logical <- (G_0 == 1) #Make logical for below
  
  #Check rowise proportion of true connections so no dangling nodes
  #sum(G_0[upper_G0]) #43 true out of 190 off-diagonals
  #mean(G_0[upper_G0]) #~22% true connections
  #(apply(G_0, 2, sum) - 1)  #Number of true connections per ROI/row, all at least 2
  
  ##.4 Generate Omega_0 from G_0 (via Lin's simu_data.R schema)
  #Sample
  set.seed(seed)
  sample <- runif(rois^2, 0.5, 1)
  set.seed(seed)
  sample <- sample * sample(c(-1, 1), rois^2, replace = TRUE)
  
  #Threshold by G_0
  init     <- diag(0.5, rois) + G_0 * matrix(sample, nrow = rois, ncol = rois)
  init_sym <- init + t(init) + diag(rowSums(G_0 + t(G_0))/2) #Ensure symmetry & p.d. X + t(X) + diag()
  
  #True Omega_0
  Omega_0 <- diag(diag(init_sym)^-0.5) %*% init_sym %*% diag(diag(init_sym)^-0.5)
  Sigma_0 <- solve(Omega_0)
  
  #Check positive definite
  if (any(eigen(Omega_0, only.values = TRUE)$values <= 0) | any(eigen(Sigma_0, only.values = TRUE)$values <= 0)) { 
    #Warning
    warning("Warning: one of Omega_0/Sigma_0 not p.d.")
  }
  
  ##B. G_k, Omega_k ~ GWish(G_k, tau_k, Omega_0/tau_k)
  #Set up storage and params
  G_k       <- list()
  Omega_k   <- list()
  Sigma_k   <- list()
  Tau_k     <- vector(mode = "numeric", length = subjects)
  n_edges   <- n_flip #Off diagonal elements to flip/change for each subject
  #alpha_tau <- subjects #Fixed shape parameter for inv_tauk
  #lambda_2  <- subjects / 50 #inv tau from inv gamma => scale 50/subjects tau, E[tau_1:k] = subjects * (1/lam2) = 50
  
  #Loop through subjects
  for (k in 1:subjects) {
    #Set seed, sample an off diagonal element / edge to flip
    set.seed(k)
    new_edges <- sample(1:(rois * (rois - 1) / 2), n_edges)
    
    #1. (G_k) Randomly add or remove n_edge # of edges in upper.tri(G_0) to generate
    G_k[[k]] <- G_0_logical
    G_k[[k]][upper.tri(G_k[[k]], diag = FALSE)][new_edges] <- !G_k[[k]][upper.tri(G_k[[k]], diag = FALSE)][new_edges] #Flip those upper tri edges
    G_k[[k]] <- as.matrix(Matrix::forceSymmetric(G_k[[k]], uplo = "U"))
    
    #Check to make sure working appropriately
    if (sum(G_k[[k]] != G_0) != 2 * n_edges) {
      warning("woops! Something's wrong, not flipping the desired number of edges.")
    }
    
    #2. (Tau_k) Randomly sim from truncated, left tailed gamma 
    set.seed(k)
    Tau_k[k] <- trunc[2] - truncdist::rtrunc(spec = "gamma", a = trunc[1], b = trunc[2],
                                             n = 1, shape = alpha_tau, rate = lambda_2)

    #3. (Omega_k) ~ GWish(G_k, Tau_k, Omega_0/Tau_k) => Omega_k^{-1} = Sigma_0 for mvtnorm
    tri_adj <- G_k[[k]]
    tri_adj[lower.tri(tri_adj, diag = TRUE)] <- FALSE #Only upper tri adj for rgwish
    
    #Set seed & sample from RGwish
    set.seed(k)
    Omega_k[[k]] <- round(BDgraph::rgwish(1, adj = tri_adj, b = Tau_k[k] + 2, D = Sigma_0 * Tau_k[k]), 4)
    Sigma_k[[k]] <- solve(Omega_k[[k]])
  }
  
  ##C. (Y_ki) Data itself (data_list of length k subjects)
  #Create data list 
  data_list <- list()
  
  #Loop through subjects and volumes to generate data
  for (k in 1:subjects) { #assumes no temporal mean trend, centered at 0
    set.seed(k)
    data_list[[k]] <- mvtnorm::rmvnorm(volumes, rep(0, rois), Sigma_k[[k]])
  }
  
  #D. True parameters
  true_params <-
    list(
      omega_0   = Omega_0,
      omega_k   = Omega_k,
      tau_k     = Tau_k,
      alpha_tau = alpha_tau,
      lambda_2  = lambda_2
    )
  
  #E. Simulation Settings
  sim_settings <-
    tibble( 
      subjects  = subjects,
      volumes   = volumes,
      rois      = rois,
      prop_true_con = prop_true_con,
      n_flip    = n_flip
    )
  
  #F. If write is TRUE then create directory and write out simulated data list as .RDS
  if(write) {
    dir.create("./sim_data", showWarnings = FALSE)
    dir.create("./sim_truth", showWarnings = FALSE)
    write_rds(data_list, paste0("./sim_data/seed_", seed, ".rds"))
    write_rds(true_params, paste0("./sim_truth/seed_", seed, ".rds"))
  }
  
  #G. Return list of data and true params
  return(list(data_list = data_list, true_params = true_params, sim_settings = sim_settings))
  
}
