library(tidyverse)
library(bayesRCM)
library(gt)
library(MASS)
library(Rcpp)
library(RcppArmadillo)
library(glasso)
library(GIGrvg)
#library(doParallel)
library(Matrix)

#Generate covariance structure for multivariate gaussian covariance matrices 
volumes  <- 250 #Up the sample size to make sure it recovers 
subjects <- 10  #Number of subject
#Test alter the alpha_tau to induce outliers
rois     <- 25 #keep it small to start testing 
#true_con <- 30 #true connections or no. of edges in G
prop_true_con <- 1/5 #Proportion of true connections, easier to control than exact number
n_edges <- 1


sim_data <- function(subjects = 20, rois = 10, volumes = 250, prop_tru_con = 1/5, n_edges = 1) {

#Start with G_0, Omega_0 (Sigma_0^{-1})
#A. Group overall graph G_0with true_con # of edges
##1. Start with diagonal p x p matrix
G_0 <- diag(1, rois)

##2. To induce relative sparsity, make every 3rd off diagonal an edge
#G_0[upper.tri(G_0)][seq(1, (subjects * (subjects - 1) / 2), by = 3)] <- TRUE 
#G_0[upper.tri(G_0)] <- (rbinom(length(G_0[upper.tri(G_0)]), 1, prob = true_con / length(G_0[upper.tri(G_0)])) == 1)
upper_G0 <- upper.tri(G_0, diag = FALSE)
set.seed(8)
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
set.seed(4)
sample <- runif(rois^2, 0.5, 1)
set.seed(4)
sample <- sample * sample(c(-1, 1), rois^2, replace = TRUE)

#Threshold by G_0
init     <- diag(0.5, rois) + G_0 * matrix(sample, nrow = rois, ncol = rois)
init_sym <- init + t(init) + diag(rowSums(G_0 + t(G_0))/2) #Ensure symmetry & p.d. X + t(X) + diag()

#True Omega_0
Omega_0 <- diag(diag(init_sym)^-0.5) %*% init_sym %*% diag(diag(init_sym)^-0.5)
Sigma_0 <- solve(Omega_0)

#Check positive definite
# if (any(eigen(Omega_0, only.values = TRUE)$values <= 0) | any(eigen(Sigma_0, only.values = TRUE)$values <= 0)) { 
#   #Warning
#   warning("Warning: one of Omega_0/Sigma_0 not p.d.")
# }

##B. G_k, Omega_k ~ GWish(G_k, tau_k, Omega_0/tau_k)
#Set up storage and params
G_k       <- list()
Omega_k   <- list()
Sigma_k   <- list()
inv_Tau_k <- Tau_k <- vector(mode = "numeric", length = subjects)
#n_edges   <- 1 #Off diagonal elements to flip/change for each subject
alpha_tau <- subjects #Fixed shape parameter for inv_tauk
lambda_2  <- subjects / 50 #inv tau from inv gamma => scale 50/subjects tau, E[tau_1:k] = subjects * (1/lam2) = 50

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
  
  #2. (Tau_k) Randomly sim tau_inv from inv. gamma distribution => tau is Gamma(a, 1/scale)
  set.seed(k)
  #Tau_k[k] <- sum(rexp(alpha_tau, lambda_2)) #Inverse tau~Gamma(a, lambda_2)
 # inv_Tau_k[k] <- rgamma(1, alpha_tau, 1/5) #Inverse Gamma: if X ~ G(a, scale) then 1 / X is IG(a, 1/scale)
  Tau_k[k]     <- rgamma(1, alpha_tau, lambda_2) #alpha_tau a little bigger to get normal shape
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

#Write out
#write_rds(data_list, sprintf("./sim/data/tau_sim/%i_%i_%i_data_list.RDS", subjects, rois, volumes))
#write_rds(Tau_k, sprintf("./sim/data/tau_sim/%i_%i_%i_tau_true.RDS", subjects, rois, volumes))


#Summary of difference between Omega_k and Omega_0 by subject
#Raw difference in values
#summary(map_dbl(.x = Omega_k, ~mean(.x - Omega_0)))
#summary(map_dbl(.x = Omega_k, ~mean(solve(.x) - Sigma_0)))

#Difference in zeros (should be no more than 2 * n_edge in max/min direction)
#summary(map_dbl(.x = Omega_k, ~mean(sum(abs(.x[upper.tri(.x, diag = FALSE)]) <= 0.001) - sum(Omega_0[upper.tri(Omega_0, diag = FALSE)] == 0))))
#hist(inv_Tau_k, breaks = 10)
#hist(Tau_k, breaks = 10)

  #Return
  return(tibble(data = data_list, tau = Tau_k))
}


#Simulation grid
sim.df <- 
  expand_grid(
    subjects = c(10, 25, 50, 100),
    rois = c(10, 25, 50),
    volumes = c(250, 500),
    prop_tru_con = 1/5,
    n_edges = 1
  )

#Simulate data based on grid
sim_res.df <- 
  sim.df %>%
    mutate(
      sim = map(.x = 1:n(), ~sim_data(subjects[.x], rois[.x], volumes[.x], prop_tru_con[.x], n_edges[.x]))
    )

#Write out for later
write_rds(sim_res.df, "./sim/data/sim_res.RDS")

#Read in results for testing
sim_res.df <- read_rds("./sim/data/sim_res.RDS")

#Set sampling params
y <- sim_res.df$sim[[1]]$data
subjects  <- length(y)
n_samples <- 500
n_burn    <- 0
n_updates <- 0
alpha_tau <- subjects #Fixed shape parameter for inv_tauk
lambda_2  <- subjects / 50 #Keep mean at 50

#Set lambda 1-3 penalty gamma a, b  hyperparams
alpha <- c(1, 1, 1) #1 - G_k, 2 - tau_k, 3 - Omega_0 (Bayesian sparse covariance/glasso)
beta  <- c(1, 1, 1)

#Sample sizes
K  <- length(y)
p  <- ncol(y[[1]])
Sk <- map(.x = y, ~t(.x) %*% .x)

#Initialize Omega_k estimates for initial values
#Omega_k - tune lambda by bic and then grab G and Omega_0
lambda_grid <- 10^seq(-3, 0, length.out = 10)
bic         <- vector(mode = "numeric", length = 0)

#Tune lambda for independent glasso
for (lam in lambda_grid) {
  omega_k <- ind_graphs(y, 0.1)
  bic   <- c(bic, bic_cal(y, omega_k))
}
#Compute initial est. via best bic
omega_k <- ind_graphs(y, lambda_grid[which.min(bic)])
adj_k   <- map(.x = omega_k, ~abs(.x) > 0.001)
#Omega0  <- apply(abind::abind(omega_k, along=3),1:2,mean)
omega_0 <- Reduce("+", omega_k) / K

#Loop through just to make sure PD
for (k in 1:K) {
  if (any(eigen(omega_k[[k]])$values < 0)) {
    warning(str_c("omega_", k, "was not positive definite in sim."))
    omega_k[[k]] <- 
      omega_k[[k]] |>
      (\(x) {as.matrix(Matrix::nearPD(x)$mat)})()
  }
}

#Initialize estimates for tau_k
#Tau vector of subject specific regularization param on Omega_0
tau_vec <- vector(mode = "numeric", length = K)
#Set tau's MH stepsize
step_tau  <- rep(0.5, K)

#Iterate over each subject, find optimal tau_k based on posterior
# tibble(
#   tau_grid = 1:500,
#   tau_res  = map_dbl(.x = tau_grid, ~log_tau_posterior(.x, omega_k[[1]], omega_0, alpha_tau  = alpha_tau, lambda_2 = 0))
# ) %>%
#   ggplot(aes(x = tau_grid, y = tau_res)) +
#   geom_point() +
#   geom_line() +
#   theme_minimal() -> tau_plot

for (k in 1:K) {
  #print(k)
  #Tau posterior for fixed omega_k, omega_0, and lambda_2 = 0
  f_opt <- function(tau) {
    -1 * log_tau_posterior(tau, omega_k[[k]], omega_0, alpha_tau = alpha_tau, lambda_2 = 0)
  }
  #Optimize in 1D
  tau_vec[k] <- c(optimize(f_opt, lower = 1, upper = 200, tol = 0.01)$min)
}



