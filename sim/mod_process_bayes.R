#Load the good stuff
library(bayesRCM)
library(tidyverse)
library(glasso)
library(gt)

# #Set working directory to bayesRCM
setwd("/panfs/jay/groups/1/zhangl4/nevil066/bayes/bayesRCM")

#Prelims and functions
source("./sim/sim_funcs.R")

#######################################################################################
#1. Model Fits

#Settings for sim results
#N <- nrow(sim_data.df)
threshold <- 0.001
#n_cores   <- parallel::detectCores() - 1
#n_cores   <- 7
n_cores   <- 50

#Storage lists
#res_list_bayes <- vector(mode = "list", length = N)
n <- 1
write <- TRUE
#N <- 2

in_res <- list.files("./sim/sim_res/bayes/model_results/")
in_files <- str_c("./sim/sim_res/bayes/model_results/", in_res)
names <- str_remove(in_res, ".rds")

#in_res[1]
#in_files[1]
#names[1]
            
print(length(in_files))
N <- length(in_files)

#stop("Pause")
#Loop to fit bayesRCM::rcm
for(n in 1:N) {

#Parallelized loop
#n_cores <- parallel::detectCores() - 1
# cl <- parallel::makeCluster(n_cores)
# on.exit(parallel::stopCluster(cl)) #When done stop cluster
# doParallel::registerDoParallel(cl) #Initialize clusters
# `%dopar%` <- foreach::`%dopar%`

#Loop
# fit_mods <- foreach::foreach(
#   n         = 1:N, #Subject index
#   .packages = c("bayesRCM", "glasso", "stringr") #Packages
#   # .noexport = c("graph_update","gwish_ij_update", "rgwish") #Functions necessary to export
# ) %dopar% {
  print(n)
  #Unique name for directory to store figs
  
  #Grab data & dimensions
  #y <- data_list <- sim_data.df$data_list[[n]] #Subject data list
  #p <- ncol(y[[1]]) #Dimension of data rois
  #K <- length(y)    #Subjects
  
  #Fit model (5000 post samps, 1000 burn)
  # n_samples = 100; n_burn = 10; n_cores = n_cores; n_updates = 5;
  # result <- bayesRCM::rcm(y, n_samples = 100, n_burn = 10, n_cores = n_cores, n_updates = 1)
    #result <- rcm(y, n_samples = 4000, n_burn = 1000, n_cores = n_cores, n_updates = 10)
    #write_rds(result, str_c("./sim/sim_res/bayes/model_results/", name, ".rds"))
    result <- read_rds(in_files[n])
    n_iter <- ncol(result$omega_0)
    p <- 10
    
    #Grab posterior median and mean for Omega_K
    omegak.df <-
        map_df(.x = 1:n_iter,
                      ~tibble(iteration = .x, 
                                      omega_k   = list(as.data.frame(result$omega_k[,,.x]) %>% mutate(roi = 1:n())))
        ) %>%
        unnest(cols = c(omega_k)) %>%
        rename_with(
            .cols = -iteration,
            ~str_replace(.x, "V", "Sub. ")
        ) %>%
        pivot_longer(
            cols = contains("Sub"),
            names_to = "subject",
            values_to = "value"
        ) %>%
        group_by(subject, roi) %>%
        summarise(
            post_mean = mean(value),
            post_median = median(value),
            .groups = "drop"
        ) %>%
        nest(result = -c(subject)) %>%
        mutate(
            mean_mat = map(.x = result, ~fill_mat(.x$post_mean, p)),
            med_mat  = map(.x = result, ~fill_mat(.x$post_median, p))
        )
    
    #Save for later (1 = L1 loss/posterior median; 2 = L2 loss/posterior mean)
    omega_k2 <- omegak.df$mean_mat
    adj_k2   <- map(.x = omega_k2, ~abs(.x) > threshold)
    omega_k1 <- omegak.df$med_mat
    adj_k1   <- map(.x = omega_k1, ~abs(.x) > threshold)
    
    #Posterior mean and median omega_0 matrices
    omega_02 <- apply(result$omega_0, 1, mean) %>% fill_mat(., p)
    omega_01 <- apply(result$omega_0, 1, median) %>% fill_mat(., p)
    
    res_list_bayes <- list(omega_k1 = omega_k1, adj_k1 = adj_k1, omega_01 = omega_01,
                                                  omega_k2 = omega_k2, adj_k2 = adj_k2, omega_02 = omega_02)
    
    #Store bayesRCM results
    if(write){
        readr::write_rds(res_list_bayes, str_c("./sim/sim_res/bayes/post_results/", names[n],".rds"))
    }
} #End loop

