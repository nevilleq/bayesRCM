#Load the good stuff
library(bayesRCM)
library(tidyverse)
library(glasso)
library(gt)

# #Set working directory to bayesRCM
setwd("/panfs/jay/groups/1/zhangl4/nevil066/bayes/bayesRCM")

#Prelims and functions
source("./sim_test/sim_funcs.R")

#############################################################################
#Read in sim data

#create storage directories
dir.create("./sim_test/sim_res/", showWarnings = FALSE)
dir.create("./sim_test/sim_res/bayes/", showWarnings = FALSE)
dir.create("./sim_test/sim_res/bayes/model_results", showWarnings = FALSE)
dir.create("./sim_test/sim_res/bayes/sim_results", showWarnings = FALSE)
dir.create("./sim_test/sim_res/bayes/post_results", showWarnings = FALSE)



#file_path to read in data
in_path <- "./sim_test/sim_data/"
in_list <- list.files("./sim_test/sim_data/", pattern = "setting")

#Read in sim settings
in_grid <- list.files("./sim_test/sim_data/", pattern = "grid")
sim_settings.df <- read_rds(str_c(in_path, in_grid, sep = "/"))

print("Data read")

#Read in all sim data
sim_data.df <-
  tibble(
    in_path  = in_path,
    in_list  = in_list,
    setting  = str_split_fixed(in_list, "_", 2)[,1] %>% parse_number(),
    seed     = str_split_fixed(in_list, "_", 2)[,2] %>% parse_number(),
    in_files = str_c(in_path, in_list, sep = "/")
  ) %>%
  dplyr::select(setting, seed, in_files) %>%
  left_join(., sim_settings.df, by = c("setting", "seed" = "seed_k")) %>%
  mutate(
    sim_res     = map(.x = in_files, ~read_rds(.x)),
    data_list   = map(sim_res, "data_list"),
    true_params = map(sim_res, "true_params")
  ) %>%
  dplyr::select(-sim_res) #%>%
  #group_by(setting) %>%
  #slice(1) %>%
  #ungroup() #Temporary to test just one sim at each setting

#Write out true params to sim_res
sim_data.df %>%
  dplyr::select(setting, seed, true_params) %>%
  write_rds(., "./sim_test/sim_res/true_param_df.rds")

#For testing
sim_data.df$data_list[[1]] -> y

#######################################################################################
#1. Model Fits

#Settings for sim results
N <- nrow(sim_data.df)
threshold <- 0.001
#n_cores   <- parallel::detectCores() - 1
#n_cores   <- 7
n_cores   <- 20

#Storage lists
#res_list_bayes <- vector(mode = "list", length = N)
#n <- which(sim_data.df$setting == 3 & sim_data.df$seed == 302)
n <- 1
write <- TRUE
#N <- 3

#print("Starting 3 iterations")

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
  name <- 
    with(
      sim_data.df,
      paste0("setting", setting[n],
             "_seed", seed[n])
    )
  
  #Grab data & dimensions
  y <- data_list <- sim_data.df$data_list[[n]] #Subject data list
  p <- ncol(y[[1]]) #Dimension of data rois
  K <- length(y)    #Subjects
  
  #Fit model (5000 post samps, 1000 burn)
  # n_samples = 100; n_burn = 10; n_cores = n_cores; n_updates = 5;
  # result <- bayesRCM::rcm(y, n_samples = 100, n_burn = 10, n_cores = n_cores, n_updates = 1)
   # result <- rcm(y, n_samples = 4000, n_burn = 1000, n_cores = n_cores, n_updates = 10)
   # write_rds(result, str_c("./sim_test/sim_res/bayes/model_results/", name, ".rds"))
    # result <- read_rds("./sim_test/sim_res/bayes/model_results/setting1_seed1.rds")
    result <- tryCatch({
      #if (n != 2){
      #  bayesRCM::rcm(y, n_samples = 1000, n_burn = 100, n_cores = n_cores, n_updates = 10)
      #} else {
      #  1 = 2
      #}
      rcm(y, n_samples = 4000, n_burn = 1000, n_cores = n_cores, n_updates = 10)
      #return(result)
      }, error = function(e){
        cat("Error in bayesRCM::rcm(), likely in subject updates", "\n", conditionMessage(e), "\n")
      })
       #result <- bayesRCM::rcm(y, n_samples = 100, n_burn = 10, n_cores = n_cores, n_updates = 1)
      # result <- rcm(y, n_samples = 4000, n_burn = 1000, n_cores = n_cores, n_updates = 10)
        
        #Write out model result
        write_rds(result, str_c("./sim_test/sim_res/bayes/model_results/", name, ".rds"))
        # result <- read_rds("./sim_test/sim_res/bayes/model_results/setting1_seed1.rds")
        
        if(is.null(result)) {
          res_list_bayes <- list(omega_k1 = NULL, adj_k1 = NULL, omega_01 = NULL,
                                 omega_k2 = NULL, adj_k2 = NULL, omega_02 = NULL)
          if(write){
            readr::write_rds(res_list_bayes, str_c("./sim_test/sim_res/bayes/post_results/", name,".rds"))
          }
          next()
    } else {
    n_iter <- ncol(result$omega_0)
    
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
        readr::write_rds(res_list_bayes, str_c("./sim_test/sim_res/bayes/post_results/", name,".rds"))
    }
  }
} #End loop