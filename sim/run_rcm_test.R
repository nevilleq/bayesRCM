#Load the good stuff
library(bayesRCM)
library(tidyverse)
library(glasso)
library(gt)

# #Set working directory to bayesRCM
#setwd("/panfs/jay/groups/1/zhangl4/nevil066/bayes/bayesRCM")
setwd("~/dissertation/bayesRCM/")

#Prelims and functions
source("./sim/sim_funcs.R")

#file_path to read in data
in_path <- "./sim/sim_data/"
in_list <- list.files("./sim/sim_data/", pattern = "setting")

#Read in sim settings
in_grid <- list.files("./sim/sim_data/", pattern = "grid")
sim_settings.df <- read_rds(str_c(in_path, in_grid, sep = "/"))

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

#For testing
sim_data.df$data_list[[1]] -> y

#######################################################################################
#1. Model Fits

#Settings for sim results
N <- nrow(sim_data.df)
threshold <- 0.001
#n_cores   <- parallel::detectCores() - 1
#n_cores   <- 7
n_cores   <- 50

#Storage lists
#res_list_bayes <- vector(mode = "list", length = N)
n <- 1



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
  n_samps <- 10
  n_burn <- 0
  n_updates <- 0
  

  #Fit model (5000 post samps, 1000 burn)
  # n_samples = 100; n_burn = 10; n_cores = n_cores; n_updates = 5;
  # result <- bayesRCM::rcm(y, n_samples = 100, n_burn = 10, n_cores = n_cores, n_updates = 1)
    result <- rcm(y, n_samples = n_samps, n_burn = n_burn, n_cores = n_cores, n_updates = n_updates)
    print("Minutes ellapsed...")
    print(sum(result$timer)/60)
    