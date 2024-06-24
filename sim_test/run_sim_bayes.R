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
  group_by(setting) %>%
  slice(1:10) %>%
  ungroup() %>%
  left_join(., sim_settings.df, by = c("setting", "seed" = "seed_k")) %>%
  mutate(
    sim_res     = map(.x = in_files, ~read_rds(.x)),
    data_list   = map(sim_res, "data_list"),
    true_params = map(sim_res, "true_params")
  ) %>%
  dplyr::select(-sim_res) #%>%
 # group_by(setting) %>%
 # slice(10) %>%
 # ungroup() #Temporary to test just one sim at each setting

 #Print number of sims to run
 print("Number of sims: ")
 print(nrow(sim_data.df))
 
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
n_cores   <- parallel::detectCores() - 1
print(n_cores)
#n_cores   <- 7
#n_cores   <- 50

#Storage lists
#res_list_bayes <- vector(mode = "list", length = N)
n <- 1
write <- TRUE
#N <- 2

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
   #result <- bayesRCM::rcm(y, n_samples = 100, n_burn = 10, n_cores = n_cores, n_updates = 1)
    result <- rcm(y, n_samples = 4000, n_burn = 1000, n_cores = n_cores, n_updates = 10)
    write_rds(result, str_c("./sim_test/sim_res/bayes/model_results/", name, ".rds"))
    # result <- read_rds("./sim/sim_res/bayes/model_results/setting1_seed1.rds")
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
} #End loop

#####################################################################################
#2. Read in results and true params
#file_path to read in results
in_path <- "./sim_test/sim_res/bayes/post_results/"
in_list <- list.files(in_path, pattern = "setting")

#Read in results
sim_res.df <-
    tibble(
        in_path  = in_path,
        in_list  = in_list,
        setting  = str_split_fixed(in_list, "_", 3)[ ,1] %>% parse_number(),
        seed     = str_split_fixed(in_list, "_", 3)[ ,2] %>% parse_number(),
        in_files = str_c(in_path, in_list, sep = "/")
    ) %>%
    dplyr::select(setting, seed, in_files) %>%
    left_join(., sim_settings.df, by = c("setting", "seed" = "seed_k")) %>%
    mutate(
        result = map(.x = in_files, ~read_rds(.x)),
    ) %>%
    unnest(result)

#Pivot wider
sim_res.df <- 
    sim_res.df %>% 
    mutate(
        name = names(result)
    ) %>%
    pivot_wider(
        names_from = "name",
        values_from = "result"
    )

#Read in true params
true_param.df <- read_rds("./sim_test/sim_res/true_param_df.rds")

#Join by setting and seed
sim_res.df <-
    left_join(sim_res.df, true_param.df, by = c("setting", "seed"))

#####################################################################################
#3. Model diagnostics  
N <- nrow(sim_res.df)


for(n in 1:N) {

#Unique name for directory to store figs
    name <- 
        with(
            sim_res.df,
            paste0("setting", setting[n],
                          "_seed", seed[n])
        )

#True parameters
true_params <- sim_res.df$true_params[[n]]
Omega_0 <- true_params$omega_0
Adj_0   <- abs(Omega_0) > threshold
Omega_k <- true_params$omega_k
Adj_k   <- map(.x = true_params$omega_k, ~abs(.x) > threshold)

#Diagnostic .df
sim_diag.df <-
    sim_res.df[n, ] %>%
    mutate(
        #L1 Loss (posterior median)
        diff_k1 = map(.x = omega_k1, ~get_diff(.x, Omega_k)),
        diff_01 = map(.x = omega_01, ~Omega_0 - .x),
        norm_01 = map(.x = diff_01, ~tibble(
            L1        = norm(.x, "1"),
            Frobenius = norm(.x, "F"),
            Spectral  = norm(.x, "2")
        )),
        norm_k1 = map(.x = diff_k1, ~get_norms(.x)),
        k1_diag = map(.x = adj_k1, ~get_bin_diag_k(.x, Adj_k)),
        adj_01  = map(.x = omega_01, ~abs(.x) > threshold),
        O1_diag = map(.x = adj_01, ~get_bin_diag_0(.x, Adj_0)),
        
        #L2 Loss (posterior mean)
        diff_k2 = map(.x = omega_k2, ~get_diff(.x, Omega_k)),
        diff_02 = map(.x = omega_02, ~Omega_0 - .x),
        norm_02 = map(.x = diff_02, ~tibble(
            L1        = norm(.x, "1"),
            Frobenius = norm(.x, "F"),
            Spectral  = norm(.x, "2")
        )),
        norm_k2 = map(.x = diff_k2, ~get_norms(.x)),
        k2_diag = map(.x = adj_k2, ~get_bin_diag_k(.x, Adj_k)),
        adj_02  = map(.x = omega_02, ~abs(.x) > threshold),
        O2_diag = map(.x = adj_02, ~get_bin_diag_0(.x, Adj_0))
    ) %>%
dplyr::select(setting, seed, subjects:n_flip, 
                            all_of(starts_with("norm")), all_of(ends_with("diag")))

    #Write out result per iteration
    write_rds(sim_diag.df, str_c("./sim_test/sim_res/bayes/sim_results/", name, ".rds"))

}

###############################################################################
#Visualize diagnostics
# in_path  <- "./sim/sim_res/bayes/sim_results/"
# in_files <- list.files(in_path)
# 
# #Read in each result
# sim_res.df <-
#   tibble(
#     in_path = str_c(in_path, in_files),
#     result  = map(.x = in_path, ~read_rds(.x))
#   ) %>%
#   unnest(result) %>%
#   dplyr::select(-in_path)
# 
# #Omega_0 norm results L1/median
# sim_res.df %>%
#   dplyr::select(setting, seed, norm_01) %>%
#   unnest(norm_01)
# 
# #Omega_0 norm results L2/mean
# sim_res.df %>%
#   dplyr::select(setting, seed, norm_02) %>%
#   unnest(norm_02)
# 
# #Omega_0 diagnostics / L1 post loss median  
# sim_res.df %>%
#   dplyr::select(setting, seed, O1_diag) %>%
#   unnest(O1_diag)
# 
# #Omega_k norm results / L1 post loss median
# sim_res.df %>%
#   dplyr::select(setting, seed, norm_k1) %>%
#   unnest(norm_k1)
# 
# #Omega_k norm results / L2 post loss mean
# sim_res.df %>%
#   dplyr::select(setting, seed, norm_k2) %>%
#   unnest(norm_k2)
# 
# #Omega_k diagnostics L1 post loss / median
# sim_res.df %>%
#   dplyr::select(setting, seed, k1_diag) %>%
#   unnest(k1_diag)
# 
# #Omega_k diagnostics L2 post loss / mean
# sim_res.df %>%
#   dplyr::select(setting, seed, k2_diag) %>%
#   unnest(k2_diag)