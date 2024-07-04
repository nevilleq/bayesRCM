#Load the good stuff
library(glasso) #individual glasso
library(JGL)    #Joint fused/group glasso
library(rcm)  
library(bayesRCM)
library(tidyverse)
library(gt)

#WD
setwd("/panfs/jay/groups/1/zhangl4/nevil066/bayes/bayesRCM")

#Prelims and functions
source("./sim/sim_funcs.R")

#############################################################################
#Read in sim data

#create storage directories
dir.create("./sim/sim_res/", showWarnings = FALSE)
dir.create("./sim/sim_res/freq/", showWarnings = FALSE)
dir.create("./sim/sim_res/freq/model_results/", showWarnings = FALSE)
dir.create("./sim/sim_res/freq/sim_results/", showWarnings = FALSE)

#file_path to read in data
in_path <- "./sim/sim_data/"
in_list <- list.files("./sim/sim_data/", pattern = "setting")

#Read in sim settings
in_grid <- list.files("./sim/sim_data/", pattern = "grid")
sim_settings.df <- read_rds(str_c(in_path, in_grid, sep = "/"))
print("Data read....")
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
  dplyr::select(-sim_res)

#Write out true params to sim_res
sim_data.df %>%
  dplyr::select(setting, seed, true_params) %>%
  write_rds(., "./sim/sim_res/true_param_df.rds")

#Take only 1 sim per setting for now
# sim_data.df <-
#   sim_data.df %>%
#   group_by(setting) %>%
#   slice(1) %>%
#   ungroup() #Temporary to test just one sim at each setting

#####################################################################################
#2. Read in results and true params
#file_path to read in results
in_path <- "./sim/sim_res/freq/model_results/"
in_list <- list.files(in_path, pattern = "setting")

#Read in results
sim_res.df <-
  tibble(
    in_path  = in_path,
    in_list  = in_list,
    model    = str_split_fixed(in_list, "_", 3)[ ,1],
    setting  = str_split_fixed(in_list, "_", 3)[ ,2] %>% parse_number(),
    seed     = str_split_fixed(in_list, "_", 3)[ ,3] %>% parse_number(),
    in_files = str_c(in_path, in_list, sep = "/")
  ) %>%
  dplyr::select(model, setting, seed, in_files) %>%
  left_join(., sim_settings.df, by = c("setting", "seed" = "seed_k")) %>%
  mutate(
    result = map(.x = in_files, ~read_rds(.x)),
  ) %>%
  unnest(result)

#Take out lambdas
lam.df <- sim_res.df %>% filter(model == "lam")

#Pivot wider
sim_res.df <- 
  sim_res.df %>% 
  filter(model != "lam") %>%
  mutate(
    name = names(result)
  ) %>%
  pivot_wider(
    names_from = "name",
    values_from = "result"
  )

#Read in true params
true_param.df <- read_rds("./sim/sim_res/true_param_df.rds")

#Join by setting and seed
sim_res.df <-
  left_join(sim_res.df, true_param.df, by = c("setting", "seed"))

#####################################################################################
#3. Model diagnostics  
n_start <- 1
N <- nrow(sim_res.df)
n_start <- 2152
n_end <- N

print("Starting model diagnostic processing...")

for(n in n_start:n_end) {
  
  #Print
  print("Processing file:")
  print(sprintf("%i/%i", n, N))
  
  #Unique name for directory to store figs
  name <- 
    with(
      sim_res.df,
      paste0(model[n], "_",
             "setting", setting[n],
             "_seed", seed[n])
    )

#True parameters
threshold <- 0.001
true_params <- sim_res.df$true_params[[n]]
Omega_0 <- true_params$omega_0
Adj_0   <- abs(Omega_0) > threshold
Omega_k <- true_params$omega_k
Adj_k   <- map(.x = Omega_k, ~abs(.x) > threshold)

#Diagnostic .df
sim_diag.df <-
  sim_res.df[n, ] %>%
  mutate(
    diff_k = map(.x = omega_k, ~get_diff(.x, Omega_k)),
    diff_0 = map(.x = omega_0, ~Omega_0 - .x),
    norm_0 = map(.x = diff_0, ~tibble(
      L1        = norm(.x, "1"),
      Frobenius = norm(.x, "F"),
      Spectral  = norm(.x, "2")
    )),
    norm_k = map(.x = diff_k, ~get_norms(.x)),
    k_diag = map(.x = adj_k, ~get_bin_diag_k(.x, Adj_k)),
    adj_0  = map(.x = omega_0, ~abs(.x) > threshold),
    O_diag = map(.x = adj_0, ~get_bin_diag_0(.x, Adj_0))
  ) %>%
  dplyr::select(model, setting, seed, subjects:n_flip, 
                all_of(starts_with("norm")), all_of(ends_with("diag")))

  #Write out result per iteration
  write_rds(sim_diag.df, str_c("./sim/sim_res/freq/sim_results/", name, ".rds"))

} #End loop
#############################################################################
#Visualize final results
#Read in frequentist
# in_path  <- "./sim/sim_res/freq/sim_results/"
# in_files <- list.files(in_path)
# 
# #Read in each result
# sim_res_freq.df <-
#   tibble(
#     in_path = str_c(in_path, in_files),
#     result  = map(.x = in_path, ~read_rds(.x))
#   ) %>%
#   unnest(result) %>%
#   dplyr::select(-in_path)

# #Read in bayesian results
# in_path  <- "./sim/sim_res/bayes/sim_results/"
# in_files <- list.files(in_path)
# 
# #Read in each result
# sim_res_bayes.df <-
#   tibble(
#     in_path = str_c(in_path, in_files),
#     result  = map(.x = in_path, ~read_rds(.x))
#   ) %>%
#   unnest(result) %>%
#   mutate(model = "bayesRCM") %>%
#   dplyr::select(model, everything(), -c(in_path))
# 
# #Final result data frame
# sim_res.df <-
#   bind_rows(
#     sim_res_freq.df,
#     sim_res_bayes.df
#   ) %>%
#   arrange(model, setting)

# #If no bayes results just use this
# sim_res.df <- sim_res_freq.df
# 
# #Display sim settings
# sim_res.df %>%
#   dplyr::select(setting, subjects:lambda_2) %>%
#   distinct() %>%
#   gt() %>%
#   tab_header("Simulation Settings")

# #Omega_0 norm results
# sim_res.df %>%
#   dplyr::select(model:seed, norm_0) %>%
#   unnest(norm_0) %>%
#   group_by(setting) %>%
#   gt() %>%
#   tab_header("Omega_0 Difference Norms by Setting")
# 
# #Omega_0 diagnostics  
# sim_res.df %>%
#   dplyr::select(model:seed, O_diag) %>%
#   unnest(O_diag) %>%
#   group_by(setting) %>%
#   gt() %>%
#   tab_header("Omega_0 Diagnostics by Setting")

# #Omega_k norm results  
# omega_k_norm.gt <-
#   sim_res.df %>%
#   dplyr::select(model:seed, norm_k) %>%
#   unnest(norm_k) %>%
#   mutate(setting = str_c("Setting ", setting)) %>%
#   group_by(subject, setting) %>%
#   gt() %>%
#   tab_header("Omega_K Difference Norm by Setting & Subject")
# omega_k_norm.gt

# omega_k_norm_sum.gt <-
#   sim_res.df %>%
#   dplyr::select(model:seed, norm_k) %>%
#   unnest(norm_k) %>%
#   mutate(setting = str_c("Setting ", setting)) %>%
#   dplyr::select(-seed) %>%
#   group_by(model, setting) %>%
#   summarise(
#     across(
#       .cols = where(is.numeric),
#       .fns  = list(mean = mean, sd = sd),
#       .names = "{.col}.{.fn}"
#     ),
#     .groups = "drop"
#   ) %>%
#   group_by(setting) %>%
#   gt() %>%
#   tab_header("Omega_K Diff-Norm Summarised over Subjects by Setting")
# 
# omega_k_norm_sum.gt

# #Omega_k diagnostics
# omega_k_diag.gt <-
#   sim_res.df %>%
#   dplyr::select(model:seed, k_diag) %>%
#   unnest(k_diag) %>%
#   mutate(setting = str_c("Setting ", setting)) %>%
#   group_by(subject, setting) %>%
#   gt() %>%
#   tab_header("Omega_K Diagnostics by Setting & Subject")

# omega_k_diag_sum.gt <-
#   sim_res.df %>%
#   dplyr::select(model:seed, k_diag) %>%
#   unnest(k_diag) %>%
#   mutate(setting = str_c("Setting ", setting)) %>%
#   dplyr::select(model, setting, mcc, accuracy, kappa) %>%
#   group_by(model, setting) %>%
#   summarise(
#     across(
#       .cols = where(is.numeric),
#       .fns  = list(mean = mean, sd = sd),
#       .names = "{.col}.{.fn}"
#     ),
#     .groups = "drop"
#   ) %>%
#   group_by(setting) %>%
#   gt() %>%
#   tab_header("Omega_K Diagnostics Summarised over Subjects by Setting")
# 
# omega_k_diag_sum.gt
