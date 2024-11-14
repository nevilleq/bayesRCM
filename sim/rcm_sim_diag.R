#Load the good stuff
library(glasso) #individual glasso
library(JGL)    #Joint fused/group glasso
library(rcm)  
library(bayesRCM)
library(tidyverse)
library(gt)

#Set wd
#setwd("/panfs/jay/groups/1/zhangl4/nevil066/bayes/bayesRCM")
#setwd("~/dissertation/bayesRCM/")

#Prelims and functions
source("./sim/sim_funcs.R")

###############################################################################
#Read in results (freq & bayes)
#file_path to read in results
in_path <- "./sim/sim_res/freq/model_results/"
in_list <- list.files(in_path, pattern = "setting")

#Read in sim settings
in_grid <- list.files("./sim/sim_data/", pattern = "grid")
sim_settings.df <- read_rds(str_c("./sim/sim_data/", in_grid, sep = "/"))

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
  filter(model == "rcm") %>% #Filter which model here
  left_join(., sim_settings.df, by = c("setting", "seed" = "seed_k")) %>%
  mutate(
    result = map(.x = in_files, ~read_rds(.x)),
  ) %>%
  unnest(result)

#Take out lambdas
lam.df <- sim_res.df %>% filter(model == "rcmlam")

#Pivot wider
sim_res.df <- 
  sim_res.df %>% 
  filter(model != "lam", model != "rcmlam") %>%
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
N <- nrow(sim_res.df)
print("Model to process: ")
sim_res.df %>% pull(model) %>% unique() %>% print()

#Loop through and grab model diagnostics (Diff Norms, MCC, Accuracy, etc.)
for(n in 1:N) {
  
  #Display iteration
  #print(n)
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