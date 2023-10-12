#Load the good stuff
library(glasso) #individual glasso
library(JGL)    #Joint fused/group glasso
library(rcm)  
library(bayesRCM)
library(tidyverse)
library(gt)

#Prelims and functions
source("./sim/sim_funcs.R")

###############################################################################
#Read in results (freq & bayes)
#Read in frequentist
in_path  <- "./sim/sim_res/freq/sim_results/"
in_files <- list.files(in_path)

#Read in each result
sim_res_freq.df <-
  tibble(
    in_path = str_c(in_path, in_files),
    result  = map(.x = in_path, ~read_rds(.x))
  ) %>%
  unnest(result) %>%
  dplyr::select(-in_path)

#Read in bayesian results
in_path  <- "./sim/sim_res/bayes/sim_results/"
in_files <- list.files(in_path)

#Read in each result
sim_res_bayes.df <-
  tibble(
    in_path = str_c(in_path, in_files),
    result  = map(.x = in_path, ~read_rds(.x))
  ) %>%
  unnest(result) %>%
  mutate(model = "bayesRCM") %>%
  dplyr::select(model, everything(), -c(in_path)) %>%
  rename(norm_0 = norm_01, norm_k = norm_k1, k_diag = k1_diag, O_diag = O1_diag) %>% #Use L1 loss/post median
  dplyr::select(-c(norm_02, norm_k2, k2_diag, O2_diag))

#Final result data frame
sim_res.df <-
  bind_rows(
    sim_res_freq.df,
    sim_res_bayes.df
  ) %>%
  arrange(model, setting)

#If no bayes results just use this
#sim_res.df <- sim_res_freq.df

#Display sim settings
sim_res.df %>%
  dplyr::select(setting, subjects:lambda_2) %>%
  distinct() %>%
  gt() %>%
  tab_header("Simulation Settings")

#Omega_0 norm results
sim_res.df %>%
  dplyr::select(model:seed, norm_0) %>%
  unnest(norm_0) %>%
  dplyr::select(-seed) %>%
  mutate(setting = str_c("Setting ", setting)) %>%
  group_by(setting) %>%
  arrange(rev(desc(setting)), rev(desc(L1))) %>%
  gt() %>%
  tab_header("Omega_0 Difference Norms by Setting")

#Omega_0 diagnostics  
sim_res.df %>%
  dplyr::select(model:seed, O_diag) %>%
  unnest(O_diag) %>%
  dplyr::select(-seed) %>%
  mutate(
    setting_str = as.factor(str_c("Setting ", setting)) %>%
                  fct_reorder(setting, .desc = FALSE)
  ) %>%
  dplyr::select(model, setting_str, mcc, accuracy, kappa) %>%
  group_by(setting_str) %>%
  arrange(setting_str, desc(mcc)) %>%
  gt() %>%
  tab_header("Omega_0 Diagnostics by Setting")

#Omega_k norm results  
omega_k_norm.gt <-
  sim_res.df %>%
  dplyr::select(model:seed, norm_k) %>%
  unnest(norm_k) %>%
  mutate(setting = str_c("Setting ", setting)) %>%
  group_by(subject, setting) %>%
  gt() %>%
  tab_header("Omega_K Difference Norm by Setting & Subject")
omega_k_norm.gt

omega_k_norm_sum.gt <-
  sim_res.df %>%
  dplyr::select(model:seed, norm_k) %>%
  unnest(norm_k) %>%
  mutate(setting = str_c("Setting ", setting)) %>%
  dplyr::select(-seed) %>%
  group_by(model, setting) %>%
  summarise(
    across(
      .cols = where(is.numeric),
      .fns  = list(mean = mean, sd = sd),
      .names = "{.col}.{.fn}"
    ),
    .groups = "drop"
  ) %>%
  group_by(setting) %>%
  gt() %>%
  tab_header("Omega_K Diff-Norm Summarised over Subjects by Setting")

omega_k_norm_sum.gt

#Omega_k diagnostics
omega_k_diag.gt <-
  sim_res.df %>%
  dplyr::select(model:seed, k_diag) %>%
  unnest(k_diag) %>%
  dplyr::select(-seed) %>%
  mutate(
    setting_str = as.factor(str_c("Setting ", setting)) %>%
      fct_reorder(setting, .desc = FALSE)
  ) %>%
  dplyr::select(model, subject, setting_str, mcc, accuracy, kappa) %>%
  group_by(setting_str, subject) %>%
  arrange(setting_str, subject, desc(mcc)) %>%
  gt() %>%
  tab_header("Omega_K Diagnostics by Setting & Subject")
omega_k_diag.gt

omega_k_diag_sum.gt <-
  sim_res.df %>%
  dplyr::select(model:seed, k_diag) %>%
  unnest(k_diag) %>%
  dplyr::select(-seed) %>%
  mutate(
    setting_str = as.factor(str_c("Setting ", setting)) %>%
      fct_reorder(setting, .desc = FALSE)
  ) %>%
  dplyr::select(model, setting_str, mcc, accuracy, kappa) %>%
  group_by(model, setting_str) %>%
  summarise(
    across(
      .cols = where(is.numeric),
      .fns  = list(mean = mean, sd = sd),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  ) %>%
  group_by(setting_str) %>%
  arrange(setting_str, desc(accuracy_mean)) %>%
  gt() %>%
  tab_header("Omega_K Diagnostics Summarised over Subjects by Setting")

omega_k_diag_sum.gt
