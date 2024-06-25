#Load the good stuff
library(bayesRCM)
library(tidyverse)
library(glasso)
library(gt)

# #Set working directory to bayesRCM
setwd("/panfs/jay/groups/1/zhangl4/nevil066/bayes/bayesRCM")

#Prelims and functions
source("./sim/sim_funcs.R")

#############################################################################
#Read in sim data

#create storage directories
dir.create("./sim/sim_res/", showWarnings = FALSE)
dir.create("./sim/sim_res/bayes/", showWarnings = FALSE)
dir.create("./sim/sim_res/bayes/model_results", showWarnings = FALSE)
dir.create("./sim/sim_res/bayes/sim_results", showWarnings = FALSE)
dir.create("./sim/sim_res/bayes/post_results", showWarnings = FALSE)



#file_path to read in data
in_path <- "./sim/sim_data/"
in_list <- list.files("./sim/sim_data/", pattern = "setting")

#Read in sim settings
in_grid <- list.files("./sim/sim_data/", pattern = "grid")
sim_settings.df <- read_rds(str_c(in_path, in_grid, sep = "/"))

print("reading sims")

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

print("write true param df")  

#Write out true params to sim_res
sim_data.df %>%
  dplyr::select(setting, seed, true_params) %>%
  write_rds(., "./sim/sim_res/true_param_df.rds")
