#Load the good stuff
library(rcm)  
library(tidyverse)
library(gt)

#Set wd
#setwd("/panfs/jay/groups/1/zhangl4/nevil066/bayes/bayesRCM")
setwd("~/dissertation/bayesRCM/")

#Prelims and functions
source("./sim/sim_funcs.R")

###############################################################################
#Read in results (freq & bayes)
#Read in frequentist
print("Reading freq sims....")
in_path  <- "./sim/sim_res/freq/model_results/"
in_files <- list.files(in_path)
in_files <- in_files[str_detect(in_files, "lam")]

#Read in each result
rcm_res_lam.df <-
  tibble(
    in_path = str_c(in_path, in_files),
    result  = map(.x = in_path, ~as.data.frame(read_rds(.x)))
  ) %>%
  unnest(result) %>%
  dplyr::select(-in_path)

#Print
print(head(rcm_res_lam.df))