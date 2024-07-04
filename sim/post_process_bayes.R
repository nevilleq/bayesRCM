#Load the good stuff
library(bayesRCM)
library(tidyverse)
library(glasso)
library(gt)

# #Set working directory to bayesRCM
setwd("/panfs/jay/groups/1/zhangl4/nevil066/bayes/bayesRCM")

#Prelims and functions
source("./sim/sim_funcs.R")

#file_path to read in data
in_path <- "./sim/sim_data/"

#Read in sim settings
in_grid <- list.files("./sim/sim_data/", pattern = "grid")
sim_settings.df <- read_rds(str_c(in_path, in_grid, sep = "/"))

threshold <- 0.001

#####################################################################################
#2. Read in results and true params
#file_path to read in results
in_path <- "./sim/sim_res/bayes/post_results/"
in_list <- list.files(in_path, pattern = "setting")

print("Read results")

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
true_param.df <- read_rds("./sim/sim_res/true_param_df.rds")

#Join by setting and seed
sim_res.df <-
    left_join(sim_res.df, true_param.df, by = c("setting", "seed"))

#####################################################################################
#3. Model diagnostics  
n_start <- 1
N <- nrow(sim_res.df)
n_start <- 568
n_end <- N

for(n in n_start:n_end) {

#Print
print("Processing file:")
print(sprintf("%i/%i", n, N))

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
    write_rds(sim_diag.df, str_c("./sim/sim_res/bayes/sim_results/", name, ".rds"))

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
