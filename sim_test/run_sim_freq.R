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
source("./sim_test/sim_funcs.R")

#############################################################################
#Read in sim data

#create storage directories
dir.create("./sim_test/sim_res/", showWarnings = FALSE)
dir.create("./sim_test/sim_res/freq/", showWarnings = FALSE)
dir.create("./sim_test/sim_res/freq/model_results/", showWarnings = FALSE)
dir.create("./sim_test/sim_res/freq/sim_results/", showWarnings = FALSE)

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
  write_rds(., "./sim_test/sim_res/true_param_df.rds")

#Take only 1 sim per setting for now
# sim_data.df <-
#   sim_data.df %>%
#   group_by(setting) %>%
#   slice(1) %>%
#   ungroup() #Temporary to test just one sim at each setting


#######################################################################################
#1. Model Fits

#Settings for sim results
N <- nrow(sim_data.df)
threshold <- 0.001

#Storage lists
res_list_ig <- res_list_jgl <- res_list_rcm <- vector(mode = "list", length = N)
#n <- 1
write <- TRUE

#Loop to compute subject/group networks for each model type
#for(n in 1:N) {
  
#Parallelized loop
#n_cores <- parallel::detectCores() - 1
n_cores <- 50
cl <- parallel::makeCluster(n_cores)
on.exit(parallel::stopCluster(cl)) #When done stop cluster
doParallel::registerDoParallel(cl) #Initialize clusters
`%dopar%` <- foreach::`%dopar%`

#Loop
fit_mods <- foreach::foreach(
  n         = 1:N, #Subject index
  .packages = c("bayesRCM", "glasso", "stringr") #Packages
 # .noexport = c("graph_update","gwish_ij_update", "rgwish") #Functions necessary to export
) %dopar% {
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
#n_iter <- ncol(result$omega_0)
K <- length(y)    #Subjects

## 2.1 Independent glasso  
#Omega_k ind glasso est - tune lambda by bic and then grab G and Omega_0
lambda_grid <- 10^seq(-3, 1, length.out = 40)
bic         <- vector(mode = "numeric", length = length(lambda_grid))

#Tune lambda for independent glasso
for (i in 1:length(lambda_grid)) {
  omega_k <- bayesRCM::ind_graphs(y, lambda_grid[i])
  omega_k <- array(unlist(omega_k), dim = c(p, p, K))
  bic[i]  <- bayesRCM::bic_cal(y, omega_k)
}

#Compute initial est. for omega_k, G_k/adj_k via best mBIC, and Omega_0
lam1_glasso <- lambda_grid[which.min(bic)]
omega_k <- bayesRCM::ind_graphs(y, lam1_glasso)
adj_k   <- purrr::map(.x = omega_k, ~abs(.x) > threshold)
omega_0 <- Reduce("+", omega_k) / K

#Result list
res_list_ig[[n]] <- list(omega_k = omega_k, adj_k = adj_k, omega_0 = omega_0)

#Store ind_glasso results
if(write){
  readr::write_rds(res_list_ig[[n]], str_c("./sim_test/sim_res/freq/model_results/iglasso_", name,".rds"))
}

#######################################
## 2.2 Fused group glasso 
#Tuning
lambda_grid <- 
  tidyr::expand_grid(
    lam1 = lam1_glasso,
    lam2 = 10^seq(-3, 0.5, length.out = 20)
  )
bic <- vector(mode = "numeric", length = nrow(lambda_grid))

#Tune lambda for independent glasso
for (i in 1:nrow(lambda_grid)) {
  omega_k <- 
    JGL::JGL(Y = y, penalty = "group",
             lambda1 = lambda_grid$lam1[i],
             lambda2 = lambda_grid$lam2[i],
             return.whole.theta = TRUE,
             warm = res_list_ig$omega_k)$theta
  omega_k <- array(unlist(omega_k), dim = c(p, p, K))
  bic[i]  <- bayesRCM::bic_cal(y, omega_k)
}

#Grab best tune
lam2_glasso <- as.numeric(lambda_grid[which.min(bic), "lam2"])

#Refit
jgl_out <- JGL::JGL(Y = y, penalty = "group",
                    lambda1 = lam1_glasso, 
                    lambda2 = lam2_glasso,
                    return.whole.theta = TRUE,
                    warm = res_list_ig$omega_k)

#Store and save results
omega_k <- jgl_out$theta
adj_k   <- purrr::map(.x = omega_k, ~abs(.x) > threshold)
omega_0 <- Reduce("+", omega_k) / K

#Result list
res_list_jgl[[n]] <- list(omega_k = omega_k, adj_k = adj_k, omega_0 = omega_0)

#Store jgl results
if(write){
  readr::write_rds(res_list_jgl[[n]], str_c("./sim_test/sim_res/freq/model_results/jgl_", name,".rds"))
}
#########################################
## 2.3 Frequentist RCM
#Tuning
lambda_grid <- 
  tidyr::expand_grid(
    lam1 = lam1_glasso,
    lam2 = 10^seq(-1, 2, length.out = 10),
    lam3 = 10^seq(-3, 0, length.out = 10)
  )
bic <- vector(mode = "numeric", length = nrow(lambda_grid))

#Tune lambda for independent glasso
for (i in 1:nrow(lambda_grid)) {
  res <- rcm::randCov(y, 
                 lambda1 = lambda_grid$lam1[i], 
                 lambda2 = lambda_grid$lam2[i], 
                 lambda3 = lambda_grid$lam3[i])
  omega_k <- res$Omegas
  omega_0 <- res$Omega0
  bic[i] <- bayesRCM::mbic_cal(y, omega_0, omega_k, lambda_grid$lam2[i])
}

#Grab best tune
lam_rcm <- lambda_grid[which.min(bic), ]

#Refit
rcm_out <- rcm::randCov(y, 
                   lambda1 = lam_rcm$lam1, 
                   lambda2 = lam_rcm$lam2, 
                   lambda3 = lam_rcm$lam3)

#Store and save results
omega_k <- rcm_out$Omegas
omega_k <- purrr::map(.x = seq(dim(omega_k)[3]), ~omega_k[ , , .x])
adj_k   <- purrr::map(.x = omega_k, ~abs(.x) > threshold)
omega_0 <- rcm_out$Omega0

#Result list
res_list_rcm[[n]] <- list(omega_k = omega_k, adj_k = adj_k, omega_0 = omega_0)

#If write out
if(write){
  readr::write_rds(res_list_rcm[[n]], str_c("./sim_test/sim_res/freq/model_results/rcm_", name,".rds"))
}

#Write out lambda tuning results
lambda_res <- list(lam1_glasso, lam2_glasso, lam_rcm)
readr::write_rds(lambda_res, str_c("./sim_test/sim_res/freq/model_results/lam_", name, ".rds"))
} #End loop

#Close parallel connection
parallel::stopCluster(cl)

#####################################################################################
#2. Read in results and true params
#file_path to read in results
in_path <- "./sim_test/sim_res/freq/model_results/"
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
true_param.df <- read_rds("./sim_test/sim_res/true_param_df.rds")

#Join by setting and seed
sim_res.df <-
  left_join(sim_res.df, true_param.df, by = c("setting", "seed"))

#####################################################################################
#3. Model diagnostics  
N <- nrow(sim_res.df)

#Loop through and grab model diagnostics (Diff Norms, MCC, Accuracy, etc.)
for(n in 1:N) {
  
  #Display iteration
  print(n)
  
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
  write_rds(sim_diag.df, str_c("./sim_test/sim_res/freq/sim_results/", name, ".rds"))

} #End loop
#############################################################################
#Visualize final results
#Read in frequentist
# in_path  <- "./sim_test/sim_res/freq/sim_results/"
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
# in_path  <- "./sim_test/sim_res/bayes/sim_results/"
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
