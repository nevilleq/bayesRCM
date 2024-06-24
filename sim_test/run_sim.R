#Load the good stuff!
#For rcm(s)
library(Rcpp)
library(RcppArmadillo)
library(GIGrvg) #Generalized inverse gaussian
library(glasso) #individual glasso
library(JGL)    #Joint fused/group glasso
library(rcm)    #previous freq. bi-level graphical model
library(bayesRCM) #our method
library(doParallel)

#Display
library(tidyverse)
library(gt)
library(igraph)
library(ggraph)
#library(brainconn) figure this out
library(tidygraph)
#library(Matrix)

#My Colours (from viridis)
my_purple <- "#440154FF"
my_yellow <- "#FDE725FF"

#Set Theme for ggplot2
theme_set(theme_minimal() + 
            theme(plot.title = element_text(hjust = 0.5),
                  plot.subtitle = element_text(hjust = 0.5),
                  legend.position = "bottom",
                  plot.background  = element_rect(fill = "white", colour = "white")
            )
)

#function for tidy date
tidy_date <- function() {
  Sys.Date() %>%
    str_replace_all("-", "_")
}

#Function to fill matrix w/upper triangle and make symettric (tidy results)
fill_mat <- function(res_vec, p) {
  null.mat <- matrix(data = 0, nrow = p, ncol = p)
  null.mat[upper.tri(null.mat, diag = TRUE)] <- res_vec
  null_temp.mat <- null.mat
  diag(null_temp.mat) <- 0
  full_precision <- null.mat + t(null_temp.mat)
  return(full_precision)
}

#Get Binary Class Diagnostics for Adj_k
get_bin_diag_k <- function(pred_list, truth_list) {
  # 
  # pred_list  <- a$adj_k[[1]]
  # truth_list <- Adj_k
  #Grab list of upper triangle vectors
  pred_up_tri   <- map(.x = pred_list, ~.x[upper.tri(.x, diag = FALSE)])
  truth_up_tri  <- map(.x = truth_list, ~.x[upper.tri(.x, diag = FALSE)])
  
  #Make table and get confusion matrix results & diagnostics
  table         <- map2(.x = pred_up_tri, .y = truth_up_tri, 
                        ~table(.x, .y))
  conf_list     <- map(.x = table, ~caret::confusionMatrix(.x, positive = "TRUE"))
  
  #Put together final diagnostic tibble rows = subjects columns = which sub & diagnostic stats
  diagnostic.df <- 
    map_df(.x = conf_list, ~as_tibble(as.data.frame(t(.x$overall)))) %>%
    janitor::clean_names() %>%
    mutate(
      subject = str_c("Sub. ", 1:length(pred_list)), 
      mcc     = map2_dbl(.x = pred_up_tri, .y = truth_up_tri,
                         ~mltools::mcc(preds = .x, actuals = .y))
    ) %>%
    dplyr::select(subject, mcc, everything())
  
  #Return diag .df to unnest
  return(diagnostic.df)
}

get_bin_diag_0 <- function(omega, truth) {
  
  #Upper triangle vec
  pred_up_tri   <- omega[upper.tri(omega, diag = FALSE)]
  truth_up_tri  <- truth[upper.tri(truth, diag = FALSE)]
  
  #Make table and get confusion matrix results & diagnostics
  table <- table(pred_up_tri, truth_up_tri)
  conf  <- caret::confusionMatrix(table, positive = "TRUE")
  
  #Put together final diagnostic tibble rows = subjects columns = which sub & diagnostic stats
  diagnostic.df <- 
    as_tibble(as.data.frame(t(conf$overall))) %>%
    janitor::clean_names() %>%
    mutate(
      mcc = mltools::mcc(preds = pred_up_tri, actuals = truth_up_tri)
    ) %>%
    dplyr::select(mcc, everything())
  
  return(diagnostic.df)
}

#create storage directories
dir.create("./sim/sim_res", showWarnings = FALSE)
dir.create("./sim/figures", showWarnings = FALSE)
#dir.create("./sim/sim_data")
#dir.create("./sim/sim_truth")

#######################################################################################

# 1. Generate some data  
#Set parameters
subjects  <- 20
volumes   <- 500
rois      <- 10
alpha_tau <- 25
lambda_2  <- alpha_tau / 50  
prop_true_con <- 1/5
n_flip    <- 1
seed      <- 4
write     <- TRUE
threshold <- 0.001

run_sim <- function(subjects = 20, volumes = 500, rois = 10, alpha_tau = 25,
                    lambda_2 = 1/2, prop_true_con = 1/5, n_flip = 1, write = FALSE,
                    seed = 4, threshold = 0.001){

#Simulate data
sim_res <- sim_data(subjects = subjects, volumes = volumes, rois = rois, alpha_tau = alpha_tau,
                    lambda_2 = lambda_2, prop_true_con = prop_true_con, n_flip = n_flip, write = FALSE, seed = seed)

#Unique name for directory to store figs
name  <- str_c("sub", subjects,
               "_vol", volumes,
               "_atau", alpha_tau,
               "_lam", lambda_2,
               "_prop", prop_true_con,
               "_nf", n_flip,
               "_seed", seed)

########################################################################################

# 2. Fit the models  
## 2.1 bayesRCM
#result   <- rcm(sim_res$data_list, n_samples = 2500, n_burn = 500, n_cores = 4, n_updates = 5)
result      <- read_rds("./sim/sim_data/result_04_12_23.rds") #Temporary results for testing
#res_list <- list(sim_settings = sim_res$sim_settings, result = result, true_params = sim_res$true_params, data = sim_res$data_list) 

#Write out true params and data
if(write) {
write_rds(sim_res$true_params, str_c("./sim/sim_truth/", name,".rds"))
write_rds(sim_res$data_list, str_c("./sim/sim_data/", name,".rds"))
}

#Read in bayesRCM result
#result         <- read_rds(str_c("./sim/sim_res/bayesrcm_", name,".rds")) #Uncomment for full run
sim_param.df   <- sim_res$sim_settings
true_params    <- sim_res$true_params
#result         <- sim_res$result

#Data and dimensions
y <- data_list <- sim_res$data_list #Subject data list
p <- ncol(y[[1]]) #Dimension of data rois
n_iter <- ncol(result$omega_0)
K <- length(y)

#Posterior mean and median omega_k matrices
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

omega_k <- omegak.df$mean_mat
adj_k   <- map(.x = omega_k, ~abs(.x) > threshold)
#omega_k <- omegak.df$med_mat
#adj_k   <- map(.x = omega_k, ~abs(.x) > threshold)

#Posterior mean and median omega_0 matrices
omega_0 <- apply(result$omega_0, 1, mean) %>% fill_mat(., p)
#omega_0 <- apply(result$omega_0, 1, median) %>% fill_mat(., p)

res_list_bayesrcm <- list(omega_k = omega_k, adj_k = adj_k, omega_0 = omega_0)

if(write) {
write_rds(res_list_bayesrcm, str_c("./sim/sim_res/bayesrcm_", name,".rds"))
}
#################################
## 2.2 Independent glasso  
#Omega_k ind glasso est - tune lambda by bic and then grab G and Omega_0
lambda_grid <- 10^seq(-3, 1, length.out = 40)
bic         <- vector(mode = "numeric", length = length(lambda_grid))

#Tune lambda for independent glasso
for (i in 1:length(lambda_grid)) {
  omega_k <- ind_graphs(y, lambda_grid[i])
  bic[i]     <- bic_cal(y, omega_k)
}

#Compute initial est. for omega_k, G_k/adj_k via best mBIC, and Omega_0
lam1_glasso <- lambda_grid[which.min(bic)]
omega_k <- ind_graphs(y, lam1_glasso)
adj_k   <- map(.x = omega_k, ~abs(.x) > threshold)
omega_0 <- Reduce("+", omega_k) / K

#Result list
res_list_ig <- list(omega_k = omega_k, adj_k = adj_k, omega_0 = omega_0)

#Store ind_glasso results
if(write){
write_rds(res_list_ig, str_c("./sim/sim_res/iglasso_", name,".rds"))
}
#######################################
## 2.3 Fused group glasso 
#Tuning
lambda_grid <- 
  expand_grid(
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
  bic[i] <- bic_cal(y, omega_k)
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
adj_k   <- map(.x = omega_k, ~abs(.x) > threshold)
omega_0 <- Reduce("+", omega_k) / K

#Result list
res_list_jgl <- list(omega_k = omega_k, adj_k = adj_k, omega_0 = omega_0)

#Store jgl results
if(write){
write_rds(res_list_jgl, str_c("./sim/sim_res/jgl_", name,".rds"))
}
#########################################
## 2.4 Frequentist RCM
#Tuning
lambda_grid <- 
  expand_grid(
    lam1 = lam1_glasso,
    lam2 =  10^seq(-3, 1, length.out = 10),
    lam3 = 10^seq(-3, 1, length.out = 10)
  )
bic <- vector(mode = "numeric", length = nrow(lambda_grid))

#Tune lambda for independent glasso
for (i in 1:nrow(lambda_grid)) {
  res <- randCov(y, 
                 lambda1 = lambda_grid$lam1[i], 
                 lambda2 = lambda_grid$lam2[i], 
                 lambda3 = lambda_grid$lam3[i])
  omega_k <- res$Omegas
  omega_0 <- res$Omega0
  bic[i] <- mbic_cal(y, omega_0, omega_k, lambda_grid$lam2[i])
}

#Grab best tune
lam_rcm <- lambda_grid[which.min(bic), ]

#Refit
rcm_out <- randCov(y, 
                   lambda1 = lam_rcm$lam1, 
                   lambda2 = lam_rcm$lam2, 
                   lambda3 = lam_rcm$lam3)

#Store and save results
omega_k <- rcm_out$Omegas
omega_k <- map(.x = seq(dim(omega_k)[3]), ~omega_k[ , , .x])
adj_k   <- map(.x = omega_k, ~abs(.x) > threshold)
omega_0 <- rcm_out$Omega0

#Result list
res_list_rcm <- list(omega_k = omega_k, adj_k = adj_k, omega_0 = omega_0)

#Store jgl results
  if(write){
  write_rds(res_list_rcm, str_c("./sim/sim_res/rcm_", name,".rds"))
  }

#########################################################################################

# 3. Diagnostics  
##Read in saved results
root_dir  <- "./sim/sim_res/"
in_files  <- list.files(root_dir)
meta_data <- str_remove(in_files, ".rds") 
var_names <- c("model", "subjects", "volumes", "alpha_tau", 
               "lambda_2", "prop_true_conn", "n_flip", "seed")

#Read in sim results and tidy
sim_result.df <-
  tibble(
    meta_data = str_remove(in_files, ".rds"),
    in_path   = str_c(root_dir, in_files)
  ) %>% 
  tidyr::separate(col = meta_data, sep = "_", into = var_names) %>% #Separate params into columns
  mutate(
    across(
      .cols = subjects:seed,
      ~parse_number(.x) #Go across params and grab numeric
    )
  ) %>%
  mutate(
    result = map(.x = in_path, ~read_rds(.x))
  ) %>%
  unnest(result) %>%
  mutate(
    name = names(result)
  ) %>%
  pivot_wider(
    names_from = "name",
    values_from = "result"
  )

#True parameters
Omega_0 <- true_params$omega_0
Adj_0   <- abs(Omega_0) > threshold
Omega_k <- true_params$omega_k
Adj_k   <- map(.x = true_params$omega_k, ~abs(.x) > threshold)

#Diff func
get_diff <- function(x_list, y_list) {
  map2(.x = x_list, .y = y_list, ~.y - .x)
}

#Get norms
get_norms <- function(x_list) {
  tibble(
    subject   = str_c("Sub. ", 1:length(x_list)),
    L1        = map_dbl(.x = x_list, ~norm(.x, "1")),
    Frobenius = map_dbl(.x = x_list, ~norm(.x, "F")),
    Spectral  = map_dbl(.x = x_list, ~norm(.x, "2"))
  )
}

#Diagnostic .df
sim_diag.df <-
  sim_result.df %>%
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
  dplyr::select(-c(in_path:omega_0, contains(c("diff", "adj"))))

  #Return final diag data frame for plotting/summary
  return(sim_diag.df)  

} #End function

#Run sim (test)
sim_diag.df <- map_df(.x = 4:5, ~run_sim(seed = .x))

#Omega_0 norm results
sim_diag.df %>%
  dplyr::select(model:seed, norm_0) %>%
  unnest(norm_0)

#Omega_0 diagnostics  
sim_diag.df %>%
  dplyr::select(model:seed, O_diag) %>%
  unnest(O_diag)

#Omega_k norm results  
sim_diag.df %>%
  dplyr::select(model:seed, norm_k) %>%
  unnest(norm_k)

#Omega_k diagnostics
sim_diag.df %>%
  dplyr::select(model:seed, k_diag) %>%
  unnest(k_diag)
