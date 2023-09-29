#Load libraries
library(tidyverse)
library(bayesRCM)

#Fix Omega_0 for all iterations of each sim setting
n_sim <- 10

#Grid of simulation settings (length 10 for now)
sim_grid.df <-
  expand_grid(
    subjects  = c(20, 50),
    volumes   = c(250, 500),
    rois      = 10,
    alpha_tau = c(20, 30), #Vary alpha (three different distributions of tau) <- this is key 
    lambda_2  = c(0.5, 0.75), #Narrow + high values, moderate + most high some low, wider + more low
    prop_true_con = c(1/5),
    n_flip    = 1, #2, 4, 6, 8, --> let this vary by subject then no need to vary
    seed_0    = 4,
    seed_k    = 1:n_sim
  ) %>%
  nest(data = seed_k) %>%
  mutate(setting = 1:nrow(.)) %>%
  unnest(data) %>%
  mutate(seed_k = seed_k + (n_sim * setting)) %>% #Make sure each seed_k/subject sim is unique
  dplyr::select(setting, everything())

#For now filter to 20 subjects only (for now)
sim_grid.df <-
  sim_grid.df %>%
  filter(subjects == 20)

#Dir.create
dir.create("./sim/sim_data/", showWarnings = FALSE)

#Simulate data
N          <- nrow(sim_grid.df) #number of total iterations, number of settings is 
N_settings <- N / length(unique(sim_grid.df$seed_k))

for(n in 1:N) {
  #print(n)
  #Sim data from bayesRCM::sim_data()
  sim_res <- 
    with(
    sim_grid.df,
    sim_data(subjects = subjects[n], volumes = volumes[n], rois = rois[n], alpha_tau = alpha_tau[n],
             lambda_2 = lambda_2[n], prop_true_con = prop_true_con[n], n_flip = n_flip[n], write = FALSE,
             seed_0   = seed_0[n], seed_k = seed_k[n])
    )
  
  #Unique name for directory to store figs
  name  <- 
    with(
      sim_grid.df,
      paste0("setting", setting[n],
            "_seed", seed_k[n])
    )
  
  #Write out
  write_rds(sim_res, str_c("./sim/sim_data/", name, ".rds"))
  
  #return
  #return(NULL)
}

#Finally, write out sim settings to be attached later
write_rds(sim_grid.df, "./sim/sim_data/sim_grid.rds")



