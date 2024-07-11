#Load libraries
library(tidyverse)
library(bayesRCM)
library(latex2exp)

#Fix Omega_0 for all iterations of each sim setting
n_sim <- 100

#Set alpha, lambda pairs
alpha = c(10, 30)
lambda = c(1/2, 1/4)

#Plot what tau_k looks like
set.seed(44)
tau_dist.gg <-
  expand_grid(
    alpha_tau = alpha,
    lambda_2  = lambda
  ) %>%
  mutate(
    setting      = 1:nrow(.),
    setting      = str_c("alpha=", alpha_tau, ", lambda=", lambda_2) %>%
                   as.factor() %>%
                   fct_reorder(setting),
    distribution = map2(.x = alpha_tau, .y = lambda_2,
                        ~100 - truncdist::rtrunc(spec = "gamma", a = 0, b = 100,
                                           n = 1000, shape = .x, rate = .y)),
  ) %>%
  unnest(distribution) %>%
  filter(!(alpha_tau == 10 & lambda_2 == 0.5), !(alpha_tau == 30 & lambda_2 == 0.25)) %>% 
  mutate(
    type = ifelse(setting == "alpha=10, lambda=0.25", 
                     "Fewer Outliers - Lower Sub. Variability",
                     "More Outliers - Higher Sub. Variability") %>%
           as_factor() %>%
           fct_reorder(., distribution, mean, .desc = TRUE)
  ) %>%
  ggplot(aes(x = distribution, y = after_stat(density), fill = setting, colour = setting)) +
  geom_density(alpha = 0.36, adjust = 1.6) +
  geom_histogram(binwidth = 2, alpha = 0.2) +
  facet_wrap(~type, scales = "fixed") +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 25)) +
  labs(
    y = "Density", 
    x = " ",
    title = TeX(r'(Simulated Distribution of $\tau_k \sim 100 - \Gamma_{(0, 100)}(\alpha_\tau, \lambda_2)$)'),
#    subtitle = "Subject Specific D.F. & Regularization"
  ) +
  scale_fill_viridis_d(" ") +
  scale_colour_viridis_d(" ") +
  theme_minimal() +
  theme(
    strip.placement  = "outside",
    legend.position = "bottom",
    plot.background = element_rect(fill = "white", colour = NA),
    plot.title      = element_text(hjust = 0.5, margin = margin(t = 0, r = 0, b = 5, l = 0)),
    #axis.text.x     = element_text(angle = 45, vjust = 1, hjust = 1),
    panel.spacing.x = unit(0.4, "lines"),
    panel.spacing.y = unit(1.2, "lines"),
    strip.text.x    = element_text(margin = margin(t = 0, r = 0, b = 10, l = 0)),
    plot.subtitle   = ggtext::element_markdown(hjust = 0.5)
  )

#Display
tau_dist.gg

ggsave("./sim/figures/tau_dist.png", tau_dist.gg, height = 4, width = 6, dpi = 300)

#Grid of simulation settings (length 10 for now)
sim_grid.df <-              
  expand_grid(               #(1) Design pipeline to run sim & get summary tables/figures & (2) put up on overleaf
    subjects  = c(20, 50),   #Keep at 20 for now run 50
    volumes   = c(100, 500), #Try 100, 250, 500
    rois      = 10, #Try up to 20
    alpha_tau = alpha, #Vary alpha (three different distributions of tau) <- this is key 
    lambda_2  = lambda, #Narrow + high values, moderate + most high some low, wider + more low
    prop_true_con = c(1/5),
    n_flip    = 1, #2, 4, 6, 8, --> let this vary by subject then no need to vary
    seed_0    = 4,
    seed_k    = 1:n_sim
  ) %>%
  filter(!(alpha_tau == 10 & lambda_2 == 0.5), !(alpha_tau == 30 & lambda_2 == 0.25)) %>% 
  nest(data = seed_k) %>%
  mutate(setting = 1:nrow(.)) %>%
  unnest(data) %>%
  mutate(seed_k = seed_k + (n_sim * setting)) %>% #Make sure each seed_k/subject sim is unique
  dplyr::select(setting, everything())

#For now filter to 20 subjects only (for now)
# sim_grid.df <-
#   sim_grid.df %>%
#   filter(subjects == 20)

#Dir.create
dir.create("./sim/sim_data/", showWarnings = FALSE)

#Simulate data
N          <- nrow(sim_grid.df) #number of total iterations, number of settings is 
N_settings <- N / length(unique(sim_grid.df$seed_k))

for(n in 1:N) {
  print(n)
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



