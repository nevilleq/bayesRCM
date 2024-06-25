#Load the good stuff
library(glasso) #individual glasso
library(JGL)    #Joint fused/group glasso
library(rcm)  
library(bayesRCM)
library(tidyverse)
library(gt)

#Set wd
setwd("/panfs/jay/groups/1/zhangl4/nevil066/bayes/bayesRCM")

#Prelims and functions
source("./sim/sim_funcs.R")

###############################################################################
#Read in results (freq & bayes)
#Read in frequentist
print("Reading freq sims....")
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
print("Reading bayes sims....")
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

#Grab sim settings
sim_settings.df <-
  sim_res.df %>%
  dplyr::select(setting, subjects:lambda_2) %>%
  mutate(setting = str_c("Setting ", setting)) %>%
  distinct() 

#Display sim settings
#sim_settings.df %>%
#  gt() %>%
#  tab_header("Simulation Settings")


#####################################################################
# OmegaO Dif Norms and Diagnostics  

print("Omega0 norm box gg")
# O Norm Viz
O_norm_box.gg <-
  sim_res.df %>%
  dplyr::select(model:seed, norm_0) %>%
  unnest(norm_0) %>%
  dplyr::select(-seed) %>%
  mutate(setting = str_c("Setting ", setting)) %>%
  group_by(setting) %>%
  #arrange(rev(desc(setting)), rev(desc(L1))) %>%
  pivot_longer(
    cols = L1:Spectral,
    names_to = "type",
    values_to = "norm"
  ) %>%
  mutate(type = as.factor(type) %>% fct_reorder(norm, mean, .desc = FALSE)) %>%
  group_by(type) %>%
  mutate(model = as.factor(model) %>% fct_reorder(norm, mean, .desc = FALSE)) %>%
  left_join(sim_settings.df, ., by = c("setting")) %>%
  mutate(volumes = str_c(volumes, " Volumes"), subjects = str_c(subjects, " Subjects")) %>%
  ggplot(aes(x = model, y = norm, colour = type, fill = type)) +
  geom_boxplot(colour = "grey80", alpha = 0.66) +
  labs(
    x = "Model",
    y = "||Truth - Estimate||",
    title = "Simulation Norm Difference Diagnostics: Omega_0"
  ) +
  scale_colour_viridis_d("Norm") +
  scale_fill_viridis_d("Norm") +
  facet_grid(subjects ~ volumes, scales = "free_y") +
  theme_minimal() +
  theme(
    strip.placement  = "outside",
    legend.position = "bottom",
    plot.background = element_rect(fill = "white", colour = NA),
    plot.title      = element_text(hjust = 0.5, margin = margin(t = 0, r = 0, b = 5, l = 0)),
    #axis.text.x     = element_text(angle = 45, vjust = 1, hjust = 1),
    panel.spacing.x = unit(0.5, "lines"),
    panel.spacing.y = unit(2, "lines"),
    strip.text.x    = element_text(margin = margin(t = 0, r = 0, b = 10, l = 0)),
    plot.subtitle   = ggtext::element_markdown(hjust = 0.5)
  )

O_norm_box.gg
ggsave("./sim/figures/O_norm_box.png", O_norm_box.gg, height = 4, width = 6)


# O norm viz ridges/density
print("Omega0 norm ridges gg")

colour <- "grey80"
O_norm.gg <-
  sim_res.df %>%
  dplyr::select(model:seed, norm_0) %>%
  unnest(norm_0) %>%
  dplyr::select(-seed) %>%
  mutate(setting = str_c("Setting ", setting)) %>%
  group_by(setting) %>%
  #arrange(rev(desc(setting)), rev(desc(L1))) %>%
  pivot_longer(
    cols = L1:Spectral,
    names_to = "type",
    values_to = "norm"
  ) %>%
  mutate(type = as.factor(type) %>% fct_reorder(norm, mean, .desc = FALSE)) %>%
  group_by(type) %>%
  mutate(model = as.factor(model) %>% fct_reorder(norm, mean, .desc = TRUE)) %>%
  left_join(sim_settings.df, ., by = c("setting")) %>%
  mutate(volumes = str_c(volumes, " Volumes"), subjects = str_c(subjects, " Subjects")) %>%
  ggplot(aes(x = norm, y = model, colour = type, fill = type)) +
  ggridges::stat_density_ridges(quantile_lines = FALSE, colour = NA,
                                alpha = 0.4, size = 0.6, jittered_points = FALSE,
                                bandwith = 0.2) +
  facet_wrap(~type, ncol = 3) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.5), minor_breaks = seq(0, 1, by = 0.5)) +
  #geom_point(position = position_jitter(width = 0.1), alpha = 0.16) +
  stat_summary(
    geom = "point",
    fun = "mean",
    size = 6,
    shape = "|",
    alpha = 1
  ) +
  labs(
    y = "Model",
    x = "||Truth - Estimate||",
    title = "Simulation Norm Difference Diagnostics: Omega_0"
  ) +
  scale_colour_viridis_d("Norm") +
  scale_fill_viridis_d("Norm") +
  facet_grid(subjects ~ volumes, scales = "free_y") +
  theme_minimal() +
  theme(
    strip.placement  = "outside",
    legend.position = "bottom",
    plot.background = element_rect(fill = "white", colour = NA),
    plot.title      = element_text(hjust = 0.5, margin = margin(t = 0, r = 0, b = 5, l = 0)),
    #axis.text.x     = element_text(angle = 45, vjust = 1, hjust = 1),
    panel.spacing.x = unit(0.5, "lines"),
    panel.spacing.y = unit(2, "lines"),
    strip.text.x    = element_text(margin = margin(t = 0, r = 0, b = 10, l = 0)),
    plot.subtitle   = ggtext::element_markdown(hjust = 0.5)
  )

O_norm.gg
ggsave("./sim/figures/O_norm_ridges.png", O_norm.gg, height = 6, width = 8)

# O-Diag Viz
print("Omega0 diag box gg")

O_diag_box.gg <-
  sim_res.df %>%
  dplyr::select(model:seed, O_diag) %>%
  unnest(O_diag) %>%
  dplyr::select(-seed) %>%
  mutate(
    setting = as.factor(str_c("Setting ", setting)) %>%
      fct_reorder(setting, .desc = FALSE)
  ) %>%
  dplyr::select(model, setting, mcc, accuracy, kappa) %>%
  pivot_longer(
    cols = mcc:kappa,
    names_to = "type",
    values_to = "stat"
  ) %>%
  mutate(type = as.factor(type) %>% fct_reorder(stat, mean, .desc = TRUE)) %>%
  group_by(type) %>%
  mutate(model = as.factor(model) %>% fct_reorder(stat, mean, .desc = TRUE)) %>%
  left_join(sim_settings.df, ., by = c("setting")) %>%
  mutate(volumes = str_c(volumes, " Volumes"), subjects = str_c(subjects, " Subjects")) %>%
  ggplot(aes(x = model, y = stat, colour = type, fill = type)) +
  geom_boxplot(colour = "grey80", alpha = 0.66) +
  labs(
    x = "Model",
    y = "Statistic",
    title = "Simulation Diagnostics: Omega_0"
  ) +
  scale_colour_viridis_d("Diagnostic") +
  scale_fill_viridis_d("Diagnostic") +
  scale_y_continuous(breaks = seq(0, 1, by = 0.25), minor_breaks = seq(0, 1, by = 0.25), limits = c(0, 1)) +
  facet_grid(subjects ~ volumes, scales = "free_y") +
  theme_minimal() +
  theme(
    strip.placement  = "outside",
    legend.position = "bottom",
    plot.background = element_rect(fill = "white", colour = NA),
    plot.title      = element_text(hjust = 0.5, margin = margin(t = 0, r = 0, b = 5, l = 0)),
    #axis.text.x     = element_text(angle = 45, vjust = 1, hjust = 1),
    panel.spacing.x = unit(0.5, "lines"),
    panel.spacing.y = unit(2, "lines"),
    strip.text.x    = element_text(margin = margin(t = 0, r = 0, b = 10, l = 0)),
    plot.subtitle   = ggtext::element_markdown(hjust = 0.5)
  )

O_diag_box.gg
ggsave("./sim/figures/O_diag_box.png", O_diag_box.gg, height = 6, width = 8)

# O-diag ridges
print("Omega0 diag ridges gg")

O_diag_ridges.gg <-
  sim_res.df %>%
  dplyr::select(model:seed, O_diag) %>%
  unnest(O_diag) %>%
  dplyr::select(-seed) %>%
  mutate(
    setting = as.factor(str_c("Setting ", setting)) %>%
      fct_reorder(setting, .desc = FALSE)
  ) %>%
  dplyr::select(model, setting, mcc, accuracy, kappa) %>%
  pivot_longer(
    cols = mcc:kappa,
    names_to = "type",
    values_to = "stat"
  ) %>%
  mutate(type = as.factor(str_to_title(type)) %>% fct_reorder(stat, mean, .desc = TRUE)) %>%
  group_by(type) %>%
  mutate(model = as.factor(model) %>% fct_reorder(stat, mean, .desc = FALSE)) %>%
  left_join(sim_settings.df, ., by = c("setting")) %>%
  mutate(volumes = str_c(volumes, " Volumes"), subjects = str_c(subjects, " Subjects")) %>%
  ggplot(aes(x = stat, y = model, colour = type, fill = type)) +
  ggridges::stat_density_ridges(quantile_lines = FALSE, colour = NA,
                                alpha = 0.4, size = 0.6, jittered_points = TRUE,
                                bandwith = 0.2) +
  facet_wrap(~type, ncol = 3) +
  scale_x_continuous(breaks = seq(0, 1, by = 1), minor_breaks = seq(0, 1, by = 1)) +
  #geom_point(position = position_jitter(width = 0.1), alpha = 0.16) +
  stat_summary(
    geom = "point",
    fun = "mean",
 #   col = colour,
    size = 6,
    shape = "|",
    alpha = 1
  ) +
  labs(
    y = "Model",
    x = "Statistic",
    title = "Simulation Diagnostics: Omega_0"
  ) +
  scale_colour_viridis_d("Diagnostic") +
  scale_fill_viridis_d("Diagnostic") +
  facet_grid(subjects ~ volumes, scales = "free_y") +
  theme_minimal() +
  theme(
    strip.placement  = "outside",
    legend.position = "bottom",
    plot.background = element_rect(fill = "white", colour = NA),
    plot.title      = element_text(hjust = 0.5, margin = margin(t = 0, r = 0, b = 5, l = 0)),
    #axis.text.x     = element_text(angle = 45, vjust = 1, hjust = 1),
    panel.spacing.x = unit(0.5, "lines"),
    panel.spacing.y = unit(2, "lines"),
    strip.text.x    = element_text(margin = margin(t = 0, r = 0, b = 10, l = 0)),
    plot.subtitle   = ggtext::element_markdown(hjust = 0.5)
  )

O_diag_ridges.gg
ggsave("./sim/figures/O_diag_ridges.png", O_diag_ridges.gg, height = 6, width = 8)

stop("Done for now")
#################################################################################
#OmegaK Norm
# K-norm diagnostics
print("OmegaK Norms")

k_norm_box <- function(alpha, lambda) {
k_norm_box.gg <-
  sim_res.df %>%
  dplyr::select(model:seed, norm_k) %>%
  unnest(norm_k) %>%
  mutate(setting = str_c("Setting ", setting)) %>%
  dplyr::select(-seed) %>%
  pivot_longer(
    cols = L1:Spectral,
    names_to = "type",
    values_to = "norm"
  ) %>%
  mutate(type = as.factor(str_remove(type, str_c("_", stat))) %>% fct_reorder(norm, mean, .desc = FALSE)) %>%
  group_by(type) %>%
  mutate(model = as.factor(model) %>% fct_reorder(norm, mean, .desc = FALSE)) %>%
  left_join(sim_settings.df, ., by = c("setting")) %>%
  filter(alpha == alpha, lambda_2 == lambda) %>%
  mutate(volumes = str_c(volumes, " Volumes"), subjects = str_c(subjects, " Subjects")) %>%
  ggplot(aes(x = model, y = norm, colour = type, fill = type)) +
  geom_boxplot(colour = "grey80", alpha = 0.66) +
  labs(
    x = "Model",
    y = "Avg. ||Truth - Estimate||",
    title = "Simulation Norm Difference Diagnostics: Omega_k",
    subtitle = sprintf("Alpha = %g, Lambda2 = %g", alpha, lambda),
  #  caption = "| denotes the average norm."
  ) +
  scale_colour_viridis_d("Norm") +
  scale_fill_viridis_d("Norm") +
  facet_grid(subjects ~ volumes, scales = "free_y") +
  theme_minimal() +
  theme(
    strip.placement  = "outside",
    legend.position = "bottom",
    plot.background = element_rect(fill = "white", colour = NA),
    plot.title      = element_text(hjust = 0.5, margin = margin(t = 0, r = 0, b = 5, l = 0)),
    #axis.text.x     = element_text(angle = 45, vjust = 1, hjust = 1),
    panel.spacing.x = unit(0.5, "lines"),
    panel.spacing.y = unit(2, "lines"),
    strip.text.x    = element_text(margin = margin(t = 0, r = 0, b = 10, l = 0)),
    plot.subtitle   = ggtext::element_markdown(hjust = 0.5)
  )
  return(k_norm_box.gg)
}

#Alpha = 20, 40; lambda_2 = 0.5, 0.75
alpha = c(20, 40, 20, 40)
lambda = c(0.5, 0.75, 0.75, 0.5)

#Test
k_norm_box.gg <- k_norm_box(alpha = 20, lambda = 0.5)
k_norm_box.gg 

#Iterate over all alpha, lambda & save
k_box_list <- map2(.x = alpha, .y = lambda,
                   ~k_norm_box(.x, .y) %>%
                     ggsave(sprintf("./sim/figures/k_norm_box_alpha%g_lambda%g.png", .x, .y), ., height= 6, width = 8))


# K-norm ridges
print("OmegaK Norms ridges")

k_norm_ridge <- function(alpha, lambda) {

k_norm_ridges.gg <-
  sim_res.df %>%
  dplyr::select(model:seed, norm_k) %>%
  unnest(norm_k) %>%
  mutate(setting = str_c("Setting ", setting)) %>%
  dplyr::select(-seed) %>%
  group_by(model, setting) %>%
  # summarise(
  #   across(
  #     .cols = where(is.numeric),
  #     .fns  = list(mean = mean, median = median, sd = sd),
  #     .names = "{.col}_{.fn}"
  #   ),
  #   .groups = "drop"
  # ) %>%
  pivot_longer(
    cols = L1:Spectral,
    names_to = "type",
    values_to = "norm"
  ) %>%
  mutate(type = as.factor(type) %>% fct_reorder(norm, mean, .desc = FALSE)) %>%
  group_by(type) %>%
  mutate(model = as.factor(model) %>% fct_reorder(norm, mean, .desc = FALSE)) %>%
  left_join(sim_settings.df, ., by = c("setting")) %>%
  filter(alpha == alpha, lambda_2 == lambda) %>% #Filter statement for alpha lambda
  mutate(volumes = str_c(volumes, " Volumes"), subjects = str_c(subjects, " Subjects")) %>%
  ggplot(aes(x = norm, y = model, colour = type, fill = type)) +
  ggridges::stat_density_ridges(quantile_lines = FALSE, colour = NA,
                                alpha = 0.6, size = 0.6, jittered_points = FALSE,
                                bandwith = 0.2) +
  #facet_wrap(~type, ncol = 3) +
  scale_x_continuous(limits = c(0, 5)) +
  #geom_point(position = position_jitter(width = 0.1), alpha = 0.16) +
  stat_summary(
    geom = "point",
    fun = "mean",
   # col = colour,
    size = 6,
    shape = "|",
    alpha = 1
  ) +
  labs(
    y = "Model",
    x = "Statistic",
    title = "Simulation Diagnostics: Omega_k",
    subtitle = sprintf("Alpha = %g, Lambda2 = %g", alpha, lambda),
    caption = "| denotes the average norm."
  ) +
  scale_colour_viridis_d("Norm") +
  scale_fill_viridis_d("Norm") +
  facet_grid(subjects ~ volumes, scales = "free_y") +
  theme_minimal() +
  theme(
    strip.placement  = "outside",
    legend.position = "bottom",
    plot.background = element_rect(fill = "white", colour = NA),
    plot.title      = element_text(hjust = 0.5, margin = margin(t = 0, r = 0, b = 5, l = 0)),
    #axis.text.x     = element_text(angle = 45, vjust = 1, hjust = 1),
    panel.spacing.x = unit(0.5, "lines"),
    panel.spacing.y = unit(2, "lines"),
    strip.text.x    = element_text(margin = margin(t = 0, r = 0, b = 10, l = 0)),
    plot.subtitle   = ggtext::element_markdown(hjust = 0.5)
  )
  
  return(k_norm_ridges.gg)
}

#Test
k_norm_ridges.gg <- k_norm_ridge(alpha = 20, lambda = 0.5)
k_norm_ridges.gg

#Iterate over all pairs and save
#Iterate over all alpha, lambda & save
k_box_list <- map2(.x = alpha, .y = lambda,
                   ~k_norm_ridge(.x, .y) %>%
                     ggsave(sprintf("./sim/figures/k_norm_ridges_alpha%g_lambda%g.png", .x, .y), ., height= 6, width = 8))

#ggsave("./sim/figures/k_norm_ridges.png", k_norm_ridges.gg, height = 4, width = 6)

#################################################################################
#OmegaK Diag
print("OmegaK diag")

# K-norm diagnostics
k_diag_box <- function(alpha, lambda) {
k_diag_box.gg <-
  sim_res.df %>%
  dplyr::select(model:seed, k_diag) %>%
  unnest(k_diag) %>%
  mutate(setting = str_c("Setting ", setting)) %>%
  dplyr::select(model:kappa, -c(seed)) %>%
  group_by(model, setting) %>%
  summarise(
    across(
      .cols = where(is.numeric),
      .fns  = list(mean = mean, median = median, sd = sd),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = contains(stat),
    names_to = "type",
    values_to = "norm"
  ) %>%
  mutate(type = as.factor(str_remove(type, str_c("_", stat))) %>% fct_reorder(norm, mean, .desc = FALSE)) %>%
  group_by(type) %>%
  mutate(model = as.factor(model) %>% fct_reorder(norm, mean, .desc = FALSE)) %>%
  left_join(sim_settings.df, ., by = c("setting")) %>%
  filter(alpha == alpha, lambda_2 == lambda) %>%
  mutate(volumes = str_c(volumes, " Volumes"), subjects = str_c(subjects, " Subjects")) %>%
  ggplot(aes(x = model, y = norm, colour = type, fill = type)) +
  geom_boxplot(colour = "grey80", alpha = 0.66) +
  labs(
    x = "Model",
    y = "Avgerage Across Subjects",
    title = "Simulation Diagnostics: Omega_k",
    subtitle = sprintf("Alpha = %g, Lambda2 = %g", alpha, lambda),
   # caption = "| denotes the average."
  ) +
  scale_colour_viridis_d("Norm") +
  scale_fill_viridis_d("Norm") +
  facet_grid(subjects ~ volumes, scales = "free_y") +
  theme_minimal() +
  theme(
    strip.placement  = "outside",
    #legend.position = "none",
    plot.background = element_rect(fill = "white", colour = NA),
    plot.title      = element_text(hjust = 0.5, margin = margin(t = 0, r = 0, b = 5, l = 0)),
    #axis.text.x     = element_text(angle = 45, vjust = 1, hjust = 1),
    panel.spacing.x = unit(0.5, "lines"),
    panel.spacing.y = unit(2, "lines"),
    strip.text.x    = element_text(margin = margin(t = 0, r = 0, b = 10, l = 0)),
    plot.subtitle   = ggtext::element_markdown(hjust = 0.5)
  )

  return(k_diag_box.gg)

}
k_box_list <- map2(.x = alpha, .y = lambda,
                   ~k_diag_box(.x, .y) %>%
                     ggsave(sprintf("./sim/figures/k_diag_box_alpha%g_lambda%g.png", .x, .y), ., height= 6, width = 8))
# K-diag ridges
k_diag_ridge <- function(alpha, lambda) {
k_diag_ridges.gg <-
  sim_res.df %>%
  dplyr::select(model:seed, k_diag) %>%
  unnest(k_diag) %>%
  mutate(setting = str_c("Setting ", setting)) %>%
  dplyr::select(model:kappa, -c(seed)) %>%
  group_by(model, setting) %>%
  summarise(
    across(
      .cols = where(is.numeric),
      .fns  = list(mean = mean, median = median, sd = sd),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = contains(stat),
    names_to = "type",
    values_to = "value"
  ) %>%
  mutate(type = as.factor(str_to_title(str_remove(type, str_c("_", stat)))) %>% fct_reorder(value, mean, .desc = TRUE)) %>%
  group_by(type) %>%
  mutate(model = as.factor(model) %>% fct_reorder(value, mean, .desc = FALSE)) %>%
  left_join(sim_settings.df, ., by = c("setting")) %>%
  filter(alpha == alpha, lambda_2 == lambda) %>%
  mutate(volumes = str_c(volumes, " Volumes"), subjects = str_c(subjects, " Subjects")) %>%
  ggplot(aes(x = value, y = model, colour = type, fill = type)) +
  ggridges::stat_density_ridges(quantile_lines = FALSE, colour = colour,
                                alpha = 0.6, size = 0.6, jittered_points = TRUE,
                                bandwith = 0.2) +
  facet_wrap(~type, ncol = 3) +
  scale_x_continuous(breaks = seq(0, 1, by = 1), minor_breaks = seq(0, 1, by = 1), limits = c(0, 1)) +
  #geom_point(position = position_jitter(width = 0.1), alpha = 0.16) +
  stat_summary(
    geom = "point",
    fun = "mean",
    col = colour,
    size = 6,
    shape = "|",
    alpha = 1
  ) +
  labs(
    y = "Model",
    x = "Statistic",
    title = "Simulation Diagnostics: Omega_k",
    subtitle = sprintf("Alpha = %g, Lambda2 = %g", alpha, lambda),
    caption = "| denotes the average."
  ) +
  scale_colour_viridis_d("Diagnostic") +
  scale_fill_viridis_d("Diagnostic") +
  facet_grid(subjects ~ volumes, scales = "free_y") +
  theme_minimal() +
  theme(
    strip.placement  = "outside",
    legend.position = "none",
    plot.background = element_rect(fill = "white", colour = NA),
    plot.title      = element_text(hjust = 0.5, margin = margin(t = 0, r = 0, b = 5, l = 0)),
    #axis.text.x     = element_text(angle = 45, vjust = 1, hjust = 1),
    panel.spacing.x = unit(0.5, "lines"),
    panel.spacing.y = unit(2, "lines"),
    strip.text.x    = element_text(margin = margin(t = 0, r = 0, b = 10, l = 0)),
    plot.subtitle   = ggtext::element_markdown(hjust = 0.5)
  )

  return(k_diag_ridges.gg)
}

k_diag_list <- map2(.x = alpha, .y = lambda,
                   ~k_diag_ridge(.x, .y) %>%
                     ggsave(sprintf("./sim/figures/k_diag_ridges_alpha%g_lambda%g.png", .x, .y), ., height= 6, width = 8))


#############################################################################################
#K-tables
print("OmegaK Norms gt")

omega_k_norm_sum.gt <-
  sim_res.df %>%
  dplyr::select(model:seed, k_diag) %>%
  unnest(k_diag) %>%
  mutate(setting = str_c("Setting ", setting)) %>%
  dplyr::select(model:kappa, -c(seed)) %>%
  group_by(model, setting) %>%
  summarise(
    across(
      .cols = where(is.numeric),
      .fns  = list(mean = mean, median = median, sd = sd),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  ) %>%
  group_by(setting) %>%
  gt() %>%
  tab_header("Omega_K Diff-Norm Summarised over Subjects by Setting")

print("OmegaK Norms")
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

#omega_k_diag_sum.gt
