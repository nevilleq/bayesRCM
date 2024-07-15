#Load the good stuff
library(glasso) #individual glasso
library(JGL)    #Joint fused/group glasso
library(rcm)  
library(bayesRCM)
library(tidyverse)
library(gt)

#Set wd
#setwd("/panfs/jay/groups/1/zhangl4/nevil066/bayes/bayesRCM")

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

mod_order <- c("bayesRCM", "rcm", "jgl", "iglasso")
diag_order <- c("MCC", "TPR", "TNR", "FDR")
norm_order <- c("Spectral", "L1", "Frobenius")
height <- 6
width <- 8


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
  mutate(type = as.factor(type) %>% fct_relevel(norm_order)) %>%
  left_join(sim_settings.df, ., by = c("setting")) %>%
  mutate(
    volumes = str_c(volumes, " Volumes"), 
    subjects = str_c(subjects, " Subjects"),
    outliers = ifelse(alpha_tau == 10 & lambda_2 == 0.25, 
                      "Fewer Outliers",
                      "More Outliers") %>%
               as_factor() %>%
               fct_relevel("Fewer Outliers")
  ) %>%
  group_by(type) %>%
  mutate(model = as.factor(model) %>% fct_relevel(mod_order)) %>%
  ggplot(aes(x = type, y = norm, colour = model, fill = model)) +
  geom_boxplot(colour = "grey80", alpha = 0.72, size = 0.4) +
  labs(
    x = "Precision Norm Difference",
    y = "Group-level ||Truth - Estimate||",
    #title = "Simulation Norm Difference Diagnostics: Omega_0"
  ) +
  scale_colour_viridis_d("Model") +
  scale_fill_viridis_d("Model", option = "magma") +
  facet_grid(subjects ~ outliers + volumes, scales = "fixed") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25)) +
  theme_minimal() +
  theme(
    strip.placement  = "outside",
    legend.position = "top",
    plot.background = element_rect(fill = "white", colour = NA),
    plot.title      = element_text(hjust = 0.5, margin = margin(t = 0, r = 0, b = 5, l = 0)),
    #axis.text.x     = element_text(angle = 45, vjust = 1, hjust = 1),
    panel.spacing.x = unit(0.4, "lines"),
    panel.spacing.y = unit(1.2, "lines"),
    strip.text.x    = element_text(margin = margin(t = 0, r = 0, b = 10, l = 0)),
    plot.subtitle   = ggtext::element_markdown(hjust = 0.5)
  )

O_norm_box.gg
ggsave("./sim/figures/O_norm_box.png", O_norm_box.gg, height = height, width = width, dpi = 300)


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
  mutate(type = as.factor(type) %>% fct_relevel(norm_order)) %>%
  left_join(sim_settings.df, ., by = c("setting")) %>%
  mutate(
    volumes = str_c(volumes, " Volumes"), 
    subjects = str_c(subjects, " Subjects"),
    outliers = ifelse(alpha_tau == 10 & lambda_2 == 0.25, 
                      "Fewer Outliers",
                      "More Outliers") %>%
      as_factor() %>%
      fct_relevel("Fewer Outliers")
  ) %>%
  group_by(type) %>%
  mutate(model = as.factor(model) %>% fct_relevel(mod_order)) %>%
  ggplot(aes(x = norm, y = type, colour = model, fill = model)) +
  ggridges::stat_density_ridges(quantile_lines = FALSE,
                                alpha = 0.32, jittered_points = FALSE,
                                bandwidth = 0.08) +
 # facet_wrap(~type, ncol = 3) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.5), minor_breaks = seq(0, 1, by = 0.5), limits = c(0, 1)) +
  #geom_point(position = position_jitter(width = 0.1), alpha = 0.16) +
  stat_summary(
    geom = "point",
    fun = "median",
    size = 6,
    shape = "|",
    alpha = 1
  ) +
  labs(
    x = "Norm",
    y = "Group-level Precision ||Truth - Estimate||",
    #title = "Simulation Norm Difference Diagnostics: Omega_0"
  ) +
  scale_colour_viridis_d("Model", option = "magma") +
  scale_fill_viridis_d("Model", option = "magma") +
  facet_grid(subjects ~ outliers + volumes, scales = "free_y") +
  theme_minimal() +
  theme(
    strip.placement  = "outside",
    legend.position = "top",
    plot.background = element_rect(fill = "white", colour = NA),
    plot.title      = element_text(hjust = 0.5, margin = margin(t = 0, r = 0, b = 5, l = 0)),
    #axis.text.x     = element_text(angle = 45, vjust = 1, hjust = 1),
    panel.spacing.x = unit(1.2, "lines"),
    panel.spacing.y = unit(0.4, "lines"),
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
  dplyr::select(model, setting, FDR, TPR, TNR, MCC) %>%
  pivot_longer(
    cols = FDR:MCC,
    names_to = "type",
    values_to = "stat"
  ) %>%
  mutate(type = as.factor(type) %>% fct_relevel(diag_order)) %>%
  left_join(sim_settings.df, ., by = c("setting")) %>%
  mutate(
    volumes = str_c(volumes, " Volumes"), 
    subjects = str_c(subjects, " Subjects"),
    outliers = ifelse(alpha_tau == 10 & lambda_2 == 0.25, 
                      "Fewer Outliers",
                      "More Outliers") %>%
      as_factor() %>%
      fct_relevel("Fewer Outliers")
  ) %>%
  mutate(model = as.factor(model) %>% fct_relevel(mod_order)) %>%
  ggplot(aes(x = type, y = stat, colour = model, fill = model)) +
  geom_boxplot(colour = "grey80", alpha = 0.8, size = 0.4, width = 0.76) +
  labs(
    x = "Binary Diagnostics of Group Graph/Network",
    y = "Score",
 #   title = "Simulation Diagnostics: Omega_0",
  ) +
  scale_colour_viridis_d("Model") +
  scale_fill_viridis_d("Model", option = "magma") +
  scale_y_continuous(breaks = seq(0, 1, by = 0.25), minor_breaks = seq(0, 1, by = 0.25), limits = c(0, 1)) +
  facet_grid(subjects ~ outliers + volumes, scales = "free_y") +
  theme_minimal() +
  theme(
    strip.placement  = "outside",
    legend.position = "top",
    plot.background = element_rect(fill = "white", colour = NA),
    plot.title      = element_text(hjust = 0.5, margin = margin(t = 0, r = 0, b = 5, l = 0)),
    #axis.text.x     = element_text(angle = 45, vjust = 1, hjust = 1),
    panel.spacing.x = unit(2, "lines"),
    panel.spacing.y = unit(1, "lines"),
    strip.text.x    = element_text(margin = margin(t = 0, r = 0, b = 10, l = 0)),
    plot.subtitle   = ggtext::element_markdown(hjust = 0.5)
  )

O_diag_box.gg
ggsave("./sim/figures/O_diag_box.png", O_diag_box.gg, height = height, width = width, dpi = 300)

# # O-diag ridges
# print("Omega0 diag ridges gg")
# 
# O_diag_ridges.gg <-
#   sim_res.df %>%
#   dplyr::select(model:seed, O_diag) %>%
#   unnest(O_diag) %>%
#   dplyr::select(-seed) %>%
#   mutate(
#     setting = as.factor(str_c("Setting ", setting)) %>%
#       fct_reorder(setting, .desc = FALSE)
#   ) %>%
#   dplyr::select(model, setting, mcc, accuracy, kappa) %>%
#   pivot_longer(
#     cols = mcc:kappa,
#     names_to = "type",
#     values_to = "stat"
#   ) %>%
#   mutate(type = as.factor(type) %>% fct_reorder(stat, mean, .desc = FALSE)) %>%
#   left_join(sim_settings.df, ., by = c("setting")) %>%
#   mutate(
#     volumes = str_c(volumes, " Volumes"), 
#     subjects = str_c(subjects, " Subjects"),
#     outliers = ifelse(alpha_tau == 10 & lambda_2 == 0.25, 
#                       "Fewer Outliers",
#                       "More Outliers") %>%
#       as_factor() %>%
#       fct_relevel("Fewer Outliers")
#   ) %>%
#   group_by(type) %>%
#   mutate(model = as.factor(model) %>% fct_reorder(stat, median, .desc = FALSE)) %>%
#   ggplot(aes(x = stat, y = model, colour = type, fill = type)) +
#   ggridges::stat_density_ridges(quantile_lines = FALSE, colour = NA,
#                                 alpha = 0.4, size = 0.6, jittered_points = TRUE,
#                                 bandwith = 0.2) +
#   facet_wrap(~type, ncol = 3) +
#   scale_x_continuous(breaks = seq(0, 1, by = 1), minor_breaks = seq(0, 1, by = 1)) +
#   #geom_point(position = position_jitter(width = 0.1), alpha = 0.16) +
#   stat_summary(
#     geom = "point",
#     fun = "mean",
#  #   col = colour,
#     size = 6,
#     shape = "|",
#     alpha = 1
#   ) +
#   labs(
#     y = "Model",
#     x = "Statistic",
#     title = "Simulation Diagnostics: Omega_0"
#   ) +
#   scale_colour_viridis_d("Diagnostic") +
#   scale_fill_viridis_d("Diagnostic") +
#   facet_grid(subjects ~ volumes, scales = "free_y") +
#   theme_minimal() +
#   theme(
#     strip.placement  = "outside",
#     legend.position = "bottom",
#     plot.background = element_rect(fill = "white", colour = NA),
#     plot.title      = element_text(hjust = 0.5, margin = margin(t = 0, r = 0, b = 5, l = 0)),
#     #axis.text.x     = element_text(angle = 45, vjust = 1, hjust = 1),
#     panel.spacing.x = unit(0.5, "lines"),
#     panel.spacing.y = unit(2, "lines"),
#     strip.text.x    = element_text(margin = margin(t = 0, r = 0, b = 10, l = 0)),
#     plot.subtitle   = ggtext::element_markdown(hjust = 0.5)
#   )
# 
# O_diag_ridges.gg
# ggsave("./sim/figures/O_diag_ridges.png", O_diag_ridges.gg, height = 6, width = 8)

#stop("Done for now")
#################################################################################
#OmegaK Norm
# K-norm diagnostics
print("OmegaK Norms")

k_norm_box.gg <-
  sim_res.df %>%
  dplyr::select(model:seed, norm_k) %>%
  unnest(norm_k) %>%
  mutate(setting = str_c("Setting ", setting)) %>%
#  dplyr::select(-seed) %>%
  pivot_longer(
    cols = L1:Spectral,
    names_to = "type",
    values_to = "norm"
  ) %>%
  mutate(type = as.factor(str_remove(type, str_c("_", norm))) %>% fct_relevel(norm_order)) %>%
  mutate(model = as.factor(model) %>% fct_relevel(mod_order)) %>%
  left_join(sim_settings.df, ., by = c("setting")) %>%
  mutate(
    volumes = str_c(volumes, " Volumes"), 
    subjects = str_c(subjects, " Subjects"),
    outliers = ifelse(alpha_tau == 10 & lambda_2 == 0.25, 
                      "Fewer Outliers",
                      "More Outliers") %>%
               as_factor() %>%
               fct_relevel("Fewer Outliers")
  ) %>%
  dplyr::select(subjects, volumes, outliers, model, subject, type, seed, norm) %>%
  group_by(subjects, volumes, outliers, model, seed, type) %>%
  summarise(norm = mean(norm)) %>% 
  ungroup() %>%
  ggplot(aes(x = type, y = norm, colour = model, fill = model)) +
  geom_boxplot(colour = "grey80", alpha = 0.72, size = 0.4) +
  labs(
    x = "Precision Norm Difference",
    y = "Avg. Subject-level ||Truth - Estimate||",
  #  title = "Simulation Norm Difference Diagnostics: Omega_k",
  #  subtitle = sprintf("Alpha = %g, Lambda2 = %g", alpha, lambda),
  #  caption = "| denotes the average norm."
  ) +
  scale_colour_viridis_d("Model") +
  scale_fill_viridis_d("Model", option = "magma") +
  facet_grid(subjects ~ outliers + volumes, scales = "fixed") +
  scale_y_continuous(limits = c(0, 2)) +
  theme_minimal() +
  theme(
    strip.placement  = "outside",
    legend.position = "top",
    plot.background = element_rect(fill = "white", colour = NA),
    plot.title      = element_text(hjust = 0.5, margin = margin(t = 0, r = 0, b = 5, l = 0)),
    #axis.text.x     = element_text(angle = 45, vjust = 1, hjust = 1),
    panel.spacing.x = unit(0.5, "lines"),
    panel.spacing.y = unit(2, "lines"),
    strip.text.x    = element_text(margin = margin(t = 0, r = 0, b = 10, l = 0)),
    plot.subtitle   = ggtext::element_markdown(hjust = 0.5)
  )
k_norm_box.gg
ggsave("./sim/figures/k_norm_box.png", k_norm_box.gg, height = height, width = width, dpi = 300)

#################################################################################
#OmegaK Diag
print("OmegaK diag")

k_diag_box.gg <-
  sim_res.df %>%
  dplyr::select(model:seed, k_diag) %>%
  unnest(k_diag) %>%
  mutate(setting = str_c("Setting ", setting)) %>%
  group_by(model, setting, seed) %>%
  dplyr::select(everything(), -c(FPR, FNR)) %>% 
  summarise(
    across(
      .cols = MCC:FDR,
      ~mean(.x, na.rm = TRUE),
      .names = "{.col}"
    ),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = MCC:FDR,
    names_to = "type",
    values_to = "stat"
  ) %>%
  mutate(type = as.factor(type) %>% fct_relevel(diag_order)) %>%
  left_join(sim_settings.df, ., by = c("setting")) %>%
  mutate(
    volumes = str_c(volumes, " Volumes"), 
    subjects = str_c(subjects, " Subjects"),
    outliers = ifelse(alpha_tau == 10 & lambda_2 == 0.25, 
                      "Fewer Outliers",
                      "More Outliers") %>%
      as_factor() %>%
      fct_relevel("Fewer Outliers")
  ) %>%
  mutate(model = as.factor(model) %>% fct_relevel(mod_order)) %>%
  ggplot(aes(x = type, y = stat, colour = model, fill = model)) +
  geom_boxplot(colour = "grey80", alpha = 0.8, size = 0.4, width = 0.76) +
  labs(
    x = "Binary Diagnostics of Group Graph/Network",
    y = "Avg. Score Across Subjects",
  #  title = "Simulation Diagnostics: Omega_k",
  #  subtitle = sprintf("Alpha = %g, Lambda2 = %g", alpha, lambda),
   # caption = "| denotes the average."
  ) +
  scale_colour_viridis_d("Model") +
  scale_fill_viridis_d("Model", option = "magma") +
  facet_grid(subjects ~ outliers + volumes, scales = "fixed") +
  scale_y_continuous(limits = c(0, 1)) +
  theme_minimal() +
  theme(
    strip.placement  = "outside",
    legend.position = "top",
    plot.background = element_rect(fill = "white", colour = NA),
    plot.title      = element_text(hjust = 0.5, margin = margin(t = 0, r = 0, b = 5, l = 0)),
    #axis.text.x     = element_text(angle = 45, vjust = 1, hjust = 1),
    panel.spacing.x = unit(1, "lines"),
    panel.spacing.y = unit(3, "lines"),
    strip.text.x    = element_text(margin = margin(t = 0, r = 0, b = 10, l = 0)),
    plot.subtitle   = ggtext::element_markdown(hjust = 0.5)
  )

  k_diag_box.gg
  ggsave("./sim/figures/k_diag_box.png", k_diag_box.gg, height = 8, width = 12, dpi = 300)


#############################################################################################
#K-tables
print("Sim tables (gt)...")

  #Prep data
O_norm.df <-
    sim_res.df %>%
    dplyr::select(model:seed, norm_0) %>%
    unnest(norm_0) %>%
    dplyr::select(-seed) %>%
    mutate(setting = str_c("Setting ", setting)) %>%
    group_by(model, setting) %>%
    summarise(
      across(
        .cols = L1:Spectral,
        ~mean(.x, na.rm = TRUE),
        .names = "{.col}"
      ),
      .groups = "drop"
    ) %>%
  #  mutate(type = as.factor(type) %>% fct_relevel(norm_order)) %>%
    left_join(sim_settings.df, ., by = c("setting")) %>%
    mutate(
      volumes = str_c(volumes, " Volumes"), 
      subjects = str_c(subjects, " Subjects"),
      outliers = ifelse(alpha_tau == 10 & lambda_2 == 0.25, 
                        "Fewer Outliers",
                        "More Outliers") %>%
        as_factor() %>%
        fct_relevel("Fewer Outliers")
    ) %>%
    mutate(model = as.factor(model) %>% fct_relevel(mod_order)) %>%
    dplyr::select(subjects, volumes, outliers, model, L1, Spectral, Frobenius) %>%
    rename(FL2 = Frobenius, SL2 = Spectral) %>%
    rename_with(.cols = where(is.numeric), ~str_c("G", .x))
    

O_diag.df <-
  sim_res.df %>%
  dplyr::select(model:seed, O_diag) %>%
  unnest(O_diag) %>%
  dplyr::select(-seed) %>%
  mutate(setting = str_c("Setting ", setting)) %>%
  dplyr::select(model, setting, FDR, TPR, TNR, MCC) %>%
  group_by(model, setting) %>%
  summarise(
    across(
      .cols = where(is.numeric),
      ~mean(.x, na.rm = TRUE),
      .names = "{.col}"
    ),
    .groups = "drop"
  ) %>%
  left_join(sim_settings.df, ., by = c("setting")) %>%
  mutate(
    volumes = str_c(volumes, " Volumes"), 
    subjects = str_c(subjects, " Subjects"),
    outliers = ifelse(alpha_tau == 10 & lambda_2 == 0.25, 
                      "Fewer Outliers",
                      "More Outliers") %>%
      as_factor() %>%
      fct_relevel("Fewer Outliers")
  ) %>%
  mutate(model = as.factor(model) %>% fct_relevel(mod_order)) %>%
  dplyr::select(subjects, volumes, outliers, model, MCC, TPR, TNR, FDR) %>%
#  rename(ACC = Accuracy) %>%
  rename_with(.cols = where(is.numeric), ~str_c("G", .x))

k_norm.df <-
  sim_res.df %>%
  dplyr::select(model:seed, norm_k) %>%
  unnest(norm_k) %>%
  mutate(setting = str_c("Setting ", setting)) %>%
  group_by(model, setting) %>%
  summarise(
    across(
      .cols = L1:Spectral,
      ~mean(.x, na.rm = TRUE),
      .names = "{.col}"
    ),
    .groups = "drop"
  ) %>%
  mutate(model = as.factor(model) %>% fct_relevel(mod_order)) %>%
  left_join(sim_settings.df, ., by = c("setting")) %>%
  mutate(
    volumes = str_c(volumes, " Volumes"), 
    subjects = str_c(subjects, " Subjects"),
    outliers = ifelse(alpha_tau == 10 & lambda_2 == 0.25, 
                      "Fewer Outliers",
                      "More Outliers") %>%
      as_factor() %>%
      fct_relevel("Fewer Outliers")
  ) %>%
  dplyr::select(subjects, volumes, outliers, model, L1, Spectral, Frobenius) %>%
  rename(FL2 = Frobenius, SL2 = Spectral) %>%
  rename_with(.cols = where(is.numeric), ~str_c("I", .x))

k_diag.df <-
  sim_res.df %>%
  dplyr::select(model:seed, k_diag) %>%
  unnest(k_diag) %>%
  mutate(setting = str_c("Setting ", setting)) %>%
  group_by(model, setting) %>%
  dplyr::select(everything(), -c(FPR, FNR)) %>% 
  summarise(
    across(
      .cols = MCC:FDR,
      ~mean(.x, na.rm = TRUE),
      .names = "{.col}"
    ),
    .groups = "drop"
  ) %>%
  left_join(sim_settings.df, ., by = c("setting")) %>%
  mutate(
    volumes = str_c(volumes, " Volumes"), 
    subjects = str_c(subjects, " Subjects"),
    outliers = ifelse(alpha_tau == 10 & lambda_2 == 0.25, 
                      "Fewer Outliers",
                      "More Outliers") %>%
      as_factor() %>%
      fct_relevel("Fewer Outliers")
  ) %>%
  mutate(model = as.factor(model) %>% fct_relevel(mod_order)) %>%
  dplyr::select(subjects, volumes, outliers, model, MCC, TPR, TNR, FDR) %>%
#  rename(ACC = Accuracy) %>%
  rename_with(.cols = where(is.numeric), ~str_c("I", .x))

#######################################################
#Table

diag.df <-
  left_join(
    O_diag.df,
    k_diag.df,
    by = c("subjects", "volumes", "outliers", "model")
  )

norm.df <-
  left_join(
    O_norm.df,
    k_norm.df,
    by = c("subjects", "volumes", "outliers", "model")
  )

table.df <-
  left_join(
    diag.df,
    norm.df,
    by = c("subjects", "volumes", "outliers", "model")
  )

max_cols <- c("MCC", "TPR", "TNR")
min_cols <- c("FDR", "L1", "FL2")
max_to_bold <- c(str_c("G", max_cols), str_c("I", max_cols))
min_to_bold <- c(str_c("G", min_cols), str_c("I", min_cols))

#20 Sub
table_20sub.gt <-
  table.df %>%
  filter(subjects == "20 Subjects") %>%
  dplyr::select(-c(subjects), -all_of(contains("SL2"))) %>%
  mutate(
    across(
      where(is.numeric),
      ~round(.x, 3)
    )
  ) %>%
  group_by(volumes, outliers) %>%
  mutate(
    across(
      .cols = all_of(min_to_bold),
      ~ifelse(.x == min(.x), str_c("<b>", .x, "</b>"), as.character(.x))
    )
  ) %>%
  mutate(
    across(
      .cols = all_of(max_to_bold),
      ~ifelse(.x == max(.x), str_c("<b>", .x, "</b>"), as.character(.x))
    )
  ) %>%
  mutate(
    across(
      .cols = -c(model),
      ~ifelse(.x == "<b>0</b>", "<b>0.000</b>", .x)
    )
  ) %>%
  mutate(
    across(
      .cols = -c(model),
      ~ifelse(.x == "<b>1</b>", "<b>1.000</b>", .x)
    )
  ) %>%
  rename(` ` = model) %>%
  gt() %>%
  fmt_markdown(columns = everything()) %>%
  # data_color(
  #   columns = GMCC:IFDR,
  #   colors = scales::col_numeric(
  #     palette = c("blue", "white", "red"),
  #     domain = c(-0.05, 1.05))
  # ) %>%
  # data_color(
  #   columns = GL1:GFL2,
  #   colors = scales::col_numeric(
  #     palette = c("blue", "white"),
  #     domain = c(0.3, 1))
  # ) %>%
  # data_color(
  #   columns = IL1:IFL2,
  #   colors = scales::col_numeric(
  #     palette = c("blue", "white"),
  #     domain = c(0.4, 2))
  # ) %>%
  tab_spanner(
    label = md("__Graph Network Diagnostics__"),
    columns = GMCC:IFDR
  ) %>%
  tab_spanner(
    label = md("__Precision Matrix Estimation__"),
    columns = GL1:IFL2
  ) %>%
  tab_footnote(
    footnote = gt::html("G: Group, I: Individual"),
    locations = cells_column_labels(
      columns = c(GMCC, GL1)
    )
  ) %>%
  tab_footnote(
    footnote = gt::html("MCC: Mathews Correlation Coef., TPR: True Positive, TNR: True Negative, FDR: False Discovery"),
    locations = cells_column_labels(
      columns = GMCC
    )
  ) %>%
  tab_footnote(
    footnote = gt::html("L1: L1 Matrix Norm, FL2: Frobenius L2 Matrix Norm"),
    locations = cells_column_labels(
      columns = GL1
    )
  )

table_20sub.gt
gtsave(table_20sub.gt, "./sim/figures/table_20.png", expand = 20)


#50 sub


table_50sub.gt <-
  table.df %>%
  filter(subjects == "50 Subjects") %>%
  dplyr::select(-c(subjects), -all_of(contains("SL2"))) %>%
  mutate(
    across(
      where(is.numeric),
      ~round(.x, 3)
    )
  ) %>%
  group_by(volumes, outliers) %>%
  mutate(
    across(
      .cols = all_of(min_to_bold),
      ~ifelse(.x == min(.x), str_c("<b>", .x, "</b>"), as.character(.x))
    )
  ) %>%
  mutate(
    across(
      .cols = all_of(max_to_bold),
      ~ifelse(.x == max(.x), str_c("<b>", .x, "</b>"), as.character(.x))
    )
  ) %>%
  mutate(
    across(
      .cols = -c(model),
      ~ifelse(.x == "<b>0</b>", "<b>0.000</b>", .x)
    )
  ) %>%
  mutate(
    across(
      .cols = -c(model),
      ~ifelse(.x == "<b>1</b>", "<b>1.000</b>", .x)
    )
  ) %>%
  rename(` ` = model) %>%
  gt() %>%
  fmt_markdown(columns = everything()) %>%
  # data_color(
  #   columns = GMCC:IFDR,
  #   colors = scales::col_numeric(
  #     palette = c("blue", "white", "red"),
  #     domain = c(-0.05, 1.05))
  # ) %>%
  # data_color(
  #   columns = GL1:GFL2,
  #   colors = scales::col_numeric(
  #     palette = c("blue", "white"),
  #     domain = c(0.25, 1))
  # ) %>%
  # data_color(
  #   columns = IL1:IFL2,
  #   colors = scales::col_numeric(
  #     palette = c("blue", "white"),
  #     domain = c(0.4, 2))
  # ) %>%
  tab_spanner(
    label = md("__Graph Network Diagnostics__"),
    columns = GMCC:IFDR
  ) %>%
  tab_spanner(
    label = md("__Precision Matrix Estimation__"),
    columns = GL1:IFL2
  ) %>%
  tab_footnote(
    footnote = gt::html("G: Group, I: Individual"),
    locations = cells_column_labels(
      columns = c(GMCC, GL1)
    )
  ) %>%
  tab_footnote(
    footnote = gt::html("MCC: Mathews Correlation Coef., TPR: True Positive, TNR: True Negative, FDR: False Discovery"),
    locations = cells_column_labels(
      columns = GMCC
    )
  ) %>%
  tab_footnote(
    footnote = gt::html("L1: L1 Matrix Norm, FL2: Frobenius L2 Matrix Norm"),
    locations = cells_column_labels(
      columns = GL1
    )
  )
table_50sub.gt
gtsave(table_50sub.gt, "./sim/figures/table_50.png", expand = 20)

