library(tidyverse)

rcm.lam <- list.files("./sim/sim_res/freq/model_results/", pattern = "rcmlam")
print(rcm.lam)
rcm.lam.list <- map_df(.x = str_c("./sim/sim_res/freq/model_results/", rcm.lam), ~read_rds(.x))
head(rcm.lam.list)

rcm.lam.df <- 
rcm.lam.list %>% 
mutate(
file = rcm.lam,
setting = str_split_fixed(file, "_", 3)[,2] %>% parse_number(),
seed = str_split_fixed(file, "_", 3)[,3] %>% parse_number()
) %>%
dplyr::select(-file) %>%
group_by(setting) %>%
summarise(
 across(
 all_of(contains("lam")),
 ~mean(.x)
 )
)

rcm.lam.df