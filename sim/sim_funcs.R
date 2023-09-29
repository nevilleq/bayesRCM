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