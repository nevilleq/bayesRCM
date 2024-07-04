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
   
  #pred_list  <- sim_diag.df$adj_k[[1]]
  #truth_list <- Adj_k
  #Grab list of upper triangle vectors
  pred_up_tri   <- map(.x = pred_list, ~factor(.x[upper.tri(.x, diag = FALSE)], levels = c("FALSE", "TRUE")))
  truth_up_tri  <- map(.x = truth_list, ~factor(.x[upper.tri(.x, diag = FALSE)], levels = c("FALSE", "TRUE")))
  
  #Make table and get confusion matrix results & diagnostics
  table         <- map2(.x = pred_up_tri, .y = truth_up_tri, 
                        ~table(.x, .y))
  conf_list     <- map(.x = table, ~caret::confusionMatrix(.x, positive = "TRUE"))
  
  #Put together final diagnostic tibble rows = subjects columns = which sub & diagnostic stats
  diagnostic.df <- 
    map_df(.x = conf_list, ~as_tibble(as.data.frame(t(.x$overall)))) %>%
    janitor::clean_names() %>%
    rename(Accuracy = accuracy) %>%
    mutate(
      subject = str_c("Sub. ", 1:length(pred_list)), 
      MCC     = map2_dbl(.x = pred_up_tri, .y = truth_up_tri,
                         ~mltools::mcc(preds = .x, actuals = .y)),
      tp      = map_dbl(.x = conf_list, ~.x$table[2,2]),
      fp      = map_dbl(.x = conf_list, ~.x$table[2,1]),
      tn      = map_dbl(.x = conf_list, ~.x$table[1,1]),
      fn      = map_dbl(.x = conf_list, ~.x$table[1,2]),
      TPR     = map2_dbl(.x = tp, .y = fn, ~.x / (.x + .y)), #tp / (tp + fn),
      FPR     = map2_dbl(.x = fp, .y = tn, ~.x / (.x + .y)), #fp / (fp + tn),
      TNR     = map2_dbl(.x = tn, .y = fp, ~.x / (.x + .y)), #tn / (tn + fp),
      FNR     = map2_dbl(.x = fn, .y = tp, ~.x / (.x + .y)), #fn / (fn + tp),
      FDR     = map2_dbl(.x = fp, .y = tp, ~.x / (.x + .y))  #fp / (fp + tp)
    ) %>%
    dplyr::select(Accuracy, MCC, TPR:FDR)
  
  #Return diag .df to unnest
  return(diagnostic.df)
}

get_bin_diag_0 <- function(omega, truth) {
  
  #Upper triangle vec
  pred_up_tri   <- factor(omega[upper.tri(omega, diag = FALSE)], levels = c("FALSE", "TRUE"))
  truth_up_tri  <- factor(truth[upper.tri(truth, diag = FALSE)], levels = c("FALSE", "TRUE"))
  
  #Make table and get confusion matrix results & diagnostics
  table <- table(pred_up_tri, truth_up_tri)
  conf  <- caret::confusionMatrix(table, positive = "TRUE")
  
  #Put together final diagnostic tibble rows = subjects columns = which sub & diagnostic stats
  diagnostic.df <- 
    as_tibble(as.data.frame(t(conf$overall))) %>%
    janitor::clean_names() %>%
    rename(Accuracy = accuracy) %>%
    mutate(
      MCC = mltools::mcc(preds = pred_up_tri, actuals = truth_up_tri),
      tp  = conf$table[2,2],
      fp  = conf$table[2,1],
      tn  = conf$table[1,1],
      fn  = conf$table[1,2],
      TPR = tp / (tp + fn),
      FPR = fp / (fp + tn),
      TNR = tn / (tn + fp),
      FNR = fn / (fn + tp),
      FDR = fp / (fp + tp)
    ) %>%
    dplyr::select(Accuracy, MCC, TPR:FDR)
  
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