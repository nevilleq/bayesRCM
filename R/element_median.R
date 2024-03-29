#' Element-wise Median
#'
#' @param matrix_list A list of matrices of the same size
#'
#' @return matrix of the element-wise median of a list of matrices
#' @export
#'
#' @examples
elementwise_median <- function(matrix_list) {
  # Get the dimensions of the first matrix in the list
  matrix_dims <- dim(matrix_list[[1]])
  
  # Initialize a result matrix with NA values
  result_matrix <- matrix(NA, nrow = matrix_dims[1], ncol = matrix_dims[2])
  
  #Loop through elements i, j
  for (i in 1:matrix_dims[1]) {
    for (j in 1:matrix_dims[2]) {
      # Extract the corresponding element from each matrix and calculate the median
      element_values <- sapply(matrix_list, function(matrix) matrix[i, j])
      result_matrix[i, j] <- median(element_values)
    }
  }
  
  return(result_matrix)
}

#' Element-wise SD
#'
#' @param matrix_list A list of matrices of the same size
#'
#' @return matrix of the element-wise median of a list of matrices
#' @export
#'
#' @examples
elementwise_sd <- function(matrix_list) {
  # Get the dimensions of the first matrix in the list
  matrix_dims <- dim(matrix_list[[1]])
  
  # Initialize a result matrix with NA values
  result_matrix <- matrix(NA, nrow = matrix_dims[1], ncol = matrix_dims[2])
  
  #Loop through elements i, j
  for (i in 1:matrix_dims[1]) {
    for (j in 1:matrix_dims[2]) {
      # Extract the corresponding element from each matrix and calculate the median
      element_values <- sapply(matrix_list, function(matrix) matrix[i, j])
      result_matrix[i, j] <- sd(element_values)
    }
  }
  
  return(result_matrix)
}

#' Element-wise Variance
#'
#' @param matrix_list A list of matrices of the same size
#'
#' @return matrix of the element-wise median of a list of matrices
#' @export
#'
#' @examples
elementwise_var <- function(matrix_list) {
  # Get the dimensions of the first matrix in the list
  matrix_dims <- dim(matrix_list[[1]])
  
  # Initialize a result matrix with NA values
  result_matrix <- matrix(NA, nrow = matrix_dims[1], ncol = matrix_dims[2])
  
  #Loop through elements i, j
  for (i in 1:matrix_dims[1]) {
    for (j in 1:matrix_dims[2]) {
      # Extract the corresponding element from each matrix and calculate the median
      element_values <- sapply(matrix_list, function(matrix) matrix[i, j])
      result_matrix[i, j] <- sd(element_values)^2
    }
  }
  
  return(result_matrix)
}

#' Element-wise Mean
#'
#' @param matrix_list A list of matrices of the same size
#'
#' @return matrix of the element-wise median of a list of matrices
#' @export
#'
#' @examples
elementwise_mean <- function(matrix_list) {
  # Get the dimensions of the first matrix in the list
  matrix_dims <- dim(matrix_list[[1]])
  
  # Initialize a result matrix with NA values
  result_matrix <- matrix(NA, nrow = matrix_dims[1], ncol = matrix_dims[2])
  
  #Loop through elements i, j
  for (i in 1:matrix_dims[1]) {
    for (j in 1:matrix_dims[2]) {
      # Extract the corresponding element from each matrix and calculate the median
      element_values <- sapply(matrix_list, function(matrix) matrix[i, j])
      result_matrix[i, j] <- mean(element_values)
    }
  }
  
  return(result_matrix)
}

#' Element-wise SD
#'
#' @param matrix_list A list of matrices of the same size
#'
#' @return matrix of the element-wise median of a list of matrices
#' @export
#'
#' @examples
elementwise_sd <- function(matrix_list) {
  # Get the dimensions of the first matrix in the list
  matrix_dims <- dim(matrix_list[[1]])
  
  # Initialize a result matrix with NA values
  result_matrix <- matrix(NA, nrow = matrix_dims[1], ncol = matrix_dims[2])
  
  #Loop through elements i, j
  for (i in 1:matrix_dims[1]) {
    for (j in 1:matrix_dims[2]) {
      # Extract the corresponding element from each matrix and calculate the median
      element_values <- sapply(matrix_list, function(matrix) matrix[i, j])
      result_matrix[i, j] <- sd(element_values)
    }
  }
  
  return(result_matrix)
}

#' Element-wise Variance
#'
#' @param matrix_list A list of matrices of the same size
#'
#' @return matrix of the element-wise median of a list of matrices
#' @export
#'
#' @examples
elementwise_var <- function(matrix_list) {
  # Get the dimensions of the first matrix in the list
  matrix_dims <- dim(matrix_list[[1]])
  
  # Initialize a result matrix with NA values
  result_matrix <- matrix(NA, nrow = matrix_dims[1], ncol = matrix_dims[2])
  
  #Loop through elements i, j
  for (i in 1:matrix_dims[1]) {
    for (j in 1:matrix_dims[2]) {
      # Extract the corresponding element from each matrix and calculate the median
      element_values <- sapply(matrix_list, function(matrix) matrix[i, j])
      result_matrix[i, j] <- sd(element_values)^2
    }
  }
  
  return(result_matrix)
}














