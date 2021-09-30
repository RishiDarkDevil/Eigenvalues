compute_covariance <- function(mat){
  cov_matrix <- mat
  for(i in 1:ncol(mat)){
    cov_matrix[, i] = cov_matrix[, i] - mean(cov_matrix[, i])
  }
  return((t(cov_matrix) %*% cov_matrix)/(nrow(mat)-1))
}
