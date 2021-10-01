# Practical QR-Algorithm
# Solve Linear Upper Triangular System
Solve.Eigenvect <- function(T.mat, evals, P){
  eigenvects <- matrix(rep(NA, ncol(T.mat)*nrow(T.mat)), ncol = ncol(T.mat), nrow = nrow(T.mat))
  v <- rep(NA, ncol(T.mat))
  for (k in 2:ncol(T.mat)) {
    v <- rep(NA, length(v))
    if(k < ncol(T.mat))
      v[(k+1):length(v)] <- 0
    v[k] <- 1
    for(j in (k-1):1){
      v[j] <- -sum(T.mat[j,(j+1):k]*v[(j+1):k]) / (T.mat[j,j] - T.mat[k,k])
    }
    eigenvects[,k] <- v
    
  }
  eigenvects[,1] <- c(1, rep(0, ncol(T.mat)-1))
  eigenvects <- P %*% eigenvects
  colnames(eigenvects) <- evals
  return(eigenvects)
}

# Finds Eigenvalues of a Matrix
eigenvalues <- function(M, n = 1000, tol = 1e-3){ # P in case of symmetric M contains the eigenvectors 
  P <- diag(nrow(M))
  if(ncol(M) > 2){
    hess <- hessenberg(M)
    M <- hess$H
  }
  for(i in 1:n){
    QR <- Eff.QR(M)
    old_M <- diag(M)
    M <- QR$R %*% QR$Q
    new_M <- diag(M)
    P <- P %*% QR$Q
    if(sum((new_M - old_M)^2, na.rm = TRUE) <= tol)
      break
  }
  if(i < n){
    cat(paste("Convergence after", i, "steps\n"))
  }
  else{
    cat(paste("Max. iterations exceeded\n"))
  }
  eigenvals <- rep(NA, nrow(M))
  for(i in 1:nrow(M)){
    eigenvals[i] <- as.numeric(M[i,i])
  }
  names(eigenvals) <- NULL
  if(ncol(M) > 2){
    P <- t(hess$Q) %*% P
  }
  eigenvects <- Solve.Eigenvect(M, eigenvals, P)
  return(list("eigenvalues" = eigenvals, "eigenvectors" = eigenvects))
}

# Some Examples:
# Ex-1
(A <- matrix(c(6,-1,4,1), 2, 2))
eigenvalues(A, tol = 1e-20)

# Ex-2
(A <- matrix(sample(1:10, 9, replace = TRUE), 3, 3))
(mat <- (A + t(A)) / 2)
(e_vals <- eigenvalues(mat))
eigen(mat)

# Ex-3
X1 <- rnorm(100, 1, 2)
X2 <- 2*X1 + rnorm(100)
X3 <- rnorm(100, 10, 3)
X4 <- 3*X3 + 5*X2
test_data <- cbind(X1, X2, X3, X4)
(cov_matrix <- cov(test_data))
(e_vals <- eigenvalues(cov_matrix, tol = 1e-10))
eigen(cov_matrix)
cov_matrix %*% e_vals$eigenvectors[,1]
e_vals$eigenvalues[1]*e_vals$eigenvectors[,1]

# Ex-4
(mat <- matrix(rnorm(100), 10, 10))
(mat <- (mat + t(mat)) / 2)
(cov_matrix <- compute_covariance(mat))
cov(mat)
eigenvalues(cov_matrix, tol = 1e-10)
eigen(cov_matrix)

# Ex-5
(mat <- matrix(rnorm(2500), 50, 50))
(cov_matrix <- compute_covariance(mat))
eigenvalues(cov_matrix, tol = 1e-5)
eigen(cov_matrix)