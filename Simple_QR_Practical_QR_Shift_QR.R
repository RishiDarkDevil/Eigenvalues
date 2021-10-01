# Simple QR-Algorithm
eigenvalues <- function(M, n = 1000, tol = 1e-3){ # P in case of symmetric M contains the eigenvectors 
  P <- diag(nrow(M))
  for(i in 1:n){
    QR <- Eff.QR(M)
    old_M <- diag(M)
    M <- QR$R %*% QR$Q
    P <- P %*% QR$Q
    new_M <- diag(M)
    if(max(abs(diag(M[-1,-ncol(M)]))) <= tol)
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
  return(list("eigenvalues" = eigenvals, "P" = P))
}


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
    if(max(abs(diag(M[-1,-ncol(M)]))) <= tol)
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


# Practical QR-Algorithm with Shifts
# Solve Linear Upper Triangular System
Solve.Eigenvect <- function(T.mat, evals, P){
  eigenvects <- matrix(rep(NA, ncol(T.mat)*nrow(T.mat)), ncol = ncol(T.mat), nrow = nrow(T.mat))
  v <- rep(NA, ncol(T.mat))
  for (k in 2:ncol(T.mat)) {
    v <- rep(NA, length(v))
    if(k != ncol(T.mat))
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
    shift <- M[nrow(M), ncol(M)]
    QR <- Eff.QR(M - diag(shift, nrow(M)))
    old_M <- diag(M)
    M <- QR$R %*% QR$Q + diag(shift, nrow(M))
    new_M <- diag(M)
    P <- P %*% QR$Q
    if(max(abs(diag(M[-1,-ncol(M)]))) <= tol)
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

