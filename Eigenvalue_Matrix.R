compute_covariance <- function(mat){
  cov_matrix <- mat
  for(i in 1:ncol(mat)){
    cov_matrix[, i] = cov_matrix[, i] - mean(cov_matrix[, i])
  }
  return((t(cov_matrix) %*% cov_matrix)/(nrow(mat)-1))
}

# Doing multiplication with Householder Matrix from the Right Side
hmult2 <- function(u, x) {x - 2*sum(x*u)*u}

# Computing H from the packed Hessenberg Decomposed Matrix
compute.H <- function(m, H.diag){
  for (i in 1:(ncol(m)-2)) {
    m[(i+1):nrow(m),i] <- 0
    m[(i+1),i] = H.diag[i]
  }
  return(list("H" = m))
}

# Computing Q from the packed Hessenberg Decomposed Matrix
compute.QH <- function(m){
  Q <- diag(ncol(m))
  for(i in 1:(ncol(m)-2)){
    for (j in 2:ncol(m)) {
      Q[(i+1):nrow(Q), j] <- hmult(m[(i+1):nrow(m), i], Q[(i+1):nrow(Q), j])
    }
  }
  return(list("Q"=Q))
}

# Calculates Hessenberg Reduction of a Matrix
hessenberg <- function(mat){
  H.diags = numeric((ncol(mat)-2))
  for (i in 1:(ncol(mat)-2)) {
    u = shaver(mat[(i+1):nrow(mat),i])
    for (j in i:ncol(mat)) {
      mat[(i+1):nrow(mat),j] <- hmult(u, mat[(i+1):nrow(mat),j])
    }
    for (j in 1:nrow(mat)) {
      mat[j, (i+1):ncol(mat)] <- hmult2(u, mat[j, (i+1):ncol(mat)])
    }
    H.diags[i] = mat[i+1,i]
    mat[(i+1):nrow(mat),i] = u
  }
  computeQ <- compute.QH(mat)
  return(list("Packed_Matrix" = mat, "H_Diag" = H.diags, "Q" = computeQ$Q, "H" = compute.H(mat, H.diags)$H))
}

# Making Unit vectors
unit <- function(v) {v/sqrt(sum(v*v))}

# Doing multiplication with Householder Matrix
hmult <- function(u, x) {x - 2*sum(u*x)*u}

# Creating shaver for shaving off the cols to concentrate norm on top 
shaver <- function(x){
  x[1] <- x[1] - sqrt(sum(x*x))
  if(sum(abs(x)) == 0){
    return(x)
  }
  unit(x)
}

# Computing R from the packed QR Decomposed Matrix
compute.R <- function(m,R.diag){
  for (i in 1:ncol(m)) {
    m[i:nrow(m),i] = numeric(nrow(m)-i+1)
    m[i,i] = R.diag[i]
  }
  return(list("R" = m))
}

# Computing Qx when the packed QR Decomposed Matrix is passed along with x vector
multiply.Q <- function(m,x){
  res.vect <- x
  ulist <- list()
  for (i in ncol(m):1) {
    u = m[i:nrow(m),i]
    v <- numeric(nrow(m))
    v[i:nrow(m)] = u
    name <- paste("u",i,":",sep = "")
    ulist[[name]] <- v
    res.vect <- hmult(v,res.vect)
  }
  return(list("Qx"=res.vect, "ulist"=ulist))
}

# Computing Q from the packed QR Decomposed Matrix - Not Needed Did extra
compute.Q <- function(m){
  e <- numeric(nrow(m))
  e[1] = 1
  A <- multiply.Q(m,e)
  Q <- A$Qx
  
  for (i in 2:nrow(m)){
    e <- numeric(nrow(m))
    e[i] = 1
    Q <- rbind(Q,multiply.Q(m,e)$Qx)
  }
  return(list("Q" = t(Q), "U" = rev(A$ulist)))
}

# Efficient QR decomposition
Eff.QR <- function(mat){
  R.diags = numeric(ncol(mat))
  
  for (i in 1:ncol(mat)) {
    if(i != nrow(mat)){
      u = shaver(mat[i:nrow(mat),i])
    }
    else{
      u = 1
    }
    for (j in i:ncol(mat)) {
      mat[i:nrow(mat),j] <- hmult(u, mat[i:nrow(mat),j])
    }
    R.diags[i] = mat[i,i]
    mat[i:nrow(mat),i] = u
  }
  computeQ <- compute.Q(mat)
  return(list("Packed_Matrix" = mat, "R_Diag" = R.diags, "Q" = computeQ$Q, "R" = compute.R(mat, R.diags)$R, "U" = computeQ$U))
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

# Finds Eigenvalues of a Matrix with QR-Algo with Shifts
eigenvalues.shift <- function(M, n = 1000, tol = 1e-8){ # P in case of symmetric M contains the eigenvectors 
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

# Finds Eigenvalues of a Matrix with QR-Algo without Shifts
eigenvalues.wo.shifts <- function(M, n = 1000, tol = 1e-8){ # P in case of symmetric M contains the eigenvectors 
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

# Simple QR-Algorithm
eigenvalues <- function(M, n = 1000, tol = 1e-8){ # P in case of symmetric M contains the eigenvectors 
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

eigenval <- function(dataset, method = "noshift", ...){
  tryCatch({
    if(method == "shift")
      eigenvalues.shift(dataset, ...)
    else if(method == "noshift")
      eigenvalues.wo.shifts(dataset, ...)
    else
      eigenvalues(dataset)}, error = function(e) print("Error!Try some Other Method/tol/n..."))
}

# If it's a csv file which contains the Data Matrix
dataset <- as.matrix(read.csv(file.choose(), header = FALSE))

# If it's already loaded in R(eg. as data) it should be in a form of a matrix, then you can simply run the bottom command
# dataset <- data

# Eigenvalue calculates the eigen value of a matrix
eigenval(dataset)
# eigen(dataset)

# Examples:
# eigenval(dataset, method = "noshift", tol = 1e-10, n = 1500)
# eigenval(dataset, method = "shift", tol = 1e-10, n = 1500)
# eigenval(dataset, method = "simple", tol = 1e-10, n = 1500)
