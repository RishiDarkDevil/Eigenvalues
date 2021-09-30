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

# Some Examples: 
# Ex-1
(mat <- matrix(c(1,2,3,4,2,3,4,5,3,4,5,6,4,5,6,7), 4, 4))
(hessen <- hessenberg(mat))
zapsmall(hessen$Q %*% t(hessen$Q)) # Verifies Orthogonal Matrix
t(hessen$Q) %*% hessen$H %*% hessen$Q # Verifies it actally worked
mat

# Ex-2
(mat <- matrix(c(1.1,1.1,1.1, 1.1, 2.7, 2.7, 1.1, 2.7, 5), 3, 3))
(hessen <- hessenberg(mat))
zapsmall(hessen$Q %*% t(hessen$Q))
t(hessen$Q) %*% hessen$H %*% hessen$Q
mat

# Ex-2
(mat <- matrix(sample(1:5, 9, replace = TRUE), 3, 3))
(hessen <- hessenberg(mat))
zapsmall(hessen$Q %*% t(hessen$Q))
t(hessen$Q) %*% hessen$H %*% hessen$Q
mat

# Ex-3
(mat <- matrix(sample(1:10, 100, replace = TRUE), 10, 10))
(hessen <- hessenberg(mat))
zapsmall(hessen$Q %*% t(hessen$Q))
t(hessen$Q) %*% hessen$H %*% hessen$Q
mat

# Ex-4
(mat <- matrix(sample(1:10, 16, replace = TRUE), 4, 4))
(hessen <- hessenberg(mat))
zapsmall(hessen$Q %*% t(hessen$Q))
t(hessen$Q) %*% hessen$H %*% hessen$Q
mat

# Ex-5
(mat <- matrix(rnorm(100), 10, 10))
(hessen <- hessenberg(mat))
zapsmall(hessen$Q %*% t(hessen$Q))
t(hessen$Q) %*% hessen$H %*% hessen$Q
mat

# Ex-6
(mat <- matrix(rnorm(100), 10, 10))
(mat <- (mat + t(mat)) / 2)
(hessen <- hessenberg(mat))
zapsmall(hessen$Q %*% t(hessen$Q))
zapsmall(hessen$H)
t(hessen$Q) %*% hessen$H %*% hessen$Q
mat

