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
Eff.QR <- function(mat, check_rank = FALSE){
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

# Some Example Cases
# Ex-1
(outMatrix <- matrix(data = 1:12,4,3))
Eff.QR(outMatrix)
Eff.QR(outMatrix)$Q %*% Eff.QR(outMatrix)$R

# Ex-2
(mat <- matrix(c(1,1,1,2,2,2,3,3,3), nrow = 3, ncol = 3))
Eff.QR(mat)
Eff.QR(mat)$Q %*% Eff.QR(mat)$R

# Ex-3
(mat <- matrix(c(1,1,1,7,-2,2,3,-2,2,8,3,2,5,1,5,4,4,0,3,4), nrow = 5, ncol = 4))
Eff.QR(mat)
Eff.QR(mat)$Q %*% Eff.QR(mat)$R