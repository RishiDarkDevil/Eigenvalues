# Eigenvalue Calculation of a Matrix
Here I have written the code to calculate the Eigenvalues of a Matrix, using the QR-Algorithm.
I have mainly used 3 variants:
 - Theoretical QR-Algorithm: Which Fails for Large Matrices and is Slower.
 - Practical QR-Algorithm: Which is Better than the Theoretical QR-Algorithm, which involves reducing Matrix to Upper-Hessenberg(Tridiagonal for Symmetric Matrix) Form and then Applying QR-Algorithm.
 - Practical QR-Algorithm with Rayleigh Shifts: Which helps in convergence little faster.

I have added the code to calculate the Sample Variance-Covariance Matrix also.\
The Hessenberg Reduction of a Matrix and the QR-Decomposition of a Matrix are also been performed and the respective codes are present.\
For Both the above Reduction and Decompositions the Householder Matrices are used.
