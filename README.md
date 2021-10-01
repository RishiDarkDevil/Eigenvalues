# Eigenvalue Calculation of a Matrix
Here I have written the code to calculate the Eigenvalues of a Matrix, using the QR-Algorithm.
I have mainly used 3 variants:
 - Theoretical QR-Algorithm: Which Fails for Large Matrices and is Slower.
 - Practical QR-Algorithm: Which is Better than the Theoretical QR-Algorithm, which involves reducing Matrix to Upper-Hessenberg(Tridiagonal for Symmetric Matrix) Form and then Applying QR-Algorithm.
 - Practical QR-Algorithm with Rayleigh Shifts: Which helps in convergence little faster.

I have added the code to calculate the Sample Variance-Covariance Matrix also.\
The Hessenberg Reduction of a Matrix and the QR-Decomposition of a Matrix are also been performed and the respective codes are present.\
For Both the above Reduction and Decompositions the Householder Matrices are used.

Usage:
 - Run the R script file using "Rscript Eigenvalues_of_Covariance_Matrix.R" from command line.
 - A GUI Prompt to select the .csv file containing data matrix is presented. Select the .csv file containing data matrix.
 - Prints out the number of steps after which convergence is obtained, along with the Eigenvalues and Eigenvectors.(If Convergence is not attained, then change the inputs of tol, n, method in the end of the file where function is called in the Script).

Tip:
 - Sometimes the path to the RScript.exe in the bin folder of the installed directory of R is not added to path. Add it if it is already not added in order to use the Rscript command.
 - Alternatively, and a better method is to run the code from within RStudio or other R IDEs.
