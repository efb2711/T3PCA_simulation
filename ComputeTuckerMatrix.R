ComputeTuckerMatrix <- function ( X , Y )
  # Computes the Tucker congruence between two matrices
  # X (n x p): first matrix
  # Y (n x p): second matrix
{
  source("ComputePhiMatrix.R")
  PhiMatrix = ComputePhiMatrix( X , Y )
  Tucker = mean( abs( diag(PhiMatrix) ) )
  return(Tucker)
}