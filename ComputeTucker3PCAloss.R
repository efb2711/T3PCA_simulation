ComputeTucker3PCAloss <- function ( A , B1 , C , Hmatrix , B2 , Xdata , Ydata , Alfa )
{
  # A (i x r1): component matrix for 1st mode of 3D (common mode)
  # B1 (j x r2): component matrix for 2nd mode of 3D
  # C (k x r3): component matrix for 3th mode of 3D
  # Hmatrix (r1 x r2r3): CoreArray (in matrix form)
  # B2 (l x r1): component matrix for 2nd mode of 2D
  # Xdata (i x jk): 3D data matrix
  # Ydata (i x l): 2D data matrix
  # Alfa (0-1): weight for the 3D part
  
  Model3D = A %*% Hmatrix %*% kronecker( t(C) , t(B1) )
  Model2D = A %*% t( B2 )
  Loss = ( Alfa * sum( ( Xdata - Model3D ) ^ 2 ) ) + ( ( 1 - Alfa ) * sum( ( Ydata - Model2D ) ^ 2 ) )
  return( Loss )
}