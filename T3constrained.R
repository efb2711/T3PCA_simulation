T3constrained <- function( X , FixedA , n , m , p , r1 , r2 , r3 , start , conv = .0000001 , StartSeed , StartB , StartC , StartH )
{
  # fit a Tucker3 model (with orthonormal B and C) with the constraint that A should equal the given BestA (ideally FixedA is also orthonormal)
  
  # X (i x jk): 3D-data matrix (displayed in matrix form)
  # FixedA (i x r1): fixed A component matrix (for the row mode) [ideally A is orthonormal but is not necessary]
  # n: number of rows
  # m: number of columns
  # p: number of slices
  # r1: number of component for the row-mode A
  # r2: number of component for the column-mode B
  # r3: number of component for the slice-mode C
  # start: type of starting values for B and C
  #   0 = generalized eigenvalue decomposition
  #   1 = random starting point
  #   2 = user specified component matrices
  # conv: tolerance value (default: .0000001)
  # StartSeed: seed to start the analysis
  # [ optional: specify StartB, StartC and StartH ]
  #  StartB (j x r2): [optional] initial values for the component matrix for the column-mode B (should be orthonormal)
  #  StartC (k x r3): [optional] initial values for the component matrix for the slice-mode C (should be orthonormal)
  #  StartH (r1 x r2r3): [optional] initial values for the core array H
  
  set.seed( StartSeed )
  library( ThreeWay )
  
  checkinput = 1
  
  if ( sum( ( ( t(FixedA) %*% FixedA ) - diag(r1) ) ^ 2 ) > .0000001 )
  {
    cat(" ",fill=TRUE)
    cat("A is not orthonormal (which is not necessary but recommended)",fill=TRUE)
    cat(" ",fill=TRUE)
  }
  
  if ( (start != 0) & (start != 1) & (start != 2) )
  {
    cat(" ",fill=TRUE)
    cat("start should be 0 (rational start), 1 (random start) or 2 (user specified start)",fill=TRUE)
    cat(" ",fill=TRUE)
    checkinput=0
  }
  
  if (checkinput == 1)
  {
    if ( length( dim(X) ) != 2 )
    {
      X = matrix( X , n , m*p ) # P slices are concatenated horizontally
    }
    X = as.matrix( X )
    
    cputime = system.time(
    {
      ss = sum( X ^ 2 )
      
      # determine initial B and C
      if ( start == 0 )
      {
        Z = permnew( X , n , m , p )
        EIG = eigen( Z %*% t( Z ) )
        B = EIG$vectors[, 1:r2]
        Z = permnew( Z , m , p , n )
        EIG = eigen( Z %*% t( Z ) )
        C = EIG$vectors[, 1:r3]
      }
      else
      {
        if ( start == 1 )
        {
          if ( m >= r2 )
          {
            B = orth( matrix( runif( m * r2 , 0 , 1 ) , m , r2 ) - 0.5 )
          }
          else
          {
            B = orth( matrix( runif( r2 * r2 , 0 , 1 ) , r2 , r2 ) - 0.5 )
            B = B[1:m, ]
          }
          
          if ( p >= r3 )
          {
            C = orth( matrix( runif( p * r3 , 0 , 1 ) , p , r3 ) - 0.5 )
          }
          else
          {
            C = orth( matrix( runif( r3 * r3 , 0 , 1 ) , r3 , r3 ) - 0.5 )
            C = C[1:p, ]
          }
        }
        else
        {
          if ( sum( ( ( t(StartB) %*% StartB ) - diag(r2) ) ^ 2 ) > .0000001 )
          {
            B = orth( StartB )
          }
          else
          {
            B = StartB
          }
          
          if ( sum( ( ( t(StartC) %*% StartC ) - diag(r3) ) ^ 2 ) > .0000001 )
          {
            C = orth( StartC )
          }
          else
          {
            C = StartC
          }
          
          H = StartH
        }
      }
      
      # compute loss (and eventually initial H)
      if ( start != 2 )
      {
        Z = permnew( t( FixedA ) %*% X , r1 , m , p )
        Z = permnew( t( B ) %*% Z , r2 , p , r1 )
        H = permnew( t( C ) %*% Z , r3 , r1 , r2 )
        f = ss - sum( H ^ 2 )
      }
      else
      {
        Z = B %*% permnew( FixedA %*% H , n , r2 , r3 )
        Z = C %*% permnew( Z , m , r3 , n )
        Z = permnew( Z , p , n , m )
        f = sum( ( X - Z ) ^ 2 )
      }
      
      iter = 0
      fold = f + ( 2 * conv * f )
      while ( ( fold - f ) > ( f * conv ) )
      {
        iter = iter + 1
        fold = f
        
        Z = permnew( X , n , m , p )
        Z = permnew( Z , m , p , n )
        Z = permnew( t(C) %*% Z , r3 , n , m )
        Z = permnew( t(FixedA) %*% Z , r1 , m , r3 )
        B = qr.Q( qr( Z %*% ( t(Z) %*% B ) ) , complete = FALSE )
        rm(Z)
        
        Z = permnew( t(FixedA) %*% X , r1 , m , p )
        Z = permnew( t(B) %*% Z , r2 , p , r1 )
        C = qr.Q( qr( Z %*% ( t(Z) %*% C ) ) , complete = FALSE )
        rm(Z)
        
        Z = permnew( t(FixedA) %*% X , r1 , m , p )
        Z = permnew( t(B) %*% Z , r2 , p , r1 )
        H = permnew( t(C) %*% Z , r3 , r1 , r2 )
        rm(Z)
        
        f = ss - sum( H ^ 2 )
      }
    })
    
    fp = 100 * ( ( ss - f ) / ss )
    
    Out = list()
    Out$B = B
    Out$C = C
    Out$H = H
    Out$f = f
    Out$fp = fp
    Out$iter = iter
    Out$cputime = cputime[1]
    return(Out)
  }
}