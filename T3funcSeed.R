T3funcSeed <- function ( X , n , m , p , r1 , r2 , r3 , start , conv , StartSeed , StartA , StartB , StartC , StartH ) 
{
  # fits a Tucker3 model with orthonormal A, B and C
  
  # X (n x mp): matricized 3D-data array (frontal slices)
  # n: number of elements of first mode (the common mode: rows)
  # m: number of elements of second mode (columns)
  # p: number of elements of third mode (slabs)
  # r1,r2,r3: rank of the model Tucker3-model (rank1 is the rank of the PCA-model)
  # start: type of start for the analysis
  #        0 = rational start (from generalized eigenvalue decomposition: orthonormal component matrices)
  #        1 = random start (random orthonormalized component matrices)
  #        2 = user-specified initial component matrices (best to take orthonormal component matrices)
  # conv: convergence criterion (tolerance value)
  # StartSeed: seed to start the analysis
  #  [ StartA, StartB, StartC, StartH: eventually one may specify initial A,B, C and H matrices (necessary if start=2) ]
  #  [  StartA (n x r1): initial row component matrix                                                   ]
  #  [  StartB (m x r2): initial column component matrix                                                ]
  #  [  StartC (p x r3): initial slab component matrix                                                  ]
  #  [  StartH (r1 x r2r3): initial (matricized) core array (frontal slices)                            ]
  
  set.seed( StartSeed ) # only for runif
  library( ThreeWay )
  
  checkinput = 1
  
  if ( (start != 0) & (start != 1) & (start != 2) )
  {
    cat(" ",fill=TRUE)
    cat("start should be 0 (rational start), 1 (random start) or 2 (user specified start)",fill=TRUE)
    cat(" ",fill=TRUE)
    checkinput=0
  }
  
  if ( length( dim(X) ) != 2 )
  {
    X = matrix( X , n , m*p ) # P slices are concatenated horizontally
  }
  
  X = as.matrix(X)
  
  if (checkinput == 1)
  {
    cputime = system.time(
    {
      ss = sum( X ^ 2 )
      
      if (start == 0)
      {
        EIG = eigen(X %*% t(X))
        A = EIG$vectors[, 1:r1]
        Z = permnew(X, n, m, p)
        EIG = eigen(Z %*% t(Z))
        B = EIG$vectors[, 1:r2]
        Z = permnew(Z, m, p, n)
        EIG = eigen(Z %*% t(Z))
        C = EIG$vectors[, 1:r3]
      }
      else
      {
        if (start == 1)
        {
          if (n >= r1)
          {
            A = orth( matrix( runif( n * r1 , 0 , 1 ) , n , r1 ) - 0.5 )
          }
          else
          {
            A = orth( matrix( runif( r1 * r1 , 0 , 1 ) , r1 , r1 ) - 0.5 )
            A = A[1:n, ]
          }
          
          if (m >= r2)
          {
            B = orth( matrix( runif( m * r2 , 0 , 1 ) , m , r2 ) - 0.5 )
          }
          else
          {
            B = orth( matrix( runif( r2 * r2 , 0 , 1 ) , r2 , r2 ) - 0.5 )
            B = B[1:m, ]
          }
          
          if (p >= r3)
          {
            C = orth( matrix( runif( p * r3 , 0 , 1 ) , p , r3 ) - 0.5 )
          }
          else
          {
            C = orth( matrix( runif( r3 * r3 , 0 , 1 ) , r3 , r3 ) - 0.5 )
            C = C[1:p, ]
          }
        }
        else #user-specified value (should be orthonormal A, B and C)
        {
          if ( sum( ( ( t(StartA) %*% StartA ) - diag(r1) ) ^ 2 ) > .0000001 )
          {
            A = orth( StartA )
          }
          else
          {
            A = StartA
          }
          
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
      
      # Compute loss function value
      if (start != 2)
      {
        Z = permnew(t(A) %*% X, r1, m, p)
        Z = permnew(t(B) %*% Z, r2, p, r1)
        H = permnew(t(C) %*% Z, r3, r1, r2)
        f = ss - sum( H ^ 2 )
      }
      else
      {
        Z = B %*% permnew(A %*% H, n, r2, r3)
        Z = C %*% permnew(Z, m, r3, n)
        Z = permnew(Z, p, n, m)
        f = sum( (X - Z) ^ 2 )
      }
      
      iter = 0
      fold = f + (2 * conv * f)
      while (fold - f > f * conv)
      {
        iter = iter + 1
        fold = f
        
        Z = permnew(X, n, m, p)
        Z = permnew(t(B) %*% Z, r2, p, n)
        Z = permnew(t(C) %*% Z, r3, n, r2)
        A = qr.Q(qr(Z %*% (t(Z) %*% A)), complete = FALSE)
        
        Z = permnew(X, n, m, p)
        Z = permnew(Z, m, p, n)
        Z = permnew(t(C) %*% Z, r3, n, m)
        Z = permnew(t(A) %*% Z, r1, m, r3)
        B = qr.Q(qr(Z %*% (t(Z) %*% B)), complete = FALSE)
        
        Z = permnew(t(A) %*% X, r1, m, p)
        Z = permnew(t(B) %*% Z, r2, p, r1)
        C = qr.Q(qr(Z %*% (t(Z) %*% C)), complete = FALSE)
        
        Z = permnew(t(A) %*% X, r1, m, p)
        Z = permnew(t(B) %*% Z, r2, p, r1)
        H = permnew(t(C) %*% Z, r3, r1, r2)
        
        f = ss - sum( H ^ 2 )
      }
    })
    
    fp = 100 * (ss - f) / ss
    
    La = H %*% t(H)
    Y = permnew(H, r1, r2, r3)
    Lb = Y %*% t(Y)
    Y = permnew(Y, r2, r3, r1)
    Lc = Y %*% t(Y)
    
    out = list()
    out$A = A
    out$B = B
    out$C = C
    out$H = H
    out$f = f
    out$fp = fp
    out$iter = iter
    out$cputime = cputime[1]
    out$La = La
    out$Lb = Lb
    out$Lc = Lc
    return(out)
  }
  else
  {
    out = list()
  }
}