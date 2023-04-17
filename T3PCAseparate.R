T3PCAseparate <- function( X3D , Y , r1 , r2 , r3 , conv , OriginalAlfa , AlternativeLossF , nRuns , StartSeed )
{
  # estimates a T3 for X3D (with orthonormal A, B1 and C) and a PCA for Y (orthonormal A and orthogonal B2)
  
  # INPUT
  #    X3D (n x m x p): 3D-data array
  #    Y (n x q): 2D-data matrix
  #    r1,r2,r3: rank of the model Tucker3-model (rank1 is the rank of the PCA-model)
  #    conv: value for convergence (tolerance value)
  #    OriginalAlfa (0-1): weight for the 3D-block (1-alfa is weight for the
  #         2D-block; .50 implies the unweighted version) between 0 and 1
  #         alfa = 1 (only take 3D-block into account)
  #         alfa = .50 (regular sum of both parts of the loss function)
  #         alfa = 0 (only take 2D-block into account)
  #         0 and 1 are not possible as input values
  #    AlternativeLossF: using the alternative loss function?
  #       0 = no (use original loss function: weighted SSQ; weighted met alfa)
  #       1 = yes (use weighted loss function with scaled SSQ: scaled by the SSQ in X and y )
  #    nRuns: number of runs
  #    StartSeed: seed to start the analysis
  #
  # OUTPUT
  #   Sol: structure for the solution containing the following information
  #        - Info
  #            * AnalysisType
  #            * nRows
  #            * nColumns3D
  #            * nSlices
  #            * nColumns2D
  #            * RankVector
  #            * TolPercentage
  #            * OriginalAlfa
  #            * Alfa
  #            * AlternativeLossF
  #            * nRuns
  #            * StartSeed
  #        - A1
  #        - B1
  #        - C
  #        - H
  #        - A2
  #        - B2
  #        - LossWeighted
  #        - LossUnweighted
  #        - FitPercentage
  #        - FitPercentage3D
  #        - nIter
  #        - FitValues
  #        - nIterValues
  #        - La
  #        - Lb
  #        - Lc
  #        - Fit3D
  #        - Fit2D
  #        - Fit
  #        - Fit3Dsize
  #        - Fit2Dsize
  #        - Fitsize
  #        - CpuTime
  #        - TimeSeconds
  
  set.seed( StartSeed )
  library( ThreeWay )
  source("T3funcSeed.R")
  
  AnalysisSeeds = GenerateSeed( nRuns+1 , StartSeed )
  
  temp = dim( X3D )
  i = temp[1]
  j = temp[2]
  k = temp[3]
  rm( temp )
  
  temp = dim( Y )
  l = temp[2]
  rm( temp )
  
  X = matrix( X3D , i , j*k ) # P slices are concatenated horizontally
  ssq3D = sum( X ^ 2 )
  ssq2D = sum( Y ^ 2 )
  
  if ( AlternativeLossF == 1 )
  {
    Alfa = ( OriginalAlfa * ssq2D ) / ( ( OriginalAlfa * ssq2D ) + ( ( 1 - OriginalAlfa ) * ssq3D ) )
  }
  else
  {
    Alfa = OriginalAlfa
  }
  
  FitValues = matrix( 0 , 1 , nRuns+1 )
  nIterValues = matrix( 0 , 1 , nRuns+1 )
  
  cputime = system.time(
  {
    # Tucker3 on X
    BestSol3D = T3funcSeed( X , i , j , k , r1 , r2 , r3 , 0 , .0000001 , AnalysisSeeds[1] )
    FitValues[1] = BestSol3D$f
    nIterValues[1] = BestSol3D$iter
    for ( runtel in 1:nRuns )
    {
      tempSol3D = T3funcSeed( X , i , j , k , r1 , r2 , r3 , 1 , .0000001 , AnalysisSeeds[ runtel+1 ] )
      FitValues[ runtel+1 ] = tempSol3D$f
      nIterValues[ runtel+1 ] = tempSol3D$iter
      if( tempSol3D$f < BestSol3D$f )
      {
        BestSol3D = tempSol3D
      }
      rm( tempSol3D )
    }
    BestA1 = BestSol3D$A #orthonormal
    BestB1 = BestSol3D$B #orthonormal
    BestC = BestSol3D$C #orthonorma
    BestH = BestSol3D$H
    BestIter = BestSol3D$iter
    BestFitPercentage3D = BestSol3D$fp
    
    # PCA on Y
    tempSol2D = svd( Y , r1 , r1 )
    EigValues = tempSol2D$d
    BestA2 = tempSol2D$u #orthonormal
    BestB2 = tempSol2D$v %*% diag( EigValues[1:r1] ) # geeft orthogonal B2 (niet zo voor SegmentedType == 1)
    rm( tempSol2D , EigValues )    
    Model3D = BestA1 %*% BestH %*% kronecker( t( BestC ) , t( BestB1 ) )
    Model2D = BestA2 %*% t( BestB2 )
    Loss3D = sum( ( Model3D - X ) ^ 2 ) 
    Loss2D = sum( ( Model2D - Y ) ^ 2 )
    FitPercentage = ( ( Alfa * ( sum(BestH^2) / ssq3D ) ) + ( ( 1 - Alfa ) * ( sum(Model2D^2) / ssq2D ) ) ) * 100
  })
  
  # compute intrinsic eigenvalues
  La = BestH %*% t( BestH )
  J = permnew( BestH , r1 , r2 , r3 )
  Lb = J %*% t( J )
  J = permnew( J , r2 , r3 , r1 )
  Lc = J %*% t( J )
  
  Out = list()
  Out$Info = list()
  Out$Info$AnalysisType = "separate"
  Out$Info$nRows = i
  Out$Info$nColumns3D = j
  Out$Info$nSlices = k
  Out$Info$nColumns2D = l
  Out$Info$RankVector = cbind( r1 , r2 , r3 )
  Out$Info$TolPercentage = conv
  Out$Info$OriginalAlfa = OriginalAlfa
  Out$Info$Alfa = Alfa
  Out$Info$AlternativeLossF = AlternativeLossF
  Out$Info$nRuns = nRuns
  Out$Info$StartSeed = StartSeed
  Out$A1 = BestA1
  Out$B1 = BestB1
  Out$C = BestC
  Out$H = BestH
  Out$A2 = BestA2
  Out$B2 = BestB2
  Out$LossWeighted = ( Alfa * Loss3D ) + ( ( 1 - Alfa ) * Loss2D )
  Out$LossUnweighted = Loss3D + Loss2D
  Out$FitPercentage = FitPercentage
  Out$FitPercentage3D = BestFitPercentage3D
  Out$nIter = BestIter
  Out$FitValues = FitValues
  Out$nIterValues = nIterValues
  Out$La = La
  Out$Lb = Lb
  Out$Lc = Lc
  Out$Fit3D = Loss3D
  Out$Fit2D = Loss2D
  Out$Fit = ( Alfa * Loss3D ) + ( ( 1 - Alfa ) * Loss2D )
  Out$Fit3Dsize = Loss3D / ( i * j * k )
  Out$Fit2Dsize = Loss2D / ( i * l )
  Out$Fitsize = Out$Fit / ( ( i * j * k ) + ( i * l ) )
  Out$CpuTime = cputime[1]
  Out$TimeSeconds = round( cputime[1] , 2 )
  return( Out )
}