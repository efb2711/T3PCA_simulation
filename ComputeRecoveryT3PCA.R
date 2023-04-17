ComputeRecoveryT3PCA <- function( TrueSol , Sol , SeparateAnalysis , ProxyFit )
{
  #compute congruences between true (without extra component) and obtained T3PCA-solution
  
  # SeparateAnalysis: Sol obtained from the separate analysis strategy?
  #    0 = no (Sol from simultaneous or segmented strategy)
  #    1 = yes (Sol from separate analysis)
  # Extracomp: solution with extra component
  #            0 = no (no extra component)
  #            1 = yes
  
  library( MASS )
  source( "procr.R" )
  source( "ComputeTuckerMatrix.R" )
  
  checkinput = 1
  
  if ( ( SeparateAnalysis != 0 ) & ( SeparateAnalysis != 1 ) )
  {
    cat(" ",fill=TRUE)
    cat("SeparateAnalysis should be 0 (simultaneous/segmented) or 1 (separate)", fill=TRUE )
    cat(" ",fill=TRUE)
    checkinput=0
  }
  
  Out = list()
  
  if ( checkinput == 1 )
  {
    #Comparing badness-of-fit and badness-of-data
    Out$BOFBOD = Sol$Fit - TrueSol$BOD
    Out$BOFBOD3D = Sol$Fit3D - TrueSol$BOD3D
    Out$BOFBOD2D = Sol$Fit2D - TrueSol$BOD2D
    Out$BOFBODsize = Sol$Fitsize - TrueSol$BODsize
    Out$BOFBOD3Dsize = Sol$Fit3Dsize - TrueSol$BOD3Dsize
    Out$BOFBOD2Dsize = Sol$Fit2Dsize - TrueSol$BOD2Dsize
    
    Out$BODreached = Sol$Fit <= ProxyFit
    Out$BODdiff = Sol$Fit - ProxyFit
    Out$nRunsOptimalFit = sum( Sol$FitValues == min( Sol$FitValues ) )
    
    #Recovery of M (BOR)
    r1 = Sol$Info$RankVector[1]
    r2 = Sol$Info$RankVector[2]
    r3 = Sol$Info$RankVector[3]
    
    if ( SeparateAnalysis == 1 )
    {
      Model3D = Sol$A1 %*% matrix( Sol$H , r1 , r2*r3 ) %*% kronecker( t(Sol$C) , t(Sol$B1) )
      Model2D = Sol$A2 %*% t( Sol$B2 )
    }
    else
    {
      Model3D = Sol$A %*% matrix( Sol$H , r1 , r2*r3 ) %*% kronecker( t(Sol$C) , t(Sol$B1) )
      Model2D = Sol$A %*% t( Sol$B2 )
    }
    
    Out$Recov = 0
    Out$Recov3D = sum( ( Model3D - TrueSol$TrueWay3DMatrixNoExtraComp ) ^ 2 )
    Out$Recov2D = sum( ( Model2D - TrueSol$TrueWay2DnoExtraComp ) ^ 2 )
    Out$Recov = ( Sol$Info$Alfa * Out$Recov3D ) + ( ( 1 - Sol$Info$Alfa ) * Out$Recov2D )
    Out$Recovsize = Out$Recov / ( ( Sol$Info$nRows * Sol$Info$nColumns3D * Sol$Info$nSlices ) + ( Sol$Info$nRows * Sol$Info$nColumns2D  ) )
    Out$Recov3Dsize = Out$Recov3D / ( Sol$Info$nRows * Sol$Info$nColumns3D * Sol$Info$nSlices )
    Out$Recov2Dsize = Out$Recov2D / ( Sol$Info$nRows * Sol$Info$nColumns2D  )
    
    Out$RecovExtra = 0
    Out$Recov3DExtra = sum( ( Model3D - TrueSol$TrueWay3DMatrix ) ^ 2 )
    Out$Recov2DExtra = sum( ( Model2D - TrueSol$TrueWay2D ) ^ 2 )
    Out$RecovExtra = ( Sol$Info$Alfa * Out$Recov3DExtra ) + ( ( 1 - Sol$Info$Alfa ) * Out$Recov2DExtra )
    Out$RecovsizeExtra = Out$RecovExtra / ( ( Sol$Info$nRows * Sol$Info$nColumns3D * Sol$Info$nSlices ) + ( Sol$Info$nRows * Sol$Info$nColumns2D  ) )
    Out$Recov3DsizeExtra = Out$Recov3DExtra / ( Sol$Info$nRows * Sol$Info$nColumns3D * Sol$Info$nSlices )
    Out$Recov2DsizeExtra = Out$Recov2DExtra / ( Sol$Info$nRows * Sol$Info$nColumns2D  )
    
    if ( SeparateAnalysis == 1 )
    {
      #Recovery measures for component matrices
      optrotationA1 = procr( Sol$A1 , TrueSol$A[,1:r1] )
      optrotationB1 = procr( Sol$B1 , TrueSol$B1 )
      optrotationC = procr( Sol$C , TrueSol$C )
      #counterrotation for the core (based on optimal rotations for A, B and C)
      InvRotatA1 = ginv( optrotationA1$RotationMatrix )
      rotatedCore = InvRotatA1 %*% Sol$H %*% kronecker( t(ginv( optrotationC$RotationMatrix )) , t(ginv( optrotationB1$RotationMatrix )) )
      #counterrotation for E (based on optimal rotation for A)
      optrotationA2 = procr( Sol$A2 , TrueSol$A[,1:r1] )
      InvRotatA2 = ginv( optrotationA2$RotationMatrix )
      rotationB2 = t( InvRotatA2 %*% t( Sol$B2 ) )   # moet dit niet Sol$B2 zijn ???????????????????????
      optrotationB2 = procr( Sol$B2 , TrueSol$B2[,1:r1] )
      InvRotatB2 = ginv( optrotationB2$RotationMatrix )
      rotationA2 = Sol$A2 %*% InvRotatB2
    }
    else
    {
      #Recovery measures for component matrices
      optrotationA = procr( Sol$A , TrueSol$A[,1:r1] )
      optrotationB1 = procr( Sol$B1 , TrueSol$B1 )
      optrotationC = procr( Sol$C , TrueSol$C )
      #counterrotation for the core (based on optimal rotations for A, B and C)
      InvRotatA = ginv( optrotationA$RotationMatrix )
      rotatedCore = InvRotatA %*% Sol$H %*% kronecker( t(ginv( optrotationC$RotationMatrix )) , t(ginv( optrotationB1$RotationMatrix )) )
      #counterrotation for E (based on optimal rotation for A)
      rotationB2 = t( InvRotatA %*% t( Sol$B2 ) )   # moet dit niet Sol$B2 zijn ???????????????????????
      optrotationB2 = procr( Sol$B2 , TrueSol$B2[,1:r1] )  #to find optimal E rotation
      rotationA = Sol$A %*% ginv( optrotationB2$RotationMatrix ) #although this is not very relevant
    }
    
    #compute Tucker congruence coefficients
    if ( SeparateAnalysis == 1 )
    {
      Out$TuckerA1 = optrotationA1$TuckerCongruence
    }
    else
    {
      Out$TuckerA = optrotationA$TuckerCongruence
    }
    Out$TuckerB1 = optrotationB1$TuckerCongruence
    Out$TuckerC = optrotationC$TuckerCongruence
    Out$TuckerCore = sum( ( rotatedCore - TrueSol$CoreMatrix[1:r1,] ) ^ 2 )
    Out$TuckerCoreSize = Out$TuckerCore / ( r1 * r2 * r3 )
    
    if ( SeparateAnalysis == 1 )
    {
      Out$TuckerA2optimal = optrotationA2$TuckerCongruence
      Out$TuckerB2 = ComputeTuckerMatrix( rotationB2 , TrueSol$B2[,1:r1] )
      #Out$TuckerB2size = Out$TuckerB2 / ( Sol$Info$nColumns2D * r1 )
      Out$TuckerB2optimal = optrotationB2$TuckerCongruence
      Out$TuckerA2 = ComputeTuckerMatrix( rotationA2 , TrueSol$A2[,1:r1] )
      #Out$TuckerA2size = Out$TuckerA2 / ( Sol$Info$nRows * r1 )
    }
    else
    {
      Out$TuckerB2 = ComputeTuckerMatrix( rotationB2 , TrueSol$B2[,1:r1] )
      #Out$TuckerB2size = Out$TuckerB2 / ( Sol$Info$nRows * Sol$Info$nColumns2D )
      Out$TuckerB2optimal = optrotationB2$TuckerCongruence
    }
  }
  
  return(Out)
}