GenerateTucker3_PCA <- function( i , j , k , l , rank1 , rank2 , rank3 , noise , ExplVarVector , OriginalAlfa , AlternativeLossF , nExtraComponents , GenerateSeed )
{
  # generates a Tucker3_PCA true model and data
  
  # A1, B1, C en A2 zijn orthonormal en B2 orthogonal (A1 and A2 have 'rank1' components in common and 'nExtraComponents' distinctive components )
  
  # als A (ook A1 en A2), B1 en C orthonormal zijn dan klopt manipulatie van explained variance
  # klopt manipulatie van explained variance per component when B1, C of B2 niet orthonormal/gonal zijn (dit wordt niet geïmplementeerd) ??
 
  #    i: number of elements of first mode (the common mode: rows)
  #    j: number of elements of second mode (columns)
  #    k: number of elements of third mode (slabs)
  #    l: number of covariates
  #    rank1,rank2,rank3: rank of the model Tucker3-model (rank1 is the rank of the PCA-model)
  #    noise (0-1): amount of noise added to the model
  #    ExplVarVector(1 x rank1 + 1): elements denote percentage of explained variance per component for the COMMON mode (should sum to one)
  #         [[ percentage of explained variance is the ssq of the model with only this component divided by the ssq of the total model (with all components) ]]
  #    OriginalAlfa: weight for the 3D-block (1-alfa is weight for the
  #         2D-block; .50 implies the unweighted version) between 0 and 1
  #         alfa = 1 (only take 3D-block into account)
  #         alfa = .50 (regular sum of both parts of the loss function)
  #         alfa = 0 (only take 2D-block into account)
  #         0 and 1 are not possible as input values
  #    AlternativeLossF: using the alternative loss function?
  #       0 = no (use original loss function: weighted SSQ)
  #       1 = yes (use weighted loss function with scaled SSQ: scaled by the SSQ in X and y )
  #    nExtraComponents: number of extra components for the common mode
  #    GenerateSeed: seed for the generation
  #
  # OUTPUT
  #   Sol: structure for the solution containing the following information
  #        - nRow
  #        - nCol1
  #        - nSlice
  #        - nCol2
  #        - RankVector
  #        - Error
  #        - ExplVarVector
  #        - OriginalAlfa
  #        - AlternativeLossF
  #        - nExtraComponents
  #        - GenerateSeed
  #        - A
  #        - B1
  #        - C
  #        - B2orth
  #        - B2
  #        - A1
  #        - A2
  #        - Core
  #        - CoreMatrix
  #        - TrueWay3DMatrix
  #        - TrueWay3DMatrixNoExtraComp
  #        - TrueWay3D
  #        - TrueWay3DnoExtraComp
  #        - TrueWay2D
  #        - TrueWay2DnoExtraComp
  #        - Way3D
  #        - Way3DMatrix
  #        - Way2D
  #        - BOD3D
  #        - BOD2D
  #        - BOD3Dsize
  #        - BOD2Dsize
  #        - BOD
  #        - BODsize
  
  library( ThreeWay , MASS )
  
  set.seed( GenerateSeed ) # rnorm and runif used
  
  # check input  
  checkinput = 1
  
  if( rank1 > min(i,j*k,l) )
  {
    cat(" ",fill=TRUE)
    cat("rank1 should be an integer between 1 and " , min(i,l,j*k) , fill=TRUE )
    cat(" ",fill=TRUE)
    checkinput=0
  }
  
  if( rank2 > min(j,i*k) )
  {
    cat(" ",fill=TRUE)
    cat("rank2 should be an integer between 1 and " , min(j,i*k) , fill=TRUE)
    cat(" ",fill=TRUE)
    checkinput=0
  }
  
  if( rank3 > min(k,i*j) )
  {
    cat(" ",fill=TRUE)
    cat("rank3 should be an integer between 1 and " , min(k,i*j) , fill=TRUE)
    cat(" ",fill=TRUE)
    checkinput=0
  }
  
  if ( (rank1 > rank2*rank3) | (rank2 > rank1*rank3) | (rank3 > rank1*rank2) )
  {
    cat(" ",fill=TRUE)
    cat("None of the ranks can be larger than the products of the other two (e.g., rank1 > rank2*rank3 is not allowed",fill=TRUE)
    cat(" ",fill=TRUE)
    checkinput=0
  }
  
  if ( (noise < 0) || (noise >= 1) )
  {
    cat(" ",fill=TRUE)
    cat("Noise is wrong. Should be a value between 0 and 1 (0 is possible, 1 not)",fill=TRUE)
    cat(" ",fill=TRUE)
    checkinput=0
  }
  
  if( length(ExplVarVector) != (rank1+nExtraComponents) )
  {
    cat(" ",fill=TRUE)
    cat("ExplVarVector should contain ", rank1+1 ," elements",fill=TRUE)
    cat(" ",fill=TRUE)
    checkinput=0
  }
  else
  {
    if ( abs( sum(ExplVarVector) - 1 ) > .0000000001 )
    {
      cat(" ",fill=TRUE)
      cat("ExplVarVector should sum to one",fill=TRUE)
      cat(" ",fill=TRUE)
      checkinput=0
    }
    else
    {
      if( sum( sign(ExplVarVector) == -1 ) > 0 )
      {
        cat(" ",fill=TRUE)
        cat("ExplVarVector should only contain positive values (between 0 and 1)",fill=TRUE)
        cat(" ",fill=TRUE)
        checkinput=0
      }
    }
  }
  
  if( ( nExtraComponents < 0 ) | ( (nExtraComponents%%1) != 0 ) )
  {
    cat(" ",fill=TRUE)
    cat("nExtraComponents should be a positive integer (or zero)",fill=TRUE)
    cat(" ",fill=TRUE)
    checkinput=0
  }
  
  if ( (OriginalAlfa < 0) || (OriginalAlfa > 1) )
  {
    cat(" ",fill=TRUE)
    cat("OriginalAlfa should be between 0 and 1",fill=TRUE)
    cat(" ",fill=TRUE)
    checkinput=0
  }
  
  if ( (AlternativeLossF !=0) && (AlternativeLossF != 1) )
  {
    cat(" ",fill=TRUE)
    cat("AlternativeLossF should be 0 or 1",fill=TRUE)
    cat(" ",fill=TRUE)
    checkinput=0
  }
  
  # if ( length(OrthonVector) != 4 )
  # {
  #   cat(" ",fill=TRUE)
  #   cat("OrthonVector should contain 4 elements that are only 0 (no orthogonality) or 1 (orthogonality)",fill=TRUE)
  #   cat(" ",fill=TRUE)
  #   checkinput=0
  # }
  # else
  # {
  #   if ( sum( (OrthonVector !=0) & (OrthonVector != 1) ) != 0 )
  #   {
  #     cat(" ",fill=TRUE)
  #     cat("OrthonVector should contain 4 elements that are only 0 (no orthogonality) or 1 (orthogonality)",fill=TRUE)
  #     cat(" ",fill=TRUE)
  #     checkinput=0
  #   }
  #   else
  #   {
  #     if( ( OrthonVector[1]==1 ) & ( OrthonVector[4]==1 ) )
  #     {
  #       cat(" ",fill=TRUE)
  #       cat("OrthonVector for A and B2 cannot be both 1 at the same time",fill=TRUE)
  #       cat(" ",fill=TRUE)
  #       checkinput=0
  #     }
  #   }
  # }
  
  if (checkinput == 1)
  {
    # if I > JK (or J > IK or K > IJ) then Kiers has a more effcient algortihm for T3
    
    # generate matrices
    Sol = list()
    Sol$nRow = i
    Sol$nCol1 = j
    Sol$nSlice = k
    Sol$nCol2 = l
    Sol$RankVector = cbind( rank1 , rank2 , rank3 )
    Sol$Error = noise
    # Sol$OrthonVector = OrthonVector
    Sol$ExplVarVector = ExplVarVector
    Sol$OriginalAlfa = OriginalAlfa
    Sol$AlternativeLossF = AlternativeLossF
    Sol$nExtraComponents = nExtraComponents
    Sol$GenerateSeed = GenerateSeed
    
    # centered component matrices (A, B1 and C orthonormal; B2 is taken orthogonal)
    Sol$A = orth( scale( matrix( rnorm( i * (rank1+(2*nExtraComponents)) , 0 , 1 ) , i , rank1+(2*nExtraComponents) ) , T , F ) )
    Sol$B1 = orth( scale( matrix( rnorm( j * rank2 , 0 , 1 ) , j , rank2 ) , T , F ) )
    Sol$C = orth( scale( matrix( rnorm( k * rank3 , 0 , 1 ) , k , rank3 ) , T , F ) ) 
    Sol$B2orth = orth( scale( matrix( rnorm( l * (rank1+nExtraComponents) , 0 , 1 ) , l , rank1+nExtraComponents ) , T , F ) )
    Sol$B2 = matrix( 0 , l , rank1+nExtraComponents )
    
    # if ( OrthonVector[1] == 1 ) # orthonormal (centered) component matrices
    # {
    #   if ( dim(Sol$A)[1] >= dim(Sol$A)[2]  )
    #   {
    #      Sol$A = orth( Sol$A )
    #   }
    #   else
    #   {
    #     cat(" ",fill=TRUE)
    #     cat("A has more components than elements (rows)",fill=TRUE)
    #   }
    # }
        
    # if ( OrthonVector[2] == 1 ) # orthonormal (centered) component matrices
    # {
    #   if ( dim(Sol$B1)[1] >= dim(Sol$B1)[2]  )
    #   {
    #     Sol$B1 = orth( Sol$B1 )
    #   }
    #   else
    #   {
    #     cat(" ",fill=TRUE)
    #     cat("B1 has more components than elements (rows)",fill=TRUE)
    #   }
    # }
    
    # if ( OrthonVector[3] == 1 ) # orthonormal (centered) component matrices
    # {
    #   if ( dim(Sol$C)[1] >= dim(Sol$C)[2]  )
    #   {
    #     Sol$C = orth( Sol$C )
    #   }
    #   else
    #   {
    #     cat(" ",fill=TRUE)
    #     cat("C has more components than elements (rows)",fill=TRUE)
    #   }
    # }
    
    # if ( OrthonVector[4] == 1 ) # orthonormal (centered) component matrices
    # {
    #   if ( dim(Sol$B2orth)[1] >= dim(Sol$B2orth)[2]  )
    #   {
    #     Sol$B2orth = orth( Sol$B2orth )
    #   }
    #   else
    #   {
    #     cat(" ",fill=TRUE)
    #     cat("B2 has more components than elements (rows)",fill=TRUE)
    #   }
    # }
    
    if ( nExtraComponents > 0 )
    {
      Sol$A1 = Sol$A[ , cbind( t( 1:rank1 ) , t( (rank1+1):(rank1+nExtraComponents) ) ) ]
      Sol$A2 = Sol$A[ , cbind( t( 1:rank1 ) , t( (rank1+1+nExtraComponents):(rank1+nExtraComponents+nExtraComponents) ) ) ]
      #Sol$A2 = Sol$A1 ## in case you only want common components
    }
    else
    {
      Sol$A1 = Sol$A
      Sol$A2 = Sol$A
    }
    
    # generate core structure (adapted to incorporate additional component)
    Sol$Core = array( rnorm( (rank1+nExtraComponents)*rank2*rank3 , 0 , 1 ) , cbind( rank1+nExtraComponents , rank2 , rank3 ) )
    CoreSq = Sol$Core * Sol$Core
    TotalSqCore = sum(CoreSq)
    
    for (rowtel in 1:dim(Sol$Core)[1])
    {
      tempweight = ( ExplVarVector[rowtel] / sum( CoreSq[rowtel,,] ) ) * TotalSqCore
      Sol$Core[rowtel,,] = Sol$Core[rowtel,,] * sqrt(tempweight)
      rm(tempweight)
    }
    
    # convert three-way core array (q1 x q2 x q3) to a core matrix (q1 x q2q3)
    Sol$CoreMatrix = matrix( Sol$Core , rank1+nExtraComponents , rank2*rank3 )
    
    # Computing the true model
    Sol$TrueWay3DMatrix = matrix( 0 , i , j*k )
    Sol$TrueWay3DMatrixNoExtraComp = matrix( 0 , i , j*k )
    for ( rowtel in 1:(rank1+nExtraComponents) )
    {
      for ( coltel in 1:rank2 )
      {
        for ( slicetel in 1:rank3 )
        {
          temparray = Sol$Core[rowtel,coltel,slicetel] * Sol$A1[,rowtel] %*% kronecker( t( Sol$C[,slicetel] ) , t( Sol$B1[,coltel] ) ) # is i-jk
          Sol$TrueWay3DMatrix = Sol$TrueWay3DMatrix + temparray
          if ( rowtel <= rank1 )
          {
            Sol$TrueWay3DMatrixNoExtraComp = Sol$TrueWay3DMatrixNoExtraComp + temparray
          }
          rm(temparray)
        }
      }
    }
    
    # convert matricized data matrix to a 3D-data array
    Sol$TrueWay3D = array( Sol$TrueWay3DMatrix , cbind(i,j,k) )
    Sol$TrueWay3DnoExtraComp = array( Sol$TrueWay3DMatrixNoExtraComp , cbind(i,j,k) )
    
    EigVals2D = runif( rank1+nExtraComponents , 0 , 1 )
    TotalArraySq = sum( EigVals2D * EigVals2D )
    Sol$TrueWay2D = matrix( 0 , i , l )
    Sol$TrueWay2DnoExtraComp = matrix( 0 , i , l )
    for( rowtel in 1:(rank1+nExtraComponents) )
    {
      tempweight = ExplVarVector[rowtel] * TotalArraySq
      tempmatrix = sqrt(tempweight) * Sol$A2[,rowtel] %*% t( Sol$B2orth[,rowtel] )
      Sol$TrueWay2D = Sol$TrueWay2D + tempmatrix
      if ( rowtel <= rank1 )
      {
        Sol$TrueWay2DnoExtraComp = Sol$TrueWay2DnoExtraComp + tempmatrix
      }
      rm(tempmatrix)
      rm(tempweight)
    }
    
    EigVal = sqrt( ExplVarVector * TotalArraySq )
    tempIdentity = diag( rank1+nExtraComponents )
    diag(tempIdentity) = EigVal
    Sol$B2 = Sol$B2orth %*% tempIdentity
    Sol$EigenValues2D = EigVal
    rm(EigVal,tempIdentity)
    
    # adding error
    error3D = array( rnorm( i * j * k , 0 , 1 ), cbind(i,j,k) )
    error2D = matrix( rnorm( i * l , 0 , 1 ) , i , l )
    error3D = error3D * sqrt( sum( Sol$TrueWay3D ^ 2 ) / sum( error3D ^ 2 ) )
    error2D = error2D * sqrt( sum( Sol$TrueWay2D ^ 2 ) / sum( error2D ^ 2 ) )
    errorlevel = noise / ( 1 - noise );
    Sol$Way3D = Sol$TrueWay3D + ( error3D * sqrt(errorlevel) )
    Sol$Way3DMatrix = matrix( Sol$Way3D , i , j*k )
    Sol$Way2D = Sol$TrueWay2D + ( error2D * sqrt(errorlevel) )
    
    # computing badness-of-data  
    if ( AlternativeLossF == 1 )
    {
      ssq3D = sum( Sol$Way3D ^ 2 )
      ssq2D = sum( Sol$Way2D ^ 2 )
      Alfa = ( OriginalAlfa * ssq2D ) / ( ( OriginalAlfa * ssq2D ) + ( ( 1 - OriginalAlfa ) * ssq3D ) )
    }
    else
    {    
      Alfa = OriginalAlfa
    }
        
    Sol$BOD3D = sum( ( Sol$Way3D - Sol$TrueWay3DnoExtraComp ) ^ 2 )
    Sol$BOD2D = sum( ( Sol$Way2D - Sol$TrueWay2DnoExtraComp ) ^ 2 )
    Sol$BOD3Dsize = Sol$BOD3D / ( i * j * k )
    Sol$BOD2Dsize = Sol$BOD2D / ( i * l )  
    Sol$BOD = ( Alfa * Sol$BOD3D ) + ( ( 1 - Alfa ) * Sol$BOD2D )
    Sol$BODsize = Sol$BOD / ( ( i * j * k ) + ( i * l ) )
    
    Sol$BOD3DExtra = sum( ( Sol$Way3D - Sol$TrueWay3D ) ^ 2 )
    Sol$BOD2DExtra = sum( ( Sol$Way2D - Sol$TrueWay2D ) ^ 2 )
    Sol$BOD3DsizeExtra = Sol$BOD3DExtra / ( i * j * k )
    Sol$BOD2DsizeExtra = Sol$BOD2DExtra / ( i * l )  
    Sol$BODExtra = ( Alfa * Sol$BOD3DExtra ) + ( ( 1 - Alfa ) * Sol$BOD2DExtra )
    Sol$BODsizeExtra = Sol$BODExtra / ( ( i * j * k ) + ( i * l ) )
    return(Sol)
  }
}