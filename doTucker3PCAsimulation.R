doTucker3PCAsimulation <- function( nRow , nCol3D , nSlice , nCol2D , r1 , r2 , r3 , noise , Orthon , ExplVarVector , OriginalAlfa = .5 , AlternativeLossF = 1 , nExtraComponents = 0 , conv = .0000001 , nRuns = 100 , GenerateSeed , AnalSeed )
{
  source("GenerateSeed.R")
  source("GenerateTucker3_PCA.R")
  source("T3PCAfunc.R")
  source("T3PCAsegmented.R")
  source("T3PCAseparate.R")
  source("ComputeRecoveryT3PCA.R")
  
  TrueSol = GenerateTucker3_PCA( nRow , nCol3D , nSlice , nCol2D , r1 , r2 , r3 , noise , Orthon , ExplVarVector , OriginalAlfa , AlternativeLossF , nExtraComponents , GenerateSeed )
  
  AnalysisSeeds = GenerateSeed( 4 , AnalSeed )  
  SolSimul = T3PCAfunc( TrueSol$Way3D , TrueSol$Way2D , nRow , nCol3D , nSlice , nCol2D , r1 , r2 , r3 , conv , OriginalAlfa , AlternativeLossF , nRuns , AnalysisSeeds[1] )
  SolSegmented3D = T3PCAsegmented( TrueSol$Way3D , TrueSol$Way2D , r1 , r2 , r3 , conv , OriginalAlfa , AlternativeLossF , nRuns , 1 , AnalysisSeeds[2] )
  SolSegmented2D = T3PCAsegmented( TrueSol$Way3D , TrueSol$Way2D , r1 , r2 , r3 , conv , OriginalAlfa , AlternativeLossF , nRuns , 2 , AnalysisSeeds[3] )
  SolSeparate = T3PCAseparate( TrueSol$Way3D , TrueSol$Way2D , r1 , r2 , r3 , conv , OriginalAlfa , AlternativeLossF , nRuns , AnalysisSeeds[4] )
  
  RecovSimul = ComputeRecoveryT3PCA( TrueSol , SolSimul , 0 )
  RecovSegmented3D = ComputeRecoveryT3PCA( TrueSol , SolSegmented3D , 0 )
  RecovSegmented2D = ComputeRecoveryT3PCA( TrueSol , SolSegmented2D , 0 )
  RecovSeparate = ComputeRecoveryT3PCA( TrueSol , SolSeparate , 1 )
  
  labs = cbind( "" , "A1" , "B1" , "C" , "Core" , "A2" , "B2" )
  rec1 = cbind( "simulult" , RecovSimul$TuckerA , RecovSimul$TuckerB1 , RecovSimul$TuckerC , RecovSimul$TuckerCoreSize , RecovSimul$TuckerA , RecovSimul$TuckerB2optimal )
  rec2 = cbind( "segmen3D" , RecovSegmented3D$TuckerA , RecovSegmented3D$TuckerB1 , RecovSegmented3D$TuckerC , RecovSegmented3D$TuckerCoreSize , RecovSegmented3D$TuckerA , RecovSegmented3D$TuckerB2optimal )
  rec3 = cbind( "segmen2D" , RecovSegmented2D$TuckerA , RecovSegmented2D$TuckerB1 , RecovSegmented2D$TuckerC , RecovSegmented2D$TuckerCoreSize , RecovSegmented2D$TuckerA , RecovSegmented2D$TuckerB2optimal )
  rec4 = cbind( "separate" , RecovSeparate$TuckerA1 , RecovSeparate$TuckerB1 , RecovSeparate$TuckerC , RecovSeparate$TuckerCoreSize , RecovSeparate$TuckerA2optimal , RecovSeparate$TuckerB2optimal )
  RecovSummary = t( cbind( t(labs) , t(rec1) , t(rec2) , t(rec3) , t(rec4) ) )
  rm( labs , rec1 , rec2 , rec3 , rec4 )
  
  Out = list()
  Out$TrueSol = TrueSol
  Out$SolSimul = SolSimul
  Out$SolSegmented3D = SolSegmented3D
  Out$SolSegmented2D = SolSegmented2D
  Out$SolSeparate = SolSeparate
  Out$RecovSimul = RecovSimul
  Out$RecovSegmented3D = RecovSegmented3D
  Out$RecovSegmented2D = RecovSegmented2D
  Out$RecovSeparate = RecovSeparate
  Out$RecovSummary = RecovSummary
  return( Out )
}