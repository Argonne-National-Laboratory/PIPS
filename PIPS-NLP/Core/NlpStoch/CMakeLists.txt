add_library(nlpstoch
  sFactory.C    
  sFactoryAug.C 
  sFactoryAugSparseRowMajSC.C
  sData.C
  sLinsys.C   
  sLinsysRoot.C 
  sLinsysRootAug.C
  sLinsysRootAugSparseRowMaj.C
  sLinsysLeaf.C  
  sVars.C 
  StochMonitor.C 
  sResiduals.C
  sTree.C 
  sTreeImpl.C 
  sInfo.C
)

add_library(aggregationstoch sFactoryAugAggregationPrecond.C sLinsysRootAggregation.C ) 
