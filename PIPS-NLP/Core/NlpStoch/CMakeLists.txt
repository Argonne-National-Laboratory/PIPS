add_library(nlpstoch
  sFactory.C    
  sFactoryAug.C 
  sData.C
  sLinsys.C   
  sLinsysRoot.C 
  sLinsysRootAug.C
  sLinsysLeaf.C  
  sVars.C 
  StochMonitor.C 
  sResiduals.C
  sTree.C 
  sTreeImpl.C 
  sInfo.C
)

add_library(aggregationstoch sFactoryAugAggregationPrecond.C sLinsysRootAggregation.C ) 