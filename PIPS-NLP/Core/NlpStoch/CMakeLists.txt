if(HAVE_MUMPS)
	set(sparseSCFiles sFactoryAugSparseRowMajSC.C sLinsysRootAugSparseRowMaj.C)
	set(sparseSCFiles ${sparseSCFiles})
endif(HAVE_MUMPS)

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
  ${sparseSCFiles}
)

add_library(aggregationstoch sFactoryAugAggregationPrecond.C sLinsysRootAggregation.C ) 
