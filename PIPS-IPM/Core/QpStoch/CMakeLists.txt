if(HAVE_PARDISO)
  add_library(ooqpstoch sFactory.C sFactoryAug.C 
	sFactoryAugSchurLeaf.C sFactoryAugComm2SchurLeaf.C
	sData.C
	sLinsys.C sLinsysRoot.C sLinsysRootAug.C 
	sLinsysRootComm2.C sLinsysRootAugComm2.C 
	sLinsysLeaf.C 
	sLinsysLeafSchurSlv.C 
	sVars.C StochMonitor.C sResiduals.C
	sTree.C sTreeImpl.C sTreeCallbacks.C 
	sInterfaceCallbacks.C)
else(HAVE_PARDISO)
  add_library(ooqpstoch sFactory.C sFactoryAug.C 
	sData.C
	sLinsys.C sLinsysRoot.C sLinsysRootAug.C 
	sLinsysRootComm2.C sLinsysRootAugComm2.C 
	sLinsysLeaf.C 
	sVars.C StochMonitor.C sResiduals.C
	sTree.C sTreeImpl.C sTreeCallbacks.C 
	sInterfaceCallbacks.C)
endif(HAVE_PARDISO)