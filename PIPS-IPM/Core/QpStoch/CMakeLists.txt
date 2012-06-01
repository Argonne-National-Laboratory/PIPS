add_library(ooqpstoch sFactory.C sFactoryAug.C 
  sData.C
  sLinsys.C sLinsysRoot.C sLinsysRootAug.C 
  sLinsysLeaf.C sVars.C StochMonitor.C sResiduals.C
  sTree.C sTreeImpl.C sTreeCallbacks.C 
  PIPSIPMInterface.C)
