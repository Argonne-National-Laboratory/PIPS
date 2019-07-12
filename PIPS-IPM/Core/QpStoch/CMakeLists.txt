
set(LIBRARY_OOQPSTOCH_SOURCES sFactory.C
   sFactoryAug.C
   sFactoryAugSchurLeaf.C
   sFactoryAugComm2SchurLeaf.C
   sData.C
   sLinsys.C
   sLinsysRoot.C
   sLinsysRootAug.C
   sLinsysRootComm2.C
   sLinsysRootAugComm2.C
   sLinsysLeaf.C
   sLinsysLeafSchurSlv.C
   sVars.C
   StochMonitor.C
   sResiduals.C
   sTree.C
   sTreeImpl.C
   sTreeCallbacks.C
   sInterfaceCallbacks.C)

if(HAVE_MUMPS AND WITH_MUMPS_LEAF)
   set(LIBRARY_OOQPSTOCH_SOURCES ${LIBRARY_OOQPSTOCH_SOURCES}
   sFactoryAugMumpsLeaf.C
   sLinsysLeafMumps.C)
endif(HAVE_MUMPS AND WITH_MUMPS_LEAF)

add_library(ooqpstoch ${LIBRARY_OOQPSTOCH_SOURCES})
