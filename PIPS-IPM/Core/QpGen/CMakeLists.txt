
# libooqpgendense_OBJ = $(QPGENDENSEOBJ) $(MPSREADEROBJ) \
# 	$(CQPGENDENSEOBJ) \
# 	$(QPGENDIR)/QpGenDense.o 


# QPGENCOMMONSRC = $(QPGENDIR)/QpGenVars.C $(QPGENDIR)/QpGenData.C \
# 	$(QPGENDIR)/QpGenResiduals.C $(QPGENDIR)/QpGen.C \
# 	$(QPGENDIR)/QpGenLinsys.C 

# QPGENSPARSESRC = $(QPGENCOMMONSRC) $(QPGENDIR)/QpGenSparseSeq.C \
# 	$(QPGENDIR)/QpGenSparseLinsys.C
# QPGENSPARSEOBJ = $(QPGENSPARSESRC:.C=.o)

# QPGENDENSESRC = $(QPGENCOMMONSRC) \
# 	$(QPGENDIR)/QpGenDenseLinsys.C $(QPGENDIR)/QpGenDense.C
# QPGENDENSEOBJ = $(QPGENDENSESRC:.C=.o)

# QPGENSTOCHCOMMONSRC = $(QPGENCOMMONSRC) \
# 	      $(QPGENDIR)/QpGenStoch.C \
# 	      $(QPGENDIR)/sFactory.C \
# 	      $(QPGENDIR)/QpGenStochData.C \
# 	      $(QPGENDIR)/QpGenStochLinsys.C \
# 	      $(QPGENDIR)/sLinsys.C \
# 	      $(QPGENDIR)/QpGenStochLinsysRoot.C \
# 	      $(QPGENDIR)/sLinsysRoot.C \
# 	      $(QPGENDIR)/QpGenStochLinsysLeaf.C \
# 	      $(QPGENDIR)/sLinsysLeaf.C \
# 	      $(QPGENDIR)/QpGenStochVars.C \
# 	      $(QPGENDIR)/StochMonitor.C \
# 	      $(QPGENDIR)/QpGenResiduals2.C 

# QPGENSTOCHCOMMONOBJ = $(QPGENSTOCHCOMMONSRC:.C=.o)

# QPGENSTOCHAUGEXTSRC = $(QPGENDIR)/QpGenStochAugExt.C \
# 	      $(QPGENDIR)/QpGenStochLinsysRootAug.C \
# 	      $(QPGENDIR)/QpGenStochLinsysRootAugExt.C
# QPGENSTOCHAUGREDSRC = $(QPGENDIR)/QpGenStochAugRed.C \
# 	      $(QPGENDIR)/sFactoryAug.C \
# 	      $(QPGENDIR)/QpGenStochLinsysRootAug.C \
# 	      $(QPGENDIR)/QpGenStochLinsysRootAugRed.C \
# 	      $(QPGENDIR)/sLinsysRootAug.C
# QPGENSTOCHAUGREDPRSRC = $(QPGENDIR)/QpGenStochAugRedPr.C \
#               $(QPGENDIR)/QpGenStochLinsysRootAugRedPrecond.C \
# 	      $(QPGENDIR)/QpStochAugRedPrCG.C \
# 	      $(QPGENDIR)/QpStochLinsysRootAugRedPrCG.C \
# 	      $(QPGENDIR)/QpStochAugRedPrPCG.C \
# 	      $(QPGENDIR)/QpStochLinsysRootAugRedPrPCG.C \
# 	      $(QPGENDIR)/RemoteMatTimesVec.C
# QPGENSTOCHNRMEQNSRC = $(QPGENDIR)/QpGenStochNrmEqn.C \
# 	      $(QPGENDIR)/QpGenStochLinsysRootAug.C \
# 	      $(QPGENDIR)/QpGenStochLinsysRootAugRed.C \
# 	      $(QPGENDIR)/QpGenStochLinsysRootNrmEqn.C
# QPGENSTOCHNRMEQNPRSRC = $(QPGENDIR)/QpStochNrmEqnPr.C \
# 	      $(QPGENDIR)/QpGenStochLinsysRootAug.C \
# 	      $(QPGENDIR)/QpGenStochLinsysRootAugRed.C \
# 	      $(QPGENDIR)/QpGenStochLinsysRootNrmEqn.C \
# 	      $(QPGENDIR)/QpGenStochLinsysRootNrmEqnPrecond.C \
# 	      $(QPGENDIR)/RemoteMatTimesVec.C
# QPGENSTOCHSCASRC = $(QPGENDIR)/sFactoryAugSca.C \
# 	      $(QPGENDIR)/sLinsysRootAugSca.C
# QPGENSTOCHEMTLSRC = $(QPGENDIR)/sFactoryAugEmtl.C \
# 	      $(QPGENDIR)/sLinsysRootAugEmtl.C \
#         $(QPGENDIR)/sFactoryAugEmtlSym.C \
#         $(QPGENDIR)/sLinsysRootAugEmtlSym.C

# QPGENSTOCHAUGEXTOBJ = $(QPGENSTOCHAUGEXTSRC:.C=.o)
# QPGENSTOCHAUGREDOBJ = $(QPGENSTOCHAUGREDSRC:.C=.o)
# QPGENSTOCHAUGREDPROBJ = $(QPGENSTOCHAUGREDPRSRC:.C=.o)
# QPGENSTOCHNRMEQNOBJ = $(QPGENSTOCHNRMEQNSRC:.C=.o)
# QPGENSTOCHNRMEQNPROBJ = $(QPGENSTOCHNRMEQNPRSRC:.C=.o)
# QPGENSTOCHSCAOBJ = $(QPGENSTOCHSCASRC:.C=.o)
# QPGENSTOCHEMTLOBJ = $(QPGENSTOCHEMTLSRC:.C=.o)


# QPGENSTOCHOBJ = $(QPGENSTOCHSRC:.C=.o) $(QPGENDIR)/QpGenStochDriver.o 

# QPGENSCASRC = $(QPGENCOMMONSRC) \
# 	$(QPGENDIR)/QpGenScaLinsys.C $(QPGENDIR)/QpGenSca.C \
# 	$(QPGENDIR)/QpGenScaMehrotraDriver.C

# QPGENSCAOBJ = $(QPGENSCASRC:.C=.o)


# libooqpgensca_OBJ = $(QPGENSCAOBJ) $(MPSREADEROBJ) 
# libooqpgensca_STATIC = $(locallibdir)/libooqpgensca.a

# libstoch_OBJ = $(QPGENSTOCHCOMMONOBJ)
# libstoch_STATIC = $(locallibdir)/libstoch.a

# libstochaugext_OBJ = $(QPGENSTOCHAUGEXTOBJ) 
# libstochaugext_STATIC = $(locallibdir)/libstochaugext.a
# libstochaugred_OBJ = $(QPGENSTOCHAUGREDOBJ)
# libstochaugred_STATIC = $(locallibdir)/libstochaugred.a
# libstochaugredpr_OBJ = $(QPGENSTOCHAUGREDPROBJ) $(QPGENSTOCHAUGREDOBJ) $(BICGSTABSOLVEROBJ) $(CGSOLVEROBJ)  $(PCGSOLVEROBJ)
# libstochaugredpr_STATIC = $(locallibdir)/libstochaugredpr.a
# libstochnrmeqn_OBJ = $(QPGENSTOCHNRMEQNOBJ)
# libstochnrmeqn_STATIC = $(locallibdir)/libstochnrmeqn.a
# libstochnrmeqnpr_OBJ = $(QPGENSTOCHNRMEQNPROBJ)  $(CGSOLVEROBJ) 
# libstochnrmeqnpr_STATIC = $(locallibdir)/libstochnrmeqnpr.a 

# libstochscaaug_OBJ = $(QPGENSTOCHSCAOBJ)
# libstochscaaug_STATIC = $(locallibdir)/libstochscaaug.a

# libstochemtlaug_OBJ = $(QPGENSTOCHEMTLOBJ)
# libstochemtlaug_STATIC = $(locallibdir)/libstochemtlaug.a

#- libooqpgensparse_OBJ = $(QPGENSPARSEOBJ) $(MPSREADEROBJ) \
#- 	$(CQPGENSPARSEOBJ) \
#- 	$(QPGENDIR)/QpGenSparseMa27.o $(MA27SOLVEROBJ) 
#- libooqpgensparse_STATIC = $(locallibdir)/libooqpgensparse.a

# libooqpgendense_OBJ = $(QPGENDENSEOBJ) $(MPSREADEROBJ) \
# 	$(CQPGENDENSEOBJ) \
# 	$(QPGENDIR)/QpGenDense.o 
# libooqpgendense_STATIC = $(locallibdir)/libooqpgendense.a

# libooqpgenstoch_OBJ = $(QPGENSTOCHOBJ) 
# libooqpgenstoch_STATIC = $(locallibdir)/libooqpgenstoch.a


# qpgen_sparse_gondzio_OBJECTS = $(QPGENDIR)/QpGenSparseGondzioDriver.o \
# 	$(libooqpgondzio_STATIC) $(libooqpgensparse_STATIC) \
# 	$(libooqpsparse_STATIC) $(libooqpbase_STATIC)

# qpgen_sparse_mehrotra_OBJECTS = $(QPGENDIR)/QpGenSparseMehrotraDriver.o \
# 	$(libooqpmehrotra_STATIC) $(libooqpgensparse_STATIC) \
# 	$(libooqpsparse_STATIC) $(libooqpbase_STATIC)


# # qpgen_stochunitc_OBJECTS = $(QPGENDIR)/QpGenStochDriverUC.o \
# # 	$(QPGENDIR)/QpGenStochDriver.o  \
# # 	$(MA57SOLVEROBJ) $(QPGENDIR)/QpGenSparseMa57.o \
# # 	$(BICGSTABSOLVEROBJ) $(CGSOLVEROBJ) \
# # 	$(libooqpgenstoch_STATIC) \
# # 	$(libooqpmehrotrastoch_STATIC) $(libooqpgondzio_STATIC) \
# # 	$(libooqpsparse_STATIC) $(libooqpdense_STATIC) \
# # 	$(libooqpgensparse_STATIC)  \
# # 	$(libooqpbase_STATIC) $(libooqpstochla_STATIC) 



# qpgen_dense_gondzio_OBJECTS = $(QPGENDIR)/QpGenDenseGondzioDriver.o \
# 	$(libooqpgondzio_STATIC) $(libooqpgendense_STATIC) \
# 	$(libooqpdense_STATIC) $(libooqpbase_STATIC)

# qpgen_sparse_ma57_gondzio_OBJECTS = \
# 	$(QPGENDIR)/QpGenSparseMa57GondzioDriver.o \
# 	$(QPGENDIR)/QpGenSparseMa57.o $(MA57SOLVEROBJ)  \
#   $(WSMPSOLVEROBJ) \
# 	$(libooqpgondzio_STATIC) $(libooqpgensparse_STATIC) \
# 	$(libooqpsparse_STATIC) \
#   $(libooqpdense_STATIC) $(libooqpbase_STATIC)

# qpstoch_driver_OBJECTS = $(MA57SOLVEROBJ) $(libstochscaaug_STATIC) \
#   $(libstochemtlaug_STATIC) $(libstoch_STATIC) $(libooqpstochla_STATIC) \
# 		$(libooqpmehrotrastoch_STATIC) $(libstochscaaug_STATIC) \
#     $(libooqpsca_STATIC) $(libooqpemtl_STATIC) \
# 		$(libooqpsparse_STATIC) $(libooqpdense_STATIC) $(libooqpbase_STATIC)

# qpstoch_dualdriver_OBJECTS = $(libooqpmehrotra_STATIC) $(WSMPSOLVEROBJ) $(QPGENDIR)/QpGenSparseMa57.o $(libooqpgensparse_STATIC) \
# 			    $(qpstoch_driver_OBJECTS)

# qpgen_sca_mehrotra_OBJECTS = \
# 	$(QPGENDIR)/QpGenScaMehrotraDriver.o \
# 	$(libooqpgensca_STATIC) $(libooqpmehrotra_STATIC) \
# 	$(libooqpgen_STATIC) $(libooqpsca_STATIC) \
# 	$(libooqpsparse_STATIC) $(libooqpbase_STATIC)
# qpgen-sparse-gondzio.exe: $(qpgen_sparse_gondzio_OBJECTS)
# 	$(CXX) -o $@ $(CXXFLAGS) $(LDFLAGS) \
# 		$(qpgen_sparse_gondzio_OBJECTS) $(MA27LIB) $(BLAS) $(FLIBS)

# qpgen-sparse-mehrotra.exe: $(qpgen_sparse_mehrotra_OBJECTS)
# 	$(CXX) -o $@ $(CXXFLAGS) $(LDFLAGS) \
# 		$(qpgen_sparse_mehrotra_OBJECTS) $(MA27LIB) $(BLAS) $(FLIBS)

# qpgen-dense-gondzio.exe: $(qpgen_dense_gondzio_OBJECTS)
# 	$(CXX) -o $@ $(CXXFLAGS) $(LDFLAGS) \
# 		$(qpgen_dense_gondzio_OBJECTS) $(LAPACK) $(BLAS)  $(FLIBS)

# qpgen-comparisontest.exe: $(qpgen_comparisontest_OBJECTS)
# 	$(CXX) -o $@ $(CXXFLAGS) $(LDFLAGS) \
#  	        $(qpgen_comparisontest_OBJECTS) \
# 		$(MA57LIB) $(METISLIB) $(LAPACK) $(BLAS)

# qpstoch-unitcommit.exe: $(QPGENDIR)/QpGenStochDriverUC.o $(qpstoch_dualdriver_OBJECTS) \
# 			$(libstochaugext_STATIC) \
# 			$(libstochaugredpr_STATIC) $(libstochaugred_STATIC) \
# 			$(libstochnrmeqnpr_STATIC) $(libstochnrmeqn_STATIC) 
# 	$(CXX) -o $@ $(CXXFLAGS) $(LDFLAGS) $(QPGENDIR)/QpGenStochDriverUC.o \
# 	       $(libstochaugext_STATIC) $(libstochaugred_STATIC) $(libstochaugredpr_STATIC) \
# 	       $(libstochnrmeqnpr_STATIC) $(libstochnrmeqn_STATIC) \
# 	       $(qpstoch_dualdriver_OBJECTS)  $(STOCH_LIBS) $(EMTL_LIB) $(SCA_LIBS)

# qpstoch-unitcommit-mkl.exe: $(QPGENDIR)/QpGenStochDriverUC.o $(qpstoch_dualdriver_OBJECTS) \
# 			$(libstochaugext_STATIC) \
# 			$(libstochaugredpr_STATIC) $(libstochaugred_STATIC) \
# 			$(libstochnrmeqnpr_STATIC) $(libstochnrmeqn_STATIC) 
# 	$(CXX) -o $@ $(CXXFLAGS) $(LDFLAGS) $(QPGENDIR)/QpGenStochDriverUC.o \
# 	       $(libstochaugext_STATIC) $(libstochaugred_STATIC) $(libstochaugredpr_STATIC) \
# 	       $(libstochnrmeqnpr_STATIC) $(libstochnrmeqn_STATIC) \
# 	       $(qpstoch_dualdriver_OBJECTS)  $(STOCH_LIBS)\
#          $(EMTL_LIB) $(SCA_LIBS_MKL)

# qpstoch-unitcommit_nrmeqn.exe: $(QPGENDIR)/QpGenStochDriverUC_nrmeqn.o $(qpstoch_dualdriver_OBJECTS) $(libstochnrmeqn_STATIC)
# 	$(CXX) -o $@ $(CXXFLAGS) $(LDFLAGS) $(QPGENDIR)/QpGenStochDriverUC_nrmeqn.o \
# 	       $(libstochnrmeqn_STATIC) $(qpstoch_dualdriver_OBJECTS)  $(STOCH_LIBS) $(FLIBS)


# qpstoch-unitcommit_large.exe: $(QPGENDIR)/QpGenStochDriverUC_large.o $(qpstoch_dualdriver_OBJECTS) \
# 			      $(libstochaugext_STATIC) \
# 			      $(libstochaugred_STATIC)  $(libstochaugredpr_STATIC) \
# 			      $(libstochnrmeqnpr_STATIC) $(libstochnrmeqn_STATIC) 
# 	$(CXX) -o $@ $(CXXFLAGS) $(LDFLAGS) $(QPGENDIR)/QpGenStochDriverUC_large.o \
# 	       $(libstochaugext_STATIC) $(libstochaugred_STATIC) $(libstochaugredpr_STATIC) \
# 	       $(libstochnrmeqnpr_STATIC) $(libstochnrmeqn_STATIC) \
# 	       $(qpstoch_dualdriver_OBJECTS)  $(STOCH_LIBS) $(LAPACK)  $(BLAS) -lgfortran
# #         $(EMTL_LIB) $(SCA_LIBS)

# qpstoch-unitcommitQuadr_large.exe: $(QPGENDIR)/sDriverUCQuadr_large.o $(qpstoch_dualdriver_OBJECTS) \
# 			      $(libstochaugext_STATIC) \
# 			      $(libstochaugred_STATIC)  $(libstochaugredpr_STATIC) \
# 			      $(libstochnrmeqnpr_STATIC) $(libstochnrmeqn_STATIC) 
# 	$(CXX) -o $@ $(CXXFLAGS) $(LDFLAGS) $(QPGENDIR)/sDriverUCQuadr_large.o \
# 	       $(libstochaugext_STATIC) $(libstochaugred_STATIC) $(libstochaugredpr_STATIC) \
# 	       $(libstochnrmeqnpr_STATIC) $(libstochnrmeqn_STATIC) \
# 	       $(qpstoch_dualdriver_OBJECTS)  $(STOCH_LIBS) $(LAPACK)  $(BLAS) -lgfortran

# qpstoch-unitcommit_hybr.exe: $(QPGENDIR)/QpGenStochDriverUC_hybr.o $(qpstoch_dualdriver_OBJECTS) \
# 			      $(libstochaugext_STATIC) \
# 			      $(libstochaugred_STATIC)  $(libstochaugredpr_STATIC) \
# 			      $(libstochnrmeqnpr_STATIC) $(libstochnrmeqn_STATIC) 
# 	$(CXX) -o $@ $(CXXFLAGS) $(LDFLAGS) $(QPGENDIR)/QpGenStochDriverUC_hybr.o \
# 	       $(libstochaugext_STATIC) $(libstochaugred_STATIC) $(libstochaugredpr_STATIC) \
# 	       $(libstochnrmeqnpr_STATIC) $(libstochnrmeqn_STATIC) \
# 	       $(qpstoch_dualdriver_OBJECTS)  $(STOCH_LIBS) \
#          $(EMTL_LIB) $(SCA_LIBS)

# qpstoch-unitcommit_hybr-mkl.exe: $(QPGENDIR)/QpGenStochDriverUC_hybr.o $(qpstoch_dualdriver_OBJECTS) \
# 			      $(libstochaugext_STATIC) \
# 			      $(libstochaugred_STATIC)  $(libstochaugredpr_STATIC) \
# 			      $(libstochnrmeqnpr_STATIC) $(libstochnrmeqn_STATIC) 
# 	$(CXX) -o $@ $(CXXFLAGS) $(LDFLAGS) $(QPGENDIR)/QpGenStochDriverUC_hybr.o \
# 	       $(libstochaugext_STATIC) $(libstochaugred_STATIC) $(libstochaugredpr_STATIC) \
# 	       $(libstochnrmeqnpr_STATIC) $(libstochnrmeqn_STATIC) \
# 	       $(qpstoch_dualdriver_OBJECTS)  $(STOCH_LIBS) \
#          $(EMTL_LIB) $(SCA_LIBS_MKL)


# qpstoch-building.exe: $(QPGENDIR)/QpGenStochDriverBuilding.o $(qpstoch_dualdriver_OBJECTS) $(libstochaugext_STATIC)
# 	$(CXX) -o $@ $(CXXFLAGS) $(LDFLAGS) $(QPGENDIR)/QpGenStochDriverBuilding.o \
# 	       $(libstochaugext_STATIC) $(qpstoch_dualdriver_OBJECTS)  $(STOCH_LIBS) $(FLIBS)

# qpgen-sca-mehrotra.exe: $(qpgen_sca_mehrotra_OBJECTS)
# 	$(CXX) -o $@ $(CXXFLAGS) $(LDFLAGS) $(LIBS) \
# 		$(qpgen_sca_mehrotra_OBJECTS) $(SCA_LIBS)



# all_ma57: qpgen-sparse-ma57-gondzio.exe
# qpgen-sparse-ma57-gondzio.exe: $(qpgen_sparse_ma57_gondzio_OBJECTS)
# 	$(CXX) -o $@ $(CXXFLAGS) $(LDFLAGS) \
# 		$(qpgen_sparse_ma57_gondzio_OBJECTS) $(MA57LIB) $(METISLIB) $(WSMPLIB) $(BLAS) $(FLIBS)

# $(libooqpgensparse_STATIC) : $(libooqpgensparse_OBJ)
# 	ar ru $@ $(libooqpgensparse_OBJ)
# 	$(RANLIB) $@

# $(libooqpgendense_STATIC) : $(libooqpgendense_OBJ)
# 	ar ru $@ $(libooqpgendense_OBJ)
# 	$(RANLIB) $@

# $(libooqpgenstoch_STATIC) : $(libooqpgenstoch_OBJ)
# 	ar ru $@ $(libooqpgenstoch_OBJ)
# 	$(RANLIB) $@

# $(libstoch_STATIC) :  $(libstoch_OBJ)
# 	ar ru $@ $(libstoch_OBJ)
# 	$(RANLIB) $@

# $(libstochaugext_STATIC) : $(libstochaugext_OBJ)
# 	ar ru $@ $(libstochaugext_OBJ)
# 	$(RANLIB) $@
# $(libstochaugred_STATIC) : $(libstochaugred_OBJ)
# 	ar ru $@ $(libstochaugred_OBJ)
# 	$(RANLIB) $@
# $(libstochaugredpr_STATIC) : $(libstochaugredpr_OBJ)
# 	ar ru $@ $(libstochaugredpr_OBJ)
# 	$(RANLIB) $@
# $(libstochnrmeqn_STATIC) : $(libstochnrmeqn_OBJ)
# 	ar ru $@ $(libstochnrmeqn_OBJ)
# 	$(RANLIB) $@
# $(libstochnrmeqnpr_STATIC) : $(libstochnrmeqnpr_OBJ)
# 	ar ru $@ $(libstochnrmeqnpr_OBJ)
# 	$(RANLIB) $@
# $(libstochscaaug_STATIC) : $(libstochscaaug_OBJ)
# 	ar ru $@ $(libstochscaaug_OBJ)
# 	$(RANLIB) $@
# $(libstochemtlaug_STATIC) : $(libstochemtlaug_OBJ)
# 	ar ru $@ $(libstochemtlaug_OBJ)
# 	$(RANLIB) $@
# $(libooqpgensca_STATIC) : $(libooqpgensca_OBJ)
# 	ar ru $@ $(libooqpgensca_OBJ)
# 	$(RANLIB) $@

# all: qpgen-sparse-mehrotra.exe
# all: qpgen-sparse-gondzio.exe
# all: qpgen-dense-gondzio.exe

# all_stoch: qpstoch-unitcommit.exe qpstoch-unitcommit_large.exe qpstoch-unitcommit_nrmeqn.exe
# all_stoch: qpstoch-building.exe

# all_libs: $(libooqpgensparse_STATIC) $(libooqpgendense_STATIC) 
# all_libs: $(libstoch_STATIC) $(libstochaugext_STATIC)

