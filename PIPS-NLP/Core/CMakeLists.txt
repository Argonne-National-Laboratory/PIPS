add_subdirectory(SparseLinearAlgebra)
add_subdirectory(DenseLinearAlgebra)
add_subdirectory(StochLinearAlgebra)

add_subdirectory(NlpSolvers)
add_subdirectory(NlpInfo)
add_subdirectory(CInterface)
add_subdirectory(NlpStoch)
#add_subdirectory(MtxSchurDecomp)
add_subdirectory(StructureSolver)


if(HAVE_PETSC)
add_subdirectory(PetscLinearAlgebra)
endif(HAVE_PETSC)


add_library(nlpbase Abstract/OoqpVersion.C Abstract/Variables.C Abstract/Data.C Abstract/Solver.C Abstract/Status.C 
  Abstract/OoqpMonitor.C Abstract/IotrRefCount.C Abstract/DoubleLinearSolver.C 
  Abstract/SolverOption.C Abstract/RegularizationAlg.C Abstract/pipsOptions.C
  Vector/OoqpVector.C Vector/SimpleVector.C Vector/VectorUtilities.C
  Utilities/drand.C Utilities/sort.C)

if(HAVE_MA27)
	set(solvers LinearSolvers/Ma27Solver/Ma27Solver.C )
endif(HAVE_MA27)

if(HAVE_MA57)
	set(solvers ${solvers} LinearSolvers/Ma57Solver/Ma57Solver.C )
endif(HAVE_MA57)

if(HAVE_MA86)
	set(solvers ${solvers} LinearSolvers/Ma86Solver/Ma86Solver.C )
endif(HAVE_MA86)

if(HAVE_PARDISO)
	set(solvers ${solvers} LinearSolvers/PardisoSolver/PardisoSolver.C 
	LinearSolvers/PardisoSolver/PardisoSchurSolver.C )
endif(HAVE_PARDISO)


if(HAVE_UMFPACK)
	add_library(LU_solvers LinearSolvers/UmfPackSolver/UmfPackSolver.C)
else(HAVE_UMFPACK)
	add_library(LU_solvers Utilities/empty.cpp)
endif(HAVE_UMFPACK)

#if(HAVE_PARDISO OR HAVE_MA57)
add_library(nlpwithsolver NlpGen/NlpGenSparseWithSolver.C)
#endif(HAVE_PARDISO OR HAVE_MA57)


add_library(nlpgensparse
  NlpGen/NlpGenVars.C NlpGen/NlpGenData.C NlpGen/NlpGenResiduals.C 
  NlpGen/NlpGen.C NlpGen/NlpGenLinsys.C NlpInfo/NlpInfo.C #NlpGen
  NlpGen/NlpGenSparse.C NlpGen/NlpGenSparseLinsys.C #NlpGenSparseLinSys
  Readers/MpsReader.C Readers/hash.C ${solvers}
  )


