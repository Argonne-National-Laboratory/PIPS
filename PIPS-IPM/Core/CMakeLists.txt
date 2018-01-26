add_subdirectory(QpSolvers)
add_subdirectory(QpPreprocess)
add_subdirectory(SparseLinearAlgebra)
add_subdirectory(DenseLinearAlgebra)
add_subdirectory(StochLinearAlgebra)
add_subdirectory(QpStoch)

add_library(ooqpbase Abstract/OoqpVersion.C Abstract/Variables.C Abstract/Data.C Abstract/Solver.C Abstract/Status.C Abstract/Scaler.C
  Abstract/Presolver.C Abstract/OoqpMonitor.C Abstract/IotrRefCount.C Abstract/DoubleLinearSolver.C 
  Vector/OoqpVector.C Vector/SimpleVector.C Vector/VectorUtilities.C
  Utilities/drand.C Utilities/sort.C)

if(HAVE_MA27)
	set(solvers LinearSolvers/Ma27Solver/Ma27Solver.C QpGen/QpGenSparseMa27.C)
endif(HAVE_MA27)
if(HAVE_MA57)
	set(solvers ${solvers} LinearSolvers/Ma57Solver/Ma57Solver.C QpGen/QpGenSparseMa57.C)
endif(HAVE_MA57)
#if(HAVE_MA86)
#	set(solvers ${solvers} LinearSolvers/Ma86Solver/Ma86Solver.C QpGen/QpGenSparseMa86.C)
#endif(HAVE_MA86)
if(HAVE_PARDISO)
	set(solvers ${solvers} LinearSolvers/PardisoSolver/PardisoSolver.C LinearSolvers/PardisoSolver/PardisoSchurSolver.C)
endif(HAVE_PARDISO)

add_library(ooqpgensparse 
  QpGen/QpGenVars.C QpGen/QpGenData.C QpGen/QpGenResiduals.C QpGen/QpGen.C QpGen/QpGenLinsys.C #QpGen
  QpGen/QpGenSparseSeq.C QpGen/QpGenSparseLinsys.C #QpGenSparse
  Readers/MpsReader.C Readers/hash.C #Readers
  ${solvers})
