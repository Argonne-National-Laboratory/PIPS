add_subdirectory(QpSolvers)
add_subdirectory(QpPreprocess)
add_subdirectory(SparseLinearAlgebra)
add_subdirectory(DenseLinearAlgebra)
add_subdirectory(StochLinearAlgebra)
add_subdirectory(QpStoch)

add_library(ooqpbase Abstract/OoqpVersion.C Abstract/Options.C Abstract/Observer.C Abstract/Variables.C Abstract/Data.C Abstract/Solver.C Abstract/Status.C Abstract/Scaler.C
  Abstract/OoqpMonitor.C Abstract/IotrRefCount.C Abstract/DoubleLinearSolver.C 
  Vector/OoqpVector.C Vector/SimpleVector.C Vector/VectorUtilities.C
  Utilities/drand.C Utilities/sort.C Utilities/pipschecks.C)

if(HAVE_MA27)
	set(solvers LinearSolvers/Ma27Solver/Ma27Solver.C QpGen/QpGenSparseMa27.C)
endif(HAVE_MA27)
if(HAVE_MA57)
	set(solvers ${solvers} LinearSolvers/Ma57Solver/Ma57Solver.C QpGen/QpGenSparseMa57.C)
endif(HAVE_MA57)
#if(HAVE_MA86)
#	set(solvers ${solvers} LinearSolvers/Ma86Solver/Ma86Solver.C QpGen/QpGenSparseMa86.C)
#endif(HAVE_MA86)

if(HAVE_MUMPS)
	set(solvers ${solvers} LinearSolvers/MumpsSolver/MumpsSolverBase.C)
endif(HAVE_MUMPS)

if(WITH_MUMPS_LEAF)
	set(solvers ${solvers} LinearSolvers/MumpsSolver/MumpsSolverLeaf.C)
endif(WITH_MUMPS_LEAF)

if(WITH_MUMPS_ROOT)
	set(solvers ${solvers} LinearSolvers/MumpsSolver/MumpsSolverRoot.C)
endif(WITH_MUMPS_ROOT)

if(HAVE_PARDISO)
	set(solvers ${solvers} LinearSolvers/PardisoSolver/PardisoSolver.C LinearSolvers/PardisoSolver/PardisoSchurSolver.C LinearSolvers/PardisoSolver/PardisoIndefSolver.C)
endif(HAVE_PARDISO)

set(solvers ${solvers} LinearSolvers/Preconditioners/SCsparsifier.C)

add_library(ooqpgensparse 
  QpGen/QpGenVars.C QpGen/QpGenData.C QpGen/QpGenResiduals.C QpGen/QpGen.C QpGen/QpGenLinsys.C QpGen/QpGenOptions.C #QpGen
  QpGen/QpGenSparseSeq.C QpGen/QpGenSparseLinsys.C #QpGenSparse
  Readers/MpsReader.C Readers/hash.C #Readers
  ${solvers})
