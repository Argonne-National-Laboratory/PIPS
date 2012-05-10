add_subdirectory(QpSolvers)
add_subdirectory(SparseLinearAlgebra)
add_subdirectory(DenseLinearAlgebra)

add_library(ooqpbase Abstract/OoqpVersion.C Abstract/Variables.C Abstract/Data.C Abstract/Solver.C Abstract/Status.C 
  Abstract/OoqpMonitor.C Abstract/IotrRefCount.C Abstract/DoubleLinearSolver.C 
  Vector/OoqpVector.C Vector/SimpleVector.C Vector/VectorUtilities.C
  Utilities/drand.C Utilities/sort.C)

add_library(ooqpgensparse 
  QpGen/QpGenVars.C QpGen/QpGenData.C QpGen/QpGenResiduals.C QpGen/QpGen.C QpGen/QpGenLinsys.C #QpGen
  QpGen/QpGenSparseSeq.C QpGen/QpGenSparseLinsys.C #QpGenSparse
  Readers/MpsReader.C Readers/hash.C #Readers
  LinearSolvers/Ma27Solver/Ma27Solver.C QpGen/QpGenSparseMa27.C) #Linear solvers
