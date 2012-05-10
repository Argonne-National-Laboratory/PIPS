include_directories(Core/Abstract Core/Vector Core/Utilities Core/QpSolvers Core/QpGen
  Core/SparseLinearAlgebra Core/DenseLinearAlgebra Core/Readers 
  Core/LinearSolvers/Ma27Solver Core/LinearSolvers/Ma57Solver)
add_subdirectory(Core)


set(EXECUTABLE_OUTPUT_PATH "${PROJECT_BINARY_DIR}")

add_executable(qpgen-sparse-mehrotra.exe Core/QpGen/QpGenSparseMehrotraDriver.C)
target_link_libraries(qpgen-sparse-mehrotra.exe ooqpgensparse ooqpsparse ooqpbase ooqpmehrotra ${MA27_LIBRARY} ${MATH_LIBS})
