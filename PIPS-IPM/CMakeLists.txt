include_directories(Core/Abstract Core/Vector Core/Utilities Core/QpSolvers Core/QpGen
  Core/SparseLinearAlgebra Core/DenseLinearAlgebra Core/Readers 
  Core/LinearSolvers/Ma27Solver Core/LinearSolvers/Ma57Solver)
include_directories(Core/StochLinearAlgebra Core/QpStoch)
add_subdirectory(Core)


set(EXECUTABLE_OUTPUT_PATH "${PROJECT_BINARY_DIR}")

if(HAVE_MA27)
  add_executable(qpgen-sparse-mehrotra.exe Core/QpGen/QpGenSparseMehrotraDriver.C)
  target_link_libraries(qpgen-sparse-mehrotra.exe 
    ooqpgensparse ooqpsparse ooqpbase ooqpmehrotra 
    ${MA27_LIBRARY} ${MATH_LIBS})
endif(HAVE_MA27)

if(HAVE_MA57)
  add_executable(qpgen-sparse-ma57-gondzio.exe Core/QpGen/QpGenSparseGondzioDriver.C)
  target_link_libraries(qpgen-sparse-ma57-gondzio.exe 
    ooqpgensparse ooqpsparse ooqpbase ooqpgondzio 
    ${MA57_LIBRARY} ${METIS_LIBRARY} ${MATH_LIBS})
endif(HAVE_MA57)

if(HAVE_MA57 AND HAVE_METIS)
  add_executable(pipsipm-unitcommit.exe Drivers/unitCommitment.cpp)
  target_link_libraries(pipsipm-unitcommit.exe
    ooqpstoch ooqpstochla ooqpmehrotrastoch
    ooqpgensparse ooqpbase ooqpsparse ooqpdense
    ${MA57_LIBRARY} ${METIS_LIBRARY} ${MATH_LIBS})

  add_executable(pipsipmFromRaw Drivers/pipsipmFromRaw.cpp)
  target_link_libraries(pipsipmFromRaw
    stochInput ${COIN_LIBS}
    ooqpstoch ooqpstochla ooqpmehrotrastoch
    ooqpgensparse ooqpbase ooqpsparse ooqpdense
    ${MA57_LIBRARY} ${METIS_LIBRARY} ${MATH_LIBS})

  add_executable(ooqpFromRaw Drivers/ooqpFromRaw.cpp)
  target_link_libraries(ooqpFromRaw
    stochInput ${COIN_LIBS}
    ooqpgensparse ooqpbase ooqpmehrotra ooqpsparse ooqpdense
    ${MA57_LIBRARY} ${METIS_LIBRARY} ${MATH_LIBS})

endif(HAVE_MA57 AND HAVE_METIS)
