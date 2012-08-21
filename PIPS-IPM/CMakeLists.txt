#add_definitions(-DTIMING -DSTOCH_TESTING) # timing output
include_directories(Core/Abstract Core/Vector Core/Utilities Core/QpSolvers Core/QpGen
  Core/SparseLinearAlgebra Core/DenseLinearAlgebra Core/Readers 
  Core/LinearSolvers/Ma27Solver Core/LinearSolvers/Ma57Solver 
  Core/LinearSolvers/Ma86Solver Core/LinearSolvers/PardisoSolver)
include_directories(Core/StochLinearAlgebra Core/QpStoch)
add_subdirectory(Core)


set(EXECUTABLE_OUTPUT_PATH "${PROJECT_BINARY_DIR}")

if(HAVE_MA27)
  add_executable(qpgen-sparse-mehrotra.exe Core/QpGen/QpGenSparseMehrotraDriver.C)
  target_link_libraries(qpgen-sparse-mehrotra.exe 
    ooqpgensparse ooqpsparse ooqpbase ooqpmehrotra ${PARDISO_LIBRARY}
    ${MA27_LIBRARY} ${MATH_LIBS})
endif(HAVE_MA27)

if(HAVE_MA57)
  add_executable(qpgen-sparse-ma57-gondzio.exe Core/QpGen/QpGenSparseGondzioDriver.C)
  target_link_libraries(qpgen-sparse-ma57-gondzio.exe 
    ooqpgensparse ooqpsparse ooqpbase ooqpgondzio ooqpdense  ${PARDISO_LIBRARY}}
    ${MA57_LIBRARY} ${METIS_LIBRARY} ${MATH_LIBS})
endif(HAVE_MA57)

if(HAVE_MA57 AND HAVE_METIS)
  # add_executable(pipsipm-unitcommit.exe Drivers/unitCommitment.cpp)
  # target_link_libraries(pipsipm-unitcommit.exe
  #   ooqpstoch ooqpstochla ooqpmehrotrastoch
  #   ooqpgensparse ooqpbase ooqpsparse ooqpdense
  #   ${MA57_LIBRARY} ${METIS_LIBRARY} ${MATH_LIBS})

  add_executable(pipsipmFromRaw Drivers/pipsipmFromRaw.cpp)
  target_link_libraries(pipsipmFromRaw
    stochInput ${COIN_LIBS}
    ooqpstoch ooqpstochla ooqpmehrotrastoch
    ooqpgensparse ooqpbase ooqpsparse ooqpdense ${PARDISO_LIBRARY}}
    ${MA57_LIBRARY} ${METIS_LIBRARY} ${MATH_LIBS})

  add_executable(ooqpFromRaw Drivers/ooqpFromRaw.cpp)
  target_link_libraries(ooqpFromRaw
    stochInput ${COIN_LIBS}
    ooqpgensparse ooqpbase ooqpmehrotra ooqpsparse ooqpdense ${PARDISO_LIBRARY}}
    ${MA57_LIBRARY} ${METIS_LIBRARY} ${MATH_LIBS})
  
  add_executable(qpgen-sparse-ma57-mehrotra.exe Core/QpGen/QpGenSparseMa57MehrotraDriver.C)
  target_link_libraries(qpgen-sparse-ma57-mehrotra.exe 
    ooqpgensparse ooqpsparse ooqpbase ooqpmehrotra ooqpdense 
    ${MA57_LIBRARY} ${METIS_LIBRARY} ${MATH_LIBS})
  
endif(HAVE_MA57 AND HAVE_METIS)

if (HAVE_MA86 AND HAVE_METIS)
  add_executable(qpgen-sparse-ma86-mehrotra.exe Core/QpGen/QpGenSparseMa86MehrotraDriver.C)
  target_link_libraries(qpgen-sparse-ma86-mehrotra.exe
    ooqpgensparse ooqpsparse ooqpbase ooqpmehrotra ooqpdense
    ${MA86_LIB} ${METIS_LIBRARY} ${MATH_LIBS})
endif(HAVE_MA86 AND HAVE_METIS)

if(HAVE_PARDISO)
  add_executable(pipsipmFromRaw_schur Drivers/pipsipmFromRaw_schur.cpp)
  target_link_libraries(pipsipmFromRaw_schur
    stochInput ${COIN_LIBS}
    ooqpstoch ooqpstochla ooqpmehrotrastoch
    ooqpgensparse ooqpbase ooqpsparse ooqpdense ${PARDISO_LIBRARY}}
    ${MA57_LIBRARY} ${METIS_LIBRARY} ${MATH_LIBS})
endif(HAVE_PARDISO)

