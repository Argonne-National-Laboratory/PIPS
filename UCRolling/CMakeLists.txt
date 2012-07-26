

include_directories(../PIPS-S/Basic)
include_directories(../PIPS-S/Core)
include_directories(../PIPS-S/CoinBALPFactorization)
include_directories(../Lagrange/LagrangeSubproblemSolver)
include_directories(../Lagrange/RecourseSubproblemSolver)
include_directories(../Lagrange/Drivers)
include_directories(../Lagrange/Bundle)

include_directories(../PIPS-IPM/Core/Abstract ../PIPS-IPM/Core/Vector ../PIPS-IPM/Core/Utilities ../PIPS-IPM/Core/QpSolvers ../PIPS-IPM/Core/QpGen
  ../PIPS-IPM/Core/SparseLinearAlgebra ../PIPS-IPM/Core/DenseLinearAlgebra ../PIPS-IPM/Core/Readers
  ../PIPS-IPM/Core/LinearSolvers/Ma57Solver)
include_directories(../PIPS-IPM/Core/StochLinearAlgebra ../PIPS-IPM/Core/QpStoch)

set(EXECUTABLE_OUTPUT_PATH "${PROJECT_BINARY_DIR}")

add_library(ucRolling ucRollingModel.cpp)


add_executable(ucRollingTest ucRollingTest.cpp)
add_executable(ucRollingCbcExtensive ucRollingCbcExtensive.cpp)
add_executable(ucRollingOOQPExtensive ucRollingOOQPExtensive.cpp)


target_link_libraries(ucRollingTest ucRolling pipss stochInput CbcBALPInterface ClpBALPInterface ${CBC_LIBS} ${CLP_LIB} ${COIN_LIBS} ${MATH_LIBS})
target_link_libraries(ucRollingCbcExtensive ucRolling pipss stochInput CbcBALPInterface ClpBALPInterface ${CBC_LIBS} ${CLP_LIB} ${COIN_LIBS} ${MATH_LIBS})
target_link_libraries(ucRollingOOQPExtensive ucRolling stochInput ${COIN_LIBS}
    ooqpgensparse ooqpbase ooqpmehrotra ooqpsparse ooqpdense
    ${MA57_LIBRARY} ${METIS_LIBRARY} ${MATH_LIBS} ${PARDISO_LIBRARY})


if(HAVE_SCIP)
	add_executable(ucRollingProxBundle ucRollingProxBundle.cpp)
	target_link_libraries(ucRollingProxBundle ucRolling bundle ClpBALPInterface ScipLagrangeSolver ${SCIP_LIBS} CbcRecourseSolver ClpRecourseSolver pipss stochInput ${CBC_LIBS} ${CLP_LIB} ${COIN_LIBS} ooqpgensparse ooqpbase ooqpmehrotra ooqpsparse ooqpdense ${MA57_LIBRARY} ${METIS_LIBRARY} ${MATH_LIBS})

endif(HAVE_SCIP)

if(HAVE_CB)
	include_directories("${CB_DIR}/include")
	add_executable(ucRollingBundle ucRollingBundle.cpp)
	target_link_libraries(ucRollingBundle ucRolling pipss stochInput CbcLagrangeSolver ClpRecourseSolver ${CBC_LIBS} ${CLP_LIB} ${COIN_LIBS} ${CB_LIB} ${MATH_LIBS})
endif(HAVE_CB)

