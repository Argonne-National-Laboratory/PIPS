

include_directories(../PIPS-S/Basic)
include_directories(../PIPS-S/Core)
include_directories(../PIPS-S/CoinBALPFactorization)
include_directories(../LagrangeBnC/LagrangeSubproblemSolver)
include_directories(../LagrangeBnC/RecourseSubproblemSolver)
include_directories(../LagrangeBnC/Drivers)
set(EXECUTABLE_OUTPUT_PATH "${PROJECT_BINARY_DIR}")

add_library(ucRolling ucRollingModel.cpp)


add_executable(ucRollingTest ucRollingTest.cpp)
add_executable(ucRollingCbcExtensive ucRollingCbcExtensive.cpp)


target_link_libraries(ucRollingTest ucRolling pipss stochInput CbcBALPInterface ClpBALPInterface ${CBC_LIBS} ${CLP_LIB} ${COIN_LIBS} ${MATH_LIBS})
target_link_libraries(ucRollingCbcExtensive ucRolling pipss stochInput CbcBALPInterface ClpBALPInterface ${CBC_LIBS} ${CLP_LIB} ${COIN_LIBS} ${MATH_LIBS})
if(HAVE_CB)
	include_directories("${CB_DIR}/include")
	add_executable(ucRollingBundle ucRollingBundle.cpp)
	target_link_libraries(ucRollingBundle ucRolling pipss stochInput CbcLagrangeSolver ClpRecourseSolver ${CBC_LIBS} ${CLP_LIB} ${COIN_LIBS} ${CB_LIB} ${MATH_LIBS})
endif(HAVE_CB)

