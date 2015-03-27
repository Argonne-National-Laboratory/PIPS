include_directories(../PIPS-S/Basic)
include_directories(../Lagrange/RecourseSubproblemSolver)
include_directories(../Lagrange/LagrangeSubproblemSolver)

set(EXECUTABLE_OUTPUT_PATH "${PROJECT_BINARY_DIR}")

add_executable(ucCoefficientDiving ucCoefficientDiving.cpp)

target_link_libraries(ucCoefficientDiving pipss stochInput ClpBALPInterface ClpRecourseSolver ${CLP_LIB} ${COIN_LIBS} ${MATH_LIBS})


add_executable(cutsTest cutsTest.cpp OsiSubproblemWrapper.cpp)

target_link_libraries(cutsTest pipss stochInput ClpBALPInterface ClpRecourseSolver CbcLagrangeSolver ${CBC_LIBS} ${CLP_LIB} ${COIN_LIBS} ${MATH_LIBS})


add_executable(ucRootNodeClp ucRootNodeClp.cpp OsiSubproblemWrapper.cpp)
target_link_libraries(ucRootNodeClp pipss stochInput ClpBALPInterface ClpRecourseSolver ${CBC_LIBS} ${CLP_LIB} ${COIN_LIBS} ${MATH_LIBS})


add_executable(ucTestSol ucTestSol.cpp)
target_link_libraries(ucTestSol pipss stochInput ClpRecourseSolver ${CBC_LIBS} ${CLP_LIB} ${COIN_LIBS} ${MATH_LIBS})

include_directories(../PIPS-S/Core)
include_directories(../PIPS-S/CoinBALPFactorization)
add_executable(ucRootNodePIPSS ucRootNodePIPSS.cpp OsiSubproblemWrapper.cpp)
target_link_libraries(ucRootNodePIPSS pipss stochInput ClpBALPInterface ClpRecourseSolver ${CBC_LIBS} ${CLP_LIB} ${COIN_LIBS} ${MATH_LIBS} ${Boost_LIBRARIES})

