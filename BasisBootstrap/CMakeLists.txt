include_directories(../PIPS-S/Basic)
include_directories(../PIPS-S/Core)
include_directories(../PIPS-S/CoinBALPFactorization)
include_directories(../Lagrange/RecourseSubproblemSolver)
include_directories(ScenarioReduction)

set(EXECUTABLE_OUTPUT_PATH "${PROJECT_BINARY_DIR}")

#add_executable(simpleBasisBootstrapClp simpleBasisBootstrapClp.cpp)
#add_executable(scenRedBasisBootstrapClp basisBootstrapScenRedClp.cpp)
add_library(scenred ScenarioReduction/fastForwardSelection.cpp ScenarioReduction/scenarioReductionUtilities.cpp)

#target_link_libraries(simpleBasisBootstrapClp pipss stochInput ClpBALPInterface ClpRecourseSolver ${CBC_LIBS} ${CLP_LIB} ${COIN_LIBS} ${MATH_LIBS})
#target_link_libraries(scenRedBasisBootstrapClp pipss scenred stochInput ClpBALPInterface ClpRecourseSolver ${CBC_LIBS} ${CLP_LIB} ${COIN_LIBS} ${MATH_LIBS})
