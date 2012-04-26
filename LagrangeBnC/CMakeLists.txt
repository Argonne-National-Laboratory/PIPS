include_directories(LagrangeSubproblemSolver)
include_directories(RecourseSubproblemSolver)
include_directories(../PIPS-S/Basic)
include_directories(../BasisBootstrap/ScenarioReduction)

add_library(CbcLagrangeSolver LagrangeSubproblemSolver/CbcLagrangeSolver.cpp)
add_library(ClpRecourseSolver RecourseSubproblemSolver/ClpRecourseSolver.cpp)
add_library(CbcRecourseSolver RecourseSubproblemSolver/CbcRecourseSolver.cpp)


set(EXECUTABLE_OUTPUT_PATH "${PROJECT_BINARY_DIR}")

add_executable(lagrangeRootNode Drivers/lagrangeRootNode.cpp)
add_executable(lagrangeCombinedScenRedRootNode Drivers/lagrangeCombinedScenRedRootNode.cpp)
add_executable(lagrangeCombinedSMPS Drivers/lagrangeCombinedScenRedRootNodeSMPS.cpp)
add_executable(lagrangeRootNodeSMPS Drivers/lagrangeRootNodeSMPS.cpp)

target_link_libraries(lagrangeRootNode pipss stochInput ClpBALPInterface CbcLagrangeSolver ClpRecourseSolver ${CBC_LIBS} ${CLP_LIB} ${COIN_LIBS} ${MATH_LIBS})
target_link_libraries(lagrangeCombinedScenRedRootNode pipss stochInput scenred ClpBALPInterface CbcLagrangeSolver ClpRecourseSolver ${CBC_LIBS} ${CLP_LIB} ${COIN_LIBS} ${MATH_LIBS})
target_link_libraries(lagrangeCombinedSMPS pipss stochInput scenred ClpBALPInterface CbcLagrangeSolver CbcRecourseSolver ${CBC_LIBS} ${CLP_LIB} ${COIN_LIBS} ${MATH_LIBS})
target_link_libraries(lagrangeRootNodeSMPS pipss stochInput scenred ClpBALPInterface CbcLagrangeSolver CbcRecourseSolver ${CBC_LIBS} ${CLP_LIB} ${COIN_LIBS} ${MATH_LIBS})
