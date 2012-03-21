include_directories(LagrangeSubproblemSolver)
include_directories(RecourseSubproblemSolver)
include_directories(../PIPS-S/Basic)

add_library(CbcLagrangeSolver LagrangeSubproblemSolver/CbcLagrangeSolver.cpp)
add_library(ClpRecourseSolver RecourseSubproblemSolver/ClpRecourseSolver.cpp)


set(EXECUTABLE_OUTPUT_PATH "${PROJECT_BINARY_DIR}")

add_executable(lagrangeRootNodeSerial Drivers/lagrangeRootNodeSerial.cpp)

target_link_libraries(lagrangeRootNodeSerial pipss stochInput ClpBALPInterface CbcLagrangeSolver ClpRecourseSolver ${CBC_LIBS} ${CLP_LIB} ${COIN_LIBS} ${MATH_LIBS})
