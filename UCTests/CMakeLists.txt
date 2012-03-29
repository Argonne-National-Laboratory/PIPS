include_directories(../PIPS-S/Basic)
include_directories(../LagrangeBnC/RecourseSubproblemSolver)
include_directories(../LagrangeBnC/LagrangeSubproblemSolver)

set(EXECUTABLE_OUTPUT_PATH "${PROJECT_BINARY_DIR}")

add_executable(ucCoefficientDiving ucCoefficientDiving.cpp)

target_link_libraries(ucCoefficientDiving pipss stochInput ClpBALPInterface ClpRecourseSolver ${CLP_LIB} ${COIN_LIBS} ${MATH_LIBS})


add_executable(cutsTest cutsTest.cpp OsiSubproblemWrapper.cpp)

target_link_libraries(cutsTest pipss stochInput ClpBALPInterface ClpRecourseSolver CbcLagrangeSolver ${CBC_LIBS} ${CLP_LIB} ${COIN_LIBS} ${MATH_LIBS})
