include_directories(../PIPS-S/Basic)
include_directories(../PIPS-S/Core)
include_directories(../PIPS-S/CoinBALPFactorization)
include_directories(../LagrangeBnC/RecourseSubproblemSolver)

set(EXECUTABLE_OUTPUT_PATH "${PROJECT_BINARY_DIR}")

add_executable(simpleBasisBootstrapClp simpleBasisBootstrapClp.cpp)

target_link_libraries(simpleBasisBootstrapClp pipss stochInput ClpBALPInterface ClpRecourseSolver ${CBC_LIBS} ${CLP_LIB} ${COIN_LIBS} ${MATH_LIBS})
