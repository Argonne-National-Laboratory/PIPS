

include_directories(../PIPS-S/Basic)
set(EXECUTABLE_OUTPUT_PATH "${PROJECT_BINARY_DIR}")

add_library(ucRolling ucRollingModel.cpp)


add_executable(ucRollingTest ucRollingTest.cpp)


target_link_libraries(ucRollingTest ucRolling pipss stochInput CbcBALPInterface ClpBALPInterface ${CBC_LIBS} ${CLP_LIB} ${COIN_LIBS} ${MATH_LIBS})
