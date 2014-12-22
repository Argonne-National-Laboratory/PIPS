#project(PIPS-S) # this doesn't work for some reason, cmake starts barfing

include_directories(CoinBALPFactorization)
include_directories(Basic)
include_directories(Core)
include_directories(Drivers)

add_library(pipss Core/BALinearAlgebra.cpp Core/BALPSolverBase.cpp Core/BALPSolverDual.cpp Core/BALPSolverPrimal.cpp Core/BAPFIPar.cpp Core/PIPSSInterface.cpp
	Basic/denseVector.cpp Basic/sparseVector.cpp Basic/BA.cpp Basic/BAVector.cpp Basic/BAData.cpp
	CoinBALPFactorization/CoinBALPFactorization1.cpp CoinBALPFactorization/CoinBALPFactorization2.cpp 
	CoinBALPFactorization/CoinBALPFactorization3.cpp CoinBALPFactorization/CoinBALPFactorization4.cpp
	Basic/PIPSLogging.cpp)	
	
set(EXECUTABLE_OUTPUT_PATH "${PROJECT_BINARY_DIR}")

add_executable(pipssFromRaw Drivers/pipssFromRaw.cpp)
target_link_libraries(pipssFromRaw pipss stochInput ${COIN_LIBS} ${MATH_LIBS} ${Boost_LIBRARIES})

add_executable(pipssSMPS Drivers/pipssSMPS.cpp)
target_link_libraries(pipssSMPS pipss stochInput ${COIN_LIBS} ${MATH_LIBS} ${Boost_LIBRARIES})

add_executable(pipssmemleak Drivers/memleak.cpp)
target_link_libraries(pipssmemleak pipss stochInput ${COIN_LIBS} ${MATH_LIBS} ${Boost_LIBRARIES})

add_executable(clpFromRaw Drivers/clpFromRaw.cpp)
target_link_libraries(clpFromRaw pipss stochInput ClpBALPInterface ${CLP_LIB} ${COIN_LIBS} ${MATH_LIBS} ${Boost_LIBRARIES})

add_executable(clpSMPS Drivers/clpSMPS.cpp)
target_link_libraries(clpSMPS pipss stochInput ClpBALPInterface ${CLP_LIB} ${COIN_LIBS} ${MATH_LIBS} ${Boost_LIBRARIES})

add_executable(pipssBootstrap Drivers/pipssBootstrap.cpp)
target_link_libraries(pipssBootstrap pipss stochInput ${CLP_LIB} ${COIN_LIBS} ${MATH_LIBS} ${Boost_LIBRARIES})

add_executable(clpBootstrap Drivers/clpBootstrap.cpp)
target_link_libraries(clpBootstrap pipss stochInput ClpBALPInterface ${CLP_LIB} ${COIN_LIBS} ${MATH_LIBS})

#add_subdirectory(CoinBALPFactorization)
#add_subdirectory(Basic)
#add_subdirectory(Core)


