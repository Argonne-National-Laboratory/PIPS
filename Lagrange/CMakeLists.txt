include_directories(LagrangeSubproblemSolver)
include_directories(RecourseSubproblemSolver)
include_directories(Bundle)
include_directories(../PIPS-S/Basic)
include_directories(../PIPS-S/Core)
include_directories(../PIPS-S/CoinBALPFactorization)
include_directories(../BasisBootstrap/ScenarioReduction)
include_directories(../PIPS-IPM/Core/Abstract ../PIPS-IPM/Core/Vector ../PIPS-IPM/Core/Utilities ../PIPS-IPM/Core/QpSolvers ../PIPS-IPM/Core/QpGen
  ../PIPS-IPM/Core/SparseLinearAlgebra ../PIPS-IPM/Core/DenseLinearAlgebra ../PIPS-IPM/Core/Readers 
  ../PIPS-IPM/Core/LinearSolvers/Ma57Solver)
include_directories(../PIPS-IPM/Core/StochLinearAlgebra ../PIPS-IPM/Core/QpStoch)

add_library(CbcLagrangeSolver LagrangeSubproblemSolver/CbcLagrangeSolver.cpp)
add_library(ClpRecourseSolver RecourseSubproblemSolver/ClpRecourseSolver.cpp)
add_library(CbcRecourseSolver RecourseSubproblemSolver/CbcRecourseSolver.cpp)

add_library(combine Drivers/combineScenarios.cpp)
add_library(bundle Bundle/cuttingPlaneBALP.cpp Bundle/levelBALP.cpp
	Bundle/l1bundleBALP.cpp Bundle/lInfTrustBALP.cpp Bundle/l1TrustBALP.cpp
	Bundle/proximalBAQP.cpp)

set(EXECUTABLE_OUTPUT_PATH "${PROJECT_BINARY_DIR}")

add_executable(lagrangeRootNode Drivers/lagrangeRootNode.cpp)
add_executable(lagrangeCombinedScenRedRootNode Drivers/lagrangeCombinedScenRedRootNode.cpp)
add_executable(lagrangeCombinedSMPS Drivers/lagrangeCombinedScenRedRootNodeSMPS.cpp)
add_executable(lagrangeRootNodeSMPS Drivers/lagrangeRootNodeSMPS.cpp)
add_executable(cpmSMPS Drivers/cpmSMPS.cpp)
add_executable(cpmRaw Drivers/cpmRaw.cpp)
add_executable(levelSMPS Drivers/levelSMPS.cpp)
add_executable(proxLevelSMPS Drivers/proxLevelSMPS.cpp)
add_executable(proxL1SMPS Drivers/l1bundleSMPS.cpp)
add_executable(proxLinfTrustSMPS Drivers/linftrustSMPS.cpp)
add_executable(proxL1TrustSMPS Drivers/l1trustSMPS.cpp)
add_executable(proxLinfTrustRaw Drivers/linftrustRaw.cpp)
add_executable(proxBundleSMPS Drivers/proxBundleSMPS.cpp)
add_executable(proxBundleRaw Drivers/proxBundleRaw.cpp)
add_executable(testSols Drivers/testSols.cpp)

target_link_libraries(lagrangeRootNode pipss stochInput ClpBALPInterface CbcLagrangeSolver ClpRecourseSolver ${CBC_LIBS} ${CLP_LIB} ${COIN_LIBS} ${MATH_LIBS})
target_link_libraries(lagrangeCombinedScenRedRootNode pipss combine stochInput scenred ClpBALPInterface CbcLagrangeSolver ClpRecourseSolver ${CBC_LIBS} ${CLP_LIB} ${COIN_LIBS} ${MATH_LIBS})
target_link_libraries(lagrangeCombinedSMPS pipss combine stochInput scenred ClpBALPInterface CbcLagrangeSolver CbcRecourseSolver ${CBC_LIBS} ${CLP_LIB} ${COIN_LIBS} ${MATH_LIBS})
target_link_libraries(lagrangeRootNodeSMPS pipss stochInput scenred ClpBALPInterface CbcLagrangeSolver CbcRecourseSolver ${CBC_LIBS} ${CLP_LIB} ${COIN_LIBS} ${MATH_LIBS})
target_link_libraries(cpmSMPS bundle pipss stochInput scenred ClpBALPInterface CbcLagrangeSolver CbcRecourseSolver ClpRecourseSolver ${CBC_LIBS} ${CLP_LIB} ${COIN_LIBS} ${MATH_LIBS})
target_link_libraries(cpmRaw bundle pipss stochInput scenred ClpBALPInterface CbcLagrangeSolver CbcRecourseSolver ClpRecourseSolver ${CBC_LIBS} ${CLP_LIB} ${COIN_LIBS} ${MATH_LIBS})
target_link_libraries(levelSMPS bundle pipss stochInput scenred ClpBALPInterface CbcLagrangeSolver CbcRecourseSolver ClpRecourseSolver ${CBC_LIBS} ${CLP_LIB} ${COIN_LIBS} ${MATH_LIBS})
target_link_libraries(proxLevelSMPS bundle pipss stochInput scenred ClpBALPInterface CbcLagrangeSolver CbcRecourseSolver ClpRecourseSolver ${CBC_LIBS} ${CLP_LIB} ${COIN_LIBS} ${MATH_LIBS})
target_link_libraries(proxL1SMPS bundle pipss stochInput scenred ClpBALPInterface CbcLagrangeSolver CbcRecourseSolver ClpRecourseSolver ${CBC_LIBS} ${CLP_LIB} ${COIN_LIBS} ${MATH_LIBS})
target_link_libraries(proxLinfTrustSMPS bundle pipss stochInput scenred ClpBALPInterface CbcLagrangeSolver CbcRecourseSolver ClpRecourseSolver ${CBC_LIBS} ${CLP_LIB} ${COIN_LIBS} ${MATH_LIBS})
target_link_libraries(proxL1TrustSMPS bundle pipss stochInput scenred ClpBALPInterface CbcLagrangeSolver CbcRecourseSolver ClpRecourseSolver ${CBC_LIBS} ${CLP_LIB} ${COIN_LIBS} ${MATH_LIBS})
target_link_libraries(proxLinfTrustRaw bundle pipss stochInput scenred ClpBALPInterface CbcLagrangeSolver CbcRecourseSolver ClpRecourseSolver ${CBC_LIBS} ${CLP_LIB} ${COIN_LIBS} ${MATH_LIBS})
target_link_libraries(proxBundleSMPS bundle pipss ooqpgensparse ooqpbase ooqpmehrotra ooqpsparse ooqpdense
    ${MA57_LIBRARY} ${METIS_LIBRARY} stochInput scenred ClpBALPInterface CbcLagrangeSolver CbcRecourseSolver ClpRecourseSolver ${CBC_LIBS} ${CLP_LIB} ${COIN_LIBS} ${MATH_LIBS})
target_link_libraries(proxBundleRaw bundle pipss ooqpgensparse ooqpbase ooqpmehrotra ooqpsparse ooqpdense
    ${MA57_LIBRARY} ${METIS_LIBRARY} stochInput scenred ClpBALPInterface CbcLagrangeSolver CbcRecourseSolver ClpRecourseSolver ${CBC_LIBS} ${CLP_LIB} ${COIN_LIBS} ${MATH_LIBS})
target_link_libraries(testSols bundle pipss stochInput scenred ClpBALPInterface CbcLagrangeSolver CbcRecourseSolver ClpRecourseSolver ${CBC_LIBS} ${CLP_LIB} ${COIN_LIBS} ${MATH_LIBS})

if(HAVE_SCIP)
	include_directories("${SCIP_INCDIR}")
	add_executable(testSolsHeur1 Drivers/testSolsHeur1.cpp)
	add_library(ScipLagrangeSolver LagrangeSubproblemSolver/ScipLagrangeSolver.cpp)
	add_executable(proxBundleScipOOQPSMPS Drivers/proxBundleScipOOQPSMPS.cpp)
	add_executable(proxBundleScipPIPSSMPS Drivers/proxBundleScipPIPSSMPS.cpp)
	add_executable(proxBundleScipRaw Drivers/proxBundleScipRaw.cpp)
	target_link_libraries(proxBundleScipPIPSSMPS bundle pipss ooqpstoch ooqpstochla ooqpmehrotrastoch ooqpgensparse ooqpbase ooqpmehrotra ooqpsparse ooqpdense
    ${MA57_LIBRARY} ${METIS_LIBRARY} stochInput scenred ClpBALPInterface ScipLagrangeSolver ${SCIP_LIBS} CbcRecourseSolver ClpRecourseSolver ${CBC_LIBS} ${CLP_LIB} ${COIN_LIBS} ${MATH_LIBS})
	target_link_libraries(proxBundleScipOOQPSMPS bundle pipss ooqpgensparse ooqpbase ooqpmehrotra ooqpsparse ooqpdense
    ${MA57_LIBRARY} ${METIS_LIBRARY} stochInput scenred ClpBALPInterface ScipLagrangeSolver ${SCIP_LIBS} CbcRecourseSolver ClpRecourseSolver ${CBC_LIBS} ${CLP_LIB} ${COIN_LIBS} ${MATH_LIBS})
	target_link_libraries(proxBundleScipRaw bundle pipss ooqpstoch ooqpstochla ooqpmehrotrastoch ooqpgensparse ooqpbase ooqpmehrotra ooqpsparse ooqpdense
    ${MA57_LIBRARY} ${METIS_LIBRARY} stochInput scenred ClpBALPInterface ScipLagrangeSolver ${SCIP_LIBS} CbcRecourseSolver ClpRecourseSolver ${CBC_LIBS} ${CLP_LIB} ${COIN_LIBS} ${MATH_LIBS})
	target_link_libraries(testSolsHeur1 bundle pipss ooqpgensparse ooqpbase ooqpmehrotra ooqpsparse ooqpdense
    ${MA57_LIBRARY} ${METIS_LIBRARY} stochInput scenred ClpBALPInterface ScipLagrangeSolver ${SCIP_LIBS} CbcRecourseSolver ClpRecourseSolver ${CBC_LIBS} ${CLP_LIB} ${COIN_LIBS} ${MATH_LIBS})
endif(HAVE_SCIP)

if(HAVE_CB AND HAVE_SCIP)
	include_directories("${CB_DIR}/include")
	add_executable(conicBundleSMPS Drivers/conicBundleDriverSMPS.cpp)
	add_executable(conicBundleRaw Drivers/conicBundleDriverRaw.cpp)
	add_executable(conicBundleParRaw Drivers/conicBundleParDriver1Raw.cpp)
	add_executable(conicBundleParSMPS Drivers/conicBundleParDriver1SMPS.cpp)
	add_executable(conicBundleCombinedSMPS Drivers/conicBundleDriverCombinedSMPS.cpp)
	target_link_libraries(conicBundleSMPS pipss stochInput scenred ClpBALPInterface ScipLagrangeSolver CbcRecourseSolver ${SCIP_LIBS} ${CLP_LIB} ${COIN_LIBS} ${CB_LIB} ${MATH_LIBS})
	target_link_libraries(conicBundleRaw pipss stochInput scenred ClpBALPInterface ScipLagrangeSolver ${SCIP_LIBS} CbcRecourseSolver ${CLP_LIB} ${COIN_LIBS} ${CB_LIB} ${MATH_LIBS})
	target_link_libraries(conicBundleParRaw pipss stochInput scenred ClpBALPInterface CbcLagrangeSolver CbcRecourseSolver ${CBC_LIBS} ${CLP_LIB} ${COIN_LIBS} ${CB_LIB} ${MATH_LIBS})
	target_link_libraries(conicBundleParSMPS pipss stochInput scenred ClpBALPInterface ScipLagrangeSolver ${SCIP_LIBS} CbcRecourseSolver ${CBC_LIBS} ${CLP_LIB} ${COIN_LIBS} ${CB_LIB} ${MATH_LIBS})
	target_link_libraries(conicBundleCombinedSMPS pipss combine stochInput scenred ClpBALPInterface CbcLagrangeSolver CbcRecourseSolver ${CBC_LIBS} ${CLP_LIB} ${COIN_LIBS} ${CB_LIB} ${MATH_LIBS})
endif(HAVE_CB AND HAVE_SCIP)

