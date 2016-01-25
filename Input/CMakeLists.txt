add_library(stochInput rawInput.cpp SMPSInput.cpp combinedInput.cpp stochasticInput.cpp)

add_library(nlpstochInput rawInput.cpp combinedInput.cpp stochasticInput.cpp)

add_library(genNLparallel amplGenStochInput.cpp amplGenStochInput_AddSlack.cpp stochasticInput.cpp AmplData_NL.cpp ${AMPL_LIBRARY})
target_include_directories(genNLparallel PUBLIC ${AMPL_INCLUDE_DIR})

add_library(genStructparallel StructJuMPInput.cpp stochasticInput.cpp)

add_library(genNLserial  AmplData_NL.cpp ${AMPL_LIBRARY})
target_include_directories(genNLserial PUBLIC ${AMPL_INCLUDE_DIR})

#add_library(multiStageInput multiStageInputTree.cpp stochasticInput.cpp)

add_library(nlpstochInputBase stochasticInput.cpp)
add_subdirectory(OPF_Matpower)

