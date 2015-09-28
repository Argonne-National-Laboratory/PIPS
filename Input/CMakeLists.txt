add_library(stochInput rawInput.cpp SMPSInput.cpp combinedInput.cpp stochasticInput.cpp)

if(HAVE_SML)
	add_library(SMLStochInput SMLInput.cpp)
endif(HAVE_SML)

add_library(genNLparallel amplGenStochInput.cpp amplGenStochInput_AddSlack.cpp stochasticInput.cpp AmplData_NL.cpp) 
add_library(genNLserial  AmplData_NL.cpp)


include_directories(OPF_Matpower) 
add_subdirectory(OPF_Matpower)
add_library(dcopfcinput dcopflowInput.cpp)