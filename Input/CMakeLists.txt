add_library(stochInput rawInput.cpp SMPSInput.cpp combinedInput.cpp)

if(HAVE_SML)
	add_library(SMLStochInput SMLInput.cpp)
endif(HAVE_SML)
