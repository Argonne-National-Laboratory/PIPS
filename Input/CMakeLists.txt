add_library(stochInput rawInput.cpp SMPSInput.cpp)

if(HAVE_SML)
	add_library(SMLStochInput SMLInput.cpp)
endif(HAVE_SML)
