
if(HAVE_UMFPACK) 
	add_library(reducedSolverSpecial ReducedSpaceSolverStateOnly.C) 
else(HAVE_UMFPACK) 
	add_library(reducedSolverSpecial empty.cpp) 
endif(HAVE_UMFPACK)


