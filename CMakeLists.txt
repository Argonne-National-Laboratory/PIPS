cmake_minimum_required(VERSION 2.8)

project(PIPSAll)

#find_package(MPI REQUIRED)
if(NOT ${CMAKE_SYSTEM_NAME} MATCHES "BlueGeneP-static") 
set(CMAKE_CXX_COMPILER "mpicxx")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0") # need this explicitly for mpicxx
endif(NOT ${CMAKE_SYSTEM_NAME} MATCHES "BlueGeneP-static") 

# inspired by elemental-bgp
if(MATH_LIBS)
	set(NEED_MATH FALSE)
	message(STATUS "Using user-defined MATH_LIBS=${MATH_LIBS}")
else(MATH_LIBS)
	set(NEED_MATH TRUE)
endif(MATH_LIBS)

if(NEED_MATH)
	if(IS_DIRECTORY $ENV{MKLROOT})
		set(MATH_LIBS "-Wl,--start-group  $ENV{MKLROOT}/lib/intel64/libmkl_intel_lp64.a $ENV{MKLROOT}/lib/intel64/libmkl_gnu_thread.a $ENV{MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -fopenmp -lpthread")
	else(IS_DIRECTORY $ENV{MKLROOT})
		enable_language( Fortran )
		find_package(LAPACK REQUIRED)
		set(MATH_LIBS "${LAPACK_LIBRARIES}")
	endif(IS_DIRECTORY $ENV{MKLROOT})
endif(NEED_MATH)





set(SHARED_DIR ${CMAKE_CURRENT_SOURCE_DIR}/SharedLibraries)
set(SML_DIR ${SHARED_DIR}/sml-0.7.2)
set(AMPL_DIR ${SHARED_DIR}/amplsolvers)

# include CoinUtils
#set(COIN_DIR ${Simplex_SOURCE_DIR}/CoinUtils-2.8.0)
#set(COIN_DIR ${SHARED_DIR}/Clp-1.14)
set(COIN_DIR ${SHARED_DIR}/Cbc-2.7.6)

if(NOT EXISTS "${COIN_DIR}/lib/pkgconfig/coinutils.pc")
	message(FATAL_ERROR "Please run \"./configure --enable-static; make install\" in Clp-1.14")
endif(NOT EXISTS "${COIN_DIR}/lib/pkgconfig/coinutils.pc")

set(ENV{PKG_CONFIG_PATH} ${COIN_DIR}/lib/pkgconfig)
find_package(PkgConfig REQUIRED)
pkg_check_modules(COINUTILS REQUIRED coinutils)
pkg_check_modules(CLP REQUIRED clp)
find_library(CLP_LIB NAMES libClp.a Clp HINTS ${CLP_LIBRARY_DIRS})
set(COIN_LIBS "")
foreach(f ${COINUTILS_LIBRARIES})
	if(NOT ${f} MATCHES "m") # skip libm
		find_library(X NAMES ${f} HINTS ${COINUTILS_LIBRARY_DIRS})
		set(COIN_LIBS ${COIN_LIBS} ${X})
		unset(X CACHE)
	endif(NOT ${f} MATCHES "m")
endforeach(f)
#message(STATUS "${COIN_LIBS}")



find_library(AMPL_LIBRARY NAMES amplsolver.a PATHS ${AMPL_DIR})
find_library(SML_LIBRARY NAMES libsml.a PATHS ${SML_DIR}/lib)
find_library(DL_LIBRARY NAMES dl)
find_library(FL_LIBRARY NAMES fl)
set(HAVE_SML TRUE)
if(${AMPL_LIBRARY} MATCHES "AMPL_LIBRARY-NOTFOUND")
	message(STATUS "amplsolver library not found. Won't build code that requires SML.")
	set(HAVE_SML FALSE)
endif(${AMPL_LIBRARY} MATCHES "AMPL_LIBRARY-NOTFOUND")
if(${SML_LIBRARY} MATCHES "SML_LIBRARY-NOTFOUND" AND ${HAVE_SML})
	message(STATUS "SML library not found. Won't build code that requires SML.")
	set(HAVE_SML FALSE)
endif(${SML_LIBRARY} MATCHES "SML_LIBRARY-NOTFOUND" AND ${HAVE_SML})
if(${DL_LIBRARY} MATCHES "DL_LIBRARY-NOTFOUND" AND ${HAVE_SML})
	message(STATUS "dl library not found. Won't build code that requires SML.")
	set(HAVE_SML FALSE)
endif(${DL_LIBRARY} MATCHES "DL_LIBRARY-NOTFOUND" AND ${HAVE_SML})
if(${FL_LIBRARY} MATCHES "FL_LIBRARY-NOTFOUND" AND ${HAVE_SML})
	message(STATUS "flex library not found. Won't build code that requires SML.")
	set(HAVE_SML FALSE)
endif(${FL_LIBRARY} MATCHES "FL_LIBRARY-NOTFOUND" AND ${HAVE_SML})

if(HAVE_SML)
	include_directories(${SML_DIR}/include)
endif(HAVE_SML)


# see if we have gurobi
# must set GUROBI_DIR environmental variable with root of gurobi architecture directory
# e.g. gurobi461/linux64
if(IS_DIRECTORY $ENV{GUROBI_DIR})
	find_library(GUROBIC NAMES gurobi46 HINTS "$ENV{GUROBI_DIR}/lib")
	find_library(GUROBICPP NAMES "gurobi_c++" HINTS "$ENV{GUROBI_DIR}/lib")
	find_package (Threads)
	set(HAVE_GUROBI TRUE)
	set(GUROBI_LIBS ${GUROBIC} ${GUROBICPP} ${CMAKE_THREAD_LIBS_INIT})
	message(STATUS "Building with Gurobi")
else(IS_DIRECTORY $ENV{GUROBI_DIR})
	set(HAVE_GUROBI FALSE)
	message(STATUS "Building without Gurobi")
endif(IS_DIRECTORY $ENV{GUROBI_DIR})



find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIR})


#set(CMAKE_BUILD_TYPE DEBUG)
set(CMAKE_BUILD_TYPE RELWITHDEBINFO)
#set(CMAKE_BUILD_TYPE RELEASE)
#message("FLAGS: ${CMAKE_CXX_FLAGS_RELEASE}")

# must set TAU_MAKEFILE, TAU_OPTIONS environment variables
#set(CMAKE_CXX_COMPILER "tau_cxx.sh") 

#IF(CMAKE_COMPILER_IS_GNUCXX)
#	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
#ENDIF(CMAKE_COMPILER_IS_GNUCXX)

include_directories(Input)
include_directories(SolverInterface)

include_directories(${COINUTILS_INCLUDE_DIRS})

add_subdirectory(Input)
add_subdirectory(SolverInterface)
add_subdirectory(PIPS-S)


#if(NOT ${CMAKE_SYSTEM_NAME} MATCHES "BlueGeneP-static") 
#add_subdirectory(src/Tests)
#add_subdirectory(src/SuboptimizationSerial)
#add_subdirectory(src/Utilities)
#add_subdirectory(src/SMPS)
#endif(NOT ${CMAKE_SYSTEM_NAME} MATCHES "BlueGeneP-static") 


#if(HAVE_SML)
#	add_subdirectory(src/SML)
#endif(HAVE_SML)
#if(HAVE_GUROBI)
#	add_subdirectory(src/GUROBI)
#endif(HAVE_GUROBI)

