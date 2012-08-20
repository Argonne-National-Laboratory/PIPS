cmake_minimum_required(VERSION 2.8)

project(PIPSAll)

if (BUILD_SHARED_LIBS)
message(STATUS "Using shared libraries")
set(CMAKE_FIND_LIBRARY_SUFFIXES ".so;.a")
endif (BUILD_SHARED_LIBS)

#find_package(MPI REQUIRED)
if(NOT ${CMAKE_SYSTEM_NAME} MATCHES "BlueGeneP-static")
set(CMAKE_CXX_COMPILER "mpicxx")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0") # need this explicitly for mpicxx
endif(NOT ${CMAKE_SYSTEM_NAME} MATCHES "BlueGeneP-static")

add_definitions(-DTIMING -DSTOCH_TESTING)

# inspired by elemental-bgp
if(MATH_LIBS)
	set(NEED_MATH FALSE)
	message(STATUS "Using user-defined MATH_LIBS=${MATH_LIBS}")
else(MATH_LIBS)
	set(NEED_MATH TRUE)
endif(MATH_LIBS)

if(NEED_MATH)
	if(IS_DIRECTORY $ENV{MKLROOT})
		set(MATH_LIBS "-Wl,--start-group  $ENV{MKLROOT}/lib/intel64/libmkl_intel_lp64.a $ENV{MKLROOT}/lib/intel64/libmkl_sequential.a $ENV{MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lgfortran")
	else(IS_DIRECTORY $ENV{MKLROOT})
		enable_language( Fortran )
		find_package(LAPACK REQUIRED)
		set(MATH_LIBS ${LAPACK_LIBRARIES};${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES})
	endif(IS_DIRECTORY $ENV{MKLROOT})
endif(NEED_MATH)




set(SHARED_DIR ${CMAKE_CURRENT_SOURCE_DIR}/SharedLibraries)
set(SML_DIR ${SHARED_DIR}/sml-0.7.2)
set(AMPL_DIR ${SHARED_DIR}/amplsolvers)
set(MA27_DIR ${SHARED_DIR}/MA27)
set(MA57_DIR ${SHARED_DIR}/MA57)
set(MA86_DIR ${SHARED_DIR}/MA86/src)
set(METIS_DIR ${SHARED_DIR}/metis-4.0.1)
# include CoinUtils
#set(COIN_DIR ${Simplex_SOURCE_DIR}/CoinUtils-2.8.0)
#set(COIN_DIR ${SHARED_DIR}/Clp-1.14)
set(COIN_DIR ${SHARED_DIR}/Cbc-2.7.6)
set(PARDISO_DIR ${SHARED_DIR}/PARDISO)

if(NOT EXISTS "${COIN_DIR}/lib/pkgconfig/coinutils.pc")
	message(FATAL_ERROR "Please run \"./configure --enable-static; make install\" in SharedLibraries/Cbc-2.7.6")
endif(NOT EXISTS "${COIN_DIR}/lib/pkgconfig/coinutils.pc")

# if(NOT EXISTS "${MA27_DIR}/libMA27.a" AND NOT EXISTS "${MA57_DIR}/libMA57.a")
# 	message(FATAL_ERROR "Could not find MA27 or MA57. Please run \"make\" in SharedLibraries/MA27 or in SharedLibraries/MA57")
# endif(NOT EXISTS "${MA27_DIR}/libMA27.a" AND NOT EXISTS "${MA57_DIR}/libMA57.a")

if(EXISTS "${MA57_DIR}/libMA57.a" AND NOT EXISTS "${METIS_DIR}/libmetis.a")
  message(FATAL_ERROR "Could not find metis needed by MA57. Please run \"make\" in SharedLibraries/metis-4.0.1")
endif(EXISTS "${MA57_DIR}/libMA57.a" AND NOT EXISTS "${METIS_DIR}/libmetis.a")

set(ENV{PKG_CONFIG_PATH} ${COIN_DIR}/lib/pkgconfig)
find_package(PkgConfig REQUIRED)
pkg_check_modules(COINUTILS REQUIRED coinutils)
pkg_check_modules(CLP REQUIRED clp)
pkg_check_modules(CBC REQUIRED cbc)
find_library(CLP_LIB NAMES libClp.a Clp HINTS ${CLP_LIBRARY_DIRS})
set(COIN_LIBS "")
foreach(f ${COINUTILS_LIBRARIES})
	if(NOT ${f} MATCHES "m" AND NOT ${f} MATCHES "blas" AND NOT ${f} MATCHES "lapack") # skip libm
		find_library(X NAMES ${f} HINTS ${COINUTILS_LIBRARY_DIRS})
		set(COIN_LIBS ${COIN_LIBS} ${X})
		unset(X CACHE)
	endif(NOT ${f} MATCHES "m" AND NOT ${f} MATCHES "blas" AND NOT ${f} MATCHES "lapack")
endforeach(f)
#message(STATUS "${COIN_LIBS}")


find_library(CBCSolver_LIB NAMES libCbcSolver.a CbcSolver HINTS ${CBC_LIBRARY_DIRS})
find_library(CBC_LIB NAMES libCbc.a Cbc HINTS ${CBC_LIBRARY_DIRS})
find_library(OSICLP_LIB NAMES libOsiClp.a OsiClp HINTS ${CBC_LIBRARY_DIRS})
find_library(OSI_LIB NAMES libOsi.a Osi HINTS ${CBC_LIBRARY_DIRS})
find_library(CGL_LIB NAMES libCgl.a Cgl HINTS ${CBC_LIBRARY_DIRS})

set(CBC_LIBS ${CBCSolver_LIB} ${CBC_LIB} ${CGL_LIB} ${OSICLP_LIB} ${OSI_LIB})

find_library(AMPL_LIBRARY NAMES amplsolver.a PATHS ${AMPL_DIR})
find_library(SML_LIBRARY NAMES libsml.a PATHS ${SML_DIR}/lib)
find_library(DL_LIBRARY NAMES dl)
find_library(FL_LIBRARY NAMES fl)

find_library(MA27_LIBRARY MA27 PATHS ${MA27_DIR})
find_library(MA57_LIBRARY MA57 PATHS ${MA57_DIR})
find_library(METIS_LIBRARY NAME libmetis.a PATHS ${METIS_DIR})
find_library(PARDISO_LIBRARY name libpardiso490-GNU430-X86-64.so PATHS ${PARDISO_DIR})

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

#see if we have METIS
#if we don't, also set HAVE_MA57 to false.
set(HAVE_METIS TRUE)
if(${METIS_LIBRARY} MATCHES "METIS_LIBRARY-NOTFOUND")
  message(STATUS "metis library not found. Won't build code that requires METIS (i.e. QpStoch and code that requires MA57).")
  set(HAVE_METIS FALSE)
endif(${METIS_LIBRARY} MATCHES "METIS_LIBRARY-NOTFOUND")
if(HAVE_METIS)
  include_directories(${SHARED_DIR}/metis-4.0.1/Lib)
endif(HAVE_METIS)

set(HAVE_MA86 HAVE_METIS)
# Check if we have all libraries that Ma86 depends on exists.
if(${MA86_LIB} MATCHES "MA86_LIB-NOTFOUND")
  message(STATUS "MA86 library not found. Won't build code that requires MA86.")
  set(HAVE_MA86 FALSE)
endif(${MA86_LIB} MATCHES "MA86_LIB-NOTFOUND") 


# see if we have MA27 or MA57
set(HAVE_MA27 TRUE)
set(HAVE_MA57 HAVE_METIS)
#message(STATUS ${MA27_LIBRARY})
if(${MA27_LIBRARY} MATCHES "MA27_LIBRARY-NOTFOUND")
  message(STATUS "MA27 library not found. Won't build code that requires MA27.")
  set(HAVE_MA27 FALSE)
endif(${MA27_LIBRARY} MATCHES "MA27_LIBRARY-NOTFOUND")
if(HAVE_MA57 AND ${MA57_LIBRARY} MATCHES "MA57_LIBRARY-NOTFOUND")
  message(STATUS "MA57 library not found. Won't build code that requires MA57.")
  set(HAVE_MA57 FALSE)
endif(HAVE_MA57 AND ${MA57_LIBRARY} MATCHES "MA57_LIBRARY-NOTFOUND")


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

set(CB_DIR "${SHARED_DIR}/ConicBundle")
set(CB_LIB "${CB_DIR}/lib/libcb.a")
if(EXISTS ${CB_LIB})
	set(HAVE_CB TRUE)
else(EXISTS ${CB_LIB})
	set(HAVE_CB FALSE)
	message(STATUS "Building without ConicBundle library")
endif(EXISTS ${CB_LIB})

set(HAVE_PARDISO FALSE)
if(${CMAKE_SYSTEM_NAME} MATCHES "BlueGeneP-static")
	set(HAVE_PARDISO TRUE)
	set(PARDISO_LIB ${SHARED_DIR}/PARDISO/libpardiso.a ${SHARED_DIR}/PARDISO/libmetis41_pardiso.a ${SHARED_DIR}/PARDISO/libmetis41-P_pardiso.a ${SHARED_DIR}/PARDISO/libpils_pardiso.a ${SHARED_DIR}/PARDISO/libpardiso.a)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
else(${CMAKE_SYSTEM_NAME} MATCHES "BlueGeneP-static")
	set(PARDISO_LIB ${SHARED_DIR}/PARDISO/libpardiso490-GNU430-X86-64.so)
	if(EXISTS ${PARDISO_LIB})
		FIND_PACKAGE( OpenMP REQUIRED )
		set(HAVE_PARDISO TRUE)
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
		set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
	endif(EXISTS ${PARDISO_LIB})
endif(${CMAKE_SYSTEM_NAME} MATCHES "BlueGeneP-static")


set(MA86_DIR "${SHARED_DIR}/MA86")
set(MA86_LIB "${MA86_DIR}/lib/libhsl_ma86.a")
set(MA86_INC "${MA86_DIR}/include")

set(SCIP_DIR "${SHARED_DIR}/scip-2.1.1")
set(SCIP_LIBDIR "${SCIP_DIR}/lib")
set(SCIP_INCDIR "${SCIP_DIR}/src")
if(EXISTS "${SCIP_LIBDIR}/libscip.a")
	set(HAVE_SCIP TRUE)
	set(SCIP_LIBS ${SCIP_LIBDIR}/libscip.a;${SCIP_LIBDIR}/libobjscip.a;${SCIP_LIBDIR}/liblpiclp.a;${SCIP_LIBDIR}/libnlpi.a)
	if (NOT ${CMAKE_SYSTEM_NAME} MATCHES "BlueGeneP-static")
		set(SCIP_LIBS ${SCIP_LIBS};"readline";"ncurses")
	endif(NOT ${CMAKE_SYSTEM_NAME} MATCHES "BlueGeneP-static")
else(EXISTS "${SCIP_LIBDIR}/libscip.a")
	set(HAVE_SCIP FALSE)
	message(STATUS "Building without SCIP")
endif(EXISTS "${SCIP_LIBDIR}/libscip.a")


find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIR})


#set(CMAKE_BUILD_TYPE DEBUG)
set(CMAKE_BUILD_TYPE RELWITHDEBINFO)
#set(CMAKE_BUILD_TYPE RELEASE)
#message("FLAGS: ${CMAKE_CXX_FLAGS_RELEASE}")

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")

# must set TAU_MAKEFILE, TAU_OPTIONS environment variables
#set(CMAKE_CXX_COMPILER "tau_cxx.sh")

#IF(CMAKE_COMPILER_IS_GNUCXX)
#	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
#ENDIF(CMAKE_COMPILER_IS_GNUCXX)

include_directories(Input)
include_directories(SolverInterface)

include_directories(${COINUTILS_INCLUDE_DIRS})

add_subdirectory(PIPS-IPM)
add_subdirectory(Input)
add_subdirectory(SolverInterface)
add_subdirectory(PIPS-S)
add_subdirectory(Lagrange)
add_subdirectory(Utilities)
add_subdirectory(UCTests)
add_subdirectory(UCRolling)
add_subdirectory(BasisBootstrap)

# this will print a warning message, should figure out how to silence it
#add_library(pipsForSwift SHARED)
# add any libraries needed by swift app here
#target_link_libraries(pipsForSwift pipss stochInput ClpRecourseSolver ${CLP_LIB} ${COIN_LIBS} ${MATH_LIBS})

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

