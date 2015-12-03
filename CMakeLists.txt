project(PIPSAll)

# include different "whole archive" linking options depending on compiler
if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  set(WHOLE_ARCHIVE "-Wl,-all_load")
  set(NO_WHOLE_ARCHIVE "-Wl,-noall_load")
  message(STATUS "SETTING HAVE_CLANG")
  set(HAVE_CLANG 1)
else ()
  set(WHOLE_ARCHIVE "-Wl,--whole-archive")
  set(NO_WHOLE_ARCHIVE "-Wl,--no-whole-archive")
endif (CMAKE_CXX_COMPILER_ID MATCHES "Clang")

cmake_minimum_required(VERSION 2.8)
if(COMMAND cmake_policy)
	cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

if (BUILD_SHARED_LIBS)
	message(STATUS "Using shared libraries")
	set(CMAKE_FIND_LIBRARY_SUFFIXES ".so;.a")
endif (BUILD_SHARED_LIBS)

find_package(MPI REQUIRED)
message(STATUS "CMAKE_SYSTEM_NAME ${CMAKE_SYSTEM_NAME}")
if(${CMAKE_SYSTEM_NAME} MATCHES "BlueGeneP-static")
	set(CMAKE_FIND_LIBRARY_SUFFIXES ".so;.a") #need this to find PARDISO
	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,-dy")
else(${CMAKE_SYSTEM_NAME} MATCHES "BlueGeneP-static")
  if(${CMAKE_SYSTEM_NAME} MATCHES "BlueGeneQ-static")
	set(CMAKE_FIND_LIBRARY_SUFFIXES ".so;.a") #need this to find PARDISO
	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,-dy")
  else(${CMAKE_SYSTEM_NAME} MATCHES "BlueGeneQ-static")
        if(IS_XE6 OR IS_XK7 OR IS_XC30)
          #set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -dynamic")
			else(IS_XE6 OR IS_XK7 OR IS_XC30)
	    	set(CMAKE_CXX_COMPILER "mpicxx")
	    	set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O02") # need this explicitly for mpicxx
	    	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -fPIC -fpermissive")
        endif(IS_XE6 OR IS_XK7 OR IS_XC30)
  endif(${CMAKE_SYSTEM_NAME} MATCHES "BlueGeneQ-static")
endif(${CMAKE_SYSTEM_NAME} MATCHES "BlueGeneP-static")


#add_definitions(-DTIMING -DSTOCH_TESTING)

# inspired by elemental-bgp
if(MATH_LIBS)
	set(NEED_MATH FALSE)
	message(STATUS "Using user-defined MATH_LIBS=${MATH_LIBS}")
else(MATH_LIBS)
	set(NEED_MATH TRUE)
	message(STATUS "MATH_LIBS=${MATH_LIBS}")
endif(MATH_LIBS)

if(NEED_MATH)
	if(IS_DIRECTORY $ENV{MKLROOT})
		set(MATH_LIBS "-Wl,--start-group  $ENV{MKLROOT}/lib/intel64/libmkl_intel_lp64.a $ENV{MKLROOT}/lib/intel64/libmkl_sequential.a $ENV{MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lgfortran -ldl")
	else(IS_DIRECTORY $ENV{MKLROOT})
		enable_language( Fortran )
		find_package(LAPACK REQUIRED)
		message(STATUS " LAPACK_LIBRARIE:  ${LAPACK_LIBRARIES}")
		message(STATUS " CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES:  ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES}")
		set(MATH_LIBS ${LAPACK_LIBRARIES} -ldl -lfl -lgomp ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES})
		message(STATUS " MATH_LIBS:  ${MATH_LIBS}")
	endif(IS_DIRECTORY $ENV{MKLROOT})
endif(NEED_MATH)

set(SHARED_DIR ${CMAKE_CURRENT_SOURCE_DIR}/ThirdPartyLibs)
set(AMPL_DIR ${SHARED_DIR}/ASL/src)
set(MA27_DIR ${SHARED_DIR}/MA27/src)
set(MA57_DIR ${SHARED_DIR}/MA57/src)
set(MA86_DIR ${SHARED_DIR}/MA86/src)
set(METIS_DIR ${SHARED_DIR}/METIS/src)
set(UMFPACK_DIR ${SHARED_DIR}/UMFPACK/src)

# include CoinUtils
set(COIN_DIR ${SHARED_DIR}/CBC/src)

if(NOT EXISTS "${COIN_DIR}/lib/pkgconfig/coinutils.pc")
	message(FATAL_ERROR "Please run \"./wgetCBC.sh\" in ${COIN_DIR}")
endif(NOT EXISTS "${COIN_DIR}/lib/pkgconfig/coinutils.pc")

if(NOT EXISTS "${MA27_DIR}/lib/libma27.a" AND NOT EXISTS "${MA57_DIR}/lib/libma57.a")
 	message(FATAL_ERROR "Could not find MA27 nor MA57. We need one HSL solver")
endif(NOT EXISTS "${MA27_DIR}/lib/libma27.a" AND NOT EXISTS "${MA57_DIR}/lib/libma57.a")

if(EXISTS "${MA57_DIR}/lib/libma57.a" AND NOT EXISTS "${METIS_DIR}/libmetis.a")
  message(FATAL_ERROR "Could not find METIS needed by MA57.")
endif(EXISTS "${MA57_DIR}/lib/libma57.a" AND NOT EXISTS "${METIS_DIR}/libmetis.a")

set(ENV{PKG_CONFIG_PATH} ${COIN_DIR}/lib/pkgconfig)
find_package(PkgConfig REQUIRED)
pkg_check_modules(COINUTILS REQUIRED coinutils)
pkg_check_modules(CLP REQUIRED clp)
pkg_check_modules(CBC REQUIRED cbc)
#find_library(CLP_LIB NAMES libClp.a Clp HINTS ${CLP_LIBRARY_DIRS})
find_library(CLP_LIB NAMES libClp.a Clp HINTS ${COIN_DIR}/lib)

set(COIN_LIBS "")
#message(STATUS "ZLIB_LIBRARY ${ZLIB_LIBRARY}")
#message(STATUS "COINUTILS_LIBRARIES ${COINUTILS_LIBRARIES}")
foreach(f ${COINUTILS_LIBRARIES})
	if(NOT ${f} MATCHES "m" AND NOT ${f} MATCHES "blas" AND NOT ${f} MATCHES "lapack") # skip libm
	    message(STATUS "looking after ${f}")
		find_library(X NAMES ${f} HINTS "${COIN_DIR}/lib" PATHS /usr/lib /usr/local/lib /usr/lib64/ /usr/lib/x86_64-linux-gnu/)
		message(STATUS "found ${X}")
		set(COIN_LIBS ${COIN_LIBS} ${X})
		unset(X CACHE)
	endif(NOT ${f} MATCHES "m" AND NOT ${f} MATCHES "blas" AND NOT ${f} MATCHES "lapack")
endforeach(f)
message(STATUS "COIN_LIBS: ${COIN_LIBS}")

message(STATUS "CBC_LIBRARY_DIRS: ${CBC_LIBRARY_DIRS}")
find_library(CBCSolver_LIB NAMES libCbcSolver.a CbcSolver HINTS ${CBC_LIBRARY_DIRS})
find_library(CBC_LIB NAMES libCbc.a Cbc HINTS ${CBC_LIBRARY_DIRS})
find_library(OSICLP_LIB NAMES libOsiClp.a OsiClp HINTS ${CBC_LIBRARY_DIRS})
find_library(OSI_LIB NAMES libOsi.a Osi HINTS ${CBC_LIBRARY_DIRS})
find_library(CGL_LIB NAMES libCgl.a Cgl HINTS ${CBC_LIBRARY_DIRS})

set(CBC_LIBS ${CBCSolver_LIB} ${CBC_LIB} ${CGL_LIB} ${OSICLP_LIB} ${OSI_LIB})
message(STATUS "CBC_LIBS: ${CBC_LIBS}")

find_library(AMPL_LIBRARY NAMES amplsolver.a asl PATHS ${AMPL_DIR} /usr/lib /usr/local/lib /usr/lib64)
find_path(AMPL_INCLUDE_DIR asl_pfgh.h HINTS ${AMPL_DIR} /usr/include/asl /usr/local/include/asl)
message(STATUS "AMPL_LIBRARY = ${AMPL_LIBRARY}")
message(STATUS "AMPL_INCLUDE_DIR = ${AMPL_INCLUDE_DIR}")

find_library(DL_LIBRARY NAMES dl)
find_library(FL_LIBRARY NAMES fl PATHS /usr/lib /usr/local/lib /usr/local/opt/flex)

find_library(MA27_LIBRARY ma27 PATHS ${MA27_DIR}/lib)
find_library(MA57_LIBRARY ma57 PATHS ${MA57_DIR}/lib)
find_library(MA86_LIBRARY hsl_ma86 PATHS ${MA86_DIR}/lib)
find_library(METIS_LIBRARY NAME libmetis.a HINTS ${METIS_DIR} PATHS ${METIS_DIR})

set(HAVE_AMPL TRUE)
if(${AMPL_LIBRARY} MATCHES "AMPL_LIBRARY-NOTFOUND")
	message(STATUS "amplsolver library not found. Won't build code that requires ASL.")
	set(HAVE_AMPL FALSE)
else(${AMPL_LIBRARY} MATCHES "AMPL_LIBRARY-NOTFOUND")
  include_directories(${AMPL_DIR}) 
endif(${AMPL_LIBRARY} MATCHES "AMPL_LIBRARY-NOTFOUND")


#see if we have METIS
#if we don't, also set HAVE_MA57 to false.
set(HAVE_METIS TRUE)
if(${METIS_LIBRARY} MATCHES "METIS_LIBRARY-NOTFOUND")
  message(STATUS "metis library not found. Won't build code that requires METIS (i.e. QpStoch and code that requires MA57).")
  set(HAVE_METIS FALSE)
endif(${METIS_LIBRARY} MATCHES "METIS_LIBRARY-NOTFOUND")
if(HAVE_METIS)
  include_directories(${SHARED_DIR}/METIS/src/Lib)
endif(HAVE_METIS)

set(HAVE_MA86 HAVE_METIS)
# Check if we have all libraries that Ma86 depends on exists.
if(${MA86_LIBRARY} MATCHES "MA86_LIBRARY-NOTFOUND")
  message(STATUS "MA86 library not found. Won't build code that requires MA86.")
  set(HAVE_MA86 FALSE)
endif(${MA86_LIBRARY} MATCHES "MA86_LIBRARY-NOTFOUND") 


# see if we have MA27 or MA57
set(HAVE_MA27 TRUE)
if(${MA27_LIBRARY} MATCHES "MA27_LIBRARY-NOTFOUND")
  message(STATUS "MA27 library not found. Won't build code that requires MA27.")
  set(HAVE_MA27 FALSE)
endif(${MA27_LIBRARY} MATCHES "MA27_LIBRARY-NOTFOUND")

set(HAVE_MA57 HAVE_METIS)
if(HAVE_MA57 AND ${MA57_LIBRARY} MATCHES "MA57_LIBRARY-NOTFOUND")
  message(STATUS "MA57 library not found. Won't build code that requires MA57.")
  set(HAVE_MA57 FALSE)
endif(HAVE_MA57 AND ${MA57_LIBRARY} MATCHES "MA57_LIBRARY-NOTFOUND")

# see if we have Umfpack
message(STATUS "${UMFPACK_DIR}")
message(STATUS "${UMFPACK_LIBRARY_PIPS}")
set(HAVE_UMFPACK TRUE)
set(HAVE_AMD FALSE)
if(${UMFPACK_LIBRARY_PIPS} MATCHES "UMFPACK_LIBRARY_PIPS-NOTFOUND")
  message(STATUS "UMFPACK library not found. Won't build code that requires UMFPACK.")
  set(HAVE_UMFPACK FALSE)
else (${UMFPACK_LIBRARY_PIPS} MATCHES "UMFPACK_LIBRARY_PIPS-NOTFOUND")
  message(STATUS "A library with UMFPACK found.")
  include_directories(${UMFPACK_DIR}/UMFPACK/Include)
  find_library(AMD_LIBRARY NAMES libamd.a PATHS ${UMFPACK_DIR}/AMD/Lib)
  message(STATUS "AMD library: ${AMD_LIBRARY}")
  if(${AMD_LIBRARY} MATCHES "AMD_LIBRARY-NOTFOUND")
    message(STATUS "AMD library not found. We require AMD to use UMFPACK.")
    set(HAVE_UMFPACK FALSE)
  else (${AMD_LIBRARY} MATCHES "AMD_LIBRARY-NOTFOUND")
    set(HAVE_AMD TRUE)
    include_directories(${UMFPACK_DIR}/AMD/Include)
    include_directories(${UMFPACK_DIR}/SuiteSparse_config)
    find_library(SUITESPARSE_LIBRARY NAMES libsuitesparseconfig.a PATHS ${UMFPACK_DIR}/SuiteSparse_config)
    if(${SUITESPARSE_LIBRARY} MATCHES "SUITESPARSE_LIBRARY-NOTFOUND")
      message(STATUS "SUITESPARSE_LIBRARY not found. We require SUITESPARSE_LIBRARY to use UMFPACK.")
      set(HAVE_UMFPACK FALSE)
    else (${SUITESPARSE_LIBRARY} MATCHES "SUITESPARSE_LIBRARY-NOTFOUND")
      set(UMF_ALLINONE_LIBRARY  ${UMFPACK_LIBRARY_PIPS} ${AMD_LIBRARY} ${SUITESPARSE_LIBRARY})
      message(STATUS "Find Umfpack: ${UMF_ALLINONE_LIBRARY}")
    endif(${SUITESPARSE_LIBRARY} MATCHES "SUITESPARSE_LIBRARY-NOTFOUND")
  endif(${AMD_LIBRARY} MATCHES "AMD_LIBRARY-NOTFOUND")
endif(${UMFPACK_LIBRARY_PIPS} MATCHES "UMFPACK_LIBRARY_PIPS-NOTFOUND")



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

set(CB_DIR "${SHARED_DIR}/ConicBundle/src")
set(CB_LIB "${CB_DIR}/lib/libcb.a")
if(EXISTS ${CB_LIB})
	set(HAVE_CB TRUE)
else(EXISTS ${CB_LIB})
	set(HAVE_CB FALSE)
	message(STATUS "Building without ConicBundle library")
endif(EXISTS ${CB_LIB})


set(PARDISO_DIR ${SHARED_DIR}/PARDISO/src)
set(HAVE_PARDISO FALSE)
if(${CMAKE_SYSTEM_NAME} MATCHES "BlueGeneP-static")
	set(HAVE_PARDISO TRUE)
	#find_library(PARDISO_LIBRARY NAMES libpardiso490-BGP_32.so PATHS ${PARDISO_DIR})
	#message(STATUS ${PARDISO_LIBRARY})
	set(PARDISO_LIBRARY "${PARDISO_DIR}/libpardiso490-BGP_32.so")
	set(PARDISO_LIBRARY32 "${PARDISO_DIR}/libpardiso490-BGP_32_32fact.so")
	#set(PARDISO_LIBRARY "libpardiso490-BGP_32.so")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    message(STATUS "BGP")
else(${CMAKE_SYSTEM_NAME} MATCHES "BlueGeneP-static")
    if(${CMAKE_SYSTEM_NAME} MATCHES "BlueGeneQ-static")
	    set(HAVE_PARDISO TRUE)
	    set(PARDISO_LIBRARY "${PARDISO_DIR}/libpardiso490-BGQ_64.so")
	    set(PARDISO_LIBRARY32 "${PARDISO_DIR}/libpardiso490-BGQ_64.so")
	    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS) ${OpenMP_EXE_LINKER_FLAGS}")
        message(STATUS "BGQ")
    else(${CMAKE_SYSTEM_NAME} MATCHES "BlueGeneQ-static")
        if(IS_XE6)
          	#find_library(PARDISO_LIBRARY NAMES "libpardiso495-GNU470-XE6.so"	  PATHS "/scratch/rosa/petra/shared_libs/")
	  		set(PARDISO_LIBRARY "/scratch/rosa/petra/shared_libs/libpardiso495-GNU430-XE6.so")
	  		message(STATUS "PARDISO library found, using ${PARDISO_LIBRARY}")
        else(IS_XE6)
	  		if(IS_XK7)
	    		set(PARDISO_LIBRARY "/scratch/todi/petra/shared_libs/libpardiso495-GNU430-XK7.so")
	    	else(IS_XK7)
	      		if(IS_XC30)
					set(PARDISO_LIBRARY "${PARDISO_LIBS_STATIC}")
	      		else(IS_XC30)
	        		set(PARDISO_LIBRARY ${PARDISO_DIR}/libpardiso490-GNU430-X86-64.so)
	      		endif(IS_XC30)
	    	endif(IS_XK7)
		endif(IS_XE6)

		if(EXISTS ${PARDISO_LIBRARY} OR IS_XC30)
			FIND_PACKAGE( OpenMP REQUIRED )
			set(HAVE_PARDISO TRUE)
			set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
			set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
		else(EXISTS ${PARDISO_LIBRARY} OR IS_XC30)
			set(PARDISO_LIBRARY "")
		endif(EXISTS ${PARDISO_LIBRARY} OR IS_XC30)
  endif(${CMAKE_SYSTEM_NAME} MATCHES "BlueGeneQ-static")
endif(${CMAKE_SYSTEM_NAME} MATCHES "BlueGeneP-static")

set(MA86_DIR "${SHARED_DIR}/MA86/src")
set(MA86_LIB "${MA86_DIR}/lib/libhsl_ma86.a")
set(MA86_INC "${MA86_DIR}/include")

set(SCIP_DIR "${SHARED_DIR}/SCIP/src")
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


find_package(Boost COMPONENTS system filesystem thread log log_setup)
include_directories(${Boost_INCLUDE_DIR})
if (${Boost_FOUND}$ MATCHES "FALSE")
  message("! A Boost version containing BOOST.LOG has NOT been found.")
  message("! Please specify the location of the Boost library containing 'log' package.")
  message("! For this create a 'Toolchain.cmake' file or edit the one from the root ")
  message("! directory and use 'set(BOOST_ROOT \"/path/to/boost\")'.")
  message("! Then run  'cmake -DCMAKE_TOOLCHAIN_FILE=../Toolchain.cmake ..' in the build directory.")
  message(FATAL_ERROR "Configuration aborted due to missing BOOST libraries!")
endif(${Boost_FOUND}$ MATCHES "FALSE")


###############################################################################

if(HAVE_MA27) 
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DWITH_MA27") 
else(HAVE_MA27) 
  set(MA27_LIBRARY "") 
endif(HAVE_MA27)

if(HAVE_MA57) 
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DWITH_MA57") 
else(HAVE_MA57) 
  set(MA57_LIBRARY "") 
endif(HAVE_MA57)

if(HAVE_MA86) 
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DWITH_MA86") 
else(HAVE_MA57) 
  set(MA86_LIBRARY "") 
endif(HAVE_MA86) 

if(HAVE_PARDISO) 
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DWITH_PARDISO") 
else(HAVE_PARDISO) 
  set(PARDISO_LIBRARY "") 
endif(HAVE_PARDISO)

if(HAVE_UMFPACK) 
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DWITH_UMFPACK") 
else(HAVE_UMFPACK) 
  set(UMF_ALLINONE_LIBRARY "") 
endif(HAVE_UMFPACK)

###############################################################################


#set(CMAKE_BUILD_TYPE DEBUG)
#set(CMAKE_BUILD_TYPE RELWITHDEBINFO)
#set(CMAKE_BUILD_TYPE RELEASE)
#message("FLAGS: ${CMAKE_CXX_FLAGS_RELEASE}")

# must set TAU_MAKEFILE, TAU_OPTIONS environment variables
#set(CMAKE_CXX_COMPILER "tau_cxx.sh")

#IF(CMAKE_COMPILER_IS_GNUCXX)
#	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
#ENDIF(CMAKE_COMPILER_IS_GNUCXX)

include_directories(Input)
include_directories(Input/OPF_Matpower)

include_directories(SolverInterface)

#message(STATUS "${COINUTILS_INCLODE_DIRS}")
include_directories(${COINUTILS_INCLUDE_DIRS})

add_subdirectory(PIPS-IPM)
add_subdirectory(Input)
add_subdirectory(SolverInterface)
add_subdirectory(PIPS-S)
add_subdirectory(Lagrange)
add_subdirectory(Utilities)
add_subdirectory(UCTests)
#add_subdirectory(UCRolling)
add_subdirectory(BasisBootstrap)
add_subdirectory(PIPS-NLP)

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


#if(HAVE_GUROBI)
#	add_subdirectory(src/GUROBI)
#endif(HAVE_GUROBI)
message(STATUS "CMAKE_CXX_FLAGS  ${CMAKE_CXX_FLAGS}")

