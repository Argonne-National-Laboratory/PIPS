set(CMAKE_SYSTEM_NAME BlueGeneP-static)

# in build directory:
# cmake -DCMAKE_TOOLCHAIN_FILE=../BGPToolchain.cmake ..
# for Cbc:
# CXX="bgxlc++_r" ./configure --enable-static=yes --enable-shared=no --host=powerpc-bgp-linux
set(CMAKE_C_COMPILER mpixlc_r)
set(CMAKE_CXX_COMPILER mpixlcxx_r)
set(CMAKE_Fortran_COMPILER bgxlf_r)

set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")
#set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O4")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g -O2 -qfloat=maf")

#set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} /home/mlubin/flex-2.5.35)

set(BGP_LAPACK "/soft/apps/LAPACK/liblapack_bgp.a")
set(PURE_ESSL "/soft/apps/ESSL-4.4.1-1/lib/libesslbg.a")

# set this version to match the compiler used
# below corresponds to @ibm-compilers-apr2012 in softenv
set(IBMCMP_BASE "/soft/apps/ibmcmp-apr2012")
set(XLF_BASE "${IBMCMP_BASE}/xlf/bg/11.1/bglib")
set(XLF_LIBS "-L${XLF_BASE} -lxlfmath -lxlf90_r")

set(XLSMP_BASE "${IBMCMP_BASE}/xlsmp/bg/1.7/bglib")
set(XLOMP_SER "-L${XLSMP_BASE} -lxlomp_ser")
set(XLSMP "-L${XLSMP_BASE} -lxlsmp")

set(MATH_LIBS "${BGP_LAPACK};${PURE_ESSL};${XLF_LIBS};${XLSMP}")

set(BOOST_ROOT "/home/mlubin/boost_1_47_0")
set(OpenMP_CXX_FLAGS "-qsmp=omp")
