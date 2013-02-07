set(CMAKE_SYSTEM_NAME BlueGeneQ-static)

# in build directory:
# cmake -DCMAKE_TOOLCHAIN_FILE=../BGQToolchain.cmake ..
# for Cbc:
# CXX="bgxlc++_r" ./configure --enable-static=yes --enable-shared=no --host=powerpc64-bgq-linux
set(CMAKE_C_COMPILER mpixlc_r)
set(CMAKE_CXX_COMPILER mpixlcxx_r)
set(CMAKE_Fortran_COMPILER bgxlf_r)

set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")
#set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O4")
#set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g -O2 -qfloat=maf")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g -O2 -qfloat=nofold")

set(BGP_LAPACK "/soft/libraries/alcf/current/xl/LAPACK/lib/liblapack.a")
set(PURE_ESSL "/soft/libraries/essl/5.1.1-0/lib64/libesslbg.a")

# set this version to match the compiler used
# below corresponds to @ibm-compilers-apr2012 in softenv
set(IBMCMP_BASE "/soft/compilers/ibmcmp-nov2012")
set(XLF_BASE "${IBMCMP_BASE}/xlf/bg/14.1/bglib64")
set(XLF_LIBS "-L${XLF_BASE} -lxlfmath -lxlf90_r")

set(XLSMP_BASE "${IBMCMP_BASE}/xlsmp/bg/3.1/bglib64")
set(XLOMP_SER "-L${XLSMP_BASE} -lxlomp_ser")
set(XLSMP "-L${XLSMP_BASE} -lxlsmp")

set(MATH_LIBS "${BGP_LAPACK};${PURE_ESSL};${XLF_LIBS};${XLSMP}")

set(BOOST_ROOT "/home/cpetra/boost_1_47_0")
set(OpenMP_CXX_FLAGS "-qsmp=omp -qnoipa")
