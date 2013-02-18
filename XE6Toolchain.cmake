set(CMAKE_SYSTEM_NAME Linux)

# in build directory:
# cmake -DCMAKE_TOOLCHAIN_FILE=../XE6Toolchain.cmake ..
# for Cbc:
# CXX="CC" ./configure --enable-static=yes --enable-shared=no --host=x86_64-unknown-linux-gnu
set(CMAKE_C_COMPILER cc)
set(CMAKE_CXX_COMPILER CC)
set(CMAKE_Fortran_COMPILER ftn)

set(IS_XE6 TRUE)

#set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")
#set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O4")
#set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g -O2 -dynamic") # -qfloat=maf 

set(ROSA_LAPACK_BLAS "/opt/cray/libsci/12.0.00/GNU/47/interlagos/lib/libsci_gnu_mp.a")

set(MATH_LIBS "${ROSA_LAPACK_BLAS}")

set(BOOST_ROOT "/users/petra/boost_1_47_0")
#set(PARDISO_LIBRARY32 "/scratch/rosa/petra/shared_libs/libpardiso491-GNU430-X86-64.so")

