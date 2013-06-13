set(CMAKE_SYSTEM_NAME Linux)

# in build directory:
# cmake -DCMAKE_TOOLCHAIN_FILE=../XC30Toolchain.cmake ..
# for Cbc:
# CXX="CC" ./configure --enable-static=yes --enable-shared=no --host=x86_64-unknown-linux-gnu
set(CMAKE_C_COMPILER cc)
set(CMAKE_CXX_COMPILER CC)
set(CMAKE_Fortran_COMPILER ftn)

set(IS_XC30 TRUE)

#set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")
#set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O4")
#set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g -O2") # -qfloat=maf 
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2") # -qfloat=maf 


#set(ROSA_LAPACK_BLAS "/opt/cray/libsci/12.0.00/GNU/47/sandybridge/lib/libsci_gnu_mp.a")
set(ROSA_LAPACK_BLAS "-L/opt/cray/libsci/12.0.00/GNU/47/sandybridge/lib/ -lsci_gnu_mp  -lsci_gnu_mp -lsci_gnu_mp")

#set(ROSA_LAPACK_BLAS "-L/opt/intel/13.0.1.117/mkl/lib/intel64/  -lifcore -lifport -lsvml -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lmkl_intel_lp64 -lmkl_sequential -lmkl_core")
set(MATH_LIBS "${ROSA_LAPACK_BLAS}")

#set(PARDISO_LIBS_STATIC "-L/scratch/daint/petra/shared_libs/ -lmetis41_pardiso -lmetis41-P_pardiso -lpardiso -lpils_pardiso -lmetis41_pardiso -lmetis41-P_pardiso -lpardiso -lpils_pardiso")
set(PARDISO_LIBS_STATIC "-L/scratch/daint/petra/shared_libs/ -lmetis -lpardiso -lmetis -lpardiso")
#set(MATH_LIBS "${MATH_LIBS} ${PARDISO_LIBS_STATIC}")

set(BOOST_ROOT "/users/petra/boost_1_47_0")

