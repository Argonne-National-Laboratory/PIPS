#add_definitions(-DTIMING -DSTOCH_TESTING) # timing output
include_directories(Core/Abstract Core/Vector Core/Utilities Core/QpSolvers Core/QpPreprocess Core/QpGen
  Core/SparseLinearAlgebra Core/DenseLinearAlgebra Core/Readers
  Core/LinearSolvers/Ma27Solver Core/LinearSolvers/Ma57Solver
  Core/LinearSolvers/Ma86Solver Core/LinearSolvers/PardisoSolver Core/LinearSolvers/BiCGStabSolver Core/LinearSolvers/MumpsSolver)
include_directories(Core/StochLinearAlgebra Core/QpStoch)
add_subdirectory(Core)

set(EXECUTABLE_OUTPUT_PATH "${PROJECT_BINARY_DIR}")

if (HAVE_MA57 AND HAVE_METIS)
    add_executable(qpgen-sparse-ma57-mehrotra.exe Core/QpGen/QpGenSparseMa57MehrotraDriver.C)
    target_link_libraries(qpgen-sparse-ma57-mehrotra.exe
      ooqpgensparse ooqpsparse ooqpbase ooqpmehrotra ooqpdense
      ${MA57_LIBRARY} ${METIS_LIBRARY} ${MATH_LIBS})

    # add_executable(pipsipm-unitcommit.exe Drivers/unitCommitment.cpp)
    # target_link_libraries(pipsipm-unitcommit.exe
    #   ooqpstoch ooqpstochla ooqpmehrotrastoch
    #   ooqpgensparse ooqpbase ooqpsparse ooqpdense
    #   ${MA57_LIBRARY} ${METIS_LIBRARY} ${MATH_LIBS})

  if (HAVE_CLANG)
  else()
    add_library(pipsipm-shared SHARED Drivers/pipsipm_C_callbacks.cpp)
    target_link_libraries(pipsipm-shared
      ${WHOLE_ARCHIVE}
      stochInput ${COIN_LIBS}
      ooqpstoch ooqpstochla ooqpmehrotrastoch
      ooqpgensparse ooqpbase ooqpsparse ooqpdense
      ${MA57_LIBRARY} ${METIS_LIBRARY}
      ${NO_WHOLE_ARCHIVE}
      ${MATH_LIBS})
  endif (HAVE_CLANG)

endif(HAVE_MA57 AND HAVE_METIS)

#missing linkage mc68 using metis.
#if (HAVE_MA86 AND HAVE_METIS)
#  add_executable(qpgen-sparse-ma86-mehrotra.exe Core/QpGen/QpGenSparseMa86MehrotraDriver.C)
#  target_link_libraries(qpgen-sparse-ma86-mehrotra.exe
#    ooqpgensparse ooqpsparse ooqpbase ooqpmehrotra ooqpdense
#    ${MA86_LIB} ${METIS_LIBRARY} ${MATH_LIBS})
#endif(HAVE_MA86 AND HAVE_METIS)

#for PIPS_IPM there is the option of linking MKL_PARDISO instead of PARDISO
if (HAVE_PARDISO OR (HAVE_MKL_PARDISO AND WITH_MKL_PARDISO))
  if(HAVE_MKL_PARDISO AND WITH_MKL_PARDISO)
    set(PARDISO_LIBRARY_CHOSEN "${MKL_PARDISO_LIBS}")
  else(HAVE_MKL_PARDISO AND WITH_MKL_PARDISO)
    set(PARDISO_LIBRARY_CHOSEN "${PARDISO_LIBRARY}")
  endif(HAVE_MKL_PARDISO AND WITH_MKL_PARDISO)

  if(HAVE_MA27)
    add_executable(qpgen-sparse-mehrotra.exe Core/QpGen/QpGenSparseMehrotraDriver.C)
    target_link_libraries(qpgen-sparse-mehrotra.exe
      ooqpgensparse ooqpsparse ooqpbase ooqpmehrotra
    ${MA27_LIBRARY} ${PARDISO_LIBRARY_CHOSEN} ${MATH_LIBS})
  endif(HAVE_MA27)

  if(HAVE_MA57 AND HAVE_METIS)
    add_executable(qpgen-sparse-ma57-gondzio.exe Core/QpGen/QpGenSparseGondzioDriver.C)
    target_link_libraries(qpgen-sparse-ma57-gondzio.exe
      ooqpgensparse ooqpsparse ooqpbase ooqpgondzio ooqpdense
      ${MA57_LIBRARY} ${METIS_LIBRARY} ${PARDISO_LIBRARY_CHOSEN} ${MATH_LIBS})

  if( 0 )
    add_executable(pipsipmFromRaw_comm2_schur Drivers/pipsipmFromRaw_comm2_schur.cpp)
    target_link_libraries(pipsipmFromRaw_comm2_schur
      stochInput ${COIN_LIBS}
      ooqpstoch ooqpstochla ooqpmehrotrastoch
      ooqpgensparse ooqpbase ooqpsparse ooqpdense
      ${MA57_LIBRARY} ${METIS_LIBRARY} ${PARDISO_LIBRARY_CHOSEN} ${MATH_LIBS})

    add_executable(pipsipmBatchFromRaw_schur Drivers/pipsipmBatchFromRaw_schur.cpp)
    target_link_libraries(pipsipmBatchFromRaw_schur
      stochInput ${COIN_LIBS}
      ooqpstoch ooqpstochla ooqpmehrotrastoch
      ooqpgensparse ooqpbase ooqpsparse ooqpdense
      ${MA57_LIBRARY} ${METIS_LIBRARY} ${PARDISO_LIBRARY_CHOSEN} ${MATH_LIBS})

    add_executable(pipsipmBatchFromRaw Drivers/pipsipmBatchFromRaw.cpp)
    target_link_libraries(pipsipmBatchFromRaw
      stochInput ${COIN_LIBS}
      ooqpstoch ooqpstochla ooqpmehrotrastoch
      ooqpgensparse ooqpbase ooqpsparse ooqpdense
      ${MA57_LIBRARY} ${METIS_LIBRARY} ${PARDISO_LIBRARY_CHOSEN} ${MATH_LIBS})
  endif( 0 )   
   
    # add_executable(pipsipmFromRaw_schur32 Drivers/pipsipmFromRaw_schur32.cpp)
    # target_link_libraries(pipsipmFromRaw_schur32
    #   stochInput ${COIN_LIBS}
    #   ooqpstoch ooqpstochla ooqpmehrotrastoch
    #   ooqpgensparse ooqpbase ooqpsparse ooqpdense
    #   ${MA57_LIBRARY} ${METIS_LIBRARY} ${PARDISO_LIBRARY32} ${MATH_LIBS})

    if(HAVE_MA27)
#    add_executable(pipsipmFromRaw Drivers/pipsipmFromRaw.cpp)
#    target_link_libraries(pipsipmFromRaw
#      stochInput ${COIN_LIBS}
#      ooqpstoch ooqpstochla ooqpmehrotrastoch
#      ooqpgensparse ooqpbase ooqpsparse ooqpdense
#      ${MA27_LIBRARY} ${MA57_LIBRARY} ${METIS_LIBRARY} ${PARDISO_LIBRARY_CHOSEN} ${MATH_LIBS})

    add_executable(ooqpFromRaw Drivers/ooqpFromRaw.cpp)
    target_link_libraries(ooqpFromRaw
      stochInput ${COIN_LIBS}
      ooqpgensparse ooqpbase ooqpmehrotra ooqpsparse ooqpdense
      ${MA57_LIBRARY} ${METIS_LIBRARY} ${MA27_LIBRARY} ${MATH_LIBS})

    add_executable(pipsipmFromRaw_schur Drivers/pipsipmFromRaw_schur.cpp)
    target_link_libraries(pipsipmFromRaw_schur
      stochInput ${COIN_LIBS}
      ooqpstoch ooqpstochla ooqpmehrotrastoch
      ooqpgensparse ooqpbase ooqpsparse ooqpdense
      ${MA27_LIBRARY} ${MA57_LIBRARY} ${METIS_LIBRARY} ${PARDISO_LIBRARY_CHOSEN} ${MATH_LIBS})

      add_executable(meanSolveAndRecourseEval Drivers/meanSolveAndRecourseEval.cpp)
      target_link_libraries(meanSolveAndRecourseEval 
        stochInput ${COIN_LIBS}
	      ooqpstoch ooqpstochla ooqpmehrotrastoch
	     ooqpgensparse ooqpbase ooqpmehrotra ooqpsparse ooqpdense
	     ${MA27_LIBRARY} ${MA57_LIBRARY} ${METIS_LIBRARY} ${PARDISO_LIBRARY_CHOSEN}
	     ${MATH_LIBS})
    endif(HAVE_MA27)

#    if(HAVE_AMPL)
#      add_executable(pipsipmPyomo_schur Drivers/pipsipmPyomo_schur.cpp)
#      target_link_libraries(pipsipmPyomo_schur
#	pyomoStochInput stochInput ${COIN_LIBS}
#	ooqpstoch ooqpstochla ooqpmehrotrastoch
#	ooqpgensparse ooqpbase ooqpsparse ooqpdense
#	${AMPL_LIBRARY} ${MA57_LIBRARY} ${METIS_LIBRARY} ${PARDISO_LIBRARY_CHOSEN} ${MATH_LIBS})
#    endif(HAVE_AMPL)
  endif(HAVE_MA57 AND HAVE_METIS)
endif(HAVE_PARDISO OR (HAVE_MKL_PARDISO AND WITH_MKL_PARDISO))

#MA27 or MA57 and METIS are prerequisites, this should be fine without additional checks.
#It still works if only one of MA27 or MA57 is available.
add_executable(pipsipmFromRaw Drivers/pipsipmFromRaw.cpp)
target_link_libraries(pipsipmFromRaw
    stochInput ${COIN_LIBS}
    ooqpstoch ooqpstochla ooqpmehrotrastoch
    ooqpgensparse ooqpbase ooqpsparse ooqpdense
    ${MA27_LIBRARY} ${MA57_LIBRARY} ${METIS_LIBRARY} ${MATH_LIBS} ${PARDISO_LIBRARY_CHOSEN})
    

add_executable(pipsipmCallbackExample Drivers/callbackExample.cpp)
target_link_libraries(pipsipmCallbackExample
    stochInput ${COIN_LIBS}
    ooqpstoch ooqpstochla ooqpgondziostoch ooqpqpscaler ooqpequistochscaler
    ooqpgensparse ooqpbase ooqpsparse ooqpdense ooqpqppresolver ooqpqpstochpresolver
    ${MA27_LIBRARY} ${MA57_LIBRARY} ${MUMPS_LIBRARY} ${PARMETIS_LIBRARY} ${METIS_LIBRARY} ${PARDISO_LIBRARY_CHOSEN} ${MATH_LIBS} ${COIN_LIBS})
if (BUILD_GDX_SOURCE)
   add_definitions(-D_P3_DllExport= -DSTDCALL= -DGDX_CALLCONV=unwanted -DLEG)
   add_library(statgdx
     Drivers/statgdx/clibtypes.cpp
     Drivers/statgdx/datastorage.cpp
     Drivers/statgdx/gdlaudit.cpp
     Drivers/statgdx/gdxdclib.cpp
     Drivers/statgdx/gdxstatic.cpp
     Drivers/statgdx/gmsdata.cpp
     Drivers/statgdx/gmsgen.cpp
     Drivers/statgdx/gmsglob.cpp
     Drivers/statgdx/gmsglobx.cpp
     Drivers/statgdx/gmsheapnew.cpp
     Drivers/statgdx/gmslibname.cpp
     Drivers/statgdx/gmsobj.cpp
     Drivers/statgdx/gmsspecs.cpp
     Drivers/statgdx/gmsstrm.cpp
     Drivers/statgdx/gxdefs.cpp
     Drivers/statgdx/gxfile.cpp
     Drivers/statgdx/idglobal_p3.cpp
     Drivers/statgdx/math_p3.cpp
     Drivers/statgdx/p3Custom2.cpp
     Drivers/statgdx/p3ieeefp.cpp
     Drivers/statgdx/p3io.cpp
     Drivers/statgdx/p3library.cpp
     Drivers/statgdx/p3platform.cpp
     Drivers/statgdx/p3private.cpp
     Drivers/statgdx/p3process.cpp
     Drivers/statgdx/p3threads.cpp
     Drivers/statgdx/p3utils.cpp
     Drivers/statgdx/paldoorg.cpp
     Drivers/statgdx/pchutil.cpp
     Drivers/statgdx/runner.cpp
     Drivers/statgdx/strhash.cpp
     Drivers/statgdx/strutilx.cpp
     Drivers/statgdx/system_p3.cpp
     Drivers/statgdx/sysutils_p3.cpp
     Drivers/statgdx/utils2.cpp
     Drivers/statgdx/xcompress.cpp
     Drivers/statgdx/gclgms.c )
     
   link_directories(PIPS-IPM)
   add_executable(gmspips Drivers/gmspips/gmspips.cpp Drivers/gmspips/gmspipsio.c)
   add_definitions(-DGMS_PIPS -DGMS_MPI -DLINKCONSTR -DGDXSOURCE -fno-strict-aliasing -fwrapv ) # -DGMS_LOG
   target_include_directories(gmspips PUBLIC Drivers/statgdx Drivers/gmspips/)
   target_link_libraries(gmspips
       dl 
       statgdx
       stochInput ${COIN_LIBS}
       ooqpstoch ooqpstochla ooqpgondziostoch ooqpequistochscaler
       ooqpgensparse ooqpbase ooqpsparse ooqpdense ooqpqppresolver ooqpqpstochpresolver
       ${MA27_LIBRARY} ${MA57_LIBRARY} ${METIS_LIBRARY} ${MUMPS_LIBRARY}  ${PARDISO_LIBRARY_CHOSEN} ${MATH_LIBS} ${COIN_LIBS})
   
   add_executable(gmschk Drivers/gmspips/gmspipschk.cpp Drivers/gmspips/gmspipsio.c)
   target_include_directories(gmschk PUBLIC Drivers/statgdx Drivers/gmspips/)
   target_link_libraries(gmschk dl statgdx)
else()
   set(GAMSDIR Drivers/gams)
   add_executable(gmspips Drivers/gmspips/gmspips.cpp Drivers/gmspips/gmspipsio.c ${GAMSDIR}/apifiles/C/api/gmomcc.c ${GAMSDIR}/apifiles/C/api/gevmcc.c ${GAMSDIR}/apifiles/C/api/gdxcc.c)
   add_definitions(-DGMS_PIPS -DGMS_MPI -DLINKCONSTR) # -DGMS_LOG
   target_include_directories(gmspips PUBLIC ${GAMSDIR}/apifiles/C/api/ Drivers/gmspips/)
   target_link_libraries(gmspips
       dl
       stochInput ${COIN_LIBS}
       ooqpstoch ooqpstochla ooqpgondziostoch ooqpequistochscaler
       ooqpgensparse ooqpbase ooqpsparse ooqpdense ooqpqppresolver ooqpqpstochpresolver
       ${MA27_LIBRARY} ${MA57_LIBRARY} ${METIS_LIBRARY}  ${PARDISO_LIBRARY_CHOSEN} ${MATH_LIBS} ${COIN_LIBS})
   
   add_executable(gmschk Drivers/gmspips/gmspipschk.cpp Drivers/gmspips/gmspipsio.c ${GAMSDIR}/apifiles/C/api/gmomcc.c ${GAMSDIR}/apifiles/C/api/gevmcc.c ${GAMSDIR}/apifiles/C/api/gdxcc.c)
   target_include_directories(gmschk PUBLIC ${GAMSDIR}/apifiles/C/api/ Drivers/gmspips/)
   target_link_libraries(gmschk dl)
endif(BUILD_GDX_SOURCE)
