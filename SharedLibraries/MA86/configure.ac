AC_INIT(hsl_ma86, 1.3.0, hsl@stfc.ac.uk)
AC_CONFIG_SRCDIR([src/hsl_ma86s.f90])
AM_INIT_AUTOMAKE(hsl_ma86, 1.3.0)

# Check for required programs
AC_PROG_F77
AC_PROG_FC
AC_PROG_CC
AC_PROG_RANLIB

# Establish some required flags
AC_FC_LIBRARY_LDFLAGS # Linking C to Fortran

# Establish correct module search path flag
AX_F90_MODULE_FLAG

# Check for external libraries (e.g. BLAS, LAPACK, MeTiS)
ACX_BLAS(,[AC_MSG_WARN([No BLAS library found. make check will not work.])])
HSL_METIS(,
   [AC_MSG_WARN([No MeTiS library found. make check will not work.])]
   )

# Output data
AC_CONFIG_FILES(Makefile src/Makefile example/Fortran/Makefile example/C/Makefile)
AC_OUTPUT
