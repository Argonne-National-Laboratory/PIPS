# Note: loosely based on immxdb_lib_metis.m4 and acx_blas.m4
AC_DEFUN([HSL_METIS], [
AC_MSG_CHECKING(for METIS library)
AC_REQUIRE([AC_PROG_CC])
#
# User hints...
#
AC_ARG_WITH([metis],
   [AC_HELP_STRING([--with-metis=<lib>],
   [user METIS library <lib>])])
case $with_metis in
   -* | */* | *.a | *.so | *.so.* | *.o) METIS_LIBS="$with_metis" ;;
    *) METIS_LIBS="-L$with_metis -lmetis" ;;
esac

# Get fortran linker names for function of interest
AC_F77_FUNC(metis_nodend)

hsl_metis_ok=no

# Check supplied location
if test $hsl_metis_ok = no; then
if test "x$METIS_LIBS" != x; then
   save_LIBS="$LIBS"; LIBS="$METIS_LIBS $LIBS -lm"
   AC_MSG_CHECKING([for $metis_nodend in $METIS_LIBS])
   AC_TRY_LINK_FUNC($metis_nodend, [hsl_metis_ok=yes], [METIS_LIBS=""])
   AC_MSG_RESULT($hsl_metis_ok)
   LIBS="$save_LIBS"
fi
fi

AC_SUBST(METIS_LIBS)

# Try just -lmetis
if test $hsl_metis_ok = no; then
   AC_CHECK_LIB(metis, $metis_nodend, [hsl_metis_ok=yes; METIS_LIBS="-lmetis"], [], [-lm])
fi

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$hsl_metis_ok" = xyes; then
        ifelse([$1],,AC_DEFINE(HAVE_METIS,1,[Define if you have a MeTiS library.]),[$1])
        :
else
        hsl_metis_ok
        $2
fi
])dnl HSL_METIS
