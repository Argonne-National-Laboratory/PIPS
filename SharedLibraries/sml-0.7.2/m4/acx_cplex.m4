# SYNOPSIS
#
#   ACX_CPLEX([ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND])
#
# DESCRIPTION
#
#   This macro searches for the Cplex library in the location specified by
#   the user by using the --with-cplex option to configure. If no location
#   was specified then the macro fails. Upon sucessful completion the
#   variables CPLEX_INCLUDE and CPLEX_LIBS are set.
#
#   ACTION-IF-FOUND is a list of shell commands to run if the Cplex library is
#   found, and ACTION-IF-NOT-FOUND is a list of commands to run it if it is
#   not found. If ACTION-IF-FOUND is not specified, the default action will
#   define HAVE_CPLEX. If ACTION-IF-NOT-FOUND is not specified then an error
#   will be generated halting configure.
#
# LAST MODIFICATION
#
#   2009-12-10
#
# COPYLEFT
#
#   Copyright (c) 2009 Marco Colombo <m.colombo@ed.ac.uk>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved.

AC_DEFUN([ACX_CPLEX], [

AC_ARG_WITH([cplex],
            [AC_HELP_STRING([--with-cplex=<dir>],
                            [path to the Cplex installation])])

case $with_cplex in
    no | "") ;;
    yes) AC_MSG_WARN([Cplex installation path missing in --with-cplex]) ;;
    *) CPLEX_DIR=$with_cplex ;;
esac

#
# Locate the Cplex headers
#
if test -n "${CPLEX_DIR}" ; then

    INCLUDE_DIR="${CPLEX_DIR}/include"
    AC_MSG_CHECKING(CPLEX include directory)
    if test [! -d "${INCLUDE_DIR}"] ; then
        AC_MSG_RESULT(failed)
        AC_MSG_WARN([directory ${INCLUDE_DIR} does not exist])
    else
        AC_MSG_RESULT([${INCLUDE_DIR}])
        CPLEX_INCLUDE="-I${INCLUDE_DIR}"
    fi
fi

#
# Locate the Cplex libraries
#
if test -n "$CPLEX_INCLUDE" ; then
    LIB_DIR="${CPLEX_DIR}/lib"
    SYSTEM=`basename $(ls -d ${LIB_DIR}/*/ | head -n1)`
    LIB_DIR="${LIB_DIR}/${SYSTEM}"
    LIBFORMAT=`basename $(ls -d ${LIB_DIR}/*/ | head -n1)`
    if test [-z "$LIBFORMAT" -o ! -d "${LIB_DIR}/${LIBFORMAT}"] ; then
        AC_MSG_WARN([Cannot find the location of the CPLEX library])
     else
        CPLEX_LIB_DIR="-L$CPLEX_DIR/lib/$SYSTEM/$LIBFORMAT"
    fi
fi

#
# Check headers and libraries
#
if test -n "${CPLEX_LIB_DIR}" ; then
    old_CFLAGS=$CFLAGS
    old_LDFLAGS=$LDFLAGS
    CFLAGS=$CPLEX_INCLUDE
    LDFLAGS=$CPLEX_LIB_DIR

    AC_CHECK_HEADER(ilcplex/cplex.h, [cplex_header_ok=yes], , [/* quiet */])
    AC_CHECK_LIB(cplex, CPXversion, [cplex_library_ok=yes], , -lm -lpthread)

    CFLAGS=$old_CFLAGS
    LDFLAGS=$old_LDFLAGS
fi

#
# Execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND
#
if test [ x"$cplex_header_ok" = xyes -a x"$cplex_library_ok" = xyes ] ; then
    CPLEX_LIBS="$CPLEX_LIB_DIR -lilocplex -lcplex -lpthread"
    AC_SUBST(CPLEX_INCLUDE, [$CPLEX_INCLUDE])
    AC_SUBST(CPLEX_LIBS, [$CPLEX_LIBS])
    ifelse([$1], , AC_DEFINE(HAVE_CPLEX, 1,
                             [Define if you have the CPLEX library.]), [$1])
else
    $2
    :
fi

]) # ACX_CPLEX
