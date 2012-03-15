# ===========================================================================
#                http://autoconf-archive.cryp.to/acx_ampl.html
# ===========================================================================
#
# SYNOPSIS
#
#   ACX_AMPL([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
#
# DESCRIPTION
#
#   This macro looks for a library that implements the amplsolver reader
#   routines for interacting with AMPL generated .nl files. On success, it
#   sets the AMPL_LIBS output variable to hold the requisite library linkages,
#   and the AMPL_INCLUDE output variable to hold the neccessary include flags.
#
#   ACTION-IF-FOUND is a list of shell commands to run if an AMPL library is
#   found, and ACTION-IF-NOT-FOUND is a list of commands to run it if it is
#   not found. If ACTION-IF-FOUND is not specified, the default action will
#   define HAVE_AMPL.
#
#   This macro requires autoconf 2.50 or later.
#
# LAST MODIFICATION
#
#   2009-12-10
#
# COPYLEFT
#
#   Copyright (c) 2009 Jonathan D. Hogg <J.Hogg@ed.ac.uk>
#   based on ACX_AMPL:
#   Copyright (c) 2008 Steven G. Johnson <stevenj@alum.mit.edu>
#
#   This program is free software: you can redistribute it and/or modify it
#   under the terms of the GNU General Public License as published by the
#   Free Software Foundation, either version 3 of the License, or (at your
#   option) any later version.
#
#   This program is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
#   Public License for more details.
#
#   You should have received a copy of the GNU General Public License along
#   with this program. If not, see <http://www.gnu.org/licenses/>.
#
#   As a special exception, the respective Autoconf Macro's copyright owner
#   gives unlimited permission to copy, distribute and modify the configure
#   scripts that are the output of Autoconf when processing the Macro. You
#   need not follow the terms of the GNU General Public License when using
#   or distributing such scripts, even though portions of the text of the
#   Macro appear in them. The GNU General Public License (GPL) does govern
#   all other use of the material that constitutes the Autoconf Macro.
#
#   This special exception to the GPL applies to versions of the Autoconf
#   Macro released by the Autoconf Macro Archive. When you make and
#   distribute a modified version of the Autoconf Macro, you may extend this
#   special exception to the GPL to apply to your modified version as well.

AC_DEFUN([ACX_AMPL], [
AC_PREREQ(2.50)
acx_ampl_ok=no

AC_ARG_WITH(ampl,
	[AC_HELP_STRING([--with-ampl=<dir>], [path to the AmplSolver library])])
case $with_ampl in
	yes | "") ;;
	no) acx_ampl_ok=disable ;;
	*) AMPL_LIBS="$with_ampl/amplsolver.a"; AMPL_INCLUDE="-I$with_ampl" ;;
esac

# Add -ldl to LIBS
AC_SEARCH_LIBS(dlopen, dl)

# First, check AMPL_LIBS environment variable
if test $acx_ampl_ok = no; then
if test "x$AMPL_LIBS" != x; then
	save_LIBS="$LIBS"; LIBS="$AMPL_LIBS $LIBS"
	save_CFLAGS="$CFLAGS"
	CFLAGS="$AMPL_INCLUDE"
	AC_CHECK_HEADER(asl_pfgh.h, , [AMPL_INCLUDE=""], [/* quiet */])
	AC_MSG_CHECKING([for ASL_alloc in $AMPL_LIBS])
	AC_TRY_LINK_FUNC(ASL_alloc, [acx_ampl_ok=yes], [AMPL_LIBS=""])
	AC_MSG_RESULT($acx_ampl_ok)
	LIBS="$save_LIBS"
	CFLAGS="$save_CFLAGS"
fi
fi

# Generic AMPL library?
if test $acx_ampl_ok = no; then
	AC_CHECK_LIB(ampl, ASL_alloc, [acx_ampl_ok=yes; AMPL_LIBS="-lamplsolver"])
fi

AC_SUBST(AMPL_LIBS)
AC_SUBST(AMPL_INCLUDE)

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_ampl_ok" = xyes; then
        ifelse([$1],,AC_DEFINE(HAVE_AMPL,1,[Define if you have a AMPL library.]),[$1])
        :
else
        acx_ampl_ok=no
        $2
fi
])dnl ACX_AMPL
