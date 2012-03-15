# (c) Jonathan Hogg 2009
#
# This file is part of SML.
#
# SML is free software: you can redistribute it and/or modify it under the
# terms of the GNU Lesser General Public License as published by the Free
# Software Foundation, using version 3 of the License.
#
# SML is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program. If not, see http://www.gnu.org/licenses/.

AC_INIT([SML],[0.7.2],[sml-support@maths.ed.ac.uk])
AC_CONFIG_AUX_DIR([config])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_SRCDIR([src/sml.tab.ypp])
AM_INIT_AUTOMAKE
AC_PROG_CXX
AC_PROG_LEX
AC_PROG_YACC
AC_PROG_RANLIB
AC_PROG_INSTALL

# generate nicer looking output from the makefiles with automake 1.11
# (run make V=1 to see all the details)
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])

# check for headers and functions not available under MinGW
AC_CHECK_HEADERS(direct.h)
AC_CHECK_FUNCS(strndup)

# amplsolver is required
ACX_AMPL( [], AC_ERROR([Ampl Solver Library Required]) )

# OOPS interface requires OOPS and dependant libraries
ACX_BLAS(HAVE_BLAS=yes,HAVE_BLAS=no)
ACX_LAPACK(HAVE_LAPACK=yes,HAVE_LAPACK=no)
IMMDX_LIB_METIS( [HAVE_METIS=yes], [] )
ACX_OOPS( [HAVE_OOPS=yes], [] )
if test "$HAVE_OOPS" = "yes" ; then
  AC_SUBST(OOPS_INTERFACE, [oops])
fi

# Cplex interface
ACX_CPLEX( [HAVE_CPLEX=yes], [] )
if test "$HAVE_CPLEX" = "yes" ; then
  AC_SUBST(CPLEX_INTERFACE, [cplex])
fi

AC_ARG_WITH(
   [mps],
   [AC_HELP_STRING([--with-mps],
      [enable SML-MPS interface])
   ],
   [
      if test "$withval" == yes ; then
         AC_SUBST(MPS_INTERFACE, [mps])
      fi
   ],
   [
      AC_SUBST(MPS_INTERFACE, [mps])
   ])


AC_CONFIG_FILES([Makefile src/Makefile interfaces/Makefile
                 interfaces/oops/Makefile interfaces/mps/Makefile
                 interfaces/cplex/Makefile])
AC_OUTPUT

echo "==========================================="
echo "Summary:"
echo "Interfaces to build: $MPS_INTERFACE $OOPS_INTERFACE $CPLEX_INTERFACE"
