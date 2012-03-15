# (c) Marco Colombo 2009
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

AM_CXXFLAGS = -I$(srcdir)/../../src @CPLEX_INCLUDE@

bin_PROGRAMS = smlcplex
smlcplex_SOURCES = sml-cplex.cpp sml-cplex.h wrapper.cpp
smlcplex_LDADD = ../../src/libsml.a @AMPL_LIBS@ -lfl @CPLEX_LIBS@
