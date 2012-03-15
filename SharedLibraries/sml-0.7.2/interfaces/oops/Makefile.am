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

AM_CXXFLAGS = -I$(srcdir)/../../src @OOPS_INCLUDE@

bin_PROGRAMS = smloops
smloops_SOURCES = wrapper.cpp OOPSBlock.cpp sml-oops.cpp OOPSBlock.h \
						sml-oops.h
smloops_LDADD = ../../src/libsml.a @AMPL_LIBS@ \
		@OOPS_LIB@ @METIS_LIB@ @LAPACK_LIBS@ @BLAS_LIBS@ -lfl
