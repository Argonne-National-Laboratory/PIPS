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
 
SUBDIRS = src interfaces
EXTRA_DIST = $(DOC_FILES) $(TEST_FILES) $(THIRDPARTY_FILES) \
				 gpl.txt lgpl.txt

DOC_FILES  = doc/userguide.pdf doc/userguide.bib doc/userguide.tex \
				 doc/fig/Interface_Example.fig doc/fig/Interface_Example.jpg \
				 doc/fig/MSND_struct.fig doc/fig/MSND_struct.jpg \
				 doc/Makefile
TEST_FILES = tests/alm_xi.dat tests/alm_xi.mod \
		tests/alm.mod tests/alm.dat \
		tests/newsvendor.mod tests/newsvendor.dat \
		tests/ampl_msnd_test_slck.mod tests/ampl_msnd_test_slck.dat
# scopfIEEE model does not converge under OOPS due to numeric
# issues in the dual space. Not distributed as a good example.
#tests/scopfIEEE_ampl.mod tests/scopfIEEE.dat tests/scopfIEEE.mod
THIRDPARTY_FILES = \
				ThirdParty/oops/COPYING \
				ThirdParty/oops/lib/liboops.a \
				ThirdParty/oops/include/oops/Algebra.h \
				ThirdParty/oops/include/oops/BlockDenseAlg.h \
				ThirdParty/oops/include/oops/BlockDiagSimpleAlg.h \
				ThirdParty/oops/include/oops/CallBack.h \
				ThirdParty/oops/include/oops/DblBordDiagSimpleAlg.h \
				ThirdParty/oops/include/oops/DenseVector.h \
				ThirdParty/oops/include/oops/GlobalOpt.h \
				ThirdParty/oops/include/oops/hopdm.h \
				ThirdParty/oops/include/oops/LogiVector.h \
				ThirdParty/oops/include/oops/MatrixDense.h \
				ThirdParty/oops/include/oops/MatrixSparseSimple.h \
				ThirdParty/oops/include/oops/OopsInterface.h \
				ThirdParty/oops/include/oops/ParAlgList.h \
				ThirdParty/oops/include/oops/parutil.h \
				ThirdParty/oops/include/oops/SparseVector.h \
				ThirdParty/oops/include/oops/StrSparseVector.h \
				ThirdParty/oops/include/oops/Tree.h \
				ThirdParty/oops/include/oops/Vector.h \
				ThirdParty/oops/include/oops/WriteMps.h
