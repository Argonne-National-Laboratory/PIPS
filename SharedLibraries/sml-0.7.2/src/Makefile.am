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

AM_YFLAGS = -d
#AM_YFLAGS += -t -v
AM_LEXFLAGS = 
#AM_LEXFLAGS += -d
AM_CXXFLAGS = @AMPL_INCLUDE@

lib_LIBRARIES = libsml.a

libsml_a_SOURCES = AmplModel.cpp AmplModel.h AmplsolverCalls.cpp 		 	\
						 AmplsolverCalls.h backend.cpp backend.h       		 	\
						 CompDescr.h CompDescrParam.cpp CompDescrParam.h    	\
						 DataNodes.h ExpandedModel.cpp ExpandedModel.h 		 	\
						 GlobalVariables.h misc.h ModelComp.cpp ModelComp.h 	\
						 ExpandedModelInterface.h nodes.cpp nodes.h  			\
						 SetElement.cpp SetElement.h Set.cpp Set.h            \
						 SetNode.cpp SetNode.h sml.cpp sml.h  						\
						 StochModelComp.cpp StochModelComp.h StochModel.cpp   \
						 StochModel.h symtab.cpp symtab.h                     \
						 sml.tab.ypp sml.yy.lpp data.tab.ypp data.yy.lpp

noinst_HEADERS = doxygen.h
include_HEADERS = ExpandedModelInterface.h sml.h

# Force these to get built first
BUILT_SOURCES = sml.tab.h data.tab.h

########################################################################
# YACC and LEX:
# The proper automake way is confusing and annoying, lets do it our way.
########################################################################
data.yy.cpp: data.yy.lpp
	$(SHELL) $(YLWRAP) $< lex.data.c $@ -- $(LEXCOMPILE) -Pdata

data.tab.cpp: data.tab.ypp
	$(SHELL) $(YLWRAP) $< y.tab.c $@ y.tab.h $*.h y.output $*.output -- $(YACCCOMPILE) -pdata
