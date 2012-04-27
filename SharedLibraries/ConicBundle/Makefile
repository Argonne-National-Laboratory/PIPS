#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*   Type....: Makefile                                                      *
#*   File....: makefile                                                      *
#*   Name....: DFN Pre makefile                                              *
#*   Author..: Thorsten Koch, modified by Christoph Helmberg                 *
#*   Copyright by Authors, All rights reserved                               *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

ARCH		:=	$(shell uname -m | sed -e s/sun../sparc/)
OSTYPE		:=	$(shell uname -s | tr A-Z a-z)
CXX		=	g++
CC		=	gcc

MODE            =       DEBU
#MODE           =       OPTI

CONICBUNDLE	=	.
CPPFLAGS	=	-I$(CONICBUNDLE)/include -I$(CONICBUNDLE)/CBsources \
			-I$(CONICBUNDLE)/Matrix -I$(CONICBUNDLE)/Tools 

CBLIBOBJECT	=	memarray.o \
			CBSolver.o MatCBSolver.o \
			CFunction.o CB_CSolver.o \
                        MatNBSolver.o Nsumproblem.o Nfunproblem.o Nbundle.o \
                        MatFCBSolver.o problem.o sumproblem.o funproblem.o \
			bundle.o hkweight.o qp_solver.o qp_sdpblock.o \
			idscaling.o diagscaling.o fullscaling.o lowrankscaling.o \
			diagtrscaling.o fulltrscaling.o lowranktrscaling.o \
			lowrankSMscaling.o lowranktrSMscaling.o \
			MatConcavefun.o \
			MatLPfun.o MatLPBCfun.o \
			MatConefun.o coneproblem.o \
			MatSDPfun.o lmaxproblem.o coeffmat.o bigmat.o \
			MatSOCfun.o socproblem.o \
			indexmat.o matrix.o symmat.o  eigval.o ldl.o chol.o \
			qr.o trisolve.o nnls.o sparssym.o sparsmat.o lanczpol.o \

CTESTOBJECT	=	c_main.o 


CXXTESTOBJECT	=	cxx_main.o

MATTESTOBJECT	=	mat_main.o

TARGET		=	lib/libcb.a t_c t_cxx t_mat

#-----------------------------------------------------------------------------

GCCWARN		=	-W -Wall -pedantic -Wcast-qual -Wwrite-strings \
			-Wnon-virtual-dtor -Wcast-align -Wconversion \
			-Wno-char-subscripts -Wpointer-arith -Wundef

#-----------------------------------------------------------------------------

# Echo Settings for backslash escapes on different ostypes, most of them with -e
# only darwin on MAC without -e
ECHO.darwin  = 
ECHO.linux   = -e
ECHO.irix    = -e
ECHO.sunos   = -e

ECHOFLAGS    =  $(ECHO.$(OSTYPE))

#--- Linux.i686.g++ settings ---------------------------------------------------
DEBU.linux.i686.g++ = 	-g  
OPTI.linux.i686.g++ =   -DNDEBUG -O3
WARN.linux.i686.g++ =	$(GCCWARN)
DEPD.linux.i686.g++ =	-MM
LINK.linux.i686.g++ =	-lm
AR.linux.i686.g++   =	ar
ARFLAGS.linux.i686.g++ =	cr
RANLIB.linux.i686.g++ =	ranlib
OPTI.linux.i686.gcc =	-DNDEBUG -O3
WARN.linux.i686.gcc =	$(GCCWARN)
DEBU.linux.i686.gcc = 	-g 
DEPD.linux.i686.gcc =	-MM
LINK.linux.i686.gcc =	-lm
#--- linux.x86_64.g++ settings ---------------------------------------------------
DEBU.linux.x86_64.g++ = -g 
OPTI.linux.x86_64.g++ = -DNDEBUG -O3  
WARN.linux.x86_64.g++ =	$(GCCWARN)
DEPD.linux.x86_64.g++ =	-MM
LINK.linux.x86_64.g++ =	-lm 
AR.linux.x86_64.g++   =	ar
ARFLAGS.linux.x86_64.g++ =	cr
RANLIB.linux.x86_64.g++ =	ranlib
OPTI.linux.x86_64.gcc =	-DNDEBUG -O3
WARN.linux.x86_64.gcc =	$(GCCWARN)
DEBU.linux.x86_64.gcc = 	-g 
DEPD.linux.x86_64.gcc =	-MM
LINK.linux.x86_64.gcc =	-lm
#--- MAC Darwin.i386.g++ settings ---------------------------------------------------
DEBU.darwin.i386.g++ = 	-g  
OPTI.darwin.i386.g++ =   -DNDEBUG -O3
WARN.darwin.i386.g++ =	$(GCCWARN)
DEPD.darwin.i386.g++ =	-MM
LINK.darwin.i386.g++ =	-lm
AR.darwin.i386.g++   =	ar
ARFLAGS.darwin.i386.g++ =	cr
RANLIB.darwin.i386.g++ =	ranlib
OPTI.darwin.i386.gcc =	-DNDEBUG -O3
WARN.darwin.i386.gcc =	$(GCCWARN)
DEBU.darwin.i386.gcc = 	-g 
DEPD.darwin.i386.gcc =	-MM
LINK.darwin.i386.gcc =	-lm
#--- linux.alpha.g++ settings ---------------------------------------------------
OPTI.linux.alpha.g++ =	-DNDEBUG -pipe -O6 -fomit-frame-pointer -fschedule-insns2 -funroll-loops -felide-constructors 
WARN.linux.alpha.g++ =	$(GCCWARN)
DEBU.linux.alpha.g++ =	-g 
DEPD.linux.alpha.g++ =	-MM
LINK.linux.alpha.g++ =	-lm
AR.linux.alpha.g++   =	ar
ARFLAGS.linux.alpha.g++ =	cr
RANLIB.linux.alpha.g++ =	ranlib
OPTI.linux.alpha.gcc =	-DNDEBUG -pipe -O6 -fomit-frame-pointer -fschedule-insns2 -funroll-loops -felide-constructors 
WARN.linux.alpha.gcc =	$(GCCWARN)
DEBU.linux.alpha.gcc =	 -g 
DEPD.linux.alpha.gcc =	-MM
LINK.linux.alpha.gcc =	-lm
#--- linux.alpha.cxx settings ---------------------------------------------------
OPTI.linux.alpha.cxx =	-DNDEBUG -fast 
WARN.linux.alpha.cxx =	-std strict_ansi
DEBU.linux.alpha.cxx =	 -g 
DEPD.linux.alpha.cxx =	-M -std strict_ansi -noimplicit_include
LINK.linux.alpha.cxx =	-lm
AR.linux.alpha.cxx   =	ar
ARFLAGS.linux.alpha.cxx =	cr
RANLIB.linux.alpha.cxx =	ranlib
OPTI.linux.alpha.ccc =	-DNDEBUG -fast 
WARN.linux.alpha.ccc =	
DEBU.linux.alpha.ccc =	 -g 
DEPD.linux.alpha.ccc =	-M
LINK.linux.alpha.ccc =	-lm
#--- Sun.Solaris.g++ settings ------------------------------------------------
DEBU.sunos.sparc.g++ =	-g 
OPTI.sunos.sparc.g++ =	-DNDEBUG -fPIC -pipe -O2 -fomit-frame-pointer -fschedule-insns2 -funroll-loops -felide-constructors 
WARN.sunos.sparc.g++ =	$(GCCWARN)
DEPD.sunos.sparc.g++ =	-MM
LINK.sunos.sparc.g++ =	-lm
AR.sunos.sparc.g++   =	ar
ARFLAGS.sunos.sparc.g++ =	cr
RANLIB.sunos.sparc.g++ =	ranlib
OPTI.sunos.sparc.gcc =	-DNDEBUG -fPIC -pipe -O2 -fomit-frame-pointer -fschedule-insns2 -funroll-loops -mcpu=supersparc
WARN.sunos.sparc.gcc =	$(GCCWARN)
DEBU.sunos.sparc.gcc =	-g 
DEPD.sunos.sparc.gcc =	-MM
LINK.sunos.sparc.gcc =	-lm 
#--- Sun.Solaris.CC settings -------------------------------------------------
DEBU.sunos.sparc.CC=	
OPTI.sunos.sparc.CC=	-DNDEBUG -fast 
WARN.sunos.sparc.CC=
DEPD.sunos.sparc.CC=	-xM1
LINK.sunos.sparc.CC=	-lm
AR.sunos.sparc.CC =	CC
ARFLAGS.sunos.sparc.CC =	-xar -o
RANLIB.sunos.sparc.CC =	ranlib
OPTI.sunos.sparc.cc=	-DNDEBUG -fast 
WARN.sunos.sparc.cc=
DEBU.sunos.sparc.cc=	-g 
DEPD.sunos.sparc.cc=	-xM1
LINK.sunos.sparc.cc=	-lm
#--- O.irix.IP32.CC settings ---------------------------------------------------
OPTI.irix.IP32.CC=	-DNDEBUG -Ofast -OPT:Olimit=3000
WARN.irix.IP32.CC=	-LANG:std:ansi-for-init-scope
DEBU.irix.IP32.CC=	-g 
DEPD.irix.IP32.CC=	-M
LINK.irix.IP32.CC=	-lCio -lm
AR.irix.IP32.CC   =	CC
ARFLAGS.irix.IP32.CC =	-ar -o
RANLIB.irix.IP32.CC =	true
OPTI.irix.IP32.cc=	-DNDEBUG -Ofast -OPT:Olimit=3000
WARN.irix.IP32.cc=	
DEBU.irix.IP32.cc=	-g 
DEPD.irix.IP32.cc=	-M
LINK.irix.IP32.cc=	-lm
#-----------------------------------------------------------------------------

CXXFLAGS        =       $($(MODE).$(OSTYPE).$(ARCH).$(CXX))\
			$(WARN.$(OSTYPE).$(ARCH).$(CXX))

CCFLAGS		= 	$($(MODE).$(OSTYPE).$(ARCH).$(CC))\
			$(WARN.$(OSTYPE).$(ARCH).$(CC))

LDFLAGS		=	$(LINK.$(OSTYPE).$(ARCH).$(CXX))
DFLAGS		=	$(DEPD.$(OSTYPE).$(ARCH).$(CXX)) 

AR		=       $(AR.$(OSTYPE).$(ARCH).$(CXX))
ARFLAGS		=	$(ARFLAGS.$(OSTYPE).$(ARCH).$(CXX))
RANLIB		=	$(RANLIB.$(OSTYPE).$(ARCH).$(CXX))

OBJDIR		=	$(MODE).$(OSTYPE).$(ARCH).$(CXX)
OBJCTEST	=	$(addprefix $(OBJDIR)/,$(CTESTOBJECT))
OBJCXXTEST	=	$(addprefix $(OBJDIR)/,$(CXXTESTOBJECT))
OBJMATTEST	=	$(addprefix $(OBJDIR)/,$(MATTESTOBJECT))
OBJCBLIB	=	$(addprefix $(OBJDIR)/,$(CBLIBOBJECT))

VPATH	        =       . $(CONICBUNDLE)/Matrix $(CONICBUNDLE)/CBsources

all:		$(TARGET)

t_c:		$(OBJCTEST) lib/libcb.a
		$(CXX) $(CXXFLAGS) $(OBJCTEST) -Llib -lcb $(LDFLAGS) -o $@

t_cxx:		$(OBJCXXTEST) lib/libcb.a
		$(CXX) $(CXXFLAGS) $(OBJCXXTEST) -Llib -lcb $(LDFLAGS)  -o $@

t_mat:		$(OBJMATTEST) lib/libcb.a
		$(CXX) $(CXXFLAGS) $(OBJMATTEST) -Llib -lcb $(LDFLAGS)  -o $@

lib/libcb.a:   	include/CBconfig.hxx $(OBJCBLIB)
		@if [ ! -d lib ]; then mkdir lib; fi
	        $(AR) $(ARFLAGS) lib/libcb.a $(OBJCBLIB) 
		$(RANLIB) lib/libcb.a

clean:
		-rm -rf OPTI.* DEBU.* $(TARGET)

$(OBJDIR)/%.o:	%.cxx
		@if [ ! -d $(OBJDIR) ]; then mkdir $(OBJDIR); fi
		$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

$(OBJDIR)/%.o:	%.c
		@if [ ! -d $(OBJDIR) ]; then mkdir $(OBJDIR); fi
		$(CC) $(CPPFLAGS) $(CCFLAGS) -c $< -o $@

$(OBJDIR)/%.d: %.c
		@if [ ! -d $(OBJDIR) ]; then mkdir $(OBJDIR); fi
		$(SHELL) -ec '$(CC) $(DFLAGS) $(CPPFLAGS) $< \
		| sed '\''s/\($*\)\.o[ :]*/$(OBJDIR:/=\/)\/\1.o $(OBJDIR:/=\/)\/\1.d : /g'\'' >$@; \
		[ -s $@ ] || rm -f $@'

$(OBJDIR)/%.d: %.cxx
		@if [ ! -d $(OBJDIR) ]; then mkdir $(OBJDIR); fi
		$(SHELL) -ec '$(CXX) $(DFLAGS) $(CPPFLAGS) $< \
		| sed '\''s/\($*\)\.o[ :]*/$(OBJDIR:/=\/)\/\1.o $(OBJDIR:/=\/)\/\1.d : /g'\'' >$@; \
		[ -s $@ ] || rm -f $@'

include/CBconfig.hxx: Makefile
ifeq ($(MODE),OPTI)
		rm -f include/CBconfig.hxx
		echo $(ECHOFLAGS) "#ifndef __CBCONFIG_HXX__"> include/CBconfig.hxx
		echo $(ECHOFLAGS) "#define __CBCONFIG_HXX__">> include/CBconfig.hxx
		echo $(ECHOFLAGS) "#define CONICBUNDLE_DEBUG 0">> include/CBconfig.hxx
		echo $(ECHOFLAGS) "#endif" >> include/CBconfig.hxx
else
		rm -f include/CBconfig.hxx
		echo $(ECHOFLAGS) "#ifndef __CBCONFIG_HXX__"> include/CBconfig.hxx
		echo $(ECHOFLAGS) "#define __CBCONFIG_HXX__">> include/CBconfig.hxx
		echo $(ECHOFLAGS) "#define CONICBUNDLE_DEBUG 1">> include/CBconfig.hxx
		echo $(ECHOFLAGS) "#endif" >> include/CBconfig.hxx
endif


#--- for some compilers it is not helpful to generate dependencies
#    automatically, therefore the dependencies are given explicitly

include depend

#include		$(OBJCTEST:.o=.d)
#include		$(OBJCXXTEST:.o=.d)
#include		$(OBJMATTEST:.o=.d)
#include		$(OBJCBLIB:.o=.d)
