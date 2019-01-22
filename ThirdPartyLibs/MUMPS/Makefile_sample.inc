#
#  This file is part of MUMPS 5.1.2, released
#  on Mon Oct  2 07:37:01 UTC 2017
#
# These settings for a PC under Debian/linux with standard packages :
# metis (parmetis), scotch (ptscotch), openmpi, gfortran

# packages installation : 
# apt-get install libmetis-dev libparmetis-dev libscotch-dev libptscotch-dev libatlas-base-dev openmpi-bin libopenmpi-dev lapack-dev
PIPS_DIR=$(HOME)/git/PIPS

# Begin orderings
#LSCOTCHDIR = /usr/lib
#ISCOTCH   = -I/usr/include/scotch

#LSCOTCH   = -L$(LSCOTCHDIR) -lptesmumps -lptscotch -lptscotcherr
#LSCOTCH   = -L$(LSCOTCHDIR) -lesmumps -lscotch -lscotcherr

LPORDDIR = $(topdir)/PORD/lib/
IPORD    = -I$(topdir)/PORD/include/
LPORD    = -L$(LPORDDIR) -lpord

#LMETISDIR = /usr/lib 
LMETISDIR = $(PIPS_DIR)/ThirdPartyLibs/METIS/lib
#IMETIS    = -I/usr/include/parmetis
IMETIS    = -I$(PIPS_DIR)/ThirdPartyLibs/METIS/include -I$(PIPS_DIR)/ThirdPartyLibs/METIS/parmetis-4.0.3/metis/include

# LMETIS    = -L$(LMETISDIR) -lparmetis -lmetis
LMETIS    = -L$(LMETISDIR) -lparmetis -lmetis -lscalapack-openmpi -lblacs-openmpi  -lblacsF77init-openmpi -lblacsCinit-openmpi

# Corresponding variables reused later
#ORDERINGSF = -Dmetis -Dpord -Dparmetis #-Dscotch -Dptscotch
ORDERINGSF = -Dparmetis -Dpord #-Dscotch
ORDERINGSC  = $(ORDERINGSF)

LORDERINGS = $(LMETIS) $(LPORD) #$(LSCOTCH)
IORDERINGSF = #$(ISCOTCH)
IORDERINGSC = $(IMETIS) $(IPORD) #$(ISCOTCH)
# End orderings
################################################################################

PLAT    =
LIBEXT  = .a
OUTC    = -o 
OUTF    = -o 
RM = /bin/rm -f
CC = mpicc
FC = mpif90
FL = mpif90
AR = ar vr 
RANLIB = ranlib
LAPACK = -llapack
SCALAP  = -lscalapack-openmpi -lblacs-openmpi  -lblacsF77init-openmpi -lblacsCinit-openmpi

LIBPAR = $(SCALAP) $(LAPACK)

INCSEQ = -I$(topdir)/libseq
LIBSEQ  = $(LAPACK) -L$(topdir)/libseq -lmpiseq

LIBBLAS = -lblas
LIBOTHERS = -lpthread

#Preprocessor defs for calling Fortran from C (-DAdd_ or -DAdd__ or -DUPPER)
CDEFS   = -DAdd_

#Begin Optimized options
# uncomment -fopenmp in lines below to benefit from OpenMP
OPTF    = -fPIC -O3  -fopenmp
OPTL    = -fPIC -O3  -fopenmp
OPTC    = -fPIC -O3  -fopenmp
#End Optimized options
 
INCS = $(INCPAR)
LIBS = $(LIBPAR)
LIBSEQNEEDED = 
