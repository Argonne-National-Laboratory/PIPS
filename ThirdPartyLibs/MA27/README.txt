# 
# Please Get MA27 from HSL (http://www.hsl.rl.ac.uk)
# Note that we need double precision FORTRAN source code.
#
# After you have the MA27 source code ma27.f, please put it in the folder src,
# and use the attached Makefile to install MA27
#

F77      = mpif77
FC       = $(F77)
FFLAGS   = -O2 -fPIC

OBJS = ma27.o 

libMA27.a:  $(OBJS) Makefile
	ar ru $@ $(OBJS)

$(OBJS) : Makefile

