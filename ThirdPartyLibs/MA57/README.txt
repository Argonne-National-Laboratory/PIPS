# 
# Please Get MA57 from HSL (http://www.hsl.rl.ac.uk)
# Note that we need double precision FORTRAN source code.
#
# After you have the MA57 source code ma57.f, please put it in the folder src,
# and use the attached Makefile to install MA57

F77      = mpif77 
FC       = $(F77)
FFLAGS   = -O2 -fPIC

OBJS = ma57.o

libMA57.a:  $(OBJS) Makefile
	ar ru $@ $(OBJS)

$(OBJS) : Makefile
