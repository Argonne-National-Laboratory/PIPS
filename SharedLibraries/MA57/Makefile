F77	= f95
#F77	 = mpixlf77_r
#F77     = /soft/mvapich2/1.4-gcc-4.4.2/bin/mpif77

FC= $(F77)

FFLAGS   = -O2

OBJS = ma57.o 

libMA57.a:  $(OBJS) Makefile
	ar ru $@ $(OBJS)

$(OBJS) : Makefile
