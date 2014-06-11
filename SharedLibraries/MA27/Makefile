F77      = mpif77
FC       = $(F77)
FFLAGS   = -O

OBJS =   ma27.o 

libMA27.a:  $(OBJS) Makefile
	ar ru $@ $(OBJS)

$(OBJS) : Makefile

