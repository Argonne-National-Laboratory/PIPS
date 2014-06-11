
# To build shared library use:
#     make ENABLE_SHARED=1

ifeq ($(ENABLE_SHARED),1)

FPIC = -fPIC

all: libMA57.so

else

FPIC =

all: libMA57.a

endif

#F77	= f95
#F77	 = mpixlf77_r
F77 = mpif77 -fPIC
#F77     = /soft/mvapich2/1.4-gcc-4.4.2/bin/mpif77

FC= $(F77)

FFLAGS   = -O2 $(FPIC)

OBJS = ma57.o

libMA57.so:  $(OBJS) Makefile
	gcc -shared $(FFLAGS) -o $@ $(OBJS)

libMA57.a:  $(OBJS) Makefile
	ar ru $@ $(OBJS)

$(OBJS) : Makefile
