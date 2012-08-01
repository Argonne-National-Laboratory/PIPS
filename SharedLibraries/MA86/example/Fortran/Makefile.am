AM_FCFLAGS = -I $(top_builddir)/src

TESTS_ENVIRONMENT = $(srcdir)/runtest.sh
EXTRA_DIST = runtest.sh
TESTS = \
	hsl_ma86ss \
	hsl_ma86ds \
	hsl_ma86cs \
	hsl_ma86zs \
	hsl_ma86ss1 \
	hsl_ma86ds1 \
	hsl_ma86cs1 \
	hsl_ma86zs1

check_PROGRAMS = \
	hsl_ma86ss \
	hsl_ma86ds \
	hsl_ma86cs \
	hsl_ma86zs \
	hsl_ma86ss1 \
	hsl_ma86ds1 \
	hsl_ma86cs1 \
	hsl_ma86zs1

hsl_ma86ss_SOURCES = hsl_ma86ss.f90
hsl_ma86ss_LDADD = $(top_builddir)/src/libhsl_ma86.a $(LAPACK_LIBS) $(BLAS_LIBS) $(METIS_LIBS) $(FCLIBS)
EXTRA_DIST += hsl_ma86ss.output hsl_ma86ss.data
hsl_ma86ds_SOURCES = hsl_ma86ds.f90
hsl_ma86ds_LDADD = $(top_builddir)/src/libhsl_ma86.a $(LAPACK_LIBS) $(BLAS_LIBS) $(METIS_LIBS) $(FCLIBS)
EXTRA_DIST += hsl_ma86ds.output hsl_ma86ds.data
hsl_ma86cs_SOURCES = hsl_ma86cs.f90
hsl_ma86cs_LDADD = $(top_builddir)/src/libhsl_ma86.a $(LAPACK_LIBS) $(BLAS_LIBS) $(METIS_LIBS) $(FCLIBS)
EXTRA_DIST += hsl_ma86cs.output hsl_ma86cs.data
hsl_ma86zs_SOURCES = hsl_ma86zs.f90
hsl_ma86zs_LDADD = $(top_builddir)/src/libhsl_ma86.a $(LAPACK_LIBS) $(BLAS_LIBS) $(METIS_LIBS) $(FCLIBS)
EXTRA_DIST += hsl_ma86zs.output hsl_ma86zs.data
hsl_ma86ss1_SOURCES = hsl_ma86ss1.f90
hsl_ma86ss1_LDADD = $(top_builddir)/src/libhsl_ma86.a $(LAPACK_LIBS) $(BLAS_LIBS) $(METIS_LIBS) $(FCLIBS)
EXTRA_DIST += hsl_ma86ss1.output hsl_ma86ss1.data
hsl_ma86ds1_SOURCES = hsl_ma86ds1.f90
hsl_ma86ds1_LDADD = $(top_builddir)/src/libhsl_ma86.a $(LAPACK_LIBS) $(BLAS_LIBS) $(METIS_LIBS) $(FCLIBS)
EXTRA_DIST += hsl_ma86ds1.output hsl_ma86ds1.data
hsl_ma86cs1_SOURCES = hsl_ma86cs1.f90
hsl_ma86cs1_LDADD = $(top_builddir)/src/libhsl_ma86.a $(LAPACK_LIBS) $(BLAS_LIBS) $(METIS_LIBS) $(FCLIBS)
EXTRA_DIST += hsl_ma86cs1.output hsl_ma86cs1.data
hsl_ma86zs1_SOURCES = hsl_ma86zs1.f90
hsl_ma86zs1_LDADD = $(top_builddir)/src/libhsl_ma86.a $(LAPACK_LIBS) $(BLAS_LIBS) $(METIS_LIBS) $(FCLIBS)
EXTRA_DIST += hsl_ma86zs1.output hsl_ma86zs1.data
