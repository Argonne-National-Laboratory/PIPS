# Copyright (C) 2006, 2008 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

## $Id: Makefile.am 2313 2011-06-19 14:33:08Z stefan $

AUTOMAKE_OPTIONS = foreign

EXTRA_DIST = get.lapack INSTALL.LAPACK

lib_LTLIBRARIES = libcoinlapack.la

LAPACK_FILES = \
	dbdsqr.f \
	dgebd2.f \
	dgebrd.f \
	dgelq2.f \
	dgelqf.f \
	dgeqr2.f \
	dgeqrf.f \
	dgesvd.f \
	dgesv.f \
	dgetf2.f \
	dgetrf.f \
	dgetrs.f \
	disnan.f \
	dlabrd.f \
	dlacpy.f \
	dlae2.f \
	dlaev2.f \
	dlaisnan.f \
	dlamch.f \
	dlange.f \
	dlanst.f \
	dlansy.f \
	dlapy2.f \
	dlarf.f \
	dlarfb.f \
	dlarfg.f \
	dlarfp.f \
	dlarft.f \
	dlartg.f \
	dlas2.f \
	dlascl.f \
	dlaset.f \
	dlasq1.f \
	dlasq2.f \
	dlasq3.f \
	dlasq4.f \
	dlasq5.f \
	dlasq6.f \
	dlasr.f \
	dlasrt.f \
	dlaswp.f \
	dlassq.f \
	dlasv2.f \
	dlatrd.f \
	dorg2l.f \
	dorg2r.f \
	dorgbr.f \
	dorglq.f \
	dorgl2.f \
	dorgql.f \
	dorgqr.f \
	dorgtr.f \
	dorm2r.f \
	dormbr.f \
	dorml2.f \
	dormlq.f \
	dormqr.f \
	dpotf2.f \
	dpotrf.f \
	dpotrs.f \
	dsteqr.f \
	dsterf.f \
	dsyev.f \
	dsytd2.f \
	dsytrd.f \
	ieeeck.f \
	iladlc.f \
	iladlr.f \
	ilaenv.f \
	iparmq.f

libcoinlapack_la_SOURCES = $(LAPACK_FILES)

libcoinlapack_la_LDFLAGS = $(LT_LDFLAGS)

dlamch.o: dlamch.f
	$(F77) $(AM_FFLAGS) $(DLAMCH_FFLAGS) -c -o $@ $<

dlamch.obj: dlamch.f
	$(F77) $(AM_FFLAGS) $(DLAMCH_FFLAGS) -c -o $@ `$(CYGPATH_W) '$<'`

dlamch.lo: dlamch.f
	$(LIBTOOL) --tag=F77 --mode=compile $(F77) $(AM_FFLAGS) $(DLAMCH_FFLAGS) -c -o $@ $<

pkgconfiglibdir = $(libdir)/pkgconfig
pkgconfiglib_DATA = coinlapack.pc

test:
	@echo "No test available for Lapack."

########################################################################
#                         Maintainer Stuff                             #
########################################################################

CLEANFILES =

# Files that are generated and should be cleaned with make distclean
DISTCLEANFILES = $(LAPACK_FILES) $(MAKEOKFILE)

include BuildTools/Makemain.inc
