# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

## $Id: Makefile.am 255 2011-07-16 13:54:07Z stefan $

# Author:  Andreas Waechter           IBM    2006-04-13

AUTOMAKE_OPTIONS = foreign

########################################################################
#                List files that should be distributed                 #
########################################################################

EXTRA_DIST = $(EXAMPLE_FILES) README AUTHORS LICENSE

DISTCLEANFILES = $(EXAMPLE_CLEAN_FILES)

datacoindir = $(datadir)/coin/Data/Sample
datacoin_DATA = $(EXAMPLE_FILES)

coindocdir = $(datadir)/coin/doc/Data/Sample
coindoc_DATA = README AUTHORS LICENSE

pkgconfiglibdir = $(libdir)/pkgconfig
pkgconfiglib_DATA = coindatasample.pc

include BuildTools/Makemain.inc

test:
	@echo "No test available."
