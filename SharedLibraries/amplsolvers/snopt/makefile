# /****************************************************************
# Copyright (C) 1997-1998 Lucent Technologies
# All Rights Reserved
#
# Permission to use, copy, modify, and distribute this software and
# its documentation for any purpose and without fee is hereby
# granted, provided that the above copyright notice appear in all
# copies and that both that the copyright notice and this
# permission notice and warranty disclaimer appear in supporting
# documentation, and that the name of Lucent or any of its entities
# not be used in advertising or publicity pertaining to
# distribution of the software without specific, written prior
# permission.
#
# LUCENT DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
# INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS.
# IN NO EVENT SHALL LUCENT OR ANY OF ITS ENTITIES BE LIABLE FOR ANY
# SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
# WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER
# IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION,
# ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF
# THIS SOFTWARE.
# ****************************************************************/

.SUFFIXES: .c .o
CC = cc
CFLAGS = -O
S = ..
# "S = .." assumes this directory is solvers/snopt .

# See $S/makefile for comments about CFLAGS.

.c.o:
	$(CC) -c $(CFLAGS) -I$S $*.c

# The following assumes compilation of Fortran 77 with f2c + cc.
# To use another Fortran compiler for libamplsnopt.a, see $S/README.f77
# and adjust the snopt rule below appropriately, e.g., to
#	$(F77) -o snopt $(snopt)
# If your Fortran compiler does not append underscores to external names,
# change snopt.o to snopt_.o below.

snopt = snopt.o snoptver.o libamplsnopt.a $S/amplsolver.a
snopt: $(snopt)
	$(CC) -o snopt $(snopt) -lf2c -lm

# If necessary, change `uname -sr` in the following sed command to a
# phrase summarizing your machine and operating system.
snoptver.c: snoptver.c0
	sed "s/Machine_details/`uname -sr`/" snoptver.c0 >snoptver.c

# make xsum.out to check for transmission errors.
# This assumes you have the xsum program, whose source
# you can get by asking netlib@netlib.bell-labs.com to
#	send xsum.c from f2c/src

xs0 = README makefile snopt.c snopt_.c snoptver.c0

xsum.out: $(xs0)
	xsum $(xs0) >xsum1.out
	cmp xsum0.out xsum1.out && mv xsum1.out xsum.out || diff xsum[01].out
