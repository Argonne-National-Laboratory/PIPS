# Specify extensions of files to delete when cleaning
CLEANEXTS   = o a 

# Specify the target file 
OUTPUTFILE  = libstatgdx.a
CXX = mpicxx
CXXFLAGS = -D_GCL_RHACK_ -D_P3_DllExport= -DSTDCALL= -DGDX_CALLCONV=unwanted -DLEG -Wall -Wno-long-long -fpermissive
CFLAGS = -D_GCL_RHACK_ -D_P3_DllExport= -DSTDCALL= -DGDX_CALLCONV=unwanted -DLEG -Wall -Wno-long-long  
# Default target
.PHONY: all
all: $(OUTPUTFILE) xp1

RUNNER = runner.o
RUNNER =

$(OUTPUTFILE): clibtypes.o gmsheapnew.o p3Custom2.o paldoorg.o datastorage.o gmslibname.o p3ieeefp.o pchutil.o gdlaudit.o gmsobj.o p3io.o gdxdclib.o gmsspecs.o p3library.o strhash.o gdxstatic.o gmsstrm.o p3platform.o strutilx.o gmsdata.o gxdefs.o p3private.o system_p3.o gmsgen.o gxfile.o p3process.o sysutils_p3.o gmsglob.o idglobal_p3.o p3threads.o gmsglobx.o math_p3.o p3utils.o xcompress.o gclgms.o utils2.o $(RUNNER)
	ar ru $@ $^
	ranlib $@

xp1: $(OUTPUTFILE) xp1.o
	g++ -m64 -o xp1 xp1.o libstatgdx.a -ldl

.PHONY: clean 
clean:
	for file in $(CLEANEXTS); do rm -f *.$$file; done


