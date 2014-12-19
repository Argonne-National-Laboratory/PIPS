# Build instructions for Sierra:

1. Assumptions:
- we're going to use GNU compilers for now. In theory, it's
  straightforward to swap out GCC for {Intel, PGI}

2. Load MPI modules from MVAPICH. Only MPICH works with PIPS. For some reason,
   the versions of OpenMPI on LC do not cooperate with PIPS (but these could
   work, perhaps, given enough effort).

    # Loading MPICH for Infiniband (i.e., MVAPICH)
    use mvapich2-gnu

The `use mvapich2-gnu` command loads the proper compiler wrapper
scripts and environment variables. This step is important because
build scripts use the `mpicxx` command to build MPI-based C++
libraries.

3. Build shared libraries
   a) Build Cbc (within Cbc directory, `./configure --enable-static && make install`
   b) Build MA57/MA27 (within MA57/MA27) directory, just `make`.
   c) Build METIS (within Metis directory), 

4. Remove any existing `build` directory (if present) and make a
   `build` directory inside the root PIPS directory, `cd` into it, and
   then call CMake. *Assuming you are in the root directory of the
   PIPS repo* (this caveat is very important!), the following compound
   command will work:
   
    # 1. Remove old build directory
	# 2. Create new build directory
	# 3. Change to new build directory
	# 4. Call CMake
	rm -rf build && mkdir build && cd build && cmake \
	-DBoost_NO_BOOST_CMAKE=TRUE \
	-DBOOST_ROOT:PATHNAME=/usr/gapps/ppp/boost-nompi-1.57.0 \
	-DBoost_LIBRARY_DIRS:FILEPATH=/usr/gapps/ppp/boost-nompi-1.57.0/lib \
	-DBoost_INCLUDE_DIRS:FILEPATH=/usr/gapps/ppp/boost-nompi-1.57.0/include \
	..

 A few remarks are in order:
 - you must be a member of group `ppp`; the easiest way to check is to run
   the `groups` command at the command line.
 - these instructions assume Boost 1.57.0 is installed in
   `/usr/gapps/ppp/boost-nompi-1.57.0`; access requires `ppp` group
   membership
 - PIPS requires Boost 1.56.0 or later; the LC has installed Boost 1.57.0,
   but it currently does not work with PIPS
 - All of the complexity in the CMake command stems from disabling the
   CMake default path checks for Boost installations, followed by
   manually specifying paths for Boost.
 - The instructions above are for an "out-of-source build", _i.e._, a
   build in a directory separate from the source directories. CMake
   recommends this practice because it makes rebuilding easier: just
   delete the build directory, make a new build directory, and rerun
   cmake.  CMake caches information in many different places, and an
   in-source build would require deleting all of the various
   information cached by CMake, which is often a pain to debug. For
   that reason, please build PIPS out-of-source.

5. Build executables with `make`.

# Running PIPS:

All of the PIPS executables should be run either at the command line
(in an `mxterm`) with `srun`, or in a Moab submission script submitted
using `msub`. For details on any of these commands, see the [LC tutorials](https://computing.llnl.gov/?set=training&page=index).

What follows are some notes on the various executables built in the
`build` directory:

## will run with just "raw input file" (special internal text-based format):
- cbcExtensive: solve using Cbc/Clp
- clpFromRaw: solve using Clp
    + writes to output files
	+ basis info contained in output files
	+ optional basis file as input
	+ looks like basis file used for warm starts
- cpmRaw: solve using PIPS-S
- lagrangeCombinedScenRedRootNode: solve using Clp
- lagrangeRootNode: solve using Clp
- pipssFromRaw: solve using PIPS-S
- testSols: check solution (using Clp?); also needs # of iterations

## will not run with just "raw input file"/not applicable:
- cbcExtensiveSMPS (requires SMPS input file)
- clpBootstrap (requires basis input file, basis output file)
- clpSMPS (requires SMPS input file)
- cutsTest (requires LP relaxation file)
- lagrangeCombinedSMPS (requires SMPS input file, # per subproblem)
- lagrangeRootNodeSMPS
- ooqpFromRaw: solves QPs
- pipsipmBatchFromRaw:
- pipsipmBatchFromRaw_schur:
- pipsipmFromRaw:
- pipsipmFromRaw_comm2_schur:
- pipsipmFromRaw_schur:
- pipssBootstrap: (requires basis input, basis output)
- qpgen-sparse-ma57-gondzio: some QP solver
- qpgen-sparse-ma57-mehrotra: some QP solver
- qpgen-sparse-mehrotra: some QP solver
- ucCoefficientDiving: requires LP relaxation basis
- ucRollingCbcExtensive: requires # scenarios, time horizon...no clue here
- ucRollingOOQPExtensive: solves some QP
- ucRollingTest: no idea what input is
- ucRootNodeClp: requires LP relaxation basis
- ucRootNodePIPSS: requires LP relaxation basis
- ucTestSol: requires fractional solution file, rounding cutoff

