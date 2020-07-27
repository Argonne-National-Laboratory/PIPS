/* PIPS-IPM                                                             
 * Authors: Cosmin G. Petra, Miles Lubin
 * (C) 2012 Argonne National Laboratory, see documentation for copyright
 */
#include <stdlib.h>
#include <iostream>
#include <iomanip>
using namespace std;
#include <unistd.h>

#include "pipschecks.h"
#include "PardisoSchurSolver.h"
#include "SparseStorage.h"
#include "SparseSymMatrix.h"
#include "SimpleVector.h"
#include "SimpleVectorHandle.h"
#include "DenseGenMatrix.h"
#include "pipsdef.h"
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include "pipsport.h"

#include "Ma57Solver.h"

#include "mpi.h"
#include "omp.h"

#ifdef HAVE_GETRUSAGE
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#endif

#ifdef WITH_MKL_PARDISO
#include "mkl_pardiso.h"
#include "mkl_types.h"
#endif

extern int gOoqpPrintLevel;
extern double g_iterNumber;
extern int gOuterBiCGIter;
#ifdef STOCH_TESTING
extern double g_scenNum;
static int rhsCount=0;
#endif
using namespace std;

#ifndef WITH_MKL_PARDISO
extern "C" void pardisoinit (void   *, int    *,   int *, int *, double *, int *);
extern "C" void pardiso     (void   *, int    *,   int *, int *,    int *, int *, 
                  double *, int    *,    int *, int *,   int *, int *,
                     int *, double *, double *, int *, double *);
extern "C" void pardiso_get_schur(void*, int*, int*, int*, double*, int*, int*);

extern "C" void pardiso_chkmatrix  (int *, int *, double *, int *, int *, int *);
extern "C" void pardiso_chkvec     (int *, int *, double *, int *);
extern "C" void pardiso_printstats (int *, int *, double *, int *, int *, int *,
                           double *, int *);
#endif

int dumpAugMatrix(int n, int nnz, int nSys, //size, nnz and size of the (1,1) block
		  double* eltsA, 
		  int* rowptr, 
		  int* colidx, 
		  const char* fname=nullptr);
int dumpSysMatrix(SparseSymMatrix* Msys,
                  const char* fname=nullptr);
int dumpRhs(SimpleVector& v);

const static int pardiso_verbosity = 0;


#ifdef TIMING_FLOPS
extern "C" {
    void HPM_Init(void);
    void HPM_Start(char *);
    void HPM_Stop(char *);
    void HPM_Print(void);
    void HPM_Print_Flops(void);
    void HPM_Print_Flops_Agg(void);
    void HPM_Terminate(char*);
}
#endif

#define SHRINK_SC

PardisoSchurSolver::PardisoSchurSolver( SparseSymMatrix * sgm )
{
  Msys = sgm;
  n = -1;
  nnz = -1;
  nSC = -1;
  rowptrAug = nullptr;
  colidxAug = nullptr;
  eltsAug = nullptr;
  shrinked2orgSC = nullptr;
  nvec = nullptr;
  nvec2 = nullptr;
  nvec_size = -1;
  // - we do not have the augmented system yet; most of initialization done during the
  // first solve call
  useSparseRhs = false;
  first = true; firstSolve = true;

#ifndef WITH_MKL_PARDISO
  num_threads = PIPSgetnOMPthreads();
#endif

  maxfct = 1;
  mnum = 1;
  phase = 0;
  msglvl = pardiso_verbosity;
  solver = 0;
  mtype = -2;
  useSparseRhs = false;

  int myRank; MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

  // todo proper parameter
  char* var = getenv("PARDISO_SPARSE_RHS_LEAF");
  assert(!useSparseRhs);
  if( var != nullptr )
  {
     int use;
     sscanf(var, "%d", &use);
     if( use == 1 )
        useSparseRhs = true;
  }

  // todo proper parameter
  var = getenv("PARDISO_SYMB_INTERVAL");
  symbFactorInterval = symbFactorIntervalDefault;

  if( var != nullptr )
  {
     int interval;
     sscanf(var, "%d", &interval);
     if( interval >= 1 )
        symbFactorInterval = interval;
  }

  // todo proper parameter
  var = getenv("PARDISO_PIVOT_PERTURBATION");
  pivotPerturbationExp = pivotPerturbationExpDefault;

  if( var != nullptr )
  {
     int exp;
     sscanf(var, "%d", &exp);
     if( exp >= 1 )
        pivotPerturbationExp = exp;
  }

  // todo proper parameter
  var = getenv("PARDISO_NITERATIVE_REFINS");
  nIterativeRefins = nIterativeRefinsDefault;

  if( var != nullptr )
  {
     int n;
     sscanf(var, "%d", &n);
     assert(n >= 0);

     if( n >= 0 )
        nIterativeRefins = n;
  }

  parallelForwardBackward = parallelForwardBackwardDefault,

  // todo proper parameter
  var = getenv("PARDISO_PARALLEL_SOLVE");
  if( var != nullptr )
  {
     int n;
     sscanf(var, "%d", &n);
     if( n == 0 )
        parallelForwardBackward = false;
     else if( n == 1 )
        parallelForwardBackward = true;
  }

  factorizationTwoLevel = factorizationTwoLevelDefault,

  // todo proper parameter
  var = getenv("PARDISO_FACTORIZE_TWOLEVEL");
  if( var != nullptr )
  {
     int n;
     sscanf(var, "%d", &n);
     if( n == 0 )
        factorizationTwoLevel = false;
     else if( n == 1 )
        factorizationTwoLevel = true;
  }


  if( myRank == 0 )
  {
     printf(" using pivot perturbation 10^-%d \n", pivotPerturbationExp);

     printf(" using maximum of %d iterative refinements  \n", nIterativeRefins);

     if( useSparseRhs )
        printf(" using PARDISO_SPARSE_RHS_LEAF \n");
     else
        printf(" NOT using PARDISO_SPARSE_RHS_LEAF \n");
  }

  nrhs = 1;
}

PardisoSchur32Solver::PardisoSchur32Solver( SparseSymMatrix * sgm )
    : PardisoSchurSolver( sgm )
{}

void PardisoSchurSolver::firstCall()
{
   iparm[0] = 0; /* make init set iparm to default values */

#ifndef WITH_MKL_PARDISO
   int error = 0;
   pardisoinit(pt, &mtype, &solver, iparm, dparm, &error);

   if( error != 0 )
   {
      cout << "PardisoSolver ERROR during pardisoinit:" << error << "." << endl;
      exit(1);
   }
#else
   pardisoinit(pt, &mtype, iparm);
#endif

   setIparm(iparm);
} 

void PardisoSchur32Solver::firstCall()
{
#ifndef WITH_MKL_PARDISO
   int error = 0;
   pardisoinit(pt, &mtype, &solver, iparm, dparm, &error);

   if( error != 0 )
   {
      cout << "PardisoSolver ERROR during pardisoinit:" << error << "." << endl;
      assert(false);
   }
#else
   pardisoinit(pt, &mtype, iparm);
#endif

   setIparm(iparm);

#ifndef WITH_MKL_PARDISO
   iparm[28] = 1; //32-bit factorization
#else
   iparm[27] = 1; //input must be in single precision as well
#endif

} 

// this function is called only once and creates the augmented system
void PardisoSchurSolver::firstSolveCall(SparseGenMatrix& R, 
					SparseGenMatrix& A,
					SparseGenMatrix& C,
					SparseGenMatrix& F,
					SparseGenMatrix& G,
					int nSC0)
{
  int nR,nA,nC,nF,nG,nx;
  nnz=0;

  F.getSize(nF,nx); nnz += F.numberOfNonZeros();
  G.getSize(nG,nx); nnz += G.numberOfNonZeros();
  R.getSize(nR,nx); nnz += R.numberOfNonZeros();
  A.getSize(nA,nx); nnz += A.numberOfNonZeros();
  C.getSize(nC,nx); nnz += C.numberOfNonZeros();
  const int Msize = Msys->size();

  if( nR == 0 )
     assert(R.numberOfNonZeros() == 0);
  if( nA == 0 )
     assert(A.numberOfNonZeros() == 0);
  if( nC == 0 )
     assert(C.numberOfNonZeros() == 0);
  if( nF == 0 )
     assert(F.numberOfNonZeros() == 0);
  if( nG == 0 )
     assert(G.numberOfNonZeros() == 0);

  assert( F.getStorageRef().isValid() );
  assert( G.getStorageRef().isValid() );
  assert( R.getStorageRef().isValid() );
  assert( A.getStorageRef().isValid() );
  assert( C.getStorageRef().isValid() );
  assert( Msys->getStorageRef().isValid() );

  // todo not implemented yet
  assert(R.numberOfNonZeros() == 0);

  if( nF > 0 || nG > 0 )
    nSC = nSC0;
  else
    nSC = nx;

#ifdef TIMING
  cout << "firstSolveCall: nR=" << nR << " nA=" << nA << " nC=" << nC << " nF=" << nF << " nG=" << nG << " nSC=" << nSC << " sizeKi=" << (nR+nA+nC)<< endl
      << " nnzR=" << R.numberOfNonZeros()
      << " nnzA=" << A.numberOfNonZeros()
      << " nnzC=" << C.numberOfNonZeros()
      << " nnzF=" << F.numberOfNonZeros()
      << " nnzG=" << C.numberOfNonZeros() << endl;
#endif

  n = nR+nA+nC+nSC;

  assert( Msize == nR+nA+nC );
  nnz += Msys->numberOfNonZeros();
  nnz += nSC; //space for the 0 diagonal of 2x2 block

  // the lower triangular part of the augmented system in row-major
  SparseSymMatrix augSys(n, nnz);
  
  // pointer for augmented system
  int* krowAug = augSys.getStorageRef().krowM;
  int* jcolAug = augSys.getStorageRef().jcolM;
  double* MAug = augSys.getStorageRef().M;

  //
  //put (1,1) block in the augmented system
  //
  memcpy(krowAug, Msys->getStorageRef().krowM, sizeof(int)*Msize);
  memcpy(jcolAug, Msys->getStorageRef().jcolM, sizeof(int)*Msys->numberOfNonZeros());
  memcpy(MAug,    Msys->getStorageRef().M,     sizeof(double)*Msys->numberOfNonZeros());


  int nnzIt = Msys->numberOfNonZeros();
  //
  //put A and C block in the augmented system as At and Ct in the lower triangular part
  //

  if( nA > 0 || nC > 0 || nF > 0 || nG > 0 )
  {
    const bool putA = A.numberOfNonZeros() > 0;
    const bool putC = C.numberOfNonZeros() > 0;

    // putA = TRUE => nA > 0
    assert(putA == (putA && (nA > 0)));
    assert(putC == (putC && (nC > 0)));

    // initialize variables for At
    SparseGenMatrix At(putA ? nx : 0, putA ? nA : 0, putA ? A.numberOfNonZeros() : 0);
    int* krowAt=At.getStorageRef().krowM;
    int* jcolAt=At.getStorageRef().jcolM;
    double *MAt=At.getStorageRef().M;

    if( putA )
       A.getStorageRef().transpose(krowAt, jcolAt, MAt);

    const int colShiftA=nR;

    // initialize variables for Ct
    SparseGenMatrix Ct(putC ? nx : 0, putC ? nC : 0, putC ? C.numberOfNonZeros() : 0);
    int* krowCt=Ct.getStorageRef().krowM;
    int* jcolCt=Ct.getStorageRef().jcolM;
    double *MCt=Ct.getStorageRef().M;

    if( putC )
       C.getStorageRef().transpose(krowCt, jcolCt, MCt);

    const int colShiftC=nR+nA;

    int row=Msize;
    for( ; row < Msize + nx; row++ ) {
      krowAug[row]=nnzIt;

      if( putA ) {
         for(int c=krowAt[row-Msize]; c< krowAt[row-Msize+1]; c++) {
            const int j=jcolAt[c];
            jcolAug[nnzIt]=j+colShiftA;
            MAug[nnzIt]   =MAt[c];
            nnzIt++;
         }
      }

      if( putC )
      {
         for(int c=krowCt[row-Msize]; c< krowCt[row-Msize+1]; c++) {
            const int j=jcolCt[c];
            jcolAug[nnzIt]=j+colShiftC;
            MAug[nnzIt]   =MCt[c];
            nnzIt++;
         }
      }
      //add the zero from the diagonal
      jcolAug[nnzIt]=row;
      MAug[nnzIt]=0.0;
      nnzIt++;

    }
    krowAug[row]=nnzIt;
  }
  assert( nnzIt = Msys->numberOfNonZeros() + A.numberOfNonZeros() + C.numberOfNonZeros() + nx );

  //
  // add linking constraint matrices F and G
  //
  if( nF > 0 || nG > 0 )
  {
     // put diagonal in zero block
     int row = Msize + nx;
     assert(row <= n - nF - nG);
     for( ; row < n - nF - nG; row++ ) {
        krowAug[row]=nnzIt;
        jcolAug[nnzIt] = row;
        MAug[nnzIt] = 0.0;
        nnzIt++;
     }

     // are there linking equality constraints?
     if( nF > 0 )
     {
        // put F in the lower triangular part below R (and below 0 block)

        int* krowF=F.getStorageRef().krowM;
        int* jcolF=F.getStorageRef().jcolM;
        double *MF=F.getStorageRef().M;

        const bool putF = F.numberOfNonZeros() > 0;

        const int row0 = n - nF - nG;
        int row = row0;
        assert(row0 >= Msize+nx);

        for( ; row < n - nG; row++ ) {
          krowAug[row]=nnzIt;

          if( putF ) {
             for(int c=krowF[row-row0]; c< krowF[row-row0+1]; c++) {
                jcolAug[nnzIt]=jcolF[c];
                MAug[nnzIt]   =MF[c];
                nnzIt++;
             }
          }
          //add the zero from the diagonal
          jcolAug[nnzIt]=row;
          MAug[nnzIt]=0.0;
          nnzIt++;
        }
        krowAug[row]=nnzIt;
     }

     // are there linking equality constraints?
     if( nG > 0 )
     {
        // put G in the lower triangular part below F

        int* krowG=G.getStorageRef().krowM;
        int* jcolG=G.getStorageRef().jcolM;
        double *MG=G.getStorageRef().M;

        const bool putG = G.numberOfNonZeros() > 0;

        const int row0 = n - nG;
        int row = row0;

        for( ; row < n; row++ ) {
          krowAug[row]=nnzIt;

          if( putG ) {
             for(int c=krowG[row-row0]; c< krowG[row-row0+1]; c++) {
                jcolAug[nnzIt]=jcolG[c];
                MAug[nnzIt]   =MG[c];
                nnzIt++;
             }
          }
          //add the zero from the diagonal
          jcolAug[nnzIt]=row;
          MAug[nnzIt]=0.0;
          nnzIt++;
        }
        krowAug[row]=nnzIt;
     }
  }
  assert( nnzIt = Msys->numberOfNonZeros() + A.numberOfNonZeros() + C.numberOfNonZeros() + F.numberOfNonZeros() + G.numberOfNonZeros() + nSC );

#ifdef SHRINK_SC

  // remove empty or zero rows from augSys and check that in Msys none are removed!
  int* shrinked2orgAug = nullptr;

  augSys.deleteZeroRowsCols(shrinked2orgAug);

  n = augSys.size();

#ifndef NDEBUG
  for( int i = 0; i < Msize; i++ )
  {
     if( i != shrinked2orgAug[i] )
     {
        std::cout << "zero row in (1,1) block of Schur complement!" << std::endl;
        std::cout << "i=" << i << " shrinked2orgAug[i]=" << shrinked2orgAug[i] << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
     }

     if( Msys->getStorageRef().krowM[i] == Msys->getStorageRef().krowM[i + 1] )
     {
        std::cout << "(2) zero row in (1,1) block of Schur complement!"<< std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
     }
   }
#endif

  nSC = augSys.size() - Msize;

  assert(shrinked2orgSC == nullptr);
  shrinked2orgSC = new int[nSC];

  // shrinked2orgSC should only map Schur complement part of augmented saddle-point system
  for( int i = 0; i < nSC; i++ )
  {
     shrinked2orgSC[i] = shrinked2orgAug[i + Msize] - Msize;
     assert(shrinked2orgSC[i] >= i);
  }

  delete[] shrinked2orgAug;
#endif

  nnz = augSys.numberOfNonZeros();

  // we need to transpose to get the augmented system in the row-major upper triangular format of PARDISO
  rowptrAug = new int[n+1];
  colidxAug = new int[nnz];
  eltsAug   = new double[nnz];

  augSys.getStorageRef().transpose(rowptrAug,colidxAug,eltsAug);

  assert(rowptrAug[n] == nnz);

  //save the indices for diagonal entries for a streamlined later update
  int* krowMsys = Msys->getStorageRef().krowM;
  int* jcolMsys = Msys->getStorageRef().jcolM;

  for(int r=0; r<Msize; r++) {

    // Msys - find the index in jcol for the diagonal (r,r)
    int idxDiagMsys=-1;
    for(int idx=krowMsys[r]; idx<krowMsys[r+1]; idx++)
      if(jcolMsys[idx]==r) {idxDiagMsys=idx; break;}

    assert(idxDiagMsys>=0);

    // aug  - find the index in jcol for the diagonal (r,r)
    int idxDiagAug=-1;
    for(int idx=rowptrAug[r]; idx<rowptrAug[r+1]; idx++)
      if(colidxAug[idx]==r) {idxDiagAug=idx; break;}

    assert(idxDiagAug>=0);

    diagMap.insert( pair<int,int>(idxDiagMsys,idxDiagAug) );
  }

  //convert to Fortran indexing
  for(int it=0; it<n+1; it++)   rowptrAug[it]++;
  for(int it=0; it<nnz; it++)   colidxAug[it]++;

  //allocate temp vector(s)
  nvec=new double[n];
  nvec2=new double[n];
  nvec_size = n;
} 

void PardisoSchurSolver::setIparm(int* iparm){

   /* common parameters */
   iparm[1] = 2; // 2 and 3 are for METIS - 2 is METIS 4.1, 3 is METIS 5.1, 0 for min degree ordering

   /* NOTE: if iparm[9] is less than 13 mkl_pardiso will not consistently produce the same schur complement as the other pardiso (on some examples)
    * this might not be an issue should be kept in mind though
    */
   iparm[30] = 0; // do not specify sparse rhs at this point ! MKL_PARDISO can either set iparm[35] or iparm[30]

#ifndef WITH_MKL_PARDISO
   iparm[2] = PIPSgetnOMPthreads();

   iparm[7] = nIterativeRefins; // max number of iterative refinement steps
   iparm[9] = pivotPerturbationExp; // pivot perturbation 10^{-x} * |A|_\{\inf}
   iparm[10] = 1; // default, scaling for IPM KKT used with either mtype=11/13 or mtype=-2/-4/6 and iparm[12]=1
   iparm[12] = 2;// 0 disable matching, 1 enable matching, no other settings

   if( factorizationTwoLevel )
      iparm[23] = 1; // parallel Numerical Factorization (0=used in the last years, 1=two-level scheduling)
   else
      iparm[23] = 0;

   if( parallelForwardBackward )
      iparm[24] = 1; // parallelization for the forward and backward solve. 0=sequential, 1=parallel solve.
   else
      iparm[24] = 0;

#else
   /* From INTEL (instead of iparm[2] which is not defined there):
    *  You can control the parallel execution of the solver by explicitly setting the MKL_NUM_THREADS environment variable.
    *  If fewer OpenMP threads are available than specified, the execution may slow down instead of speeding up.
    *  If MKL_NUM_THREADS is not defined, then the solver uses all available processors.
    */
   iparm[7] = 0; // MKL_PARDISO runs into troubles otherwise
   iparm[9] = 13; // MKL_PARDISO need this in order to compute same schur decomposition as schenk pardiso
   iparm[10] = 0; // scaling for IPM KKT; used with IPARM(13)=1 or 2
   iparm[12] = 0; // improved accuracy for IPM KKT; used with IPARM(11)=1; use 2 for advanced matchings and higher accuracy.

   /* NOTE: requires iparm[23] = 1 which in return requires iparm[10] = iparm[12] = 0 */
   /* even though the documentation does not tell so setting iparm[23] = 10 is highly unstable and might result in segmentation faults */
   iparm[35] = -2; // compute the schur complement

   /* mkl_pardiso has no chkmatrix method - instead one can set iparm[26] */
   #ifndef NDEBUG
   iparm[26] = 1;
   #endif
#endif

}

bool PardisoSchurSolver::iparmUnchanged(){

   /* put all Parameters that should stay be checked against init into this array */
   static const int check_iparm[] = { 1, 30, 10, 12, 23, 24 };

   bool unchanged = true;
   bool print = false;

   int iparm_compare[64];
   setIparm(iparm_compare);


   vector<int> to_compare(check_iparm, check_iparm + sizeof(check_iparm) / sizeof(check_iparm[0]) );

   for(int i = 0; i < 64; ++i)
   {
      // if entry should be compared
      if(std::find(to_compare.begin(), to_compare.end(), i) != to_compare.end())
      {
         if(iparm[i] != iparm_compare[i])
         {
            if(print)
               std::cout << "ERROR - PardisoSolver: elements in iparm changed at " << i << ": "
                  << iparm[i] << " != " << iparm_compare[i] << "(new)" << std::endl;
            unchanged = false;
         }
      }
   }
   return unchanged;
}

void PardisoSchurSolver::diagonalChanged( int /* idiag */, int /* extent */ )
{
  this->matrixChanged();
}

void PardisoSchurSolver::matrixChanged()
{
  if(first) { firstCall(); first=false;}

  // we don't have the right hand-size, therefore we can't (re)factorize 
  // the augmented system at this point.

  //dumpSysMatrix(Msys);
}
 
void PardisoSchurSolver::schur_solve(SparseGenMatrix& R, 
				     SparseGenMatrix& A,
				     SparseGenMatrix& C,
				     SparseGenMatrix& F,
				     SparseGenMatrix& G,
				     DenseSymMatrix& SC0)
{
  int* rowptrSC = nullptr;
  int* colidxSC = nullptr;
  double* eltsSC = nullptr;

  computeSC(SC0.size(), R, A, C, F, G, rowptrSC, colidxSC, eltsSC);

#ifdef SHRINK_SC
  for( int r = 0; r < nSC; r++ )
  {
     const int r_org = shrinked2orgSC[r];
     assert(r_org >= r && r_org < SC0.size());

     for( int ci = rowptrSC[r]; ci < rowptrSC[r + 1]; ci++ )
     {
        const int c = colidxSC[ci];
        const int c_org = shrinked2orgSC[c];

        assert(c >= r);
        assert(c_org >= c && c_org < SC0.size());

        SC0[c_org][r_org] += eltsSC[ci];

        // NOTE: we only save half of the matrix, so we don't need
        // if( r_org != c_org ) SC0[r_org][c_org] += eltsSC[ci];
     }
  }
#else
  for( int r = 0; r < nSC; r++ )
  {
     for( int ci = rowptrSC[r]; ci < rowptrSC[r + 1]; ci++ )
     {
        const int c = colidxSC[ci];
        assert(c >= r);

        SC0[c][r] += eltsSC[ci];

     }
  }
#endif

  delete[] rowptrSC; delete[] colidxSC; delete[] eltsSC;
}


void PardisoSchurSolver::schur_solve_sparse(SparseGenMatrix& R,
                 SparseGenMatrix& A,
                 SparseGenMatrix& C,
                 SparseGenMatrix& F,
                 SparseGenMatrix& G,
                 SparseSymMatrix& SC0)
{
  int* rowptrSC = nullptr;
  int* colidxSC = nullptr;
  double* eltsSC = nullptr;

  computeSC(SC0.size(), R, A, C, F, G, rowptrSC, colidxSC, eltsSC);

  int* rowptrBase = SC0.krowM();
  int* colidxBase = SC0.jcolM();
  double* eltsBase = SC0.M();

#ifdef SHRINK_SC
  assert(SC0.size() >= nSC);

  // add to summed Schur complement
  for( int r = 0; r < nSC; r++ )
  {
     const int r_org = shrinked2orgSC[r];
     assert(r_org >= r && r_org < SC0.size());

     int cbase = rowptrBase[r_org];

     // catch empty diagonal
     if( rowptrSC[r + 1] - rowptrSC[r] == 1 && eltsSC[rowptrSC[r]] == 0.0 )
        continue;

     for( int j = rowptrSC[r]; j < rowptrSC[r + 1]; j++ )
     {
        const int c_org = shrinked2orgSC[colidxSC[j]];
        assert(c_org >= colidxSC[j] && c_org < SC0.size());

        while( colidxBase[cbase] != c_org )
        {
           cbase++;
           assert(cbase < rowptrBase[r_org + 1]);
        }

        eltsBase[cbase] += eltsSC[j];
     }
  }
#else
  assert(SC0.size() == nSC);


  // add to summed Schur complement
  for( int r = 0; r < nSC; r++ )
  {
     int cbase = rowptrBase[r];

     // catch empty diagonal
     if( rowptrSC[r + 1] - rowptrSC[r] == 1 && eltsSC[rowptrSC[r]] == 0.0 )
        continue;

     for( int j = rowptrSC[r]; j < rowptrSC[r + 1]; j++ )
     {
        const int c = colidxSC[j];

        while( colidxBase[cbase] != c )
        {
           cbase++;
           assert(cbase < rowptrBase[r + 1]);
        }

        eltsBase[cbase] += eltsSC[j];
     }
  }
#endif

  delete[] rowptrSC; delete[] colidxSC; delete[] eltsSC;
}


void PardisoSchurSolver::computeSC(int nSCO,
/*const*/SparseGenMatrix& R,
/*const*/SparseGenMatrix& A,
/*const*/SparseGenMatrix& C,
/*const*/SparseGenMatrix& F,
/*const*/SparseGenMatrix& G, int*& rowptrSC, int*& colidxSC, double*& eltsSC)
{
   assert(!rowptrSC && !colidxSC && !eltsSC);

   bool doSymbFact = false;
   if( firstSolve )
   {
      firstSolveCall(R, A, C, F, G, nSCO);
      firstSolve = false;
      doSymbFact = true;
   }
   else
   {
      //update diagonal entries in the PARDISO aug sys
      const double* eltsMsys = Msys->getStorageRef().M;
      map<int, int>::iterator it;

#if 0
      double max = -1e20;
      double min = 1e20;
      double minAbs = 1e20;

      for(it=diagMap.begin(); it!=diagMap.end(); it++)
      {
         const double elem = eltsMsys[it->first];
         if(elem > max)
         max = elem;
         if(elem < min)
         min = elem;
         if(std::fabs(elem) < minAbs && elem > 0.0 )
         minAbs = std::fabs(elem);

      }
      std::cout << "local Schur diag: min/max/minAbs  " << min << " " << max << " " << minAbs << std::endl;
#endif

      for( it = diagMap.begin(); it != diagMap.end(); it++ )
         eltsAug[it->second] = eltsMsys[it->first];
   }

   // call PARDISO
   int error = 0;

#ifndef NDEBUG
#ifndef WITH_MKL_PARDISO
   pardiso_chkmatrix(&mtype, &n, eltsAug, rowptrAug, colidxAug, &error);
   if( error != 0 )
   {
      cout << "PARDISO matrix error " << error << endl;
      exit(1);
   }
#endif
#endif

   const int nIter = (int) g_iterNumber;
   int myRank; MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

   if( (nIter % symbFactorInterval) == 0 )
   {
      doSymbFact = true;
      if( myRank == 0 )
         printf("PardisoSchur: starting symbolic analysis ... ");
   }

   int phase = 22; // numerical factorization
   if( doSymbFact )
      phase = 12; // numerical factorization & symb analysis

   assert(iparmUnchanged());

   /* compute schur complement */
#ifndef WITH_MKL_PARDISO
   iparm[37] = nSC; // compute Schur-complement

   #ifdef TIMING_FLOPS
   HPM_Start("PARDISOFact");
   #endif

   pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, eltsAug, rowptrAug,
         colidxAug, nullptr, &nrhs, iparm, &msglvl, nullptr, nullptr, &error ,dparm);

 #ifdef TIMING_FLOPS
   HPM_Stop("PARDISOFact");
 #endif

   if( doSymbFact && myRank == 0 )
      printf("finished \n");

   const int nnzSC = iparm[38];

#ifdef TIMING
   if(1001*(myRank/1001)==myRank)
   printf("rank %d perturbPiv %d peakmem %d\n", myRank, iparm[13], iparm[14]); // same for MKL and pardiso
   //cout << "NNZ(SCHUR) " << nnzSC << "    SPARSITY " << nnzSC/(1.0*nSC*nSC) << endl;
 #endif

   if( error != 0 )
   {
      printf("PardisoSolver - ERROR during factorization: %d. Phase param=%d\n",
            error, phase);
      exit(1);
   }
   rowptrSC = new int[nSC + 1];
   colidxSC = new int[nnzSC];
   eltsSC = new double[nnzSC];

   pardiso_get_schur(pt, &maxfct, &mnum, &mtype, eltsSC, rowptrSC, colidxSC);

   //convert back to C/C++ indexing
   for( int it = 0; it < nSC + 1; it++ )
      rowptrSC[it]--;
   for( int it = 0; it < nnzSC; it++ )
      colidxSC[it]--;
#else
   assert(iparm[35] == -2);
   assert(iparm[23] == 1);

   /* perm array has to be defined - it is of size of augmented system (dense) and inicates rows/colums we want to have in the schur complement with 1
    * rest set to 0
    * if specified like {1, 0, 0, 1, 0, 0} mkl_pardiso will first reorder row 0 and 3 to the bottom (hopefully stable sort) and then use the last two rows for computing
    * the schur complement
    * setting all entries to 1 will return the input matrix as schur complement
    */
   int perm[n] = {0};
   for(int i = n - nSC; i < n; ++i)
      perm[i] = 1;

   /* reordering and symbolic factorization */
   phase = 11;

   pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, eltsAug, rowptrAug,
   	   	   colidxAug, perm, &nrhs, iparm, &msglvl, nullptr, nullptr, &error);

   /* iparm[35] should now contain the number of non-zero entries in the Schur-complement */
   /* if it does not most likely the general iparm setup is faulty */
   assert(error == 0); assert(iparm[35] >= 0);

   /* preallocation of schur matrix arrays */
   int nnzSC = iparm[35];

   /* get deleted inside of destructor of schur_transposed */
   int* rowptrSCtransp = new int[nSC + 1];
   int* colidxSCtransp = new int[nnzSC];
   double* eltsSCtransp = new double[nnzSC];

   phase = 22;

   /* reset iparm[35] to -2 necessary! */
   iparm[35] = -2;
   int step = 1;

   // mkl pardiso returns arrays with zero-based index via export
   pardiso_export(pt, eltsSCtransp, rowptrSCtransp, colidxSCtransp, &step, iparm, &error);

   #ifdef TIMING_FLOPS
   HPM_Start("PARDISOFact");
   #endif

   /* factorization call */
   pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, eltsAug, rowptrAug,
         colidxAug, perm, &nrhs, iparm, &msglvl, nullptr, nullptr, &error);

   #ifdef TIMING_FLOPS
   HPM_Stop("PARDISOFact");
   #endif

   #ifdef TIMING
   if(1001*(myRank/1001)==myRank)
      printf("rank %d perturbPiv %d peakmem %d\n", myRank, iparm[13], iparm[14]); // same for MKL and pardiso
   //cout << "NNZ(SCHUR) " << nnzSC << "    SPARSITY " << nnzSC/(1.0*nSC*nSC) << endl;
   #endif
   if( error != 0 )
   {
      printf("PardisoSolver - ERROR during factorization: %d. Phase param=%d\n",
            error, phase);
      exit(1);
   }

   /////////////////////////////////////////////////////
   // transpose the matrix since it is given in lower triangular form and we are using upper triangular form
   /////////////////////////////////////////////////////
   SparseStorage schur_transposed( nSC, nSC, nnzSC, rowptrSCtransp, colidxSCtransp, eltsSCtransp, true);

   rowptrSC = new int[nSC + 1];
   colidxSC = new int[nnzSC];
   eltsSC = new double[nnzSC];

   schur_transposed.transpose(rowptrSC, colidxSC, eltsSC);

#endif

   assert(subMatrixIsOrdered(rowptrSC, colidxSC, 0, nSC));
}


void PardisoSchurSolver::solve( OoqpVector& rhs_in )
{ 
  SimpleVector& rhs=dynamic_cast<SimpleVector&>(rhs_in);

  int error = 0;
  assert(iparmUnchanged());

  assert(nvec_size == n);
  double* const x_n = nvec;
  double* const rhs_n = nvec2;

  const int dim=rhs.length();
  assert(dim >= 0 && dim <= n);

  memset(x_n, 0, dim * sizeof(double));
  memcpy(rhs_n, rhs.elements(), dim * sizeof(double));

  if( n > dim )
     memset(&rhs_n[dim], 0, (n - dim) * sizeof(double));

#ifdef TIMING_FLOPS
  HPM_Start("PARDISOSolve");
#endif

  // solving phase
#ifndef WITH_MKL_PARDISO
  phase = 33; /* solve - iterative refinement */

  int* rhsSparsity = nullptr;
  if( useSparseRhs )
  {
     iparm[30] = 1; //sparse rhs
     rhsSparsity = new int[n]();

     for( int i = 0; i < dim; i++  )
        if( !PIPSisZero(rhs_n[i]) )
           rhsSparsity[i] = 1;
  }
  else
  {
     iparm[30] = 0;
  }

  pardiso (pt , &maxfct , &mnum, &mtype, &phase,
	   &n, eltsAug, rowptrAug, colidxAug,
	   rhsSparsity, &nrhs,
	   iparm , &msglvl, rhs_n, x_n, &error
	   ,dparm
  );

  iparm[30] = 0;
  delete[] rhsSparsity;

  assert(error == 0);
#else
  /* pardiso from mkl does not support same functionality as pardiso-project
   *
   * pardiso project:
   * when computing the schur complement S with factorization matrices we will get
   *
   * [A11 A12]   [L11 0] [I 0] [U11 U12]
   * [A21 A22] = [L12 I] [0 S] [0     I]
   *
   * a subsequent solve call will then only solve for A11 x1 = b1 instead of the full
   * system.
   *
   * pardiso mkl:
   * while the schur complement is the same, the factorization computed, stored and
   * used for solve calls is a full factorization. thus pardiso from intel will always
   * solve the full system
   *
   * workaround is to solve
   *
   * (phase 331)
   * [L11   0] [z1] = [b1]
   * [L12   I] [z2] = [b2]
   *
   * (phase 332)
   * [I 0] [y1]   [z1]
   * [0 S] [y2] = [z2]
   *
   * (phase 333)
   * [U11 U12] [x1]   [y1]
   * [0     I] [x2] = [y2]
   *
   */


  // forward substitution
   double* z_n = new double[nvec_size];
   assert(iparm[7] == 0);
   assert(iparm[35] = -2);

   // this is necessary for usage of stage = 331/332/333
   iparm[9] = 0;

   // HACK: keeping iparm[35] = -2 will, for some reason, not compute the correct result
   // iparm[35] will be set to -2 after stage 333
   iparm[35] = 2;

   phase = 331;
   pardiso (pt , &maxfct , &mnum, &mtype, &phase,
         &n, eltsAug, rowptrAug, colidxAug,
         nullptr, &nrhs, iparm, &msglvl, rhs_n, z_n, &error);
   assert(error == 0);

   // diagonal substitution
   phase = 332;
   double* y_n = new double[nvec_size];

   pardiso (pt , &maxfct , &mnum, &mtype, &phase,
         &n, eltsAug, rowptrAug, colidxAug,
         nullptr, &nrhs, iparm, &msglvl, z_n, y_n, &error);
   assert(error == 0);

   // backward substitution
   for(int i = dim; i < nvec_size; ++i)
      y_n[i] = 0.0;
   phase = 333;

   pardiso (pt , &maxfct , &mnum, &mtype, &phase,
         &n, eltsAug, rowptrAug, colidxAug,
         nullptr, &nrhs, iparm, &msglvl, y_n, x_n, &error);
   assert(error == 0);

   iparm[35] = -2;

   delete[] z_n;
   delete[] y_n;
#endif


#ifdef TIMING_FLOPS
  HPM_Stop("PARDISOSolve");
#endif


#ifdef PRINT_SC_RESIDUAL
  //compute residual (alternative)  
  double* tmp_resid=new double[dim];
  memcpy(tmp_resid, rhs.elements(), dim*sizeof(double));
  double mat_max=0;
  for( int i = 0; i < dim; i++ )
  {
     for( int p = rowptrAug[i]; p < rowptrAug[i + 1]; p++ )
     {
        const int j = colidxAug[p - 1] - 1;
        if( j + 1 <= dim )
        {
           //r[i] = r[i] + M(i,j)*x(j)
           tmp_resid[i] -= eltsAug[p - 1] * x_n[j];

           if( abs(eltsAug[p - 1]) > mat_max )
              mat_max = abs(eltsAug[p - 1]);

           if( j != i )
           {
              //r[j] = r[j] + M(j,i)*x(i)
              tmp_resid[j] -= eltsAug[p - 1] * x_n[i];
           }
        }
     }
  }

  double res_norm2=0.0, res_nrmInf=0, sol_inf=0.;
  for( int i = 0; i < dim; i++ )
  {
     res_norm2 += tmp_resid[i] * tmp_resid[i];
     if( res_nrmInf < fabs(tmp_resid[i]) )
        res_nrmInf = tmp_resid[i];
     if( abs(x_n[i]) > sol_inf )
        sol_inf = abs(x_n[i]);
  }
  res_norm2 = sqrt(res_norm2);

  const double rhsNorm=rhs.twonorm();
  //if(min(res_nrmInf/(mat_max*sol_inf),res_norm2/(mat_max*sol_inf))>1e-3) {
  if(min(res_nrmInf/rhsNorm,res_norm2/rhsNorm)>1e-6) {
    cout << "PardisoSchurSolve large residual - norms resid="<< res_norm2 << ":" << res_nrmInf 
	 << " rhs=" << rhsNorm << " sol="<<sol_inf << " mat="<< mat_max 
	 << " #refin.=" << iparm[6] 
	 << " rel.res.nrmInf=" << res_nrmInf/rhsNorm 
	 << " bicgiter=" << gOuterBiCGIter<< endl;

  }
  delete[] tmp_resid;
#endif
  
  memcpy(&rhs[0], x_n, dim * sizeof(double));
}

void PardisoSchur32Solver::solve( OoqpVector& rhs_in )
{ 
  SimpleVector& rhs=dynamic_cast<SimpleVector&>(rhs_in);

  /* same for mkl_pardiso and pardiso */
  int error = 0;
  int phase = 33; // solve and iterative refinement

  SimpleVector x_n(n);
  SimpleVector rhs_n(n);
  int dim=rhs.length();
  memcpy(&rhs_n[0], rhs.elements(), dim*sizeof(double));
  for(int i=dim; i<n; i++) rhs_n[i]=0.0;

  //double start = MPI_Wtime();
  pardiso (pt , &maxfct , &mnum, &mtype, &phase,
	   &n, eltsAug, rowptrAug, colidxAug, 
	   nullptr, &nrhs,
	   iparm , &msglvl, rhs_n.elements(), x_n.elements(), &error
#ifndef WITH_MKL_PARDISO
	   ,dparm
#endif
  );
  //cout << "---pardiso:" << MPI_Wtime()-start << endl;


  //compute residual (alternative)
#ifdef PRINT_SC_RESIDUAL
  double* tmp_resid=new double[dim];
  memcpy(tmp_resid, rhs.elements(), dim*sizeof(double));
  for(int i=0; i<dim; i++) {
    for(int p=rowptrAug[i]; p<rowptrAug[i+1]; p++) {
      int j=colidxAug[p-1]-1;
      if(j+1<=dim) {
	//r[i] = r[i] + M(i,j)*x(j)
	tmp_resid[i] -= eltsAug[p-1]*x_n[j];
	
	if(j!=i) {
	  //r[j] = r[j] + M(j,i)*x(i)
	  tmp_resid[j] -=  eltsAug[p-1]*x_n[i];
	}
      }
    }
  }
  double res_norm2=0.0, res_nrmInf=0; 
  for(int i=0; i<dim; i++) {
      res_norm2 += tmp_resid[i]*tmp_resid[i];
      if(res_nrmInf<fabs(tmp_resid[i]))
	 res_nrmInf=tmp_resid[i];
  }
  res_norm2 = sqrt(res_norm2);

  double rhsNorm=rhs.twonorm();
  if(res_norm2/rhsNorm>1e-3) {
    cout << "PardisoSchur32Solve::solve big residual --- rhs.nrm=" << rhsNorm 
	 << " rel.res.nrm2=" << res_norm2/rhsNorm
	 << " rel.res.nrmInf=" << res_nrmInf/rhsNorm << endl;
  }
  delete[] tmp_resid;
#endif
  

  memcpy(&rhs[0], x_n.elements(), dim*sizeof(double));
}

void PardisoSchurSolver::solve(GenMatrix& rhs_in)
{
  assert(false && "Function not supported. Use PardisoSolver for this functionality.");
}

PardisoSchurSolver::~PardisoSchurSolver()
{
  int phase = -1; /* Release internal memory . */
  msglvl = 0;
  int error = 0;
  
  pardiso (pt, &maxfct, &mnum, &mtype, &phase,
	   &n, nullptr, rowptrAug, colidxAug, nullptr, &nrhs,
	   iparm, &msglvl, nullptr, nullptr, &error
#ifndef WITH_MKL_PARDISO
	   , dparm
#endif
  );
  if ( error != 0) {
    printf ("PardisoSchurSolver - ERROR in pardiso release: %d", error ); 
  }
  
  delete[] rowptrAug;
  delete[] colidxAug;
  delete[] eltsAug;
  delete[] shrinked2orgSC;
  delete[] nvec;
  delete[] nvec2;
}

#ifdef STOCH_TESTING
int dumpAugMatrix(int n, int nnz, int nSys,
		  double* elts, int* rowptr, int* colidx, 
		  const char* fname)
{
  char filename[1024];
  int myRank; MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
  if(fname==nullptr) 
    sprintf(filename, "augMat-r%d-i%g-s%g.dat", myRank, g_iterNumber,  g_scenNum+1);
  else 
    sprintf(filename, "%s-r%d-i%g-s%g.dat", fname, myRank, g_iterNumber,  g_scenNum+1);

  cout << "saving to:" << filename << " ...";

  ofstream fd(filename);
  fd << scientific;
  fd.precision(16);

  fd << n << endl << nSys << endl << nnz << endl;
  for(int it=0; it<n+1; it++) {fd << rowptr[it] << " ";}  fd << endl;
  for(int it=0; it<nnz; it++) {fd << colidx[it] << " ";}  fd << endl;
  for(int it=0; it<nnz; it++) {fd << elts[it]   << " ";}  fd << endl;

  cout << filename << " done!" << endl;
  return 0;
}

int dumpSysMatrix(SparseSymMatrix* Msys, const char* fname)
{
  char filename[1024];
  if(fname==nullptr) 
    sprintf(filename, "Qdump-%g-s%g.dat", g_iterNumber,  g_scenNum+1);
  else 
    sprintf(filename, "%s-%g-s%g.dat", fname, g_iterNumber,  g_scenNum+1);

  cout << "saving to:" << filename << " ...";

  int n  =Msys->size();
  int nnz=Msys->numberOfNonZeros();

  // we need to transpose to get the augmented system in the row-major upper triangular format of  PARDISO 
  int* rowptr  = new int[n+1];
  int* colidx  = new int[nnz];
  double* elts = new double[nnz];
  Msys->getStorageRef().transpose(rowptr,colidx,elts);


  ofstream fd(filename);
  fd << scientific;
  fd.precision(16);

  fd << n << endl << nnz << endl;
  for(int it=0; it<n+1; it++) {fd << rowptr[it]+1 << " ";}  fd << endl;
  for(int it=0; it<nnz; it++) {fd << colidx[it]+1 << " ";}  fd << endl;
  for(int it=0; it<nnz; it++) {fd << elts[it]   << " ";}  fd << endl;

  delete[] rowptr; delete[] colidx; delete[] elts;

  cout << " Done!" << endl;
  

  return 0;
}
int dumpRhs(SimpleVector& v)
{
  rhsCount++;
  char filename[1024];
  sprintf(filename, "rhsDump-%g-%d.dat", g_iterNumber,   rhsCount);
  cout << "saving to: " << filename << " ...";

  ofstream fd(filename);
  fd << scientific;
  fd.precision(16);

  fd << v.length() << endl;
  for(int i=0; i<v.length(); i++) fd << v[i]   << " ";  
  fd << endl;

  cout << "done!" << endl;
  return 0;
}
#endif
