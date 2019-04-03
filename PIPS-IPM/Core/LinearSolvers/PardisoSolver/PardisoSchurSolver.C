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

int pardiso_stuff(int n, int nnz, int n0, int* rowptr, int* colidx, double* elts, const vector< vector<double> >& vRhs);

int dumpAugMatrix(int n, int nnz, int nSys, //size, nnz and size of the (1,1) block
		  double* eltsA, 
		  int* rowptr, 
		  int* colidx, 
		  const char* fname=NULL);
int dumpSysMatrix(SparseSymMatrix* Msys,
                  const char* fname=NULL);
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


PardisoSchurSolver::PardisoSchurSolver( SparseSymMatrix * sgm )
{
  Msys = sgm;
  n = -1;
  nnz = -1;
  nSC = -1;
  rowptrAug = NULL;
  colidxAug = NULL;
  eltsAug = NULL;

  nvec = NULL;
  nvec2 = NULL;
  nvec_size = -1;
  // - we do not have the augmented system yet; most of initialization done during the
  // first solve call

  first = true; firstSolve = true;
#ifndef WITH_MKL_PARDISO
  num_threads = PIPSgetnOMPthreads();
#endif
}

PardisoSchur32Solver::PardisoSchur32Solver( SparseSymMatrix * sgm )
    : PardisoSchurSolver( sgm )
{}

void PardisoSchurSolver::firstCall()
{
   int mtype = -2;
#ifndef WITH_MKL_PARDISO
   int solver = 0, error;
   pardisoinit(pt, &mtype, &solver, iparm, dparm, &error);

   if( error != 0 )
   {
      cout << "PardisoSolver ERROR during pardisoinit:" << error << "." << endl;
      exit(1);
   }
#else
   /* enable matrix checker (default disabled) - mkl pardiso does not have chkmatrix */
   //iparm[26] = 1;
   pardisoinit(pt, &mtype, iparm);
#endif
} 

void PardisoSchur32Solver::firstCall()
{
   int mtype = -2;
#ifndef WITH_MKL_PARDISO
   int solver = 0, error;
   pardisoinit(pt, &mtype, &solver, iparm, dparm, &error);

   if( error != 0 )
   {
      cout << "PardisoSolver ERROR during pardisoinit:" << error << "." << endl;
      assert(false);
   }
#else
   /* enable matrix checker (default disabled) - mkl pardiso does not have chkmatrix */
   //iparm[26] = 1;
   pardisoinit(pt, &mtype, iparm);
#endif

#ifndef WITH_MKL_PARDISO
   iparm[28]=1; //32-bit factorization
#else 
   iparm[27]=1; //input must be in single precision as well
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
  int Msize=Msys->size();

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


  int nnzIt=Msys->numberOfNonZeros();

  //
  //put A and C block in the augmented system as At and Ct in the lower triangular part
  //

  if( nA > 0 || nC > 0 )
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

  nnz=augSys.numberOfNonZeros();

  // we need to transpose to get the augmented system in the row-major upper triangular format of  PARDISO 
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

   //
   // symbolic analysis
   //
   /* same for mkl_pardiso and pardiso */
//   mtype = -2;
//   int phase = 11; // analysis
//   int maxfct = 1, mnum = 1, nrhs = 1;
//
//   //iparm[1] = 2; // 2 is for metis, 0 for min degree
//   iparm[7] = 8; // # iterative refinements
//   //iparm[9] =10; // pivot perturbation 10^{-xxx}
//   iparm[23] = 1; //Parallel Numerical Factorization (0 = used in the last years/INTEL calls it classical algorithm , 1 = two-level scheduling)
//   int msglvl = 0; // with statistical information
//
//#ifndef WITH_MKL_PARDISO
//   iparm[2] = num_threads;
//   iparm[10] = 1; // scaling for IPM KKT; used with IPARM(13)=1 or 2
//   iparm[12] = 2; // improved accuracy for IPM KKT; used with IPARM(11)=1;
//                  // if needed, use 2 for advanced matchings and higer accuracy.
//   iparm[32] = 1; // compute determinant - no equivalent for MKL_PARDISO
//   iparm[37] = Msys->size(); //compute Schur-complement
//#else
//   /* From INTEL (instead of iparm[2] which is not defined there):
//   *  You can control the parallel execution of the solver by explicitly setting the MKL_NUM_THREADS environment variable.
//   *  If fewer OpenMP threads are available than specified, the execution may slow down instead of speeding up.
//   *  If MKL_NUM_THREADS is not defined, then the solver uses all available processors.
//   */
//   iparm[10] = 1; // default, scaling for IPM KKT used with either mtype=11/13 or mtype=-2/-4/6 and iparm[12]=1
//   iparm[12] = 1; // 0 disable matching, 1 enable matching, no other settings
//
//   // schur complement stuff
//   iparm[35] and perm todo
//#endif
//
//
//   /* no dparm in MKL PARDISO anymore */
//   pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, eltsAug, rowptrAug, colidxAug, 
//	   NULL, &nrhs,
//         iparm, &msglvl, NULL, NULL, &error
//#ifndef WITH_MKL_PARDISO
//         ,dparm
//#endif
//   );
//
//   if( error != 0 )
//   {
//      printf("PardisoSolver - ERROR during symbolic factorization: %d\n",
//            error);
//      assert(false);
//   }

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
  int* rowptrSC;
  int* colidxSC;
  double* eltsSC;

  computeSC(SC0.size(), R, A, C, F, G, rowptrSC, colidxSC, eltsSC);

  for( int r = 0; r < nSC; r++ )
  {
     for( int ci = rowptrSC[r]; ci < rowptrSC[r + 1]; ci++ )
     {
        const int c = colidxSC[ci];
        assert(c >= r);

        SC0[c][r] += eltsSC[ci];

#ifndef DENSE_USE_HALF
        if( r != c )
           SC0[r][c] += eltsSC[ci];
#endif
     }
  }

  delete[] rowptrSC; delete[] colidxSC; delete[] eltsSC;
}


void PardisoSchurSolver::schur_solve_sparse(SparseGenMatrix& R,
                 SparseGenMatrix& A,
                 SparseGenMatrix& C,
                 SparseGenMatrix& F,
                 SparseGenMatrix& G,
                 SparseSymMatrix& SC0)
{
  int* rowptrSC;
  int* colidxSC;
  double* eltsSC;

  computeSC(SC0.size(), R, A, C, F, G, rowptrSC, colidxSC, eltsSC);

  int* rowptrBase = SC0.krowM();
  int* colidxBase = SC0.jcolM();
  double* eltsBase = SC0.M();

  assert(SC0.size() == nSC);

  // add to summed Schur complement todo: exploit block structure, get start and end of block
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

  delete[] rowptrSC; delete[] colidxSC; delete[] eltsSC;
}


void PardisoSchurSolver::computeSC(int nSCO,
/*const*/SparseGenMatrix& R,
/*const*/SparseGenMatrix& A,
/*const*/SparseGenMatrix& C,
/*const*/SparseGenMatrix& F,
/*const*/SparseGenMatrix& G, int*& rowptrSC, int*& colidxSC, double*& eltsSC)
{
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
   int mtype = -2, error;

#ifndef NDEBUG
#ifndef WITH_MKL_PARDISO
   pardiso_chkmatrix(&mtype, &n, eltsAug, rowptrAug, colidxAug, &error);
   if( error != 0 )
   {
      cout << "PARDISO matrix error " << error << endl;
      exit(1);
   }
#else
   /* enable matrix checker (default disabled) - mkl pardiso does not have chkmatrix */
   iparm[26] = 1;
#endif
#endif

   const int nIter = (int) g_iterNumber;
   const int symbEvery = 5;
   if( (nIter % symbEvery) == 0 )
      doSymbFact = true;

   /* same for mkl_pardiso and pardiso */
   int phase = 22; // numerical factorization
   if( doSymbFact )
      phase = 12; // numerical factorization & symb analysis

   int maxfct = 1; // max number of fact having same sparsity pattern to keep at the same time
   int mnum = 1; // actual matrix (as in index from 1 to maxfct)
   int nrhs = 1;
   int msglvl = pardiso_verbosity;
   //int myRankp; MPI_Comm_rank(MPI_COMM_WORLD, &myRankp);
   //if (myRankp==0) msglvl=1;

   iparm[7] = 8; // max number of iterative refinement steps
   iparm[9] = 13;//todo was 1 // pivot perturbation 10^{-xxx} todo better 7 or higher?
   iparm[30] = 0; // do not specify sparse rhs at this point

#ifndef WITH_MKL_PARDISO
   iparm[2] = num_threads;

   iparm[10] = 1; // scaling for IPM KKT; used with IPARM(13)=1 or 2
   iparm[12] = 2; // improved accuracy for IPM KKT; used with IPARM(11)=1; use 2 for advanced matchings and higher accuracy.
   // iparm[32] = 1; // compute determinant - not available in MKL_PARDISO
#else
   /* From INTEL (instead of iparm[2] which is not defined there):
    *  You can control the parallel execution of the solver by explicitly setting the MKL_NUM_THREADS environment variable.
    *  If fewer OpenMP threads are available than specified, the execution may slow down instead of speeding up.
    *  If MKL_NUM_THREADS is not defined, then the solver uses all available processors.
    */
   iparm[10] = 1; // default, scaling for IPM KKT used with either mtype=11/13 or mtype=-2/-4/6 and iparm[12]=1
   iparm[12] = 1;// 0 disable matching, 1 enable matching, no other settings
#endif


#ifdef PARDISO_PARALLEL_AGGRESSIVE
   iparm[23] = 1; // parallel Numerical Factorization (0=used in the last years, 1=two-level scheduling)

   #ifndef WITH_MKL_PARDISO
   // iparm[1] = 3; // 3 Metis 5.1 (only for PARDISO >= 6.0) - not available in MKL_PARDISO
   iparm[24] = 1;// parallelization for the forward and backward solve. 0=sequential, 1=parallel solve.
   // iparm[27] = 1; // Parallel metis - not for MKL_PARDISO
   #else
   iparm[24] = 2; // not sure but two seems to be the appropriate equivalent here
                  // one rhs -> parallelization, multiple rhs -> parallel forward backward subst
   #endif
#else
   iparm[1] = 2; // 2 is for metis, 0 for min degree
   iparm[23] = 0; // parallel Numerical Factorization (0=used in the last years, 1=two-level scheduling)
   iparm[24] = 0; // parallelization for the forward and backward solve. 0=sequential
#endif

   /* compute schur complement */
#ifndef WITH_MKL_PARDISO
   iparm[37] = nSC; // Msys->size(); // compute Schur-complement

   #ifdef TIMING
   //dumpAugMatrix(n,nnz,iparm[37], eltsAug, rowptrAug, colidxAug);
   //double o=MPI_Wtime();
   #endif
   #ifdef TIMING_FLOPS
   HPM_Start("PARDISOFact");
   #endif

   pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, eltsAug, rowptrAug,
         colidxAug, NULL, &nrhs, iparm, &msglvl, NULL, NULL, &error ,dparm);

 #ifdef TIMING_FLOPS
   HPM_Stop("PARDISOFact");
 #endif
   int nnzSC=iparm[38];
 #ifdef TIMING
   int myRank; MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
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
   // perm array has to be defined is of size of augmented system (dense) and inicates rows/colums we want to have in the schur complement with 1
   // rest set to 0
   // setting all entries to 1 will return the input matrix as schur complement

   // iparm[35] < 0 only possible if iparm[23] != 0 and since iparm[23]=1 disables scaling and matching we have to set it to iparm[23]=10
   iparm[23] = 1; // 10 will lead to a seg fault thrown in pardiso_export or weird behaviour of te first pardiso call (depending on iparm[35])
   iparm[35] = -2;

   int perm[n] = {0};
   for(int i = n - nSC; i < n; ++i)
      perm[i] = 1;

   /* reordering and symbolic factorization */
   phase = 11;

   pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, eltsAug, rowptrAug,
   	   	   colidxAug, perm, &nrhs, iparm, &msglvl, NULL, NULL, &error);

   /* iparm[35] should now contain the number of non-zero entries in the Schur-complement */
   assert(error == 0); assert(iparm[35] >= 0);

   /* preallocation of schur matrix arrays */
   int nnzSC = iparm[35];
   rowptrSC = new int[nSC + 1];
   colidxSC = new int[iparm[35]];
   eltsSC = new double[iparm[35]];

   phase = 22;
   /* reset iparm[35] to -2 */
   iparm[35] = -2;
   int step = 1;

   // mkl pardiso returns arrays with zero-based index via export
   pardiso_export(pt, eltsSC, rowptrSC, colidxSC, &step, iparm, &error);

   #ifdef TIMING
   //dumpAugMatrix(n,nnz,iparm[37], eltsAug, rowptrAug, colidxAug);
   //double o=MPI_Wtime();
   #endif

   #ifdef TIMING_FLOPS
   HPM_Start("PARDISOFact");
   #endif

   /* factorization call */
   pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, eltsAug, rowptrAug,
         colidxAug, perm, &nrhs, iparm, &msglvl, NULL, NULL, &error);

   #ifdef TIMING_FLOPS
   HPM_Stop("PARDISOFact");
   #endif

   #ifdef TIMING
   int myRank; MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
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


   int ind, pend;

   /////////////////////////////////////////////////////
   // transpose the matrix since it is given in lower triangular form and we are using upper triangular form
   /////////////////////////////////////////////////////
   // todo move somewhere else

   //cummulative sum to find the number of elements in each row of At, ie column of A.
   int w[n] = {0};

   int* rowptrSCTP = new int[nSC + 1];
   int* colidxSCTP = new int[nnzSC];
   double* eltsSCTP = new double[nnzSC];

   for(int i = 0; i < nnzSC; i++)
      w[ colidxSC[i] ]++;

   rowptrSCTP[0] = 0;
   for(int i = 1; i <= nSC; i++)
   {
      rowptrSCTP[i] = rowptrSCTP[i-1] + w[i-1];
      w[i-1] = rowptrSCTP[i-1];
   }

   // populate colidxSCTP and eltsSCTP now
   for(int i = 0; i < nSC; i++)
   {
      pend = rowptrSC[ i + 1 ];
      for(int j = rowptrSC[i]; j < pend; j++)
      {
         ind = w[ colidxSC[j] ];
         colidxSCTP[ ind ] = i;
         eltsSCTP[ ind ] = eltsSC[j];
         w[ colidxSC[j] ] = ind + 1;
      }
   }

   for(int i = 0; i < nSC + 1; ++i)
      rowptrSC[i] = rowptrSCTP[i];
   for(int i = 0; i < nnzSC; ++i)
   {
      colidxSC[i] = colidxSCTP[i];
      eltsSC[i] = eltsSCTP[i];
   }

   delete[] colidxSCTP;
   delete[] eltsSCTP;
   delete[] rowptrSCTP;
#endif

// TODO remove when finished debugging mkl - print schur matrix
//   for( int it = 0; it < nSC; it++ )
//   {
//      int k = rowptrSC[it];
//      for( int jt = 0; jt < nSC && k < rowptrSC[it + 1]; jt++ )
//      {
//         if( colidxSC[k] == jt )
//         {
//            std::cout << std::fixed << std::setprecision(2) << eltsSC[k]
//                  << "\t";
//            ++k;
//         }
//         else
//            std::cout << "0\t";
//      }
//      std::cout << std::endl;
//   }
//   std::cout << "-----------------------------------------------------------------------------------" << std::endl;

   assert(subMatrixIsOrdered(rowptrSC, colidxSC, 0, nSC));

}


void PardisoSchurSolver::solve( OoqpVector& rhs_in )
{ 
  SimpleVector& rhs=dynamic_cast<SimpleVector&>(rhs_in);

  /* same for mkl_pardiso and pardiso */
  int mtype=-2, error;
  int phase = 33; // solve and iterative refinement
  int maxfct = 1; // max number of fact having same sparsity pattern to keep at the same time
  int mnum = 1; // actual matrix (as in index from 1 to maxfct)
  int nrhs = 1;

  int msglvl = pardiso_verbosity;
  //int myRankp; MPI_Comm_rank(MPI_COMM_WORLD, &myRankp);
  //if (myRankp==0) msglvl=1;

  iparm[7] = 8; // max number of iterative refinement steps

#ifndef WITH_MKL_PARDISO
  iparm[2] = num_threads;

  iparm[10] = 1; // scaling for IPM KKT; used with IPARM(13)=1 or 2
  iparm[12] = 2; // improved accuracy for IPM KKT; used with IPARM(11)=1; use 2 for advanced matchings and higher accuracy.
#else
  /* From INTEL (instead of iparm[2] which is not defined there):
   *  You can control the parallel execution of the solver by explicitly setting the MKL_NUM_THREADS environment variable.
   *  If fewer OpenMP threads are available than specified, the execution may slow down instead of speeding up.
   *  If MKL_NUM_THREADS is not defined, then the solver uses all available processors.
   */
  iparm[10] = 1; // default, scaling for IPM KKT used with either mtype=11/13 or mtype=-2/-4/6 and iparm[12]=1
  iparm[12] = 1;// 0 disable matching, 1 enable matching, no other settings
#endif


#ifdef PARDISO_PARALLEL_AGGRESSIVE
  iparm[23] = 1; // parallel Numerical Factorization (0=used in the last years, 1=two-level scheduling)

  #ifndef WITH_MKL_PARDISO
  iparm[24] = 1;// parallelization for the forward and backward solve. 0=sequential, 1=parallel solve.
  #else
  iparm[24] = 2; // not sure but two seems to be the appropriate equivalent here; one rhs -> parallelization, multiple rhs -> parallel forward backward subst
  #endif
#else
  iparm[23] = 0; // parallel Numerical Factorization (0=used in the last years, 1=two-level scheduling)
#endif

  assert(nvec_size == n);
  double* const x_n = nvec;
  double* const rhs_n = nvec2;

  const int dim=rhs.length();
  assert(dim >= 0 && dim <= n);

  memset(x_n, 0, dim * sizeof(double));
  memcpy(rhs_n, rhs.elements(), dim * sizeof(double));

  // todo sparsity vector?
  if( n > dim )
     memset(&rhs_n[dim], 0, (n - dim) * sizeof(double));

#ifdef TIMING_FLOPS
  HPM_Start("PARDISOSolve");
#endif

  // todo remove when done debugging
  //////////////////////////////////////
  iparm[1-1] = 1;         /* No solver default */
  iparm[2-1] = 2;         /* Fill-in reordering from METIS */
  iparm[10-1] = 8;        /* Perturb the pivot elements with 1E-13 */
  iparm[11-1] = 0;        /* Use nonsymmetric permutation and scaling MPS */
  iparm[13-1] = 0;        /* Maximum weighted matching algorithm is switched-off (default for symmetric). Try iparm[12] = 1 in case of inappropriate accuracy */
  iparm[14-1] = 0;        /* Output: Number of perturbed pivots */
  iparm[18-1] = -1;       /* Output: Number of nonzeros in the factor LU */
  iparm[19-1] = -1;       /* Output: Mflops for LU factorization */
  iparm[36-1] = 1;        /* Use Schur complement */
  // added
  iparm[5] = 0;
  iparm[23] = 1;
  iparm[30] = 0;
  iparm[35] = -2;
  iparm[7] = 0;
  iparm[20] = 0;
  iparm[24] = 0;
  iparm[26] = 0;

  for( int it = 0; it < n; it++ )
  {
     int k = rowptrAug[it] - 1;
     for( int jt = 0; jt < n && k < rowptrAug[it + 1] - 1; jt++ )
     {
        if( colidxAug[k] - 1 == jt )
        {
           std::cout << std::fixed << std::setprecision(4) << eltsAug[k]
                 << "\t";
           ++k;
        }
        else
           std::cout << "0\t";
     }
     std::cout << std::endl;
  }
  //////////////////////////////////////

  // solving phase
  phase = 33;

  pardiso (pt , &maxfct , &mnum, &mtype, &phase,
	   &n, eltsAug, rowptrAug, colidxAug, 
	   NULL, &nrhs,
	   iparm , &msglvl, rhs_n, x_n, &error
#ifndef WITH_MKL_PARDISO
	   ,dparm
#endif

  );
  assert(error == 0);


  /////TODO remove when done debugging
  //////////////////////////////////////

  // compute the full sparse representation of the sparse symmetric input
  // cout elems per row
  int nelems[n] = {0};
  for(int i = 0; i < n; ++i)
  {
     // diag elem
     for(int j = rowptrAug[i] - 1; j < rowptrAug[i + 1] - 1; ++j)
     {
        if(i == colidxAug[j] - 1)
           nelems[i]++;
        else
        {
           nelems[colidxAug[j] - 1]++;
           nelems[i]++;
        }
     }
  }

  // fill rowptr array
  int rowptrAugFull[n + 1];

  rowptrAugFull[0] = 0;
  for(int i = 0; i < n; ++i)
     rowptrAugFull[i + 1] = rowptrAugFull[i] + nelems[i];

  int colidxAugFull[rowptrAugFull[n]];
  for( int i = 0; i < rowptrAugFull[n]; ++i)
     colidxAugFull[i] = -1;

  double eltsAugFull[rowptrAugFull[n]];

  // fill in col and value
  for(int i = 0; i < n; ++i)
  {
     int rowstart = rowptrAug[i] - 1;
     int rowend = rowptrAug[i+1] - 1;

     for(int k = rowstart; k < rowend; ++k)
     {
        double value = eltsAug[k];
        int col = colidxAug[k] - 1;

        int kk = rowptrAugFull[i];
        int colfull = colidxAugFull[kk];
        while(colfull != -1)
        {
           kk++;
           colfull = colidxAugFull[kk];
        }

        colidxAugFull[kk] = col;
        eltsAugFull[kk] = value;

        if(col != i )
        {
           kk = rowptrAugFull[col];
           int colfull = colidxAugFull[kk];

           while(colfull != -1)
           {
              kk++;
              colfull = colidxAugFull[kk];
           }

           colidxAugFull[kk] = i;
           eltsAugFull[kk] = value;
        }
     }
   }

  // print full matrix
   std::cout << std::endl << "dim: " << dim << "\tn: " << n << std::endl;
   for( int it = 0; it < n; it++ )
   {
      int k = rowptrAugFull[it];
      for( int jt = 0; jt < n && k < rowptrAugFull[it + 1]; jt++ )
      {
         if( colidxAugFull[k] == jt )
         {
            std::cout << std::fixed << std::setprecision(4) << eltsAugFull[k]
                  << "\t";
            ++k;
         }
         else
            std::cout << "0\t";
      }
      std::cout << std::endl;
   }

   std::cout << std::endl;

   // print x
  for(int i = 0; i < n; ++i)
     std::cout << x_n[i] << " ";
   std::cout << std::endl;

   // print Ax == b
   double res[n];
   for( int i = 0; i < n; ++i )
   {
      res[i] = 0.0;
      for( int j = rowptrAugFull[i]; j < rowptrAugFull[i + 1]; ++j )
         res[i] += eltsAugFull[j] * x_n[colidxAugFull[j]];
      std::cout << res[i] << " == " << rhs_n[i] << std::endl;
   }

   exit(1);
   //////////////////////////////////////

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

    //if(res_norm2/rhsNorm>1e-4 && 0==gOuterBiCGIter) {
    //  int myRank; MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    //  char tmp[1024];
    //  sprintf(tmp, "augMat-r%d-i%g-s%g-b%d.dat", myRank, g_iterNumber,  g_scenNum+1, gOuterBiCGIter);
      
    //  cout << "saving matrix to " << tmp << endl;
      
    //  dumpAugMatrix(n,nnz,iparm[37], eltsAug, rowptrAug, colidxAug);
    //  dumpRhs(rhs);
    //}
  }
  delete[] tmp_resid;
#endif
  
  memcpy(&rhs[0], x_n, dim * sizeof(double));
}

void PardisoSchur32Solver::solve( OoqpVector& rhs_in )
{ 
  SimpleVector& rhs=dynamic_cast<SimpleVector&>(rhs_in);

  /* same for mkl_pardiso and pardiso */
  int mtype=-2, error;
  int phase = 33; // solve and iterative refinement
  int maxfct = 1; // max number of fact having same sparsity pattern to keep at the same time
  int mnum = 1; // actual matrix (as in index from 1 to maxfct)
  int nrhs = 1;

  int msglvl = pardiso_verbosity;
  //int myRankp; MPI_Comm_rank(MPI_COMM_WORLD, &myRankp);
  //if (myRankp==0) msglvl=1;

  iparm[7] = 8; // max number of iterative refinement steps


  // todo this is not working at the moment
#ifndef WITH_MKL_PARDISO
  iparm[2] = num_threads;

  iparm[10] = 1; // scaling for IPM KKT; used with IPARM(13)=1 or 2
  iparm[12] = 2; // improved accuracy for IPM KKT; used with IPARM(11)=1; use 2 for advanced matchings and higher accuracy.
#else
  /* From INTEL (instead of iparm[2] which is not defined there):
   *  You can control the parallel execution of the solver by explicitly setting the MKL_NUM_THREADS environment variable.
   *  If fewer OpenMP threads are available than specified, the execution may slow down instead of speeding up.
   *  If MKL_NUM_THREADS is not defined, then the solver uses all available processors.
   */
  iparm[10] = 1; // default, scaling for IPM KKT used with either mtype=11/13 or mtype=-2/-4/6 and iparm[12]=1
  iparm[12] = 1;// 0 disable matching, 1 enable matching, no other settings
#endif

  iparm[23] = 0; //Parallel Numerical Factorization (0=used in the last years, 1=two-level scheduling)
  iparm[24] = 0; //Parallel Numerical Factorization (0=used in the last years, 1=two-level scheduling)
  
  SimpleVector x_n(n);
  SimpleVector rhs_n(n);
  int dim=rhs.length();
  memcpy(&rhs_n[0], rhs.elements(), dim*sizeof(double));
  for(int i=dim; i<n; i++) rhs_n[i]=0.0;

  //double start = MPI_Wtime();
  pardiso (pt , &maxfct , &mnum, &mtype, &phase,
	   &n, eltsAug, rowptrAug, colidxAug, 
	   NULL, &nrhs,
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



// void PardisoSchurSolver::solve( OoqpVector& rhs_in )
// { 
//   SimpleVector& rhs=dynamic_cast<SimpleVector&>(rhs_in);

//   dumpRhs(rhs);

//   int mtype=-2, error;
//   int phase=33;      //solve and iterative refinement
//   int maxfct=1, mnum=1, nrhs=1;
//   iparm[2]=num_threads;
//   //iparm[5]=1;    //replace rhs with sol 
//   iparm[7]=0;    // # of iterative refinements
//   //iparm[1] = 2;// 2 is for metis, 0 for min degree 
//   //iparm[ 9] =12; // pivot perturbation 10^{-xxx} 
//   iparm[10] = 1; // scaling for IPM KKT; used with IPARM(13)=1 or 2
//   iparm[12] = 2; // improved accuracy for IPM KKT; used with IPARM(11)=1; 
//                  // if needed, use 2 for advanced matchings and higer accuracy.
//   iparm[23] = 1; //Parallel Numerical Factorization (0=used in the last years, 1=two-level scheduling)
//   int msglvl=0;  // with statistical information

//   int dim=rhs.length();
//   SimpleVector x(dim); x.setToZero();
//   SimpleVector dx(nvec, dim);  
//   SimpleVector b(dim); b.copyFrom(rhs);
//   double rhsNorm = b.twonorm(); double relResNorm;
//   int refinSteps=0;
//   SimpleVector rhs_n(n);

//   vector<double> histRelRes;
//   SimpleVector x_last(dim);
//   do {

//     memcpy(&rhs_n[0], rhs.elements(), dim*sizeof(double));
//     for(int i=dim; i<n; i++) rhs_n[i]=0.0;

//     double start = MPI_Wtime();
//     pardiso (pt , &maxfct , &mnum, &mtype, &phase,
// 	     &n, eltsAug, rowptrAug, colidxAug, 
// 	     NULL, &nrhs,
// 	     iparm , &msglvl, rhs_n.elements(), nvec, &error, dparm );   

//     cout << "---pardiso:" << MPI_Wtime()-start << endl;

//     if ( error != 0) {
//       printf ("PardisoSchurSolver - ERROR during single rhs: %d\n", error );
//       exit(0);
//     }
//     x.axpy(1.0,dx);
//     //residual
//     rhs.copyFrom(b);
//     Msys->mult(1.0, rhs.elements(), 1,  -1.0, x.elements(),1);
//     relResNorm = rhs.twonorm()/rhsNorm;

//      //compute residual (alternative)
//     double* tmp_resid=new double[dim];
//     memcpy(tmp_resid, b.elements(), dim*sizeof(double));
//     for(int i=0; i<dim; i++) {
//       for(int p=rowptrAug[i]; p<rowptrAug[i+1]; p++) {
// 	int j=colidxAug[p-1]-1;
// 	if(j+1<=dim) {
// 	  //r[i] = r[i] + M(i,j)*x(j)
// 	  tmp_resid[i] -= eltsAug[p-1]*x[j];
	
// 	  if(j!=i) {
// 	    //r[j] = r[j] + M(j,i)*x(i)
// 	    tmp_resid[j] -=  eltsAug[p-1]*x[i];
// 	  }
// 	}
//       }
//     }

//     double res_norm2=0.0; 
//     for(int i=0; i<dim; i++) res_norm2 += tmp_resid[i]*tmp_resid[i];
//     res_norm2 = sqrt(res_norm2);

//     cout << "--- rhs.nrm=" << rhsNorm << " rel.res.nrm=" << relResNorm 
// 	 << " rel.res.nrm=" << res_norm2/rhsNorm
// 	 << endl;
//     delete[] tmp_resid;

//     if(histRelRes.size()>0) {
//       if(histRelRes[histRelRes.size()-1] < relResNorm) {
// 	cout << "PardisoSchurSolver::residual norm increase in iter refin, last="
// 	     << histRelRes[histRelRes.size()-1] << " new=" << relResNorm 
// 	     << ". Restoring iterate..." << endl;
// 	x.copyFrom(x_last);
// 	relResNorm=histRelRes[histRelRes.size()-1];
// 	break;
//       }
//     }
//     x_last.copyFrom(x);
//     histRelRes.push_back(relResNorm);
//     if(relResNorm<1e-9) {
//       break;
//     }
//     refinSteps++;
//   } while(refinSteps<=5);
  
//   if(relResNorm>1e-8) {

//     Ma57Solver slv(Msys);
//     slv.matrixChanged();
//     SimpleVector x_ma(dim); x_ma.copyFrom(b);
//     slv.solve(x_ma);
    
//     Msys->mult(1.0, b.elements(), 1,  -1.0, x_ma.elements(),1);
//     double relResNorm_ma57 = b.twonorm()/rhsNorm;
//     //cout <<"Ma57 rel resid norm:" << relResNorm_ma57 << endl << endl;
//     if(relResNorm_ma57 < relResNorm) {
//       x.copyFrom(x_ma);
//       relResNorm=relResNorm_ma57;
//     }
//     if(relResNorm>1e-8) {

//       cout << "Pardiso iter refinements: " << refinSteps 
// 	   << "   rel resid nrm=" << relResNorm 
// 	   << ". Ma57 rel resid=" << relResNorm_ma57 << endl;
//       //cout << "Pardiso conv history:";
//       //for(size_t it=0; it<histRelRes.size(); it++)
//       //cout << histRelRes[it] << " ";
//       //cout << endl;
//       // cout <<"Ma57 rel resid norm:" << relResNorm_ma57 << endl << endl; 
//     }
//   }

//   rhs.copyFrom(x);
// }



void PardisoSchurSolver::solve(GenMatrix& rhs_in)
{
  assert(false && "Function not supported. Use PardisoSolver for this functionality.");
}

PardisoSchurSolver::~PardisoSchurSolver()
{
  int phase = -1; /* Release internal memory . */
  int mtype = -2;
  int maxfct=  1, mnum=1, nrhs=1, msglvl=0, error;
  
  pardiso (pt, &maxfct, &mnum, &mtype, &phase,
	   &n, NULL, rowptrAug, colidxAug, NULL, &nrhs,
	   iparm, &msglvl, NULL, NULL, &error
#ifndef WITH_MKL_PARDISO
	   , dparm
#endif
  );
  if ( error != 0) {
    printf ("PardisoSchurSolver - ERROR in pardiso release: %d", error ); 
  }
  
  if(rowptrAug) delete[] rowptrAug;
  if(colidxAug) delete[] colidxAug;
  if(eltsAug)   delete[] eltsAug;
  
  if(nvec) delete[] nvec;
  if(nvec2) delete[] nvec2;
}

#ifdef STOCH_TESTING
int dumpAugMatrix(int n, int nnz, int nSys,
		  double* elts, int* rowptr, int* colidx, 
		  const char* fname)
{
  char filename[1024];
  int myRank; MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
  if(fname==NULL) 
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
  if(fname==NULL) 
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
// int dumpSysMatrix(SparseSymMatrix& M)
// {
//   int n=M.size(), nnz=M.numberOfNonZeros();

//   int*  rowptr = new int[n+1];
//   int*  colidx = new int[nnz];
//   double* elts = new double[nnz];

//   M.getStorageRef().transpose(rowptr,colidx,elts);

//   //FORTRAN indexes
//   for(int it=0; it<n+1; it++) rowptr[it]++;
//   for(int it=0; it<nnz; it++) colidx[it]++;

//   char fname[128];
//   strcpy(fname, "matDump.dat");
//   dumpAugMatrix(n, nnz, elts, rowptr, colidx, fname);
//   return 0;
// }

/*
int pardiso_stuff(int n, int nnz, int n0, int* rowptr, int* colidx, double* elts, const vector< vector<double> >& vRhs)
{
  void  *pt[64];
  int iparm[64];
  double dparm[64];
  int num_threads;

  char *var = getenv("OMP_NUM_THREADS");
  if(var != NULL) {
    sscanf( var, "%d", &num_threads );
    cout << "Using " << num_threads << " threads." << endl;
  }
  else {
    printf("Set environment OMP_NUM_THREADS\n");
    exit(1);
  }
  int solver=0, mtype=-2, error;
  pardisoinit (pt,  &mtype, &solver, iparm, dparm, &error);
  if (error!=0) {
    cout << "Pardiso solver ERROR during pardisoinit:" << error << "." << endl;
    exit(1);
  }

  pardiso_chkmatrix(&mtype, &n, elts, rowptr, colidx, &error);
  if(0!=error)
    cout << "Error in matrix detected: " << error << endl;
  //int zero=0;
  //pardiso_printstats (&mtype, &n, elts, rowptr, colidx, &zero, NULL, &error);

  //
  // symbolic analysis & numerical factorization
  //
  int phase=12; //analysis & fact
  int maxfct=1, mnum=1, nrhs=1;
  iparm[2]=num_threads;
  iparm[7]=8;     //# iterative refinements
  //iparm[1] = 2; // 2 is for metis, 0 for min degree
  //iparm[ 9] =10; // pivot perturbation 10^{-xxx}
  iparm[10] = 1; // scaling for IPM KKT; used with IPARM(13)=1 or 2
  iparm[12] = 2; // improved accuracy for IPM KKT; used with IPARM(11)=1;
                 // if needed, use 2 for advanced matchings and higer accuracy.
  iparm[23] = 1; //Parallel Numerical Factorization (0=used in the last years, 1=two-level scheduling)

  int msglvl=0;  // with statistical information
  iparm[32] = 0; // compute determinant
  if(n0!=n)  iparm[37] = n0; //compute Schur-complement

  double start=omp_get_wtime();
  pardiso (pt, &maxfct , &mnum, &mtype, &phase,
           &n, elts, rowptr, colidx,
           NULL, &nrhs,
           iparm , &msglvl, NULL, NULL, &error, dparm );
  cout << "Factorization: " << omp_get_wtime()-start << endl;
  if ( error != 0) {
    printf ("PardisoSolver - ERROR during factorization: %d\n", error );
    exit(1);
  }


  ////////////////////////////////////////////////////////////
  int nSC=n-n0;
  if(nSC>0) {
    int nnzSC=iparm[38];
    int* rowptrSC =new int[nSC+1];
    int* colidxSC =new int[nnzSC];
    double* eltsSC=new double[nnzSC];
    
    pardiso_get_schur(pt, &maxfct, &mnum, &mtype, eltsSC, rowptrSC, colidxSC);
    //!cout << "Schur complement (2nd stage) n=" << nSC << "   nnz=" << nnzSC << endl;
    
    //convert back to C/C++ indexing
    for(int it=0; it<nSC+1; it++) rowptrSC[it]--;
    for(int it=0; it<nnzSC; it++) colidxSC[it]--;
    
    
    //for(int r=0; r<nSC; r++) {
    //  for(int ci=rowptrSC[r]; ci<rowptrSC[r+1]; ci++) {
    //	int c=colidxSC[ci];
    //	SC0[r][c] += eltsSC[ci];
    //	if(r!=c)
    //	  SC0[c][r] += eltsSC[ci];
    //}
    //}
    delete[] rowptrSC; delete[] colidxSC; delete[] eltsSC;
  }

  



  ////////////////////////////////////////////////////////////


  double* rhs=new double[n];
  double* sol=new double[n];
  double* resid=new double[n0];
  //solve with each rhs
  for(int r=0; r<vRhs.size(); r++) {
    memcpy(rhs,   &vRhs[r][0], n0*sizeof(double));
    memcpy(resid, &vRhs[r][0], n0*sizeof(double));


    //
    // solve
    //
    phase=33;
    iparm[7]=8;    // # of iterative refinements
    int nrhs=1;
    msglvl=0;

    start = omp_get_wtime();
    pardiso (pt , &maxfct , &mnum, &mtype, &phase,
	     &n, elts, rowptr, colidx, 
	     NULL, &nrhs,
	     iparm , &msglvl, rhs, sol, &error, dparm );   
    double t=omp_get_wtime()-start; 

    //cout << "Sol:";
    //for(int i=2000; i<2003; i++) cout << sol[i] << " ";
    //cout << endl;

    //compute residual
    for(int i=0; i<n0; i++) {
      for(int p=rowptr[i]; p<rowptr[i+1]; p++) {
	int j=colidx[p-1]-1;
	if(j+1<=n0) {
	  //r[i] = r[i] + M(i,j)*x(j)
	  resid[i] -= elts[p-1]*sol[j];
	
	  if(j!=i) {
	    //r[j] = r[j] + M(j,i)*x(i)
	    resid[j] -=  elts[p-1]*sol[i];
	  }
	}
      }
    }
    double res_norm=0.0, rhs_norm=0.0;
    for(int i=0; i<n0; i++) res_norm += resid[i]*resid[i];
    for(int i=0; i<n0; i++) rhs_norm += rhs[i]*rhs[i];
    res_norm=sqrt(res_norm); rhs_norm=sqrt(rhs_norm);
    cout << "Solve: " << t << " sec."
         << " Rhs.nrm: " << rhs_norm 
	 << " Rel.resid.nrm=" << res_norm/rhs_norm
	 << endl;
  }
  delete[] rhs; delete[] sol; delete[] resid;

  //
  // clean-up
  //
  phase = -1; // Release internal memory .

  pardiso (pt, &maxfct, &mnum, &mtype, &phase,
           &n, NULL, rowptr, colidx, NULL, &nrhs,
           iparm, &msglvl, NULL, NULL, &error, dparm );
  if ( error != 0) {
    printf ("Pardiso solver - ERROR in pardiso release: %d", error );
    exit(1);
  }

  return 0; 
}
*/
