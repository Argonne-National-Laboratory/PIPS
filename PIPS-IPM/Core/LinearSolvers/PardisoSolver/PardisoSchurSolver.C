/* PIPS-IPM                                                             
 * Authors: Cosmin G. Petra, Miles Lubin
 * (C) 2012 Argonne National Laboratory, see documentation for copyright
 */
#include <stdlib.h>
#include <iostream>
using namespace std;

#include "PardisoSchurSolver.h"
#include "SparseStorage.h"
#include "SparseSymMatrix.h"
#include "SimpleVector.h"
#include "SimpleVectorHandle.h"
#include "DenseGenMatrix.h"
#include <cstdlib>

#include "Ma57Solver.h"

#include "mpi.h"

#ifdef HAVE_GETRUSAGE
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#endif

extern int gOoqpPrintLevel;

using namespace std;

extern "C" void pardisoinit (void   *, int    *,   int *, int *, double *, int *);
extern "C" void pardiso     (void   *, int    *,   int *, int *,    int *, int *, 
                  double *, int    *,    int *, int *,   int *, int *,
                     int *, double *, double *, int *, double *);
extern "C" void pardiso_schur(void*, int*, int*, int*, double*, int*, int*);

extern "C" void pardiso_chkmatrix  (int *, int *, double *, int *, int *, int *);
extern "C" void pardiso_chkvec     (int *, int *, double *, int *);
extern "C" void pardiso_printstats (int *, int *, double *, int *, int *, int *,
                           double *, int *);


PardisoSchurSolver::PardisoSchurSolver( SparseSymMatrix * sgm )
{
  Msys = sgm;
  n = -1;
  nvec=NULL;
  // - we do not have the augmented system yet; most of intialization done during the 
  // first solve call

  first = true; firstSolve = true;
   
  /* Numbers of processors, value of OMP_NUM_THREADS */
  char *var = getenv("OMP_NUM_THREADS");
  if(var != NULL)
    sscanf( var, "%d", &num_threads );
  else {
    printf("Set environment OMP_NUM_THREADS");
    exit(1);
  }
}

void PardisoSchurSolver::firstCall()
{
  int solver=0, mtype=-2, error;
  pardisoinit (pt,  &mtype, &solver, iparm, dparm, &error); 
  if (error!=0) {
    cout << "PardisoSchurSolver ERROR during pardisoinit:" << error << "." << endl;
    exit(1);
  }
} 

// this function is called only once and creates the augmented system
void PardisoSchurSolver::firstSolveCall(SparseGenMatrix& R, 
					SparseGenMatrix& A,
					SparseGenMatrix& C)
{
  int nR,nA,nC;
  int nnz=0;
  R.getSize(nR,nSC); nnz += R.numberOfNonZeros();
  A.getSize(nA,nSC); nnz += A.numberOfNonZeros();
  C.getSize(nC,nSC); nnz += C.numberOfNonZeros();
  int Msize=Msys->size();

  n = nR+nA+nC+nSC;
  assert( Msize == nR+nA+nC );
  nnz += Msys->numberOfNonZeros();
  nnz += nSC; //space for the 0 diagonal of 2x2 block

  // the lower triangular part of the augmented system in row-major
  SparseSymMatrix augSys( n, nnz);
  
  //
  //put (1,1) block in the augmented system
  //
  memcpy(augSys.getStorageRef().krowM, Msys->getStorageRef().krowM, sizeof(int)*Msize);
  memcpy(augSys.getStorageRef().jcolM, Msys->getStorageRef().jcolM, sizeof(int)*Msys->numberOfNonZeros());
  memcpy(augSys.getStorageRef().M,     Msys->getStorageRef().M,     sizeof(double)*Msys->numberOfNonZeros());

  //
  //put C block in the augmented system as C^T in the lower triangular part
  //
  if(C.numberOfNonZeros()>0 && nC>0) {
    int nnzIt=Msys->numberOfNonZeros();
    
    SparseGenMatrix Ct(nSC,nC,C.numberOfNonZeros());
    int* krowCt=Ct.getStorageRef().krowM;
    int* jcolCt=Ct.getStorageRef().jcolM;
    double *MCt=Ct.getStorageRef().M;
    
    C.getStorageRef().transpose(krowCt, jcolCt, MCt);
    
    int colShift=nR+nA;
    
    int* krowAug= augSys.getStorageRef().krowM;
    int* jcolAug = augSys.getStorageRef().jcolM;
    double* MAug = augSys.getStorageRef().M;
    
    int row=Msize;
    for(; row<n; row++) {
      krowAug[row]=nnzIt;
      
      for(int c=krowCt[row-Msize]; c< krowCt[row-Msize+1]; c++) {
	
	int j=jcolCt[c];
	  
	jcolAug[nnzIt]=j+colShift;
	MAug[nnzIt]   =MCt[c];
	nnzIt++;
      }
      //add the zero from the diagonal
      jcolAug[nnzIt]=row;
      MAug[nnzIt]=0.0;
      nnzIt++;
      
    }
    krowAug[row]=nnzIt;
  }

  // -- to do
  //put A/R block in the augmented system as A^T/R^T in the lower triangular part
  assert(A.numberOfNonZeros()==0);
  //

  nnz=augSys.numberOfNonZeros();
  // we need to transpose to get the augmented system in the row-major upper triangular format of  PARDISO 
  rowptrAug = new int[n+1];
  colidxAug = new int[nnz];
  eltsAug   = new double[nnz];

  augSys.getStorageRef().transpose(rowptrAug,colidxAug,eltsAug);
      

  //save the indeces for diagonal entries for a streamlined later update
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

  //
  // symbolic analysis
  //
  int mtype=-2, error;
  int phase=11; //analysis
  int maxfct=1, mnum=1, nrhs=1;
  iparm[2]=num_threads;
  iparm[7]=3;     //# iterative refinements
  //iparm[1] = 2; // 2 is for metis, 0 for min degree 
  //iparm[ 9] =10; // pivot perturbation 10^{-xxx} 
  iparm[10] = 1; // scaling for IPM KKT; used with IPARM(13)=1 or 2
  iparm[12] = 2; // improved accuracy for IPM KKT; used with IPARM(11)=1; 
                 // if needed, use 2 for advanced matchings and higer accuracy.
  iparm[23] = 1; //Parallel Numerical Factorization (0=used in the last years, 1=two-level scheduling)

  int msglvl=0;  // with statistical information
  iparm[32] = 1; // compute determinant
  iparm[37] = Msys->size(); //compute Schur-complement

  pardiso (pt , &maxfct , &mnum, &mtype, &phase,
	   &n, eltsAug, rowptrAug, colidxAug, 
	   NULL, &nrhs,
	   iparm , &msglvl, NULL, NULL, &error, dparm );

  if ( error != 0) {
    printf ("PardisoSolver - ERROR during symbolic factorization: %d\n", error );
    assert(false);
  }
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
}
 
void PardisoSchurSolver::schur_solve(SparseGenMatrix& R, 
				     SparseGenMatrix& A,
				     SparseGenMatrix& C,
				     DenseSymMatrix& SC0)
{
  if(firstSolve) { 
    firstSolveCall(R,A,C); firstSolve=false; 
  } else {

    //update diagonal entries in the PARDISO aug sys
    double* eltsMsys = Msys->getStorageRef().M;
    map<int,int>::iterator it;
    for(it=diagMap.begin(); it!=diagMap.end(); it++)
      eltsAug[it->second] = eltsMsys[it->first];
  }

  // call PARDISO
  int mtype=-2, error;
  pardiso_chkmatrix(&mtype,&n, eltsAug, rowptrAug, colidxAug, &error);
  if(error != 0) {
    cout << "PARDISO check matrix error" << error << endl;
    exit(1);
  }
  

  int phase=22; //Numerical factorization
  int maxfct=1, mnum=1, nrhs=1;
  iparm[2]=num_threads;
  iparm[7]=3;     //# iterative refinements
  //iparm[1] = 2; // 2 is for metis, 0 for min degree 
  //iparm[ 9] = 10; // pivot perturbation 10^{-xxx} 
  iparm[10] = 1; // scaling for IPM KKT; used with IPARM(13)=1 or 2
  iparm[12] = 2; // improved accuracy for IPM KKT; used with IPARM(11)=1; 
                 // if needed, use 2 for advanced matchings and higer accuracy.
  iparm[23] = 1; //Parallel Numerical Factorization (0=used in the last years, 1=two-level scheduling)

  int msglvl=0;  // with statistical information
  //iparm[32] = 1; // compute determinant
  iparm[37] = Msys->size(); //compute Schur-complement

#ifdef TIMING
  double o=MPI_Wtime();
#endif

  pardiso (pt , &maxfct , &mnum, &mtype, &phase,
	   &n, eltsAug, rowptrAug, colidxAug, 
	   NULL, &nrhs,
	   iparm , &msglvl, NULL, NULL, &error, dparm );


  int nnzSC=iparm[38];

#ifdef TIMING
  cout << "PARDISOSCHUR FACT(AUGMAT) " << MPI_Wtime()-o << endl;
  //cout << "NNZ(SCHUR) " << nnzSC << "    SPARSITY " << nnzSC/(1.0*nSC*nSC) << endl;
#endif
  if ( error != 0) {
    printf ("PardisoSolver - ERROR during factorization: %d\n", error );
    assert(false);
  }

  int* rowptrSC =new int[nSC+1];
  int* colidxSC =new int[nnzSC];
  double* eltsSC=new double[nnzSC];

  pardiso_schur(pt, &maxfct, &mnum, &mtype, eltsSC, rowptrSC, colidxSC);
  //!cout << "Schur complement (2nd stage) n=" << nSC << "   nnz=" << nnzSC << endl;

  //convert back to C/C++ indexing
  for(int it=0; it<nSC+1; it++) rowptrSC[it]--;
  for(int it=0; it<nnzSC; it++) colidxSC[it]--;


  for(int r=0; r<nSC; r++) {
    for(int ci=rowptrSC[r]; ci<rowptrSC[r+1]; ci++) {
      int c=colidxSC[ci];
      SC0[r][c] += eltsSC[ci];
      if(r!=c)
	SC0[c][r] += eltsSC[ci];
    }
  }

  delete[] rowptrSC; delete[] colidxSC; delete[] eltsSC;
}


void PardisoSchurSolver::solve( OoqpVector& rhs_in )
{ 

  //Ma57Solver slv(Msys);
  //slv.matrixChanged();
  //slv.solve(rhs_in);
  //return;
  
  SimpleVector& rhs=dynamic_cast<SimpleVector&>(rhs_in);

  int mtype=-2, error;
  int phase=33;      //solve and iterative refinement
  int maxfct=1, mnum=1, nrhs=1;
  iparm[2]=num_threads;
  //iparm[5]=1;    //replace rhs with sol 
  iparm[7]=0;    // # of iterative refinements
  //iparm[1] = 2;// 2 is for metis, 0 for min degree 
  //iparm[ 9] =12; // pivot perturbation 10^{-xxx} 
  iparm[10] = 1; // scaling for IPM KKT; used with IPARM(13)=1 or 2
  iparm[12] = 2; // improved accuracy for IPM KKT; used with IPARM(11)=1; 
                 // if needed, use 2 for advanced matchings and higer accuracy.
  iparm[23] = 1; //Parallel Numerical Factorization (0=used in the last years, 1=two-level scheduling)
  int msglvl=0;  // with statistical information

  int dim=rhs.length();
  SimpleVector x(dim); x.setToZero();
  SimpleVector dx(nvec, dim);  
  SimpleVector b(dim); b.copyFrom(rhs);
  double rhsNorm = b.twonorm(); double relResNorm;
  int refinSteps=0;
  SimpleVector rhs_n(n);

  vector<double> histRelRes;
  SimpleVector x_last(dim);
  do {

    memcpy(&rhs_n[0], rhs.elements(), dim*sizeof(double));
    for(int i=dim; i<n; i++) rhs_n[i]=0.0;

    pardiso (pt , &maxfct , &mnum, &mtype, &phase,
	     &n, eltsAug, rowptrAug, colidxAug, 
	     NULL, &nrhs,
	     iparm , &msglvl, rhs_n.elements(), nvec, &error, dparm );   
    if ( error != 0) {
      printf ("PardisoSchurSolver - ERROR during single rhs: %d\n", error );
      exit(0);
    }
    x.axpy(1.0,dx);
    //residual
    rhs.copyFrom(b);
    Msys->mult(1.0, rhs.elements(), 1,  -1.0, x.elements(),1);
    relResNorm = rhs.twonorm()/rhsNorm;

    if(histRelRes.size()>0) {
      if(histRelRes[histRelRes.size()-1] < relResNorm) {
	cout << "PardisoSchurSolver::residual norm increase in iter refin, last="
	     << histRelRes[histRelRes.size()-1] << " new=" << relResNorm 
	     << ". Restoring iterate..." << endl;
	x.copyFrom(x_last);
	relResNorm=histRelRes[histRelRes.size()-1];
	break;
      }
    }
    x_last.copyFrom(x);
    histRelRes.push_back(relResNorm);
    if(relResNorm<1e-9) {
      break;
    }
    refinSteps++;
  } while(refinSteps<=9);
  
  if(relResNorm>1e-8) {

    Ma57Solver slv(Msys);
    slv.matrixChanged();
    SimpleVector x_ma(dim); x_ma.copyFrom(b);
    slv.solve(x_ma);
    
    Msys->mult(1.0, b.elements(), 1,  -1.0, x_ma.elements(),1);
    double relResNorm_ma57 = b.twonorm()/rhsNorm;
    //cout <<"Ma57 rel resid norm:" << relResNorm_ma57 << endl << endl;
    if(relResNorm_ma57 < relResNorm) {
      x.copyFrom(x_ma);
      relResNorm=relResNorm_ma57;
    }
    if(relResNorm>1e-8) {

      cout << "Pardiso iter refinements: " << refinSteps 
	   << "   rel resid nrm=" << relResNorm 
	   << "Ma57 rel resid =" << relResNorm_ma57 << endl;
      //cout << "Pardiso conv history:";
      //for(size_t it=0; it<histRelRes.size(); it++)
      //cout << histRelRes[it] << " ";
      //cout << endl;
      // cout <<"Ma57 rel resid norm:" << relResNorm_ma57 << endl << endl; 
    }
  }

  rhs.copyFrom(x);
}

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
	   iparm, &msglvl, NULL, NULL, &error, dparm );
  if ( error != 0) {
    printf ("PardisoSchurSolver - ERROR in pardiso release: %d", error ); 
  }
  if(rowptrAug) delete[] rowptrAug;
  if(colidxAug) delete[] colidxAug;
  if(eltsAug)   delete[] eltsAug;

  if(nvec) delete[] nvec;
}

