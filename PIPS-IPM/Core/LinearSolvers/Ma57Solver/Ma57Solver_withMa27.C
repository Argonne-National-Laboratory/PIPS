/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP
 * Modified by Cosmin Petra to perform solves with the factors.
 */

#include "Ma57Solver.h"
#include "SparseStorage.h"
#include "SparseSymMatrix.h"
#include "SimpleVector.h"
#include "SimpleVectorHandle.h"
#include "DenseGenMatrix.h"
#include "pipsport.h"

#ifdef HAVE_GETRUSAGE
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#endif

///!
#include <cmath>
#include "Ma27Solver.h"

extern int gOoqpPrintLevel;

#ifndef MIN
#define MIN(a,b) ((a > b) ? b : a)
#endif

#include <mpi.h>

void dumpdata(int* irow, int* jcol, double*M, int n, int nnz)
{
  printf("======================================================\n");
  for(int i=0; i<nnz; i++)  printf("%6d %6d %10.2f\n", irow[i], jcol[i], M[i]);
  printf("\n");
  /*  for(int i=0; i<nnz; i++)  printf("%10d ", jcol[i]);
  printf("\n");
  for(int i=0; i<nnz; i++)  printf("%10.2f ", M[i]);
  printf("\n");
  */
  printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
}


Ma27Solver* ma27=nullptr;
SparseSymMatrix * sgm_=nullptr;

Ma57Solver::Ma57Solver( SparseSymMatrix * sgm )
{
  irowM = 0;
  jcolM = 0;
  fact  = 0;
  ifact = 0;
  keep  = 0;

  iworkn  = nullptr; dworkn = nullptr;
  niworkn = ndworkn = 0;

  ipessimism = 1.4;
  rpessimism = 1.4;

  FNAME(ma57id)( cntl, icntl );
  //icntl[1] = -1; // don't print warning messages
  icntl[8] = 10; // up to 10 steps of iterative refinement
  icntl[5] = 5; // 4 use Metis; 5 automatic choice(MA47 or Metis); 3 min
	      // degree ordering as in MA27; 2 use MC47;

  // set initial value of "Treat As Zero" parameter
  kTreatAsZero = 1.e-10;      this->setTreatAsZero();

  // set initial value of Threshold parameter
  kThresholdPivoting = 1.e-5; this->setThresholdPivoting();

  // set the largest value of ThresholdPivoting parameter we are
  // willing to tolerate.
  kThresholdPivotingMax = 1.e-1;

  // set the increase factor for ThresholdPivoting parameter
  kThresholdPivotingFactor = 10.0;

  // set the required precision for each linear system solve
  kPrecision = 1.e-9;

  mStorage = sgm->getStorageHandle();
  n        = mStorage->n;
  M        = mStorage->M;

  nnz = mStorage->krowM[n];

  //!
  ma27=new Ma27Solver(sgm);
  sgm_=sgm;

}

void Ma57Solver::firstCall()
{
  irowM = new int[nnz];
  jcolM = new int[nnz];

  int * krowM = mStorage->krowM;
  for( int i = 0; i < n; i++ ) {
    for( int k = krowM[i]; k < krowM[i+1]; k++ ) {
      irowM[k] = i + 1;
    }
  }
  for( int k = 0; k < nnz; k++ ) {
    jcolM[k] = mStorage->jcolM[k] + 1;
  }

  lkeep = ( nnz > n ) ? (5 * n + 2 *nnz + 42) : (6 * n + nnz + 42);
  keep = new int[lkeep];

  int * iwork = new int[5 * n];
  FNAME(ma57ad)( &n, &nnz, irowM, jcolM, &lkeep, keep, iwork, icntl,
	   info, rinfo );

  delete [] iwork;

  lfact = info[8];
  lfact = 2*(int) (rpessimism * lfact);
  fact  = new double[lfact];
  lifact = info[9];
  lifact = (int) (ipessimism * lifact);
  ifact  = new int[lifact];

}
void Ma57Solver::diagonalChanged( int /* idiag */, int /* extent */ )
{
  cout << "aaaaaaaaaa" << endl;
  this->matrixChanged();

  ma27->matrixChanged();

}

void Ma57Solver::matrixChanged()
{
  if( !keep ) this->firstCall();

  int * iwork = new_iworkn(n);

  int done = 0, tries = 0;;
  do {
#ifdef HAVE_GETRUSAGE
    rusage before;
    if( gOoqpPrintLevel >= 100 ) {
      getrusage( RUSAGE_SELF, &before );
    }
#endif

    //!log
    //dumpdata(irowM, jcolM, M, n, nnz);
    //assert(false);

    FNAME(ma57bd)( &n,       &nnz,    M,     fact,  &lfact,  ifact,
	     &lifact,  &lkeep,  keep,  iwork,  icntl,  cntl,
	     info,     rinfo );
#ifdef HAVE_GETRUSAGE
    rusage  after;
    if( gOoqpPrintLevel >= 100 ) {
      getrusage( RUSAGE_SELF, &after );
      cout << "For try " << tries + 1
	   << " the factorization took "
	   << (double) (after.ru_utime.tv_sec - before.ru_utime.tv_sec)
	+ (after.ru_utime.tv_usec - before.ru_utime.tv_usec) / 1000000.0
	   << " seconds.\n";
    }

    //cout << "MA57 factorization: " << info[21] << " 2x2 pivots were used in the factorization" <<  endl;
    //cout << "MA57 factorization: " << info[23] << " negative eigenvalues were found." <<  endl;
#endif
    cout << "Ma57: numerical rank: " << info[24] << " (mat size: " << n << ")" << endl;

    if( info[0] != 0 ) cout << "ma57bd: Factorization: info[0]=: " << info[0] << endl;
    //assert(false);
    switch( info[0] ) {
    case 0: done = 1;
      break;
    case -3: {
      int ic = 0;
      int lnfact = (int) (info[16] * rpessimism);
      double * newfact = new double[lnfact];
      FNAME(ma57ed)( &n, &ic, keep, fact, &lfact, newfact, &lnfact,
	       ifact, &lifact, ifact, &lifact, info );
      delete [] fact;
      fact = newfact;
      lfact = lnfact;
      rpessimism *= 1.1;
      cout << "Resizing real part. pessimism = " << rpessimism << endl;
    }; break;
    case -4: {
      int ic = 1;
      int lnifact = (int) (info[17] * ipessimism);
      int * nifact = new int[ lnifact ];
      FNAME(ma57ed)( &n, &ic, keep, fact, &lfact, fact, &lfact,
	       ifact, &lifact, nifact, &lnifact, info );
      delete [] ifact;
      ifact = nifact;
      lifact = lnifact;
      ipessimism *= 1.1;
      cout << "Resizing int part. pessimism = " << ipessimism << endl;
    }; break;
    default:
      if( info[0] >= 0 ) done = 1;
      assert( info[0] >= 0 );
    } // end switch
    tries++;
  } while( !done );
  freshFactor = 1;

  //delete [] iwork;
}

void Ma57Solver::solve( OoqpVector& rhs_in )
{

  SimpleVector x_ma27(n), rhs_sav(n);
  x_ma27.copyFrom(rhs_in);
  rhs_sav.copyFrom(rhs_in);

	//    int job = 0; // Solve using A
	//    int one = 1;

	//    SimpleVectorHandle work( new SimpleVector(n) );
	//    SimpleVector & rhs = dynamic_cast<SimpleVector &>(rhs_in);

	//    double * drhs = rhs.elements();
	//    double * dwork  = work->elements();

	//    int * iwork = new int[n];

	//    rusage before;
	//    getrusage( RUSAGE_SELF, &before );

	//    FNAME(ma57cd)( &job,       &n,
	//  	   fact,       &lfact,    ifact,  &lifact,
	//  	   &one,       drhs,      &n,
	//  	   dwork,      &n,        iwork,
	//  	   icntl,      info );


	//    rusage after;
	//    getrusage( RUSAGE_SELF, &after );
	//    cout << "Solution with the factored matrix took "
	//         << (double) (after.ru_utime.tv_sec - before.ru_utime.tv_sec)
	//        + (after.ru_utime.tv_usec - before.ru_utime.tv_usec) / 1000000.0
	//  	 << " seconds.\n";

	//    delete [] iwork;
	int job = 0;
	if( freshFactor ) {
		icntl[8] = 1; // No iterative refinement
	} else {
		icntl[8] = 10; // Iterative refinement
	}
	// MIKE: are these structure ever released??

	SimpleVectorHandle x( new SimpleVector(n) );
	SimpleVectorHandle resid( new SimpleVector(n) );
	SimpleVectorHandle work( new SimpleVector(5 * n) );
	SimpleVector & rhs = dynamic_cast<SimpleVector &>(rhs_in);

	double * drhs = rhs.elements();
	double * dx   = x->elements();
	double * dresid = resid->elements();
	double * dwork  = work->elements();

	int * iwork = new_iworkn(n);


  /*static int s = 0;
  int mype; MPI_Comm_rank(MPI_COMM_WORLD,&mype);
  if (mype == 0 && s==105) {
    printf("RHS IN\n\n");
    for (int i = 0; i < n; i++) {
      printf("%d: %.10E\n", i, rhs[i]);
    }
  }*/

#ifdef HAVE_GETRUSAGE
	rusage before;
	if( gOoqpPrintLevel >= 100 ) {
		getrusage( RUSAGE_SELF, &before );
	}
#endif

	int done = 0;
	int refactorizations = 0;
	int dontRefactor =  (kThresholdPivoting > kThresholdPivotingMax);
  while( !done && refactorizations < 10 ) {

    icntl[9]=1;//condition number

    FNAME(ma57dd)( &job,       &n,        &nnz,   M,        irowM,   jcolM,
        fact,       &lfact,    ifact,  &lifact,  drhs,    dx,
        dresid,      dwork,    iwork,  icntl,    cntl,    info,
        rinfo );
    if( resid->infnorm() < kPrecision*( 1 + rhs.infnorm() ) ) {
      // resids are fine, use them
      done = 1;
      cout << "Ma57: relative norm of residuals for linear system: "
	   <<  resid->infnorm()/rhs.infnorm() << endl;
      cout << "Ma57: condition number: " << rinfo[10] << " " << rinfo[11] << " " << rinfo[12] << endl;
    } else {
      // resids aren't good enough.
      if( freshFactor ) { // We weren't doing iterative refinement,
        // let's do so
        job = 2;
        icntl[8] = 10;
        // Mark this factorization as stale
        freshFactor = 0;
        // And grow more pessimistic about the next factorization
        if( kThresholdPivoting >= kThresholdPivotingMax ) {
          // We have already refactored as with a high a pivtol as we
          // are willing to use
          dontRefactor = 1;
        } else {
          // refactor with a higher Threshold Pivoting parameter
          kThresholdPivoting *= kThresholdPivotingFactor;
          if( kThresholdPivoting > kThresholdPivotingMax )
            kThresholdPivoting = kThresholdPivotingMax;
          this->setThresholdPivoting();
          cout << "Setting ThresholdPivoting parameter to "
            << kThresholdPivoting << " for future factorizations" << endl;
        }
      } else if ( dontRefactor ) {
        // We might have tried a refactor, but the pivtol is above our
        // limit.
        done = 1;
      } else {
        // Otherwise, we have already tried iterative refinement, and
        // have already increased the ThresholdPivoting parameter
        cout << "Refactoring with Threshold Pivoting parameter"
          << kThresholdPivoting << endl;
        this->matrixChanged();
        refactorizations++;
        // be optimistic about the next factorization
        job = 0;
        icntl[8] = 1;
      } // end else we hava already tried iterative refinement
    } // end else resids aren't good enough
  } // end while not done

  //!
  SimpleVector& x_ma57 = *x;
  ma27->matrixChanged();
  ma27->solve(x_ma27);

  for(int i=0; i<x_ma57.length(); i++)
    if(fabs((x_ma57[i]-x_ma27[i])/(1+x_ma27[i])) >1) {
      printf("[%6d]  ma57:%22.14e  ma27:%22.14e \n", i, x_ma57[i], x_ma27[i]);
    }

  SimpleVector resid_ma57(n);
  resid_ma57.copyFrom(rhs_sav);
  sgm_->mult(1.0, resid_ma57.elements(), 1, -1.0, x_ma57.elements(), 1);
  cout << "MA57 resid.nrm=" << resid_ma57.infnorm() << "   rhs.nrm=" << rhs_sav.infnorm()
       << "   n=" << n << endl;


  SimpleVector resid_ma27(n);
  resid_ma27.copyFrom(rhs_sav);
  sgm_->mult(-1.0, resid_ma27.elements(), 1, 1.0, x_ma27.elements(), 1);
  cout << "MA27 resid.nrm=" << resid_ma27.infnorm() << endl;


  rhs.copyFrom( *x );
  //rhs.copyFrom(x_ma27);

  //delete [] iwork; //it is cached now
}

/*void Ma57Solver::Refine( OoqpVector& x_in, OoqpVector& rhs_in )
{
  int job=2; //calculate r=b-Ax, solve A(dx)=r, update solution and exit.

  SimpleVector & rhs = dynamic_cast<SimpleVector &>(rhs_in);
  SimpleVector & x   = dynamic_cast<SimpleVector &>(x_in);

  double * drhs   = rhs.elements();
  double * dx     = x.elements();
  double * dresid = new double[n];
  double * dwork  = new double[5*n];

  icntl[8]=2;//steps of iterative refinement

  int * iwork = new_iworkn(n);
  FNAME(ma57dd)( &job,       &n,        &nnz,   M,        irowM,   jcolM,
	   fact,       &lfact,    ifact,  &lifact,  drhs,    dx,
	   dresid,      dwork,    iwork,  icntl,    cntl,    info, rinfo );

  if(info[0]!=0) cout << "ma57dd: info[0]=: " << info[0] << endl;

  delete[] dwork; delete[] dresid;

  }*/

Ma57Solver::~Ma57Solver()
{
  delete [] irowM;
  delete [] jcolM;
  delete [] fact;
  delete [] ifact;
  delete [] keep;
  if(iworkn) delete[] iworkn;
  if(dworkn) delete[] dworkn;
}

void Ma57Solver::Lsolve( OoqpVector& x )
{
  solve(2,x);
}

void Ma57Solver::Dsolve( OoqpVector& x )
{
  solve(3,x);
}

void Ma57Solver::Ltsolve( OoqpVector& x )
{
  solve(4,x);
}


void Ma57Solver::solve(int solveType, OoqpVector& rhs_in)
{
  //dumpdata(irowM, jcolM, M, n, nnz);
  if(solveType < 1 || solveType > 4) {
    assert("Unknown JOB assigned for use in MA57CD!" && 0);
  } else if(solveType==1) {
    // we prefer iterative refinement.
    solve(rhs_in);
    return;
  } /*else */

  int job = solveType; // Solve using A
  int one = 1;

  SimpleVector & rhs = dynamic_cast<SimpleVector &>(rhs_in);
  //!SimpleVectorHandle work( new SimpleVector(n) );

  double * drhs   = rhs.elements();
  double * dwork  = new_dworkn(n);
  int * iwork     = new_iworkn(n);

#ifdef HAVE_GETRUSAGE
  rusage before;
  getrusage( RUSAGE_SELF, &before );
#endif

  FNAME(ma57cd)( &job,       &n,
	   fact,       &lfact,    ifact,  &lifact,
	   &one,       drhs,      &n,
	   dwork,      &n,        iwork,
	   icntl,      info );
  /*
#ifdef HAVE_GETRUSAGE
  rusage after;
  getrusage( RUSAGE_SELF, &after );
  cout << "Solution with the factored matrix took "
       << (double) (after.ru_utime.tv_sec - before.ru_utime.tv_sec)
    + (after.ru_utime.tv_usec - before.ru_utime.tv_usec) / 1000000.0
       << " seconds.\n";
#endif
  */
  //delete [] iwork; //!using a sort of cache now
}

void Ma57Solver::solve(GenMatrix& rhs_in)
{
  DenseGenMatrix &rhs = dynamic_cast<DenseGenMatrix&>(rhs_in);
  int N,NRHS;
  // rhs vectors are on the "rows", for continuous memory
  rhs.getSize(NRHS,N);
  assert(n==N);


  // we need checks on the residuals, can't do that with multiple RHS
  for (int i = 0; i < NRHS; i++) {
    SimpleVector v(rhs[i],N);
    solve(v);
  }


//   int job = 1;

//   const int BLOCKSIZE = 20;

//   double * dwork  = new_dworkn(n*BLOCKSIZE);
//   int dworksize = n*BLOCKSIZE;
//   int * iwork     = new_iworkn(n);

//   for (int startcol = 0; startcol < NRHS; startcol += BLOCKSIZE) {
//     double *drhs = rhs[startcol];
//     int endcol = MIN(startcol+BLOCKSIZE,NRHS);
//     int numcols = endcol-startcol;
//     //cout << "MA57 multiple RHS" << endl;
//     FNAME(ma57cd)( &job,       &n,
// 		   fact,       &lfact,    ifact,  &lifact,
// 		   &numcols,       drhs,      &n,
// 		   dwork,      &dworksize,        iwork,
// 		   icntl,      info );
//     assert(info[0] >= 0);
//     if (info[0] > 0) {
//       printf("warning from ma57cd, info[0]=%d\n",info[0]);
//     }
//   }

}


int* Ma57Solver::new_iworkn(int dim)
{
  if(niworkn != dim) {
    if(iworkn) delete[] iworkn;
    iworkn = new int[dim];
    niworkn=dim;
  }else {
    if(nullptr == iworkn)
      iworkn = new int[dim];
  }
  return iworkn;
}

double* Ma57Solver::new_dworkn(int dim)
{
  if(ndworkn != dim) {
    if(dworkn) delete[] dworkn;
    dworkn = new double[dim];
    ndworkn = dim;
  }else {
    if(nullptr == dworkn)
      dworkn = new double[dim];
  }
  return dworkn;
}
