/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP 
 * Modefied by Cosmin Petra to perform solves with the factors.
 */

#include <iostream>
using namespace std;

#include "PardisoSolver.h"
#include "SparseStorage.h"
#include "SparseSymMatrix.h"
#include "SimpleVector.h"
#include "SimpleVectorHandle.h"
#include "DenseGenMatrix.h"

#ifdef HAVE_GETRUSAGE
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#endif

extern int gOoqpPrintLevel;


extern "C" void pardisoinit (void   *, int    *,   int *, int *, double *, int *);
extern "C" void pardiso     (void   *, int    *,   int *, int *,    int *, int *, 
                  double *, int    *,    int *, int *,   int *, int *,
                     int *, double *, double *, int *, double *);


extern "C" void pardiso_chkmatrix  (int *, int *, double *, int *, int *, int *);
extern "C" void pardiso_chkvec     (int *, int *, double *, int *);
extern "C" void pardiso_printstats (int *, int *, double *, int *, int *, int *,
                           double *, int *);



//PardisoSolver::PardisoSolver( SparseSymMatrix * sgm )
//{
//	first = true;
// store sgm and call pardisoinit
//}

PardisoSolver::PardisoSolver( SparseSymMatrix * sgm_ )
{
  cout << "Murat Mut Pardiso Solver ***************************************************" << endl;
  first = true;
  // initialize each data element.
  
   pardisoInitMurat(sgm_);
  
  // store sgm and call pardisoinit
  
 
}


// We initialize Pardiso.
void PardisoSolver::pardisoInitMurat(SparseSymMatrix * sgm)
{
    

   sgm_ = sgm;
   
  // We need to initialize pt
  //            mtype
  
  
    n = sgm_->size();
    nnz = sgm_->numberOfNonZeros();
    mStorage = SparseStorageHandle( sgm_->getStorage() ); 
    
    cout << "what is n?" << n<< endl;  //Murat  
    cout << "what is nnz?" << nnz << endl;  //Murat  
    
    
    //krowMt = new int[n+1];
	//jcolMt = new int[nnz];
	//Mt     = new double[nnz];
	
	/* RHS and solution vectors. */
    //  double   b[8], x[8];
    // int      nrhs = 1;          /* Number of right hand sides. */

    /* Internal solver memory pointer pt,                  */
    /* 32-bit: int pt[64]; 64-bit: long int pt[64]         */
    /* or void *pt[64] should be OK on both architectures  */ 
    
    // void    *pt[64]; 
    
    /* Pardiso control parameters. */
    
    
    int      maxfct, mnum, phase, error, msglvl, solver, numRows,  numCols;

    /* Number of processors. */
    int      num_procs;

    /* Auxiliary variables. */
    char    *var;
    int      i;

    double   ddum;              /* Double dummy */
    int      idum;              /* Integer dummy. */
    
    sgm_-> getSize(numRows, numCols);  //Murat
   
    cout << "numRCols = " <<  numCols << endl;  //Murat  
   
    cout << "numRows = "  <<  numRows << endl;            //Murat  
    cout << "what is krow " <<   (mStorage->krowM)[5] << endl;  //Murat
    cout << "what is jrow " <<   mStorage->jcolM[2] << endl;  //Murat
    
   
    error = 0;
    solver = 1;  /* use sparse direct solver */
 
}

// We need to solve Pardiso. 
void PardisoSolver::pardisoMurat()



{
 // pardiso();
}


void PardisoSolver::firstCall()
{

pardisoinit (pt,  &mtype, &solver, iparm, dparm, &error); 

cout << "enters firstCall??? " << endl;  //Murat  


} 

 
void PardisoSolver::diagonalChanged( int /* idiag */, int /* extent */ )
{
  this->matrixChanged();
}

void PardisoSolver::matrixChanged()
{
	if (first) { firstCall(); first = false; }
// compute numerical factorization
}
 
void PardisoSolver::solve( OoqpVector& rhs_in )
{
    
	SimpleVector & rhs = dynamic_cast<SimpleVector &>(rhs_in);
	double * drhs = rhs.elements();
	
    double * rhscpy = new double[n];        //These two lines from Miles in WSMP
    memcpy(rhscpy,drhs,n*sizeof(double));	//Miles in WSMP
	
// solve linear system with rhs_in
}


PardisoSolver::~PardisoSolver()
{

}



void PardisoSolver::solve(GenMatrix& rhs_in)
 {
 
 }



 /*void PardisoSolver::Lsolve( OoqpVector& x )
   {
  
   }  
 */

 /*  void PardisoSolver::Dsolve( OoqpVector& x )
  {
  
  }  */

 /* void PardisoSolver::Ltsolve( OoqpVector& x )
    {
  
    } */



