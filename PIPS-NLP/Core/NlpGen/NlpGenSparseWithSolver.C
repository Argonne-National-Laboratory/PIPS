/* PIPS-NLP                                                         	*
 * Authors: Nai-Yuan Chiang                      		*
 * (C) 2015 Argonne National Laboratory			*/

#include "NlpGenSparseWithSolver.h"
#include "NlpGenSparseLinsys.h"

#include "NlpGenData.h"

#include "SparseLinearAlgebraPackage.h"
#include "SparseSymMatrix.h"
#include "SimpleVector.h"
#include "SimpleVectorHandle.h"

#include "SparseLinearAlgebraPackage.h"

#include "DoubleLinearSolver.h"


#ifdef WITH_MA27
#include "Ma27Solver.h"
#endif
#ifdef WITH_MA57
#include "Ma57Solver.h"
#endif
#ifdef WITH_PARDISO
#include "PardisoSolver.h"
#endif 
#ifdef WITH_UMFPACK
#include "UmfPackSolver.h"
#endif

//#include "MtxSchurDecompSolver.h"
//#include "ReducedSpaceSolver.h"

#ifdef WITH_PETSC
#include "PetscIterativeSolver_Sparse.h"
#endif


extern int gOuterSolve;
extern int separateHandDiag;
extern int gSymLinearSolver;

extern int gBuildSchurComp;
extern int gUseReducedSpace;

extern int gUsePetsc;
extern int gUsePetscOuter;
extern int gUser_Defined_PC;


NlpGenSparseWithSolver::NlpGenSparseWithSolver( int nx, int my, int mz,
					  int nnzQ, int nnzA, int nnzC ) :
  NlpGenSparse( nx, my, mz, nnzQ, nnzA, nnzC )
{
  la = SparseLinearAlgebraPackage::soleInstance();
  fullQ = ((gUsePetsc==0) || (gUsePetscOuter==0) || (gUser_Defined_PC!=2));
}

LinearSystem * NlpGenSparseWithSolver::makeLinsys( Data * prob_in,int setMulti)
{
  assert("Not Implemented" &&0);
}


LinearSystem * NlpGenSparseWithSolver::makeLinsys( Data * prob_in )
{
  NlpGenData * prob = (NlpGenData *) prob_in;
  int n, nnzTotal;
  int ns = mz;
  DoubleLinearSolver * solver = NULL;
  LinearSystem *resultLS = NULL;

	if(gOuterSolve >= 3 && separateHandDiag==0){
      n = nx + ns + my + mz;    // for x s y z
	  nnzTotal = n + nnzQ + nnzA + nnzC + ns;
	}
	else if(gOuterSolve >= 3 && separateHandDiag==1){
	  n = nx + ns + my + mz;	  // for x s y z
	  nnzTotal = nnzQ + nnzA + nnzC + ns;
	}
	else {
	  n = nx + my + mz;			// only x y z, we need compress the linear system later
	  nnzTotal = n + nnzQ + nnzA + nnzC;
	}
	
	SparseSymMatrixHandle Mat( new SparseSymMatrix( n, nnzTotal ) );
	
	if(gOuterSolve >= 3 ){
	  if(separateHandDiag==0){
	    SimpleVectorHandle v( new SimpleVector(n) );
	    v->setToZero();
	    Mat->setToDiagonal(*v);
	  }

	  if(fullQ) prob->putQIntoAt( *Mat, 0, 0 );
	  prob->putAIntoAt( *Mat, nx + ns, 0);
	  	  
	  if(ns>0){
		prob->putCIntoAt( *Mat, nx + ns + my, 0 ); 
		SparseSymMatrixHandle tempDiagMat( new SparseSymMatrix( ns, ns ) );
		int *tempDiagRowId = new int[ns];
		int *tempDiagColId = new int[ns];
		double *tempDiagEleId = new double[ns];
		int info;
		
		for(int i=0;i<ns;i++){
		  tempDiagRowId[i]=i;
		  tempDiagColId[i]=i;
		  tempDiagEleId[i]=-1;
		}
		tempDiagMat->putSparseTriple( tempDiagRowId,ns, tempDiagColId, tempDiagEleId,info );
		Mat->symAtPutSubmatrix( nx + ns + my, nx, *tempDiagMat, 0, 0, ns, ns);
		
        delete[] tempDiagRowId;
     	delete[] tempDiagColId;
     	delete[] tempDiagEleId;
      }
	}
	else{
	  SimpleVectorHandle v( new SimpleVector(n) );
	  v->setToZero();
	  Mat->setToDiagonal(*v);

	  prob->putQIntoAt( *Mat, 0, 0 );
	  prob->putAIntoAt( *Mat, nx, 0);
	  prob->putCIntoAt( *Mat, nx + my, 0 );	
	}



  if(0==gUseReducedSpace){
	if(1==gBuildSchurComp){
	  if(1==gUsePetsc && gUsePetscOuter==0){
#ifdef WITH_PETSC
		solver = new PetscIterativeSolver_Sparse( Mat,my+mz);
#endif
	  }
      else if(0==gSymLinearSolver){
#ifdef WITH_MA27
		solver = new Ma27Solver( Mat);
#endif 
      }
      else if(1==gSymLinearSolver){
#ifdef WITH_MA57
		solver = new Ma57Solver( Mat);
#endif 
      }  
      else if(2==gSymLinearSolver){
#ifdef WITH_PARDISO
		solver = new PardisoSolver( Mat,my+mz);
#endif 
      }
	  else if(3==gSymLinearSolver){
#ifdef WITH_UMFPACK
		solver = new UmfPackSolver( Mat,my+mz);
#endif 
      }	
    }
    else{
	    // we only support dense schur comp
	}

  }
  else if (1==gUseReducedSpace){
	assert(0==gBuildSchurComp);
	if(0==gBuildSchurComp){
	  int decisionVarSize = prob->schurSize;
	  int *decisionVarID = prob->schurVarConID;

	  int fullVarXSize = nx;
	  int fullVarYSize = my;
	  int fullVarSSize = mz;

//	  solver = new ReducedSpaceSolver(Mat,decisionVarSize,decisionVarID,fullVarXSize,fullVarYSize,fullVarSSize);
	}
  }	

  if(solver == NULL)		
    assert("solver invalid" && 0);

  resultLS = new NlpGenSparseLinsys( this, prob, la, Mat, solver );


  assert(resultLS != NULL);

  return resultLS;

}




