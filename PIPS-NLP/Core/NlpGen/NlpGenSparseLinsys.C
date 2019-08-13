/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

/* 2015. Modified by Nai-Yuan Chiang for NLP*/

#include "NlpGenSparseLinsys.h"
#include "DoubleLinearSolver.h"
#include "SparseSymMatrix.h"
#include "SparseGenMatrix.h"

#include "RegularizationAlg.h"

#include "NlpGenData.h"
#include "NlpGenVars.h"
#include "SparseGenMatrixHandle.h"

extern int gOuterSolve;
extern int separateHandDiag;
//extern int gAssumeMatSingular;

extern int gUsePetsc;
extern int gUser_Defined_PC;
extern int gUsePetscOuter;
extern int gisNLP;
extern int gPipsPrtLV;

NlpGenSparseLinsys::NlpGenSparseLinsys(  NlpGen * factory_in,
				       NlpGenData * data,
				       LinearAlgebraPackage * la,
				       SparseSymMatrix * Mat_in,
				       DoubleLinearSolver * solver_in )
  : NlpGenLinsys( factory_in, data, la ), solver(solver_in)
{
  SpReferTo( Mat, Mat_in );

  firstXDiagUpdate = true;
  firstSDiagUpdate = true;
  firstYDiagUpdate = true;
  firstZDiagUpdate = true;

  firstQUpdate = true;
  firstAUpdate = true;
  firstCUpdate = true;

}


void NlpGenSparseLinsys::putXDiagonal( OoqpVector& xdiag )
{
  Mat->atPutDiagonal( 0, xdiag );
}

void NlpGenSparseLinsys::putSDiagonal( OoqpVector& sdiag )
{
  Mat->atPutDiagonal( nx, sdiag );
}

void NlpGenSparseLinsys::putYDualDiagonal( OoqpVector& ydiag )
{
  if(gOuterSolve < 3){
  	Mat->atPutDiagonal( nx, ydiag );
  }else if(gOuterSolve >= 3 && separateHandDiag==0){
	Mat->atPutDiagonal( nx + mz, ydiag );
  }else{
	assert(0);
  }
  
}

void NlpGenSparseLinsys::putZDiagonal( OoqpVector& zdiag )
{
  if(gOuterSolve < 3){
  	Mat->atPutDiagonal( nx + my, zdiag );
  }else if(gOuterSolve >= 3 && separateHandDiag==0){
	Mat->atPutDiagonal( nx + mz + my, zdiag );
  }else{
	assert(0);
  }
}

void NlpGenSparseLinsys::setXDiagonal( OoqpVector& xdiag )
{
  Mat->copyDiagonalVal_From( 0, xdiag, firstXDiagUpdate, xDiagIdxMap);
  firstXDiagUpdate = false;
}

void NlpGenSparseLinsys::setSDiagonal( OoqpVector& sdiag )
{
  Mat->copyDiagonalVal_From( nx, sdiag, firstSDiagUpdate, sDiagIdxMap);
  firstSDiagUpdate = false;
}

void NlpGenSparseLinsys::setYDiagonal( OoqpVector& ydiag )
{
  if(gOuterSolve < 3){
  	Mat->copyDiagonalVal_From( nx, ydiag, firstYDiagUpdate, yDiagIdxMap);
  }else if(gOuterSolve >= 3 && separateHandDiag==0){
	Mat->copyDiagonalVal_From( nx + mz, ydiag, firstYDiagUpdate, yDiagIdxMap );
  }else{
	assert(0);
  }
  firstYDiagUpdate = false;
}

void NlpGenSparseLinsys::setZDiagonal( OoqpVector& zdiag )
{
  if(gOuterSolve < 3){
  	Mat->copyDiagonalVal_From( nx + my, zdiag, firstZDiagUpdate, zDiagIdxMap );
  }else if(gOuterSolve >= 3 && separateHandDiag==0){
	Mat->copyDiagonalVal_From( nx + mz + my, zdiag, firstZDiagUpdate, zDiagIdxMap );
  }else{
	assert(0);
  }
  firstZDiagUpdate = false;
}


void NlpGenSparseLinsys::setAdditiveDiagonal()
{
  assert(gOuterSolve >= 3 && separateHandDiag==1);
  Mat->setAdditiveDiagonal(*additiveDiag);
}


void NlpGenSparseLinsys::solveCompressed( OoqpVector & arhs )
{
  solver->solve( arhs );
  KryIter=solver->KryIter;
}  

NlpGenSparseLinsys::~NlpGenSparseLinsys()
{
  delete solver;
}


void NlpGenSparseLinsys::factor(Data *prob, Variables *vars)
{
  factorNoMatChange(prob, vars,NULL);
  solver->matrixChanged();  
}

void NlpGenSparseLinsys::factor(Data *prob_in, Variables *vars,RegularizationAlg *RegInfo)
{
  int Num_NegEVal = -1;
  bool skipUpdateReg = false;
  double priReg=0.0, dualReg=0.0;

  NlpGenData * prob = (NlpGenData *) prob_in;
  
  // DoEvalReg =   1 -> when factorizing the matrix, add regularizations to correct inertia and singularity 
  //				  	(we ONLY calling this routine once in IBR)
  //			  2 -> when factorizing the matrix, add regularizations to correct singularity only 
  //					(this is always the 1st call of this routine when IFR is used)
  //			  0 -> when factorizing the matrix, force to use primal regularizaion. called iff xWx tests fail  
  //					(the other calls of this routine when IFR is used, now matrix is nonsingular for sure)
  if(RegInfo->DoEvalReg >= 1){
    RegInfo->newLinearSystem(); 
    
    if(RegInfo->ForceReg)
      factorNoMatChange(prob, vars,RegInfo);
    else
      factorNoMatChange(prob, vars,NULL);
    
    Num_NegEVal = solver->matrixChanged();
    
    if(gPipsPrtLV>=3) 
      printf("NlpGenSparseLinsys (serial): Num_NegEVal is %d and  my+mz is %lld\n", Num_NegEVal, my+mz);

    // check if matrix is singular
    if( Num_NegEVal < 0 || (Num_NegEVal < my+mz && RegInfo->DoEvalReg == 1) )
      RegInfo->MatrixSingular = 1;
    else
      RegInfo->MatrixSingular = 0;
    
    // skip update regularization if: 	1) have correct inertia and mat is nonsingular 
    //				OR 	2) mat is nonsingular and we will do inertia-free test later
    if( (RegInfo->DoEvalReg == 1 && Num_NegEVal == my+mz) || (RegInfo->DoEvalReg == 2 && Num_NegEVal != -1) ) 
      skipUpdateReg = true;	
  }
  
  // update regularization
  while(!skipUpdateReg){
    RegInfo->computeRegularization(priReg,dualReg,prob->currMu);
	
    factorNoMatChange(prob, vars,RegInfo);
    Num_NegEVal = solver->matrixChanged();	

    if(gPipsPrtLV>=3) 
      printf("NlpGenSparseLinsys (serial): Num_NegEVal is %d and  my+mz is %lld (on refactorization with regularizations: %g %g)\n", Num_NegEVal, my+mz, priReg, dualReg);

    // check if matrix is singular
    if(Num_NegEVal < 0)
      RegInfo->MatrixSingular = 1;
    else
      RegInfo->MatrixSingular = 0;
    
    // skip update regularization if: 	1) have correct inertia and mat is nonsingular 
    //				OR 	2) mat is nonsingular and we will do inertia-free test later	
    //				OR  	3) we are doing inertia-free test now (mat is definitely nonsingular)
    if( (RegInfo->DoEvalReg == 1 && Num_NegEVal == my+mz) 
	|| (RegInfo->DoEvalReg == 2 && Num_NegEVal != -1) || RegInfo->DoEvalReg == 0 )
      {
	skipUpdateReg = true;
      }	
  } 
}

//FIXME_ME: if we do reduced space, we do not need to updae kkt matrices
void
NlpGenSparseLinsys::UpdateMatrices( Data * prob_in, int const updateLevel)
{
	int useUpdate=updateLevel;
	if(!gisNLP) useUpdate=1;

	NlpGenData * prob = (NlpGenData *) prob_in;

	if(useUpdate>=2){
	  if(gOuterSolve < 3){
	    if(fullQ)
		  prob->setQIntoAt( *Mat, 0, 0, firstQUpdate, QmatIdxMap);
	    prob->setAIntoAt( *Mat, nx, 0, firstAUpdate, AmatIdxMap);
	    prob->setCIntoAt( *Mat, nx + my, 0, firstCUpdate, CmatIdxMap);
	  }else{
		if(fullQ)
		  prob->setQIntoAt( *Mat, 0, 0,  firstQUpdate, QmatIdxMap);
		prob->setAIntoAt( *Mat, nx + mz, 0, firstAUpdate, AmatIdxMap);
		prob->setCIntoAt( *Mat, nx + mz +my, 0,firstCUpdate, CmatIdxMap);
	  }
	  firstQUpdate = false;
	  firstAUpdate = false;
	  firstCUpdate = false;
	}

	if(useUpdate>=1){
	  if(gOuterSolve < 3){
	    if( nxlow + nxupp >= 0 ) setXDiagonal( *dd );
	    if( my > 0) setYDiagonal( *temp_diagY);	  	
	    if( mclow + mcupp > 0 ) setZDiagonal( *nomegaInv );
	  }else if (gOuterSolve >= 3 && separateHandDiag==1){
		joinRHSXSYZ(*additiveDiag,*dd,*temp_diagS,*temp_diagY ,*temp_diagZ);
		setAdditiveDiagonal();
	  }else if(gOuterSolve >= 3 && separateHandDiag==0){
		if( nxlow + nxupp >= 0 ) setXDiagonal( *dd );
	    if( my > 0) setYDiagonal( *temp_diagY);	  
		if( mclow + mcupp > 0 ){
			setSDiagonal( *temp_diagS );
			setZDiagonal( *temp_diagZ );			
		}
	  }else{
		assert(0);
	  }

	}
	
}


