/* PIPS-NLP                                                         	*
 * Authors: Nai-Yuan Chiang                      		*
 * (C) 2015 Argonne National Laboratory			*/

#include "MtxSchurDecompSolver.h"

#include "assert.h"
#include "stdlib.h"

#include "DoubleMatrix.h"
#include "SparseSymMatrix.h"

#include "SparseSymMatrix.h"
#include "SparseGenMatrix.h"

#include "RegularizationAlg.h"
#include "NlpGenLinsys.h"

#ifdef WITH_MA27
#include "Ma27Solver.h"
#endif 
#ifdef WITH_MA57
#include "Ma57Solver.h"
#endif 
#ifdef WITH_PARDISO
#include "PardisoSolver.h"
#endif 

#include "DeSymIndefSolver.h"
#include "sLinsysRoot.h"

#include "OoqpVector.h"
#include "SimpleVector.h"

extern int gBuildSchurComp;
extern int gSymLinearSolver;
extern int gSolveSchurScheme;

#ifndef MIN
#define MIN(a,b) ((a > b) ? b : a)
#endif

#ifdef TIMING
  #include "mpi.h"
  extern double probGenTime;
#endif

MtxSchurDecompSolver::MtxSchurDecompSolver()
	: 
	 requireUpdate(1), firstCallFlag(1), mtxCase(0),
	 borderDim_m(0), borderDim_n(0), diag0MatDim(0), diag1MatDim(0), fullMatDim(0),
	 borderMatNNz(0), diag0MatNNz(0), diag1MatNNz(0), fullMatNNz(0),	 
	 borderMat_ele(NULL), diag0Mat_ele(NULL), diag1Mat_ele(NULL), 
	 borderMat_rowBeg(NULL), diag0Mat_rowBeg(NULL), diag1Mat_rowBeg(NULL), 
	 borderMat_colIdx(NULL), diag0Mat_colIdx(NULL), diag1Mat_colIdx(NULL), 
	 schurMatDim(0), schurMatNNz(0),	 
	 schurMat_ele(NULL), Msys(NULL), solver_schur(NULL),solver_diag1(NULL),
	 //kktSC(NULL),
	 dkktSC(NULL)
{}

MtxSchurDecompSolver::MtxSchurDecompSolver(DoubleMatrix* MatIn, const int schurDim,
					int *inputMatType, int* schurVarIDX, const int localNegaEigVal_in)
	: 
	 requireUpdate(1), firstCallFlag(1), mtxCase(0), 
	 borderDim_m(0), borderDim_n(0), diag0MatDim(0), diag1MatDim(0), fullMatDim(0),
	 borderMatNNz(0), diag0MatNNz(0), diag1MatNNz(0), fullMatNNz(0),	 
	 borderMat_ele(NULL), diag0Mat_ele(NULL), diag1Mat_ele(NULL), 
	 borderMat_rowBeg(NULL), diag0Mat_rowBeg(NULL), diag1Mat_rowBeg(NULL), 
	 borderMat_colIdx(NULL), diag0Mat_colIdx(NULL), diag1Mat_colIdx(NULL), 
	 schurMatDim(0), schurMatNNz(0),	 
	 schurMat_ele(NULL), Msys(NULL), solver_schur(NULL),solver_diag1(NULL),
	 //kktSC(NULL),
	 dkktSC(NULL)
{
  newMtxSchurDecompSolver(MatIn, schurDim,inputMatType);
  schurVarIDinFull  = schurVarIDX;
  localNegaEigVal = localNegaEigVal_in;
}

void MtxSchurDecompSolver::newMtxSchurDecompSolver(DoubleMatrix* MatIn, const int schurDim,
					int *inputMatType)
{
  MatIn->getSize(fullMatDim,fullMatDim);
  
  borderDim_m = fullMatDim - schurDim;
  borderDim_n = schurDim;
  
  diag1MatDim = fullMatDim - schurDim;

  diag0MatDim = schurDim;  
  
  schurMatDim = schurDim;
  


  
  mtxCase = inputMatType[0]+inputMatType[1]+inputMatType[2]+inputMatType[3];
  assert(mtxCase == INPUT_RLSS);

//  if(gBuildSchurComp == 1)
  {
  	schurMatNNz  = schurMatDim*schurMatDim;
  }
  
  if(inputMatType[1] == SCHURSOLVER_INPUT_SPARSE){
	if(inputMatType[0] == SCHURSOLVER_INPUT_SYMMETRIC){
	  fullMatNNz = ((SparseSymMatrix*)MatIn)->numberOfNonZeros();
	  Msys = dynamic_cast<SparseSymMatrix*>(MatIn);
	}else{
	  assert("impossible!" && 0);
	}
  }else if(inputMatType[1] == SCHURSOLVER_INPUT_DENSE) {
  	if(inputMatType[2] == SCHURSOLVER_INPUT_FULL){
	  fullMatNNz = fullMatDim*fullMatDim;
	}else{
	  assert("impossible!" && 0);
	}
	assert("impossible!" && 0);
  }else{
	  assert("impossible!" && 0);
  }

  dkktSC =  new DenseSymMatrix(schurMatDim);
  solver_schur =  dynamic_cast<DeSymIndefSolver*> (new DeSymIndefSolver(dkktSC));

}



int _findIDX(const int idinFullMat, const int vlength, int *_vector)
{
  int findIDX = -1;
  for(int i=0;i<vlength;i++){
  	if(_vector[i]==idinFullMat){
	  findIDX = i;
	  break;
	}
  }
  return findIDX; 	  
}

int _findIDXinDiag1(const int idinFullMat, const int vlength, int *_vector)
{
  int findIDX = idinFullMat;
  for(int i=0;i<vlength;i++)
  	if(_vector[i]<idinFullMat)
	  findIDX--;
  return findIDX; 	  
}

void MtxSchurDecompSolver::_firstCall_RLSS()
{
  firstCallFlag = 0;

  //save the indeces for diagonal entries for a streamlined later update
  int* krowbegFullMat 	= Msys->getStorageRef().krowM;
  int* jcolFullMat 		= Msys->getStorageRef().jcolM;
  double* eleFullMat 	= Msys->getStorageRef().M;

  diag0Mat_rowBeg  = (int*) malloc((diag0MatDim+1)*sizeof(int));
  diag1Mat_rowBeg  = (int*) malloc((diag1MatDim+1)*sizeof(int));
  borderMat_rowBeg = (int*) malloc((borderDim_m+1)*sizeof(int));
  for(int jj=0;jj<diag0MatDim+1;jj++) diag0Mat_rowBeg[jj]=0;
  for(int jj=0;jj<diag1MatDim+1;jj++) diag1Mat_rowBeg[jj]=0;
  for(int jj=0;jj<borderDim_m+1;jj++) borderMat_rowBeg[jj]=0;

  diag1VarIDinFull = (int*) malloc(diag1MatDim*sizeof(int));
  FullVarIDinDiag01 = (int*) malloc(fullMatDim*sizeof(int));
  for(int currVar=0; currVar<fullMatDim;currVar++){
  	FullVarIDinDiag01[currVar]=0;
  }

  // build var local-full map and build var  full-local map
  for(int currVar=0; currVar<fullMatDim;currVar++){
	int localVarIDInSChur = _findIDX(currVar,schurMatDim,schurVarIDinFull);
	if(localVarIDInSChur == -1){
	  //this var is in diag1
	  localVarIDInSChur = _findIDXinDiag1(currVar,schurMatDim,schurVarIDinFull);
	  assert(localVarIDInSChur<=currVar);
	  diag1VarIDinFull[localVarIDInSChur] = currVar;
	  FullVarIDinDiag01[currVar] =  (localVarIDInSChur+1);
	}else{
	  //this var is in diag0
	  FullVarIDinDiag01[currVar] = -(localVarIDInSChur+1);
	}
  }


  
  for(int currRow=0; currRow<fullMatDim;currRow++){
	//count 1st, 2nd and link size
	
	for( int k=krowbegFullMat[currRow];k<krowbegFullMat[currRow+1];k++){	
	  int currCol = jcolFullMat[k];

//	  int localColIDInSChur = _findIDX(currCol,schurMatDim,schurVarIDinFull);
//	  int localRowIDInSChur = _findIDX(currRow,schurMatDim,schurVarIDinFull);
	  int localColIDInSChur = FullVarIDinDiag01[currCol];
	  int localRowIDInSChur = FullVarIDinDiag01[currRow];
	  
	  if(localColIDInSChur < 0 && localRowIDInSChur < 0){
		//belong to Diag0, correct the neg index
		localColIDInSChur = -(localColIDInSChur+1);
		localRowIDInSChur = -(localRowIDInSChur+1);
		diag0MatNNz++;
		if(localRowIDInSChur >= localColIDInSChur)
		  diag0Mat_rowBeg[localRowIDInSChur+1]++;
		else
		  diag0Mat_rowBeg[localColIDInSChur+1]++;
	  }else if(localColIDInSChur > 0 && localRowIDInSChur >0){
		//belong to Diag1,correct the pos index
		localColIDInSChur = (localColIDInSChur-1);
		localRowIDInSChur = (localRowIDInSChur-1);
		diag1MatNNz++;
		if(localRowIDInSChur >= localColIDInSChur)
		  diag1Mat_rowBeg[localRowIDInSChur+1]++;
		else
		  diag1Mat_rowBeg[localColIDInSChur+1]++;	
	  }else if(localColIDInSChur < 0 && localRowIDInSChur > 0 ){
		//belong to borderMat
		localRowIDInSChur = (localRowIDInSChur-1);
		borderMatNNz++;
		borderMat_rowBeg[localRowIDInSChur+1]++;	  
	  }else if(localColIDInSChur > 0 && localRowIDInSChur < 0){
		//belong to borderMat
		localColIDInSChur = (localColIDInSChur-1);
		borderMatNNz++;
		borderMat_rowBeg[localColIDInSChur+1]++;	  
	  }else{assert("impossible"&&0);}
	}
  }

  assert(borderMatNNz+diag1MatNNz+diag0MatNNz==fullMatNNz);
  
  for(int j=1; j < diag0MatDim+1; j++){
	diag0Mat_rowBeg[j] += diag0Mat_rowBeg[j-1];
  }  
  for(int j=1; j < diag1MatDim+1; j++){
	diag1Mat_rowBeg[j] += diag1Mat_rowBeg[j-1];
  }
  for(int j=1; j < borderDim_m+1; j++){
	borderMat_rowBeg[j] += borderMat_rowBeg[j-1];
  }

  // alocate space
  borderMat_colIdx 	= (int*)malloc(borderMatNNz*sizeof(int));
  diag0Mat_colIdx  	= (int*)malloc(diag0MatNNz*sizeof(int));  
  diag1Mat_colIdx  	= (int*)malloc(diag1MatNNz*sizeof(int));
  
  borderMat_ele 	= (double*)malloc(borderMatNNz*sizeof(double));
  diag0Mat_ele 		= (double*)malloc(diag0MatNNz*sizeof(double));
  diag1Mat_ele 		= (double*)malloc(diag1MatNNz*sizeof(double));  

  Border_Full_eleMap 	= (int*)malloc(borderMatNNz*sizeof(int));
  Diag0_Full_eleMap  	= (int*)malloc(diag0MatNNz*sizeof(int));  
  Diag1_Full_eleMap  	= (int*)malloc(diag1MatNNz*sizeof(int));




  int *nextdiag0InRow  = (int*)malloc(diag0MatDim*sizeof(int));
  int *nextdiag1InRow  = (int*)malloc(diag1MatDim*sizeof(int));
  int *nextborderInRow = (int*)malloc(borderDim_m*sizeof(int));

  for(int jj=0;jj<diag0MatDim;jj++) nextdiag0InRow[jj] = diag0Mat_rowBeg[jj];
  for(int jj=0;jj<diag1MatDim;jj++) nextdiag1InRow[jj] = diag1Mat_rowBeg[jj];
  for(int jj=0;jj<borderDim_m;jj++) nextborderInRow[jj]= borderMat_rowBeg[jj];

  int diag0nnz_wrk=0, diag1nnz_wrk=0, bordernnz_wrk=0; 
  int wrkGoff=0;


  for(int currRow=0; currRow<fullMatDim;currRow++){
	for( int k=krowbegFullMat[currRow];k<krowbegFullMat[currRow+1];k++){	  
	  int currCol = jcolFullMat[k];
//	  int localColIDInSChur = _findIDX(currCol,schurMatDim,schurVarIDinFull);
//	  int localRowIDInSChur = _findIDX(currRow,schurMatDim,schurVarIDinFull);
	  int localColIDInSChur = FullVarIDinDiag01[currCol];
	  int localRowIDInSChur = FullVarIDinDiag01[currRow];
	  
	  if(localColIDInSChur < 0 && localRowIDInSChur < 0){
		//belong to Diag0
		localColIDInSChur = -(localColIDInSChur+1);
		localRowIDInSChur = -(localRowIDInSChur+1);
		if(localRowIDInSChur >= localColIDInSChur){
		  wrkGoff = nextdiag0InRow[localRowIDInSChur]++;
		  diag0Mat_colIdx[wrkGoff] 	= localColIDInSChur;
		  diag0Mat_ele[wrkGoff]		= eleFullMat[k];		  
		}
		else{
		  wrkGoff = nextdiag0InRow[localColIDInSChur]++;
		  diag0Mat_colIdx[wrkGoff] = localRowIDInSChur;
		  diag0Mat_ele[wrkGoff]		= eleFullMat[k];
		}
		Diag0_Full_eleMap[wrkGoff] = k;
		diag0nnz_wrk++;		
	  }
	  else if(localColIDInSChur > 0 && localRowIDInSChur > 0){
		//belong to Diag1
//		localColIDInSChur = _findIDXinDiag1(currCol,schurMatDim,schurVarIDinFull);
//		localRowIDInSChur = _findIDXinDiag1(currRow,schurMatDim,schurVarIDinFull);
//		localColIDInSChur = _findIDX(currCol,diag1MatDim,diag1VarIDinFull);
//		localRowIDInSChur = _findIDX(currRow,diag1MatDim,diag1VarIDinFull);
		localColIDInSChur = (localColIDInSChur-1);
		localRowIDInSChur = (localRowIDInSChur-1);

		if(localRowIDInSChur >= localColIDInSChur){
		  wrkGoff = nextdiag1InRow[localRowIDInSChur]++;
		  diag1Mat_colIdx[wrkGoff] = localColIDInSChur;
		  diag1Mat_ele[wrkGoff]		= eleFullMat[k];
		}	
		else{
		  wrkGoff = nextdiag1InRow[localColIDInSChur]++;
		  diag1Mat_colIdx[wrkGoff] = localRowIDInSChur;
		  diag1Mat_ele[wrkGoff]		= eleFullMat[k];
		}
		Diag1_Full_eleMap[wrkGoff] = k;
		diag1nnz_wrk++;
	  }
	  else if(localColIDInSChur < 0 && localRowIDInSChur > 0){
		//belong to borderMat
//		localRowIDInSChur = _findIDXinDiag1(currRow,schurMatDim,schurVarIDinFull);
//		localRowIDInSChur = _findIDX(currRow,diag1MatDim,diag1VarIDinFull);		
		localColIDInSChur = -(localColIDInSChur+1);
		localRowIDInSChur = (localRowIDInSChur-1);

		wrkGoff = nextborderInRow[localRowIDInSChur]++;
		borderMat_colIdx[wrkGoff] 	= localColIDInSChur;
		borderMat_ele[wrkGoff]		= eleFullMat[k];
		Border_Full_eleMap[wrkGoff] = k;
		bordernnz_wrk++;	  
	  }else if(localColIDInSChur > 0 && localRowIDInSChur < 0){
		//belong to borderMat
//		localColIDInSChur = _findIDXinDiag1(currCol,schurMatDim,schurVarIDinFull);
//		localColIDInSChur = _findIDX(currCol,diag1MatDim,diag1VarIDinFull);
		localColIDInSChur = (localColIDInSChur-1);
		localRowIDInSChur = -(localRowIDInSChur+1);
		
		wrkGoff = nextborderInRow[localColIDInSChur]++;
		borderMat_colIdx[wrkGoff] 	= localRowIDInSChur;
		borderMat_ele[wrkGoff]		= eleFullMat[k];
		Border_Full_eleMap[wrkGoff] = k;
		bordernnz_wrk++;		  
	  }else{assert("impossible"&&0);}
	}
  }
  
  for(int jj=0;jj<diag0MatDim;jj++) 
	assert(nextdiag0InRow[jj]==diag0Mat_rowBeg[jj+1]);
  for(int jj=0;jj<diag1MatDim;jj++) 
	assert(nextdiag1InRow[jj]==diag1Mat_rowBeg[jj+1]);	
  for(int jj=0;jj<borderDim_m;jj++) 
    assert(nextborderInRow[jj]==borderMat_rowBeg[jj+1]);
  assert(bordernnz_wrk == borderMatNNz &&  diag1nnz_wrk == diag1MatNNz && diag0nnz_wrk == diag0MatNNz);	


  diag0Mat = new SparseSymMatrix( diag0MatDim, diag0MatNNz, 
  					diag0Mat_rowBeg, diag0Mat_colIdx, diag0Mat_ele);

  diag1Mat = new SparseSymMatrix( diag1MatDim, diag1MatNNz, 
				  diag1Mat_rowBeg, diag1Mat_colIdx, diag1Mat_ele);

  borderMat = new SparseGenMatrix( borderDim_m, borderDim_n, borderMatNNz, 
				  borderMat_rowBeg, borderMat_colIdx, borderMat_ele);


  solver_diag1 = NULL;
  if(0==gSymLinearSolver){
#ifdef WITH_MA27  	
	solver_diag1	= new Ma27Solver( diag1Mat );
#endif
  }
  else if(1==gSymLinearSolver){
#ifdef WITH_MA57  	
	solver_diag1	= new Ma57Solver( diag1Mat );
#endif
  }
  else if(2==gSymLinearSolver){
#ifdef WITH_PARDISO
	solver_diag1	= new PardisoSolver( diag1Mat, localNegaEigVal );
#endif 	
  } 
  assert(solver_diag1);

  free(nextdiag0InRow);
  free(nextdiag1InRow);
  free(nextborderInRow);

  
}


//faster than DenseSymMatrix::atPutZeros
void MtxSchurDecompSolver::initializeKKT_Dense(DenseSymMatrix* dkktSC)
{
  int nn = dkktSC->size();

  double ** M = dkktSC->getStorageRef().M;

  for(int j=0; j<nn; j++) {
      M[0][j] = 0.0;
  }

  int nToCopy = nn*sizeof(double);

  for(int i=1; i<nn; i++) {
    memcpy(M[i], M[0], nToCopy);
  }
}

#if 0
//faster than DenseSymMatrix::atPutZeros
void MtxSchurDecompSolver::myAtPutZeros(DenseSymMatrix* mat, 
			       int row, int col, 
			       int rowExtent, int colExtent)
{
  int nn = mat->size();
  assert( row >= 0 && row + rowExtent <= nn );
  assert( col >= 0 && col + colExtent <= nn );

  double ** M = mat->getStorageRef().M;

  for(int j=col; j<col+colExtent; j++) {
      M[row][j] = 0.0;
  }

  int nToCopy = colExtent*sizeof(double);

  for(int i=row+1; i<row+rowExtent; i++) {
    memcpy(M[i]+col, M[row]+col, nToCopy);
  }
}

void MtxSchurDecompSolver::myAtPutZeros(DenseSymMatrix* mat)
{
  int n = mat->size();
  myAtPutZeros(mat, 0, 0, n, n);
}
#endif

/**
 * Computes C = Border^T * inv(Diag_1) * Border.
 */

void MtxSchurDecompSolver::addTermToDenseSchurCompl(SparseSymMatrix *DiagMat, 
										SparseGenMatrix* BordMat, DenseSymMatrix& SC) 
{

  int bord_m, bord_n, sc_mn;
  
  BordMat->getSize(bord_m, bord_n); 
  sc_mn = SC.size(); assert(sc_mn>=bord_n);

  if(bord_n==-1) bord_n = sc_mn;
  
//  N = borderDim_m;
  assert(bord_n==sc_mn);
  
  int blocksize = 64;
  DenseGenMatrix cols(blocksize,bord_m);


  for (int it=0; it < bord_n; it += blocksize) {
    int start=it;
    int end = MIN(it+blocksize,bord_n);
    int numcols = end-start;
    cols.getStorageRef().m = numcols; // avoid extra solves


    bool allzero = true;
    memset(&cols[0][0],0,bord_m*blocksize*sizeof(double));

	BordMat->getStorageRef().fromGetColBlock(start, &cols[0][0], bord_m, numcols, allzero);
//    assert(0);
    if(!allzero) {
	  solver_diag1->solve(cols);
	  BordMat->getStorageRef().transMultMat( 1.0, SC[start], numcols, bord_n, -1.0, &cols[0][0], bord_m);
    } //end !allzero
  }
}


/**
 * z0 -= Border^T  * Li\Di\ [ zi ]
 *
 * 
 */
void MtxSchurDecompSolver::addLnizi(SparseGenMatrix *borderMat,OoqpVector& z0_, OoqpVector& zi_)
{
  SimpleVector& z0 = dynamic_cast<SimpleVector&>(z0_);
  SimpleVector& zi = dynamic_cast<SimpleVector&>(zi_);

  solver_diag1->Dsolve (zi);  
  solver_diag1->Ltsolve(zi);

  //get n0= nx(parent)= #cols of A or C
  int dummy, n0;
  borderMat->getSize(dummy, n0);

  borderMat->transMult(1.0, z0, -1.0, zi);
}


/*
 *  y = alpha*Lni^T x + y
 *
 *                                  ( [ R 0 0 ]     )
 *  y = y + alpha* Di\Li\ (  [ A 0 0 ] * x )
 *                                  (  [ C 0 0 ]    )
 *
*/

void MtxSchurDecompSolver::LniTransMult(SparseGenMatrix *borderMat, 
									SimpleVector& y, double alpha, SimpleVector& x)
{
  int N, nx0;
  
  //get nx(parent) from the number of cols of A (or C). Use N as dummy
  borderMat->getSize(N,nx0);
  // a mild assert
  assert(nx0 <= x.length());

  N = diag1MatDim;
  assert( y.length() == N);

  //!memopt
  SimpleVector LniTx(N);
  SimpleVector x1(&x[0], nx0);
  
  LniTx.setToZero();
  borderMat->mult(0.0, LniTx, 1.0, x1);
  
  solver_diag1->Lsolve(LniTx); 
  solver_diag1->Dsolve(LniTx);
  solver_diag1->Ltsolve(LniTx);
  
  y.axpy(alpha,LniTx); 
 
}


void MtxSchurDecompSolver::finalizeKKT(SparseSymMatrix *diag0Mat_, DenseSymMatrix * kktd )
{
  int j, p, pend; double val;

  //alias for internal buffer of kkt
  double** dKkt = kktd->Mat();

  int SCdim = kktd->size();
  
  /////////////////////////////////////////////////////////////
  // update the KKT with Diag0 
  /////////////////////////////////////////////////////////////
  int* krowQ=diag0Mat_->krowM(); int* jcolQ=diag0Mat_->jcolM(); double* dQ=diag0Mat_->M();
  for(int i=0; i<SCdim; i++) {
    pend = krowQ[i+1];
    for(p=krowQ[i]; p<pend; p++) {
      j = jcolQ[p]; 
      val = dQ[p];
      dKkt[i][j] += val;
      if(i!=j) dKkt[j][i] += val;
    }
  }

}



int MtxSchurDecompSolver::_numericalFact_RLSS()
{
  if(1==firstCallFlag)
  	_firstCall_RLSS();


  double* eleFullMat  = Msys->getStorageRef().M;

  //update matrices with new entries
  for(int wrkGoff=0; wrkGoff < diag0MatNNz; wrkGoff++) {
	diag0Mat_ele[wrkGoff]	= eleFullMat[Diag0_Full_eleMap[wrkGoff]];
  }
  for(int wrkGoff=0; wrkGoff < diag1MatNNz; wrkGoff++) {
	diag1Mat_ele[wrkGoff]	= eleFullMat[Diag1_Full_eleMap[wrkGoff]];
  }
  for(int wrkGoff=0; wrkGoff < borderMatNNz; wrkGoff++) {
	borderMat_ele[wrkGoff]	= eleFullMat[Border_Full_eleMap[wrkGoff]];
  }  

  int negEVal = 0;

  // set kktSC to 0;
  initializeKKT_Dense(dkktSC);
  
  // First tell children to factorize. 
  negEVal += solver_diag1->matrixChanged();
 
  // build schur complement 
  addTermToDenseSchurCompl(diag1Mat, borderMat, *dkktSC);

  finalizeKKT(diag0Mat, dkktSC);

  negEVal += solver_schur->matrixChanged();

  return negEVal;

};


void MtxSchurDecompSolver::solve( OoqpVector& rhs_ )
{
  SimpleVector& rhs = dynamic_cast<SimpleVector&>(rhs_);

  double *diag0Rhs_ele 		= (double*)malloc(diag0MatDim*sizeof(double));
  double *diag1Rhs_ele 		= (double*)malloc(diag1MatDim*sizeof(double));

  int VarIDinFull;

  for(int VarIDinDiag0=0;VarIDinDiag0<diag0MatDim;VarIDinDiag0++){
	VarIDinFull = schurVarIDinFull[VarIDinDiag0];
	diag0Rhs_ele[VarIDinDiag0] = rhs.elements()[VarIDinFull];
  }
  for(int VarIDinDiag1=0;VarIDinDiag1<diag1MatDim;VarIDinDiag1++){
	VarIDinFull = diag1VarIDinFull[VarIDinDiag1];
	diag1Rhs_ele[VarIDinDiag1] = rhs.elements()[VarIDinFull];
  }

  

  SimpleVector  *diag0_Rhs =  new SimpleVector( diag0Rhs_ele, diag0MatDim );
  SimpleVector  *diag1_Rhs =  new SimpleVector( diag1Rhs_ele, diag1MatDim );


  Lsolve (diag0_Rhs,diag1_Rhs); 
  Dsolve (diag0_Rhs,diag1_Rhs);
  Ltsolve(diag0_Rhs,diag1_Rhs);


  // reorder the vector in orginal order
  for(int VarIDinDiag0=0;VarIDinDiag0<diag0MatDim;VarIDinDiag0++){
	VarIDinFull = schurVarIDinFull[VarIDinDiag0];
	rhs.elements()[VarIDinFull] = diag0Rhs_ele[VarIDinDiag0];
  }
  for(int VarIDinDiag1=0;VarIDinDiag1<diag1MatDim;VarIDinDiag1++){
	VarIDinFull = diag1VarIDinFull[VarIDinDiag1];
	rhs.elements()[VarIDinFull] = diag1Rhs_ele[VarIDinDiag1];
  }
  
  delete (diag1_Rhs);
  delete (diag0_Rhs);
  free (diag1Rhs_ele);
  free (diag0Rhs_ele);



}


void MtxSchurDecompSolver::Lsolve( OoqpVector* diag0_Rhs, OoqpVector* diag1_Rhs )
{
  SimpleVector& b0 = dynamic_cast<SimpleVector&>(*diag0_Rhs);   
  SimpleVector& b1 = dynamic_cast<SimpleVector&>(*diag1_Rhs);

  // Do Lsolve to diag1 matrix
  solver_diag1->Lsolve(b1);

  // build new rhs for schur
  addLnizi(borderMat,b0, b1);

  assert(gSolveSchurScheme==0);
  if(gSolveSchurScheme==0) {
    // Do Lsolve to Schur matrix	
	  solver_schur->Lsolve(b0);    
  }
}



void MtxSchurDecompSolver::Dsolve( OoqpVector* diag0_Rhs, OoqpVector* diag1_Rhs )
{

  SimpleVector& b0 = dynamic_cast<SimpleVector&>(*diag0_Rhs);   
  SimpleVector& b1 = dynamic_cast<SimpleVector&>(*diag1_Rhs);

  //! commented - already done in addLnizi
  //solver_diag1->Dsolve (zi);
  
  assert(gSolveSchurScheme==0);
  if(gSolveSchurScheme==0) {
    // Option 1. - solve with the factors
    solver_schur->Dsolve(b0);
  }

}

void MtxSchurDecompSolver::Ltsolve( OoqpVector* diag0_Rhs, OoqpVector* diag1_Rhs )
{
  SimpleVector& b0 = dynamic_cast<SimpleVector&>(*diag0_Rhs);   
  SimpleVector& b1 = dynamic_cast<SimpleVector&>(*diag1_Rhs);

  // Do Lsolve to Schur matrix
  assert(gSolveSchurScheme==0);
  if(gSolveSchurScheme==0) {
    // Option 1. - solve with the factors
	solver_schur->Ltsolve(b0);
  }

  SimpleVector& x0 = b0; //just another name, for clarity
  
  //! commented - already done in addLnizi
  //solver_diag1->Ltsolve (zi);

  // Li^T\bi for each child i. The backsolve needs z0
  this->LniTransMult(borderMat, b1, -1.0, x0);

}





int MtxSchurDecompSolver::matrixChanged()
{
  int negEVal;
  switch (mtxCase){
	case INPUT_RLSS: 
		negEVal = _numericalFact_RLSS();
		break;
	default:
		assert("Not implement!" &&0);
  }
  return negEVal;
}


void MtxSchurDecompSolver::buildDenseSchurMat(DoubleMatrix* MatIn, int *schurVarIDX)
{
  switch (mtxCase){
	case INPUT_RLSS: 
		_buildDenseSchurMat_RLSS(MatIn, schurVarIDX);
		break;
	default:
		assert("Not implement!" &&0);
		break;
  }
};

void MtxSchurDecompSolver::firstCall()
{

#ifdef TIMING
	double tTot=MPI_Wtime();
#endif  

  switch (mtxCase){
	case INPUT_RLSS: 
		_firstCall_RLSS();
		break;
	default:
		assert("Not implement!" &&0);
  }


#ifdef TIMING
	tTot=MPI_Wtime()-tTot;
	MPI_Barrier(MPI_COMM_WORLD);
	probGenTime += tTot;
#endif
  
};

void MtxSchurDecompSolver::solve1stVarOnly(OoqpVector& rhs_, OoqpVector& sol_SC)
{
  SimpleVector& rhs = dynamic_cast<SimpleVector&>(rhs_);

  double *diag0Rhs_ele 		= (double*)malloc(diag0MatDim*sizeof(double));
  double *diag1Rhs_ele 		= (double*)malloc(diag1MatDim*sizeof(double));

  int VarIDinFull;

  for(int VarIDinDiag0=0;VarIDinDiag0<diag0MatDim;VarIDinDiag0++){
	VarIDinFull = schurVarIDinFull[VarIDinDiag0];
	diag0Rhs_ele[VarIDinDiag0] = rhs.elements()[VarIDinFull];
  }
  for(int VarIDinDiag1=0;VarIDinDiag1<diag1MatDim;VarIDinDiag1++){
	VarIDinFull = diag1VarIDinFull[VarIDinDiag1];
	diag1Rhs_ele[VarIDinDiag1] = rhs.elements()[VarIDinFull];
  }

  SimpleVector  *diag0_Rhs =  new SimpleVector( diag0Rhs_ele, diag0MatDim );
  SimpleVector  *diag1_Rhs =  new SimpleVector( diag1Rhs_ele, diag1MatDim );

  SimpleVector& b0 = dynamic_cast<SimpleVector&>(*diag0_Rhs);   
  SimpleVector& b1 = dynamic_cast<SimpleVector&>(*diag1_Rhs);


  Lsolve (diag0_Rhs,diag1_Rhs); 
  solver_schur->Dsolve(b0);
  solver_schur->Ltsolve(b0);

  sol_SC.copyFrom(*diag0_Rhs);

  delete (diag1_Rhs);
  delete (diag0_Rhs);
  free (diag1Rhs_ele);
  free (diag0Rhs_ele);
};

