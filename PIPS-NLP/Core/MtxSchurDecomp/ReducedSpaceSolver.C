/* PIPS-NLP                                                         	*
 * Authors: Nai-Yuan Chiang                      		*
 * (C) 2015 Argonne National Laboratory			*/

#include "ReducedSpaceSolver.h"

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
#include "UmfPackSolver.h"




#include "OoqpVector.h"
#include "SimpleVector.h"
#include "SimpleVectorHandle.h"

#include "sLinsysRoot.h"
#include "MtxSchurDecompSolver.h"

#include "SaddlePointSolver.h"


extern int gSymLinearSolver;
extern int gOuterSolve;

extern int gRS_SchurSolver;

#ifndef MIN
#define MIN(a,b) ((a > b) ? b : a)
#endif
#ifndef MAX
#define MAX(a,b) ((a > b) ? a : b)
#endif


ReducedSpaceSolver::ReducedSpaceSolver()
	: Ax_solver(NULL), SC_solver(NULL), SparseAug_solver(NULL),
	  firstCallFlag(1), firstSCsolve(1),fullMatDim(0),
	  Hxx_NNz(0), Huu_NNz(0), Hss_NNz(0), Hux_NNz(0), 
	  Ax_NNz(0), Au_NNz(0), Tx_NNz(0), Tu_NNz(0), I_NNz(0), fullMatNNz(0),	
	  Hxx_ele(NULL), Huu_ele(NULL),Hss_ele(NULL),Hux_ele(NULL),
	  Ax_ele(NULL),Au_ele(NULL),Tx_ele(NULL),Tu_ele(NULL),I_ele(NULL),
	  Hxx_rowBeg(NULL),Huu_rowBeg(NULL),Hss_rowBeg(NULL),Hux_rowBeg(NULL),
	  Ax_rowBeg(NULL),Au_rowBeg(NULL), Tx_rowBeg(NULL),Tu_rowBeg(NULL),I_rowBeg(NULL),	
	  Hxx_colIdx(NULL),Huu_colIdx(NULL),Hss_colIdx(NULL),Hux_colIdx(NULL),
	  Ax_colIdx(NULL),Au_colIdx(NULL), Tx_colIdx(NULL),Tu_colIdx(NULL),I_colIdx(NULL),	
	  Hxx_Full_eleMap(NULL),Huu_Full_eleMap(NULL),Hss_Full_eleMap(NULL),Hux_Full_eleMap(NULL),
	  Ax_Full_eleMap(NULL),Au_Full_eleMap(NULL), Tx_Full_eleMap(NULL),Tu_Full_eleMap(NULL),I_Full_eleMap(NULL),
	  Hxx_Mat(NULL),Huu_Mat(NULL),Hss_Mat(NULL),Hux_Mat(NULL), Hxu_Mat(NULL),
	  Ax_Mat(NULL),Au_Mat(NULL), 
	  Tx_Mat(NULL),Tu_Mat(NULL),I_Mat(NULL),
	  Msys(NULL), localNegaEigVal(0),
	  stateVarIDinFull(NULL),
	  FullVarIDinLocal(NULL),	  
	  Hxx_Dim(0), Huu_Dim(0), Hss_Dim(0), Ax_Dim_m(0), Tx_Dim_m(0),	  
	  DecisionVarDim(0), slackVarDim(0), stateVarDim(0), dualYDim(0),	  
	  decisionVarIDinFull(NULL)
{}

ReducedSpaceSolver::ReducedSpaceSolver(DoubleMatrix* MatIn, const int DecisionVarDim_,
					int* DecisionVarIDX_, 
					const int fullVarXSize, const int fullVarYSize, const int fullVarSSize)
	: Ax_solver(NULL), SC_solver(NULL), SparseAug_solver(NULL),
	  firstCallFlag(1),firstSCsolve(1),fullMatDim(0),
	  Hxx_NNz(0), Huu_NNz(0), Hss_NNz(0), Hux_NNz(0), 
	  Ax_NNz(0), Au_NNz(0), Tx_NNz(0), Tu_NNz(0), I_NNz(0), fullMatNNz(0),
	  Hxx_ele(NULL), Huu_ele(NULL),Hss_ele(NULL),Hux_ele(NULL),
	  Ax_ele(NULL),Au_ele(NULL),Tx_ele(NULL),Tu_ele(NULL),I_ele(NULL),
	  Hxx_rowBeg(NULL),Huu_rowBeg(NULL),Hss_rowBeg(NULL),Hux_rowBeg(NULL),
	  Ax_rowBeg(NULL),Au_rowBeg(NULL), Tx_rowBeg(NULL),Tu_rowBeg(NULL),I_rowBeg(NULL),	
	  Hxx_colIdx(NULL),Huu_colIdx(NULL),Hss_colIdx(NULL),Hux_colIdx(NULL),
	  Ax_colIdx(NULL),Au_colIdx(NULL), Tx_colIdx(NULL),Tu_colIdx(NULL),I_colIdx(NULL),	
	  Hxx_Full_eleMap(NULL),Huu_Full_eleMap(NULL),Hss_Full_eleMap(NULL),Hux_Full_eleMap(NULL),
	  Ax_Full_eleMap(NULL),Au_Full_eleMap(NULL), Tx_Full_eleMap(NULL),Tu_Full_eleMap(NULL),I_Full_eleMap(NULL),
	  Hxx_Mat(NULL),Huu_Mat(NULL),Hss_Mat(NULL),Hux_Mat(NULL), Hxu_Mat(NULL),
	  Ax_Mat(NULL),Au_Mat(NULL),
	  Tx_Mat(NULL),Tu_Mat(NULL),I_Mat(NULL),
	  Msys(NULL), localNegaEigVal(0),
	  stateVarIDinFull(NULL),
	  FullVarIDinLocal(NULL),
	  Hxx_Dim(fullVarXSize-DecisionVarDim_), Huu_Dim(DecisionVarDim_), Hss_Dim(fullVarSSize), 
	  Ax_Dim_m(fullVarYSize), Tx_Dim_m(fullVarSSize),
	  DecisionVarDim(DecisionVarDim_), slackVarDim(fullVarSSize), 
	  stateVarDim(fullVarXSize-DecisionVarDim_), dualYDim(fullVarYSize),
	  decisionVarIDinFull(DecisionVarIDX_)
{
  int tempfullDim=0;
  MatIn->getSize(tempfullDim,tempfullDim);
  assert(gOuterSolve >= 3);
  assert(tempfullDim == fullVarXSize+fullVarYSize+fullVarSSize+fullVarSSize);

  fullMatDim = tempfullDim;
  fullMatNNz = ((SparseSymMatrix*)MatIn)->numberOfNonZeros();
  Msys = dynamic_cast<SparseSymMatrix*>(MatIn);

  firstHxxUpdate=true;
  firstHssUpdate=true;
  firstAxUpdate=true;
  firstTxUpdate=true;

  doBuildSc = 1;

}




int _findIDXinState(const int idinFullMat, const int vlength, int *_vector)
{
  int findIDX = idinFullMat;
  for(int i=0;i<vlength;i++)
  	if(_vector[i]<idinFullMat)
	  findIDX--;
  return findIDX; 	  
}

void ReducedSpaceSolver::firstCall()
{
  firstCallFlag = 0;

  //shortcut
  int nx = DecisionVarDim + stateVarDim;
  int ns = slackVarDim;
  int my = dualYDim;
  int mz = ns;

  int findNnz_u=0, findNnz_x=0, findNnz_xu=0;

  int* krowbegFullMat 	= Msys->getStorageRef().krowM;
  int* jcolFullMat 		= Msys->getStorageRef().jcolM;
  double* eleFullMat 	= Msys->getStorageRef().M;


  stateVarIDinFull = (int*) malloc(stateVarDim*sizeof(int));
  FullVarIDinLocal = (int*) malloc(fullMatDim*sizeof(int));
  for(int currVar=0; currVar<fullMatDim;currVar++){
  	FullVarIDinLocal[currVar]=0;
  }

  // build var local-full map and build var  full-local map
  for(int currVar=0; currVar<nx;currVar++){
	int VarID_local = _findIDX(currVar,DecisionVarDim,decisionVarIDinFull);
	if(VarID_local == -1){
	  //this var is state var
	  VarID_local = _findIDXinState(currVar,DecisionVarDim,decisionVarIDinFull);
	  stateVarIDinFull[VarID_local] = currVar;
	  FullVarIDinLocal[currVar] =  (VarID_local+1);
	}else{
	  //this var is decision
	  FullVarIDinLocal[currVar] = -(VarID_local+1);
	}
  }

  //count nozeros in Hxx, Huu, Hux
  Hxx_rowBeg  = (int*) malloc((stateVarDim+1)*sizeof(int));
  Huu_rowBeg  = (int*) malloc((DecisionVarDim+1)*sizeof(int));
  Hux_rowBeg  = (int*) malloc((DecisionVarDim+1)*sizeof(int));
  for(int jj=0;jj<stateVarDim+1;jj++) Hxx_rowBeg[jj]=0;
  for(int jj=0;jj<DecisionVarDim+1;jj++){Huu_rowBeg[jj]=0; Hux_rowBeg[jj]=0;}
  
  for(int currRow=0; currRow<nx;currRow++){
	for( int k=krowbegFullMat[currRow];k<krowbegFullMat[currRow+1];k++){	
	  int currCol = jcolFullMat[k];

	  int localColID = FullVarIDinLocal[currCol];
	  int localRowID = FullVarIDinLocal[currRow];
	  
	  if(localColID < 0 && localRowID < 0){
		// in Huu! correct the neg index
		localColID = -(localColID+1);
		localRowID = -(localRowID+1);
		Huu_NNz++;
		if(localRowID >= localColID)
		  Huu_rowBeg[localRowID+1]++;
		else
		  Huu_rowBeg[localColID+1]++;
	  }else if(localColID > 0 && localRowID >0){
		// in Hxx, correct the pos index
		localColID = (localColID-1);
		localRowID = (localRowID-1);
		Hxx_NNz++;
		if(localRowID >= localColID)
		  Hxx_rowBeg[localRowID+1]++;
		else
		  Hxx_rowBeg[localColID+1]++;	
	  }else if(localColID < 0 && localRowID > 0 ){
		//belong to Hux
		localColID = -(localColID+1);
		Hux_NNz++;
		Hux_rowBeg[localColID+1]++;	  
	  }else if(localColID > 0 && localRowID < 0){
		//belong to borderMat
		localRowID = -(localRowID+1);
		Hux_NNz++;
		Hux_rowBeg[localRowID+1]++;	  
	  }else{assert("impossible"&&0);}
	}
  }
  for(int j=1; j < stateVarDim+1; j++)Hxx_rowBeg[j] += Hxx_rowBeg[j-1]; 
  for(int j=1; j < DecisionVarDim+1; j++){
  	Huu_rowBeg[j] += Huu_rowBeg[j-1]; 
	Hux_rowBeg[j] += Hux_rowBeg[j-1];
  }

  //copy Huu, Hxx, Hux
  // alocate space
  Huu_colIdx 	= (int*)malloc(Huu_NNz*sizeof(int));
  Huu_ele 	= (double*)malloc(Huu_NNz*sizeof(double)); 
  Huu_Full_eleMap 	= (int*)malloc(Huu_NNz*sizeof(int));
  Hxx_colIdx 	= (int*)malloc(Hxx_NNz*sizeof(int));
  Hxx_ele 	= (double*)malloc(Hxx_NNz*sizeof(double)); 
  Hxx_Full_eleMap 	= (int*)malloc(Hxx_NNz*sizeof(int));
  Hux_colIdx 	= (int*)malloc(Hux_NNz*sizeof(int));
  Hux_ele 	= (double*)malloc(Hux_NNz*sizeof(double)); 
  Hux_Full_eleMap 	= (int*)malloc(Hux_NNz*sizeof(int));

  int *nextHxxInRow  = (int*)malloc(stateVarDim*sizeof(int));
  int *nextHuuInRow  = (int*)malloc(DecisionVarDim*sizeof(int));
  int *nextHuxInRow  = (int*)malloc(DecisionVarDim*sizeof(int));

  for(int jj=0;jj<stateVarDim;jj++) nextHxxInRow[jj] = Hxx_rowBeg[jj];
  for(int jj=0;jj<DecisionVarDim;jj++){
  	nextHuuInRow[jj] = Huu_rowBeg[jj];
    nextHuxInRow[jj] = Hux_rowBeg[jj];
  }
  
  findNnz_u=0;  findNnz_x=0;  findNnz_xu==0;
  int wrkGoff=0;
  for(int currRow=0; currRow<nx;currRow++){
	for( int k=krowbegFullMat[currRow];k<krowbegFullMat[currRow+1];k++){	
	  int currCol = jcolFullMat[k];

	  int localColID = FullVarIDinLocal[currCol];
	  int localRowID = FullVarIDinLocal[currRow];
	  
	  if(localColID < 0 && localRowID < 0){
		// in Huu! correct the neg index
		localColID = -(localColID+1);
		localRowID = -(localRowID+1);
		if(localRowID >= localColID){
		  wrkGoff = nextHuuInRow[localRowID]++;
		  Huu_colIdx[wrkGoff]  = localColID;
		}
		else{
		  wrkGoff = nextHuuInRow[localColID]++;
		  Huu_colIdx[wrkGoff]  = localRowID;
		}
		Huu_ele[wrkGoff] = eleFullMat[k];
	    Huu_Full_eleMap[wrkGoff] = k;
		findNnz_u++;
	  }else if(localColID > 0 && localRowID >0){
		// in Hxx, correct the pos index
		localColID = (localColID-1);
		localRowID = (localRowID-1);
		if(localRowID >= localColID){
		  wrkGoff = nextHxxInRow[localRowID]++;
		  Hxx_colIdx[wrkGoff]  = localColID;
		}
		else{
		  wrkGoff = nextHxxInRow[localColID]++;
		  Hxx_colIdx[wrkGoff]  = localRowID;
		}
		Hxx_ele[wrkGoff] = eleFullMat[k];
	    Hxx_Full_eleMap[wrkGoff] = k;
		findNnz_x++;
	  }else if(localColID < 0 && localRowID > 0 ){
		//belong to Hux
		localColID = -(localColID+1);		
		localRowID = (localRowID-1);
		wrkGoff = nextHuxInRow[localColID]++;
		Hux_colIdx[wrkGoff] 	= localRowID;
		Hux_ele[wrkGoff]		= eleFullMat[k];
		Hux_Full_eleMap[wrkGoff] = k;
		findNnz_xu++;	
	  }else if(localColID > 0 && localRowID < 0){
		//belong to borderMat
		localColID = (localColID-1);
		localRowID = -(localRowID+1);
		wrkGoff = nextHuxInRow[localRowID]++;
		Hux_colIdx[wrkGoff] 	= localColID;
		Hux_ele[wrkGoff]		= eleFullMat[k];
		Hux_Full_eleMap[wrkGoff] = k;
		findNnz_xu++;  
	  }else{assert("impossible"&&0);}
	}
  }
  free(nextHuuInRow);
  free(nextHxxInRow);
  free(nextHuxInRow);


  //copy Hss
  Hss_rowBeg  = (int*) malloc((slackVarDim+1)*sizeof(int));
  for(int jj=0;jj<slackVarDim+1;jj++) Hss_rowBeg[jj]=0;

  Hss_NNz = krowbegFullMat[nx+ns]-krowbegFullMat[nx];
  Hss_colIdx  		= (int*)malloc(Hss_NNz*sizeof(int));
  Hss_ele 			= (double*)malloc(Hss_NNz*sizeof(double));
  Hss_Full_eleMap  	= (int*)malloc(Hss_NNz*sizeof(int));  
  findNnz_u=0;
  for(int currRow=nx; currRow<nx+ns;currRow++){
	for( int k=krowbegFullMat[currRow];k<krowbegFullMat[currRow+1];k++){	
	  int currCol = jcolFullMat[k];
	  assert(currCol>=nx && currCol<=currRow);
	  Hss_rowBeg[currRow-nx+1]++;
	  Hss_colIdx[findNnz_u] = currCol-nx;
	  Hss_ele[findNnz_u] = eleFullMat[k];
	  Hss_Full_eleMap[findNnz_u] = k;
	  findNnz_u++;
	}
  }
  for(int j=1; j < ns+1; j++){Hss_rowBeg[j] += Hss_rowBeg[j-1];}
  assert(findNnz_u == Hss_NNz);


  //count nozeros in Ax, Au
  Au_rowBeg  = (int*) malloc((my+1)*sizeof(int));  
  Ax_rowBeg  = (int*) malloc((my+1)*sizeof(int));
  for(int jj=0;jj<my+1;jj++){ Au_rowBeg[jj]=0;  Ax_rowBeg[jj]=0; }
  
  for(int currRow=nx+ns; currRow<nx+ns+my;currRow++){
  	// we use krowbegFullMat[currRow+1]-1 since the last ele is given bt dual regularization
	for( int k=krowbegFullMat[currRow];k<krowbegFullMat[currRow+1]-1;k++){	
	  int currCol = jcolFullMat[k];
	  
	  int localColID = FullVarIDinLocal[currCol];
	  int localRowID = currRow-nx-ns;
	  
	  if(localColID < 0 ){
		// in Au
		Au_NNz++;
		Au_rowBeg[localRowID+1]++;
	  }else if(localColID > 0 ){
		// in Ax
		Ax_NNz++;
		Ax_rowBeg[localRowID+1]++;
	  }else{ assert("impossible"&&0);}
	}
  }
  for(int j=1; j < my+1; j++){Au_rowBeg[j] += Au_rowBeg[j-1];Ax_rowBeg[j] += Ax_rowBeg[j-1];}

  // copy Ax, Au
  Ax_colIdx  		= (int*)malloc(Ax_NNz*sizeof(int));
  Ax_ele 			= (double*)malloc(Ax_NNz*sizeof(double));
  Ax_Full_eleMap  	= (int*)malloc(Ax_NNz*sizeof(int));   
  Au_colIdx  		= (int*)malloc(Au_NNz*sizeof(int));
  Au_ele 			= (double*)malloc(Au_NNz*sizeof(double));
  Au_Full_eleMap  	= (int*)malloc(Au_NNz*sizeof(int)); 
  findNnz_u=0; findNnz_x=0;  
  for(int currRow=nx+ns; currRow<nx+ns+my;currRow++){
	for( int k=krowbegFullMat[currRow];k<krowbegFullMat[currRow+1]-1;k++){	
	  int currCol = jcolFullMat[k];
	  assert(currRow>=currCol);

	  int localColID = FullVarIDinLocal[currCol];
	  int localRowID = currRow-nx-ns;
	  
	  if(localColID < 0 ){
		// in Au! correct the neg index
		localColID = -(localColID+1);
		Au_colIdx[findNnz_u] = localColID;
		Au_ele[findNnz_u] = eleFullMat[k];
		Au_Full_eleMap[findNnz_u] = k;
		findNnz_u++;
	  }else if(localColID > 0 ){
		// in Ax, correct the pos index
		localColID = (localColID-1);
		Ax_colIdx[findNnz_x] = localColID;
		Ax_ele[findNnz_x] = eleFullMat[k];
		Ax_Full_eleMap[findNnz_x] = k;
		findNnz_x++;
	  }else{assert("impossible"&&0);}
	}
  }

  //count nozeros in Tx, Tu
  Tu_rowBeg  = (int*) malloc((mz+1)*sizeof(int));  
  Tx_rowBeg  = (int*) malloc((mz+1)*sizeof(int));
  for(int jj=0;jj<mz+1;jj++){ Tu_rowBeg[jj]=0;  Tx_rowBeg[jj]=0; }
  
  for(int currRow=nx+ns+my; currRow<nx+ns+my+mz;currRow++){
	for( int k=krowbegFullMat[currRow];k<krowbegFullMat[currRow+1]-1;k++){	
	  int currCol = jcolFullMat[k];
	  
	  int localColID = FullVarIDinLocal[currCol];
	  int localRowID = currRow-nx-ns-my;
	  
	  if(localColID < 0 ){
		// in Tu
		Tu_NNz++;
		Tu_rowBeg[localRowID+1]++;
	  }else if(localColID > 0 ){
		// in Tx
		Tx_NNz++;
		Tx_rowBeg[localRowID+1]++;
	  }else{
	    // in -I
	    assert(eleFullMat[k]==-1);
	  }
	}
  }
  for(int j=1; j < mz+1; j++){Tu_rowBeg[j] += Tu_rowBeg[j-1]; Tx_rowBeg[j] += Tx_rowBeg[j-1];} 

  // copy Tx, Tu
  Tx_colIdx  		= (int*)malloc(Tx_NNz*sizeof(int));
  Tx_ele 			= (double*)malloc(Tx_NNz*sizeof(double));
  Tx_Full_eleMap  	= (int*)malloc(Tx_NNz*sizeof(int));   
  Tu_colIdx  		= (int*)malloc(Tu_NNz*sizeof(int));
  Tu_ele 			= (double*)malloc(Tu_NNz*sizeof(double));
  Tu_Full_eleMap  	= (int*)malloc(Tu_NNz*sizeof(int)); 
  findNnz_u=0; findNnz_x=0;  
  for(int currRow=nx+ns+my; currRow<nx+ns+my+mz;currRow++){
	for( int k=krowbegFullMat[currRow];k<krowbegFullMat[currRow+1]-1;k++){	
	  int currCol = jcolFullMat[k];
	  
	  int localColID = FullVarIDinLocal[currCol];
	  int localRowID = currRow-nx-ns-my;
	  
	  if(localColID < 0 ){
		// in Tu! correct the neg index
		localColID = -(localColID+1);
		Tu_colIdx[findNnz_u] = localColID;
		Tu_ele[findNnz_u] = eleFullMat[k];
		Tu_Full_eleMap[findNnz_u] = k;
		findNnz_u++;
	  }else if(localColID > 0 ){
		// in Tx, correct the pos index
		localColID = (localColID-1);
		Tx_colIdx[findNnz_x] = localColID;
		Tx_ele[findNnz_x] = eleFullMat[k];
		Tx_Full_eleMap[findNnz_x] = k;
		findNnz_x++;
	  }
	}
  }

  // dualYDim+ns is given by dual reg
  assert(Hxx_NNz+Huu_NNz+Hss_NNz+Hux_NNz+Tx_NNz+Tu_NNz+
  		 Ax_NNz+Au_NNz+ns+dualYDim+ns==fullMatNNz);


  Hxx_Mat = new SparseSymMatrix( stateVarDim, Hxx_NNz, 
				  Hxx_rowBeg, Hxx_colIdx, Hxx_ele);
  Huu_Mat = new SparseSymMatrix( DecisionVarDim, Huu_NNz, 
				  Huu_rowBeg, Huu_colIdx, Huu_ele);
  Hux_Mat = new SparseGenMatrix( DecisionVarDim, stateVarDim, Hux_NNz, 
				  Hux_rowBeg, Hux_colIdx, Hux_ele);

  Hss_Mat = new SparseSymMatrix( ns, Hss_NNz, 
				  Hss_rowBeg, Hss_colIdx, Hss_ele);

  Au_Mat = new SparseGenMatrix( my, DecisionVarDim, Au_NNz, 
				  Au_rowBeg, Au_colIdx, Au_ele);
  Ax_Mat = new SparseGenMatrix( my, stateVarDim, Ax_NNz, 
				  Ax_rowBeg, Ax_colIdx, Ax_ele);

  Tu_Mat = new SparseGenMatrix( mz, DecisionVarDim, Tu_NNz, 
				  Tu_rowBeg, Tu_colIdx, Tu_ele);
  Tx_Mat = new SparseGenMatrix( mz, stateVarDim, Tx_NNz, 
				  Tx_rowBeg, Tx_colIdx, Tx_ele);

  Ax_solver = dynamic_cast<UmfPackSolver*> (new UmfPackSolver(Ax_Mat));

}


void ReducedSpaceSolver::solve(GenMatrix& rhs_in)
{
  DenseGenMatrix &rhs = dynamic_cast<DenseGenMatrix&>(rhs_in);
  int N,NRHS;
  // rhs vectors are on the "rows", for continuous memory
  rhs.getSize(NRHS,N);
  assert(fullMatDim==N);
    
  for (int i = 0; i < NRHS; i++) {
    SimpleVector v(rhs[i],N);
    solve(v);
  } 

}


void ReducedSpaceSolver::solve ( OoqpVector& rhs_ )
{
  int nx = stateVarDim+DecisionVarDim, my=dualYDim, ns=slackVarDim, mz=ns;
  
  SimpleVector& rhs = dynamic_cast<SimpleVector&>(rhs_);
  double *x_Rhs_ele 		= (double*)malloc(stateVarDim*sizeof(double));
  double *u_Rhs_ele 		= (double*)malloc(DecisionVarDim*sizeof(double));
  double *s_Rhs_ele, *y_Rhs_ele, *z_Rhs_ele;

  int VarIDinFull;

  for(int VarIDLocal=0;VarIDLocal<stateVarDim;VarIDLocal++){
	VarIDinFull = stateVarIDinFull[VarIDLocal];
	x_Rhs_ele[VarIDLocal] = rhs.elements()[VarIDinFull];
  }
  for(int VarIDLocal=0;VarIDLocal<DecisionVarDim;VarIDLocal++){
	VarIDinFull = decisionVarIDinFull[VarIDLocal];
	u_Rhs_ele[VarIDLocal] = rhs.elements()[VarIDinFull];
  }
  s_Rhs_ele = &(rhs.elements()[nx]);
  y_Rhs_ele = &(rhs.elements()[nx+ns]);
  z_Rhs_ele = &(rhs.elements()[nx+ns+my]);
  
  SimpleVector  *x_Rhs =  new SimpleVector( x_Rhs_ele, stateVarDim );
  SimpleVector  *u_Rhs =  new SimpleVector( u_Rhs_ele, DecisionVarDim );
  SimpleVector  *s_Rhs =  new SimpleVector( s_Rhs_ele, ns );
  SimpleVector  *y_Rhs =  new SimpleVector( y_Rhs_ele, my ); 
  SimpleVector  *z_Rhs =  new SimpleVector( z_Rhs_ele, mz );
  
  assert(my==stateVarDim);
  
  // build temp vector
  int N = stateVarDim+ns+my+mz;
  SimpleVector TempVec(N);
  SimpleVector x_temp(&TempVec[0], stateVarDim); 		x_temp.copyFrom(*x_Rhs);
  SimpleVector s_temp(&TempVec[stateVarDim], ns);		s_temp.copyFrom(*s_Rhs);
  SimpleVector y_temp(&TempVec[stateVarDim+ns], my);	y_temp.copyFrom(*y_Rhs);
  SimpleVector z_temp(&TempVec[stateVarDim+ns+my], mz);	z_temp.copyFrom(*z_Rhs);
  
  // build new u_Rhs as step (a) r_u -Ac'*C^{-T}(r_p-Hc*C^{-1}r_d)-Huc*C^{-1}r_d
  reducedSpaceJacSolve(&y_temp,&z_temp);			// C^{-1}r_d								rhs_p 
  Hxx_Mat->mult(1,x_temp,-1,y_temp); 	  			// r_p-Hc*C^{-1}r_d
  if(ns!=0) Hss_Mat->mult(1,s_temp,-1,z_temp); 	  	// r_p-Hc*C^{-1}r_d
  reducedSpaceJacTransSolve(&x_temp,&s_temp);		// C^{-T}(r_p-Hc*C^{-1}r_d)					rhs_d
  
  Hux_Mat->mult(1,*u_Rhs,-1,y_temp); 	  			// r_u = r_u - Huc*C^{-1}r_d	
  Au_Mat->transMult(1,*u_Rhs,-1,x_temp);			// r_u = r_u - Ac'*C^{-T}(r_p-Hc*C^{-1}r_d)		for x
  if(ns!=0) Tu_Mat->transMult(1,*u_Rhs,-1,s_temp);  // r_u = r_u - Ac'*C^{-T}(r_p-Hc*C^{-1}r_d)		for s

  // solve U system;
  SimpleVector *rhs_ptr = &rhs;
  solveDeltaU(*u_Rhs,rhs_ptr);

  // short cut
  SimpleVector *u_sol = u_Rhs;  

  //compute \delta_p = - C^{-1}Ac*\delta_u
  Au_Mat->mult(0,*x_Rhs,-1.0,*u_sol);					// -Ac*\delta_u						for x
  if(ns!=0) Tu_Mat->mult(0.0,*s_Rhs,-1.0,*u_sol);		// -Ac*\delta_u						for s
  reducedSpaceJacSolve(x_Rhs,s_Rhs);		  			// -C^{-1}Ac*\delta_u

  // compute \delta_d =  - C^{-T}	(Hcu*\delta_u -Hc*C^{-1}*Ac*\delta_u )
  Hxx_Mat->mult(0,*y_Rhs,-1,*x_Rhs); 	  				// Hc * C^{-1}*Ac*\delta_u			for y
  if(ns!=0)Hss_Mat->mult(0,*z_Rhs,-1,*s_Rhs); 	  		// Hc * C^{-1}*Ac*\delta_u			for z
  Hux_Mat->transMult(1,*y_Rhs,-1,*u_sol); 	  			// -Hcu*\delta_u + Hc*C^{-1}*Ac*\delta_u
  reducedSpaceJacTransSolve(y_Rhs,z_Rhs);		   		// - C^{-T}  (Hcu*\delta_u -Hc*C^{-1}*Ac*\delta_u )

  //compute \delta_p = C^{-1}(r_d-Ac*\delta_u) = rhs_p - C^{-1}Ac*\delta_u
  x_Rhs->axpy(1.0,y_temp);								// rhs_p - C^{-1}Ac*\delta_u		for x
  if(ns!=0) s_Rhs->axpy(1.0,z_temp);					// rhs_p - C^{-1}Ac*\delta_u		for s
  
  // compute \delta_d = rhs_d - C^{-T}  (Hcu*\delta_u -Hc*C^{-1}*Ac*\delta_u )  
  y_Rhs->axpy(1.0,x_temp);								// rhs_d - C^{-T}  (Hcu*\delta_u -Hc*C^{-1}*Ac*\delta_u )		for y
  if(ns!=0) z_Rhs->axpy(1.0,s_temp);					// rhs_d - C^{-T}  (Hcu*\delta_u -Hc*C^{-1}*Ac*\delta_u )		for s

	
  // reorder the vector in orginal order
  for(int VarIDLocal=0;VarIDLocal<stateVarDim;VarIDLocal++){
	VarIDinFull = stateVarIDinFull[VarIDLocal];
	rhs.elements()[VarIDinFull] = x_Rhs_ele[VarIDLocal];
  }
  for(int VarIDLocal=0;VarIDLocal<DecisionVarDim;VarIDLocal++){
	VarIDinFull = decisionVarIDinFull[VarIDLocal];
	rhs.elements()[VarIDinFull] = u_Rhs_ele[VarIDLocal];
  }


  delete (z_Rhs);
  delete (y_Rhs);
  delete (s_Rhs);
  delete (u_Rhs);
  delete (x_Rhs);
  free (u_Rhs_ele);
  free (x_Rhs_ele);  
}

int ReducedSpaceSolver::_numericalFact ( )
{
  if(1==firstCallFlag)
	firstCall();

  double* eleFullMat  = Msys->getStorageRef().M;

  //update matrices with new entries
  for(int wrkGoff=0; wrkGoff < Huu_NNz; wrkGoff++)
	Huu_ele[wrkGoff] = eleFullMat[Huu_Full_eleMap[wrkGoff]];
  for(int wrkGoff=0; wrkGoff < Hxx_NNz; wrkGoff++)
	Hxx_ele[wrkGoff] = eleFullMat[Hxx_Full_eleMap[wrkGoff]];
  for(int wrkGoff=0; wrkGoff < Hux_NNz; wrkGoff++)
	Hux_ele[wrkGoff] = eleFullMat[Hux_Full_eleMap[wrkGoff]];

  for(int wrkGoff=0; wrkGoff < Hss_NNz; wrkGoff++)
	Hss_ele[wrkGoff] = eleFullMat[Hss_Full_eleMap[wrkGoff]];

  
  for(int wrkGoff=0; wrkGoff < Ax_NNz; wrkGoff++)
	Ax_ele[wrkGoff] = eleFullMat[Ax_Full_eleMap[wrkGoff]];
  for(int wrkGoff=0; wrkGoff < Au_NNz; wrkGoff++)
	Au_ele[wrkGoff] = eleFullMat[Au_Full_eleMap[wrkGoff]];  
  for(int wrkGoff=0; wrkGoff < Tx_NNz; wrkGoff++)
	Tx_ele[wrkGoff] = eleFullMat[Tx_Full_eleMap[wrkGoff]];
  for(int wrkGoff=0; wrkGoff < Tu_NNz; wrkGoff++)
	Tu_ele[wrkGoff] = eleFullMat[Tu_Full_eleMap[wrkGoff]];  

  int negEVal = 0;

  doBuildSc =1;


  Ax_solver->matrixChanged();

  negEVal = dualYDim + slackVarDim;
  return negEVal;

}


// do C^{-1} where C is the reduced Jac	(accoring to (*))
void ReducedSpaceSolver::reducedSpaceJacSolve (OoqpVector* rhs_y_, OoqpVector* rhs_z_)
{
  SimpleVector* rhs_y = dynamic_cast<SimpleVector*>(rhs_y_);
  int A_m, A_n;
  Ax_Mat->getSize(A_m,A_n);
  assert(A_m==A_n && A_m==rhs_y->length());

  //do A solve
  Ax_solver->solve(*rhs_y);

  // if we have slack variable
  if(slackVarDim!=0){
  	SimpleVector* rhs_z = dynamic_cast<SimpleVector*>(rhs_z_);
    Tx_Mat->mult(-1,*rhs_z,1,*rhs_y);
  }
  
}

// do C^{-T} where C is the reduced Jac
void ReducedSpaceSolver::reducedSpaceJacTransSolve (OoqpVector* rhs_x_, OoqpVector* rhs_s_)
{

  SimpleVector* rhs_x = dynamic_cast<SimpleVector*>(rhs_x_);
  int A_m, A_n;
  Ax_Mat->getSize(A_m,A_n);
  assert(A_m==A_n && A_m==rhs_x->length());

  // if we have slack variable
  if(slackVarDim!=0){
  	SimpleVector* rhs_s = dynamic_cast<SimpleVector*>(rhs_s_);
	Tx_Mat->transMult(1,*rhs_x,1,*rhs_s);
	rhs_s->negate();
  }

  //do A trans solve
  Ax_solver->solveTrans(*rhs_x);
  
}


void ReducedSpaceSolver::schursolveDeltaU (OoqpVector& rhs_, OoqpVector* rhs_Full)
{
  SimpleVector& rhs = dynamic_cast<SimpleVector&>(rhs_);

  int *matType = new int[4];
  matType[0]   = SCHURSOLVER_INPUT_SYMMETRIC;
  matType[1]   = SCHURSOLVER_INPUT_SPARSE;	  
  matType[2]   = SCHURSOLVER_INPUT_LOWER;	  
  matType[3]   = SCHURSOLVER_INPUT_ROWWISE;	  

  SC_solver = new MtxSchurDecompSolver(Msys,DecisionVarDim,matType,
  					decisionVarIDinFull, dualYDim+slackVarDim);

  SC_solver->matrixChanged();

  ((MtxSchurDecompSolver*)SC_solver)->solve1stVarOnly(*rhs_Full,rhs_);
  

  
}

void ReducedSpaceSolver::_schursolveDeltaU_BuildSC_firstCall ()
{
  firstSCsolve=0;
  
  int n = stateVarDim+slackVarDim+dualYDim+slackVarDim;
  int nnzAug = Hxx_NNz + Hss_NNz + Ax_NNz + Tx_NNz + slackVarDim + n;
  //alocate the matrix and copy the data into
  kkt_diag = new SparseSymMatrix(n,nnzAug);

  SimpleVectorHandle v( new SimpleVector(n) );
  v->setToZero();
  kkt_diag->setToDiagonal(*v);
  
  kkt_diag->symAtPutSubmatrix( 0, 0, *Hxx_Mat, 0, 0, stateVarDim, stateVarDim);
  kkt_diag->symAtPutSubmatrix( stateVarDim+slackVarDim, 0, *Ax_Mat, 0, 0, dualYDim, stateVarDim);


  int len = slackVarDim;	
  int *tempDiagRowId;
  int *tempDiagColId;
  double *tempDiagEleId;
  int info;

  if(len>0){
  	kkt_diag->symAtPutSubmatrix( stateVarDim, stateVarDim, *Hss_Mat, 0, 0, slackVarDim, slackVarDim);
    kkt_diag->symAtPutSubmatrix( stateVarDim+slackVarDim+dualYDim, 0, *Tx_Mat, 0, 0, slackVarDim, stateVarDim);	  

	SparseSymMatrixHandle tempDiagMat( new SparseSymMatrix( slackVarDim, slackVarDim ) );
	tempDiagRowId = new int[slackVarDim];
	tempDiagColId = new int[slackVarDim];
	tempDiagEleId = new double[slackVarDim];
	
    for(int i=0;i<slackVarDim;i++){
	  tempDiagRowId[i]=i;
	  tempDiagColId[i]=i;
	  tempDiagEleId[i]=-1;
    }
    tempDiagMat->putSparseTriple( tempDiagRowId,len, tempDiagColId, tempDiagEleId,info );
    kkt_diag->symAtPutSubmatrix( stateVarDim+slackVarDim+dualYDim, stateVarDim, *tempDiagMat, 0, 0, 
  								slackVarDim, slackVarDim);
  }

  SparseAug_solver = NULL;
  if(0==gSymLinearSolver){
#ifdef WITH_MA27  	
  	SparseAug_solver = new Ma27Solver( kkt_diag);
#endif
  }else if(1==gSymLinearSolver){
#ifdef WITH_MA57  	
   SparseAug_solver = new Ma57Solver( kkt_diag);
#endif
  }else if(2==gSymLinearSolver){
#ifdef WITH_PARDISO
    SparseAug_solver = new PardisoSolver( kkt_diag,dualYDim+slackVarDim);
#endif 
  }else if(3==gSymLinearSolver){
#ifdef WITH_UMFPACK  
   SparseAug_solver = new UmfPackSolver( kkt_diag,dualYDim+slackVarDim);
#endif 
  }
//  else if(4==gSymLinearSolver)
//   SparseAug_solver = new SaddlePointSolver(kkt_diag,stateVarDim,dualYDim,slackVarDim);

  assert(SparseAug_solver);

  kktsc =  new DenseSymMatrix(DecisionVarDim);
  SC_solver = dynamic_cast<DeSymIndefSolver*> (new DeSymIndefSolver(kktsc));

  if(len>0){
    delete[] tempDiagRowId;
    delete[] tempDiagColId;
    delete[] tempDiagEleId;
  }
}

void ReducedSpaceSolver::schursolveDeltaU_BuildSC (OoqpVector& rhs_, OoqpVector* rhs_Full)
{
  int negEVal = 0;

  if(firstSCsolve==1){
	_schursolveDeltaU_BuildSC_firstCall();
  }else{  
	kkt_diag->symAtSetSubmatrix( 0, 0, *Hxx_Mat, 0, 0, stateVarDim, stateVarDim,
								firstHxxUpdate, LocHxxMap);
	kkt_diag->symAtSetSubmatrix( stateVarDim+slackVarDim, 0, *Ax_Mat, 0, 0, dualYDim, stateVarDim, 
								firstAxUpdate, LocAxMap);
	if(slackVarDim>0){
	  kkt_diag->symAtSetSubmatrix( stateVarDim, stateVarDim, *Hss_Mat, 0, 0, slackVarDim, slackVarDim,
	  							firstHssUpdate, LocHssMap);	
	  kkt_diag->symAtSetSubmatrix( stateVarDim+slackVarDim+dualYDim, 0, *Tx_Mat, 0, 0, slackVarDim, stateVarDim,
	  							firstTxUpdate, LocTxMap);
	}
	firstHxxUpdate=false;
	firstHssUpdate=false;
	firstAxUpdate=false;
	firstTxUpdate=false;
  }

  // set kktSC to 0;
  MtxSchurDecompSolver *tempSolver = new MtxSchurDecompSolver();

  if(doBuildSc=1){
    tempSolver->initializeKKT_Dense(kktsc);
    negEVal += SparseAug_solver->matrixChanged();
    addTermToDenseSchurCompl(*kktsc);
	tempSolver->finalizeKKT(Huu_Mat, kktsc);
	doBuildSc=0;
  }

  negEVal += SC_solver->matrixChanged();
  negEVal = dualYDim+slackVarDim;
  assert(negEVal == dualYDim+slackVarDim);
 
  SC_solver->solve(rhs_);

  delete tempSolver;
  
}

void ReducedSpaceSolver::addTermToDenseSchurCompl( DenseSymMatrix& SC) 
{
  if(Hxu_Mat == NULL)
    Hxu_Mat = new SparseGenMatrix( stateVarDim, DecisionVarDim, Hux_NNz);
  int* krowHxu=Hxu_Mat->getStorageRef().krowM;
  int* jcolHxu=Hxu_Mat->getStorageRef().jcolM;
  double *MHxu=Hxu_Mat->getStorageRef().M;

  Hux_Mat->getStorageRef().transpose(krowHxu, jcolHxu, MHxu);

  
  int N, nxP, NP,mR,nR;
  int locns = slackVarDim;
  
  Au_Mat->getSize(N, nxP); assert(N==dualYDim);
  NP = SC.size(); assert(NP>=nxP);

  if(nxP==-1) Tu_Mat->getSize(N,nxP);
  if(nxP==-1) nxP = NP;

  N = stateVarDim+locns+dualYDim+locns;
  
  int blocksize = 64;
  DenseGenMatrix cols(blocksize,N);


  for (int it=0; it < nxP; it += blocksize) {
    int start=it;
    int end = MIN(it+blocksize,nxP);
    int numcols = end-start;
    cols.getStorageRef().m = numcols; // avoid extra solves


    bool allzero = true;
    memset(&cols[0][0],0,N*blocksize*sizeof(double));

   
	Hxu_Mat->getStorageRef().fromGetColBlock(start, &cols[0][0], N, numcols, allzero);
	Au_Mat->getStorageRef().fromGetColBlock(start, &cols[0][stateVarDim+locns], N, numcols, allzero);
	Tu_Mat->getStorageRef().fromGetColBlock(start, &cols[0][stateVarDim+locns+dualYDim], N, numcols, allzero);
	  
    if(!allzero) {
	  SparseAug_solver->solve(cols);
  
      Hxu_Mat->getStorageRef().transMultMat( 1.0, SC[start], numcols, NP,  
       				      -1.0, &cols[0][0], N);
      Au_Mat->getStorageRef().transMultMat( 1.0, SC[start], numcols, NP,  
       				      -1.0, &cols[0][stateVarDim+locns], N);
      Tu_Mat->getStorageRef().transMultMat( 1.0, SC[start], numcols, NP,
       				      -1.0, &cols[0][stateVarDim+locns+dualYDim], N);
    } //end !allzero
  }

}




void ReducedSpaceSolver::solveDeltaU (OoqpVector& rhs_, OoqpVector* rhs_Full)
{
  SimpleVector& rhs = dynamic_cast<SimpleVector&>(rhs_);

  //do A trans solve
  if(1==gRS_SchurSolver)
  	schursolveDeltaU_BuildSC(rhs_,rhs_Full);
  else 
  	schursolveDeltaU(rhs_,rhs_Full);
		
}

