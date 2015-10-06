/* PIPS-NLP                                                         	*
 * Authors: Nai-Yuan Chiang                      		*
 * (C) 2015 Argonne National Laboratory			*/

#include "PartitionAugSolver.h"

#include "assert.h"
#include "stdlib.h"


#include <iostream>
#include <string>
#include <fstream>
#include <limits>
#include <sstream>


#include "DoubleMatrix.h"
#include "SparseSymMatrix.h"

#include "SparseSymMatrix.h"
#include "SparseGenMatrix.h"

#include "RegularizationAlg.h"
#include "NlpGenLinsys.h"

#include "DeSymIndefSolver.h"
#include "sLinsysRoot.h"

#include "OoqpVector.h"
#include "SimpleVector.h"

#ifdef WITH_MA27
#include "Ma27Solver.h"
#endif
#ifdef WITH_MA57
#include "Ma57Solver.h"
#endif
#ifdef WITH_PARDISO
#include "PardisoSolver.h"
#endif

using namespace std;


extern int gBuildSchurComp;
extern int gSymLinearSolver;
extern int gSolveSchurScheme;

#ifndef MIN
#define MIN(a,b) ((a > b) ? b : a)
#endif

#ifdef TIMING
  #include "mpi.h"
  extern double probGenTime;
  extern double PartSolver_GenTime;
  extern double PartSolver_SolTime;
  extern double PartSolver_FactTime;
  extern int call_sol_Times;
  extern int call_fact_Times;
#endif

int _findNoPart(int* idxVec,const int length)
{
  int maxID=0;
  for(int i=0;i<length;i++)
  	maxID = (maxID>idxVec[0])?maxID:idxVec[i];

  return maxID;
}

PartitionAugSolver::PartitionAugSolver(DoubleMatrix* MatIn, int* var_Part_idx_in,
					int *con_Part_idx_in, const int localNegaEigVal_in,
					const int fullVarSize, const int fullConSize, const int nPart)
	: 
	 var_Part_idx(var_Part_idx_in), con_Part_idx(con_Part_idx_in),
	 firstCallFlag(1), dkktSC(NULL), nb_part(nPart),
	 M_nb_Var(fullVarSize), M_nb_Con(fullConSize)
{
  if(nb_part==0){
  	nb_part = _findNoPart(var_Part_idx,M_nb_Var);
	int con_part = _findNoPart(con_Part_idx,M_nb_Con);
	nb_part = (nb_part>con_part)?nb_part:con_part;
  }

  kkt_diag.resize(nb_part,NULL);  
  diag_solver.resize(nb_part,NULL); 
  diag_Dim.resize(nb_part,0);  
  diag_Nnz.resize(nb_part,0); 
  diag_rowBeg.resize(nb_part,NULL);  
  diag_colIdx.resize(nb_part,NULL); 
  diag_ele.resize(nb_part,NULL);  
  diag_ele_Map.resize(nb_part,NULL);   
  diag_nb_Var.resize(nb_part,0);  
  diag_nb_Con.resize(nb_part,0);   

  diag_last_Dim=0;
  diag_last_Nnz=0;
  diag_last_nb_Var=0;
  diag_last_nb_Con=0;

  kkt_bord.resize(nb_part,NULL);  
  bord_Dim_m.resize(nb_part,0); 
  bord_Dim_n.resize(nb_part,0);
  bord_Nnz.resize(nb_part,0);  
  bord_rowBeg.resize(nb_part,NULL); 
  bord_colIdx.resize(nb_part,NULL);  
  bord_ele.resize(nb_part,NULL); 
  bord_ele_Map.resize(nb_part,NULL);  

  diag_IDinFull.resize(nb_part,NULL);  

  M_Dim = M_nb_Var+M_nb_Con;
  M_sys = dynamic_cast<SparseSymMatrix*>(MatIn);
  M_Nnz = M_sys->numberOfNonZeros();

  for(int i=0; i<M_nb_Var; i++){
  	int col_blockID = var_Part_idx[i];
  	if(col_blockID>0){
	  diag_Dim[col_blockID-1]++;
	  diag_nb_Var[col_blockID-1]++;
	}else{
	  diag_last_Dim++;
	  diag_last_nb_Var++;
	}
  }

  for(int i=0; i<M_nb_Con; i++){
  	int row_blockID = con_Part_idx[i];
  	if(row_blockID>0){
	  diag_Dim[row_blockID-1]++;
	  diag_nb_Con[row_blockID-1]++;
	}else{
	  diag_last_Dim++;
	  diag_last_nb_Con++;
	}  	
  }

  for(int i=0; i<nb_part; i++){
	bord_Dim_m[i] = diag_Dim[i];
	bord_Dim_n[i] = diag_last_Dim;
  }

  int sumDim=diag_last_Dim;
  for(int i=0; i<nb_part; i++){
	sumDim += diag_Dim[i];
  }
  assert(sumDim==M_Dim);

  sc_Dim = diag_last_Dim;
  assert(gBuildSchurComp==1);
  if(gBuildSchurComp == 1){
  	sc_Nnz = sc_Dim*sc_Dim;
	dkktSC = new DenseSymMatrix(sc_Dim);
	schur_solver =  dynamic_cast<DeSymIndefSolver*> (new DeSymIndefSolver(dkktSC));
  }
  assert(schur_solver);

  localNegaEigVal = localNegaEigVal_in;
}





template <typename T>
void FreeAll( T & t ) {
    T tmp;
    t.swap( tmp );
}

void PartitionAugSolver::firstCall()
{
  firstCallFlag = 0;

#ifdef TIMING
	  double tTot=MPI_Wtime();
	  MPI_Barrier(MPI_COMM_WORLD);
#endif

  std::vector<int*> nextdiag_InRow;
  std::vector<int*> nextbord_InRow;
  std::vector<int>  findID;
  nextdiag_InRow.resize(nb_part,NULL);
  nextbord_InRow.resize(nb_part,NULL); 

  int *nextdiag_last_InRow;
  std::vector<int> diagNnz_wrk, bordNnz_wrk;
  int wrkGoff; 


  //save the indeces for diagonal entries for a streamlined later update
  int* krowbegFullMat 	= M_sys->getStorageRef().krowM;
  int* jcolFullMat 		= M_sys->getStorageRef().jcolM;
  double* eleFullMat 	= M_sys->getStorageRef().M;

  for(int son=0;son<nb_part;son++){
	diag_rowBeg[son] = (int*) malloc((diag_Dim[son]+1)*sizeof(int));
	bord_rowBeg[son] = (int*) malloc((bord_Dim_m[son]+1)*sizeof(int));
	for(int jj=0;jj<bord_Dim_m[0]+1;jj++){
	  diag_rowBeg[son][jj]=0;
	  bord_rowBeg[son][jj]=0;
	} 
  }
  diag_last_rowBeg = (int*) malloc((diag_last_Dim+1)*sizeof(int));
  for(int jj=0;jj<diag_last_Dim+1;jj++) diag_last_rowBeg[jj]=0;

  
  for(int son=0;son<nb_part;son++){
    diag_IDinFull[son] = (int*) malloc(diag_Dim[son]*sizeof(int));
  } 
  diag_last_IDinFull = (int*) malloc(diag_last_Dim*sizeof(int));
  FullIDin_Diag   = (int*) malloc(M_Dim*sizeof(int));

//  for(int currVar=0; currVar<M_Dim;currVar++){
//  	FullVarIDin_Diag[currVar]=-1;
//  }

  // build var local-full map and build var full-local map
  findID.clear(); findID.resize(nb_part+1,0);
  for(int currVar=0; currVar<M_nb_Var;currVar++){
  	int blockID = var_Part_idx[currVar]-1;
	if(blockID >= 0){
	  FullIDin_Diag[currVar] = findID[blockID];
	  diag_IDinFull[blockID][findID[blockID]] = currVar;
	  findID[blockID]++;
	}else if(blockID==-1){
	  FullIDin_Diag[currVar] = findID[nb_part];
	  diag_last_IDinFull[findID[nb_part]] = currVar;
	  findID[nb_part]++;	  
	}
  }
  // build con local-full map and build con full-local map
  for(int currCon=0; currCon<M_nb_Con;currCon++){
  	int blockID = con_Part_idx[currCon]-1;
	if(blockID >= 0){
	  FullIDin_Diag[M_nb_Var+currCon] = findID[blockID];
	  diag_IDinFull[blockID][findID[blockID]] = M_nb_Var+currCon;
	  findID[blockID]++;
	}else if(blockID==-1){
	  FullIDin_Diag[M_nb_Var+currCon] = findID[nb_part];
	  diag_last_IDinFull[findID[nb_part]] = M_nb_Var+currCon;
	  findID[nb_part]++;	  
	}
  } 
  findID.clear(); FreeAll(findID);

  //count size  for Hessian part
  for(int currRow=0; currRow<M_nb_Var;currRow++){
	for( int k=krowbegFullMat[currRow];k<krowbegFullMat[currRow+1];k++){	
	  int currCol = jcolFullMat[k];
	  assert(currCol<M_nb_Var);
	  int colBlockID = var_Part_idx[currCol]-1;
	  int rowBlockID = var_Part_idx[currRow]-1;
	  int localColID = FullIDin_Diag[currCol];
	  int localRowID = FullIDin_Diag[currRow];

  	  if(colBlockID == rowBlockID && rowBlockID>=0){
		//belong to Diag[BlockID]
		int blockID = rowBlockID;
		diag_Nnz[blockID]++;
		if(localRowID >= localColID)
		  diag_rowBeg[blockID][localRowID+1]++;
		else
		  diag_rowBeg[blockID][localColID+1]++;
	  }else if(colBlockID == rowBlockID && rowBlockID==-1){
		//belong to Diag_last
		diag_last_Nnz++;
		if(localRowID >= localColID)
		  diag_last_rowBeg[localRowID+1]++;
		else
		  diag_last_rowBeg[localColID+1]++;
	  }else if(colBlockID==-1 && rowBlockID>=0){
		//belong to Bord[rowBlockID]
		int blockID = rowBlockID;
		bord_Nnz[blockID]++;
		if(localRowID >= localColID)
		  bord_rowBeg[blockID][localRowID+1]++;
		else
		  bord_rowBeg[blockID][localColID+1]++;
	  }else if(colBlockID>=0 && rowBlockID==-1){
		//belong to Bord[colBlockID]
		int blockID = colBlockID;
		bord_Nnz[blockID]++;
		if(localRowID >= localColID)
		  bord_rowBeg[blockID][localRowID+1]++;
		else
		  bord_rowBeg[blockID][localColID+1]++;
	  }else{
	    assert("impossible"&&0);
	  }
	}
  }

  int linksize =0 ;
  //count size  for Constraint part
  for(int currRow=M_nb_Var; currRow<M_Dim;currRow++){
	for( int k=krowbegFullMat[currRow];k<krowbegFullMat[currRow+1];k++){	
	  int currCol = jcolFullMat[k];
	  int consID = currRow-M_nb_Var;
	  int colBlockID;

	  int rowBlockID = con_Part_idx[consID]-1;
	  if(currCol >= M_nb_Var)
	  	colBlockID = rowBlockID;
	  else
	    colBlockID = var_Part_idx[currCol]-1;
	  
	  int localColID = FullIDin_Diag[currCol];
	  int localRowID = FullIDin_Diag[currRow];

  	  if(colBlockID == rowBlockID && rowBlockID>=0){
		//belong to Diag[BlockID]
		int blockID = rowBlockID;
		diag_Nnz[blockID]++;
		diag_rowBeg[blockID][localRowID+1]++;
	  }else if(colBlockID == rowBlockID && rowBlockID==-1){
		//belong to Diag_last
		diag_last_Nnz++;
		diag_last_rowBeg[localRowID+1]++;
	  }else if(colBlockID==-1 && rowBlockID>=0){
		//belong to Bord[rowBlockID]
		int blockID = rowBlockID;
		bord_Nnz[blockID]++;
		bord_rowBeg[blockID][localRowID+1]++;
	  }else if(colBlockID>=0 && rowBlockID==-1){
		//belong to Bord[colBlockID]
		int blockID = colBlockID;
		bord_Nnz[blockID]++;
		bord_rowBeg[blockID][localColID+1]++;
	  }else{
//	    assert("impossible"&&0);
		linksize++;
	  }
	}
  }

  int sumNNz=diag_last_Nnz;
//  assert(bord_Nnz[0]==0);
  for(int son=0;son<nb_part;son++){
    sumNNz += bord_Nnz[son]+diag_Nnz[son];
  } 
  assert(sumNNz == M_Nnz);

  for(int son=0;son<nb_part;son++){
    for(int j=0; j < diag_Dim[son]+1; j++){
	  diag_rowBeg[son][j] += diag_rowBeg[son][j-1];
    }  
    for(int j=0; j < bord_Dim_m[son]+1; j++){
	  bord_rowBeg[son][j] += bord_rowBeg[son][j-1];
    }

    // alocate space
    bord_colIdx[son] 	= (int*)malloc(bord_Nnz[son]*sizeof(int));
    bord_ele[son] 		= (double*)malloc(bord_Nnz[son]*sizeof(double));
    bord_ele_Map[son] 	= (int*)malloc(bord_Nnz[son]*sizeof(int));
    nextbord_InRow[son]  = (int*)malloc(bord_Dim_m[son]*sizeof(int));
	
    diag_colIdx[son]    = (int*)malloc(diag_Nnz[son]*sizeof(int)); 
    diag_ele[son] 		= (double*)malloc(diag_Nnz[son]*sizeof(double)); 
    diag_ele_Map[son] 	= (int*)malloc(diag_Nnz[son]*sizeof(int));  
    nextdiag_InRow[son]  = (int*)malloc(diag_Dim[son]*sizeof(int));

    for(int jj=0;jj<diag_Dim[son];jj++) nextdiag_InRow[son][jj] = diag_rowBeg[son][jj];
    for(int jj=0;jj<bord_Dim_m[son];jj++) nextbord_InRow[son][jj]= bord_rowBeg[son][jj];	
  } 


  for(int j=1; j < diag_last_Dim+1; j++){
	diag_last_rowBeg[j] += diag_last_rowBeg[j-1];
  }
  diag_last_colIdx	  = (int*)malloc(diag_last_Nnz*sizeof(int)); 
  diag_last_ele 	  = (double*)malloc(diag_last_Nnz*sizeof(double)); 
  diag_last_ele_Map   = (int*)malloc(diag_last_Nnz*sizeof(int));  
  nextdiag_last_InRow = (int*)malloc(diag_last_Dim*sizeof(int));
  for(int jj=0;jj<diag_last_Dim;jj++) nextdiag_last_InRow[jj] = diag_last_rowBeg[jj];


  
  diagNnz_wrk.resize(nb_part,0);
  bordNnz_wrk.resize(nb_part,0);
  wrkGoff = 0;
  int diag_last_Nnz_wrk=0;

  // get element of Hessian
  for(int currRow=0; currRow<M_nb_Var;currRow++){
	for( int k=krowbegFullMat[currRow];k<krowbegFullMat[currRow+1];k++){	
	  int currCol = jcolFullMat[k];
	  int colBlockID = var_Part_idx[currCol]-1;
	  int rowBlockID = var_Part_idx[currRow]-1;
	  int localColID = FullIDin_Diag[currCol];
	  int localRowID = FullIDin_Diag[currRow];

  	  if(colBlockID == rowBlockID && rowBlockID>=0){
		//belong to Diag[BlockID]
		int blockID = rowBlockID;
		assert(localRowID >= localColID);
		if(localRowID >= localColID){
		  wrkGoff = nextdiag_InRow[blockID][localRowID]++;
		  diag_colIdx[blockID][wrkGoff] = localColID;
		  diag_ele[blockID][wrkGoff]	 = eleFullMat[k];
		}	
		else{
		  wrkGoff = nextdiag_InRow[blockID][localColID]++;
		  diag_colIdx[blockID][wrkGoff] = localRowID;
		  diag_ele[blockID][wrkGoff]	= eleFullMat[k];
		}
		diag_ele_Map[blockID][wrkGoff] = k;
		diagNnz_wrk[blockID]++;
	  }else if(colBlockID == rowBlockID && rowBlockID==-1){
		//belong to Diag_last
	    assert(localRowID >= localColID);		
		if(localRowID >= localColID){
		  wrkGoff = nextdiag_last_InRow[localRowID]++;
		  diag_last_colIdx[wrkGoff] 	= localColID;
		  diag_last_ele[wrkGoff]		= eleFullMat[k];		  
		}
		else{
		  wrkGoff = nextdiag_last_InRow[localColID]++;
		  diag_last_colIdx[wrkGoff] = localRowID;
		  diag_last_ele[wrkGoff]	= eleFullMat[k];
		}
		diag_last_ele_Map[wrkGoff] = k;
		diag_last_Nnz_wrk++;		
	  }else if(colBlockID==-1 && rowBlockID>=0){
		//belong to Bord[rowBlockID]
		int blockID = rowBlockID;
		wrkGoff = nextbord_InRow[blockID][localRowID]++;
		bord_colIdx[blockID][wrkGoff]   = localColID;
		bord_ele[blockID][wrkGoff]	    = eleFullMat[k];
		bord_ele_Map[blockID][wrkGoff]  = k;
		bordNnz_wrk[blockID]++;		
	  }else if(colBlockID>=0 && rowBlockID==-1){
		//belong to Bord[colBlockID]
		int blockID = colBlockID;
		wrkGoff = nextbord_InRow[blockID][localColID]++;
		bord_colIdx[blockID][wrkGoff] 	= localRowID;
		bord_ele[blockID][wrkGoff]		= eleFullMat[k];
		bord_ele_Map[blockID][wrkGoff] = k;
		bordNnz_wrk[blockID]++;		  
	  }else{
	    assert("impossible"&&0);
	  }
	}
  }



  // get element of Jac
  for(int currRow=M_nb_Var; currRow<M_Dim;currRow++){
	for( int k=krowbegFullMat[currRow];k<krowbegFullMat[currRow+1];k++){	
	  int currCol = jcolFullMat[k];
	  int consID = currRow-M_nb_Var;
	  int colBlockID;

	  int rowBlockID = con_Part_idx[consID]-1;
	  if(currCol >= M_nb_Var)
	  	colBlockID = rowBlockID;
	  else
	    colBlockID = var_Part_idx[currCol]-1;
	  
	  int localColID = FullIDin_Diag[currCol];
	  int localRowID = FullIDin_Diag[currRow];

  	  if(colBlockID == rowBlockID && rowBlockID>=0){
		//belong to Diag[BlockID]
		int blockID = rowBlockID;
		assert(localRowID >= localColID);
		wrkGoff = nextdiag_InRow[blockID][localRowID]++;
		diag_colIdx[blockID][wrkGoff]  = localColID;
		diag_ele[blockID][wrkGoff]	   = eleFullMat[k];
		diag_ele_Map[blockID][wrkGoff] = k;
		diagNnz_wrk[blockID]++;
	  }else if(colBlockID == rowBlockID && rowBlockID==-1){
		//belong to Diag_last
  		assert(localRowID >= localColID);		
		wrkGoff = nextdiag_last_InRow[localRowID]++;
		diag_last_colIdx[wrkGoff] 	= localColID;
		diag_last_ele[wrkGoff]		= eleFullMat[k];		  
		diag_last_ele_Map[wrkGoff]  = k;
		diag_last_Nnz_wrk++;	
	  }else if(colBlockID==-1 && rowBlockID>=0){
		//belong to Bord[rowBlockID]
		int blockID = rowBlockID;
		wrkGoff = nextbord_InRow[blockID][localRowID]++;
		bord_colIdx[blockID][wrkGoff]   = localColID;
		bord_ele[blockID][wrkGoff]	    = eleFullMat[k];
		bord_ele_Map[blockID][wrkGoff]  = k;
		bordNnz_wrk[blockID]++;	
	  }else if(colBlockID>=0 && rowBlockID==-1){
		//belong to Bord[colBlockID]
		int blockID = colBlockID;
		wrkGoff = nextbord_InRow[blockID][localColID]++;
		bord_colIdx[blockID][wrkGoff] 	= localRowID;
		bord_ele[blockID][wrkGoff]		= eleFullMat[k];
		bord_ele_Map[blockID][wrkGoff]  = k;
		bordNnz_wrk[blockID]++;
	  }else{
	    assert("impossible"&&0);
	  }
	}
  }
  


  for(int jj=0;jj<diag_last_Dim;jj++) 
	assert(nextdiag_last_InRow[jj]==diag_last_rowBeg[jj+1]);
  assert (diag_last_Nnz_wrk == diag_last_Nnz);

  kkt_diag_last = new SparseSymMatrix( diag_last_Dim, diag_last_Nnz, diag_last_rowBeg, diag_last_colIdx, diag_last_ele);





  for(int son=0;son<nb_part;son++){
    for(int jj=0;jj<diag_Dim[son];jj++) 
	  assert(nextdiag_InRow[son][jj]==diag_rowBeg[son][jj+1]);	
    for(int jj=0;jj<bord_Dim_m[son];jj++) 
      assert(nextbord_InRow[son][jj]==bord_rowBeg[son][jj+1]);
    assert(bordNnz_wrk[son] == bord_Nnz[son] &&  diagNnz_wrk[son] == diag_Nnz[son]);	
  
    kkt_diag[son] = new SparseSymMatrix( diag_Dim[son], diag_Nnz[son], diag_rowBeg[son], diag_colIdx[son], diag_ele[son]);
    kkt_bord[son] = new SparseGenMatrix( bord_Dim_m[son], bord_Dim_n[son], bord_Nnz[son], 
  									bord_rowBeg[son], bord_colIdx[son], bord_ele[son]);

	diag_solver[son] = NULL;
    if(0==gSymLinearSolver){
#ifdef WITH_MA27		
	  diag_solver[son]	= new Ma27Solver( kkt_diag[son] );
#endif
    }
    else if(1==gSymLinearSolver){
#ifdef WITH_MA57		
	  diag_solver[son]	= new Ma57Solver( kkt_diag[son] );
#endif
    }
    else if(2==gSymLinearSolver){
#ifdef WITH_PARDISO
	  diag_solver[son]	= new PardisoSolver( kkt_diag[son], localNegaEigVal );
#endif	  
    } 
	assert(diag_solver[son]);
	
	free(nextdiag_InRow[son]);
    free(nextbord_InRow[son]);
  }

  free(nextdiag_last_InRow);

  FreeAll(nextdiag_InRow);
  FreeAll(nextbord_InRow);   

  diagNnz_wrk.clear(); FreeAll(diagNnz_wrk);
  bordNnz_wrk.clear(); FreeAll(bordNnz_wrk);

#ifdef TIMING
	  tTot=MPI_Wtime()-tTot;
	  MPI_Barrier(MPI_COMM_WORLD);
	  PartSolver_GenTime += tTot;
#endif

  
}



//faster than DenseSymMatrix::atPutZeros
void PartitionAugSolver::initializeKKT_Dense()
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




void PartitionAugSolver::finalizeKKT()
{
  int j, p, pend; double val;

  //alias for internal buffer of kkt
  double** dKkt = dkktSC->Mat();
  int SCdim = sc_Dim;
  
  /////////////////////////////////////////////////////////////
  // update the KKT with last diag block 
  /////////////////////////////////////////////////////////////
  int* krowQ = kkt_diag_last->krowM(); 
  int* jcolQ=kkt_diag_last->jcolM(); 
  double* dQ=kkt_diag_last->M();

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


int PartitionAugSolver::_numericalFact()
{
  if(1==firstCallFlag)
  	firstCall();

#ifdef TIMING
		  call_fact_Times++;
		  double tTot=MPI_Wtime();
		  MPI_Barrier(MPI_COMM_WORLD);
#endif


  double* eleFullMat  = M_sys->getStorageRef().M;

  //update matrices with new entries
  for(int wrkGoff=0; wrkGoff < diag_last_Nnz; wrkGoff++) {
	diag_last_ele[wrkGoff]	= eleFullMat[diag_last_ele_Map[wrkGoff]];
  }

  for(int block=0;block<nb_part;block++){
    for(int wrkGoff=0; wrkGoff < diag_Nnz[block]; wrkGoff++) {
	  diag_ele[block][wrkGoff]	= eleFullMat[diag_ele_Map[block][wrkGoff]];
    }
    for(int wrkGoff=0; wrkGoff < bord_Nnz[block]; wrkGoff++) {
	  bord_ele[block][wrkGoff]	= eleFullMat[bord_ele_Map[block][wrkGoff]];
    }  
  }
  int negEVal = 0;

  // set kktSC to 0;
  initializeKKT_Dense();
  
  // First tell children to factorize. 
  for(int block=0;block<nb_part;block++){
    negEVal += diag_solver[block]->matrixChanged();
  }
  
  // build schur complement 
  addTermToDenseSchurCompl();

  finalizeKKT();

  negEVal += schur_solver->matrixChanged();

#ifdef TIMING
		  tTot=MPI_Wtime()-tTot;
		  MPI_Barrier(MPI_COMM_WORLD);
		  PartSolver_FactTime += tTot;
#endif


  return negEVal;

};




/**
 * Computes C = Border^T * inv(Diag_1) * Border.
 */
void PartitionAugSolver::addTermToDenseSchurCompl() 
{

  int bord_m, bord_n, sc_mn;
  int blocksize = 64;
  
  sc_mn = dkktSC->size();

  for(int block=0;block<nb_part;block++){
    kkt_bord[block]->getSize(bord_m, bord_n); 
    assert(sc_mn>=bord_n);
	
	if(bord_n==-1) bord_n = sc_mn;
	assert(bord_n==sc_mn);

	DenseGenMatrix cols(blocksize,bord_m);

	for (int it=0; it < bord_n; it += blocksize) {
      int start=it;
      int end = MIN(it+blocksize,bord_n);
      int numcols = end-start;
      cols.getStorageRef().m = numcols; // avoid extra solves


      bool allzero = true;
      memset(&cols[0][0],0,bord_m*blocksize*sizeof(double));

	  kkt_bord[block]->getStorageRef().fromGetColBlock(start, &cols[0][0], bord_m, numcols, allzero);

      if(!allzero) {
	    diag_solver[block]->solve(cols);
	  	kkt_bord[block]->getStorageRef().transMultMat( 1.0, (*dkktSC)[start], numcols, bord_n, -1.0, &cols[0][0], bord_m);
      } 
    }
  }
}


/**
 * z0 -= Border^T  * Li\Di\ [ zi ]
 *
 */
void PartitionAugSolver::addLnizi(OoqpVector& z0_, OoqpVector& zi_,const int block )
{
  SimpleVector& z0 = dynamic_cast<SimpleVector&>(z0_);
  SimpleVector& zi = dynamic_cast<SimpleVector&>(zi_);

  diag_solver[block]->Dsolve (zi);  
  diag_solver[block]->Ltsolve(zi);

  //get n0= nx(parent)= #cols of A or C
  int dummy, n0;
  kkt_bord[block]->getSize(dummy, n0);

  kkt_bord[block]->transMult(1.0, z0, -1.0, zi);
}

void PartitionAugSolver::LDLtsolve_SC( OoqpVector* diag0_Rhs)
{
  SimpleVector& b0 = dynamic_cast<SimpleVector&>(*diag0_Rhs); 
  assert(gSolveSchurScheme==0);
  if(gSolveSchurScheme==0) {
    // Do Lsolve to Schur matrix	
	schur_solver->Lsolve(b0);
    // Option 1. - solve with the factors
    schur_solver->Dsolve(b0);
	schur_solver->Ltsolve(b0);
  }

}


/*
 *  y = alpha*Lni^T x + y
 *
 *                                  ( [ R 0 0 ]     )
 *  y = y + alpha* Di\Li\ (  [ A 0 0 ] * x )
 *                                  (  [ C 0 0 ]    )
 *
*/

void PartitionAugSolver::LniTransMult(SparseGenMatrix *borderMat, 
									SimpleVector& y, double alpha, SimpleVector& x, const int block)
{
  int N, nx0;
  
  //get nx(parent) from the number of cols of A (or C). Use N as dummy
  borderMat->getSize(N,nx0);
  // a mild assert
  assert(nx0 <= x.length());

  N = diag_Dim[block];
  assert( y.length() == N);

  //!memopt
  SimpleVector LniTx(N);
  SimpleVector x1(&x[0], nx0);
  
  LniTx.setToZero();
  borderMat->mult(0.0, LniTx, 1.0, x1);
  
  diag_solver[block]->Lsolve(LniTx); 
  diag_solver[block]->Dsolve(LniTx);
  diag_solver[block]->Ltsolve(LniTx);
  
  y.axpy(alpha,LniTx); 
 
}


void PartitionAugSolver::solve(GenMatrix& rhs_in)
{
  DenseGenMatrix &rhs = dynamic_cast<DenseGenMatrix&>(rhs_in);
  int N,NRHS;
  // rhs vectors are on the "rows", for continuous memory
  rhs.getSize(NRHS,N);
  assert(M_Dim==N);
    
  for (int i = 0; i < NRHS; i++) {
    SimpleVector v(rhs[i],N);
    solve(v);
  } 

}

void PartitionAugSolver::solve( OoqpVector& rhs_ )
{
  int block;
  
#ifdef TIMING
		double tTot=MPI_Wtime();
		MPI_Barrier(MPI_COMM_WORLD);
#endif

  
  SimpleVector& rhs = dynamic_cast<SimpleVector&>(rhs_);

  double *diag_last_Rhs_ele 		= (double*)malloc(diag_last_Dim*sizeof(double));
  std::vector<double*> diag_Rhs_ele;
  diag_Rhs_ele.resize(nb_part,NULL);

  for(block=0;block<nb_part;block++){
    diag_Rhs_ele[block] = (double*)malloc(diag_Dim[block]*sizeof(double));
  }


  
  int Index_Full;
  for(int Index_local=0;Index_local<diag_last_Dim;Index_local++){
	Index_Full = diag_last_IDinFull[Index_local];
	diag_last_Rhs_ele[Index_local] = rhs.elements()[Index_Full];
  }
  for(block=0;block<nb_part;block++){
	for(int Index_local=0;Index_local<diag_Dim[block];Index_local++){
	  Index_Full = diag_IDinFull[block][Index_local];
	  diag_Rhs_ele[block][Index_local] = rhs.elements()[Index_Full];
    }
  }
  
  SimpleVector  *diag_last_Rhs =  new SimpleVector( diag_last_Rhs_ele, diag_last_Dim );
  std::vector<SimpleVector*> diag_i_Rhs;
  diag_i_Rhs.resize(nb_part,NULL);

  //  Solve L
  for(block=0;block<nb_part;block++){
	diag_i_Rhs[block] =  new SimpleVector( diag_Rhs_ele[block], diag_Dim[block] );
	diag_solver[block]->Lsolve(*diag_i_Rhs[block]);
	addLnizi(*diag_last_Rhs,*diag_i_Rhs[block],block );
  }

  LDLtsolve_SC (diag_last_Rhs); 

  //  do backsolve to get x_i
  SimpleVector& x0 = *diag_last_Rhs; //just another name, for clarity
  
  // Li^T\bi for each child i. The backsolve needs z0
  for(block=0;block<nb_part;block++){
	diag_i_Rhs[block] =  new SimpleVector( diag_Rhs_ele[block], diag_Dim[block] );
    this->LniTransMult(kkt_bord[block], *diag_i_Rhs[block], -1.0, x0,block);	
  }

  // reorder the vector in orginal order
  for(int VarIDinDiag0=0;VarIDinDiag0<diag_last_Dim;VarIDinDiag0++){
	int IDinFull = diag_last_IDinFull[VarIDinDiag0];
	rhs.elements()[IDinFull] = diag_last_Rhs_ele[VarIDinDiag0];
  }

  for(block=0;block<nb_part;block++){
    for(int VarIDinDiag1=0;VarIDinDiag1<diag_Dim[block];VarIDinDiag1++){
	  int IDinFull = diag_IDinFull[block][VarIDinDiag1];
	  rhs.elements()[IDinFull] = diag_Rhs_ele[block][VarIDinDiag1];
    }
  }

  delete(diag_last_Rhs);
  for(block=0;block<nb_part;block++){
    delete(diag_i_Rhs[block]);
  }  
  FreeAll(diag_i_Rhs);
  
  free (diag_last_Rhs_ele);
  for(block=0;block<nb_part;block++){
    free(diag_Rhs_ele[block]);
  }
  FreeAll(diag_Rhs_ele);

#ifdef TIMING
		tTot=MPI_Wtime()-tTot;
		MPI_Barrier(MPI_COMM_WORLD);
		PartSolver_SolTime += tTot;
		call_sol_Times++;
#endif


}



