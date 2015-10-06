/* PIPS
   Authors: Cosmin Petra
   See license and copyright information in the documentation */

/* 2015. Modified by Nai-Yuan Chiang for NLP*/

#ifndef STOCHLEAFLINSYS
#define STOCHLEAFLINSYS

#include "sLinsys.h"
#include "sTree.h"
#include "sFactory.h"
#include "sData.h"
#include "SparseSymMatrix.h"
#include "SparseGenMatrix.h"


//#include "MtxSchurDecompSolver.h"
//#include "ReducedSpaceSolver.h"
#include "ReducedSpaceSolverStateOnly.h"
//#include "PartitionAugSolver.h"


extern int gOuterSolve;
extern int separateHandDiag;

extern int gBuildSchurComp;
extern int gUseReducedSpace;
extern int gNP_Alg;

/** This class solves the linear system corresponding to a leaf node.
 *  It just redirects the call to NlpGenSparseLinsys.
 */
class sLinsysLeaf : public sLinsys
{
 public:
  //sLinsysLeaf(NlpGenStoch * factory_, sData * prob_);
  template<class LINSOLVER>
    sLinsysLeaf(sFactory* factory,
		sData* prob_,				    
		OoqpVector* dd_, OoqpVector* dq_, OoqpVector* nomegaInv_,
		OoqpVector* rhs_, OoqpVector* additiveDiag_, LINSOLVER *linsolver=NULL);

  virtual ~sLinsysLeaf();

  virtual int factor2( sData *prob, Variables *vars);
  virtual void Lsolve ( sData *prob, OoqpVector& x );
  virtual void Dsolve ( sData *prob, OoqpVector& x );
  virtual void Ltsolve( sData *prob, OoqpVector& x );

  //virtual void Lsolve2 ( OoqpVector& x );
  //virtual void Dsolve2 ( OoqpVector& x );
  virtual void Ltsolve2( sData *prob, StochVector& x, SimpleVector& xp);


  //virtual void solveCompressed( OoqpVector& rhs );
  virtual void putXDiagonal( OoqpVector& xdiag_ );
  virtual void putSDiagonal( OoqpVector& sdiag_ );
  virtual void putYDualDiagonal( OoqpVector& ydiag_ );
  virtual void putZDiagonal( OoqpVector& zdiag );

  virtual void setAdditiveDiagonal();


  //void Ltsolve_internal(  sData *prob, StochVector& x, SimpleVector& xp);
  void sync();
  virtual void deleteChildren();
  virtual void UpdateMatrices( Data * prob_in,int const updateLevel=2);


 
 bool firstQUpdate,firstBUpdate, firstDUpdate;
 
 std::map<int,int> LocQMap;
 std::map<int,int> LocBMap;
 std::map<int,int> LocDMap;


 bool firstXDiagUpdate,firstSDiagUpdate,firstYDiagUpdate,firstZDiagUpdate;

 std::map<int,int> xDiagIdxMap;
 std::map<int,int> sDiagIdxMap;
 std::map<int,int> yDiagIdxMap;
 std::map<int,int> zDiagIdxMap;

 virtual void setXDiagonal( OoqpVector& xdiag );
 virtual void setSDiagonal( OoqpVector& sdiag );  
 virtual void setYDiagonal( OoqpVector& ydiag );  
 virtual void setZDiagonal( OoqpVector& zdiag );

  

 protected:
  sLinsysLeaf() {};

  static void mySymAtPutSubmatrix(SymMatrix& kkt, 
				  GenMatrix& B, GenMatrix& D, 
				  int locnx, int locmy, int locmz);
}; 

template<class LINSOLVER>
sLinsysLeaf::sLinsysLeaf(sFactory *factory_, sData* prob,
			 OoqpVector* dd_, 
			 OoqpVector* dq_,
			 OoqpVector* nomegaInv_,
			 OoqpVector* rhs_,
			 OoqpVector* additiveDiag_,
			 LINSOLVER* thesolver)
  : sLinsys(factory_, prob, dd_, dq_, nomegaInv_, rhs_, additiveDiag_)
{
  //int rank; MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  //double t = MPI_Wtime();

  // create the KKT system matrix
  // size = ?  nnz = ?


#ifdef NY_TIMING
	double time_Temp_1, time_Temp_2;
	time_Temp_1=MPI_Wtime();
	int mype_;
	MPI_Comm_rank(MPI_COMM_WORLD,&mype_);
#endif	  

  
  int nnzQ, nnzB, nnzD;
  int n, nnzTotal;
  int locns;

  firstQUpdate=true;
  firstBUpdate=true; 
  firstDUpdate=true;
  firstXDiagUpdate=true;
  firstSDiagUpdate=true;
  firstYDiagUpdate=true;
  firstZDiagUpdate=true;

  prob->getLocalSizes(locnx, locmy, locmz);
  locns = locmz;
  prob->getLocalNnz(nnzQ, nnzB, nnzD);

  if(gOuterSolve >= 3 && separateHandDiag==0){
	n = locnx+locns+locmy+locmz;	  // for x s y z
	nnzTotal = n + nnzQ + nnzB + nnzD + locns;
  }
  else if(gOuterSolve >= 3 && separateHandDiag==1){
	n = locnx+locns+locmy+locmz;		// for x s y z
	nnzTotal = nnzQ + nnzB + nnzD + locns;
  }
  else {
	n = locnx+locmy+locmz;		  // only x y z, we need compress the linear system later
	nnzTotal = n + nnzQ + nnzB + nnzD;
  }


  //alocate the matrix and copy the data into
  SparseSymMatrix* kktsp = new SparseSymMatrix(n,nnzTotal);
  kkt = kktsp;

  if(gOuterSolve >= 3){
	if(separateHandDiag==0){
      SimpleVectorHandle v( new SimpleVector(n) );
      v->setToZero();
      kkt->setToDiagonal(*v);
	}
	
#ifdef NY_TIMING
		 time_Temp_2=MPI_Wtime();
		 if(0==mype_){
		  cout << "\n before put matrix in sLinsysLeaf " << time_Temp_2- time_Temp_1<< "\n" << endl;
		  }
		 time_Temp_1=MPI_Wtime();
#endif	  

	kkt->symAtPutSubmatrix( 0, 0, prob->getLocalQ(), 0, 0, locnx, locnx);

#ifdef NY_TIMING
			 time_Temp_2=MPI_Wtime();
			 if(0==mype_){
			  cout << "\n put Q in sLinsysLeaf " << time_Temp_2- time_Temp_1<< "\n" << endl;
			  }
			 time_Temp_1=MPI_Wtime();
#endif	


  	kkt->symAtPutSubmatrix( locnx+locns, 0, prob->getLocalB(), 0, 0, locmy, locnx);

#ifdef NY_TIMING
			 time_Temp_2=MPI_Wtime();
			 if(0==mype_){
			  cout << "\n put A in sLinsysLeaf " << time_Temp_2- time_Temp_1<< "\n" << endl;
			  }
#endif
	  
	
	// Fix the problem in shiftRows_CorrectMap in SparseStorage, then we can replace symAtPutSubmatrix 
	// by symAtSetSubmatrix e.g:
    // kkt->symAtSetSubmatrix( 0, 0, prob->getLocalQ(), 0, 0, locnx, locnx,firstQUpdate, LocQMap);	
    // and others		

	if (locns>0){

#ifdef NY_TIMING
			 time_Temp_1=MPI_Wtime();
#endif			
      kkt->symAtPutSubmatrix( locnx+locns+locmy, 0, prob->getLocalD(), 0, 0, locmz, locnx);
#ifdef NY_TIMING
				   time_Temp_2=MPI_Wtime();
				   if(0==mype_){
					cout << "\n put C in sLinsysLeaf " << time_Temp_2- time_Temp_1<< "\n" << endl;
					}
#endif
			
#ifdef NY_TIMING
				   time_Temp_1=MPI_Wtime();
#endif	


	  SparseSymMatrixHandle tempDiagMat( new SparseSymMatrix( locns, locns ) );
	  int *tempDiagRowId = new int[locns];
	  int *tempDiagColId = new int[locns];
	  double *tempDiagEleId = new double[locns];
	  int info;

	  for(int i=0;i<locns;i++){
	    tempDiagRowId[i]=i;
	    tempDiagColId[i]=i;
	    tempDiagEleId[i]=-1;
	   }
	  tempDiagMat->putSparseTriple( tempDiagRowId,locns, tempDiagColId, tempDiagEleId,info );
	  kkt->symAtPutSubmatrix( locnx + locns + locmy, locnx, *tempDiagMat, 0, 0, locns, locns);

#ifdef NY_TIMING
				   time_Temp_2=MPI_Wtime();
				   if(0==mype_){
					cout << "\n put Slack Diag Mat in sLinsysLeaf " << time_Temp_2- time_Temp_1<< "\n" << endl;
					}
#endif
			
#ifdef NY_TIMING
				   time_Temp_1=MPI_Wtime();
#endif    
      delete[] tempDiagRowId;
      delete[] tempDiagColId;
      delete[] tempDiagEleId;
    }
  }
  else{
    SimpleVectorHandle v( new SimpleVector(n) );
    v->setToZero();
    kkt->setToDiagonal(*v);

	kkt->symAtPutSubmatrix( 0, 0, prob->getLocalQ(), 0, 0, locnx, locnx);
	if(locmz>0) {
	  kkt->symAtPutSubmatrix( locnx, 0, prob->getLocalB(), 0, 0, locmy, locnx);
	  kkt->symAtPutSubmatrix( locnx+locmy, 0, prob->getLocalD(), 0, 0, locmz, locnx);		  
	} else{
	  mySymAtPutSubmatrix(*kkt, prob->getLocalB(), prob->getLocalD(), locnx, locmy, locmz);
	}
  }
  // Fix the problem in shiftRows_CorrectMap in SparseStorage, then we can use the following codes.  
//  firstQUpdate = false;
//  firstBUpdate = false;
//  firstDUpdate = false;  

  // create the solver for the linear system

			  
#ifdef NY_TIMING
					 time_Temp_1=MPI_Wtime();
#endif 

  
  solver = NULL;
  if(0==gUseReducedSpace && gNP_Alg==0){
	if(1==gBuildSchurComp){
	  solver = new LINSOLVER( kktsp,locmy+locmz);
	}else{
	  // build sc to solve problem, we only support dense schur comp
	  assert("Not Implemented" &&0);
	}
  }
  else if (1==gUseReducedSpace && gNP_Alg==0 ){
  	// use reduced space schur solver, assume decition vars is inside kktsp
	if(0==gBuildSchurComp){
	  assert("Not Implemented" &&0);
	}else if(1==gBuildSchurComp){
	  //FIXME: prob->schursize should retunr local size
	  int decisionVarSize = prob->schurSize;
	  int *decisionVarID = prob->schurVarConID;

	  int fullVarXSize = locnx;
	  int fullVarYSize = locmy;
	  int fullVarSSize = locmz;

//	  solver = new ReducedSpaceSolver(kktsp,decisionVarSize,decisionVarID,fullVarXSize,fullVarYSize,fullVarSSize);
	}
	
  }
  else if(2==gUseReducedSpace && gNP_Alg==0){
  	// use reduced space schur solver, assume all the decition vars have already been moved
	solver = new ReducedSpaceSolverStateOnly(kktsp,locnx,locmy,locmz);
  }
  else if(gNP_Alg>0){
//	solver = new PartitionAugSolver(kktsp,  
//					prob->var_Part_idx_in,
//					prob->con_Part_idx_in, 
//					locmy+locmz, locnx+locmz, locmy+locmz, 0);
		
  }


#ifdef NY_TIMING
					 time_Temp_2=MPI_Wtime();
					 if(0==mype_){
					  cout << "\n creat solver in sLinsysLeaf " << time_Temp_2- time_Temp_1<< "\n" << endl;
					  }
#endif
			  
#ifdef NY_TIMING
					 time_Temp_1=MPI_Wtime();
#endif

  
//  solver = new LINSOLVER(kktsp);		
  assert(solver != NULL);
  
  //t = MPI_Wtime() - t;
  //if (rank == 0) printf("new sLinsysLeaf took %f sec\n",t);

  mpiComm = (dynamic_cast<StochVector*>(dd_))->mpiComm;
}



#endif
