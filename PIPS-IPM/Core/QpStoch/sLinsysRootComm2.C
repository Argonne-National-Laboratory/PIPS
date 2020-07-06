/* PIPS
   Authors: Cosmin Petra and Miles Lubin
   See license and copyright information in the documentation */


#include "sLinsysRootComm2.h"
#include "sTree.h"
#include "sFactory.h"
#include "sData.h"
#include "sDummyLinsys.h"
#include "sLinsysLeaf.h"
#include "pipsport.h"

/*********************************************************************/
/************************** ROOT *************************************/
/*********************************************************************/

#ifdef STOCH_TESTING
extern double g_iterNumber;
extern double g_scenNum;
#endif

sLinsysRootComm2::sLinsysRootComm2(sFactory * factory_, sData * prob_)
  : sLinsysRoot(factory_, prob_)
{}

sLinsysRootComm2::sLinsysRootComm2(sFactory* factory_,
			 sData* prob_,
			 OoqpVector* dd_, 
			 OoqpVector* dq_,
			 OoqpVector* nomegaInv_,
			 OoqpVector* rhs_)
  : sLinsysRoot(factory_, prob_, dd_, dq_, nomegaInv_, rhs_)
{}

sLinsysRootComm2::~sLinsysRootComm2()
{ }

//this variable is just reset in this file; children will default to the "safe" linear solver
//extern int gLackOfAccuracy; unused

void sLinsysRootComm2::factor2(sData *prob, Variables *vars)
{
  DenseSymMatrix& kktd = dynamic_cast<DenseSymMatrix&>(*kkt);
  initializeKKT(prob, vars);

  // First tell children to factorize. 
  for(size_t c=0; c<children.size(); c++) {
    children[c]->factor2(prob->children[c], vars);
  }

  for(size_t c=0; c<children.size(); c++) {
#ifdef STOCH_TESTING
    g_scenNum=c;
#endif
    if(children[c]->mpiComm == MPI_COMM_NULL)
      continue;

    children[c]->stochNode->resMon.recFactTmChildren_start();    
    //---------------------------------------------
    children[c]->addTermToDenseSchurCompl(prob->children[c], kktd);
    //---------------------------------------------
    children[c]->stochNode->resMon.recFactTmChildren_stop();
  }

#ifdef TIMING
  MPI_Barrier(MPI_COMM_WORLD);
  stochNode->resMon.recReduceTmLocal_start();
#endif 
  reduceKKT();
 #ifdef TIMING
  stochNode->resMon.recReduceTmLocal_stop();
#endif  
  finalizeKKT(prob, vars);

  factorizeKKT();

#ifdef TIMING
  afterFactor();
#endif

  //gLackOfAccuracy=0;
}

void sLinsysRootComm2::Lsolve(sData *prob, OoqpVector& x)
{
  StochVector& b = dynamic_cast<StochVector&>(x);
  assert(children.size() == b.children.size() );
  int myRank; MPI_Comm_rank(mpiComm, &myRank);

  // children compute their part
  for(size_t it=0; it<children.size(); it++) {
    children[it]->Lsolve(prob->children[it], *b.children[it]);  
  }

  // Since a depth-first traversal is used, Li\bi is already done. 
  // Do the Schur compl and L0\b0

  SimpleVector& b0 = dynamic_cast<SimpleVector&>(*b.vec);

  //this code actually works on a single CPU too :)
  if (iAmDistrib) {
    //only one process add b0
    if(myRank>0) {
      b0.setToZero();
    }
  } //else b0.writeToStream(cout);


  for(size_t it=0; it<children.size(); it++) {
#ifdef TIMING
    children[it]->stochNode->resMon.eLsolve.clear();
    children[it]->stochNode->resMon.recLsolveTmChildren_start();
#endif
    SimpleVector& zi = dynamic_cast<SimpleVector&>(*b.children[it]->vec);

    //!memopt here
    //SimpleVector tmp(zi.length());
    //tmp.copyFromArray(zi.elements());
    //children[it]->addLnizi(prob->children[it], b0, tmp);
    children[it]->addLnizi(prob->children[it], b0, zi);
#ifdef TIMING
    children[it]->stochNode->resMon.recLsolveTmChildren_stop();
#endif
  }

#ifdef TIMING  
  MPI_Barrier(MPI_COMM_WORLD);
  stochNode->resMon.eReduce.clear();//reset
  stochNode->resMon.recReduceTmLocal_start();
#endif
  if (iAmDistrib) {
    if(0==myRank) {
      double* buffer = new double[b0.length()];
      MPI_Reduce(b0.elements(), buffer, b0.length(),
		 MPI_DOUBLE, MPI_SUM, 0, mpiComm);
      b0.copyFromArray(buffer);
      delete[] buffer;
    } else 
      MPI_Reduce(b0.elements(),nullptr,b0.length(),
		 MPI_DOUBLE, MPI_SUM, 0, mpiComm);
  }
#ifdef TIMING 
  stochNode->resMon.recReduceTmLocal_stop();
#endif

#ifdef TIMING
  stochNode->resMon.eLsolve.clear();
  stochNode->resMon.recLsolveTmLocal_start();
#endif

  //no reduce/broadcast  -  this function does not do anything
  if(0==myRank)
    solver->Lsolve(b0);
#ifdef TIMING
  stochNode->resMon.recLsolveTmLocal_stop();
#endif

}

void sLinsysRootComm2::Dsolve( sData *prob, OoqpVector& x )
{
  StochVector& b = dynamic_cast<StochVector&>(x);
  int myRank; MPI_Comm_rank(mpiComm,&myRank);

  //! commented - already done in addLnizi - cpetra
  //  for(size_t it=0; it<children.size(); it++) {
  //  children[it]->Dsolve(prob->children[it], *b.children[it]);
  //}

  SimpleVector& b0 = dynamic_cast<SimpleVector&>(*b.vec);
#ifdef TIMING
  stochNode->resMon.eDsolve.clear();
  stochNode->resMon.recDsolveTmLocal_start();
#endif
  //reduce/broadcast - broadcast b0 from rank 0
  if(0==myRank)
    solveReduced(prob, b0);
#ifdef TIMING
  stochNode->resMon.recDsolveTmLocal_stop();
#endif

#ifdef TIMING
  MPI_Barrier(mpiComm);
  stochNode->resMon.eBcast.clear();
  stochNode->resMon.recBcastTmLocal_start();
#endif
  if(iAmDistrib) {
    MPI_Bcast(b0.elements(), b0.length(), MPI_DOUBLE, 0, mpiComm);
  }
#ifdef TIMING
  stochNode->resMon.recBcastTmLocal_stop();
#endif
}

void sLinsysRootComm2::Ltsolve( sData *prob, OoqpVector& x )
{
  StochVector& b   = dynamic_cast<StochVector&>(x);
  SimpleVector& b0 = dynamic_cast<SimpleVector&>(*b.vec);

#ifdef TIMING
  stochNode->resMon.eLtsolve.clear();
  stochNode->resMon.recLtsolveTmLocal_start();
#endif
  solver->Ltsolve(b0);
#ifdef TIMING
  stochNode->resMon.recLtsolveTmLocal_stop();
#endif

  SimpleVector& x0 = b0; //just another name, for clarity
  
  // Li^T\bi for each child i. The backsolve needs z0

  for(size_t it=0; it<children.size(); it++) {
    children[it]->Ltsolve2(prob->children[it], *b.children[it], x0);
  }
#ifdef TIMING
  int myRank; MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
  if(128*(myRank/128) == myRank) {
    double tTotResChildren=0.0;
    for(size_t it=0; it<children.size(); it++) {
      tTotResChildren += children[it]->stochNode->resMon.eLsolve.tmChildren;
      tTotResChildren += children[it]->stochNode->resMon.eLsolve.tmLocal;
    }
    double tReduce=stochNode->resMon.eReduce.tmLocal;
    double tBcast=stochNode->resMon.eBcast.tmLocal; 
    double tStg1=stochNode->resMon.eDsolve.tmLocal;  

    double tTotStg2Children=0.0;
    for(size_t it=0; it<children.size(); it++) {
      tTotStg2Children += children[it]->stochNode->resMon.eLtsolve.tmChildren;
      tTotStg2Children += children[it]->stochNode->resMon.eLtsolve.tmLocal;
    }
    cout << "  rank " << myRank << " "
	 << "Resid comp " << tTotResChildren << " " 
	 << "reduce " << tReduce << " "
	 << "bcast " << tBcast << " " 
	 << "1stStage solve " << tStg1 << " "
	 << "2ndStage solve " << tTotStg2Children <<endl;
  }
#endif


}

///////////////////////////////////////////////////////////
// ATOMS of FACTOR 2
//////////////////////////////////////////////////////////
  /* Atoms methods of FACTOR2 for a non-leaf linear system */

void sLinsysRootComm2::reduceKKT()
{
  DenseSymMatrix* kktd = dynamic_cast<DenseSymMatrix*>(kkt); 

  //parallel communication
  if(iAmDistrib) submatrixReduce(kktd, 0, 0, locnx, locnx, mpiComm);
}


void sLinsysRootComm2::factorizeKKT()
{
  //stochNode->resMon.recFactTmLocal_start();  
#ifdef TIMING
  extern double g_iterNumber;
  double st=MPI_Wtime();
#endif

  int myRank; MPI_Comm_rank(mpiComm,&myRank);
  if(0==myRank)
    solver->matrixChanged();

  //stochNode->resMon.recFactTmLocal_stop(); 
#ifdef TIMING
  st = MPI_Wtime()-st;
  // note, this will include noop scalapack processors
  if( 0==myRank )
    printf("  rank %d 1stSTAGE FACT %g SEC ITER %d\n", myRank, st, (int)g_iterNumber);
#endif
}

 
#define CHUNK_SIZE 1024*1024*64 //doubles  = 128 MBytes (maximum)
void sLinsysRootComm2::submatrixReduce(DenseSymMatrix* A, 
				     int row, int col, int drow, int dcol,
				     MPI_Comm comm)
{
  double ** M = A->mStorage->M;
  int n = A->mStorage->n;
#ifdef DEBUG 
  assert(n >= row+drow);
  assert(n >= col+dcol);
#endif
  int chunk_size = CHUNK_SIZE / n * n; 
  chunk_size = min(chunk_size, n*n);
  double* chunk = new double[chunk_size];

  int rows_in_chunk = chunk_size/n;
  int iRow=row;
  do {

    if(iRow+rows_in_chunk > drow)
      rows_in_chunk = drow-iRow;

    //iErr=MPI_Allreduce(&M[iRow][0], 
    //		       chunk, rows_in_chunk*n, 
    //		       MPI_DOUBLE, MPI_SUM, comm);
    MPI_Reduce(&M[iRow][0], chunk, rows_in_chunk*n,
		    MPI_DOUBLE, MPI_SUM, 0, comm);

    //copy data in M
    for(int i=iRow; i<iRow+rows_in_chunk; i++) {

      int shft = (i-iRow)*n;
      for(int j=col; j<col+dcol; j++)
	M[i][j] = chunk[shft+j];
    }
    iRow += rows_in_chunk;
  
  } while(iRow<row+drow);

  delete[] chunk;
}


