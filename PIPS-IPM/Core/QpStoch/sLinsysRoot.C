/* PIPS
   Authors: Cosmin Petra and Miles Lubin
   See license and copyright information in the documentation */


#include "sLinsysRoot.h"
#include "sTree.h"
#include "sFactory.h"
#include "sData.h"
#include "sDummyLinsys.h"
#include "sLinsysLeaf.h"
/*********************************************************************/
/************************** ROOT *************************************/
/*********************************************************************/

#ifdef STOCH_TESTING
extern double g_iterNumber;
#endif

sLinsysRoot::sLinsysRoot(sFactory * factory_, sData * prob_)
  : sLinsys(factory_, prob_), iAmDistrib(0)
{
  createChildren(prob_);
}

sLinsysRoot::sLinsysRoot(sFactory* factory_,
			 sData* prob_,
			 OoqpVector* dd_, 
			 OoqpVector* dq_,
			 OoqpVector* nomegaInv_,
			 OoqpVector* rhs_)
  : sLinsys(factory_, prob_, dd_, dq_, nomegaInv_, rhs_), iAmDistrib(0)
{
  createChildren(prob_);
}

sLinsysRoot::~sLinsysRoot()
{
  for(size_t c=0; c<children.size(); c++)
    delete children[c];
}

void sLinsysRoot::factor2(sData *prob, Variables *vars)
{
  DenseSymMatrix& kktd = dynamic_cast<DenseSymMatrix&>(*kkt);
  initializeKKT(prob, vars);
  
  // First tell children to factorize.
  for(size_t c=0; c<children.size(); c++) {
    children[c]->factor2(prob->children[c], vars);
  }

  for(size_t c=0; c<children.size(); c++) {

    if(children[c]->mpiComm == MPI_COMM_NULL)
      continue;

    children[c]->stochNode->resMon.recFactTmChildren_start();    
    //---------------------------------------------
    children[c]->addTermToDenseSchurCompl(prob->children[c], kktd);
    //---------------------------------------------
    children[c]->stochNode->resMon.recFactTmChildren_stop();
  }
  stochNode->resMon.recReduceTmLocal_start();
  
  reduceKKT();


  stochNode->resMon.recReduceTmLocal_stop();
  
  finalizeKKT(prob, vars);
  
  //printf("(%d, %d) --- %f\n", PROW,PCOL, kktd[PROW][PCOL]);

  factorizeKKT();

  //int mype; MPI_Comm_rank(MPI_COMM_WORLD,&mype);
  //if (mype==0) dumpMatrix(-1, 0, "kkt", kktd); 
  //MPI_Barrier(MPI_COMM_WORLD);
#ifdef TIMING
  afterFactor();
#endif

}

#ifdef TIMING
void sLinsysRoot::afterFactor()
{
  int mype; MPI_Comm_rank(mpiComm, &mype);

  for (size_t c=0; c<children.size(); c++) {
    if (children[c]->mpiComm == MPI_COMM_NULL) continue;
    printf("NODE %4zu SPFACT %g BACKSOLVE %g SEC PROC %d ITER %d\n", c,
     children[c]->stochNode->resMon.eFact.tmLocal,
	   children[c]->stochNode->resMon.eFact.tmChildren, mype, (int)g_iterNumber);
  }
  double redall = stochNode->resMon.eReduce.tmLocal;
  double redscat = stochNode->resMon.eReduceScatter.tmLocal;
  printf("REDUCE %g SEC ON PROC %d ITER %d REDSCAT %g DIFF %g\n", redall, 
    mype, (int)g_iterNumber, redscat, redall-redscat);
}
#endif

void sLinsysRoot::Lsolve(sData *prob, OoqpVector& x)
{
  StochVector& b = dynamic_cast<StochVector&>(x);
  assert(children.size() == b.children.size() );

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
    int myRank; MPI_Comm_rank(mpiComm, &myRank);
    if(myRank>0) {
      b0.setToZero();
    }
  } //else b0.writeToStream(cout);


  for(size_t it=0; it<children.size(); it++) {
    children[it]->stochNode->resMon.recLsolveTmChildren_start();

    SimpleVector& zi = dynamic_cast<SimpleVector&>(*b.children[it]->vec);

    //!memopt here
    SimpleVector tmp(zi.length());
    tmp.copyFromArray(zi.elements());

    children[it]->addLnizi(prob->children[it], b0, tmp);

    children[it]->stochNode->resMon.recLsolveTmChildren_stop();
  }

  if (iAmDistrib) {
 
    double* buffer = new double[b0.length()];

    MPI_Allreduce(b0.elements(), buffer, b0.length(),
		  MPI_DOUBLE, MPI_SUM, mpiComm);

    b0.copyFromArray(buffer);

    delete[] buffer;
  }

  //dumpRhs(0, "rhs",  b0);

  stochNode->resMon.recLsolveTmLocal_start();
  solver->Lsolve(b0);
  stochNode->resMon.recLsolveTmLocal_stop();
}


void sLinsysRoot::Ltsolve2( sData *prob, StochVector& x, SimpleVector& xp)
{
  StochVector& b   = dynamic_cast<StochVector&>(x);
  SimpleVector& bi = dynamic_cast<SimpleVector&>(*b.vec);
  
  stochNode->resMon.recLtsolveTmLocal_start();

  //b_i -= Lni^T x0
  this->LniTransMult(prob, bi, -1.0, xp);
  solver->Ltsolve(bi);

  stochNode->resMon.recLtsolveTmLocal_stop();

  SimpleVector& xi = bi;
  //recursive call in order to get the children to do their part
  for(size_t it=0; it<children.size(); it++) {
    children[it]->Ltsolve2(prob->children[it], *b.children[it], xi);
  }
}

void sLinsysRoot::Ltsolve( sData *prob, OoqpVector& x )
{
  StochVector& b   = dynamic_cast<StochVector&>(x);
  SimpleVector& b0 = dynamic_cast<SimpleVector&>(*b.vec);

  stochNode->resMon.recLtsolveTmLocal_start();
  solver->Ltsolve(b0);
  stochNode->resMon.recLtsolveTmLocal_stop();

  //dumpRhs(0, "sol",  b0);

  SimpleVector& x0 = b0; //just another name, for clarity
  
  // Li^T\bi for each child i. The backsolve needs z0

  for(size_t it=0; it<children.size(); it++) {
    children[it]->Ltsolve2(prob->children[it], *b.children[it], x0);
  }
}

void sLinsysRoot::Dsolve( sData *prob, OoqpVector& x )
{
  StochVector& b = dynamic_cast<StochVector&>(x);

  for(size_t it=0; it<children.size(); it++) {
    children[it]->Dsolve(prob->children[it], *b.children[it]);
  }

  SimpleVector& b0 = dynamic_cast<SimpleVector&>(*b.vec);
  solveReduced(prob, b0);
}



void sLinsysRoot::createChildren(sData* prob)
{
	sLinsys* child=NULL;
  StochVector& ddst = dynamic_cast<StochVector&>(*dd);
  StochVector& dqst = dynamic_cast<StochVector&>(*dq);
  StochVector& nomegaInvst = dynamic_cast<StochVector&>(*nomegaInv);
  StochVector& rhsst = dynamic_cast<StochVector&>(*rhs);

  //get the communicator from one of the vectors
  this->mpiComm = ddst.mpiComm;
  this->iAmDistrib = ddst.iAmDistrib;
  for(size_t it=0; it<prob->children.size(); it++) {
		if(MPI_COMM_NULL == ddst.children[it]->mpiComm) {
      child = new sDummyLinsys(dynamic_cast<sFactory*>(factory), prob->children[it]);
    } else {
      sFactory* stochFactory = dynamic_cast<sFactory*>(factory);

      if(prob->children[it]->children.size() == 0) {	
	//child = new sLinsysLeaf(dynamic_cast<QpGenStoch*>(factory), 
	child = stochFactory->newLinsysLeaf(prob->children[it],
					     ddst.children[it],
					     dqst.children[it],
					     nomegaInvst.children[it],
					     rhsst.children[it]);
      } else {
	//child = new sLinsysRoot(dynamic_cast<QpGenStoch*>(factory), 
	child = stochFactory->newLinsysRoot(prob->children[it],
					     ddst.children[it],
					     dqst.children[it],
					     nomegaInvst.children[it],
					     rhsst.children[it]);
      }
    }
    AddChild(child);
  }
}

void sLinsysRoot::deleteChildren()
{
  for(size_t it=0; it<children.size(); it++) {
    children[it]->deleteChildren();
    delete children[it];
  }
  children.clear();
}

void sLinsysRoot::putXDiagonal( OoqpVector& xdiag_ )
{
  StochVector& xdiag = dynamic_cast<StochVector&>(xdiag_);
  assert(children.size() == xdiag.children.size());

  //kkt->atPutDiagonal( 0, *xdiag.vec );
  xDiag = xdiag.vec;
 
  // propagate it to the subtree
  for(size_t it=0; it<children.size(); it++)
    children[it]->putXDiagonal(*xdiag.children[it]);
}


void sLinsysRoot::putZDiagonal( OoqpVector& zdiag_ )
{
  StochVector& zdiag = dynamic_cast<StochVector&>(zdiag_);
  assert(children.size() == zdiag.children.size());

  //kkt->atPutDiagonal( locnx+locmy, *zdiag.vec );
  zDiag = zdiag.vec;

  // propagate it to the subtree
  for(size_t it=0; it<children.size(); it++)
    children[it]->putZDiagonal(*zdiag.children[it]);
}

void sLinsysRoot::AddChild(sLinsys* child)
{
  children.push_back(child);
}


void sLinsysRoot::sync()
{
  //delete children
  deleteChildren();
  //assert(false);

  //delete local stuff
  if( nxupp + nxlow > 0 ) {
    delete dd; 
    delete dq; 
  }
  delete nomegaInv;
  delete rhs;
  if (solver) delete solver;
  if (kkt)    delete kkt;


  //allocate
  if( nxupp + nxlow > 0 ) {
    //dd      = OoqpVectorHandle(stochNode->newPrimalVector());
    //dq      = OoqpVectorHandle(stochNode->newPrimalVector());
    dd = stochNode->newPrimalVector();
    dq = stochNode->newPrimalVector();
    data->getDiagonalOfQ( *dq );
  }
  nomegaInv   = stochNode->newDualZVector();
  rhs         = stochNode->newRhs();


  data->getLocalSizes(locnx, locmy, locmz);
  createChildren(data);

  kkt = createKKT(data);
  solver = createSolver(data, kkt);
}

///////////////////////////////////////////////////////////
// ATOMS of FACTOR 2
//////////////////////////////////////////////////////////
  /* Atoms methods of FACTOR2 for a non-leaf linear system */
void sLinsysRoot::initializeKKT(sData* prob, Variables* vars)
{
  DenseSymMatrix* kktd = dynamic_cast<DenseSymMatrix*>(kkt); 
  myAtPutZeros(kktd);
}

void sLinsysRoot::reduceKKT()
{
  DenseSymMatrix* kktd = dynamic_cast<DenseSymMatrix*>(kkt); 
  //parallel communication
  if(iAmDistrib) submatrixAllReduce(kktd, 0, 0, locnx, locnx, mpiComm);
}


void sLinsysRoot::factorizeKKT()
{
  //stochNode->resMon.recFactTmLocal_start();  
#ifdef TIMING
  MPI_Barrier(mpiComm);
  extern double g_iterNumber;
  double st=MPI_Wtime();
#endif

  solver->matrixChanged();

  //stochNode->resMon.recFactTmLocal_stop(); 
#ifdef TIMING
  st = MPI_Wtime()-st;
  MPI_Barrier(mpiComm);
  int mype; MPI_Comm_rank(mpiComm, &mype);
  // note, this will include noop scalapack processors
  printf("FACT %g SEC ON PROC %d ITER %d\n", st, mype, (int)g_iterNumber);
#endif
}

 

//faster than DenseSymMatrix::atPutZeros
void sLinsysRoot::myAtPutZeros(DenseSymMatrix* mat, 
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

void sLinsysRoot::myAtPutZeros(DenseSymMatrix* mat)
{
  int n = mat->size();
  myAtPutZeros(mat, 0, 0, n, n);
}


#define CHUNK_SIZE 1024*1024*16 //doubles  = 128 MBytes (maximum)
void sLinsysRoot::submatrixAllReduce(DenseSymMatrix* A, 
				     int row, int col, int drow, int dcol,
				     MPI_Comm comm)
{
  double ** M = A->mStorage->M;
  int n = A->mStorage->n;
#ifdef DEBUG 
  assert(n >= row+drow);
  assert(n >= col+dcol);
#endif
  int iErr;
  int chunk_size = CHUNK_SIZE / n * n; 
  chunk_size = min(chunk_size, n*n);
  double* chunk = new double[chunk_size];

  int rows_in_chunk = chunk_size/n;
  int iRow=row;
  do {

    if(iRow+rows_in_chunk > drow)
      rows_in_chunk = drow-iRow;

    iErr=MPI_Allreduce(&M[iRow][0], 
		       chunk, rows_in_chunk*n, 
		       MPI_DOUBLE, MPI_SUM, comm);
    assert(iErr==MPI_SUCCESS);

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



#ifdef STOCH_TESTING
void sLinsysRoot::dumpMatrix(int scen, int proc, const char* nameToken, DenseSymMatrix& M) 
{
  int n = M.size();
  char szNumber[30];
  string strBuffer="";

  //assert(false);

  int iter = g_iterNumber;

  if(iter!=1 && iter!=5 && iter!=15 && iter!=25 && iter!=35 && iter!=45) return;


  char szFilename[256];
  if(scen==-1)
    sprintf(szFilename, "%s_%d__%d.mat", nameToken, n, iter);
  else 
    sprintf(szFilename, "%s_%03d_%d__%d.mat", nameToken, scen+1, n, iter);
  FILE* file = fopen(szFilename, "w");
  assert(file);
  

  for(int j=0; j<n; j++) {
    for(int i=0; i<n; i++) {
      sprintf(szNumber, "%22.16f ", M[i][j]);
      strBuffer += szNumber;
    }
    strBuffer += "\n";
    
    if(strBuffer.length()>1250000) {
      fwrite(strBuffer.c_str(), 1, strBuffer.length(), file);
      strBuffer = "";
    }
  }

  if(strBuffer.length()>0) {
    fwrite(strBuffer.c_str(), 1, strBuffer.length(), file);
  }
  
  fclose(file);
}

void sLinsysRoot::dumpRhs(int proc, const char* nameToken,  SimpleVector& rhs) 
{
  int n = rhs.length();
  char szNumber[30];
  string strBuffer="";


  int iter = g_iterNumber;
  if(iter!=0 && iter!=2 && iter!=20 && iter!=25 && iter!=55) return;

  char ipmPhase[4];
  if(g_iterNumber-iter>0) strcpy(ipmPhase, "co");
  else strcpy(ipmPhase, "pr");

  char szFilename[256];
  sprintf(szFilename, "%s_%s_%d__%d.mat", nameToken,ipmPhase, n, iter);


  for(int i=0; i<n; i++) {
    sprintf(szNumber, "%22.16f ", rhs[i]);
    strBuffer += szNumber;
  }

  FILE* file = fopen(szFilename, "w");
  assert(file);

  fwrite(strBuffer.c_str(), 1, strBuffer.length(), file);
  fclose(file);
}

#endif
