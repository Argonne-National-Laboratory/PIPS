/* PIPS
   Authors: Cosmin Petra and Miles Lubin
   See license and copyright information in the documentation */


#include "sLinsysRoot.h"
#include "sTree.h"
#include "sFactory.h"
#include "sData.h"
#include "sDummyLinsys.h"
#include "sLinsysLeaf.h"
#include "math.h"

/*********************************************************************/
/************************** ROOT *************************************/
/*********************************************************************/

#ifdef STOCH_TESTING
extern double g_iterNumber;
double g_scenNum;
#endif

extern int gOuterSolve;

sLinsysRoot::sLinsysRoot(sFactory * factory_, sData * prob_)
  : sLinsys(factory_, prob_), iAmDistrib(0)
{
  assert(dd!=NULL);
  createChildren(prob_);

  if(gOuterSolve) {
    // stuff for iterative refimenent and BiCG
    sol  = factory_->tree->newRhs();
    res  = factory_->tree->newRhs();
    resx = factory_->tree->newPrimalVector();
    resy = factory_->tree->newDualYVector();
    resz = factory_->tree->newDualZVector();
    if(gOuterSolve==2) {
      //BiCGStab; additional vectors needed
      sol2 = factory_->tree->newRhs();
      res2 = factory_->tree->newRhs();
      res3 = factory_->tree->newRhs();
      res4 = factory_->tree->newRhs();
      res5 = factory_->tree->newRhs();
    } else {
      sol2 = res2 = res3 = res4 = res5 = NULL;
    }
  } else {
    sol  = res  = resx = resy = resz = NULL;
    sol2 = res2 = res3 = res4 = res5 = NULL;
  }
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

  if(gOuterSolve) {
      // stuff for iterative refimenent and BiCG 
      sol  = factory_->tree->newRhs();
      res  = factory_->tree->newRhs();
      resx = factory_->tree->newPrimalVector();
      resy = factory_->tree->newDualYVector();
      resz = factory_->tree->newDualZVector();
    if(gOuterSolve==2) {
      //BiCGStab; additional vectors needed
      sol2 = factory_->tree->newRhs();
      res2 = factory_->tree->newRhs();
      res3 = factory_->tree->newRhs();
      res4 = factory_->tree->newRhs();
      res5 = factory_->tree->newRhs();
    } else {
      sol2 = res2 = res3 = res4 = res5 = NULL;
    }
  } else {
      sol  = res  = resx = resy = resz = NULL;
      sol2 = res2 = res3 = res4 = res5 = NULL;
  }
}

sLinsysRoot::~sLinsysRoot()
{
  for(size_t c=0; c<children.size(); c++)
    delete children[c];
}

//this variable is just reset in this file; children will default to the "safe" linear solver
extern int gLackOfAccuracy;

void sLinsysRoot::factor2(sData *prob, Variables *vars)
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

#ifdef STOCH_TESTING
  const double epsilon = 1e-5;
  for( int k = 0; k < locnx + locmy + locmyl + locmzl; k++)
       for( int k2 = 0; k2 < locnx + locmy + locmyl + locmzl; k2++)
       {
           if( fabs(kktd[k][k2] - kktd[k2][k]) > epsilon )
              std::cout << "SYMMETRY FAIL, > eps " << fabs(kktd[k][k2] - kktd[k2][k])  << std::endl;
           assert(fabs(kktd[k][k2] - kktd[k2][k]) <= epsilon);
       }
#endif

  factorizeKKT();

  //if (mype==0) dumpMatrix(-1, 0, "kkt", kktd);

#ifdef TIMING
  afterFactor();
#endif

  //gLackOfAccuracy=0;
}

#ifdef TIMING
void sLinsysRoot::afterFactor()
{
  int mype; MPI_Comm_rank(mpiComm, &mype);

  if( (mype/256)*256==mype) {
    for (size_t c=0; c<children.size(); c++) {
      if (children[c]->mpiComm == MPI_COMM_NULL) continue;
      
      printf("  rank %d NODE %4zu SPFACT %g BACKSOLVE %g SEC ITER %d\n", mype, c,
	     children[c]->stochNode->resMon.eFact.tmLocal,
	     children[c]->stochNode->resMon.eFact.tmChildren, (int)g_iterNumber);
    }
  }
  if( (mype/1024)*1024==mype) {
    for (size_t c=0; c<children.size(); c++) {
      if (children[c]->mpiComm == MPI_COMM_NULL) continue;
      
      double redall = stochNode->resMon.eReduce.tmLocal;
      double redscat = stochNode->resMon.eReduceScatter.tmLocal;
      printf("  rank %d REDUCE %g SEC ITER %d REDSCAT %g DIFF %g\n", mype, redall, 
	     (int)g_iterNumber, redscat, redall-redscat);
    }
  }
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
  assert(!b.vecl);

  //this code actually works on a single CPU too :)
  if (iAmDistrib) {
    //only one process add b0
    int myRank; MPI_Comm_rank(mpiComm, &myRank);
    if(myRank>0) {
      b0.setToZero();
    }
  }

  // compute B_i^T rhs_i and add it up

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
    children[it]->addLniziLinkCons(prob->children[it], b0, zi, locmy, locmz);

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
 
    double* buffer = new double[b0.length()];
    MPI_Allreduce(b0.elements(), buffer, b0.length(),
		  MPI_DOUBLE, MPI_SUM, mpiComm);

    b0.copyFromArray(buffer);

    delete[] buffer;
  }
#ifdef TIMING 
  stochNode->resMon.recReduceTmLocal_stop();
#endif
  //dumpRhs(0, "rhs",  b0);

#ifdef TIMING
  stochNode->resMon.eLsolve.clear();
  stochNode->resMon.recLsolveTmLocal_start();
#endif
  solver->Lsolve(b0);
#ifdef TIMING
  stochNode->resMon.recLsolveTmLocal_stop();
#endif

}


void sLinsysRoot::Ltsolve2( sData *prob, StochVector& x, SimpleVector& xp)
{
  StochVector& b   = dynamic_cast<StochVector&>(x);
  SimpleVector& bi = dynamic_cast<SimpleVector&>(*b.vec);

#ifdef TIMING
  stochNode->resMon.eLtsolve.clear();
  stochNode->resMon.recLtsolveTmLocal_start();
#endif
  //b_i -= Lni^T x0
  this->LniTransMult(prob, bi, -1.0, xp);
  solver->Ltsolve(bi);

#ifdef TIMING
  stochNode->resMon.recLtsolveTmLocal_stop();
#endif
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

#ifdef TIMING
  stochNode->resMon.eLtsolve.clear();
  stochNode->resMon.recLtsolveTmLocal_start();
#endif
  solver->Ltsolve(b0);
#ifdef TIMING
  stochNode->resMon.recLtsolveTmLocal_stop();
#endif
  //dumpRhs(0, "sol",  b0);

  SimpleVector& x0 = b0; //just another name, for clarity
  
  // Li^T\bi for each child i. The backsolve needs z0

  for(size_t it=0; it<children.size(); it++) {
    children[it]->Ltsolve2(prob->children[it], *b.children[it], x0);
  }
#ifdef TIMING
  int myRank; MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

  if(256*(myRank/256) == myRank) {
    double tTotResChildren=0.0;
    for(size_t it=0; it<children.size(); it++) {
	if (children[it]->mpiComm == MPI_COMM_NULL) continue;
	tTotResChildren += children[it]->stochNode->resMon.eLsolve.tmChildren;
	tTotResChildren += children[it]->stochNode->resMon.eLsolve.tmLocal;
    }
    double tComm=stochNode->resMon.eReduce.tmLocal;


    //double tTotChildren=0.0; 
    //for(size_t it=0; it<children.size(); it++) {   
    //  tTotChildren += children[it]->stochNode->resMon.eDsolve.tmChildren;
    //  tTotChildren += children[it]->stochNode->resMon.eDsolve.tmLocal;
    //} 
    double tStg1=stochNode->resMon.eDsolve.tmLocal;  

    double tTotStg2Children=0.0;
    for(size_t it=0; it<children.size(); it++) {
      if (children[it]->mpiComm == MPI_COMM_NULL) continue;
      tTotStg2Children += children[it]->stochNode->resMon.eLtsolve.tmChildren;
      tTotStg2Children += children[it]->stochNode->resMon.eLtsolve.tmLocal;
    }
    cout << "  rank " << myRank << " "
	 << "Resid comp " << tTotResChildren << " " << "reduce " << tComm << " "
	 << "1stStage solve " << tStg1 << " "
	 << "2ndStage solve " << tTotStg2Children <<endl;
  }
#endif


}

void sLinsysRoot::Dsolve( sData *prob, OoqpVector& x )
{
//#ifdef TIMING
//    double tTot = MPI_Wtime();
//#endif
  StochVector& b = dynamic_cast<StochVector&>(x);

  //! commented - already done in addLnizi - cpetra
  //  for(size_t it=0; it<children.size(); it++) {
  //  children[it]->Dsolve(prob->children[it], *b.children[it]);
  //}

  SimpleVector& b0 = dynamic_cast<SimpleVector&>(*b.vec);
#ifdef TIMING
  stochNode->resMon.eDsolve.clear();
  stochNode->resMon.recDsolveTmLocal_start();
#endif

  solveReducedLinkCons(prob, b0);

#ifdef TIMING
  stochNode->resMon.recDsolveTmLocal_stop();
#endif
}



void sLinsysRoot::createChildren(sData* prob)
{
  sLinsys* child=NULL;
  assert(dd!=NULL);
  assert(dynamic_cast<StochVector*>(dd) !=NULL);
  StochVector& ddst = dynamic_cast<StochVector&>(*dd);
  StochVector& dqst = dynamic_cast<StochVector&>(*dq);
  StochVector& nomegaInvst = dynamic_cast<StochVector&>(*nomegaInv);
  StochVector& rhsst = dynamic_cast<StochVector&>(*rhs);

  //get the communicator from one of the vectors
  this->mpiComm = ddst.mpiComm;
  this->iAmDistrib = ddst.iAmDistrib;
  for(size_t it=0; it<prob->children.size(); it++) {
      assert(ddst.children[it]!=NULL); 
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
  zDiagLinkCons = zdiag.vecl;

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
  if (iAmDistrib) {
    submatrixAllReduce(kktd, 0, 0, locnx, locnx, mpiComm);
	if( locmyl > 0 || locmzl > 0 )
	{
	   int locNxMy = locnx + locmy;
	   int locNxMyMylMzl = locnx + locmy + locmyl + locmzl;

	   assert(kktd->size() == locNxMyMylMzl);

	   // reduce upper right part
	   submatrixAllReduce(kktd, 0, locNxMy, locnx, locmyl + locmzl, mpiComm);

	   // preserve symmetry todo memopt!
	   double** M = kktd->mStorage->M;
	   for( int k = locNxMy; k < locNxMyMylMzl; k++ )
		  for( int k2 = 0; k2 < locnx; k2++ )
		    M[k][k2] = M[k2][k];

	   // reduce lower diagonal part
	   submatrixAllReduce(kktd, locNxMy, locNxMy, locmyl + locmzl, locmyl + locmzl, mpiComm);
	}
  }
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
  if( (mype/512)*512==mype )
    printf("  rank %d 1stSTAGE FACT %g SEC ITER %d\n", mype, st, (int)g_iterNumber);
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


#define CHUNK_SIZE 1024*1024*64 //doubles  = 128 MBytes (maximum)
void sLinsysRoot::submatrixAllReduce(DenseSymMatrix* A,
		             int startRow, int startCol, int nRows, int nCols,
				     MPI_Comm comm)
{
  double ** M = A->mStorage->M;
  int n = A->mStorage->n;

  assert(nRows > 0);
  assert(nCols > 0);
  assert(startRow >= 0);
  assert(startCol >= 0);

  int endRow = startRow + nRows;
  int endCol = startCol + nCols;

  assert(n >= endRow);
  if( n < endCol )
	  cout << n << " " <<  endCol << endl;
  assert(n >= endCol);

  int iErr;
  int chunk_size = (CHUNK_SIZE / n) * n;
  chunk_size = min(chunk_size, n*nRows);

  double* chunk = new double[chunk_size];

  int rows_in_chunk = chunk_size/n;

  int iRow=startRow;

  // main loop
  do {

    if( iRow + rows_in_chunk > endRow )
      rows_in_chunk = endRow - iRow;

    assert(rows_in_chunk > 0);

    iErr=MPI_Allreduce(&M[iRow][0], chunk, rows_in_chunk*n, MPI_DOUBLE, MPI_SUM, comm);

    assert(iErr==MPI_SUCCESS);

    int shift = 0;

    // copy into M
    for( int i = iRow; i < iRow + rows_in_chunk; i++ ) {
      for( int j = startCol; j < endCol; j++ )
	    M[i][j] = chunk[shift+j];

      // shift one row forward
      shift += n;
    }
    iRow += rows_in_chunk;

  } while( iRow < endRow );

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
