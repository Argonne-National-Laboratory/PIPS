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

#include "pipsport.h"

/*********************************************************************/
/************************** ROOT *************************************/
/*********************************************************************/

#ifdef STOCH_TESTING
double g_scenNum;
#endif

extern double g_iterNumber;
extern bool ipStartFound;
extern int gOuterSolve;

sLinsysRoot::sLinsysRoot(sFactory * factory_, sData * prob_)
  : sLinsys(factory_, prob_), iAmDistrib(0), sparseKktBuffer(nullptr)
{
  assert(dd!=nullptr);
  xDiag = nullptr;
  zDiag = nullptr;
  zDiagLinkCons = nullptr;
  kktDist = nullptr;

#ifdef TIMING
  int myRank; MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
  if( myRank == 0 )
     std::cout << "Rank 0: create LinSys children ..." << std::endl;
#endif

  createChildren(prob_);

#ifdef TIMING
  if( myRank == 0 )
     std::cout << "Rank 0: children created" << std::endl;
#endif

  precondSC = SCsparsifier(-1.0, mpiComm);

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
      sol2 = res2 = res3 = res4 = res5 = nullptr;
    }
  } else {
    sol  = res  = resx = resy = resz = nullptr;
    sol2 = res2 = res3 = res4 = res5 = nullptr;
  }

#ifdef DIST_PRECOND
  usePrecondDist = true;
#else
  usePrecondDist = false;
#endif

  // use sparse KKT if link structure is present
  hasSparseKkt = prob_->exploitingLinkStructure();

  usePrecondDist = usePrecondDist && hasSparseKkt && iAmDistrib;
  MatrixEntryTriplet_mpi = MPI_DATATYPE_NULL;

  initProperChildrenRange();
}

sLinsysRoot::sLinsysRoot(sFactory* factory_,
			 sData* prob_,
			 OoqpVector* dd_, 
			 OoqpVector* dq_,
			 OoqpVector* nomegaInv_,
			 OoqpVector* rhs_)
  : sLinsys(factory_, prob_, dd_, dq_, nomegaInv_, rhs_), iAmDistrib(0), sparseKktBuffer(nullptr)
{
  xDiag = nullptr;
  zDiag = nullptr;
  zDiagLinkCons = nullptr;
  kktDist = nullptr;

  createChildren(prob_);

  precondSC = SCsparsifier(-1.0, mpiComm);

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
      sol2 = res2 = res3 = res4 = res5 = nullptr;
    }
  } else {
      sol  = res  = resx = resy = resz = nullptr;
      sol2 = res2 = res3 = res4 = res5 = nullptr;
  }

#ifdef DIST_PRECOND
  usePrecondDist = true;
#else
  usePrecondDist = false;
#endif

  // use sparse KKT if (enough) 2 links are present
  hasSparseKkt = prob_->exploitingLinkStructure();

  usePrecondDist = usePrecondDist && hasSparseKkt && iAmDistrib;
  MatrixEntryTriplet_mpi = MPI_DATATYPE_NULL;

  initProperChildrenRange();
}

sLinsysRoot::~sLinsysRoot()
{
  for(size_t c=0; c<children.size(); c++)
    delete children[c];

  delete kktDist;

  delete[] sparseKktBuffer;
}

//this variable is just reset in this file; children will default to the "safe" linear solver
extern int gLackOfAccuracy;

void sLinsysRoot::factor2(sData *prob, Variables *vars)
{
  initializeKKT(prob, vars);

  // First tell children to factorize. 
  for(size_t c = 0; c < children.size(); c++)
    children[c]->factor2(prob->children[c], vars);

  for(size_t c = 0; c < children.size(); c++) {
#ifdef STOCH_TESTING
    g_scenNum=c;
#endif
    if(children[c]->mpiComm == MPI_COMM_NULL)
      continue;

    children[c]->stochNode->resMon.recFactTmChildren_start();
    //---------------------------------------------
    addTermToSchurCompl(prob, c);
    //---------------------------------------------
    children[c]->stochNode->resMon.recFactTmChildren_stop();
  }

#ifdef TIMING
  MPI_Barrier(MPI_COMM_WORLD);
  stochNode->resMon.recReduceTmLocal_start();
#endif 

  reduceKKT(prob);

 #ifdef TIMING
  stochNode->resMon.recReduceTmLocal_stop();
#endif

  finalizeKKT(prob, vars);

  factorizeKKT(prob);

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
  sLinsys* child=nullptr;
  assert(dd!=nullptr);
  assert(dynamic_cast<StochVector*>(dd) !=nullptr);
  StochVector& ddst = dynamic_cast<StochVector&>(*dd);
  StochVector& dqst = dynamic_cast<StochVector&>(*dq);
  StochVector& nomegaInvst = dynamic_cast<StochVector&>(*nomegaInv);
  StochVector& rhsst = dynamic_cast<StochVector&>(*rhs);

  //get the communicator from one of the vectors
  this->mpiComm = ddst.mpiComm;
  this->iAmDistrib = ddst.iAmDistrib;
  for(size_t it=0; it<prob->children.size(); it++) {
      assert(ddst.children[it]!=nullptr); 
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

void sLinsysRoot::initProperChildrenRange()
{
   assert(children.size() > 0);

   int childStart = -1;
   int childEnd = -1;
   for( size_t it = 0; it < children.size(); it++ )
   {
      if( childEnd != -1 )
         assert(children[it]->isDummy());

      if( children[it]->isDummy() )
      {
         // end of range?
         if( childStart != -1 && childEnd == -1 )
            childEnd = int(it);

         continue;
      }

      // start of range?
      if( childStart == -1 )
         childStart = int(it);
   }

   assert(childStart >= 0);

   if( childEnd == -1 )
   {
      assert(!children[children.size() - 1]->isDummy());
      childEnd = int(children.size());
   }

    assert(childStart < childEnd && childEnd <= int(children.size()));

    childrenProperStart = childStart;
    childrenProperEnd = childEnd;
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
   if( hasSparseKkt )
   {
      SparseSymMatrix* kkts = dynamic_cast<SparseSymMatrix*>(kkt);
      kkts->symPutZeroes();
   }
   else
   {
      DenseSymMatrix* kktd = dynamic_cast<DenseSymMatrix*>(kkt);
      myAtPutZeros(kktd);
   }
}

void sLinsysRoot::reduceKKT()
{
   reduceKKT(nullptr);
}

void sLinsysRoot::reduceKKT(sData* prob)
{
   if( usePrecondDist )
      reduceKKTdist(prob);
   else if( hasSparseKkt )
      reduceKKTsparse();
   else
      reduceKKTdense();
}


// collects (reduces) dense global Schur complement
void sLinsysRoot::reduceKKTdense()
{
   DenseSymMatrix* const kktd = dynamic_cast<DenseSymMatrix*>(kkt);

   // parallel communication
   if( iAmDistrib )
   {
      if( locnx > 0 )
         submatrixAllReduceDiagLower(kktd, 0, locnx, mpiComm);

      if( locmyl > 0 || locmzl > 0 )
      {
         const int locNxMy = locnx + locmy;
         assert(kktd->size() == locnx + locmy + locmyl + locmzl);

         // reduce lower left part
         if( locnx > 0 )
         {
            submatrixAllReduceFull(kktd, locNxMy, 0, locmyl + locmzl, locnx, mpiComm);
         }

         // reduce lower diagonal linking part
         submatrixAllReduceDiagLower(kktd, locNxMy, locmyl + locmzl, mpiComm);
      }
  }
}


// collects sparse global Schur complement
void sLinsysRoot::reduceKKTsparse()
{
   if( !iAmDistrib )
      return;

   int myRank; MPI_Comm_rank(mpiComm, &myRank);

   assert(kkt);

   SparseSymMatrix& kkts = dynamic_cast<SparseSymMatrix&>(*kkt);

   int* const krowKkt = kkts.krowM();
   double* const MKkt = kkts.M();
   const int sizeKkt = locnx + locmy + locmyl + locmzl;
   const int nnzKkt = krowKkt[sizeKkt];

   assert(kkts.size() == sizeKkt);
   assert(!kkts.isLower);

   reduceToProc0(nnzKkt, MKkt);
}

#define CHUNK_SIZE (1024*1024*64) //doubles = 128 MBytes (maximum)
void sLinsysRoot::reduceToProc0(int size, double* values)
{
   assert(values && values != sparseKktBuffer);
   assert(size > 0);

   int myRank; MPI_Comm_rank(mpiComm, &myRank);

   if( myRank == 0 && sparseKktBuffer == nullptr )
      sparseKktBuffer = new double[CHUNK_SIZE];

   const int reps = size / CHUNK_SIZE;
   const int res = size - CHUNK_SIZE * reps;
   assert(res >= 0 && res < CHUNK_SIZE);

   for( int i = 0; i < reps; i++ )
   {
      double* const start = &values[i * CHUNK_SIZE];
      MPI_Reduce(start, sparseKktBuffer, CHUNK_SIZE, MPI_DOUBLE, MPI_SUM, 0, mpiComm);

      if( myRank == 0 )
         memcpy(start, sparseKktBuffer, size_t(CHUNK_SIZE) * sizeof(double));
   }

   if( res > 0 )
   {
      double* const start = &values[reps * CHUNK_SIZE];
      MPI_Reduce(start, sparseKktBuffer, res, MPI_DOUBLE, MPI_SUM, 0, mpiComm);

      if( myRank == 0 )
         memcpy(start, sparseKktBuffer, size_t(res) * sizeof(double));
   }

#if 0
   if( myRank == 0 && sparseKktBuffer == nullptr )
      sparseKktBuffer = new double[nnzKkt];

   MPI_Reduce(MKkt, sparseKktBuffer, nnzKkt, MPI_DOUBLE, MPI_SUM, 0, mpiComm);

   if( myRank == 0 )
      memcpy(MKkt, sparseKktBuffer, size_t(nnzKkt) * sizeof(double));
#endif
}


void sLinsysRoot::registerMatrixEntryTripletMPI()
{
   assert(MatrixEntryTriplet_mpi == MPI_DATATYPE_NULL);

   const int nitems = 3;
   int blocklengths[3] = { 1, 1, 1 };
   MPI_Datatype Types[3] = { MPI_DOUBLE, MPI_INT, MPI_INT };
   MPI_Aint offsets[3];

   offsets[0] = offsetof(MatrixEntryTriplet, val);
   offsets[1] = offsetof(MatrixEntryTriplet, row);
   offsets[2] = offsetof(MatrixEntryTriplet, col);

   MPI_Type_create_struct(nitems, blocklengths, offsets, Types, &MatrixEntryTriplet_mpi);
   MPI_Type_commit(&MatrixEntryTriplet_mpi);
}

void sLinsysRoot::syncKKTdistLocalEntries(sData* prob)
{
   if( !iAmDistrib )
      return;

   assert(kkt && hasSparseKkt);

   SparseSymMatrix& kkts = dynamic_cast<SparseSymMatrix&>(*kkt);

   int* const krowKkt = kkts.krowM();
   int* const jColKkt = kkts.jcolM();
   double* const MKkt = kkts.M();

   const int childStart = childrenProperStart;
   const int childEnd = childrenProperEnd;
   int myRank; MPI_Comm_rank(mpiComm, &myRank);
   int size; MPI_Comm_size(mpiComm, &size);

   assert(size > 1);

   // MPI matrix entries triplet not registered yet?
   if( MatrixEntryTriplet_mpi == MPI_DATATYPE_NULL )
      registerMatrixEntryTripletMPI();

   // pack the entries that will be send below
   std::vector<MatrixEntryTriplet> prevEntries = this->packKKTdistOutOfRangeEntries(prob, childStart, childEnd);
   std::vector<MatrixEntryTriplet> myEntries(0);

   assert(prevEntries.size() > 0 && prevEntries[0].row == - 1 && prevEntries[0].col == - 1);

   // odd processes send first (one process back)
   if( myRank % 2 != 0 )
   {
      this->sendKKTdistLocalEntries(prevEntries);
   }

   // even processes (except last) receive first
   if( myRank % 2 == 0 && myRank != size - 1 )
   {
      assert(myEntries.size() == 0);
      myEntries = this->receiveKKTdistLocalEntries();
   }

   // even processes (except first) send
   if( myRank % 2 == 0 && myRank > 0 )
   {
      this->sendKKTdistLocalEntries(prevEntries);
   }

   // odd processes (except last) receive
   if( myRank % 2 != 0 && myRank != size - 1 )
   {
      assert(myEntries.size() == 0);
      myEntries = this->receiveKKTdistLocalEntries();
   }

   assert(myEntries.size() > 0 || myRank == size - 1 );

   int lastRow = 0;
   int lastC = -1;

#ifndef NDEBUG
   const std::vector<bool>& rowIsLocal = prob->getSCrowMarkerLocal();
   const std::vector<bool>& rowIsMyLocal = prob->getSCrowMarkerMyLocal();
#endif

   // finally, put received data into Schur complement matrix
   for( size_t i = 1; i < myEntries.size(); i++ )
   {
      const double val = myEntries[i].val;
      const int row = myEntries[i].row;
      const int col = myEntries[i].col;

      assert(myRank != size - 1);
      assert(row >= 0 && row < locnx + locmy + locmyl + locmzl);
      assert(col >= row && col < locnx + locmy + locmyl + locmzl);
      assert(val == val); // catch NaNs

      assert(rowIsMyLocal[row] || (rowIsMyLocal[col] && !rowIsLocal[row]));

      int c;

      // continue from last position?
      if( row == lastRow )
      {
         for( c = lastC + 1; c < krowKkt[row + 1]; c++ )
         {
            const int colKkt = jColKkt[c];

            if( colKkt == col )
            {
               MKkt[c] += val;
               break;
            }
         }

         // found the correct entry in last row?
         if( c != krowKkt[row + 1] )
         {
            assert(c < krowKkt[row + 1]);
            lastRow = row;
            lastC = c;
            continue;
         }
      }

      c = krowKkt[row];
      assert(col >= jColKkt[c]);

      for( ; c < krowKkt[row + 1]; c++ )
      {
         const int colKkt = jColKkt[c];

         if( colKkt == col )
         {
            MKkt[c] += val;
            break;
         }
      }

      assert(c != krowKkt[row + 1]);

      lastRow = row;
      lastC = c;
   }
}


std::vector<sLinsysRoot::MatrixEntryTriplet> sLinsysRoot::receiveKKTdistLocalEntries() const
{
   assert(kkt && hasSparseKkt);
   assert(MatrixEntryTriplet_mpi != MPI_DATATYPE_NULL);

   int myRank; MPI_Comm_rank(mpiComm, &myRank);
   int size; MPI_Comm_size(mpiComm, &size);
   const int nextRank = myRank + 1;
   assert(nextRank < size);

   // receive data from next process

   MPI_Status status;
   MPI_Probe(nextRank, 0, mpiComm, &status);

   int nEntries;
   MPI_Get_count(&status, MatrixEntryTriplet_mpi, &nEntries);

   assert(nEntries >= 1);

   std::vector<MatrixEntryTriplet> entries(nEntries);

   MPI_Recv((void*) &entries[0], nEntries, MatrixEntryTriplet_mpi, nextRank, 0, mpiComm, MPI_STATUS_IGNORE);

   assert(entries[0].row == -1 && entries[0].col == -1); // dummy check

   PIPSdebugMessage("myRank=%d received %d \n", myRank, nEntries);

   return entries;
}


void sLinsysRoot::sendKKTdistLocalEntries(const std::vector<MatrixEntryTriplet>& prevEntries) const
{
   int myRank; MPI_Comm_rank(mpiComm, &myRank);
   const int prevRank = myRank - 1;
   const int nEntries = int(prevEntries.size());

   assert(myRank >= 0);
   assert(nEntries > 0);
   assert(MatrixEntryTriplet_mpi != MPI_DATATYPE_NULL);

   PIPSdebugMessage("myRank=%d sends %d \n", myRank, nEntries);
   MPI_Send(&prevEntries[0], nEntries, MatrixEntryTriplet_mpi, prevRank, 0, mpiComm);
}

std::vector<sLinsysRoot::MatrixEntryTriplet> sLinsysRoot::packKKTdistOutOfRangeEntries(sData* prob, int childStart, int childEnd) const
{
   assert(kkt && hasSparseKkt);

   int myRank; MPI_Comm_rank(mpiComm, &myRank);

   SparseSymMatrix& kkts = dynamic_cast<SparseSymMatrix&>(*kkt);
   const std::vector<bool>& rowIsLocal = prob->getSCrowMarkerLocal();
   const std::vector<bool>& rowIsMyLocal = prob->getSCrowMarkerMyLocal();
   int* const krowKkt = kkts.krowM();
   int* const jColKkt = kkts.jcolM();
   double* const MKkt = kkts.M();
   const int sizeKkt = locnx + locmy + locmyl + locmzl;

   std::vector<MatrixEntryTriplet> packedEntries(0);

   // add dummy value
   packedEntries.push_back({-1.0, -1, -1});

   if( childStart > 0 )
   {
      assert(myRank > 0);

      // pack data
      for( int r = 0; r < sizeKkt; r++ )
      {
         const bool rIsLocal = rowIsLocal[r];

         if( rIsLocal && !rowIsMyLocal[r] )
         {
            for( int c = krowKkt[r]; c < krowKkt[r + 1]; c++ )
            {
               const int col = jColKkt[c];

               if( !rowIsMyLocal[col] )
               {
                  const double val = MKkt[c];

                  if( PIPSisZero(val) )
                     continue;

                  packedEntries.push_back({val, r, col});
               }
            }
         }

         if( rIsLocal )
            continue;

         for( int c = krowKkt[r]; c < krowKkt[r + 1]; c++ )
         {
            const int col = jColKkt[c];

            if( rowIsLocal[col] && !rowIsMyLocal[col] )
            {
               const double val = MKkt[c];

               if( PIPSisZero(val) )
                  continue;

               packedEntries.push_back({val, r, col});
            }
         }
      }
   }

   return packedEntries;
}


void sLinsysRoot::reduceKKTdist(sData* prob)
{
   assert(prob);
   assert(iAmDistrib);
   assert(kkt);

   const std::vector<bool>& rowIsLocal = prob->getSCrowMarkerLocal();
   const std::vector<bool>& rowIsMyLocal = prob->getSCrowMarkerMyLocal();

   SparseSymMatrix& kkts = dynamic_cast<SparseSymMatrix&>(*kkt);

   int* const krowKkt = kkts.krowM();
   int* const jColKkt = kkts.jcolM();
   double* const MKkt = kkts.M();
   const int sizeKkt = locnx + locmy + locmyl + locmzl;
   int nnzDistMyLocal = 0;
   int nnzDistShared = 0;
   int nnzDistLocal;

   std::vector<int> rowSizeMyLocal(sizeKkt, 0);
   std::vector<int> rowSizeShared(sizeKkt, 0);
   std::vector<int> rowSizeLocal(sizeKkt, 0);
   std::vector<int> rowIndexMyLocal(0);
   std::vector<int> colIndexMyLocal(0);

   assert(int(rowIsLocal.size()) == sizeKkt);

   // add up locally owned entries
   this->syncKKTdistLocalEntries(prob);

   // add B_0, F_0, G_0 and diagonals (all scattered)
   this->finalizeKKTdist(prob);

   precondSC.unmarkDominatedSCdistLocals(*prob, kkts);

   // compute row lengths
   for( int r = 0; r < sizeKkt; r++ )
   {
      if( rowIsMyLocal[r] )
      {
         for( int c = krowKkt[r]; c < krowKkt[r + 1]; c++ )
         {
            const int col = jColKkt[c];

            if( col < 0 )
               continue;

            nnzDistMyLocal++;
            rowSizeMyLocal[r]++;
            rowIndexMyLocal.push_back(r);
            colIndexMyLocal.push_back(col);
         }

         continue;
      }

      const bool rIsLocal = rowIsLocal[r];

      for( int c = krowKkt[r]; c < krowKkt[r + 1]; c++ )
      {
         const int col = jColKkt[c];

         if( col < 0 )
         {
            assert(-col - 1 >= 0 && -col - 1 < sizeKkt);
            assert((!(!rIsLocal && !rowIsLocal[-col - 1])));
            continue;
         }

         // is (r, col) a shared entry?
         if( !rIsLocal && !rowIsLocal[col] )
         {
            nnzDistShared++;
            rowSizeShared[r]++;
            assert(!rowIsMyLocal[col]);
         }

         if( rowIsMyLocal[col] )
         {
            nnzDistMyLocal++;
            rowSizeMyLocal[r]++;
            rowIndexMyLocal.push_back(r);
            colIndexMyLocal.push_back(col);
         }
      }
   }

   assert(int(rowIndexMyLocal.size()) == nnzDistMyLocal);

   // sum up local sizes
   MPI_Allreduce(&nnzDistMyLocal, &nnzDistLocal, 1, MPI_INT, MPI_SUM, mpiComm);
   MPI_Allreduce(&rowSizeMyLocal[0], &rowSizeLocal[0], sizeKkt, MPI_INT, MPI_SUM, mpiComm);

#ifndef NDEBUG
   {
      int nnzDistSharedMax;
      std::vector<int> rowSizeSharedMax(sizeKkt, 0);

      MPI_Allreduce(&nnzDistShared, &nnzDistSharedMax, 1, MPI_INT, MPI_MAX, mpiComm);
      MPI_Allreduce(&rowSizeShared[0], &rowSizeSharedMax[0], sizeKkt, MPI_INT, MPI_MAX, mpiComm);

      assert(nnzDistSharedMax == nnzDistShared);
      for( int i = 0; i < sizeKkt; i++ )
         assert(rowSizeShared[i] == rowSizeSharedMax[i]);
   }
#endif

   int localGatheredMyStart;
   int localGatheredMyEnd;

   std::vector<int> rowIndexGathered = PIPSallgathervInt(rowIndexMyLocal, mpiComm);
   std::vector<int> colIndexGathered = PIPSallgathervInt(colIndexMyLocal, mpiComm, localGatheredMyStart, localGatheredMyEnd);

#ifndef NDEBUG
   assert(int(rowIndexGathered.size()) == nnzDistLocal);
   assert(int(colIndexGathered.size()) == nnzDistLocal);
   assert(localGatheredMyEnd - localGatheredMyStart == nnzDistMyLocal);

   for( int i = 0; i < int(rowIndexMyLocal.size()); i++ )
   {
      assert(rowIndexMyLocal[i] == rowIndexGathered[i + localGatheredMyStart]);
      assert(colIndexMyLocal[i] == colIndexGathered[i + localGatheredMyStart]);
   }
#endif

   const int nnzDist = nnzDistLocal + nnzDistShared;

   assert(!kktDist || !kktDist->isLower);

   delete kktDist;
   kktDist = new SparseSymMatrix(sizeKkt, nnzDist, false);

   int* const krowDist = kktDist->krowM();
   int* const jColDist  = kktDist->jcolM();
   double* const MDist = kktDist->M();

   assert(krowDist[0] == 0);
   assert(sizeKkt > 0 && krowDist[1] == 0);

   memset(MDist, 0, nnzDist * sizeof(double));

   for( int r = 1; r < sizeKkt; r++ )
      krowDist[r + 1] = krowDist[r] + rowSizeLocal[r - 1] + rowSizeShared[r - 1];

   // fill in global and locally owned positions and values
   for( int r = 0; r < sizeKkt; r++ )
   {
      if( rowIsMyLocal[r] )
      {
         for( int c = krowKkt[r]; c < krowKkt[r + 1]; c++ )
         {
            assert(krowDist[r + 1] < nnzDist);

            const int col = jColKkt[c];

            if( col < 0 )
               continue;

            const double val = MKkt[c];

            MDist[krowDist[r + 1]] = val;
            jColDist[krowDist[r + 1]++] = col;
         }

         continue;
      }

      for( int c = krowKkt[r]; c < krowKkt[r + 1]; c++ )
      {
         const int col = jColKkt[c];

         if( col < 0 )
         {
            assert(-col - 1 >= 0 && -col - 1 < sizeKkt);
            assert(!(!rowIsLocal[r] && !rowIsLocal[-col - 1]));
            continue;
         }

         // is (r, col) a shared entry or locally owned?
         if( (!rowIsLocal[r] && !rowIsLocal[col]) || rowIsMyLocal[col] )
         {
            assert(krowDist[r + 1] < nnzDist);

            const double val = MKkt[c];

            MDist[krowDist[r + 1]] = val;
            jColDist[krowDist[r + 1]++] = col;
         }
      }
   }

   precondSC.resetSCdistEntries(kkts);

   // fill in gathered local pairs not inserted yet
   for( int i = 0; i < nnzDistLocal; i++ )
   {
      const int row = rowIndexGathered[i];
      const int col = colIndexGathered[i];

      assert(row >= 0 && row < sizeKkt);
      assert(col >= row && col < sizeKkt);

      // pair already added?
      if( i >= localGatheredMyStart && i < localGatheredMyEnd )
         continue;

      assert(krowDist[row + 1] < nnzDist);
      assert(MDist[krowDist[row + 1]] == 0.0);
      jColDist[krowDist[row + 1]++] = col;
   }

#ifndef NDEBUG
   assert(krowDist[0] == 0);
   assert(krowDist[sizeKkt] == nnzDist);

   for( int r = 0; r < sizeKkt; r++ )
   {
      assert(krowDist[r + 1] == krowDist[r] + rowSizeLocal[r] + rowSizeShared[r]);
      assert(krowDist[r + 1] >= krowDist[r]);
   }
#endif

   kktDist->getStorageRef().sortCols();

   assert(kktDist->getStorageRef().isValid());

   reduceToProc0(nnzDist, MDist);

   assert(kktDist->getStorageRef().isValid());
   assert(kktDist->getStorageRef().isSorted());
}

void sLinsysRoot::factorizeKKT()
{
   factorizeKKT(nullptr);
}

void sLinsysRoot::factorizeKKT(sData* prob)
{
  //stochNode->resMon.recFactTmLocal_start();  
#ifdef TIMING
  MPI_Barrier(mpiComm);
  extern double g_iterNumber;
  double st=MPI_Wtime();
#endif

  if( usePrecondDist )
  {
     int myRank; MPI_Comm_rank(mpiComm, &myRank);

     assert(kktDist);
     assert(prob);

     if( myRank == 0)
        precondSC.getSparsifiedSC_fortran(*prob, *kktDist);

#if 0
      {
         ofstream myfile;
         int mype; MPI_Comm_rank(mpiComm, &mype);

         if( mype == 0 )
         {
            printf("\n\n ...WRITE OUT kktDist! \n\n");
            myfile.open("../ADist.txt");
            int* ia = kktDist->krowM(); int* ja = kktDist->jcolM(); double* a = kktDist->M();

            for( int i = 0; i < kktDist->size(); i++ )
               for( int k = ia[i]; k < ia[i + 1]; k++ )
                  myfile << i << '\t' << ja[k - 1] << '\t' << a[k - 1] << endl;

            myfile.close();
         }

         MPI_Barrier(mpiComm);
         printf("...exiting (root) \n");
         exit(1);
      }
#endif

     solver->matrixRebuild(*kktDist);
  }
  else
  {
     // in solver allocate memory once and only reallocate if more memory needed?
     solver->matrixChanged();
  }

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
  assert( row >= 0 && row + rowExtent <= mat->size() );
  assert( col >= 0 && col + colExtent <= mat->size() );

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

void sLinsysRoot::addTermToSchurCompl(sData* prob, size_t childindex)
{
   // todo bad hack, should be removed once user parameters are available (along all global variables)
   ipIterations = ipStartFound ? static_cast<int>(g_iterNumber) : -1;

   assert(childindex < prob->children.size());
#ifdef PARDISO_BLOCKSC
   children[childindex]->addTermToSchurComplBlocked(prob->children[childindex], hasSparseKkt, *kkt);
#else
   if( hasSparseKkt )
   {
      SparseSymMatrix& kkts = dynamic_cast<SparseSymMatrix&>(*kkt);
      children[childindex]->addTermToSparseSchurCompl(prob->children[childindex], kkts);
   }
   else
   {
      DenseSymMatrix& kktd = dynamic_cast<DenseSymMatrix&>(*kkt);
      children[childindex]->addTermToDenseSchurCompl(prob->children[childindex], kktd);
   }
#endif
}

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
  assert(n >= endCol);

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
#ifndef NDEBUG
    const int iErr=MPI_Allreduce(&M[iRow][0], chunk, rows_in_chunk*n, MPI_DOUBLE, MPI_SUM, comm);
    assert(iErr==MPI_SUCCESS);
#else
    MPI_Allreduce(&M[iRow][0], chunk, rows_in_chunk*n, MPI_DOUBLE, MPI_SUM, comm);
#endif

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


void sLinsysRoot::submatrixAllReduceFull(DenseSymMatrix* A,
                   int startRow, int startCol, int nRows, int nCols,
                 MPI_Comm comm)
{
   double** const M = A->mStorage->M;

   assert(nRows > 0);
   assert(nCols > 0);
   assert(startRow >= 0);
   assert(startCol >= 0);

   const int endRow = startRow + nRows;

   assert(A->mStorage->n >= endRow);
   assert(A->mStorage->n >= startCol + nCols);

   const int buffersize = nRows * nCols;

   double* const bufferSend = new double[buffersize];
   double* const bufferRecv = new double[buffersize];

   // copy into send buffer
   int counter = 0;
   const size_t nColBytes = nCols * sizeof(double);

   for( int r = startRow; r < endRow; r++ )
   {
      memcpy(&bufferSend[counter], &M[r][startCol], nColBytes);
      counter += nCols;
   }

   assert(counter == buffersize);

#ifndef NDEBUG
   const int iErr = MPI_Allreduce(bufferSend, bufferRecv, buffersize, MPI_DOUBLE, MPI_SUM, comm);
   assert(iErr == MPI_SUCCESS);
#else
   MPI_Allreduce(bufferSend, bufferRecv, buffersize, MPI_DOUBLE, MPI_SUM, comm);
#endif

   // copy back
   counter = 0;
   for( int r = startRow; r < endRow; r++ )
   {
      memcpy(&M[r][startCol], &bufferRecv[counter], nColBytes);
      counter += nCols;
   }

   assert(counter == buffersize);

   delete[] bufferRecv;
   delete[] bufferSend;
}


void sLinsysRoot::submatrixAllReduceDiagLower(DenseSymMatrix* A,
                   int substart, int subsize,
                 MPI_Comm comm)
{
   double** const M = A->mStorage->M;

   assert(subsize >= 0);
   assert(substart >= 0);

   if( subsize == 0)
      return;

   const int subend = substart + subsize;
   assert(A->mStorage->n >= subend);

   // number of elements in lower matrix triangle (including diagonal)
   const int buffersize = (subsize * subsize + subsize) / 2;
   assert(buffersize > 0);

   double* const bufferSend = new double[buffersize];
   double* const bufferRecv = new double[buffersize];

   int counter = 0;

   for( int i = substart; i < subend; i++ )
      for( int j = substart; j <= i; j++ )
      {
         assert(counter < buffersize);
         bufferSend[counter++] = M[i][j];
      }

   assert(counter == buffersize);

#ifndef NDEBUG
   const int iErr = MPI_Allreduce(bufferSend, bufferRecv, buffersize, MPI_DOUBLE, MPI_SUM, comm);
   assert(iErr == MPI_SUCCESS);
#else
   MPI_Allreduce(bufferSend, bufferRecv, buffersize, MPI_DOUBLE, MPI_SUM, comm);
#endif

   counter = 0;
   for( int i = substart; i < subend; i++ )
      for( int j = substart; j <= i; j++ )
      {
         assert(counter < buffersize);
         M[i][j] = bufferRecv[counter++];
      }

   delete[] bufferSend;
   delete[] bufferRecv;
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
