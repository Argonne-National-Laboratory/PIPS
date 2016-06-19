/* PIPS
   Authors: Cosmin Petra and Miles Lubin
   See license and copyright information in the documentation */

/* 2015. Modified by Nai-Yuan Chiang for NLP*/

#include "sLinsysRoot.h"
#include "sTree.h"
#include "sFactory.h"
#include "sData.h"
#include "sDummyLinsys.h"
#include "sLinsysLeaf.h"
#include "OoqpVector.h"

#include "RegularizationAlg.h"

using namespace std;
/*********************************************************************/
/************************** ROOT *************************************/
/*********************************************************************/

#ifdef STOCH_TESTING
extern double g_iterNumber;
double g_scenNum;
#endif

extern int gOuterSolve;
extern int separateHandDiag;

extern int gDoIR_Aug;
extern int gDoIR_Full;

extern int gMaxIR;
extern double gIRtol;

extern int gUsePetsc;
extern int gUser_Defined_PC;
extern int gUsePetscOuter;
extern int gisNLP;

#ifdef TIMING
double t_start, troot_total, taux, tchild_total, tcomm_total;
#endif


sLinsysRoot::sLinsysRoot(sFactory * factory_, sData * prob_, bool createChild)
  : sLinsys(factory_, prob_), iAmDistrib(0)
{
  assert(dd!=NULL);
  if( createChild ) 
  	createChildren(prob_);

  sol2Bicg = res2Bicg = res3Bicg = res4Bicg = res5Bicg = NULL;
  sol2 = res2 = res3 = res4 = res5 = NULL;

  if(gOuterSolve<3) { 
    // stuff for iterative refimenent and BiCG
    sol  = factory_->tree->newRhs();
    res  = factory_->tree->newRhs();
    resx = factory_->tree->newPrimalVector();
    resy = factory_->tree->newDualYVector();
    resz = factory_->tree->newDualZVector();
    if(gOuterSolve==2) {
      //BiCGStab for compressed; additional vectors needed
      sol2 = factory_->tree->newRhs();
      res2 = factory_->tree->newRhs();
      res3 = factory_->tree->newRhs();
      res4 = factory_->tree->newRhs();
      res5 = factory_->tree->newRhs();
    }
  }else if(gOuterSolve>=3) {
	sol  = factory_->tree->newRhsXSYZ();
	res  = factory_->tree->newRhsXSYZ();
	resx = factory_->tree->newPrimalVector();
	ress = factory_->tree->newDualZVector();
	resy = factory_->tree->newDualYVector();
	resz = factory_->tree->newDualZVector();
    if(gDoIR_Aug ==1 || gUsePetscOuter!=0){
	  sol2 = factory_->tree->newRhsXSYZ();
	  res2 = factory_->tree->newRhsXSYZ();
	  res3 = res4 = res5 = NULL;
    }
    if(gOuterSolve>=4) {
      //BiCGStab; additional vectors needed
      sol2Bicg  = factory_->tree->newRhsXSYZ();
      res2Bicg  = factory_->tree->newRhsXSYZ();
      res3Bicg  = factory_->tree->newRhsXSYZ();
      res4Bicg  = factory_->tree->newRhsXSYZ();
      res5Bicg  = factory_->tree->newRhsXSYZ();
    }	
  }else {
    sol  = res  = resx = resy = resz = NULL;
  }
  firstBUpdate = true; 
  firstQUpdate = true;
}

sLinsysRoot::sLinsysRoot(sFactory* factory_,
			 sData* prob_,
			 OoqpVector* dd_, 
			 OoqpVector* dq_,
			 OoqpVector* nomegaInv_,
			 OoqpVector* rhs_,
			 OoqpVector* additiveDiag_, bool createChild)
  : sLinsys(factory_, prob_, dd_, dq_, nomegaInv_, rhs_, additiveDiag_), iAmDistrib(0)
{
  if( createChild ) 
	createChildren(prob_);

  if(gOuterSolve) {
      sol2 = res2 = res3 = res4 = res5 = NULL;  	
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
    } 
  } else {
      sol  = res  = resx = resy = resz = NULL;
      sol2 = res2 = res3 = res4 = res5 = NULL;
  }
  firstBUpdate = true; 
  firstQUpdate = true;
}

sLinsysRoot::~sLinsysRoot()
{
//  for(size_t c=0; c<children.size(); c++)
//    delete children[c];

  for(size_t it=0; it<children.size(); it++) {
    children[it]->deleteChildren();
    delete children[it];
  }
  children.clear();

	if(dd) { delete dd; dd=NULL;} 
	if(dq) { delete dq; dq=NULL;} 
	if(rhs) { delete rhs; rhs=NULL;} 
	if(nomegaInv) { delete nomegaInv; nomegaInv=NULL;} 
	if(temp_diagX) { delete temp_diagX; temp_diagX=NULL;} 
	if(temp_diagZ) { delete temp_diagZ; temp_diagZ=NULL;} 
	if(temp_diagS) { delete temp_diagS; temp_diagS=NULL;} 
	if(temp_diagY) { delete temp_diagY; temp_diagY=NULL;}   
	if(sol) delete sol; if(res) delete res; if(resx) delete resx; if(ress) delete ress; if(resy) delete resy; if(resz) delete resz; 
	if(sol2) delete sol2; if(res2) delete res2; 
	if(sol2Bicg) delete sol2Bicg;
	if(res2Bicg) delete res2Bicg; if(res3Bicg) delete res3Bicg; if(res4Bicg) delete res4Bicg; if(res5Bicg) delete res5Bicg;
  
}

//this variable is just reset in this file; children will default to the "safe" linear solver
extern int gLackOfAccuracy;

// FIXME_NY: we can skip build schur if one diag is singular
int sLinsysRoot::factor2(sData *prob, Variables *vars)
{
  int negEVal=0, tempNegEVal=0;
  int return_NegEval=-1;  
  int matIsSingular=0,matIsSingularAllReduce;
  int mype; MPI_Comm_rank(mpiComm, &mype);
  	
  DenseSymMatrix& kktd = dynamic_cast<DenseSymMatrix&>(*kkt);
  initializeKKT(prob, vars);

  // First tell children to factorize. 
  for(size_t c=0; c<children.size(); c++) {
    tempNegEVal = children[c]->factor2(prob->children[c], vars);
	if(tempNegEVal<0){
	  matIsSingular = 1; 
	}else{
	  negEVal += tempNegEVal;
	}
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
  
  //printf("(%d, %d) --- %f\n", PROW,PCOL, kktd[PROW][PCOL]);

  MPI_Allreduce(&matIsSingular, &matIsSingularAllReduce, 1, MPI_INT, MPI_SUM, mpiComm);

  if(0==matIsSingularAllReduce){
  	// all the diag mat is nonsingular
  	MPI_Allreduce(&negEVal, &return_NegEval, 1, MPI_INT, MPI_SUM, mpiComm);
	negEVal = factorizeKKT();
	if(negEVal<0){ 
	  return_NegEval = -1;
	}else{
	  return_NegEval += negEVal;
	}
  }
  
#ifdef TIMING
  afterFactor();
#endif

  return return_NegEval;

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
      double redall = stochNode->resMon.eReduce.tmLocal;
      double redscat = stochNode->resMon.eReduceScatter.tmLocal;
      printf("  rank %d REDUCE %g SEC ITER %d REDSCAT %g DIFF %g\n", mype, redall, 
	     (int)g_iterNumber, redscat, redall-redscat);
  }
}
#endif

void sLinsysRoot::Lsolve(sData *prob, OoqpVector& x)
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
  solveReduced(prob, b0);
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
	  	  if(gOuterSolve >= 3 && separateHandDiag==1){
			StochVector& additiveDiagst = dynamic_cast<StochVector&>(*additiveDiag);
			//child = new sLinsysLeaf(dynamic_cast<NlpGenStoch*>(factory),
	        child = stochFactory->newLinsysLeaf(prob->children[it],
						  ddst.children[it],
						  dqst.children[it],
						  nomegaInvst.children[it],
						  rhsst.children[it],
						  additiveDiagst.children[it]);
		  }else
	        child = stochFactory->newLinsysLeaf(prob->children[it],
						  ddst.children[it],
						  dqst.children[it],
						  nomegaInvst.children[it],
						  rhsst.children[it],NULL);		  	
	  } else {
	      if(gOuterSolve >= 3 && separateHandDiag==1){
		  	StochVector& additiveDiagst = dynamic_cast<StochVector&>(*additiveDiag);
	        //child = new sLinsysRoot(dynamic_cast<NlpGenStoch*>(factory), 
	        child = stochFactory->newLinsysRoot(prob->children[it],
						  ddst.children[it],
						  dqst.children[it],
						  nomegaInvst.children[it],
						  rhsst.children[it],
						  additiveDiagst.children[it]);
	      }
		  else
	        child = stochFactory->newLinsysRoot(prob->children[it],
						  ddst.children[it],
						  dqst.children[it],
						  nomegaInvst.children[it],
						  rhsst.children[it],NULL);			  
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

void sLinsysRoot::putSDiagonal( OoqpVector& sdiag_ )
{
  StochVector& sdiag = dynamic_cast<StochVector&>(sdiag_);
  assert(children.size() == sdiag.children.size());

  //kkt->atPutDiagonal( 0, *xdiag.vec );
  sDiag = sdiag.vec;
 
  // propagate it to the subtree
  for(size_t it=0; it<children.size(); it++)
    children[it]->putSDiagonal(*sdiag.children[it]);
}

void sLinsysRoot::putYDualDiagonal( OoqpVector& ydiag_ )
{
  StochVector& ydiag = dynamic_cast<StochVector&>(ydiag_);
  assert(children.size() == ydiag.children.size());

  yDiag = ydiag.vec;

  // propagate it to the subtree
  for(size_t it=0; it<children.size(); it++)
    children[it]->putYDualDiagonal(*ydiag.children[it]);
}


void sLinsysRoot::putZDiagonal( OoqpVector& zdiag_ )
{
  StochVector& zdiag = dynamic_cast<StochVector&>(zdiag_);
  assert(children.size() == zdiag.children.size());

  zDiag = zdiag.vec;

  // propagate it to the subtree
  for(size_t it=0; it<children.size(); it++)
    children[it]->putZDiagonal(*zdiag.children[it]);
}

void sLinsysRoot::setXDiagonal( OoqpVector& xdiag_ )
{
  StochVector& xdiag = dynamic_cast<StochVector&>(xdiag_);
  assert(children.size() == xdiag.children.size());

  //kkt->atPutDiagonal( 0, *xdiag.vec );
  xDiag = xdiag.vec;
 
  // propagate it to the subtree
  for(size_t it=0; it<children.size(); it++)
    children[it]->setXDiagonal(*xdiag.children[it]);
}

void sLinsysRoot::setSDiagonal( OoqpVector& sdiag_ )
{
  StochVector& sdiag = dynamic_cast<StochVector&>(sdiag_);
  assert(children.size() == sdiag.children.size());

  //kkt->atPutDiagonal( 0, *xdiag.vec );
  sDiag = sdiag.vec;
 
  // propagate it to the subtree
  for(size_t it=0; it<children.size(); it++)
    children[it]->setSDiagonal(*sdiag.children[it]);
}

void sLinsysRoot::setYDiagonal( OoqpVector& ydiag_ )
{
  StochVector& ydiag = dynamic_cast<StochVector&>(ydiag_);
  assert(children.size() == ydiag.children.size());

  yDiag = ydiag.vec;

  // propagate it to the subtree
  for(size_t it=0; it<children.size(); it++)
    children[it]->setYDiagonal(*ydiag.children[it]);
}


void sLinsysRoot::setZDiagonal( OoqpVector& zdiag_ )
{
  StochVector& zdiag = dynamic_cast<StochVector&>(zdiag_);
  assert(children.size() == zdiag.children.size());

  zDiag = zdiag.vec;

  // propagate it to the subtree
  for(size_t it=0; it<children.size(); it++)
    children[it]->setZDiagonal(*zdiag.children[it]);
}


void sLinsysRoot::setAdditiveDiagonal()
{
  assert(gOuterSolve >= 3 && separateHandDiag==1);

  // FIXME_NY: how to do it in root node?
  assert("not finished"&&0);
  
  StochVector& additiveDiag_ = dynamic_cast<StochVector&>(*additiveDiag);
  kkt->setAdditiveDiagonal(*additiveDiag_.vec);

  // propagate it to the subtree
  for(size_t it=0; it<children.size(); it++)
    children[it]->setAdditiveDiagonal();

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
//  if( nxupp + nxlow > 0 ) {
    //dd      = OoqpVectorHandle(stochNode->newPrimalVector());
    //dq      = OoqpVectorHandle(stochNode->newPrimalVector());
    dd = stochNode->newPrimalVector();
    dq = stochNode->newPrimalVector();
    data->getDiagonalOfQ( *dq );
//  }
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
  if(iAmDistrib)
  {
    int N;
    if(gOuterSolve>=3 )
    {
        int locns = locmz;
        N = locnx+locns+locmy+locmz;
    }
    else
      N = locnx+locmy+locmz;

    if( (stochNode->mle()+stochNode->mli())>0)
      submatrixAllReduce(kktd, 0, 0, N, N, mpiComm);
    else
      submatrixAllReduce(kktd, 0, 0, locnx, locnx, mpiComm);
  }
}

int sLinsysRoot::factorizeKKT()
{
  int negEValTemp=0;
  
  //stochNode->resMon.recFactTmLocal_start();  
#ifdef TIMING
  MPI_Barrier(mpiComm);
  extern double g_iterNumber;
  double st=MPI_Wtime();
#endif

  negEValTemp = solver->matrixChanged();

  //stochNode->resMon.recFactTmLocal_stop(); 
#ifdef TIMING
  st = MPI_Wtime()-st;
  MPI_Barrier(mpiComm);
  int mype; MPI_Comm_rank(mpiComm, &mype);
  // note, this will include noop scalapack processors
  if( (mype/256)*256==mype )
    printf("  rank %d 1stSTAGE FACT %g SEC ITER %d\n", mype, st, (int)g_iterNumber);
#endif

  return negEValTemp;
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

void
sLinsysRoot::UpdateMatrices( Data * prob_in, int const updateLevel)
{
  int useUpdate=updateLevel;
  if(!gisNLP) useUpdate=1;

  sData* prob = dynamic_cast<sData*>(prob_in);

  if(useUpdate>=2){
  	if(gOuterSolve < 3){	
	  kkt->symAtSetSubmatrix( 0, 0, prob->getLocalQ(), 0, 0, locnx, locnx, firstQUpdate, LocQMap);
	  if(locmy>0)
	  	kkt->symAtSetSubmatrix( locnx, 0, prob->getLocalB(), 0, 0, locmy, locnx, firstBUpdate, LocBMap);	
	}else{
	  kkt->symAtSetSubmatrix( 0, 0, prob->getLocalQ(), 0, 0, locnx, locnx, firstQUpdate, LocQMap);
	  if(locmy>0)
	  	kkt->symAtSetSubmatrix( locnx + locmz, 0, prob->getLocalB(), 0, 0, locmy, locnx, firstBUpdate, LocBMap);
	}
	firstQUpdate = false;
	firstBUpdate = false;	
  }
  
  // propagate it to the subtree
  for(size_t it=0; it<children.size(); it++)
    children[it]->UpdateMatrices(prob->children[it],updateLevel);

  if(useUpdate>=1){

	if(gOuterSolve < 3){
	  setXDiagonal( *dd );
	  setYDiagonal( *temp_diagY);	  
	  setZDiagonal( *nomegaInv );	  
	}else if (gOuterSolve >= 3 && separateHandDiag==1){
	  joinRHSXSYZ(*additiveDiag,*dd,*temp_diagS,*temp_diagY ,*temp_diagZ);
	  setAdditiveDiagonal();
	}else if(gOuterSolve >= 3 && separateHandDiag==0){
	  setXDiagonal( *dd );
	  setSDiagonal( *temp_diagS );
	  setYDiagonal( *temp_diagY);		  
	  setZDiagonal( *temp_diagZ );		  
	}else{
	  assert(0);
	}
  }

}

