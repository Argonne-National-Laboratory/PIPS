#include "QpGenStochLinsysRoot.h"
#include "StochTree.h"
#include "QpGenStoch.h"
#include "QpGenStochData.h"
#include "QpGenStochLinsysLeaf.h"
/*********************************************************************/
/************************** ROOT *************************************/
/*********************************************************************/

#ifdef STOCH_TESTING
extern double g_iterNumber;
#endif

QpGenStochLinsysRoot::QpGenStochLinsysRoot(QpGenStoch * factory_, QpGenStochData * prob_)
  : QpGenStochLinsys(factory_, prob_), iAmDistrib(0)
{
  createChildren(prob_);
}

QpGenStochLinsysRoot::QpGenStochLinsysRoot(QpGenStoch* factory_,
					   QpGenStochData* prob_,
					   OoqpVector* dd_, 
					   OoqpVector* dq_,
					   OoqpVector* nomegaInv_,
					   OoqpVector* rhs_)
  : QpGenStochLinsys(factory_, prob_, dd_, dq_, nomegaInv_, rhs_), iAmDistrib(0)
{
  createChildren(prob_);
}

QpGenStochLinsysRoot::~QpGenStochLinsysRoot()
{

}


void QpGenStochLinsysRoot::createChildren(QpGenStochData* prob)
{

  QpGenStochLinsys* child=NULL;
  StochVector& ddst = dynamic_cast<StochVector&>(*dd);
  StochVector& dqst = dynamic_cast<StochVector&>(*dq);
  StochVector& nomegaInvst = dynamic_cast<StochVector&>(*nomegaInv);
  StochVector& rhsst = dynamic_cast<StochVector&>(*rhs);

  //get the communicator from one of the vectors
  this->mpiComm = ddst.mpiComm;
  this->iAmDistrib = ddst.iAmDistrib;

  for(int it=0; it<prob->children.size(); it++) {
    
    if(MPI_COMM_NULL == ddst.children[it]->mpiComm) {
      child = new QpGenStochDummyLinsys(dynamic_cast<QpGenStoch*>(factory), prob->children[it]);
    } else {
      QpGenStoch* stochFactory = dynamic_cast<QpGenStoch*>(factory);

      if(prob->children[it]->children.size() == 0) {	
	//child = new QpGenStochLinsysLeaf(dynamic_cast<QpGenStoch*>(factory), 
	child = stochFactory->newLinsysLeaf(prob->children[it],
					    ddst.children[it],
					    dqst.children[it],
					    nomegaInvst.children[it],
					    rhsst.children[it]);
      } else {
	//child = new QpGenStochLinsysRoot(dynamic_cast<QpGenStoch*>(factory), 
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

void QpGenStochLinsysRoot::deleteChildren()
{
  for(int it=0; it<children.size(); it++) {
    children[it]->deleteChildren();
    delete children[it];
  }
  children.clear();
}

void QpGenStochLinsysRoot::putXDiagonal( OoqpVector& xdiag_ )
{
  StochVector& xdiag = dynamic_cast<StochVector&>(xdiag_);
  assert(children.size() == xdiag.children.size());

  //kkt->atPutDiagonal( 0, *xdiag.vec );
  xDiag = xdiag.vec;
 
  // propagate it to the subtree
  for(int it=0; it<children.size(); it++)
    children[it]->putXDiagonal(*xdiag.children[it]);
}


void QpGenStochLinsysRoot::putZDiagonal( OoqpVector& zdiag_ )
{
  StochVector& zdiag = dynamic_cast<StochVector&>(zdiag_);
  assert(children.size() == zdiag.children.size());

  //kkt->atPutDiagonal( locnx+locmy, *zdiag.vec );
  zDiag = zdiag.vec;

  // propagate it to the subtree
  for(int it=0; it<children.size(); it++)
    children[it]->putZDiagonal(*zdiag.children[it]);
}

void QpGenStochLinsysRoot::AddChild(QpGenStochLinsys* child)
{
  children.push_back(child);
}


void QpGenStochLinsysRoot::sync()
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

#ifdef STOCH_TESTING
void QpGenStochLinsysRoot::dumpMatrix(int scen, int proc, const char* nameToken, DenseSymMatrix& M) 
{
  int n = M.size();
  char szNumber[30];
  string strBuffer="";

  assert(false);

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

void QpGenStochLinsysRoot::dumpRhs(int proc, const char* nameToken,  SimpleVector& rhs) 
{
  int n = rhs.length();
  char szNumber[30];
  string strBuffer="";


  int iter = g_iterNumber;
  if(iter!=2 && iter!=20 && iter!=25 && iter!=55) return;

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
