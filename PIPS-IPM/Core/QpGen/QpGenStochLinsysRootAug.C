#include "QpGenStochLinsysRootAug.h"
#include "DeSymIndefSolver.h"
#include "QpGenStochData.h"

QpGenStochLinsysRootAug::QpGenStochLinsysRootAug(QpGenStoch * factory_, QpGenStochData * prob_)
  : QpGenStochLinsysRoot(factory_, prob_), UtV(NULL)
{ 
  prob_->getLocalSizes(locnx, locmy, locmz);

  //kkt = createKKT(prob_);
  //solver = createSolver(prob_, kkt);
};

QpGenStochLinsysRootAug::QpGenStochLinsysRootAug(QpGenStoch* factory_,
					   QpGenStochData* prob_,
					   OoqpVector* dd_, 
					   OoqpVector* dq_,
					   OoqpVector* nomegaInv_,
					   OoqpVector* rhs_)
  : QpGenStochLinsysRoot(factory_, prob_, dd_, dq_, nomegaInv_, rhs_), UtV(NULL)
{ 
  prob_->getLocalSizes(locnx, locmy, locmz);

  //kkt = createKKT(prob_);
  //solver = createSolver(prob_, kkt);
};

QpGenStochLinsysRootAug::~QpGenStochLinsysRootAug()
{
  if (UtV) delete UtV;
}


void QpGenStochLinsysRootAug::Lsolve(QpGenStochData *prob, OoqpVector& x)
{
  StochVector& b = dynamic_cast<StochVector&>(x);
  assert(children.size() == b.children.size() );

  // children compute their part
  for(int it=0; it<children.size(); it++) {
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


  for(int it=0; it<children.size(); it++) {
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


void QpGenStochLinsysRootAug::Ltsolve2( QpGenStochData *prob, StochVector& x, SimpleVector& xp)
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
  for(int it=0; it<children.size(); it++) {
    children[it]->Ltsolve2(prob->children[it], *b.children[it], xi);
  }
}

void QpGenStochLinsysRootAug::Ltsolve( QpGenStochData *prob, OoqpVector& x )
{
  StochVector& b   = dynamic_cast<StochVector&>(x);
  SimpleVector& b0 = dynamic_cast<SimpleVector&>(*b.vec);

  stochNode->resMon.recLtsolveTmLocal_start();
  solver->Ltsolve(b0);
  stochNode->resMon.recLtsolveTmLocal_stop();

  //dumpRhs(0, "sol",  b0);

  SimpleVector& x0 = b0; //just another name, for clarity
  
  // Li^T\bi for each child i. The backsolve needs z0

  for(int it=0; it<children.size(); it++) {
    children[it]->Ltsolve2(prob->children[it], *b.children[it], x0);
  }
}

void QpGenStochLinsysRootAug::Dsolve( QpGenStochData *prob, OoqpVector& x )
{
  StochVector& b = dynamic_cast<StochVector&>(x);

  for(int it=0; it<children.size(); it++) {
    children[it]->Dsolve(prob->children[it], *b.children[it]);
  }

  SimpleVector& b0 = dynamic_cast<SimpleVector&>(*b.vec);

  stochNode->resMon.recDsolveTmLocal_start();

  solver->Dsolve(b0);
  
  stochNode->resMon.recDsolveTmLocal_stop();
}



void QpGenStochLinsysRootAug::allocUtV()
{
  if(UtV==NULL)
    UtV = new DenseSymMatrix(locnx);
}

void QpGenStochLinsysRootAug::initializeUtV()
{
  int special=1; //the special gets the cake
  if(iAmDistrib) {
    //am I the special one?
    int rank; MPI_Comm_rank(mpiComm, &rank);
    int size; MPI_Comm_size(mpiComm, &size);
    //if(rank!=size-1) special=0;
    if(rank!=0) special=0;
    //special=0;
  }

  DenseSymMatrix* kktd = dynamic_cast<DenseSymMatrix*>(kkt);
  int i,j;

  if(special) {
    //get the cake
    int bytesToCopy=locnx*sizeof(double);
    for(i=0;i<locnx; i++) {
      memcpy(UtV->mStorage->M[i], kktd->mStorage->M[i], bytesToCopy);
    }
  } else {
    //put zeros
    myAtPutZeros(UtV, 0,0, locnx, locnx);
  }

}
void QpGenStochLinsysRootAug::deleteUtV()
{
  delete UtV; UtV = NULL; 
}


//faster than DenseSymMatrix::atPutZeros
void QpGenStochLinsysRootAug::myAtPutZeros(DenseSymMatrix* mat, 
					      int row, int col, 
					      int rowExtent, int colExtent)
{
  int nn = mat->size();
  assert( row >= 0 && row + rowExtent <= nn );
  assert( col >= 0 && col + colExtent <= nn );

  double ** M = mat->storage().M;

  for(int j=col; j<col+colExtent; j++) {
      M[row][j] = 0.0;
  }

  int nToCopy = colExtent*sizeof(double);

  for(int i=row+1; i<row+rowExtent; i++) {
    memcpy(M[i]+col, M[row]+col, nToCopy);
  }
}
