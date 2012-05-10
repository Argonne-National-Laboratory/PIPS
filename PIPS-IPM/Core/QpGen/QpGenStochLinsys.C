#include "QpGenStochLinsys.h"
#include "StochTree.h"
#include "QpGenStoch.h"
#include "QpGenStochData.h"
#include "SparseLinearAlgebraPackage.h"


QpGenStochLinsys::QpGenStochLinsys(QpGenStoch* factory_, QpGenStochData* prob)
  : QpGenLinsys(), kkt(NULL), solver(NULL)
{
  factory = factory_;

  nx = prob->nx; my = prob->my; mz = prob->mz;
  ixlow = prob->ixlow;
  ixupp = prob->ixupp;
  iclow = prob->iclow;
  icupp = prob->icupp;

  nxlow = prob->nxlow;
  nxupp = prob->nxupp;
  mclow = prob->mclow;
  mcupp = prob->mcupp;

  if( nxupp + nxlow > 0 ) {
    dd      = factory_->tree->newPrimalVector();
    dq      = factory_->tree->newPrimalVector();
    prob->getDiagonalOfQ( *dq );
  }
  nomegaInv   = factory_->tree->newDualZVector();
  rhs         = factory_->tree->newRhs();

  useRefs=0;
  data = prob;
  stochNode = prob->stochNode;
}

QpGenStochLinsys::QpGenStochLinsys(QpGenStoch* factory_,
				   QpGenStochData* prob,				    
				   OoqpVector* dd_, 
				   OoqpVector* dq_,
				   OoqpVector* nomegaInv_,
				   OoqpVector* rhs_)
  : QpGenLinsys(), kkt(NULL), solver(NULL)
{
  factory = factory_;


  nx = prob->nx; my = prob->my; mz = prob->mz;
  ixlow = prob->ixlow;
  ixupp = prob->ixupp;
  iclow = prob->iclow;
  icupp = prob->icupp;

  nxlow = prob->nxlow;
  nxupp = prob->nxupp;
  mclow = prob->mclow;
  mcupp = prob->mcupp;

  if( nxupp + nxlow > 0 ) {
    //dd      = OoqpVectorHandle(dd_);
    dd= dd_;
    //dq      = OoqpVectorHandle(dq_);
    dq = dq_;
  }
  //nomegaInv   = OoqpVectorHandle(nomegaInv_);
  //rhs         = OoqpVectorHandle(rhs_);
  nomegaInv = nomegaInv_;
  rhs = rhs_;

  useRefs=1;
  data = prob;
  stochNode = prob->stochNode;
}


QpGenStochLinsys::~QpGenStochLinsys()
{
  if (solver) delete solver;
  if (kkt)    delete kkt;
}

void QpGenStochLinsys::joinRHS( OoqpVector& rhs_in,  OoqpVector& rhs1_in,
				OoqpVector& rhs2_in, OoqpVector& rhs3_in )
{
  StochVector& rhs  = dynamic_cast<StochVector&>(rhs_in);
  StochVector& rhs1 = dynamic_cast<StochVector&>(rhs1_in);
  StochVector& rhs2 = dynamic_cast<StochVector&>(rhs2_in);
  StochVector& rhs3 = dynamic_cast<StochVector&>(rhs3_in);

  rhs.jointCopyFrom(rhs1, rhs2, rhs3);
}

void QpGenStochLinsys::separateVars( OoqpVector& x_in, OoqpVector& y_in,
				     OoqpVector& z_in, OoqpVector& vars_in )
{
  StochVector& x    = dynamic_cast<StochVector&>(x_in);
  StochVector& y    = dynamic_cast<StochVector&>(y_in);
  StochVector& z    = dynamic_cast<StochVector&>(z_in);
  StochVector& vars = dynamic_cast<StochVector&>(vars_in);

  vars.jointCopyTo(x, y, z);
}




void QpGenStochLinsys::factor(Data *prob_, Variables *vars)
{
  // the call to the the parent's method takes care of all necessary updates
  // to the KKT system (updating diagonals mainly). This is done reccursevely,
  // we don't have to worry about it anymore. 
  QpGenLinsys::factor(prob_, vars);

  // now DO THE LINEAR ALGEBRA!
  
  QpGenStochData* prob = dynamic_cast<QpGenStochData*>(prob_);
  // in order to avoid a call to QpGenLinsys::factor, call factor2 method.
  factor2(prob, vars);
}
 

/**
 * Computes U = Li\Gi^T.
 *        [ 0 0 0 ]
 * Gi^T = [ A 0 0 ]
 *        [ C 0 0]
 *
 * We have special structure here:
 *             [ 0 ]
 *   U   = Li\ [ A ] ,   U is (nx+my+mz)-by-(np)
 *             [ C ]
 *
 *   V = Di\U
 */
void QpGenStochLinsys::computeU_V(QpGenStochData *prob, 
				  DenseGenMatrix* U, DenseGenMatrix* V)
{
  U->scalarMult(0.0);
  V->scalarMult(0.0);

  SparseGenMatrix& A = prob->getLocalA();
  SparseGenMatrix& C = prob->getLocalC();

  int N, nxP;
  A.getSize(N, nxP); assert(N==locmy);

  N = locnx+locmy+locmz;
  SimpleVector uCol(N);

  for(int it=0; it<nxP; it++) {
    
    double* p = &uCol[0];
    for(int it1=0; it1<locnx; it1++) p[it1]=0.0;

    A.fromGetDense(0, it, &uCol[locnx],  1, locmy, 1);    
    C.fromGetDense(0, it, &uCol[locnx+locmy], 1, locmz, 1);

    solver->Lsolve(uCol);
    U->atPutDense(0, it, &uCol[0], 1, N, 1);

    solver->Dsolve(uCol);
    V->atPutDense(0, it, &uCol[0], 1, N, 1);
  } 
}

void QpGenStochLinsys::allocU(DenseGenMatrix ** U, int n0)
{
  int lines,cols;
  if(*U==NULL) {
    *U = new DenseGenMatrix(locnx+locmy+locmz, n0);
  } else {
    (*U)->getSize(lines,cols);
    
    if(lines!=locnx+locmy+locmz || cols != n0) {

      delete (*U);
      *U = new DenseGenMatrix(locnx+locmy+locmz, n0);
    }
  }
}

void QpGenStochLinsys::allocV(DenseGenMatrix ** V, int n0)
{
  int lines,cols;
  if(*V==NULL)
    *V = new DenseGenMatrix(locnx+locmy+locmz, n0);
  else {
    (*V)->getSize(lines,cols);
    
    if(lines!=locnx+locmy+locmz || cols != n0) {
      
      delete (*V);
      *V = new DenseGenMatrix(locnx+locmy+locmz, n0);
    }
  }
}



/**
 *       [ 0 Ai^T Ci^T ]          [    ]
 * b0 -= [ 0   0   0   ] * Li\Di\ [ zi ]
 *       [ 0   0   0   ]          [    ]
 *
 * Changes only the first n0 entries of b0
 */
void QpGenStochLinsys::addLnizi(QpGenStochData *prob, OoqpVector& z0_, OoqpVector& zi_)
{
  SimpleVector& z0 = dynamic_cast<SimpleVector&>(z0_);
  SimpleVector& zi = dynamic_cast<SimpleVector&>(zi_);

  solver->Dsolve (zi);
  solver->Ltsolve(zi);

  SparseGenMatrix& A = prob->getLocalA();
  SparseGenMatrix& C = prob->getLocalC();

  //get n0= nx(parent)= #cols of A or C
  int dummy, n0;
  A.getSize(dummy, n0);

  // zi2 and zi3 are just references to fragements of zi
  SimpleVector zi2 (&zi[locnx], locmy);
  SimpleVector zi3 (&zi[locnx+locmy], locmz);
  // same for z01
  SimpleVector z01 (&z0[0], n0);
  
  A.transMult(1.0, z01, -1.0, zi2);
  C.transMult(1.0, z01, -1.0, zi3);


  if(stochNode->id()==-12) {
      printf("sleeping\n");
      sleep(10);
      //usleep(10000000);
  }
}



void QpGenStochLinsys::solveCompressed( OoqpVector& rhs_ )
{
  StochVector& rhs = dynamic_cast<StochVector&>(rhs_);
  //!log
  int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  //int size; MPI_Comm_size(MPI_COMM_WORLD, &size);
  

  //!opt - mix execution of the following calls

  Lsolve (data, rhs); 
  //sleep(rank);
  //if(rank>=0)
  //{printf("-----------------\n");rhs.writeToStream(cout);printf("~~~~~~~~~~~~~~~~~~~~~~~~~~\n");}
  //!log

  Dsolve(data,rhs);
  //sleep(rank);
  //if(rank>=0)
  //{printf("-----------------\n");rhs.writeToStream(cout);printf("~~~~~~~~~~~~~~~~~~~~~~~~~~\n");}

 
  Ltsolve(data,rhs);
  //if(rank==0) {
    //printf("Lsolve\n");   rhs.writeToStream(cout);
  //}
}


/*
 *  y = alpha*Lni^T x + beta*y
 *
 *                       ( [ 0 0 0 ]     )
 *  y = beta*y + Di\Li\ (  [ A 0 0 ] * x )
 *                      (  [ C 0 0 ]    )
 */
void QpGenStochLinsys::LniTransMult(QpGenStochData *prob, 
				    SimpleVector& y, 
				    double alpha, SimpleVector& x)
{
  SparseGenMatrix& A = prob->getLocalA();
  SparseGenMatrix& C = prob->getLocalC();

  int N, nx0;

  //get nx(parent) from the number of cols of A (or C). Use N as dummy
  A.getSize(N,nx0);
  // a mild assert
  assert(nx0 <= x.length());

  N = locnx+locmy+locmz;
  assert( y.length() == N);
  
  //!memopt
  SimpleVector LniTx(N);

  // shortcuts
  SimpleVector x1(&x[0], nx0);
  SimpleVector LniTx1(&LniTx[0], locnx);
  SimpleVector LniTx2(&LniTx[locnx], locmy);
  SimpleVector LniTx3(&LniTx[locnx+locmy], locmz);
  
  LniTx1.setToZero();
  A.mult(0.0, LniTx2, 1.0, x1);
  C.mult(0.0, LniTx3, 1.0, x1);
 
  solver->Lsolve(LniTx); 
  solver->Dsolve(LniTx);

  //!execopt : does LniTx have the first nx entries zero ?
  y.axpy(alpha,LniTx); 
 
}

