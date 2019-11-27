#include "QpGenStochLinsysRootAugExt.h"
#include "DeSymIndefSolver.h"
#include "QpGenStochData.h"
#include "pipsport.h"

QpGenStochLinsysRootAugExt::QpGenStochLinsysRootAugExt(QpGenStoch * factory_, QpGenStochData * prob_)
  : QpGenStochLinsysRootAug(factory_, prob_)
{ 
  kkt = createKKT(prob_);
  solver = createSolver(prob_, kkt);
};

QpGenStochLinsysRootAugExt::QpGenStochLinsysRootAugExt(QpGenStoch* factory_,
					   QpGenStochData* prob_,
					   OoqpVector* dd_, 
					   OoqpVector* dq_,
					   OoqpVector* nomegaInv_,
					   OoqpVector* rhs_)
  : QpGenStochLinsysRootAug(factory_, prob_, dd_, dq_, nomegaInv_, rhs_)
{ 
  kkt = createKKT(prob_);
  solver = createSolver(prob_, kkt);
};

QpGenStochLinsysRootAugExt::~QpGenStochLinsysRootAugExt()
{ }

SymMatrix* 
QpGenStochLinsysRootAugExt::createKKT(QpGenStochData* prob)
{
  int n = locnx+locmy+locmz;

  return new DenseSymMatrix(n);
}

DoubleLinearSolver*
QpGenStochLinsysRootAugExt::createSolver(QpGenStochData* prob, SymMatrix* kktmat_)
{

  DenseSymMatrix* kktmat = dynamic_cast<DenseSymMatrix*>(kktmat_);
  return new DeSymIndefSolver(kktmat);
}


void QpGenStochLinsysRootAugExt::factor2(QpGenStochData *prob, Variables *vars)
{
  assert( children.size() == prob->children.size() );

  stochNode->resMon.recFactTmLocal_start();
  stochNode->resMon.recSchurMultLocal_start();
  
  // Diagonals were already updated but KKT was not formed (fixme)
  DenseSymMatrix * kktd = (DenseSymMatrix*) kkt; 
  kktd->symAtPutSubmatrix( 0, 0, prob->getLocalQ(), 0, 0, locnx, locnx);
  kktd->symAtPutSubmatrix( locnx, 0, prob->getLocalB(), 0, 0, locmy, locnx,1);
  kktd->symAtPutSubmatrix( locnx+locmy, 0, prob->getLocalD(), 0, 0, locmz, locnx,1);


  //kktd->storage().atPutZeros(locnx, locnx, locmy+locmz, locmy+locmz);
  myAtPutZeros(kktd, locnx, locnx, locmy+locmz, locmy+locmz);

  kktd->atPutDiagonal( 0, *xDiag );
  if(locmz) kktd->atPutDiagonal( locnx+locmy, *zDiag );


  stochNode->resMon.recSchurMultLocal_stop();
  stochNode->resMon.recFactTmLocal_stop();

  //dumpMatrix(-1, 0, "M", *kktd);  

  // First tell children to factorize.
  for(int it=0; it<children.size(); it++) {
    children[it]->factor2(prob->children[it], vars);
  }
  
  //!parallel communication here
  
  // each child computes   Gi Li^{-T} Di^{-1} Li^{-1} Gi^T  as
  // BLAS3 mat-mat product of U^T V, where
  //          - U=Li^{-1} Gi^T
  //          - V=Di^{-1} U

  DenseGenMatrix* U = nullptr;
  DenseGenMatrix* V = nullptr;
  

  stochNode->resMon.recSchurMultLocal_start();

  allocUtV();
  initializeUtV();

  stochNode->resMon.recSchurMultLocal_stop();

  for(int it=0; it<children.size(); it++) {

    // Allocate  U if necessary. If the sizes of 
    // the previous child are good, reuse the mem.

    if(children[it]->mpiComm == MPI_COMM_NULL)
      continue;
    children[it]->stochNode->resMon.recFactTmChildren_start();    

    children[it]->allocU(&U, locnx);
    children[it]->allocV(&V, locnx);

    children[it]->computeU_V(prob->children[it], U, V);
    children[it]->stochNode->resMon.recSchurMultChildren_start();

    //!dumping data
    //DenseSymMatrix* UtV0 = new DenseSymMatrix(locnx);
    //UtV0->scalarMult(0.0);
    //UtV0->matMult(-1.0, *U, 1, *V, 0, 1.0); 
    //dumpMatrix(it, 0, "C", *UtV0);
    //delete UtV0;
    //~dumping data

    UtV->matMult(-1.0, *U, 1, *V, 0, 1.0); 
    children[it]->stochNode->resMon.recSchurMultChildren_stop();

    children[it]->stochNode->resMon.recFactTmChildren_stop();
  }

  stochNode->resMon.recSchurMultLocal_start();
  delete U; delete V;
  //deleteUtV(); reuse this
  stochNode->resMon.recSchurMultLocal_stop();

  //parallel communication
  if(iAmDistrib) {

    double* buffer=new double[locnx*locnx];
    MPI_Allreduce(&(UtV->getStorageRef().M[0][0]),
		  buffer,
		  locnx*locnx,
		  MPI_DOUBLE,
		  MPI_SUM,
		  MPI_COMM_WORLD);

    memcpy(&UtV->mStorage->M[0][0], buffer, locnx*locnx*sizeof(double));

    delete [] buffer;
  }

  stochNode->resMon.recFactTmLocal_start();
  stochNode->resMon.recSchurMultLocal_start();

  //update the upper block of the kkt with the UtV block
  kkt->symAtPutSubmatrix(0, 0, *UtV, 0, 0, locnx, locnx);

  stochNode->resMon.recSchurMultLocal_stop();  

  //dumpMatrix(-1, 0, "S", *kktd);

  
  double tt= MPI_Wtime();
  // just trigger a local refactorization 
  // (if needed, depends on the type of lin solver).
  solver->matrixChanged();

  printf("Fact took %g seconds.\n", MPI_Wtime()-tt);
  stochNode->resMon.recFactTmLocal_stop();  
}

