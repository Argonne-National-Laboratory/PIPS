#include "sData.h"
#include "sTree.h"
#include "StochSymMatrix.h"
#include "StochGenMatrix.h"
#include "StochVector.h"
#include "QpGenVars.h"
#include "SparseLinearAlgebraPackage.h"

sData::sData(sTree* tree)
//  : QpGenData(NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL)
{
  stochNode = tree;
  Q     = SymMatrixHandle ( tree->createQ() );
  g     = OoqpVectorHandle    ( tree->createc() );

  blx   = OoqpVectorHandle    ( tree->createxlow()  );
  ixlow = OoqpVectorHandle    ( tree->createixlow() );
  bux   = OoqpVectorHandle    ( tree->createxupp()  );
  ixupp = OoqpVectorHandle    ( tree->createixupp() );


  A  = GenMatrixHandle        ( tree->createA() );
  bA = OoqpVectorHandle       ( tree->createb() );


  C     = GenMatrixHandle     ( tree->createC() );  
  bl    = OoqpVectorHandle    ( tree->createclow() );
  iclow = OoqpVectorHandle    ( tree->createiclow() );
  bu    = OoqpVectorHandle    ( tree->createcupp()  );
  icupp = OoqpVectorHandle    ( tree->createicupp() );

  sc = OoqpVectorHandle ( tree->newPrimalVector() );

  nxlow = ixlow->numberOfNonzeros(); 
  nxupp = ixupp->numberOfNonzeros(); 
  mclow = iclow->numberOfNonzeros(); 
  mcupp = icupp->numberOfNonzeros();

  createChildren();
}

sData::sData(sTree* tree_, OoqpVector * c_in, SymMatrix * Q_in,
	     OoqpVector * xlow_in, OoqpVector * ixlow_in, int nxlow_,
	     OoqpVector * xupp_in, OoqpVector * ixupp_in, int nxupp_,
	     GenMatrix  * A_in, OoqpVector * bA_in,
	     GenMatrix  * C_in,
	     OoqpVector * clow_in, OoqpVector * iclow_in, int mclow_,
	     OoqpVector * cupp_in, OoqpVector * icupp_in, int mcupp_ )
  
  : QpGenData(SparseLinearAlgebraPackage::soleInstance(),
	      c_in, Q_in, 
	      xlow_in, ixlow_in, xupp_in, ixupp_in,
	      A_in, bA_in,
	      C_in,
	      clow_in, iclow_in, cupp_in, icupp_in)
{
  nxlow = nxlow_; nxupp = nxupp_; 
  mclow = mclow_; mcupp = mcupp_;
  stochNode = tree_;

  createChildren();
}


void sData::createChildren()
{
  //follow the structure of one of the tree objects and create the same
  //structure for this class, and link this object with the corresponding 
  //vectors and matrices
  StochVector& gSt     = dynamic_cast<StochVector&>(*g);
  StochSymMatrix& QSt  = dynamic_cast<StochSymMatrix&>(*Q);
  
  StochVector& xlowSt  = dynamic_cast<StochVector&>(*blx); 
  StochVector& ixlowSt = dynamic_cast<StochVector&>(*ixlow); 
  StochVector& xuppSt  = dynamic_cast<StochVector&>(*bux); 
  StochVector& ixuppSt = dynamic_cast<StochVector&>(*ixupp);
  StochGenMatrix& ASt  = dynamic_cast<StochGenMatrix&>(*A); 
  StochVector& bASt    = dynamic_cast<StochVector&>(*bA);
  StochGenMatrix& CSt  = dynamic_cast<StochGenMatrix&>(*C);
  StochVector& clowSt  = dynamic_cast<StochVector&>(*bl); 
  StochVector& iclowSt = dynamic_cast<StochVector&>(*iclow);
  StochVector& cuppSt  = dynamic_cast<StochVector&>(*bu); 
  StochVector& icuppSt = dynamic_cast<StochVector&>(*icupp); 
  

  for(size_t it=0; it<gSt.children.size(); it++) {
    AddChild(new sData(stochNode->children[it],
	       gSt.children[it], QSt.children[it],
	       xlowSt.children[it], ixlowSt.children[it], nxlow,
	       xuppSt.children[it], ixuppSt.children[it], nxupp,
	       ASt.children[it], bASt.children[it],
	       CSt.children[it],
	       clowSt.children[it], iclowSt.children[it], mclow,
	       cuppSt.children[it], icuppSt.children[it], mcupp ));
  }

}

void sData::destroyChildren()
{
  for(size_t it=0; it<children.size(); it++) {
    children[it]->destroyChildren();
    delete children[it];
  }
  children.clear();
}

void sData::AddChild(sData* child)
{
  children.push_back(child);
}

double sData::objectiveValue( QpGenVars * vars )
{
  StochVector& x = dynamic_cast<StochVector&>(*vars->x);
  OoqpVectorHandle temp( x.clone() );

  this->getg( *temp );
  this->Qmult( 1.0, *temp, 0.5, *vars->x );

  return temp->dotProductWith( *vars->x );
}

void sData::createScaleFromQ()
{

  assert("Not implemented!" && 0);

  // Stuff the diagonal elements of Q into the vector "sc"
  this->getDiagonalOfQ( *sc);

  // Modifying scVector is equivalent to modifying sc
  /*SimpleVector & scVector = dynamic_cast<SimpleVector &>(*sc);

  int scLength = scVector.length();

  for( int i = 0; i < scLength; i++){
    if( scVector[i] > 1)
        scVector[i] = 1.0/sqrt( scVector[i]);
    else
        scVector[i] = 1.0;
    }
  */
}



sData::~sData()
{
  for(size_t it=0; it<children.size(); it++)
    delete children[it];
}

int sData::getLocalnx()
{
  StochSymMatrix& Qst = dynamic_cast<StochSymMatrix&>(*Q);
  return Qst.diag->size();
}

int sData::getLocalmy()
{
  int my, nx;
  StochGenMatrix& Ast = dynamic_cast<StochGenMatrix&>(*A);
  Ast.Bmat->getSize(my, nx);
  return my;
}

int sData::getLocalmz()
{
  int mz, nx;
  StochGenMatrix& Cst = dynamic_cast<StochGenMatrix&>(*C);
  Cst.Bmat->getSize(mz, nx);
  return mz;
}

int sData::getLocalSizes(int& nx, int& my, int& mz)
{
  StochGenMatrix& Ast = dynamic_cast<StochGenMatrix&>(*A);
  Ast.Bmat->getSize(my, nx);

  StochGenMatrix& Cst = dynamic_cast<StochGenMatrix&>(*C);
  Cst.Bmat->getSize(mz, nx);
  return 0;
}


int sData::getLocalNnz(int& nnzQ, int& nnzB, int& nnzD)
{
  StochSymMatrix& Qst = dynamic_cast<StochSymMatrix&>(*Q);
  StochGenMatrix& Ast = dynamic_cast<StochGenMatrix&>(*A);
  StochGenMatrix& Cst = dynamic_cast<StochGenMatrix&>(*C);

  nnzQ = Qst.diag->getStorageRef().len + Qst.border->getStorageRef().len;
  nnzB = Ast.Bmat->getStorageRef().len;
  nnzD = Cst.Bmat->getStorageRef().len;
  return 0;
}


SparseSymMatrix& sData::getLocalQ()
{
  StochSymMatrix& Qst = dynamic_cast<StochSymMatrix&>(*Q);
  return *Qst.diag;
}

SparseGenMatrix& sData::getLocalCrossHessian()
{
  StochSymMatrix& Qst = dynamic_cast<StochSymMatrix&>(*Q);
  return *Qst.border;
}

// T_i x_0 + W_i x_i = b_i

// This is T_i
SparseGenMatrix& sData::getLocalA()
{
  StochGenMatrix& Ast = dynamic_cast<StochGenMatrix&>(*A);
  return *Ast.Amat;
}

// This is W_i:
SparseGenMatrix& sData::getLocalB()
{
  StochGenMatrix& Ast = dynamic_cast<StochGenMatrix&>(*A);
  return *Ast.Bmat;
}

// low_i <= C_i x_0 + D_i x_i <= upp_i

// This is C_i
SparseGenMatrix& sData::getLocalC()
{
  StochGenMatrix& Cst = dynamic_cast<StochGenMatrix&>(*C);
  return *Cst.Amat;
}

// This is D_i
SparseGenMatrix& sData::getLocalD()
{
  StochGenMatrix& Cst = dynamic_cast<StochGenMatrix&>(*C);
  return *Cst.Bmat;
}



void sData::sync()
{

  
  destroyChildren();
  
  stochNode->syncPrimalVector(dynamic_cast<StochVector&>(*g)); 

//   int myRank; MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
//   printf("vec -----------------------------------------------------\n");sleep(myRank+1);  
//   stochNode->displayVectorVsTreeStructure(dynamic_cast<StochVector&>(*g), myRank);
//   printf("vec done ----------------------\n"); usleep(10000);

  stochNode->syncPrimalVector(dynamic_cast<StochVector&>(*blx));
  stochNode->syncPrimalVector(dynamic_cast<StochVector&>(*ixlow));
  stochNode->syncPrimalVector(dynamic_cast<StochVector&>(*bux));
  stochNode->syncPrimalVector(dynamic_cast<StochVector&>(*ixupp));
  stochNode->syncDualYVector(dynamic_cast<StochVector&>(*bA));
  stochNode->syncDualZVector(dynamic_cast<StochVector&>(*bl));
  stochNode->syncDualZVector(dynamic_cast<StochVector&>(*bu));
  stochNode->syncDualZVector(dynamic_cast<StochVector&>(*iclow));
  stochNode->syncDualZVector(dynamic_cast<StochVector&>(*icupp));

  stochNode->syncStochSymMatrix(dynamic_cast<StochSymMatrix&>(*Q));
  stochNode->syncStochGenMatrix(dynamic_cast<StochGenMatrix&>(*A));

  //sleep(myRank);printf("A mat------------------------------------------------\n");
  //stochNode->displayMatVsTreeStructure(dynamic_cast<StochGenMatrix&>(*A), myRank);

  stochNode->syncStochGenMatrix(dynamic_cast<StochGenMatrix&>(*C));

  createChildren();
}
