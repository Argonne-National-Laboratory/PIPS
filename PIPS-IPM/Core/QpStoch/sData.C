#include "sData.h"
#include "sTree.h"
#include "StochSymMatrix.h"
#include "StochGenMatrix.h"
#include "StochVector.h"
#include "QpGenVars.h"
#include "SparseLinearAlgebraPackage.h"
#include "mpi.h"

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

  use2Links = false;
}

sData::sData(sTree* tree_, OoqpVector * c_in, SymMatrix * Q_in,
        OoqpVector * xlow_in, OoqpVector * ixlow_in, long long nxlow_,
        OoqpVector * xupp_in, OoqpVector * ixupp_in, long long nxupp_,
        GenMatrix  * A_in, OoqpVector * bA_in,
        GenMatrix  * C_in,
        OoqpVector * clow_in, OoqpVector * iclow_in, long long mclow_,
        OoqpVector * cupp_in, OoqpVector * icupp_in, long long mcupp_,
        bool exploit2Links)

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

  if( exploit2Links )
  {
     init2LinksData();

     if( use2Links )
     {
        // compute permutation vector
        const int myl = tree_->myl();
        const int mzl = tree_->mzl();

        assert(myl >= 0 && mzl >= 0 && (mzl + myl > 0));

        std::vector<unsigned int> permvecA(myl);
        std::vector<unsigned int> permvecC(mzl);

        for( int i = 0; i < myl; ++i )
           permvecA[i] = myl - i - 1;

        for( int i = 0; i < mzl; ++i )
           permvecC[i] =  mzl - i - 1;
#if 0
        ofstream myfile;
          myfile.open ("C1.txt");
        dynamic_cast<StochGenMatrix&>(*C).writeToStreamDense(myfile);
#endif


        permuteLinkingRows(permvecA, permvecC);

        permuteLinkingRows(permvecA, permvecC);

#if 0
        ofstream myfile2;
           myfile2.open ("C2.txt");
              dynamic_cast<StochGenMatrix&>(*C).writeToStreamDense(myfile2);
#endif
     }
  }
  else
     use2Links = false;
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

void sData::permuteLinkingRows(const std::vector<unsigned int>& permvecA, const std::vector<unsigned int>& permvecC)
{
   dynamic_cast<StochGenMatrix&>(*A).permuteLinkingRows(permvecA);
   dynamic_cast<StochGenMatrix&>(*C).permuteLinkingRows(permvecC);
   dynamic_cast<StochVector&>(*bA).permuteLinkingEntries(permvecA);
   dynamic_cast<StochVector&>(*bl).permuteLinkingEntries(permvecC);
   dynamic_cast<StochVector&>(*bu).permuteLinkingEntries(permvecC);
   dynamic_cast<StochVector&>(*iclow).permuteLinkingEntries(permvecC);
   dynamic_cast<StochVector&>(*icupp).permuteLinkingEntries(permvecC);
}

void sData::init2LinksData()
{
   use2Links = true;

   int myrank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   const StochGenMatrix& Astoch = dynamic_cast<const StochGenMatrix&>(*A);
   const StochGenMatrix& Cstoch = dynamic_cast<const StochGenMatrix&>(*C);

   linkIndicatorA = Astoch.get2LinkIndicator();
   linkIndicatorC = Cstoch.get2LinkIndicator();

   int nEq = 0;
   int nIneq = 0;

   for( size_t i = 0; i < linkIndicatorA.size(); i++ )
      if( linkIndicatorA[i] )
         nEq++;

   for( size_t i = 0; i < linkIndicatorC.size(); i++ )
      if( linkIndicatorC[i] )
         nIneq++;

   if( myrank == 0 )
   {
      std::cout << "number of equality 2-links: " << nEq << " (out of "
            << linkIndicatorA.size() << " equalities) " << std::endl;
      std::cout << "number of inequality 2-links: " << nIneq << " (out of "
            << linkIndicatorC.size() << " equalities) " << std::endl;

      std::cout << "ratio: "
            << (nIneq + nEq) / ((double) linkIndicatorA.size() + linkIndicatorC.size()) << std::endl;
   }

   if( nIneq + nEq / ((double) linkIndicatorA.size() + linkIndicatorC.size()) < min2LinksRatio )
   {
      if( myrank == 0 )
         std::cout << "not enough 2-links found" << std::endl;
      use2Links = false;
   }


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
  long long my, nx;
  StochGenMatrix& Ast = dynamic_cast<StochGenMatrix&>(*A);
  Ast.Bmat->getSize(my, nx);
  return my;
}

int sData::getLocalmyl()
{
  long long myl, nxl;
  StochGenMatrix& Ast = dynamic_cast<StochGenMatrix&>(*A);
  Ast.Blmat->getSize(myl, nxl);
  return myl;
}


int sData::getLocalmz()
{
  long long mz, nx;
  StochGenMatrix& Cst = dynamic_cast<StochGenMatrix&>(*C);
  Cst.Bmat->getSize(mz, nx);
  return mz;
}

int sData::getLocalmzl()
{
  long long mzl, nxl;
  StochGenMatrix& Cst = dynamic_cast<StochGenMatrix&>(*C);
  Cst.Blmat->getSize(mzl, nxl);
  return mzl;
}

int sData::getLocalSizes(int& nx, int& my, int& mz, int& myl, int& mzl)
{
  long long nxloc, myloc, mzloc, mylloc, mzlloc;

  StochGenMatrix& Ast = dynamic_cast<StochGenMatrix&>(*A);
  Ast.Blmat->getSize(mylloc, nxloc);
  Ast.Bmat->getSize(myloc, nxloc);

  StochGenMatrix& Cst = dynamic_cast<StochGenMatrix&>(*C);
  Cst.Blmat->getSize(mzlloc, nxloc);
  Cst.Bmat->getSize(mzloc, nxloc);

  nx=nxloc; my=myloc; mz=mzloc; myl=mylloc; mzl=mzlloc;
  return 0;
}

int sData::getLocalSizes(int& nx, int& my, int& mz)
{
  long long nxll, myll, mzll;

  StochGenMatrix& Ast = dynamic_cast<StochGenMatrix&>(*A);
  Ast.Bmat->getSize(myll, nxll);

  StochGenMatrix& Cst = dynamic_cast<StochGenMatrix&>(*C);
  Cst.Bmat->getSize(mzll, nxll);

  nx=nxll; my=myll; mz=mzll;
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

// This is F_i (linking equality matrix):
SparseGenMatrix& sData::getLocalF()
{
  StochGenMatrix& Ast = dynamic_cast<StochGenMatrix&>(*A);
  return *Ast.Blmat;
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

// This is G_i (linking inequality matrix):
SparseGenMatrix& sData::getLocalG()
{
  StochGenMatrix& Cst = dynamic_cast<StochGenMatrix&>(*C);
  return *Cst.Blmat;
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
