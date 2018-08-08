#include "sData.h"
#include "sTree.h"
#include "sTreeCallbacks.h"
#include "StochSymMatrix.h"
#include "StochGenMatrix.h"
#include "StochVector.h"
#include "QpGenVars.h"
#include "SparseLinearAlgebraPackage.h"
#include "mpi.h"


static inline
int nnzTriangular(int size)
{
   assert(size >= 0);
   return ((1 + size) * size) / 2;
}

static
void appendDenseBlock(int start, int end, int& nnz, int* jcolM)
{
   assert(jcolM);
   assert(nnz >= 0);
   assert(start >= 0 && end >= start);

   for( int i = start; i < end; i++ )
      jcolM[nnz++] = i;
}

static
int appendDiagBlocks(const std::vector<int>& linkStartBlocks, const std::vector<int>& linkStartBlockLengths, int borderstart, int bordersize, int rowSC,
                                          int rowBlock, int& blockStartrow, int& nnz, int* jcolM)
{
   assert(rowBlock >= blockStartrow && blockStartrow >= 0 && borderstart >= 0 && bordersize >= 0 && nnz >= 0);

   const int block = linkStartBlocks[rowBlock];
   const int currlength = (block >= 0) ? linkStartBlockLengths[block] : bordersize;

   assert(currlength >= 1);

   // add diagonal block (possibly up to the order)

   int rownnz = currlength - (rowBlock - blockStartrow);

   for( int i = 0; i < rownnz; ++i )
      jcolM[nnz++] = rowSC + i;

   // with offdiagonal blocks?
   if( block >= 0 )
   {
      // add right off-diagonal block and border part

      const int nextlength = linkStartBlockLengths[block + 1];

      assert(nextlength >= 0);
      assert(block != int(linkStartBlockLengths.size()) - 2 || nextlength == 0);

      for( int i = rownnz; i < rownnz + nextlength; ++i )
         jcolM[nnz++] = rowSC + i;

      rownnz += nextlength + bordersize;

      for( int i = borderstart; i < borderstart + bordersize; ++i )
         jcolM[nnz++] = i;

      // last row of current block?
      if( rowBlock + 1 == blockStartrow + currlength )
         blockStartrow = rowBlock + 1;
   }

   return rownnz;
}


int sData::getSchurCompMaxNnz(const std::vector<int>& linkStartBlocks, const std::vector<int>& linkStartBlockLengths)
{
   const size_t nRows = linkStartBlocks.size();
   const size_t nBlocks = linkStartBlockLengths.size();
   size_t nRowsSparse = 0;

   int nnz = 0;

   // main loop, going over all 2-link blocks
   for( size_t block = 0; block < nBlocks; ++block )
   {
      if( linkStartBlockLengths[block] == 0 )
         continue;

      const int n2links = linkStartBlockLengths[block];
      const int nextlength = linkStartBlockLengths[block + 1];

      assert(n2links > 0);
      assert(nextlength >= 0);
      assert(block != linkStartBlockLengths.size() - 2 || nextlength == 0);

      nRowsSparse += size_t(n2links);

      // diagonal block
      nnz += nnzTriangular(n2links);

      // (one) off-diagonal block
      nnz += n2links * nextlength;
   }

   // any rows left?
   if( nRowsSparse < nRows )
   {
      const size_t nRowsDense = nRows - nRowsSparse;
      nnz += nnzTriangular(nRowsDense) + nRowsDense * nRowsSparse;
   }

   return nnz;
}


int sData::n2linksRows(const std::vector<int>& linkStartBlockLengths)
{
   int n = 0;

   for( size_t i = 0; i < linkStartBlockLengths.size(); ++i )
      n += linkStartBlockLengths[i];

   return n;
}

std::vector<int> sData::get2LinkLengthsVec(const std::vector<int>& linkStartBlocks, size_t nBlocks)
{
   std::vector<int> linkStartBlockLengths(nBlocks, 0);

   const size_t nlinks = linkStartBlocks.size();

   for( size_t i = 0; i < nlinks; i++ )
   {
      const int block = linkStartBlocks[i];

      if( block >= 0 )
      {
         assert(size_t(block) < nBlocks);
         linkStartBlockLengths[block]++;
      }
   }

   return linkStartBlockLengths;
}

SparseSymMatrix* sData::createSchurCompSymbSparseUpper()
{
   assert(children.size() > 0);

   const int nx0 = getLocalnx();
   const int my0 = getLocalmy();
   const int myl = getLocalmyl();
   const int mzl = getLocalmzl();
   const int sizeSC = nx0 + my0 + myl + mzl;
   const int nnz = getSchurCompMaxNnz();

   assert(nnz > 0);
   assert(myl >= 0 && mzl >= 0);

   int* krowM = new int[sizeSC + 1];
   int* jcolM = new int[nnz];
   double* M = new double[nnz];

   krowM[0] = 0;

   // get B_0^T (resp. A_0^T)
   SparseGenMatrix& Btrans = getLocalB().getTranspose();
   int* const startRowBtrans = Btrans.krowM();
   int* const colidxBtrans = Btrans.jcolM();

#ifndef NDEBUG
      int bm, bn;
      Btrans.getSize(bm, bn);
      assert(bm == nx0 && bn == my0);
#endif

   int nnzcount = 0;

   // dense square block, B_0^T, and dense border blocks
   for( int i = 0; i < nx0; ++i )
   {
      const int blength = startRowBtrans[i + 1] - startRowBtrans[i];
      assert(blength >= 0);

      krowM[i + 1] = krowM[i] + (nx0 - i) + blength + myl + mzl;

      appendDenseBlock(i, nx0, nnzcount, jcolM);

      for( int cb = startRowBtrans[i]; cb < startRowBtrans[i + 1]; cb++ )
         jcolM[nnzcount++] = nx0 + colidxBtrans[cb];

      appendDenseBlock(nx0 + my0, nx0 + my0 + myl + mzl, nnzcount, jcolM);

      assert(nnzcount == krowM[i + 1]);
   }

   // empty rows; put diagonal for PARDISO
   for( int i = nx0; i < nx0 + my0; ++i )
   {
      const int rowStartIdx = krowM[i];

      jcolM[rowStartIdx] = i;
      krowM[i + 1] = rowStartIdx + 1;
   }

   nnzcount += my0;

   // equality linking: sparse diagonal blocks, and mixed rows
   int blockStartrow = 0;
   const int n2linksRowsEq = n2linksRows(linkStartBlockLengthsA);
   const int bordersizeEq = linkStartBlocksA.size() - n2linksRowsEq;
   const int borderstartEq = nx0 + my0 + n2linksRowsEq;

   assert(bordersizeEq >= 0 && n2linksRowsEq <= myl);

   // todo replace mzl for sparse 2-link ink linking cons (G)
   for( int i = nx0 + my0, j = 0; i < nx0 + my0 + myl; ++i, ++j )
   {
       const int blockrownnz = appendDiagBlocks(linkStartBlocksA, linkStartBlockLengthsA, borderstartEq, bordersizeEq, i, j, blockStartrow, nnzcount, jcolM);

       appendDenseBlock(nx0 + my0 + myl, nx0 + my0 + myl + mzl, nnzcount, jcolM);
       krowM[i + 1] = krowM[i] + blockrownnz + mzl;
   }

   // inequality linking: dense border block and sparse diagonal blocks
   blockStartrow = 0;
   const int n2linksRowsIneq = n2linksRows(linkStartBlockLengthsC);
   const int bordersizeIneq = linkStartBlocksC.size() - n2linksRowsIneq;
   const int borderstartIneq = nx0 + my0 + myl + n2linksRowsIneq;

   assert(bordersizeIneq >= 0 && n2linksRowsIneq <= mzl);

   for( int i = nx0 + my0 + myl, j = 0; i < nx0 + my0 + myl + mzl; ++i, ++j )
   {
       const int blockrownnz = appendDiagBlocks(linkStartBlocksC, linkStartBlockLengthsC, borderstartIneq, bordersizeIneq, i, j, blockStartrow, nnzcount, jcolM);

       krowM[i + 1] = krowM[i] + blockrownnz;
   }

   assert(nnzcount == nnz);

   return (new SparseSymMatrix(sizeSC, nnz, krowM, jcolM, M, 0, false));
}

std::vector<unsigned int> sData::getAscending2LinkPermutation(std::vector<int>& linkStartBlocks, size_t nBlocks)
{
   const size_t size = linkStartBlocks.size();
   assert(size > 0);

   std::vector<unsigned int> permvec(size, 0);
   std::vector<int> w(nBlocks + 1, 0);

   for( size_t i = 0; i < size; ++i )
   {
      assert(linkStartBlocks[i] >= - 1 && linkStartBlocks[i] < int(nBlocks));
      w[linkStartBlocks[i] + 1]++;
   }

   // initialize start pointers
   int sum = 0;
   for( size_t i = 1; i <= nBlocks; ++i )
   {
      sum += w[i];
      w[i] = sum;
   }

   assert(unsigned(sum + w[0]) == size);

   w[0] = 0;

   for( size_t i = 0; i < size; ++i )
   {
      const int startBlock = (linkStartBlocks[i] >= 0) ? linkStartBlocks[i] : int(nBlocks);

      assert(w[startBlock] <= int(size));
      assert(permvec[w[startBlock]] == 0);

      permvec[w[startBlock]] = i;
      w[startBlock]++;
   }

#ifndef NDEBUG
     for( size_t i = 1; i < permvec.size(); i++ )
        assert(linkStartBlocks[permvec[i]] == - 1 || linkStartBlocks[permvec[i - 1]] <=  linkStartBlocks[permvec[i]]);
#endif

   // permute linkStartBlocks
   std::vector<int> tmpvec(size);

   for( size_t i = 0; i < size; ++i )
      tmpvec[i] = linkStartBlocks[permvec[i]];

   linkStartBlocks = tmpvec;

   return permvec;
}

sData::sData(sTree* tree)
//  : QpGenData(NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL)
{
   stochNode = tree;
   Q = SymMatrixHandle(tree->createQ());
   g = OoqpVectorHandle(tree->createc());

   blx = OoqpVectorHandle(tree->createxlow());
   ixlow = OoqpVectorHandle(tree->createixlow());
   bux = OoqpVectorHandle(tree->createxupp());
   ixupp = OoqpVectorHandle(tree->createixupp());

   A = GenMatrixHandle(tree->createA());
   bA = OoqpVectorHandle(tree->createb());

   C = GenMatrixHandle(tree->createC());
   bl = OoqpVectorHandle(tree->createclow());
   iclow = OoqpVectorHandle(tree->createiclow());
   bu = OoqpVectorHandle(tree->createcupp());
   icupp = OoqpVectorHandle(tree->createicupp());

   sc = OoqpVectorHandle(tree->newPrimalVector());

   nxlow = ixlow->numberOfNonzeros();
   nxupp = ixupp->numberOfNonzeros();
   mclow = iclow->numberOfNonzeros();
   mcupp = icupp->numberOfNonzeros();

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

  init2LinksData(exploit2Links);

  if( use2Links )
  {
     assert(linkStartBlocksA.size() == unsigned(tree_->myl()));
     assert(linkStartBlocksC.size() == unsigned(tree_->mzl()));

     // compute permutation vector

     const size_t nBlocks = dynamic_cast<StochVector*>(c_in)->children.size();

     std::vector<unsigned int> permvecA = getAscending2LinkPermutation(linkStartBlocksA, nBlocks);
     std::vector<unsigned int> permvecC = getAscending2LinkPermutation(linkStartBlocksC, nBlocks);

#ifndef NDEBUG
     const int myl = tree_->myl();
     const int mzl = tree_->mzl();
     assert(myl >= 0 && mzl >= 0 && (mzl + myl > 0));

#if 0
     std::vector<unsigned int> permvecA(myl);
     std::vector<unsigned int> permvecC(mzl);

     for( int i = 0; i < myl; ++i )
        permvecA[i] = myl - i - 1;

     for( int i = 0; i < mzl; ++i )
        permvecC[i] =  mzl - i - 1;

     ofstream myfile;
       myfile.open ("C1.txt");
     dynamic_cast<StochGenMatrix&>(*C).writeToStreamDense(myfile);
#endif

#endif

     permuteLinkingRows(permvecA, permvecC);

#if 0
     ofstream myfile2;
        myfile2.open ("C2.txt");
           dynamic_cast<StochGenMatrix&>(*C).writeToStreamDense(myfile2);
#endif
  }
}

void sData::writeToStreamDense(ostream& out) const
{
   int myRank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

   if( myRank == 0 ) out <<  "A: " << std::endl;
   (*A).writeToStreamDense(out);
   if( myRank == 0 ) out <<  "C: " << std::endl;
   (*C).writeToStreamDense(out);
   if( myRank == 0 ) out <<  "obj: " << std::endl;
   (*g).writeToStreamAll(out);
   if( myRank == 0 ) out <<  "bA: " << std::endl;
   (*bA).writeToStreamAll(out);
   if( myRank == 0 ) out <<  "xupp: " << std::endl;
   (*bux).writeToStreamAll(out);
   if( myRank == 0 ) out <<  "ixupp: " << std::endl;
   (*ixupp).writeToStreamAll(out);
   if( myRank == 0 ) out <<  "xlow: " << std::endl;
   (*blx).writeToStreamAll(out);
   if( myRank == 0 ) out <<  "ixlow: " << std::endl;
   (*ixlow).writeToStreamAll(out);
   if( myRank == 0 ) out <<  "cupp: " << std::endl;
   (*bu).writeToStreamAll(out);
   if( myRank == 0 ) out <<  "icupp: " << std::endl;
   (*icupp).writeToStreamAll(out);
   if( myRank == 0 ) out <<  "clow: " << std::endl;
   (*bl).writeToStreamAll(out);
   if( myRank == 0 ) out <<  "iclow: " << std::endl;
   (*iclow).writeToStreamAll(out);
}

sData*
sData::cloneFull(bool switchToDynamicStorage) const
{
   // todo Q is empty!
   StochSymMatrixHandle Q_clone(dynamic_cast<const StochSymMatrix&>(*Q).clone());
   StochGenMatrixHandle A_clone(dynamic_cast<const StochGenMatrix&>(*A).cloneFull(switchToDynamicStorage));
   StochGenMatrixHandle C_clone(dynamic_cast<const StochGenMatrix&>(*C).cloneFull(switchToDynamicStorage));

   StochVectorHandle c_clone (dynamic_cast<const StochVector&>(*g).cloneFull());
   StochVectorHandle bA_clone ( dynamic_cast<const StochVector&>(*bA).cloneFull());
   StochVectorHandle xupp_clone (dynamic_cast<const StochVector&>(*bux).cloneFull());
   StochVectorHandle ixupp_clone (dynamic_cast<const StochVector&>(*ixupp).cloneFull());
   StochVectorHandle xlow_clone ( dynamic_cast<const StochVector&>(*blx).cloneFull());
   StochVectorHandle ixlow_clone ( dynamic_cast<const StochVector&>(*ixlow).cloneFull());
   StochVectorHandle cupp_clone ( dynamic_cast<const StochVector&>(*bu).cloneFull());
   StochVectorHandle icupp_clone ( dynamic_cast<const StochVector&>(*icupp).cloneFull());
   StochVectorHandle clow_clone ( dynamic_cast<const StochVector&>(*bl).cloneFull());
   StochVectorHandle iclow_clone ( dynamic_cast<const StochVector&>(*iclow).cloneFull());

   sTree* tree_clone = stochNode; // todo

   sData* clone = new sData(tree_clone, c_clone, Q_clone, xlow_clone,
         ixlow_clone, nxlow, xupp_clone, ixupp_clone, nxupp, A_clone, bA_clone,
         C_clone, clow_clone, iclow_clone, mclow, cupp_clone, icupp_clone,
         mcupp);

   return clone;
}

void
sData::createChildren()
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

void
sData::destroyChildren()
{
   for( size_t it = 0; it < children.size(); it++ )
   {
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

void sData::init2LinksData(bool exploit2links)
{
   use2Links = exploit2links;

   if( !exploit2links )
   {
      linkStartBlocksA = std::vector<int>();
      linkStartBlocksC = std::vector<int>();
      return;
   }

   int myrank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   const StochGenMatrix& Astoch = dynamic_cast<const StochGenMatrix&>(*A);
   const StochGenMatrix& Cstoch = dynamic_cast<const StochGenMatrix&>(*C);

   int n2LinksEq = 0;
   int n2LinksIneq = 0;

   linkStartBlocksA = Astoch.get2LinkStartBlocks();
   linkStartBlocksC = Cstoch.get2LinkStartBlocks();

   linkStartBlockLengthsA = get2LinkLengthsVec(linkStartBlocksA, stochNode->children.size());
   linkStartBlockLengthsC = get2LinkLengthsVec(linkStartBlocksC, stochNode->children.size());

   printLinkConsStats();
   printLinkVarsStats();

   for( size_t i = 0; i < linkStartBlocksA.size(); ++i )
      if( linkStartBlocksA[i] >= 0 )
         n2LinksEq++;

   for( size_t i = 0; i < linkStartBlocksC.size(); ++i )
      if( linkStartBlocksC[i] >= 0 )
         n2LinksIneq++;

   assert(n2LinksEq == n2linksRows(linkStartBlockLengthsA));
   assert(n2LinksIneq == n2linksRows(linkStartBlockLengthsC));

   if( myrank == 0 )
   {
      std::cout << "number of equality 2-links: " << n2LinksEq << " (out of "
            << linkStartBlocksA.size() << " equalities) " << std::endl;
      std::cout << "number of inequality 2-links: " << n2LinksIneq << " (out of "
            << linkStartBlocksC.size() << " equalities) " << std::endl;

      std::cout << "ratio: "
            << (n2LinksEq + n2LinksIneq) / ((double) linkStartBlocksA.size() + linkStartBlocksC.size()) << std::endl;
   }

   if( n2LinksEq + n2LinksIneq / ((double) linkStartBlocksA.size() + linkStartBlocksC.size()) < min2LinksRatio )
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

double
sData::objectiveValue(QpGenVars * vars)
{
   StochVector& x = dynamic_cast<StochVector&>(*vars->x);
   OoqpVectorHandle temp(x.clone());

   this->getg(*temp);
   this->Qmult(1.0, *temp, 0.5, *vars->x);

   return temp->dotProductWith(*vars->x);
}

void
sData::createScaleFromQ()
{

   assert("Not implemented!" && 0);

   // Stuff the diagonal elements of Q into the vector "sc"
   this->getDiagonalOfQ(*sc);

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

void sData::printLinkVarsStats()
{
   int n = getLocalnx();

   std::vector<int> linkCount(n, 0);
   std::vector<int> linkCount0(n, 0);
   std::vector<int> linkCountLC(n, 0);

   StochGenMatrix& Astoch = dynamic_cast<StochGenMatrix&>(*A);
   StochGenMatrix& Cstoch = dynamic_cast<StochGenMatrix&>(*C);

   Astoch.updateKLinkVarsCount(linkCount);
   Cstoch.updateKLinkVarsCount(linkCount);

   Astoch.Bmat->getTranspose().updateNonEmptyRowsCount(linkCount0);
   Astoch.Bmat->deleteTransposed();
   Cstoch.Bmat->getTranspose().updateNonEmptyRowsCount(linkCount0);
   Cstoch.Bmat->deleteTransposed();

   if( Astoch.Blmat )
   {
      Astoch.Blmat->getTranspose().updateNonEmptyRowsCount(linkCountLC);
      Astoch.Blmat->deleteTransposed();
   }

   if( Cstoch.Blmat )
   {
      Cstoch.Blmat->getTranspose().updateNonEmptyRowsCount(linkCountLC);
      Cstoch.Blmat->deleteTransposed();
   }

   int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);

   if( rank == 0 )
   {
      std::vector<int> linkSizes(nLinkStats, 0);

      int count0 = 0;
      int countLC = 0;
      int count0LC = 0;

      for( int i = 0; i < n; i++ )
      {
         assert(linkCount[i] >= 0 && linkCount0[i] >= 0 && linkCountLC[i] >= 0);
         assert(linkCount0[i] <= 2 && linkCountLC[i] <= 2);

         if( linkCount[i] < nLinkStats )
            linkSizes[size_t(linkCount[i])]++;

         if( linkCount[i] == 0 && linkCountLC[i] == 0 && linkCount0[i] != 0 )
            count0++;

         if( linkCount[i] == 0 && linkCount0[i] == 0 && linkCountLC[i] != 0 )
            countLC++;

         if( linkCount[i] == 0 && (linkCount0[i] != 0 || linkCountLC[i] != 0) )
            count0LC++;
      }

      std::cout << "total link vars " << n << std::endl;

      for( int i = 0; i < nLinkStats; i++ )
         if( linkSizes[i] != 0 )
            std::cout << i << "-links " << linkSizes[i] << std::endl;

      std::cout << "Block0 exclusive vars " << count0 << std::endl;
      std::cout << "LC exclusive vars " << countLC << std::endl;
      std::cout << "Block0 or LC vars " << count0LC  << std::endl;
   }
}

void sData::printLinkConsStats()
{
   int myl = getLocalmyl();
   int mzl = getLocalmzl();

   int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);

   if( myl > 0 )
   {
      std::vector<int> linkCount(myl, 0);

      dynamic_cast<StochGenMatrix&>(*A).updateKLinkConsCount(linkCount);

      if( rank == 0 )
      {
         std::vector<int> linkSizes(nLinkStats, 0);

         for( int i = 0; i < myl; i++ )
            if( linkCount[i] < nLinkStats )
            {
               assert(linkCount[i] >= 0);
               linkSizes[size_t(linkCount[i])]++;
            }

         std::cout << "total equality Linking Constraints " << myl << std::endl;

         for( int i = 0; i < nLinkStats; i++ )
            if( linkSizes[i] != 0 )
               std::cout << i << "-links " << linkSizes[i] << std::endl;
      }
   }

   if( mzl > 0 )
   {
      std::vector<int> linkCount(mzl, 0);

      dynamic_cast<StochGenMatrix&>(*C).updateKLinkConsCount(linkCount);

      if( rank == 0 )
      {
         std::vector<int> linkSizes(nLinkStats, 0);

         for( int i = 0; i < mzl; i++ )
            if( linkCount[i] < nLinkStats )
            {
               assert(linkCount[i] >= 0);
               linkSizes[size_t(linkCount[i])]++;
            }

         std::cout << "total inequality Linking Constraints " << mzl << std::endl;

         for( int i = 0; i < nLinkStats; i++ )
            if( linkSizes[i] != 0 )
               std::cout << i << "-links " << linkSizes[i] << std::endl;
      }
   }
}

sData::~sData()
{
   for( size_t it = 0; it < children.size(); it++ )
      delete children[it];
}

int
sData::getLocalnx()
{
   StochSymMatrix& Qst = dynamic_cast<StochSymMatrix&>(*Q);
   return Qst.diag->size();
}

int
sData::getLocalmy()
{
   long long my, nx;
   StochGenMatrix& Ast = dynamic_cast<StochGenMatrix&>(*A);
   Ast.Bmat->getSize(my, nx);
   return my;
}

int
sData::getLocalmyl()
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

int
sData::getLocalmzl()
{
   long long mzl, nxl;
   StochGenMatrix& Cst = dynamic_cast<StochGenMatrix&>(*C);
   Cst.Blmat->getSize(mzl, nxl);
   return mzl;
}

int
sData::getLocalSizes(int& nx, int& my, int& mz, int& myl, int& mzl)
{
   long long nxloc, myloc, mzloc, mylloc, mzlloc;

   StochGenMatrix& Ast = dynamic_cast<StochGenMatrix&>(*A);
   Ast.Blmat->getSize(mylloc, nxloc);
   Ast.Bmat->getSize(myloc, nxloc);

   StochGenMatrix& Cst = dynamic_cast<StochGenMatrix&>(*C);
   Cst.Blmat->getSize(mzlloc, nxloc);
   Cst.Bmat->getSize(mzloc, nxloc);

   nx = nxloc;
   my = myloc;
   mz = mzloc;
   myl = mylloc;
   mzl = mzlloc;
   return 0;
}

int
sData::getLocalSizes(int& nx, int& my, int& mz)
{
   long long nxll, myll, mzll;

   StochGenMatrix& Ast = dynamic_cast<StochGenMatrix&>(*A);
   Ast.Bmat->getSize(myll, nxll);

   StochGenMatrix& Cst = dynamic_cast<StochGenMatrix&>(*C);
   Cst.Bmat->getSize(mzll, nxll);

   nx = nxll;
   my = myll;
   mz = mzll;
   return 0;
}

int
sData::getLocalNnz(int& nnzQ, int& nnzB, int& nnzD)
{
   StochSymMatrix& Qst = dynamic_cast<StochSymMatrix&>(*Q);
   StochGenMatrix& Ast = dynamic_cast<StochGenMatrix&>(*A);
   StochGenMatrix& Cst = dynamic_cast<StochGenMatrix&>(*C);

   nnzQ = Qst.diag->getStorageRef().len + Qst.border->getStorageRef().len;
   nnzB = Ast.Bmat->getStorageRef().len;
   nnzD = Cst.Bmat->getStorageRef().len;
   return 0;
}

int sData::getSchurCompMaxNnz()
{
   assert(children.size() > 0);

   const int n0 = getLocalnx();
   const int my = getLocalmy();
   const int myl = getLocalmyl();
   const int mzl = getLocalmzl();

#ifndef NDEBUG
   {
      int mB, nB;
      getLocalB().getSize(mB, nB);
      assert(mB == my  && nB == n0);
   }
#endif

   int nnz = 0;

   // sum up half of dense square
   nnz += nnzTriangular(n0);

   // add B_0 (or A_0, depending on notation)
   nnz += getLocalB().numberOfNonZeros();

   // add borders
   nnz += myl * n0;
   nnz += mzl * n0;

   // (empty) diagonal
   nnz += my;

   // add linking equality parts
   nnz += getSchurCompMaxNnz(linkStartBlocksA, linkStartBlockLengthsA);

   // add linking inequality parts
   nnz += getSchurCompMaxNnz(linkStartBlocksC, linkStartBlockLengthsC);

   // add linking mixed parts todo
   nnz += myl * mzl;

   return nnz;
}

SparseSymMatrix& sData::getLocalQ()
{
   StochSymMatrix& Qst = dynamic_cast<StochSymMatrix&>(*Q);
   return *Qst.diag;
}

SparseGenMatrix&
sData::getLocalCrossHessian()
{
   StochSymMatrix& Qst = dynamic_cast<StochSymMatrix&>(*Q);
   return *Qst.border;
}

// T_i x_0 + W_i x_i = b_i

// This is T_i
SparseGenMatrix&
sData::getLocalA()
{
   StochGenMatrix& Ast = dynamic_cast<StochGenMatrix&>(*A);
   return *Ast.Amat;
}

// This is W_i:
SparseGenMatrix&
sData::getLocalB()
{
   StochGenMatrix& Ast = dynamic_cast<StochGenMatrix&>(*A);
   return *Ast.Bmat;

}

// This is F_i (linking equality matrix):
SparseGenMatrix&
sData::getLocalF()
{
   StochGenMatrix& Ast = dynamic_cast<StochGenMatrix&>(*A);
   return *Ast.Blmat;
}

// low_i <= C_i x_0 + D_i x_i <= upp_i

// This is C_i
SparseGenMatrix&
sData::getLocalC()
{
   StochGenMatrix& Cst = dynamic_cast<StochGenMatrix&>(*C);
   return *Cst.Amat;
}

// This is D_i
SparseGenMatrix&
sData::getLocalD()
{
   StochGenMatrix& Cst = dynamic_cast<StochGenMatrix&>(*C);
   return *Cst.Bmat;
}

// This is G_i (linking inequality matrix):
SparseGenMatrix&
sData::getLocalG()
{
   StochGenMatrix& Cst = dynamic_cast<StochGenMatrix&>(*C);
   return *Cst.Blmat;
}

void
sData::cleanUpPresolvedData(const StochVector& rowNnzVecA, const StochVector& rowNnzVecC, const StochVector& colNnzVec)
{
   StochSymMatrix& Q_stoch = dynamic_cast<StochSymMatrix&>(*Q);

   // todo only works if Q is empty
   Q_stoch.deleteEmptyRowsCols(colNnzVec);

   // clean up equality system
   StochGenMatrix& A_stoch = dynamic_cast<StochGenMatrix&>(*A);
   StochVector& b_Astoch = dynamic_cast<StochVector&>(*bA);

   A_stoch.initStaticStorageFromDynamic(rowNnzVecA, colNnzVec);
   A_stoch.freeDynamicStorage();

   b_Astoch.removeEntries(rowNnzVecA);

   // clean up inequality system and x
   StochGenMatrix& C_stoch = dynamic_cast<StochGenMatrix&>(*C);
   StochVector& g_stoch = dynamic_cast<StochVector&>(*g);

   StochVector& blx_stoch = dynamic_cast<StochVector&>(*blx);
   StochVector& ixlow_stoch = dynamic_cast<StochVector&>(*ixlow);
   StochVector& bux_stoch = dynamic_cast<StochVector&>(*bux);
   StochVector& ixupp_stoch = dynamic_cast<StochVector&>(*ixupp);

   StochVector& bl_stoch = dynamic_cast<StochVector&>(*bl);
   StochVector& iclow_stoch = dynamic_cast<StochVector&>(*iclow);
   StochVector& bu_stoch = dynamic_cast<StochVector&>(*bu);
   StochVector& icupp_stoch = dynamic_cast<StochVector&>(*icupp);

   C_stoch.initStaticStorageFromDynamic(rowNnzVecC, colNnzVec);
   C_stoch.freeDynamicStorage();

   g_stoch.removeEntries(colNnzVec);

   blx_stoch.removeEntries(colNnzVec);
   ixlow_stoch.removeEntries(colNnzVec);
   bux_stoch.removeEntries(colNnzVec);
   ixupp_stoch.removeEntries(colNnzVec);

   bl_stoch.removeEntries(rowNnzVecC);
   iclow_stoch.removeEntries(rowNnzVecC);
   bu_stoch.removeEntries(rowNnzVecC);
   icupp_stoch.removeEntries(rowNnzVecC);

   assert(stochNode != NULL);

   // adapt sizes and tree
   sTreeCallbacks& callbackTree = dynamic_cast<sTreeCallbacks&>(*stochNode);

   callbackTree.initPresolvedData(Q_stoch, A_stoch, C_stoch, g_stoch, b_Astoch, iclow_stoch);
   callbackTree.switchToPresolvedData();

   long long dummy;
   nx = g_stoch.length();
   A_stoch.getSize( my, dummy );
   C_stoch.getSize( mz, dummy );

   nxlow = ixlow_stoch.numberOfNonzeros();
   nxupp = ixupp_stoch.numberOfNonzeros();
   mclow = iclow_stoch.numberOfNonzeros();
   mcupp = icupp_stoch.numberOfNonzeros();
}

void
sData::sync()
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
