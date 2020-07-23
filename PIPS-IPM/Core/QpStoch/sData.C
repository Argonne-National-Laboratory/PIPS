#include "sData.h"
#include "sTree.h"
#include "sTreeCallbacks.h"
#include "StochSymMatrix.h"
#include "StochGenMatrix.h"
#include "StochVector.h"
#include "SparseLinearAlgebraPackage.h"
#include "mpi.h"
#include <iostream>
#include "pipsport.h"

static
std::vector<unsigned int> getInversePermutation(const std::vector<unsigned int>& perm)
{
   size_t size = perm.size();
   std::vector<unsigned int> perm_inv(size, 0);

   for( size_t i = 0; i < size; i++ )
      perm_inv[perm[i]] = i;

   return perm_inv;
}

static inline
bool blockIsInRange(int block, int blocksStart, int blocksEnd)
{
   return ((block >= (blocksStart - 1) && block < blocksEnd) || block == -1);
}

static inline
int nnzTriangular(int size)
{
   assert(size >= 0);
   return ((1 + size) * size) / 2;
}

static
void appendRowDense(int start, int end, int& nnz, int* jcolM)
{
   assert(jcolM);
   assert(nnz >= 0);
   assert(start >= 0 && end >= start);

   for( int i = start; i < end; i++ )
      jcolM[nnz++] = i;
}

static
void appendRowSparse(int startColIdx, int endColIdx, int colOffset, const int* jcolM_append, int& nnz, int* jcolM)
{
   assert(jcolM);
   assert(nnz >= 0);
   assert(startColIdx >= 0 && startColIdx <= endColIdx);

   for( int c = startColIdx; c < endColIdx; c++ )
      jcolM[nnz++] = colOffset + jcolM_append[c];
}


static
int appendDiagBlocks(const std::vector<int>& linkStartBlockId, const std::vector<int>& linkStartBlockLengths, int borderstart, int bordersize, int rowSC,
                                          int rowBlock, int& blockStartrow, int& nnz, int* jcolM)
{
   assert(rowBlock >= blockStartrow && blockStartrow >= 0 && borderstart >= 0 && bordersize >= 0 && nnz >= 0);

   const int block = linkStartBlockId[rowBlock];
   const int currlength = (block >= 0) ? linkStartBlockLengths[block] : bordersize;

   assert(currlength >= 1);

   // add diagonal block (possibly up to the border)

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

static
void appendDiagBlocksDist(const std::vector<int>& linkStartBlockId, const std::vector<int>& linkStartBlockLengths, int borderstart, int bordersize, int rowSC,
                                          int rowBlock, int blocksStart, int blocksEnd, int& blockStartrow, int& nnz, int* jcolM)
{
   assert(rowBlock >= blockStartrow && blockStartrow >= 0 && borderstart >= 0 && bordersize >= 0 && nnz >= 0);

   const int block = linkStartBlockId[rowBlock];
   const int lastBlock = blocksEnd - 1;

   if( blockIsInRange(block, blocksStart, blocksEnd) )
   {
      const int currlength = (block >= 0) ? linkStartBlockLengths[block] : bordersize;
      const int rownnz = currlength - (rowBlock - blockStartrow);

      assert(currlength >= 1);

      // add diagonal block (possibly up to the border)
      for( int i = 0; i < rownnz; ++i )
         jcolM[nnz++] = rowSC + i;

      // with off-diagonal blocks? (at sparse part)
      if( block >= 0 )
      {
         // add right off-diagonal block and border part

         const int nextlength = (block == lastBlock) ? 0 : linkStartBlockLengths[block + 1];

         assert(nextlength >= 0);
         assert(block != int(linkStartBlockLengths.size()) - 1 || nextlength == 0);

         for( int i = rownnz; i < rownnz + nextlength; ++i )
            jcolM[nnz++] = rowSC + i;

         for( int i = borderstart; i < borderstart + bordersize; ++i )
            jcolM[nnz++] = i;

         // last row of current block?
         if( rowBlock + 1 == blockStartrow + currlength )
            blockStartrow = rowBlock + 1;
      }
   }

   // at sparse part?
   if( block >= 0 )
   {
      const int currlength = linkStartBlockLengths[block];
      assert(currlength >= 1);

      // last row of current block?
      if( rowBlock + 1 == blockStartrow + currlength )
         blockStartrow = rowBlock + 1;
   }
}


static
int appendMixedBlocks(const std::vector<int>& linkStartBlockId_Left,
      const std::vector<int>& linkStartBlockId_Right,
      const std::vector<int>& linkStartBlockLengths_Left,
      const std::vector<int>& linkStartBlockLengths_Right,
      int colStartIdxSC, int bordersize_cols,
      int rowIdx, int& colIdxOffset, int& rowBlockStartIdx, int& nnz, int* jcolM)
{
   assert(rowIdx >= rowBlockStartIdx && rowBlockStartIdx >= 0 && colStartIdxSC >= 0 && nnz >= 0);
   assert(bordersize_cols >= 0 && colIdxOffset >= 0);
   assert(linkStartBlockLengths_Left.size() == linkStartBlockLengths_Right.size());

   const int nCols = int(linkStartBlockId_Right.size());
   const int block = linkStartBlockId_Left[rowIdx];

   assert(nCols >= bordersize_cols);

   int rownnz;

   // sparse row?
   if( block >= 0 )
   {
      const int length_Right = linkStartBlockLengths_Right[block];
      const int colStartIdxBorderSC = colStartIdxSC + nCols - bordersize_cols;
      int colStartIdx = colStartIdxSC + colIdxOffset;
      rownnz = 0;

      assert(length_Right >= 0);

      // 1) left off-diagonal block (not for first block)
      if( block >= 1 )
      {
         const int prevlength_Right = linkStartBlockLengths_Right[block - 1];
         assert(prevlength_Right >= 0);

         for( int i = 0; i < prevlength_Right; ++i )
            jcolM[nnz++] = (colStartIdx + i);

         rownnz += prevlength_Right;
         colStartIdx += prevlength_Right;
      }

      // 2) diagonal block
      for( int i = 0; i < length_Right; ++i )
         jcolM[nnz++] = (colStartIdx + i);

      rownnz += length_Right;
      colStartIdx += length_Right;

      // 3) right off-diagonal block (not for last block)
      if( int(linkStartBlockLengths_Left.size()) != block + 1 )
      {
         const int nextlength_Right = linkStartBlockLengths_Right[block + 1];

         for( int i = 0; i < nextlength_Right; ++i )
            jcolM[nnz++] = (colStartIdx + i);

         rownnz += nextlength_Right;
         colStartIdx += nextlength_Right;
      }

      assert(colStartIdx <= colStartIdxBorderSC);

      // 4) right border
      for( int i = 0; i < bordersize_cols; ++i )
         jcolM[nnz++] = colStartIdxBorderSC + i;

      rownnz += bordersize_cols;

      // last row of current block?
      if( rowIdx + 1 == rowBlockStartIdx + linkStartBlockLengths_Left[block] )
      {
         rowBlockStartIdx = rowIdx + 1;

         if( block >= 1 )
            colIdxOffset += linkStartBlockLengths_Right[block - 1];
      }
      else
      {
         assert(block == linkStartBlockId_Left[rowIdx + 1]);
      }
   }
   else
   {
      // append fully dense row
      for( int i = 0; i < nCols; ++i )
         jcolM[nnz++] = colStartIdxSC + i;

      rownnz = nCols;
   }

   return rownnz;
}


static
void appendMixedBlocksDist(const std::vector<int>& linkStartBlockId_Left,
      const std::vector<int>& linkStartBlockId_Right,
      const std::vector<int>& linkStartBlockLengths_Left,
      const std::vector<int>& linkStartBlockLengths_Right,
      int colStartIdxSC, int bordersize_cols,
      int rowIdx, int blocksStart, int blocksEnd,
      int& colIdxOffset, int& rowBlockStartIdx, int& nnz, int* jcolM)
{
   assert(rowIdx >= rowBlockStartIdx && rowBlockStartIdx >= 0 && colStartIdxSC >= 0 && nnz >= 0);
   assert(bordersize_cols >= 0 && colIdxOffset >= 0);
   assert(linkStartBlockLengths_Left.size() == linkStartBlockLengths_Right.size());

   const int block = linkStartBlockId_Left[rowIdx];
   const bool blockInRange = blockIsInRange(block, blocksStart, blocksEnd);
   const int nCols = int(linkStartBlockId_Right.size());

   assert(nCols >= bordersize_cols);

   // sparse row?
   if( block >= 0 )
   {
      if( blockInRange )
      {
         const int length_Right = linkStartBlockLengths_Right[block];
         const int colStartIdxBorderSC = colStartIdxSC + nCols - bordersize_cols;
         int colStartIdx = colStartIdxSC + colIdxOffset;

         assert(length_Right >= 0);

         // 1) left off-diagonal block (not for first block)
         if( block >= 1 && block != blocksStart - 1 )
         {
            const int prevlength_Right = linkStartBlockLengths_Right[block - 1];
            assert(prevlength_Right >= 0);

            for( int i = 0; i < prevlength_Right; ++i )
               jcolM[nnz++] = (colStartIdx + i);
         }

         if( block >= 1 )
         {
            const int prevlength_Right = linkStartBlockLengths_Right[block - 1];
            assert(prevlength_Right >= 0);

            colStartIdx += prevlength_Right;
         }

         // 2) diagonal block
         for( int i = 0; i < length_Right; ++i )
            jcolM[nnz++] = (colStartIdx + i);

         colStartIdx += length_Right;

         // 3) right off-diagonal block (not for last block)
         if( block != blocksEnd - 1 )
         {
            const int nextlength_Right = linkStartBlockLengths_Right[block + 1];

            for( int i = 0; i < nextlength_Right; ++i )
               jcolM[nnz++] = (colStartIdx + i);

            colStartIdx += nextlength_Right;
         }

         assert(colStartIdx <= colStartIdxBorderSC);

         // 4) right border
         for( int i = 0; i < bordersize_cols; ++i )
            jcolM[nnz++] = colStartIdxBorderSC + i;
      }

      // last row of current block?
      if( rowIdx + 1 == rowBlockStartIdx + linkStartBlockLengths_Left[block] )
      {
         rowBlockStartIdx = rowIdx + 1;

         if( block >= 1 )
            colIdxOffset += linkStartBlockLengths_Right[block - 1];
      }
      else
      {
         assert(block == linkStartBlockId_Left[rowIdx + 1]);
      }
   }
   else if( blockInRange )
   {
      assert(block == -1);

      // append dense row, but skip entries not in range
      for( int i = 0; i < nCols; ++i )
      {
         const int blockRight = linkStartBlockId_Right[i];
         const bool blockRightInRange = blockIsInRange(blockRight, blocksStart, blocksEnd);

         if( blockRightInRange )
            jcolM[nnz++] = colStartIdxSC + i;
      }
   }
}


void sData::getSCrangeMarkers(int blocksStart, int blocksEnd, int& local2linksStartEq, int& local2linksEndEq,
      int& local2linksStartIneq, int& local2linksEndIneq)
{
   const int blocksStartReal = (blocksStart > 0) ? (blocksStart - 1) : blocksStart;
   const int nx0 = getLocalnx();
   const int my0 = getLocalmy();
   const int myl = getLocalmyl();
   local2linksStartEq = nx0 + my0;
   local2linksStartIneq = nx0 + my0 + myl;

   for( int block = 0; block < blocksStartReal; ++block )
   {
      const int lengthEq = linkStartBlockLengthsA[block];
      const int lengthIneq = linkStartBlockLengthsC[block];

      assert(lengthEq >= 0 && lengthIneq >= 0);

      local2linksStartEq += lengthEq;
      local2linksStartIneq += lengthIneq;
   }

   local2linksEndEq = local2linksStartEq + getSCdiagBlocksNRows(linkStartBlockLengthsA, blocksStart, blocksEnd);
   local2linksEndIneq = local2linksStartIneq + getSCdiagBlocksNRows(linkStartBlockLengthsC, blocksStart, blocksEnd);
}

void sData::getSCrangeMarkersMy(int blocksStart, int blocksEnd, int& local2linksStartEq, int& local2linksEndEq,
      int& local2linksStartIneq, int& local2linksEndIneq)
{
   const int nx0 = getLocalnx();
   const int my0 = getLocalmy();
   const int myl = getLocalmyl();
   local2linksStartEq = nx0 + my0;
   local2linksStartIneq = nx0 + my0 + myl;

   for( int block = 0; block < blocksStart; ++block )
   {
      const int lengthEq = linkStartBlockLengthsA[block];
      const int lengthIneq = linkStartBlockLengthsC[block];

      assert(lengthEq >= 0 && lengthIneq >= 0);

      local2linksStartEq += lengthEq;
      local2linksStartIneq += lengthIneq;
   }

   local2linksEndEq = local2linksStartEq + getSCdiagBlocksNRowsMy(linkStartBlockLengthsA, blocksStart, blocksEnd);
   local2linksEndIneq = local2linksStartIneq + getSCdiagBlocksNRowsMy(linkStartBlockLengthsC, blocksStart, blocksEnd);
}


int sData::getSCdiagBlocksNRows(const std::vector<int>& linkStartBlockLengths,
      int blocksStart, int blocksEnd)
{
   assert(blocksStart >= 0 && blocksStart < blocksEnd);
   assert(blocksEnd <= int(linkStartBlockLengths.size()));

   int nRowsRange = 0;
   const int blocksStartReal = (blocksStart > 0) ? (blocksStart - 1) : blocksStart;

   // main loop, going over specified 2-link blocks
   for( int block = blocksStartReal; block < blocksEnd; ++block )
   {
      const int length = linkStartBlockLengths[block];
      assert(length >= 0);

      nRowsRange += length;
   }
   return nRowsRange;
}


int sData::getSCdiagBlocksNRowsMy(const std::vector<int>& linkStartBlockLengths,
      int blocksStart, int blocksEnd)
{
   assert(blocksStart >= 0 && blocksStart < blocksEnd);
   assert(blocksEnd <= int(linkStartBlockLengths.size()));

   int nRowsRange = 0;

   // main loop, going over specified 2-link blocks
   for( int block = blocksStart; block < blocksEnd; ++block )
   {
      const int length = linkStartBlockLengths[block];
      assert(length >= 0);

      nRowsRange += length;
   }
   return nRowsRange;
}

int sData::getSCdiagBlocksNRows(const std::vector<int>& linkStartBlockLengths)
{
   return (getSCdiagBlocksNRows(linkStartBlockLengths, 0, int(linkStartBlockLengths.size())));
}

int sData::getSCdiagBlocksMaxNnz(size_t nRows, const std::vector<int>& linkStartBlockLengths)
{
   const size_t nBlocks = linkStartBlockLengths.size();
   size_t nRowsSparse = 0;

   int nnz = 0;

   // main loop, going over all 2-link blocks
   for( size_t block = 0; block < nBlocks; ++block )
   {
      if( linkStartBlockLengths[block] == 0 )
         continue;

      const int length = linkStartBlockLengths[block];
      const int nextlength = linkStartBlockLengths[block + 1];

      assert(length > 0);
      assert(nextlength >= 0);
      assert(block != linkStartBlockLengths.size() - 2 || nextlength == 0);

      nRowsSparse += size_t(length);

      // diagonal block
      nnz += nnzTriangular(length);

      // (one) off-diagonal block
      nnz += length * nextlength;
   }

   // any rows left?
   if( nRowsSparse < nRows )
   {
      const size_t nRowsDense = nRows - nRowsSparse;
      nnz += nnzTriangular(nRowsDense) + nRowsDense * nRowsSparse;
   }

   return nnz;
}


int sData::getSCdiagBlocksMaxNnzDist(size_t nRows, const std::vector<int>& linkStartBlockLengths, int blocksStart, int blocksEnd)
{
#ifndef NDEBUG
   const int nblocks = int(linkStartBlockLengths.size());
#endif
   assert(blocksStart >= 0);
   assert(blocksStart < blocksEnd);
   assert(blocksEnd <= nblocks);
   assert(nblocks >= 2);

   const int nRowsSparse = getSCdiagBlocksNRows(linkStartBlockLengths);
   const int nRowsSparseRange = getSCdiagBlocksNRows(linkStartBlockLengths, blocksStart, blocksEnd);

   int nnz = 0;

   // main loop, going over specified 2-link blocks
   for( int block = blocksStart; block < blocksEnd; ++block )
   {
      const int length = linkStartBlockLengths[block];

      if( length == 0 )
         continue;

      const int prevlength = (block == 0) ? 0 : linkStartBlockLengths[block - 1];

      assert(length > 0);
      assert(prevlength >= 0);
      assert(block != nblocks - 1); // length should be 0 for last block

      // diagonal block
      nnz += nnzTriangular(length);

      // above off-diagonal block
      nnz += prevlength * length;
   }

   if( blocksStart > 0 )
   {
      const int prevlength = linkStartBlockLengths[blocksStart - 1];

      nnz += nnzTriangular(prevlength);
   }

   // any rows left?
   if( nRowsSparse < int(nRows) )
   {
      const int nRowsDense = int(nRows) - nRowsSparse;
      nnz += nnzTriangular(nRowsDense);
      nnz += nRowsDense * nRowsSparseRange;
   }

   return nnz;
}

int sData::getSCmixedBlocksMaxNnz(size_t nRows, size_t nCols,
      const std::vector<int>& linkStartBlockLength_Left,
      const std::vector<int>& linkStartBlockLength_Right)
{
   assert(linkStartBlockLength_Left.size() == linkStartBlockLength_Right.size());

   const size_t nBlocks = linkStartBlockLength_Left.size();
   size_t nRowsSparse = 0;
   size_t nColsSparse = 0;

   int nnz = 0;

   // main loop, going over all 2-link blocks
   for( size_t block = 0; block < nBlocks; ++block )
   {
      const int length_Left = linkStartBlockLength_Left[block];
      const int length_Right = linkStartBlockLength_Right[block];
      assert(length_Left >= 0 && length_Right >= 0);

      nRowsSparse += size_t(length_Left);
      nColsSparse += size_t(length_Right);

      // diagonal block
      nnz += length_Left * length_Right;

      if( block == 0 )
         continue;

      const int prevlength_Left = linkStartBlockLength_Left[block - 1];
      const int prevlength_Right = linkStartBlockLength_Right[block - 1];

      assert(prevlength_Left >= 0 && prevlength_Right >= 0);

      // left off-diagonal block
      nnz += length_Left * prevlength_Right;

      // upper off-diagonal block
      nnz += prevlength_Left * length_Right;
   }

   // dense border?
   if( nRowsSparse < nRows || nColsSparse < nCols )
   {
      assert(nRowsSparse <= nRows && nColsSparse <= nCols);

      const size_t nRowsDense = nRows - nRowsSparse;
      const size_t nColsDense = nCols - nColsSparse;

      nnz += nRowsDense * nColsSparse; // lower border part without right border
      nnz += nColsDense * nRows;       // complete right border
   }

   return nnz;
}


int sData::getSCmixedBlocksMaxNnzDist(size_t nRows, size_t nCols,
      const std::vector<int>& linkStartBlockLength_Left,
      const std::vector<int>& linkStartBlockLength_Right,
      int blocksStart, int blocksEnd)
{
   assert(linkStartBlockLength_Left.size() == linkStartBlockLength_Right.size());
   assert(blocksStart >= 0);
   assert(blocksStart < blocksEnd);

#ifndef NDEBUG
   const int nBlocks = int(linkStartBlockLength_Left.size());
#endif
   const int blockLast = blocksEnd - 1;
   const int blocksStartReal =  (blocksStart > 0) ? (blocksStart - 1) : 0;
   const int nRowsSparse = getSCdiagBlocksNRows(linkStartBlockLength_Left);
   const int nColsSparse = getSCdiagBlocksNRows(linkStartBlockLength_Right);
   const int nRowsSparseRange = getSCdiagBlocksNRows(linkStartBlockLength_Left, blocksStart, blocksEnd);
   const int nColsSparseRange = getSCdiagBlocksNRows(linkStartBlockLength_Right, blocksStart, blocksEnd);

   assert(nBlocks > 1 && blocksEnd <= nBlocks);
   assert(linkStartBlockLength_Left[nBlocks - 1] == 0 && linkStartBlockLength_Right[nBlocks - 1] == 0);

   int nnz = 0;

   // main loop, going over all 2-link blocks
   for( int block = blocksStartReal; block < blocksEnd; ++block )
   {
      const int length_Left = linkStartBlockLength_Left[block];
      const int length_Right = linkStartBlockLength_Right[block];

      // left off-diagonal block
      if( block != blocksStartReal )
      {
         assert(block >= 1);
         assert(linkStartBlockLength_Right[block - 1] >= 0);

         nnz += length_Left * length_Right;
      }

      // diagonal block
      nnz += length_Left * length_Right;

      // right off-diagonal block
      if( block != blockLast )
      {
         assert(block < nBlocks - 1);

         const int nextlength_Right = linkStartBlockLength_Right[block + 1];

         assert(nextlength_Right >= 0);

         nnz += length_Left * nextlength_Right;
      }
   }

   // dense (right or lower) border?
   if( nRowsSparse < int(nRows) || nColsSparse < int(nCols) )
   {
      assert(nRowsSparse <= int(nRows) && nColsSparse <= int(nCols));

      const int nRowsDense = int(nRows) - nRowsSparse;
      const int nColsDense = int(nCols) - nColsSparse;

      nnz += nRowsDense * nColsSparseRange;  // lower left border part (without right border)
      nnz += nRowsSparseRange * nColsDense;  // upper right border
      nnz += nRowsDense * nColsDense;        // lower right border

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

std::vector<int> sData::get2LinkLengthsVec(const std::vector<int>& linkStartBlockId, size_t nBlocks)
{
   std::vector<int> linkStartBlockLengths(nBlocks, 0);

   const size_t nlinks = linkStartBlockId.size();

   for( size_t i = 0; i < nlinks; i++ )
   {
      const int block = linkStartBlockId[i];

      if( block >= 0 )
      {
         assert(size_t(block) < nBlocks);
         linkStartBlockLengths[block]++;
      }
   }

   assert(nBlocks == 0 || linkStartBlockLengths[nBlocks - 1] == 0);

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

   const int nx0NonZero = nx0 - n0LinkVars;
   int nnzcount = 0;

   assert(nx0NonZero >= 0);

   // dense square block, B_0^T, and dense border blocks todo: add space for CDCt
   for( int i = 0; i < nx0NonZero; ++i )
   {
      const int blength = startRowBtrans[i + 1] - startRowBtrans[i];
      assert(blength >= 0);

      krowM[i + 1] = krowM[i] + (nx0 - i) + blength + myl + mzl;

      appendRowDense(i, nx0, nnzcount, jcolM);

      appendRowSparse(startRowBtrans[i], startRowBtrans[i + 1], nx0, colidxBtrans, nnzcount, jcolM);

      appendRowDense(nx0 + my0, nx0 + my0 + myl + mzl, nnzcount, jcolM);

      assert(nnzcount == krowM[i + 1]);
   }

   // dense square block and rest of B_0, F_0^T, G_0^T
   for( int i = nx0NonZero; i < nx0; ++i )
   {
      appendRowDense(i, nx0, nnzcount, jcolM);

      appendRowSparse(startRowBtrans[i], startRowBtrans[i + 1], nx0, colidxBtrans, nnzcount, jcolM);

      if( myl > 0 )
      {
         SparseGenMatrix& Ft = getLocalF().getTranspose();
         const int* startRowFtrans = Ft.krowM();
         const int* colidxFtrans = Ft.jcolM();

         appendRowSparse(startRowFtrans[i], startRowFtrans[i + 1], nx0 + my0, colidxFtrans, nnzcount, jcolM);
      }

      if( mzl > 0 )
      {
         SparseGenMatrix& Gt = getLocalG().getTranspose();
         const int* startRowGtrans = Gt.krowM();
         const int* colidxGtrans = Gt.jcolM();

         appendRowSparse(startRowGtrans[i], startRowGtrans[i + 1], nx0 + my0 + myl, colidxGtrans, nnzcount, jcolM);
      }

      krowM[i + 1] = nnzcount;
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
   const int n2linksRowsEq = n2linkRowsEq();
   const int bordersizeEq = linkStartBlockIdA.size() - n2linksRowsEq;
   const int borderstartEq = nx0 + my0 + n2linksRowsEq;
   const int n2linksRowsIneq = n2linkRowsIneq();
   const int bordersizeIneq = linkStartBlockIdC.size() - n2linksRowsIneq;
   const int borderstartIneq = nx0 + my0 + myl + n2linksRowsIneq;

   assert(bordersizeEq >= 0 && n2linksRowsEq <= myl);
   assert(bordersizeIneq >= 0 && n2linksRowsIneq <= mzl);

   for( int i = nx0 + my0, j = 0, colIdxOffset = 0, blockStartrowMix = 0; i < nx0 + my0 + myl; ++i, ++j )
   {
      int blockrownnz = appendDiagBlocks(linkStartBlockIdA, linkStartBlockLengthsA, borderstartEq, bordersizeEq, i, j,
            blockStartrow, nnzcount, jcolM);

      blockrownnz += appendMixedBlocks(linkStartBlockIdA, linkStartBlockIdC, linkStartBlockLengthsA, linkStartBlockLengthsC,
            (nx0 + my0 + myl), bordersizeIneq, j, colIdxOffset, blockStartrowMix, nnzcount, jcolM);

      assert(blockStartrowMix == blockStartrow);

      krowM[i + 1] = krowM[i] + blockrownnz;
   }

   // inequality linking: dense border block and sparse diagonal blocks
   blockStartrow = 0;

   for( int i = nx0 + my0 + myl, j = 0; i < nx0 + my0 + myl + mzl; ++i, ++j )
   {
       const int blockrownnz = appendDiagBlocks(linkStartBlockIdC, linkStartBlockLengthsC, borderstartIneq, bordersizeIneq, i, j, blockStartrow, nnzcount, jcolM);

       krowM[i + 1] = krowM[i] + blockrownnz;
   }

   assert(nnzcount == nnz);

   return (new SparseSymMatrix(sizeSC, nnz, krowM, jcolM, M, 1, false));
}


SparseSymMatrix* sData::createSchurCompSymbSparseUpperDist(int blocksStart, int blocksEnd)
{
   assert(children.size() > 0);

   const int nx0 = getLocalnx();
   const int my0 = getLocalmy();
   const int myl = getLocalmyl();
   const int mzl = getLocalmzl();
   const int mylLocal = myl - getSCdiagBlocksNRows(linkStartBlockLengthsA)
      + getSCdiagBlocksNRows(linkStartBlockLengthsA, blocksStart, blocksEnd);
   const int mzlLocal = mzl - getSCdiagBlocksNRows(linkStartBlockLengthsC)
      + getSCdiagBlocksNRows(linkStartBlockLengthsC, blocksStart, blocksEnd);
   const int sizeSC = nx0 + my0 + myl + mzl;
   const int nnz = getSchurCompMaxNnzDist(blocksStart, blocksEnd);

   assert(getSchurCompMaxNnzDist(0, linkStartBlockLengthsA.size()) == getSchurCompMaxNnz());
   assert(blocksStart >= 0 && blocksStart < blocksEnd);
   assert(nnz > 0);
   assert(myl >= 0 && mzl >= 0);

   int* const krowM = new int[sizeSC + 1];
   int* const jcolM = new int[nnz];
   double* const M = new double[nnz];

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

   const int nx0NonZero = nx0 - n0LinkVars;
   const int n2linksRowsEq = n2linkRowsEq();
   const int bordersizeEq = linkStartBlockIdA.size() - n2linksRowsEq;
   const int borderstartEq = nx0 + my0 + n2linksRowsEq;
   const int n2linksRowsIneq = n2linkRowsIneq();
   const int bordersizeIneq = linkStartBlockIdC.size() - n2linksRowsIneq;
   const int borderstartIneq = nx0 + my0 + myl + n2linksRowsIneq;
   int local2linksStartEq;
   int local2linksEndEq;
   int local2linksStartIneq;
   int local2linksEndIneq;

   this->getSCrangeMarkers(blocksStart, blocksEnd, local2linksStartEq, local2linksEndEq,
         local2linksStartIneq, local2linksEndIneq);

   assert(nx0NonZero >= 0);
   assert(bordersizeEq >= 0 && n2linksRowsEq <= myl);
   assert(bordersizeIneq >= 0 && n2linksRowsIneq <= mzl);
   assert(local2linksStartEq >= nx0 + my0 && local2linksEndEq <= borderstartEq);
   assert(local2linksStartIneq >= nx0 + my0 + myl && local2linksEndIneq <= borderstartIneq);

   int nnzcount = 0;

   // dense square block, B_0^T, and dense border blocks todo: add space for CDCt
   for( int i = 0; i < nx0NonZero; ++i )
   {
      const int blength = startRowBtrans[i + 1] - startRowBtrans[i];
      assert(blength >= 0);

      krowM[i + 1] = krowM[i] + (nx0 - i) + blength + mylLocal + mzlLocal;

      appendRowDense(i, nx0, nnzcount, jcolM);
      appendRowSparse(startRowBtrans[i], startRowBtrans[i + 1], nx0, colidxBtrans, nnzcount, jcolM);

      appendRowDense(local2linksStartEq, local2linksEndEq, nnzcount, jcolM);
      appendRowDense(borderstartEq, borderstartEq + bordersizeEq, nnzcount, jcolM);
      appendRowDense(local2linksStartIneq, local2linksEndIneq, nnzcount, jcolM);
      appendRowDense(borderstartIneq, borderstartIneq + bordersizeIneq, nnzcount, jcolM);

      assert(nnzcount == krowM[i + 1]);
   }

   // dense square block and rest of B_0, F_0^T, G_0^T
   for( int i = nx0NonZero; i < nx0; ++i )
   {
      appendRowDense(i, nx0, nnzcount, jcolM);

      appendRowSparse(startRowBtrans[i], startRowBtrans[i + 1], nx0, colidxBtrans, nnzcount, jcolM);

      if( myl > 0 )
      {
         SparseGenMatrix& Ft = getLocalF().getTranspose();
         const int* startRowFtrans = Ft.krowM();
         const int* colidxFtrans = Ft.jcolM();

         appendRowSparse(startRowFtrans[i], startRowFtrans[i + 1], nx0 + my0, colidxFtrans, nnzcount, jcolM);
      }

      if( mzl > 0 )
      {
         SparseGenMatrix& Gt = getLocalG().getTranspose();
         const int* startRowGtrans = Gt.krowM();
         const int* colidxGtrans = Gt.jcolM();

         appendRowSparse(startRowGtrans[i], startRowGtrans[i + 1], nx0 + my0 + myl, colidxGtrans, nnzcount, jcolM);
      }

      krowM[i + 1] = nnzcount;
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

   for( int i = nx0 + my0, j = 0, colIdxOffset = 0, blockStartrowMix = 0; i < nx0 + my0 + myl; ++i, ++j )
   {
      appendDiagBlocksDist(linkStartBlockIdA, linkStartBlockLengthsA, borderstartEq, bordersizeEq, i, j,
            blocksStart, blocksEnd, blockStartrow, nnzcount, jcolM);

      appendMixedBlocksDist(linkStartBlockIdA, linkStartBlockIdC, linkStartBlockLengthsA, linkStartBlockLengthsC, (nx0 + my0 + myl),
            bordersizeIneq, j, blocksStart, blocksEnd, colIdxOffset, blockStartrowMix, nnzcount, jcolM);

      assert(blockStartrowMix == blockStartrow);

      krowM[i + 1] = nnzcount;
   }

   // inequality linking: dense border block and sparse diagonal blocks

   blockStartrow = 0;

   for( int i = nx0 + my0 + myl, j = 0; i < nx0 + my0 + myl + mzl; ++i, ++j )
   {
      appendDiagBlocksDist(linkStartBlockIdC, linkStartBlockLengthsC, borderstartIneq, bordersizeIneq, i, j,
            blocksStart, blocksEnd, blockStartrow, nnzcount, jcolM);

      krowM[i + 1] = nnzcount;
   }

   assert(nnzcount == nnz);

   this->initDistMarker(blocksStart, blocksEnd);

   return (new SparseSymMatrix(sizeSC, nnzcount, krowM, jcolM, M, 1, false));
}

std::vector<unsigned int> sData::get0VarsRightPermutation(const std::vector<int>& linkVarsNnzCount)
{
   const int size = int(linkVarsNnzCount.size());

   if( size == 0 )
      return std::vector<unsigned int>();

   std::vector<unsigned int> permvec(size, 0);

   int count = 0;
   int backCount = size - 1;
   for( int i = 0; i < size; ++i )
   {
      assert(count <= backCount);
      assert(linkVarsNnzCount[i] >= 0);

      if( linkVarsNnzCount[i] != 0 )
         permvec[count++] = i;
      else
         permvec[backCount--] = i;
   }

   assert(count == backCount + 1);

   return permvec;
}

std::vector<unsigned int> sData::getAscending2LinkPermutation(std::vector<int>& linkStartBlockId, size_t nBlocks)
{
   const size_t size = linkStartBlockId.size();

   if( size == 0 )
      return std::vector<unsigned int>();

   std::vector<unsigned int> permvec(size, 0);
   std::vector<int> w(nBlocks + 1, 0);

   for( size_t i = 0; i < size; ++i )
   {
      assert(linkStartBlockId[i] >= - 1 && linkStartBlockId[i] < int(nBlocks));
      w[linkStartBlockId[i] + 1]++;
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
      const int startBlock = (linkStartBlockId[i] >= 0) ? linkStartBlockId[i] : int(nBlocks);

      assert(w[startBlock] <= int(size));
      assert(permvec[w[startBlock]] == 0);

      permvec[w[startBlock]] = i;
      w[startBlock]++;
   }

#ifndef NDEBUG
     for( size_t i = 1; i < permvec.size(); i++ )
        assert(linkStartBlockId[permvec[i]] == - 1 || linkStartBlockId[permvec[i - 1]] <=  linkStartBlockId[permvec[i]]);
#endif

   // permute linkStartBlockId
   std::vector<int> tmpvec(size);

   for( size_t i = 0; i < size; ++i )
      tmpvec[i] = linkStartBlockId[permvec[i]];

   linkStartBlockId = tmpvec;

   return permvec;
}

sData::sData(sTree* tree)
//  : QpGenData(nullptr,nullptr,nullptr,nullptr,nullptr,nullptr,nullptr,nullptr,nullptr,nullptr,nullptr,nullptr,nullptr,nullptr)
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

   useLinkStructure = false;
   n0LinkVars = 0;
}

sData::sData(sTree* tree_, OoqpVector * c_in, SymMatrix * Q_in,
        OoqpVector * xlow_in, OoqpVector * ixlow_in, long long nxlow_,
        OoqpVector * xupp_in, OoqpVector * ixupp_in, long long nxupp_,
        GenMatrix  * A_in, OoqpVector * bA_in,
        GenMatrix  * C_in,
        OoqpVector * clow_in, OoqpVector * iclow_in, long long mclow_,
        OoqpVector * cupp_in, OoqpVector * icupp_in, long long mcupp_
        )
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

  useLinkStructure = false;
  n0LinkVars = 0;
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

/** Write the LP in MPS format. Only works if not distributed. */
void sData::writeMPSformat(ostream& out)
{
   // Note: only write the inequalities that have a finite rhs
   // (because no specified rhs of a row implies rhs=0).
   // Also, variable coefficients with indices in inequalitites with
   // inifnite rhs are not written because these rows do not appear in the MPs model.

   int world_size;
   MPI_Comm_size(MPI_COMM_WORLD, &world_size);

   if( world_size > 1 )
   {
      cout<<"MPS format writer only available using one Process!"<<endl;
      return;
   }
   cout<<"Writing MPS format..."<<endl;

   out <<  "NAME PIPS_to_MPS " << endl;
   out << "ROWS" <<endl;
   out << " N COST" <<endl;

   // write all row names and if they are E, L or G
   (*A).writeMPSformatRows(out, 0, nullptr);
   (*C).writeMPSformatRows(out, 1, icupp);
   (*C).writeMPSformatRows(out, 2, iclow);

   // write all variable names
   out <<  "COLUMNS " << endl;
   writeMPSColumns(out);

   // write all rhs / lhs
   out <<  "RHS " << endl;

   (*bA).writeMPSformatRhs(out, 0, nullptr);
   (*bu).writeMPSformatRhs(out, 1, icupp);
   (*bl).writeMPSformatRhs(out, 2, iclow);

   // write all variable bounds
   out <<  "BOUNDS " << endl;
   (*bux).writeMPSformatBounds(out, ixupp, true);
   (*blx).writeMPSformatBounds(out, ixlow, false);

   out <<  "ENDATA " << endl;

   cout<<"Finished writing MPS format."<<endl;
}

void sData::writeMPSColumns(ostream& out)
{
   int world_size;
   MPI_Comm_size(MPI_COMM_WORLD, &world_size);
   assert( world_size == 1 );

   int n;
   string varName;
   string rowNameStub;
   string rowNameStubLT;
   string rowNameStubGT;
   StochVector& gStoch = dynamic_cast<StochVector&>(*g);
   StochVector& icuppStoch = dynamic_cast<StochVector&>(*icupp);
   StochVector& iclowStoch = dynamic_cast<StochVector&>(*iclow);
   StochGenMatrix& AStoch = dynamic_cast<StochGenMatrix&>(*A);
   SparseGenMatrix& ASparseTrans = dynamic_cast<SparseGenMatrix*>(AStoch.Bmat)->getTranspose();
   StochGenMatrix& CStoch = dynamic_cast<StochGenMatrix&>(*C);
   SparseGenMatrix& CSparseTrans = dynamic_cast<SparseGenMatrix*>(CStoch.Bmat)->getTranspose();

   SimpleVector* gSimple = dynamic_cast<SimpleVector*>(gStoch.vec);
   n = gSimple->n;

   std::stringstream sstmCol;
   std::stringstream sstmRow;


   // linking variables:
   for( int col = 0; col<n; col++ )
   {
      sstmCol.clear();
      sstmCol.str("");
      sstmCol << " var_L_" << col;
      varName = sstmCol.str();

      // cost coefficients:
      rowNameStub = "COST";
      if( gSimple->elements()[col] != 0 )
         out<<varName<< " " << rowNameStub << " " << gSimple->elements()[col] <<endl;

      // coefficients in A_0:
      rowNameStub = "row_E_R_";
      for( int k = ASparseTrans.krowM()[col]; k<ASparseTrans.krowM()[col+1]; k++ )
         out<<varName<< " " << rowNameStub << ASparseTrans.jcolM()[k] << " " << ASparseTrans.M()[k] <<endl;

      // coefficients in F_0:
      if( AStoch.Blmat )
      {
         SparseGenMatrix& ABlmatSparseTrans = dynamic_cast<SparseGenMatrix*>(AStoch.Blmat)->getTranspose();
         rowNameStub = "row_E_L_";
         for( int k = ABlmatSparseTrans.krowM()[col]; k<ABlmatSparseTrans.krowM()[col+1]; k++ )
            out<<varName<< " " << rowNameStub << ABlmatSparseTrans.jcolM()[k] << " " << ABlmatSparseTrans.M()[k] <<endl;
         dynamic_cast<SparseGenMatrix*>(AStoch.Blmat)->deleteTransposed();
      }
      // coefficients in A_i:
      for( size_t it = 0; it < children.size(); it++ )
      {
         SparseGenMatrix& AChildSparseTrans = dynamic_cast<SparseGenMatrix*>(AStoch.children[it]->Amat)->getTranspose();
         sstmRow.clear();
         sstmRow.str("");
         sstmRow << "row_E_"<<(int)it<<"_";
         rowNameStub = sstmRow.str();
         for( int k = AChildSparseTrans.krowM()[col]; k<AChildSparseTrans.krowM()[col+1]; k++ )
            out<<varName<< " " << rowNameStub << AChildSparseTrans.jcolM()[k] << " " << AChildSparseTrans.M()[k] <<endl;
         dynamic_cast<SparseGenMatrix*>(AStoch.children[it]->Amat)->deleteTransposed();
      }

      // coefficients in C_0:
      rowNameStubLT = "row_L_R_";
      rowNameStubGT = "row_G_R_";
      for( int k = CSparseTrans.krowM()[col]; k<CSparseTrans.krowM()[col+1]; k++ )
      {
         int rowIdx = CSparseTrans.jcolM()[k];
         if( dynamic_cast<SimpleVector*>(icuppStoch.vec)->elements()[rowIdx] != 0.0)
            out<<varName<< " " << rowNameStubLT << rowIdx << " " << CSparseTrans.M()[k] <<endl;
         if( dynamic_cast<SimpleVector*>(iclowStoch.vec)->elements()[rowIdx] != 0.0)
            out<<varName<< " " << rowNameStubGT << rowIdx << " " << CSparseTrans.M()[k] <<endl;
      }
      // coefficients in G_0:
      if( CStoch.Blmat )
      {
         SparseGenMatrix& CBlmatSparseTrans = dynamic_cast<SparseGenMatrix*>(CStoch.Blmat)->getTranspose();
         rowNameStubLT = "row_L_L_";
         rowNameStubGT = "row_G_L_";
         for( int k = CBlmatSparseTrans.krowM()[col]; k<CBlmatSparseTrans.krowM()[col+1]; k++ )
         {
            int rowIdx = CBlmatSparseTrans.jcolM()[k];
            if( dynamic_cast<SimpleVector*>(icuppStoch.vecl)->elements()[rowIdx] != 0.0)
               out<<varName<< " " << rowNameStubLT << rowIdx << " " << CBlmatSparseTrans.M()[k] <<endl;
            if( dynamic_cast<SimpleVector*>(iclowStoch.vecl)->elements()[rowIdx] != 0.0)
               out<<varName<< " " << rowNameStubGT << rowIdx << " " << CBlmatSparseTrans.M()[k] <<endl;
         }
         dynamic_cast<SparseGenMatrix*>(CStoch.Blmat)->deleteTransposed();
      }
      // coefficients in C_i:
      for( size_t it = 0; it < children.size(); it++ )
      {
         SparseGenMatrix& CChildSparseTrans = dynamic_cast<SparseGenMatrix*>(CStoch.children[it]->Amat)->getTranspose();
         sstmRow.clear();
         sstmRow.str("");
         sstmRow << "row_L_"<<(int)it<<"_";
         rowNameStubLT = sstmRow.str();
         sstmRow.clear();
         sstmRow.str("");
         sstmRow << "row_G_"<<(int)it<<"_";
         rowNameStubGT = sstmRow.str();
         for( int k = CChildSparseTrans.krowM()[col]; k<CChildSparseTrans.krowM()[col+1]; k++ )
         {
            int rowIdx = CChildSparseTrans.jcolM()[k];
            if( dynamic_cast<SimpleVector*>(icuppStoch.children[it]->vec)->elements()[rowIdx] != 0.0)
               out<<varName<< " " << rowNameStubLT << rowIdx << " " << CChildSparseTrans.M()[k] <<endl;
            if( dynamic_cast<SimpleVector*>(iclowStoch.children[it]->vec)->elements()[rowIdx] != 0.0)
               out<<varName<< " " << rowNameStubGT << rowIdx << " " << CChildSparseTrans.M()[k] <<endl;
         }
         dynamic_cast<SparseGenMatrix*>(CStoch.children[it]->Amat)->deleteTransposed();
      }
   }

   // non-linking variables:
   for( size_t it = 0; it < children.size(); it++ )
   {
      SimpleVector* gSimple = dynamic_cast<SimpleVector*>(gStoch.children[it]->vec);
      n = gSimple->n;

      for( int col = 0; col<n; col++ )
      {
         sstmCol.clear();
         sstmCol.str("");
         sstmCol << " var_"<<(int)it <<"_" << col;
         varName = sstmCol.str();

         // coeffs in COST:
         rowNameStub = "COST";
         if( gSimple->elements()[col] != 0 )
            out<<varName<< " " << rowNameStub << " " << gSimple->elements()[col] <<endl;

         // coeffs in A_i:
         SparseGenMatrix& AChildSparseTrans = dynamic_cast<SparseGenMatrix*>(AStoch.children[it]->Bmat)->getTranspose();
         sstmRow.clear();
         sstmRow.str("");
         sstmRow << "row_E_"<<(int)it<<"_";
         rowNameStub = sstmRow.str();
         for( int k = AChildSparseTrans.krowM()[col]; k<AChildSparseTrans.krowM()[col+1]; k++ )
            out<<varName<< " " << rowNameStub << AChildSparseTrans.jcolM()[k] << " " << AChildSparseTrans.M()[k] <<endl;
         dynamic_cast<SparseGenMatrix*>(AStoch.children[it]->Bmat)->deleteTransposed();

         // coefficients in D_i:
         SparseGenMatrix& CChildSparseTrans = dynamic_cast<SparseGenMatrix*>(CStoch.children[it]->Bmat)->getTranspose();
         sstmRow.clear();
         sstmRow.str("");
         sstmRow << "row_L_"<<(int)it<<"_";
         rowNameStubLT = sstmRow.str();
         sstmRow.clear();
         sstmRow.str("");
         sstmRow << "row_G_"<<(int)it<<"_";
         rowNameStubGT = sstmRow.str();
         for( int k = CChildSparseTrans.krowM()[col]; k<CChildSparseTrans.krowM()[col+1]; k++ )
         {
            int rowIdx = CChildSparseTrans.jcolM()[k];
            if( dynamic_cast<SimpleVector*>(icuppStoch.children[it]->vec)->elements()[rowIdx] != 0.0)
               out<<varName<< " " << rowNameStubLT << CChildSparseTrans.jcolM()[k] << " " << CChildSparseTrans.M()[k] <<endl;
            if( dynamic_cast<SimpleVector*>(iclowStoch.children[it]->vec)->elements()[rowIdx] != 0.0)
               out<<varName<< " " << rowNameStubGT << CChildSparseTrans.jcolM()[k] << " " << CChildSparseTrans.M()[k] <<endl;
         }
         dynamic_cast<SparseGenMatrix*>(CStoch.children[it]->Bmat)->deleteTransposed();

         // coefficients in F_i:
         if( dynamic_cast<StochGenMatrix*>(AStoch.children[it])->Blmat )
         {
            SparseGenMatrix& ABlmatSparseTrans = dynamic_cast<SparseGenMatrix*>(AStoch.children[it]->Blmat)->getTranspose();
            rowNameStub = "row_E_L_";
            for( int k = ABlmatSparseTrans.krowM()[col]; k<ABlmatSparseTrans.krowM()[col+1]; k++ )
               out<<varName<< " " << rowNameStub << ABlmatSparseTrans.jcolM()[k] << " " << ABlmatSparseTrans.M()[k] <<endl;
            dynamic_cast<SparseGenMatrix*>(AStoch.children[it]->Blmat)->deleteTransposed();
         }

         // coefficients in G_i:
         if( dynamic_cast<StochGenMatrix*>(CStoch.children[it])->Blmat )
         {
            SparseGenMatrix& CBlmatSparseTrans = dynamic_cast<SparseGenMatrix*>(CStoch.children[it]->Blmat)->getTranspose();
            rowNameStubLT = "row_L_L_";
            rowNameStubGT = "row_G_L_";
            for( int k = CBlmatSparseTrans.krowM()[col]; k<CBlmatSparseTrans.krowM()[col+1]; k++ )
            {
               int rowIdx = CBlmatSparseTrans.jcolM()[k];
               if( dynamic_cast<SimpleVector*>(icuppStoch.vecl)->elements()[rowIdx] != 0.0)
                  out<<varName<< " " << rowNameStubLT << rowIdx << " " << CBlmatSparseTrans.M()[k] <<endl;
               if( dynamic_cast<SimpleVector*>(iclowStoch.vecl)->elements()[rowIdx] != 0.0)
                  out<<varName<< " " << rowNameStubGT << rowIdx << " " << CBlmatSparseTrans.M()[k] <<endl;
            }
            dynamic_cast<SparseGenMatrix*>(CStoch.children[it]->Blmat)->deleteTransposed();
         }
      }
   }

   // delete transposed matrices:
   dynamic_cast<SparseGenMatrix*>(AStoch.Bmat)->deleteTransposed();
   dynamic_cast<SparseGenMatrix*>(CStoch.Bmat)->deleteTransposed();

}

sData* sData::cloneFull(bool switchToDynamicStorage) const
{
   // todo Q is empty!
   StochSymMatrixHandle Q_clone(dynamic_cast<const StochSymMatrix&>(*Q).clone());
   StochGenMatrixHandle A_clone(dynamic_cast<const StochGenMatrix&>(*A).cloneFull(switchToDynamicStorage));
   StochGenMatrixHandle C_clone(dynamic_cast<const StochGenMatrix&>(*C).cloneFull(switchToDynamicStorage));

   StochVectorHandle c_clone (dynamic_cast<StochVector*>(g->cloneFull()));
   StochVectorHandle bA_clone ( dynamic_cast<StochVector*>(bA->cloneFull()));
   StochVectorHandle xupp_clone (dynamic_cast<StochVector*>(bux->cloneFull()));
   StochVectorHandle ixupp_clone (dynamic_cast<StochVector*>(ixupp->cloneFull()));
   StochVectorHandle xlow_clone (dynamic_cast<StochVector*>(blx->cloneFull()));
   StochVectorHandle ixlow_clone (dynamic_cast<StochVector*>(ixlow->cloneFull()));
   StochVectorHandle cupp_clone (dynamic_cast<StochVector*>(bu->cloneFull()));
   StochVectorHandle icupp_clone (dynamic_cast<StochVector*>(icupp->cloneFull()));
   StochVectorHandle clow_clone (dynamic_cast<StochVector*>(bl->cloneFull()));
   StochVectorHandle iclow_clone (dynamic_cast<StochVector*>(iclow->cloneFull()));

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

void sData::permuteLinkingCons()
{
   assert(linkConsPermutationA.size() == 0);
   assert(linkConsPermutationC.size() == 0);

   const size_t nBlocks = dynamic_cast<StochVector&>(*g).children.size();

   // compute permutation vectors
   linkConsPermutationA = getAscending2LinkPermutation(linkStartBlockIdA, nBlocks);
   linkConsPermutationC = getAscending2LinkPermutation(linkStartBlockIdC, nBlocks);

   assert(permutationIsValid(linkConsPermutationA));
   assert(permutationIsValid(linkConsPermutationC));

   dynamic_cast<StochGenMatrix&>(*A).permuteLinkingCons(linkConsPermutationA);
   dynamic_cast<StochGenMatrix&>(*C).permuteLinkingCons(linkConsPermutationC);
   dynamic_cast<StochVector&>(*bA).permuteLinkingEntries(linkConsPermutationA);
   dynamic_cast<StochVector&>(*bl).permuteLinkingEntries(linkConsPermutationC);
   dynamic_cast<StochVector&>(*bu).permuteLinkingEntries(linkConsPermutationC);
   dynamic_cast<StochVector&>(*iclow).permuteLinkingEntries(linkConsPermutationC);
   dynamic_cast<StochVector&>(*icupp).permuteLinkingEntries(linkConsPermutationC);
}

void sData::permuteLinkingVars()
{
   assert(linkVarsPermutation.size() == 0);

   linkVarsPermutation = get0VarsRightPermutation(linkVarsNnz);

   assert(permutationIsValid(linkVarsPermutation));

   dynamic_cast<StochGenMatrix&>(*A).permuteLinkingVars(linkVarsPermutation);
   dynamic_cast<StochGenMatrix&>(*C).permuteLinkingVars(linkVarsPermutation);
   dynamic_cast<StochVector&>(*g).permuteVec0Entries(linkVarsPermutation);
   dynamic_cast<StochVector&>(*bux).permuteVec0Entries(linkVarsPermutation);
   dynamic_cast<StochVector&>(*blx).permuteVec0Entries(linkVarsPermutation);
   dynamic_cast<StochVector&>(*ixupp).permuteVec0Entries(linkVarsPermutation);
   dynamic_cast<StochVector&>(*ixlow).permuteVec0Entries(linkVarsPermutation);
}


sVars* sData::getVarsUnperm(const sVars& vars, const sData& unpermData) const
{
   const std::vector<unsigned int> perm_inv_link_vars = this->getLinkVarsPermInv();   
   const std::vector<unsigned int> perm_inv_link_cons_eq = this->getLinkConsEqPermInv();   
   const std::vector<unsigned int> perm_inv_link_cons_ineq = this->getLinkConsIneqPermInv();   

   sVars* unperm_vars = new sVars(vars, unpermData.ixlow, unpermData.ixupp, unpermData.iclow, unpermData.icupp);

   if( perm_inv_link_vars.size() != 0 )
   {
      dynamic_cast<StochVector&>(*unperm_vars->x).permuteVec0Entries(perm_inv_link_vars);   
      dynamic_cast<StochVector&>(*unperm_vars->v).permuteVec0Entries(perm_inv_link_vars);   
      dynamic_cast<StochVector&>(*unperm_vars->w).permuteVec0Entries(perm_inv_link_vars);   
      dynamic_cast<StochVector&>(*unperm_vars->phi).permuteVec0Entries(perm_inv_link_vars);   
      dynamic_cast<StochVector&>(*unperm_vars->gamma).permuteVec0Entries(perm_inv_link_vars);   
   }

   if( perm_inv_link_cons_eq.size() != 0 )
   {
      dynamic_cast<StochVector&>(*unperm_vars->y).permuteLinkingEntries(perm_inv_link_cons_eq);
   }

   if( perm_inv_link_cons_ineq.size() != 0 )
   {
      dynamic_cast<StochVector&>(*unperm_vars->s).permuteLinkingEntries(perm_inv_link_cons_ineq); 
      dynamic_cast<StochVector&>(*unperm_vars->z).permuteLinkingEntries(perm_inv_link_cons_ineq); 
      dynamic_cast<StochVector&>(*unperm_vars->t).permuteLinkingEntries(perm_inv_link_cons_ineq); 
      dynamic_cast<StochVector&>(*unperm_vars->u).permuteLinkingEntries(perm_inv_link_cons_ineq); 
      dynamic_cast<StochVector&>(*unperm_vars->pi).permuteLinkingEntries(perm_inv_link_cons_ineq); 
      dynamic_cast<StochVector&>(*unperm_vars->lambda).permuteLinkingEntries(perm_inv_link_cons_ineq); 
   }

   return unperm_vars;
}

sResiduals* sData::getResidsUnperm(const sResiduals& resids, const sData& unpermData) const
{
   const std::vector<unsigned int> perm_inv_link_vars = this->getLinkVarsPermInv();   
   const std::vector<unsigned int> perm_inv_link_cons_eq = this->getLinkConsEqPermInv();   
   const std::vector<unsigned int> perm_inv_link_cons_ineq = this->getLinkConsIneqPermInv();   

   sResiduals* unperm_resids = new sResiduals(resids, unpermData.ixlow, unpermData.ixupp, unpermData.iclow, unpermData.icupp );

   if( perm_inv_link_vars.size() != 0 )
   {
      dynamic_cast<StochVector&>(*unperm_resids->rQ).permuteVec0Entries(perm_inv_link_vars);   
      dynamic_cast<StochVector&>(*unperm_resids->rv).permuteVec0Entries(perm_inv_link_vars);   
      dynamic_cast<StochVector&>(*unperm_resids->rw).permuteVec0Entries(perm_inv_link_vars);   
      dynamic_cast<StochVector&>(*unperm_resids->rgamma).permuteVec0Entries(perm_inv_link_vars);   
      dynamic_cast<StochVector&>(*unperm_resids->rphi).permuteVec0Entries(perm_inv_link_vars);   
   }

   if( perm_inv_link_cons_eq.size() != 0 )
   {
      dynamic_cast<StochVector&>(*unperm_resids->rA).permuteLinkingEntries(perm_inv_link_cons_eq);
   }

   if( perm_inv_link_cons_ineq.size() != 0 )
   {
      dynamic_cast<StochVector&>(*unperm_resids->rC).permuteLinkingEntries(perm_inv_link_cons_ineq); 
      dynamic_cast<StochVector&>(*unperm_resids->rt).permuteLinkingEntries(perm_inv_link_cons_ineq);   
      dynamic_cast<StochVector&>(*unperm_resids->ru).permuteLinkingEntries(perm_inv_link_cons_ineq);   
      dynamic_cast<StochVector&>(*unperm_resids->rz).permuteLinkingEntries(perm_inv_link_cons_ineq);
      dynamic_cast<StochVector&>(*unperm_resids->rlambda).permuteLinkingEntries(perm_inv_link_cons_ineq);   
      dynamic_cast<StochVector&>(*unperm_resids->rpi).permuteLinkingEntries(perm_inv_link_cons_ineq);   
   }

   return unperm_resids;
}

void sData::activateLinkStructureExploitation()
{
   if( useLinkStructure )
      return;
   useLinkStructure = true;

   const int myrank = PIPS_MPIgetRank(MPI_COMM_WORLD);

   /* don't attempt to use linking structure when there actually is no linking constraints */
   if( stochNode->myl() == 0 && stochNode->mzl() == 0 )
   {
      useLinkStructure = false;
      if( myrank == 0 )
         std::cout << "no linking constraints so no linking structure found" << std::endl;
      return;
   }

   const int nx0 = getLocalnx();

   const StochGenMatrix& Astoch = dynamic_cast<const StochGenMatrix&>(*A);
   const StochGenMatrix& Cstoch = dynamic_cast<const StochGenMatrix&>(*C);

   linkVarsNnz = std::vector<int>(nx0, 0);

   int n2LinksEq = 0;
   int n2LinksIneq = 0;

   Astoch.getLinkVarsNnz(linkVarsNnz);
   Cstoch.getLinkVarsNnz(linkVarsNnz);

   linkStartBlockIdA = Astoch.get2LinkStartBlocks();
   linkStartBlockIdC = Cstoch.get2LinkStartBlocks();

   linkStartBlockLengthsA = get2LinkLengthsVec(linkStartBlockIdA, stochNode->children.size());
   linkStartBlockLengthsC = get2LinkLengthsVec(linkStartBlockIdC, stochNode->children.size());

   printLinkConsStats();
   printLinkVarsStats();

   for( size_t i = 0; i < linkVarsNnz.size(); ++i )
      if( linkVarsNnz[i] == 0 )
         n0LinkVars++;

   for( size_t i = 0; i < linkStartBlockIdA.size(); ++i )
      if( linkStartBlockIdA[i] >= 0 )
         n2LinksEq++;

   for( size_t i = 0; i < linkStartBlockIdC.size(); ++i )
      if( linkStartBlockIdC[i] >= 0 )
         n2LinksIneq++;

   assert(n2LinksEq == n2linkRowsEq());
   assert(n2LinksIneq == n2linkRowsIneq());

   if( myrank == 0 )
   {
      std::cout << "number of 0-link variables: " << n0LinkVars << " (out of "
            << nx0 << " link variables) " << std::endl;
      std::cout << "number of equality 2-links: " << n2LinksEq << " (out of "
            << linkStartBlockIdA.size() << " equalities) " << std::endl;
      std::cout << "number of inequality 2-links: " << n2LinksIneq << " (out of "
            << linkStartBlockIdC.size() << " equalities) " << std::endl;

      std::cout << "ratio: "
            << (n2LinksEq + n2LinksIneq) / ((double) linkStartBlockIdA.size() + linkStartBlockIdC.size()) << std::endl;
   }

   if( (n2LinksEq + n2LinksIneq + n0LinkVars) / double(linkStartBlockIdA.size() + linkStartBlockIdC.size() + linkVarsNnz.size()) < minStructuredLinksRatio )
   {
      if( myrank == 0 )
         std::cout << "not enough linking structure found" << std::endl;
      useLinkStructure = false;
   }

   if( useLinkStructure )
   {
      assert(linkStartBlockIdA.size() == unsigned(stochNode->myl()));
      assert(linkStartBlockIdC.size() == unsigned(stochNode->mzl()));

   #ifndef NDEBUG
      const int myl = stochNode->myl();
      const int mzl = stochNode->mzl();
      assert(myl >= 0 && mzl >= 0 && (mzl + myl > 0));
   #endif

      permuteLinkingCons();
      permuteLinkingVars();
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

   std::vector<int> linkCountA(n, 0);
   std::vector<int> linkCountC(n, 0);
   std::vector<int> linkCount0(n, 0);
   std::vector<int> linkCountLC(n, 0);

   StochGenMatrix& Astoch = dynamic_cast<StochGenMatrix&>(*A);
   StochGenMatrix& Cstoch = dynamic_cast<StochGenMatrix&>(*C);

   Astoch.updateKLinkVarsCount(linkCountA);
   Cstoch.updateKLinkVarsCount(linkCountC);

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
         const int linkCountAB = linkCountA[i] + linkCountC[i];
         assert(linkCountAB >= 0 && linkCount0[i] >= 0 && linkCountLC[i] >= 0);
         assert(linkCount0[i] <= 2 && linkCountLC[i] <= 2);

         if( linkCountAB < nLinkStats )
            linkSizes[size_t(linkCountAB)]++;

         if( linkCountAB == 0 && linkCountLC[i] == 0 && linkCount0[i] != 0 )
            count0++;

         if( linkCountAB == 0 && linkCount0[i] == 0 && linkCountLC[i] != 0 )
            countLC++;

         if( linkCountAB == 0 && (linkCount0[i] != 0 || linkCountLC[i] != 0) )
            count0LC++;
      }


      int nlocal = 0;
      for( int i = 0; i < nLinkStats; i++ )
         if( linkSizes[i] != 0 )
         {
            nlocal += linkSizes[i];
            std::cout << i << "-link vars: " << linkSizes[i] << std::endl;
         }

      assert(n - nlocal >= 0);
      std::cout << "---total linking variables: " << n << " (global: " << n - nlocal << ")" <<   std::endl;

      std::cout << "   Block0 exclusive vars " << count0 << std::endl;
      std::cout << "   LC exclusive vars " << countLC << std::endl;
      std::cout << "   Block0 or LC vars " << count0LC  << std::endl;
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

         int nlocal = 0;
         for( int i = 0; i < nLinkStats; i++ )
            if( linkSizes[i] != 0 )
            {
               nlocal += linkSizes[i];
               std::cout << "equality " <<  i << "-link cons: " << linkSizes[i] << std::endl;
            }
         std::cout << "---total equality linking constraints: " << myl << " (global: " << myl - nlocal << ")" <<   std::endl;

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

         int nlocal = 0;
         for( int i = 0; i < nLinkStats; i++ )
            if( linkSizes[i] != 0 )
            {
               nlocal += linkSizes[i];
               std::cout << "inequality " <<  i << "-link cons: " << linkSizes[i] << std::endl;
            }
         std::cout << "---total inequality linking constraints: " << mzl << " (global: " << mzl - nlocal << ")" <<   std::endl;
      }
   }
}

sData::~sData()
{
   for( size_t it = 0; it < children.size(); it++ )
      delete children[it];
}

std::vector<unsigned int> sData::getLinkVarsPerm() const
{
   std::vector<unsigned int> copy = linkVarsPermutation;
   return copy;
}
std::vector<unsigned int> sData::getLinkVarsPermInv() const
{
   return getInversePermutation(linkVarsPermutation);
}
std::vector<unsigned int> sData::getLinkConsEqPermInv() const
{
   return getInversePermutation(linkConsPermutationA);
}
std::vector<unsigned int> sData::getLinkConsIneqPermInv() const
{
   return getInversePermutation(linkConsPermutationC);
}

int sData::getLocalnx()
{
   long long my, nx;
   StochGenMatrix& Ast = dynamic_cast<StochGenMatrix&>(*A);
   Ast.Bmat->getSize(my, nx);
   return nx;
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

   assert(n0 >= n0LinkVars);

   // sum up half of dense square
   nnz += nnzTriangular(n0);

   // add B_0 (or A_0, depending on notation)
   nnz += getLocalB().numberOfNonZeros();

   // add borders
   nnz += myl * (n0 - n0LinkVars);
   nnz += mzl * (n0 - n0LinkVars);

   // (empty) diagonal
   nnz += my;

   // add linking equality parts
   nnz += getSCdiagBlocksMaxNnz(linkStartBlockIdA.size(), linkStartBlockLengthsA);

   // add linking inequality parts
   nnz += getSCdiagBlocksMaxNnz(linkStartBlockIdC.size(), linkStartBlockLengthsC);

   // add linking mixed parts
   nnz += getSCmixedBlocksMaxNnz(linkStartBlockIdA.size(), linkStartBlockIdC.size(),
                                 linkStartBlockLengthsA, linkStartBlockLengthsC);

   if( myl > 0 )
   {
      SparseGenMatrix& Ft = getLocalF().getTranspose();
      const int* startRowFtrans = Ft.krowM();
      nnz += startRowFtrans[n0] - startRowFtrans[n0 - n0LinkVars];
   }

   if( mzl > 0 )
   {
      SparseGenMatrix& Gt = getLocalG().getTranspose();
      const int* startRowGtrans = Gt.krowM();
      nnz += startRowGtrans[n0] - startRowGtrans[n0 - n0LinkVars];
   }

   return nnz;
}


int sData::getSchurCompMaxNnzDist(int blocksStart, int blocksEnd)
{
   assert(children.size() > 0);

   const int n0 = getLocalnx();
   const int my = getLocalmy();
   const int myl = getLocalmyl();
   const int mzl = getLocalmzl();
   const int mylLocal = myl - getSCdiagBlocksNRows(linkStartBlockLengthsA)
      + getSCdiagBlocksNRows(linkStartBlockLengthsA, blocksStart, blocksEnd);
   const int mzlLocal = mzl - getSCdiagBlocksNRows(linkStartBlockLengthsC)
      + getSCdiagBlocksNRows(linkStartBlockLengthsC, blocksStart, blocksEnd);

#ifndef NDEBUG
   {
      int mB, nB;
      getLocalB().getSize(mB, nB);
      assert(mB == my  && nB == n0);
   }
#endif

   int nnz = 0;

   assert(n0 >= n0LinkVars);

   // sum up half of dense square
   nnz += nnzTriangular(n0);

   // add B_0 (or A_0, depending on notation)
   nnz += getLocalB().numberOfNonZeros();

   // add borders
   nnz += mylLocal * (n0 - n0LinkVars);
   nnz += mzlLocal * (n0 - n0LinkVars);

   // (empty) diagonal
   nnz += my;

   // add linking equality parts
   nnz += getSCdiagBlocksMaxNnzDist(linkStartBlockIdA.size(), linkStartBlockLengthsA, blocksStart, blocksEnd);

   // add linking inequality parts
   nnz += getSCdiagBlocksMaxNnzDist(linkStartBlockIdC.size(), linkStartBlockLengthsC, blocksStart, blocksEnd);

   // add linking mixed parts
   nnz += getSCmixedBlocksMaxNnzDist(linkStartBlockIdA.size(), linkStartBlockIdC.size(),
                                 linkStartBlockLengthsA, linkStartBlockLengthsC, blocksStart, blocksEnd);

   if( myl > 0 )
   {
      SparseGenMatrix& Ft = getLocalF().getTranspose();
      const int* startRowFtrans = Ft.krowM();
      nnz += startRowFtrans[n0] - startRowFtrans[n0 - n0LinkVars];
   }

   if( mzl > 0 )
   {
      SparseGenMatrix& Gt = getLocalG().getTranspose();
      const int* startRowGtrans = Gt.krowM();
      nnz += startRowGtrans[n0] - startRowGtrans[n0 - n0LinkVars];
   }

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

void sData::cleanUpPresolvedData(const StochVectorBase<int>& rowNnzVecA, const StochVectorBase<int>& rowNnzVecC,
      const StochVectorBase<int>& colNnzVec)
{
   StochSymMatrix& Q_stoch = dynamic_cast<StochSymMatrix&>(*Q);

   // todo only works if Q is empty - not existent?
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

   assert(stochNode != nullptr);

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


void sData::initDistMarker(int blocksStart, int blocksEnd)
{
   assert(isSCrowLocal.size() == 0);
   assert(isSCrowMyLocal.size() == 0);

   assert(linkStartBlockIdA.size() > 0 || linkStartBlockIdC.size() > 0);
   assert(blocksStart >= 0 && blocksStart < blocksEnd && blocksEnd <= int(linkStartBlockLengthsA.size()));

   const int nx0 = getLocalnx();
   const int my0 = getLocalmy();
   const int myl = getLocalmyl();
   const int mzl = getLocalmzl();
   const int sizeSC = nx0 + my0 + myl + mzl;

   assert(sizeSC > 0);

   isSCrowLocal.resize(sizeSC);
   isSCrowMyLocal.resize(sizeSC);

   for( int i = 0; i < nx0; i++ )
   {
      isSCrowLocal[i] = false;
      isSCrowMyLocal[i] = false;
   }

   for( int i = nx0; i < nx0 + my0; i++ )
   {
      isSCrowLocal[i] = false;
      isSCrowMyLocal[i] = false;
   }

   // equality linking
   for( int i = nx0 + my0, j = 0; i < nx0 + my0 + myl; i++, j++ )
   {
      assert( unsigned(j) < linkStartBlockIdA.size() );
      const int block = linkStartBlockIdA[j];
      isSCrowLocal[i] = (block != -1);
      isSCrowMyLocal[i] = (block >= blocksStart && block < blocksEnd);
   }

   // inequality linking
   for( int i = nx0 + my0 + myl, j = 0; i < nx0 + my0 + myl + mzl; i++, j++ )
   {
      assert( unsigned(j) < linkStartBlockIdC.size() );
      const int block = linkStartBlockIdC[j];
      isSCrowLocal[i] = (block != -1);
      isSCrowMyLocal[i] = (block >= blocksStart && block < blocksEnd);
   }
}

const std::vector<bool>& sData::getSCrowMarkerLocal() const
{
   assert(isSCrowLocal.size() != 0);

   return isSCrowLocal;
}

const std::vector<bool>& sData::getSCrowMarkerMyLocal() const
{
   assert(isSCrowMyLocal.size() != 0);

   return isSCrowMyLocal;
}

int sData::n2linkRowsEq() const
{
   return n2linksRows(linkStartBlockLengthsA);
}

int sData::n2linkRowsIneq() const
{
   return n2linksRows(linkStartBlockLengthsC);
}

// is root node data of sData object same on all procs?
bool sData::isRootNodeInSync() const
{
   bool in_sync = true;

   /* matrix Q */
   // todo

   /* matrix A */
   if(!dynamic_cast<const StochGenMatrix&>(*A).isRootNodeInSync())
   {
      std::cout << "ERROR: matrix A corrupted!" << std::endl;
      in_sync = false;
   }

   /* matrix C */
   if( !dynamic_cast<const StochGenMatrix&>(*C).isRootNodeInSync() )
   {
      std::cout << "ERROR: matrix C corrupted!" << std::endl;
      in_sync = false;
   }

   /* objective g */
   if( !dynamic_cast<const StochVector&>(*g).isRootNodeInSync() )
   {
      std::cout << "ERROR: objective vector corrupted!" << std::endl;
      in_sync = false;
   }

   /* rhs equality bA */
   if( !dynamic_cast<const StochVector&>(*bA).isRootNodeInSync() )
   {
      std::cout << "ERROR: rhs of A corrupted!" << std::endl;
      in_sync = false;
   }

   /* upper bounds x bux */
   if( !dynamic_cast<const StochVector&>(*bux).isRootNodeInSync() )
   {
      std::cout << "ERROR: upper bounds x corrupted!" << std::endl;
      in_sync = false;
   }

   /* index for upper bounds x ixupp */
   if( !dynamic_cast<const StochVector&>(*ixupp).isRootNodeInSync() )
   {
      std::cout << "ERROR: index upper bounds x corrupted!" << std::endl;
      in_sync = false;
   }

   /* lower bounds x blx */
   if( !dynamic_cast<const StochVector&>(*blx).isRootNodeInSync() )
   {
      std::cout << "ERROR: lower bounds x corrupted!" << std::endl;
      in_sync = false;
   }

   /* index for lower bounds x ixlow */
   if( !dynamic_cast<const StochVector&>(*ixlow).isRootNodeInSync() )
   {
      std::cout << "ERROR: index lower bounds x corrupted!" << std::endl;
      in_sync = false;
   }

   /* upper bounds C bu */
   if( !dynamic_cast<const StochVector&>(*bu).isRootNodeInSync() )
   {
      std::cout << "ERROR: rhs C corrupted!" << std::endl;
      in_sync = false;
   }

   /* index upper bounds C icupp */
   if( !dynamic_cast<const StochVector&>(*icupp).isRootNodeInSync() )
   {
      std::cout << "ERROR: index rhs C corrupted!" << std::endl;
      in_sync = false;
   }

   /* lower bounds C bl */
   if( !dynamic_cast<const StochVector&>(*bl).isRootNodeInSync() )
   {
      std::cout << "ERROR: lower bounds C corrupted!" << std::endl;
      in_sync = false;
   }

   /* index for lower bounds C iclow */
   if( !dynamic_cast<const StochVector&>(*iclow).isRootNodeInSync() )
   {
      std::cout << "ERROR: index lower bounds C corrupted!" << std::endl;
      in_sync = false;
   }

   /* sacle sc */
   // todo

   return in_sync;
}
