/*
 * SparseStorageDynamic.C
 *
 *  Created on: 07.02.2018
 *      Author: bzfrehfe
 */

#include "SparseStorageDynamic.h"
#include <cassert>
#include <algorithm>
#include <vector>

int SparseStorageDynamic::instances = 0;


SparseStorageDynamic::SparseStorageDynamic(int m, int n, int len, double spareRatio)
 : spareRatio(spareRatio), m(m), n(n), len(len)
{
   assert(m >= 0 && len >= 0);

   rowptr = new ROWPTRS[m + 1];

   if( len != 0)
   {
      M = new double[len];
      jcolM = new int[len];
   }
   else
   {
      M = NULL;
      jcolM = NULL;
   }
}

SparseStorageDynamic::SparseStorageDynamic(const SparseStorage& storage, double spareRatio)
 : spareRatio(spareRatio), m(storage.m), n(storage.n)
{
   assert(spareRatio >= 0.0);
   assert(m >= 0 && storage.len >= 0);

   const int* const orgkrowM = storage.krowM;
   const int* const orgjcolM = storage.jcolM;
   const double* const orgM = storage.M;


   // compute length of new storage

   len = 0;
   for( int i = 0; i < storage.len; i++ )
      if( orgM[i] != 0.0 )
         len++;

   for( int r = 0; r < m; r++ )
   {
      len += int(spareRatio * (orgkrowM[r + 1] - orgkrowM[r]));
      assert(orgkrowM[r + 1] - orgkrowM[r] >= 0);
   }

   if( len > 0 )
   {
      jcolM = new int[len];
      M = new double[len];
   }
   else
   {
      jcolM = NULL;
      M = NULL;
   }

   rowptr = new ROWPTRS[m + 1];

   // build storage

   int shift = 0;
   for( int r = 0; r < m; r++ )
   {
      const int offset = int(spareRatio * (orgkrowM[r + 1] - orgkrowM[r]));

      rowptr[r].start = orgkrowM[r] + shift;

      for( int j = orgkrowM[r]; j < orgkrowM[r + 1]; j++ )
      {
         if( orgM[j] != 0.0 )
         {
            assert(j + shift >= 0);

            M[j + shift] = orgM[j];
            jcolM[j + shift] = orgjcolM[j];
         }
         else
         {
            shift--;
         }
      }

      rowptr[r].end = orgkrowM[r + 1] + shift;

      shift += offset;
   }

   assert(m == 0 || rowptr[m - 1].end + int(spareRatio * (orgkrowM[m] - orgkrowM[m - 1])) == len );

   rowptr[m].start = orgkrowM[m] + shift;
   rowptr[m].end = rowptr[m].start;

   assert(rowptr[m].start == len);

   SparseStorageDynamic::instances++;
}

SparseStorageDynamic::SparseStorageDynamic( const SparseStorageDynamic &dynamicStorage)
: spareRatio(dynamicStorage.spareRatio), m(dynamicStorage.m), n(dynamicStorage.n), len(dynamicStorage.len)
{
   jcolM = new int[len];
   M = new double[len];
   rowptr = new ROWPTRS[m + 1];

   memcpy(jcolM, dynamicStorage.jcolM, len * sizeof(dynamicStorage.jcolM[0]));
   memcpy(M, dynamicStorage.M, len * sizeof(dynamicStorage.M[0]));
   memcpy(rowptr, dynamicStorage.rowptr, (m + 1) * sizeof(dynamicStorage.rowptr[0]));

   // is this necessary?    SparseStorageDynamic::instances++;
}

void SparseStorageDynamic::getSize(int& m, int& n)
{
   m = this->m;
   n = this->n;
}

SparseStorage* SparseStorageDynamic::getStaticStorage(double* rowNnz, double* colNnz) const
{
   int m_static = 0;

   // empty?
   if( n <= 0 )
   {
      assert(len == 0);
      assert(colNnz == NULL);
      assert(rowNnz != NULL);

      for( int r = 0; r < m; r++ )
         if( rowNnz[r] != 0.0 )
            m_static++;

      SparseStorage* staticStorage = new SparseStorage(m_static, n, len);

      return staticStorage;
   }

   assert(rowNnz != NULL && colNnz != NULL);

   // get m, n, len for new storage

   bool* cols = new bool[n];

   for( int i = 0; i < n; i++ )
      cols[i] = false;

   int n_static = 0;
   int len_static = 0;

   for( int r = 0; r < m; r++ )
   {
      const int start = rowptr[r].start;
      const int end = rowptr[r].end;

      int rownnz = 0;

      for( int j = start; j < end; j++ )
      {
         assert(jcolM[j] < n);

         if( M[j] == 0.0 )
            continue;

         cols[jcolM[j]] = true;
         rownnz++;
      }

      assert(rownnz == 0 || rowNnz[r] != 0.0);

      if( rownnz > 0 || rowNnz[r] != 0.0 )
      {
         len_static += rownnz;
         m_static++;
      }
   }

   // now create and fill storage

   int* colsmap = new int[n];

   for( int i = 0; i < n; i++ )
   {
      if( cols[i] )
         colsmap[i] = n_static;
#ifndef NDEBUG
      else
         colsmap[i] = -1;
#endif

      if( cols[i] || colNnz[i] != 0.0 )
         n_static++;
   }

   SparseStorage* staticStorage = new SparseStorage(m_static, n_static, len_static);

   int* const krowM_static = staticStorage->krowM;
   int* const jcolM_static = staticStorage->jcolM;
   double* const M_static = staticStorage->M;

   int nnz_static = 0;
   int rowcount = 0;

   for( int r = 0; r < m; r++ )
   {
      const int start = rowptr[r].start;
      const int end = rowptr[r].end;

      const int nnz_static_old = nnz_static;

      for( int j = start; j < end; j++ )
      {
         if( M[j] == 0.0 )
            continue;

         assert(cols[jcolM[j]]);

         const int col_static = colsmap[jcolM[j]];

#ifndef NDEBUG
         assert(col_static >= 0);
         assert(col_static < n_static);
#endif

         jcolM_static[nnz_static] = col_static;
         M_static[nnz_static++] = M[j];
      }

      if( nnz_static_old != nnz_static || rowNnz[r] != 0.0 )
         krowM_static[rowcount++] = nnz_static_old;
   }

   delete[] cols;
   delete[] colsmap;

   assert(nnz_static == len_static);
   assert(rowcount == m_static);
   krowM_static[rowcount] = nnz_static;

   return staticStorage;
}


SparseStorageDynamic* SparseStorageDynamic::getTranspose() const
{
   assert(n > 0);

   // compute nnz of each row of At (column of A)

   int* w = new int[n];

   for( int i = 0; i < n; i++ )
      w[i] = 0;

   for( int r = 0; r < m; r++ )
   {
      const int start = rowptr[r].start;
      const int end = rowptr[r].end;

      for( int j = start; j < end; j++ )
      {
         assert(M[j] != 0.0);
         w[jcolM[j]]++;
      }
   }

   int translen = 0;

   assert(spareRatio >= 0.0);

   for( int i = 0; i < n; i++ )
      translen += w[i] + int(w[i] * spareRatio);

   SparseStorageDynamic* transpose = new SparseStorageDynamic(n, m, translen, spareRatio);

   // set row pointers

   ROWPTRS* const transrowptr = transpose->rowptr;
   transrowptr[0].start = 0;

   for( int i = 1; i <= n; i++ )
   {
      const int oldstart = transrowptr[i - 1].start;
      const int oldend = oldstart + w[i - 1];
      assert(oldend >= oldstart);

      transrowptr[i - 1].end = oldend;
      transrowptr[i].start = oldend + int((oldend - oldstart) * spareRatio);

      w[i - 1] = oldstart;
   }

   transrowptr[n].start = translen;
   transrowptr[n].end = translen;


   // fill jCol and M

   int* const transjcolM = transpose->jcolM;
   double* const transM = transpose->M;

   for( int r = 0; r < m; r++ )
   {
      const int start = rowptr[r].start;
      const int end = rowptr[r].end;

      for( int j = start; j < end; j++ )
      {
         const int idx = w[jcolM[j]];

         assert(idx < translen);

         transM[idx] = M[j];
         transjcolM[idx] = r;
         w[jcolM[j]] = idx + 1;
      }
   }

   delete[] w;

   return transpose;
}


void SparseStorageDynamic::addNnzPerRow(double* vec) const
{
   for( int r = 0; r < m; r++ )
      vec[r] += rowptr[r].end - rowptr[r].start;
}

void SparseStorageDynamic::writeToStreamDense( ostream& out) const
{
   int i, k;
   //todo: instead of \t, use length of longest value in M

   for( i = 0; i < m; i++ )
   { // Row i
      int j = 0; // Column j
      const int start = rowptr[i].start;
      const int end = rowptr[i].end;
      for( k = start; k < end; k++ )
      {
         while( jcolM[k] > j )
         {
            out << 0 << '\t';
            j++;
         }
         out << M[k] << '\t';
         j++;
      }
      while( j < n )
      {
         out << 0 << '\t';
         j++;
      }
      out << endl;
   }
}

void SparseStorageDynamic::writeToStreamDenseRow( stringstream& out, int rowidx) const
{
   int j = 0; // Column j

   const int start = rowptr[rowidx].start;
   const int end = rowptr[rowidx].end;

   for( int k = start; k < end; k++ )
   {
      while( jcolM[k] > j )
      {
         out << 0 << '\t';
         j++;
      }
      out << M[k] << '\t';
      j++;
   }

   while( j < n )
   {
      out << 0 << '\t';
      j++;
   }
}

void SparseStorageDynamic::restoreOrder()
{
   for(int i=0; i<m; i++)  // row i
   {
      int start = rowptr[i].start;
      int end = rowptr[i].end;

      // todo: this is too much overhead, better to implement a random access iterator
      std::vector<std::pair<int, double> > pairVector;

      for(int j=start; j<end; j++)
         pairVector.push_back(make_pair(jcolM[j], M[j]));

      std::sort(pairVector.begin(), pairVector.end(), first_is_smaller());
      for(int j=start; j<end; j++)
      {
         jcolM[j] = pairVector[j-start].first;
         M[j] = pairVector[j-start].second;
      }
   }
}

SparseStorageDynamic::~SparseStorageDynamic()
{
   delete[] jcolM;
   delete[] rowptr;
   delete[] M;

   SparseStorageDynamic::instances--;
}
