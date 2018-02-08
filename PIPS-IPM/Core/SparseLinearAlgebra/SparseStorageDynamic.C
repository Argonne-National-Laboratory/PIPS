/*
 * SparseStorageDynamic.C
 *
 *  Created on: 07.02.2018
 *      Author: bzfrehfe
 */

#include "SparseStorageDynamic.h"
#include <cassert>

int SparseStorageDynamic::instances = 0;

SparseStorageDynamic::SparseStorageDynamic(const SparseStorage& storage, double spareRatio)
 : spareRatio(spareRatio), m(storage.m), n(storage.n)
{
   assert(spareRatio >= 0.0);

   const int* const orgkrowM = storage.krowM;
   const int* const orgjcolM = storage.jcolM;
   const double* const orgM = storage.M;

   // compute length of new storage

   len = storage.len;

   for( int r = 0; r < m; r++ )
   {
      len += spareRatio * (orgkrowM[r + 1] - orgkrowM[r]) + 1;
      assert(orgkrowM[r + 1] - orgkrowM[r] >= 0);
   }

   jcolM = new int[len];
   M = new double[len];
   rowptr = new ROWPTRS[m + 1];


   // build storage

   int shift = 0;
   for( int r = 0; r < m; r++ )
   {
      const int offset = spareRatio * (orgkrowM[r + 1] - orgkrowM[r]) + 1;

      rowptr[r].start = orgkrowM[r] + shift;

      for( int j = orgkrowM[r]; j < orgkrowM[r + 1]; j++ )
      {
         M[j + shift] = orgM[j];
         jcolM[j + shift] = orgjcolM[j];
      }

      shift += offset;

      rowptr[r].end = orgkrowM[r + 1] + shift;
   }

   rowptr[m].start = orgkrowM[m] + shift;
   rowptr[m].end = rowptr[m].start;

   SparseStorageDynamic::instances++;
}

void SparseStorageDynamic::getSize(int& m, int& n)
{
   m = this->m;
   n = this->n;
}

void SparseStorageDynamic::copyFrom(SparseStorage& targetstorage)
{
   // get m, n, len for new storage

   bool* cols = new bool[n];

   for( int i = 0; i < n; i++ )
      cols[i] = false;

   int n_static = 0;
   int m_static = 0;
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

      if( rownnz > 0 )
      {
         len_static += rownnz;
         m_static++;
      }
   }

   // now create and fill storage

   int* colsmap = new int[n];

   for( int i = 0; i < n; i++ )
      if( cols[i] )
         colsmap[i] = n_static++;
#ifndef NDEBUG
      else
         colsmap[i] = -1;
#endif

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

      if( nnz_static_old != nnz_static )
         krowM_static[rowcount++] = nnz_static_old;
   }

   delete[] cols;
   delete[] colsmap;

   assert(nnz_static == len_static);
   assert(rowcount == m_static);
   krowM_static[rowcount] = nnz_static;

   return staticStorage;
}

SparseStorageDynamic::~SparseStorageDynamic()
{
   SparseStorageDynamic::instances--;
}
