/*
 * SparseStorageDynamic.C
 *
 *  Created on: 07.02.2018
 *      Author: bzfrehfe
 */ 

#include "SparseStorageDynamic.h"
#include "pipsdef.h"
#include "pipsport.h"

#include <cassert>
#include <algorithm>
#include <vector>
#include <numeric>
#include <limits>
#include <cstring>

int SparseStorageDynamic::instances = 0;


SparseStorageDynamic::SparseStorageDynamic(int m, int n, int len, double spareRatio)
 : spareRatio(spareRatio), m(m), m_len(m), n(n), len(len), len_free(len)
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
      M = nullptr;
      jcolM = nullptr;
   }
   SparseStorageDynamic::instances++;
}

SparseStorageDynamic::SparseStorageDynamic(const SparseStorage& storage, double spareRatio)
 : spareRatio(spareRatio), m(storage.m), m_len(storage.m), n(storage.n)
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

   // store actual number of entries
   len_free = len;
   for( int r = 0; r < m; r++ )
   {
      len += static_cast<int>(spareRatio * (orgkrowM[r + 1] - orgkrowM[r]));
      assert(orgkrowM[r + 1] - orgkrowM[r] >= 0);
   }

   // actual size minus actual entries
   len_free = len - len_free;

   if( len > 0 )
   {
      jcolM = new int[len];
      M = new double[len];
   }
   else
   {
      jcolM = nullptr;
      M = nullptr;
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
            assert( orgjcolM[j] < n );
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

SparseStorageDynamic::SparseStorageDynamic(const SparseStorageDynamic &dynamicStorage) :
   spareRatio(dynamicStorage.spareRatio), m(dynamicStorage.m), m_len(dynamicStorage.m_len), n(dynamicStorage.n), len(
      dynamicStorage.len), len_free(dynamicStorage.len_free)
{
   assert(m >= 0 && len >= 0);

   rowptr = new ROWPTRS[m + 1];
   memcpy(rowptr, dynamicStorage.rowptr, (m + 1) * sizeof(dynamicStorage.rowptr[0]));

   if( len > 0 )
   {
      jcolM = new int[len];
      M = new double[len];
      memcpy(jcolM, dynamicStorage.jcolM, len * sizeof(dynamicStorage.jcolM[0]));
      memcpy(M, dynamicStorage.M, len * sizeof(dynamicStorage.M[0]));
   }
   else
   {
      jcolM = nullptr;
      M = nullptr;
   }

   SparseStorageDynamic::instances++;
}

void SparseStorageDynamic::getSize(int& m, int& n) const
{
   m = this->m;
   n = this->n;
}

const ROWPTRS SparseStorageDynamic::getRowPtr(int i) const
{
   assert( 0 <= i && i < m );
   return rowptr[i];
}

const int SparseStorageDynamic::getJcolM(int i) const
{
   assert( 0 <= i && i < len );
   return jcolM[i];
}

const double& SparseStorageDynamic::getMat(int i) const
{
   assert( 0 <= i && i < len );
   return M[i];
}

void SparseStorageDynamic::setMat(int i, double val)
{
   assert( 0 <= i && i < len );
   M[i] = val;
}

SparseStorage* SparseStorageDynamic::getStaticStorage(const int* rowNnz, const int* colNnz) const
{
   int m_static = 0;

   // empty?
   if( n <= 0 )
   {
      assert(len == 0);
      assert(colNnz == nullptr);
      assert(rowNnz != nullptr);

      for( int r = 0; r < m; r++ )
         if( rowNnz[r] != 0.0 )
            m_static++;

      SparseStorage* staticStorage = new SparseStorage(m_static, n, len);

      return staticStorage;
   }

   assert(rowNnz != nullptr && colNnz != nullptr);

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

         if( PIPSisZero(M[j]) )
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
      transrowptr[i].start = oldend + static_cast<int>((oldend - oldstart) * spareRatio);

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

      transpose->len_free -= (end - start);
   }

   delete[] w;

   if(transpose->getNVals() != getNVals() )
      std::cout << getNVals() << "\t" << transpose->getNVals() << std::endl;
   assert( transpose->getNVals() == getNVals() );

   return transpose;
}


void SparseStorageDynamic::addNnzPerRow(int* vec) const
{
#ifdef PRE_CPP11
   for( int r = 0; r < m; r++ )
      vec[r] += rowptr[r].end - rowptr[r].start;
#else
   std::transform(rowptr, rowptr + m, vec, vec, [](ROWPTRS pt, double v)->double { return v + pt.end - pt.start; });
#endif
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

double SparseStorageDynamic::abmaxnorm() const
{
   double max = 0.0;
   for( int i = 0; i < m; ++i )
   {
      for( int j = rowptr[m].start; j < rowptr[m].end; ++j )
      {
         if( std::fabs(M[j]) > max )
            max = std::fabs(M[j]);
      }
   }
   return max;
}

void SparseStorageDynamic::restoreOrder()
{
   for(int i = 0; i < m; i++)  // row i
   {
      const int start = rowptr[i].start;
      const int end = rowptr[i].end;

      // todo: this is too much overhead, better to implement a random access iterator // how do you mean?
      std::vector<std::pair<int, double> > pairVector;

      for(int j = start; j < end; j++)
         pairVector.push_back(make_pair(jcolM[j], M[j]));

      std::sort(pairVector.begin(), pairVector.end(), first_is_smaller());
      for(int j = start; j < end; j++)
      {
         jcolM[j] = pairVector[j - start].first;
         M[j] = pairVector[j - start].second;
      }
   }
}

void SparseStorageDynamic::removeEntryAtIndex(int row, int col_idx)
{
   assert( 0 <= row && row <= m );
   assert( rowptr[row].start <= col_idx && col_idx < rowptr[row].end );

   int& row_end = rowptr[row].end;

   std::swap(M[col_idx], M[row_end - 1]);
   std::swap(jcolM[col_idx], jcolM[row_end - 1]);
   --row_end;
   ++len_free;

   assert(len >= len_free);
}

void SparseStorageDynamic::removeEntryAtRowCol(int row, int col)
{
   assert(0 <= row && row < m);
   assert(0 <= col && col < n);

   int i = -1;
   const int end = rowptr[row].end;
   const int start = rowptr[row].start;

   for( i = start; i < end; i++)
   {
      if( jcolM[i] == col )
         break;
   }
   assert( jcolM[i] == col);

   removeEntryAtIndex( row, i );
}

/** adds col with coeff to row - returns false if col was in row already, true otherwise */
bool SparseStorageDynamic::addColToRow( double coeff, int col, int row )
{
   assert(0 <= row && row < m);
   assert(0 <= col && col < n);

   const int end = rowptr[row].end;
   const int start = rowptr[row].start;

   assert( start <= end );

   int i = start;
   while(i < end)
   {
      if( jcolM[i] == col )
         break;
      ++i;
   }

   if( i != end )
   {
      assert( jcolM[i] == col );
      M[i] += coeff;

      return false;
   }
   else
   {
      assert( i == end );

      /* extend storage if necessary */
      if( rowptr[row].end == rowptr[row + 1].start )
         rebuildSpareStructure(1);

      /* if our row has been empty before it will also be empty after compressing the storage values */
      assert( rowptr[row].end < rowptr[row + 1].start );
      assert( len_free > 0 );

      /* insert entry */
      const int new_col_idx = rowptr[row].end;
      ++rowptr[row].end;
      jcolM[new_col_idx] = col;
      M[new_col_idx] = coeff;
      --len_free;

      return true;
   }
}

void SparseStorageDynamic::clearRow( int row )
{
   assert( 0 <= row && row < m );
   len_free += rowptr[row].end - rowptr[row].start;
   rowptr[row].end = rowptr[row].start;

   assert(len >= len_free);
}

/* slow ! */
void SparseStorageDynamic::clearCol( int col )
{
   assert( 0 <= col && col < n);
   for(int row = 0; row < m; ++row)
   {
      for(int i = rowptr[row].start; i < rowptr[row].end; ++i)
      {
         if(jcolM[i] == col)
         {
            removeEntryAtIndex(row, i);
            --i;
         }
      }
   }
}

void SparseStorageDynamic::appendRow( const SparseStorageDynamic& storage, int row )
{
   assert( storage.getN() <= n );
   if(m_len == 0)
   {
      rowptr[0].start = rowptr[0].end = 0;
   }
   
   /* extract nonzero entries from row */
   std::vector<double> val;
   std::vector<int> idx;
   const int length_row_in_storage = storage.getRowPtr(row).end - storage.getRowPtr(row).start;

   // todo: theoretically this copying is not necessary..
   val.reserve(length_row_in_storage);
   idx.reserve(length_row_in_storage);

   for(int i = storage.getRowPtr(row).start; i < storage.getRowPtr(row).end; ++i)
   {
      if( !PIPSisZero( storage.getMat(i) ) )
      {
         val.push_back(storage.getMat(i));
         idx.push_back(storage.getJcolM(i));
      }
   }

   /* try to add a new row to the storage and extend the row ptr if necessary */
   if( m == m_len )
      extendStorageRows();

   assert(m_len > m);

   /* length of row that should be appended */
   const int total_length_row = val.size() + int(val.size() * spareRatio);
   
   /* check storage space */
   assert( rowptr[m].start == rowptr[m].end );
   assert( (len - rowptr[m].start) >= 0 );
   
   bool extended = false;
   while( total_length_row > (len - rowptr[m].start) )
   {
      if( ( (len_free - len * spareRatio) > static_cast<int>(val.size()) ) && !extended)
      {
         extended = true;
         rebuildSpareStructure();
      }
      extendStorageValues();
   }

   assert( total_length_row <= (len - rowptr[m].start) );

   // actually insert row
   std::copy(val.begin(), val.end(), M + rowptr[m].start);
   std::copy(idx.begin(), idx.end(), jcolM + rowptr[m].start);

   rowptr[m].end = rowptr[m].start + val.size();
   rowptr[m + 1].start = rowptr[m + 1].end = rowptr[m].start + total_length_row;
   ++m;

   len_free -= val.size();

   assert( len >= len_free );
   assert( len_free >= 0 );
}

void SparseStorageDynamic::scaleRow( int row, double factor )
{
   assert( 0 <= row && row < m );

#ifdef PRE_CPP11
   for(int i = rowptr[row].start; i < rowptr[row].end; ++i)
      M[i] *= factor;
#else
   std::transform( M + rowptr[row].start, M + rowptr[row].end, M + rowptr[row].start, 
      [factor](double e) -> double { return e*factor ; } );
#endif
}

/* doubles the size of rowptr */
void SparseStorageDynamic::extendStorageRows()
{
   assert(m == m_len);

   /* double size of old array */
   int m_len_tmp = 2 * m_len;

   if(m_len == 0)
      m_len_tmp = 2;

   /* extend the storage */
   ROWPTRS* rowptr_tmp = new ROWPTRS[m_len_tmp + 1];

   std::copy(rowptr, rowptr + m_len + 1, rowptr_tmp);

   std::swap(rowptr_tmp, rowptr);
   std::swap(m_len_tmp, m_len);

   assert( rowptr[m].start == rowptr[m].end );

   delete[] rowptr_tmp;
}

/** extend all arrays at the end - structure of storage stays the same */
void SparseStorageDynamic::extendStorageValues()
{
   assert(len_free >= 0);
   assert(len >= len_free);

   long len_tmp_long = len;

   /* if initial size of storage was zero set it to 2 */
   if( len_tmp_long == 0 )
   {
      len_tmp_long = 2;

      assert(len_free == 0);
   }
   else
   {
      /* double the storage size */
      assert( len_tmp_long * 2 < std::numeric_limits<int>::max() );
      len_tmp_long *= 2;
   }

   int len_tmp = static_cast<int>(len_tmp_long);

   assert(len_tmp > len);
   len_free += len_tmp - len;

   /* extend the storage */
   int * jcolM_tmp = new int[ len_tmp ];
   double * M_tmp = new double[ len_tmp ];

   if( len != 0 )
   {
      std::copy(jcolM, jcolM + len, jcolM_tmp);
      std::copy(M, M + len, M_tmp);
   }

   std::swap( jcolM, jcolM_tmp );
   std::swap( M, M_tmp );
   std::swap( len, len_tmp );

   delete[] jcolM_tmp;
   delete[] M_tmp;

#ifndef NDEBUG
#ifndef PRE_CPP11
   assert( len_free == (len - std::accumulate(rowptr, rowptr + m, 0, [](const int& a, const ROWPTRS& rp)->int { return a + (rp.end - rp.start); })) );
#endif
#endif

   assert(len >= len_free);
}


/** will rebuild the spare structure adding guaranteed_spare to every rows spare part and thus guaranteeing guaranteed_spare free entries after every row */
void SparseStorageDynamic::rebuildSpareStructure(int guaranteed_spare)
{  
   /* extend storage until we can actually restore the spareRatio */
   assert(len >= len_free);
   assert(rowptr[m].end == rowptr[m].start);

   assert( (len - len_free) * spareRatio + m * guaranteed_spare < std::numeric_limits<int>::max() );

   const int size_for_spare = static_cast<int>((len - len_free) * (1.0 + spareRatio) ) + m * guaranteed_spare;

   while( size_for_spare > len_free )
      extendStorageValues();

   /* shift entries and restore the sparse storage pattern with spareRatio */

   /* copies of current arrays */
   const std::vector<double> M_copy(M, M + len);
   const std::vector<int> jcolM_copy(jcolM, jcolM + len);

   /* offset marks start of current row */
   int offset = 0;
   for( int i = 0; i <= m; ++i)
   {
      const int start = rowptr[i].start;
      const int end = rowptr[i].end;
      const int range = end - start;
      const int spare_range = static_cast<int>(range * spareRatio) + guaranteed_spare;

      assert(offset + range + spare_range < len);

      std::copy(&M_copy.at(start), &M_copy.at(end), M + offset);
      std::copy(&jcolM_copy.at(start), &jcolM_copy.at(end), jcolM + offset);

      rowptr[i].start = offset;
      rowptr[i].end = offset + range;

      offset += range + spare_range;

      assert(offset < len);
   }

   assert( rowptr[m].start == rowptr[m].end );
#ifndef NDEBUG
#ifndef PRE_CPP11
   assert( len_free == len - std::accumulate(rowptr, rowptr + m, 0, [](int a, ROWPTRS rp)->int { return a + (rp.end - rp.start); }));
#endif
#endif
}

double SparseStorageDynamic::rowTimesVec( const double* vec, int length, int row) const
{
   assert(0 <= row && row < m);
   assert(length == n);

   double res = 0.0;
   for( int i = rowptr[row].start; i < rowptr[row].end; ++i)
   {
      assert(jcolM[i] < length);
      res += vec[jcolM[i]] * M[i];
   }

   return res;
}

void SparseStorageDynamic::axpyWithRowAt( double alpha, double* y, int length, int row) const
{
   assert(y);
   assert(0 <= row && row < m);
   assert(length == n);

   for(int i = rowptr[row].start; i < rowptr[row].end; ++i)
   {
      assert(jcolM[i] < length);
      y[jcolM[i]] += alpha * M[i];
   }
}

void SparseStorageDynamic::axpyWithRowAtPosNeg( double alpha, double * y_pos, double* y_neg, int length, int row) const
{
   assert(y_pos);
   assert(0 <= row && row < m);
   assert(length == n);

   for(int i = rowptr[row].start; i < rowptr[row].end; ++i)
   {
      assert( y_pos[jcolM[i]] >= 0 );
      assert( y_neg[jcolM[i]] >= 0 );
      assert(jcolM[i] < length);
      const double val = alpha * M[i];
      const double fin_val = y_pos[jcolM[i]] - y_neg[jcolM[i]] + val;

      if( fin_val > 0 )
      {
         y_pos[jcolM[i]] = fin_val;
         y_neg[jcolM[i]] = 0.0;
      }
      else
      {
         y_pos[jcolM[i]] = 0.0;
         y_neg[jcolM[i]] = -fin_val;
      }
   }
}

void SparseStorageDynamic::getRowMinVec(const double* colScaleVec, double* vec) const
{
   const bool coscale = (colScaleVec != nullptr);

   if( coscale )
   {
      for( int r = 0; r < m; r++ )
      {
         double minval = vec[r];
         assert(minval >= 0.0);

         for( int i = rowptr[r].start; i < rowptr[r].end; i++ )
         {
            const double absval = std::abs(M[i] * colScaleVec[jcolM[i]]);

            if( absval < minval && absval > pips_eps )
               minval = absval;
         }
         vec[r] = minval;
      }
   }
   else
   {
      for( int r = 0; r < m; r++ )
      {
         double minval = vec[r];
         assert(minval >= 0.0);

         for( int i = rowptr[r].start; i < rowptr[r].end; i++ )
         {
            const double absval = std::abs(M[i]);

            if( absval < minval && absval > pips_eps )
               minval = absval;
         }
         vec[r] = minval;
      }
   }
}

void SparseStorageDynamic::getRowMaxVec(const double* colScaleVec, double* vec) const
{
   const bool coscale = (colScaleVec != nullptr);

   if( coscale )
   {
      for( int r = 0; r < m; r++ )
      {
         double maxval = vec[r];
         assert(maxval >= 0.0);

         for( int i = rowptr[r].start; i < rowptr[r].end; i++ )
         {
            const double absval = std::abs(M[i] * colScaleVec[jcolM[i]]);

            if( absval > maxval )
               maxval = absval;
         }
         vec[r] = maxval;
      }
   }
   else
   {
      for( int r = 0; r < m; r++ )
      {
         double maxval = vec[r];
         assert(maxval >= 0.0);

         for( int i = rowptr[r].start; i < rowptr[r].end; i++ )
         {
            const double absval = std::abs(M[i]);

            if( absval > maxval )
               maxval = absval;
         }
         vec[r] = maxval;
      }
   }
}

void SparseStorageDynamic::getRowMinMaxVec(bool getMin, const double* colScaleVec, double* vec) const
{
   if( getMin )
      getRowMinVec(colScaleVec, vec);
   else
      getRowMaxVec(colScaleVec, vec);
}


SparseStorageDynamic::~SparseStorageDynamic()
{
   delete[] jcolM;
   delete[] rowptr;
   delete[] M;

   SparseStorageDynamic::instances--;
}

