/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#include <cmath>
#include <cstring>
#include <cassert>
#include "SparseStorage.h"
#include "OoqpVector.h"
#include "SimpleVector.h"
#include "pipsdef.h"
#include <limits>
#include <fstream>
#include <string>
#include <algorithm>
#include "pipsport.h"

int SparseStorage::instances = 0;

SparseStorage::SparseStorage( int m_, int n_, int len_ )
{
  int i;
  
  neverDeleteElts = 0;
  m      = m_;
  n      = n_;
  len    = len_;
  jcolM  = new int[len];
  krowM  = new int[m+1];
  for( i = 0; i <= m; i++ ) {
    krowM[i] = 0;
  }
  M      = new double[len];

  isFortranIndexed = false;

  SparseStorage::instances++;
}

SparseStorage::SparseStorage( int m_, int n_, int len_,
			      int * krowM_, int * jcolM_, double * M_,
			      int deleteElts)
{
  neverDeleteElts = (!deleteElts);

  //neverDeleteElts = 1;
  m               = m_;
  n               = n_;
  len             = len_;
  jcolM           = jcolM_;
  krowM           = krowM_;
  M               = M_;

  isFortranIndexed = false;

  SparseStorage::instances++;
}

SparseStorage::~SparseStorage()
{
  if ( !neverDeleteElts ) {
    delete [] jcolM;
    delete [] krowM;
    delete [] M;
  }

  SparseStorage::instances--;
}


void SparseStorage::copyFrom(int * krowM_, int * jcolM_, double * M_) const
{
   memcpy(jcolM_, jcolM, len * sizeof(jcolM[0]));
   memcpy(M_, M, len * sizeof(M[0]));
   memcpy(krowM_, krowM, (m + 1) * sizeof(krowM[0]));
}

void SparseStorage::getSize( int& m_, int& n_ ) const
{
  m_ = m;
  n_ = n;
}


void SparseStorage::fromGetDiagonal( int idiag, OoqpVector& vec_in )
{
  SimpleVector & vec = dynamic_cast<SimpleVector &>(vec_in);
  int extent = vec.length();

  assert( idiag + extent <= m );
  assert( idiag + extent <= n );
  
  int i, j, k;

  for ( i = idiag; i < idiag + extent; i++ ) {
    // Loop over all rows in range

    // Set the diagonal element to zero in case we don't find 
    // one in the matrix
    vec[i-idiag] = 0.0;
    for( k = krowM[i]; k < krowM[i+1]; k++ ) {
      // Loop over the elements of the sparse row
      j = jcolM[k];
      if ( i == j ) {
	vec[i-idiag] = M[k];
      }
    } // End loop over the elements of the sparse row
  } // End loop over all rows in range

}

void SparseStorage::ColumnScale( OoqpVector& scale_in )
{
  SimpleVector & scale = dynamic_cast<SimpleVector &>(scale_in);

  assert( scale.length() == n );
 
  int i, j, k;

  for ( i = 0; i < m; i++ ) {
    // Loop over all rows in the sparse matrix
    for( k = krowM[i]; k < krowM[i+1]; k++ ) {
      // Loop over the elements of the sparse row
      j = jcolM[k];
      M[k] = M[k] * scale[j];
    } // End loop over the elements of the sparse row
  } // End loop over all rows in the sparse matrix

}

void SparseStorage::RowScale( OoqpVector& scale_in )
{
  SimpleVector & scale = dynamic_cast<SimpleVector &>(scale_in);

  assert( scale.length() == m );

  int i, k;

  for ( i = 0; i < m; i++ ) {
    // Loop over all rows in the sparse matrix
    for( k = krowM[i]; k < krowM[i+1]; k++ ) {
      // Loop over the elements of the sparse row
      M[k] = M[k] * scale[i];
    } // End loop over the elements of the sparse row
  } // End loop over all rows in the sparse matrix

}

void SparseStorage::SymmetricScale( OoqpVector& scale_in )
{
  SimpleVector & scale = dynamic_cast<SimpleVector &>(scale_in);

  assert( scale.length() == n );
  assert( scale.length() == m );

  int i, j, k;

  for ( i = 0; i < m; i++ ) {
    // Loop over all rows in the sparse matrix
    for( k = krowM[i]; k < krowM[i+1]; k++ ) {
      // Loop over the elements of the sparse row
      j = jcolM[k];
      M[k] = M[k] * scale[j] * scale[i];
    } // End loop over the elements of the sparse row
  } // End loop over all rows in the sparse matrix

}


void SparseStorage::scalarMult( double num )
{
  int i, k;

  for ( i = 0; i < m; i++ ) {
    // Loop over all rows in the sparse matrix
    for( k = krowM[i]; k < krowM[i+1]; k++ ) {
      // Loop over the elements of the sparse row
      M[k] = M[k] * num;
    } // End loop over the elements of the sparse row
  } // End loop over all rows in the sparse matrix

}

void SparseStorage::getDiagonal( OoqpVector& vec_in )
{
  this->fromGetDiagonal( 0, vec_in );
}

void SparseStorage::setToDiagonal( OoqpVector& vec_in )
{
  SimpleVector & vec = dynamic_cast<SimpleVector &>(vec_in);
  int diagExtent = (m <= n) ? m : n; 

  assert( diagExtent == vec.length() );

  int i;
  for( i = 0; i <= diagExtent; i++ ) {
    krowM[i] = i; // Initialize to the diagonal matrix of all zeros
  }
  for( ; i <= m; i++ ) {
    krowM[i] = diagExtent;
  }
  
  double * v = &vec[0];
  for( i = 0; i < diagExtent; i++ ) {
    jcolM[i] = i;
    M[i]     = v[i];
  }
}

bool SparseStorage::isValid(bool verbose) const
{
   assert(krowM && jcolM && M);

   if( m < 0 || n < 0 || len < 0)
   {
      printf("isValid: negative size parameter \n");
      return false;
   }

   if( krowM[0] != 0 || krowM[m] != len )
   {
      printf("isValid: krowM broken \n");
      return false;
   }

   for( int i = 0; i < len; i++ )
      if( jcolM[i] < 0 || jcolM[i] >= n )
      {
         printf("isValid: column index out of bounds \n");
         return false;
      }

   for( int i = 0; i < m; i++ )
      if( krowM[i] > krowM[i + 1] )
      {
         printf("isValid: row indices wrongly ordered \n");
         return false;
      }

   return true;
}

bool SparseStorage::isSorted() const
{
   assert(isValid(false));

   for( int i = 0; i < m; i++ )
   {
      for( int j = krowM[i] + 1; j < krowM[i + 1]; j++ )
      {
         const int col = jcolM[j];
         const int prevcol = jcolM[j - 1];

         if( col <= prevcol )
            return false;
      }
   }

   return true;
}

void SparseStorage::fromGetDense( int row, int col, double * A, int lda,
		int rowExtent, int colExtent )
{
	int i, j, k, jcurrent;

	assert( row >= 0 && row + rowExtent <= m );
	assert( col >= 0 && col + colExtent <= n );

	for ( i = row; i < row + rowExtent; i++ ) {
		// Loop over all rows in range
		jcurrent = col - 1;
		for( k = krowM[i]; k < krowM[i+1]; k++ ) {
			// Loop over the elements of the sparse row
			j = jcolM[k];
			if ( j >= col ) {
				// j is big enough to be within range
				if ( j < col + colExtent ) {
					// j is small enough to be within range
					for ( jcurrent++; jcurrent < j; jcurrent++ ) {
						//printf("A1[%d]=%g\n", (i - row) * lda + jcurrent - col, 0.0);
						A[(i - row) * lda + jcurrent - col] = 0.0;
					}
					jcurrent = j;
					//printf("A2[%d]=%g\n", (i - row) * lda + j - col, M[k]);
					A[(i - row) * lda + j - col] = M[k];
				} else { // j is too big. 
					// There will be no more interesting elements in this row
					break;
				} // End else j is too big
			} // End if j is big enough
		} // End loop over element of the sparse row
		//!for( jcurrent++; jcurrent < n; jcurrent++ ) {
		for( jcurrent++; jcurrent < col+colExtent; jcurrent++ ) {
			A[(i - row) * lda + jcurrent - col] = 0.0;
			//printf("A3[%d]=%g\n", (i - row) * lda + jcurrent - col, 0.0);
		}
	} // End loop over all rows in range
}


// used in backsolves
// get a dense block of columns *in column-major format*
// A must be zero'd on input
// allzero is true if there are actually no nonzero entries in this block
void SparseStorage::fromGetColBlock(int col, double *A, int lda, int colExtent, bool &allzero)
{
  int i,j,k;
  for (i = 0; i < m; i++) {
    for ( k = krowM[i]; k < krowM[i+1]; k++) {
      j = jcolM[k];
      if ( j >= col ) {
        if (j < col + colExtent) {
          A[i+(j-col)*lda] = M[k];
          //if (allzero) allzero = false;
	  allzero=false;
        } else {
          break;
        }
      }
    }
  }
}


// used in backsolves
// get a dense block of columns *in column-major format*
// A must be zero'd on input
// allzero is true if there are actually no nonzero entries in this block
// colSparsity contains on exit the sparsity pattern of union of columns
void SparseStorage::fromGetColBlock(int col, double *A, int lda, int colExtent, int* colSparsity, bool &allzero)
{
  int i,j,k;
  for (i = 0; i < m; i++) {
    for ( k = krowM[i]; k < krowM[i+1]; k++) {
      j = jcolM[k];
      if ( j >= col ) {
        if (j < col + colExtent) {
          A[i+(j-col)*lda] = M[k];
          //if (allzero) allzero = false;
	  allzero=false;
	  colSparsity[i]=1;
        } else {
          break;
        }
      }
    }
  }
}


void SparseStorage::fromGetRowsBlock(const int* rowIndices, int nRows, int arrayLineSize,
      int arrayLineOffset, double* rowsArrayDense, int* rowSparsity)
{
   assert(rowsArrayDense && rowIndices);
   assert(arrayLineSize >= 0 && arrayLineOffset >= 0);

   // todo use OMP?
   for( int i = 0; i < nRows; i++ )
   {
      const int r = rowIndices[i];
      assert(r >= 0 && r < m);

      // empty row?
      if( krowM[r] == krowM[r + 1] )
         continue;

      const int offset = i * arrayLineSize + arrayLineOffset;

      for( int c = krowM[r]; c < krowM[r + 1]; c++ )
      {
         const int col = jcolM[c];
         assert(offset >= 0);

         if( rowSparsity )
            rowSparsity[arrayLineOffset + col] = 1;

         rowsArrayDense[offset + col] = M[c];
      }
   }
}

void SparseStorage::getLinkVarsNnz(std::vector<int>& vec) const
{
   assert(int(vec.size()) == n);

   for( int i = 0; i < len; i++ )
   {
      const int col = jcolM[i];
      assert(col < n);

      vec[col]++;
   }
}

void SparseStorage::atPutSpRow( int row, double A[], int lenA,
				    int jcolA[], int& info )
{
  int ik;
  int ka    = lenA - 1;
  int km_f  = krowM[row + 1] - 1;
  int km    = km_f;
  int km_s  = krowM[row];
  int count = 0;

  assert( row >= 0 && row < m );

  while( ka >= 0 ) {
    if ( km < km_s ) {
      // There are no more elements in M. All the rest of A must be
      // inserted.
      count += ka + 1;
      break;
    } else if ( jcolM[km] == jcolA[ka] ) {
      // The element in A will replace an element in M
      km--; ka--;
    } else if ( jcolM[km] > jcolA[ka] ) {
      assert( jcolA[ka] >= 0 );
      // This element is in M but not in A
      km--;
    } else {
      // The element is in A, but not in M and so must be inserted.
      assert( jcolA[ka] < n );
      ka--;
      count++;
    }
  }

  if ( count > 0 ) {
    this->shiftRows( row + 1, count, info );
  } else {
    info = 0;
  }

  if ( 0 == info ) {
    ka    = lenA - 1;
    km    = km_f;
    ik    = krowM[row + 1] - 1;

    while ( ka >= 0 ) {
      if ( km < km_s ) {
	// There are no more elements in M. All the rest of A must be
	// inserted.
	for( ; ka >= 0; ka--, ik-- ) {
	  jcolM[ik] = jcolA[ka];
	  assert( jcolM[ik] >= 0 && jcolM[ik] < n );
	  M[ik]     = A[ka];
	}
	break;
      } else if ( jcolM[km] == jcolA[ka] ) {
	// The element in A will replace an element in M
	jcolM[ik] = jcolM[km];
	M[ik]     = A[ka];
	km--; ka--; ik--;
      } else if ( jcolM[km] > jcolA[ka] ) {
	// This element is in M but not in A
	jcolM[ik] = jcolM[km];
	M[ik]     = M[km];
	km--; ik--;
      } else {
	// The element is in A, but not in M.
	jcolM[ik] = jcolA[ka];
	assert( jcolM[ik] >= 0 && jcolM[ik] < n );
	M[ik]     = A[ka];
		
	ka--; ik--;
      }
    }
  }
}

void SparseStorage:: atPutDense( int row, int col, double * A, int lda,
				     int rowExtent, int colExtent )
{ 
  int info, count;
  int i, km_f, km, ka, k;

  assert( row >= 0 && row + rowExtent <= m );
  assert( col >= 0 && col + colExtent <= n );

  for ( i = row; i < row + rowExtent; i++) {
    // Loop over all rows in range.
    km_f = krowM[i + 1] - 1;
    km   = km_f;
    ka   = colExtent - 1;
    count = 0;
    while( km >= krowM[i] || ka >= 0 ) {
      // The current row in M and the current row in A are not
      // both empty
      if ( ka < 0 ) {
	// The current row of A is empty. Insert an element from M
	km--;
      } else if ( km < krowM[i] ) {
	// The current row of M is empty. Insert an element from A.
	if ( A[(i - row) * lda + ka] == 0 ) {
	  // The current element of A is zero. Skip it.
	  ka--;
	} else {
	  count++; ka --;
	}
      } else if ( ka + col > jcolM[km] ) {
	// The current element in A comes first.
	if ( A[(i - row) * lda + ka] == 0 ) {
	  // The current element of A is zero. Skip it.
	  ka--;
	} else {
	  count++; ka--;
	}
      } else if ( ka + col < jcolM[km] ) {
	// The current element in M comes first.
	km--;
      } else {
	// The current element in M is overwritten by the element in A.
	km--; ka--;
      }
    } // End while the current row...are not both empty
    if ( count > 0 ) {
      this->shiftRows( i + 1, count, info );
    } else {
      info = 0;
    }
    if ( info != 0 ) { 
      cout << "bing\n"; return;
    }
	
    km = km_f;
    k  = krowM[i + 1] - 1;
    ka = colExtent - 1;

    while( km >= krowM[i] || ka >= 0 ) {
      // The current row in M and the current row in A are not
      // both empty
      if ( ka < 0 ) {
	// The current row of A is empty. Insert an elemnt from M
	M[k]     = M[km];
	jcolM[k] = jcolM[km];
	k--; km--;
      } else if ( km < krowM[i] ) {
	// The current row of M is empty. Insert an element from A.
	if ( A[ (i - row) * lda + ka] == 0 ) {
	  // The current element of A is zero. Skip it.
	  ka--;
	} else {
	  M[k]     = A[(i - row) * lda + ka];
	  jcolM[k] = ka + col;
	  k--; ka --;
	}
      } else if ( ka + col > jcolM[km] ) {
	// The current element in A comes first.
	if ( A[ (i - row) * lda + ka] == 0 ) {
	  // The current element of A is zero. Skip it.
	  ka--;
	} else {
	  M[k]     = A[(i - row) * lda + ka];
	  jcolM[k] = ka + col;
	  k--; ka--;
	}
      } else if ( ka + col < jcolM[km] ) {
	// The current element in M comes first.
	M[k]      = M[km];
	jcolM[k]  = jcolM[km];
	k--; km--;
      } else {
	// The current element in M is overwritten by the element in A.
	M[k]       = A[(i - row) * lda + ka];
	jcolM[k]   = ka + col;  
	k--; km--; ka--;
      }
    } // End while the current row...are not both empty
  } // End loop over all rows in range.
}

void SparseStorage::shiftRows( int row, int shift, int& info )
{
  if ( shift == 0 ) {
    info = 0;
  } else if ( krowM[m] + shift > len ) {
    // Insufficient space
    info = krowM[m] + shift - len;
  } else {
    // We perform the copy
    info = 0;
    int lcopy = krowM[m] - krowM[row];
    if ( lcopy > 0 ) {
      // There is anything to copy
      // As a consequence of lcopy > 0, col !== n
      memmove( &jcolM[ krowM[row] + shift ], 
	       &jcolM[ krowM[row] ], lcopy * sizeof(int) );
      memmove( &M[ krowM[row] + shift ],
	       &M[ krowM[row] ], lcopy * sizeof(double) );
      int i;
      for ( i = row; i <= m; i++ ) {
	krowM[i] += shift;
      }    
    } else {
      // Still adjust the starts of the rows
      int i;
      int rowStart = krowM[m] + shift; 
      for ( i = row; i <= m; i++ ) {
	krowM[i] = rowStart;
      }    
    }
  } // end else we perform the copy
}

void SparseStorage::putSparseTriple( int irow[], int lenA,
					 int jcol[], double A[], 
					 int& info )
{
  if( len < lenA ) {
    info = 1;
    return;
  }
  info = 0;
  int i, k;
  krowM[0] = 0;
  i = 0;
  for( k = 0; k < lenA; k++ ) {
    while( i < irow[k] ) {
      i++;
      krowM[i] = k;
    }
    // now i == irow[k] because irow is sorted
    jcolM[k] = jcol[k];
    M[k]     = A[k];
  }
  for( i++; i <= m; i++ ) {
    krowM[i] = lenA;
  }
}

void SparseStorage::fromGetSpRow( int row, int col,
				      double A[], int lenA, int jcolA[],
				      int& nnz,
				      int colExtent, int& info )
{
  assert( col >= 0 && col < n );
  assert( row >= 0 && row < m );
  assert( col + colExtent <= n );
  int km, colm;
  int ka = 0;
  int lastCol = col + colExtent - 1;

  info = 0;
  
  for ( km = krowM[row]; km < krowM[row+1]; km++ ) {
    colm = jcolM[km];
    if ( colm >= col ) {
      if ( colm <= lastCol ) {
	if( ka < lenA ) {
	  A[ka]     = M[km];
	  jcolA[ka] = colm;
	  ka++;
	} else {
	  // Count the number of aditional elements needed in A
	  info++;
	}
      } else {
	break;
      }
    }
  }
  nnz = ka;
}

void SparseStorage::writeToStream(ostream& out) const
{
  int i, k;

  for( i = 0; i < m; i++ ) {
    for ( k = krowM[i]; k < krowM[i+1]; k++ ) {
      out << i << '\t' << jcolM[k] << '\t' << M[k] << endl;
    }
  }
}

void SparseStorage::writeToStreamDense( ostream& out) const
{
   int i, k;
   //todo: instead of \t, use length of longest value in M

   for( i = 0; i < m; i++ )
   { // Row i
      int j = 0; // Column j
      for( k = krowM[i]; k < krowM[i + 1]; k++ )
      {

#ifndef NDEBUG
         // assert that columns are ordered
         assert(k == krowM[i] || jcolM[k - 1] < jcolM[k]);
#endif

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

void SparseStorage::writeToStreamDenseRow( stringstream& out, int rowidx) const
{
   int j = 0; // Column j
   for( int k = krowM[rowidx]; k < krowM[rowidx + 1]; k++ )
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

void indexedLexSort( int first[], int n, int swapFirst,
		     int second[], int swapSecond, int index[] )
{
  int fi, se, j, k, kinc, inc, ktemp;
  const int nincs = 12;
  const int incs[]  = {1, 5, 19, 41, 109, 209, 505,
		       929, 2161, 3905, 8929, 16001};
  
  kinc = 0;
  for ( k = 0; k < nincs; k++ ) {
    kinc = k;
    if ( incs[kinc] > n/2 ) {
      kinc--;
      break;
    }
  }
  // incs[kinc] is the greatest value in the sequence that is also less
  // than n/2.

  //for( k = 0; k < n; k++ ) index[k] = k;

  for( ; kinc >= 0; kinc-- ) {
    // Loop over all increments
    inc = incs[kinc];

    if ( !swapFirst && !swapSecond ) {
      for ( k = inc; k < n; k++ ) {
	// loop over all subarrays defined by the current increment
	ktemp = index[k];
	fi = first[ ktemp ];
	se = second[ ktemp ];
	// Insert element k into the sorted subarray
	for( j = k; j >= inc; j -= inc ) {
	  // Loop over the elements in the current subarray
	  if ( fi < first[ index[j - inc] ] ||
	       ( fi == first[ index[j - inc] ] &&
		 se < second[ index[j - inc] ] ) ) {
	    // Swap elements j and j - inc, implicitly use the fact
	    // that ktemp hold element j to avoid having to assign to 
	    // element j - inc
	    index[j]        = index[j - inc];
	  } else {
	    // There are no more elements in this sorted subarray which
	    // are less than element j
	    break;
	  }
	} // End loop over the elements in the current subarray
	// Move index[j] out of temporary storage
	index[j] = ktemp;
	// The element has been inserted into the subarray.
      } // End loop over all subarrays defined by the current increment
    } else if ( swapSecond && !swapFirst ) {
      for ( k = inc; k < n; k++ ) {
	ktemp = index[k];
	fi = first[ ktemp ];
	se = second[ k ];
	for( j = k; j >= inc; j -= inc ) {
	  if ( fi < first[ index[j - inc] ] ||
	       ( fi == first[ index[j - inc] ] &&
		 se < second[ j - inc ] ) ) {
	    index[j]        = index[j - inc];
	    second[j]       = second[j - inc];
	  } else {
	    break;
	  }
	} 
	index[j]   = ktemp;
	second[j]  = se;
      } 
    } else if ( swapFirst  && !swapSecond ) {
      for ( k = inc; k < n; k++ ) {
	ktemp = index[k];
	fi = first[ k ];
	se = second[ ktemp ];
	for( j = k; j >= inc; j -= inc ) {
	  if ( fi < first[j - inc] ||
	       ( fi == first[j - inc] &&
		 se < second[ index[j - inc ] ]) ) {
	    index[j]        = index[j - inc];
	    first[j]        = first[j - inc];
	  } else {
	    break;
	  }
	} 
	index[j]   = ktemp;
	first[j]   = fi;
      } 
    } else { // Swap both
      for ( k = inc; k < n; k++ ) {
	ktemp = index[k];
	fi = first[ k ];
	se = second[ k];
	for( j = k; j >= inc; j -= inc ) {
	  if ( fi < first[j - inc] ||
	       ( fi == first[j - inc] &&
		 se < second[ j - inc ]) ) {
	    index[j]        = index[j - inc];
	    first[j]        = first[j - inc];
	    second[j]       = second[j - inc];
	  } else {
	    break;
	  }
	} 
	index[j]   = ktemp;
	first[j]   = fi;
	second[j]  = se;
      } 
    }
  } // End loop over all increments
}


void SparseStorage::mult( double beta,  double y[], int incy,
			      double alpha, const double x[], int incx ) const
{
  int i, j, k;
  double temp;
  for( i = 0; i < m; i++ ) {
    temp = 0;
    for( k = krowM[i]; k < krowM[i+1]; k++ ) {
      j = jcolM[k];
      temp += M[k] * x[j * incx];
#ifndef NDEBUG
      assert(j < n);
#endif
    }
    y[i * incy] = beta * y[i * incy] + alpha * temp;
  }
}

void SparseStorage::transMult( double beta,  double y[], int incy,
			       double alpha, const double x[], int incx ) const
{
  int i, j, k;
  if(beta!=1.0) {
    for( j = 0; j < n; j++ ) {
      y[j * incy] *= beta;
    }
  }
  for( i = 0; i < m; i++ ) {
    for( k = krowM[i]; k < krowM[i+1]; k++ ) {
      j = jcolM[k];
#ifndef NDEBUG
      assert(j < n);
#endif
      y[j * incy] += alpha * M[k] * x[i * incx];
    }
  }
}

/* Y <- alpha* M^T X + beta*Y
 * Computes only the elements in Y that are lower triangular
 * elements in the larger matrix (Y contains a subset of the ny columns of 
 * this matrix, starting at colStart)
 *
 * X  - dense matrix with columns in continuous memory**
 * that is, along *rows* in C-style, but treated as columns
 * Y - dense matrix with rows  in continuous memory that is
 * a subset of ny columns of a larger dense matrix. The first column
 * in Y is the colStart column in the larger matrix
 * M - is 'this'
 */
void SparseStorage::transMultMatLower( double beta,  double* Y, int ny, int ldy,
				       double alpha, double *X, int ldx, int colStart)
{

  //Note: Y[j+ldy*v] is Y(j,v)

  int i, j, k;
  if(beta!=1.0) {
    for( int v = 0; v < ny; v++) {
      int ldyv=ldy*v;
      int startrow=v+colStart; assert(false);
      for( j = startrow; j<n; j++ ) {
	Y[j +ldyv] *= beta;
      }
    }
  }
  for( i = 0; i < m; i++ ) {
    for( k = krowM[i]; k < krowM[i+1]; k++ ) {
      j = jcolM[k];
#ifndef NDEBUG
      assert(j<n);
#endif
      //int endcol=j+colStart;
      for (int v = j+colStart; v<ny; v++) { 
	Y[j+v*ldy] += alpha * M[k] * X[i+v*ldx];
      }
    }
  }
}


// Y <- alpha*M^TX + beta*Y
// X dense matrix with columns in continuous memory**
// that is, along *rows* in C-style, but treated as columns
void SparseStorage::transMultMat( double beta,  double* Y, int ny, int ldy,
						 double alpha, double *X, int ldx)
{
  int i, j, k;
  if(beta!=1.0) {
    for( int v = 0; v < ny; v++)
      for( j = 0; j < n; j++ ) {
	Y[j +ldy*v] *= beta;
      }
  }
  for( i = 0; i < m; i++ ) {
    for( k = krowM[i]; k < krowM[i+1]; k++ ) {
      j = jcolM[k];
#ifndef NDEBUG
      assert(j<n);
#endif
      for (int v = 0; v<ny; v++) { 
	Y[j+v*ldy] += alpha * M[k] * X[i+v*ldx];
      }
    }
  }
}

// adds this * x to the part of row yrow of y that is in the upper half of y
// y is assumed to be symmetric and sorted!
void SparseStorage::multMatSymUpper( double beta, SparseStorage& y,
              double alpha, double x[], int yrow, int ycolstart ) const
{
   assert(yrow >= 0 && yrow < y.m);
   assert(ycolstart >= 0 && ycolstart < y.n);
   assert(y.n == y.m);
   assert(y.n >= m + ycolstart);

   int* const krowM_y = y.krowM;
   int* const jcolM_y = y.jcolM;
   double* const M_y = y.M;

   // assert that yrow is sorted
#ifndef NDEBUG
   for( int ci = krowM_y[yrow] + 1; ci < krowM_y[yrow + 1]; ci++ )
      assert( jcolM_y[ci - 1] < jcolM_y[ci] );
#endif

   // scale row yrow
   if( beta != 1.0 )
      for( int c_y = krowM_y[yrow]; c_y != krowM_y[yrow + 1]; c_y++ )
         M_y[c_y] *= beta;

   // add this * x to yrow (and exploit that y is symmetric)
   int c_y = krowM_y[yrow];
   for( int r = 0; r < m; r++ )
   {
      double yrx;
      const int colplace_y = r + ycolstart; // the column in y where to place yrx

      // not in upper half of y?
      if( colplace_y < yrow )
         continue;

      // compute y_(r,.) * x
      yrx = 0.0;

      for( int c = krowM[r]; c != krowM[r + 1]; c++ )
      {
         const int col = jcolM[c];
         yrx += x[col] * M[c];
      }

      if( PIPSisZero(yrx) )
         continue;

      yrx *= alpha;

      for( ; c_y != krowM_y[yrow + 1]; c_y++ )
      {
         const int col_y = jcolM_y[c_y];
         if( col_y == colplace_y )
         {
            M_y[c_y] += yrx;
            break;
         }
      }
      assert(c_y != krowM_y[yrow + 1]);
   }
}

// y only contains the elements starting at firstrow
// this is useful in the symmetric reduce, because we can
// avoid a temporary array, and also save some computation
void SparseStorage::transMultLower( double beta,  double y[],
			       double alpha, double x[], int firstrow )
{
  int i, j, k;
  if(beta!=1.0) {
    for( j = 0; j < n; j++ ) {
      y[j] *= beta;
    }
  }
  for( i = 0; i < m; i++ ) {
    for( k = krowM[i]; k < krowM[i+1]; k++ ) {
      j = jcolM[k];
      if (j < firstrow) continue;
#ifndef NDEBUG
      assert(j<n);
#endif
      // nonzero element at (i,j)
      // multiply with x[i] and add to y[j]
      y[j - firstrow] += alpha * M[k] * x[i];
    }
  }
}

#ifndef MIN
#define MIN(a,b) ( (a > b) ? b : a)
#endif

// beta = 1.0
// Y <- alpha*M^TX + Y
void SparseStorage::transMultMatLower( double* Y, int ny, int firstrow,
						 double alpha, double *X, int ldx)
{
  int i, j, k, v;
	// map from col number to offset in Y
	int *offset = new int[ny];
	offset[0] = 0;
	for ( k = 1; k < ny; k++) {
		offset[k] = offset[k-1] + (n-firstrow-k+1);
	}
  for( i = 0; i < m; i++ ) {
    for( k = krowM[i]; k < krowM[i+1]; k++ ) {
      j = jcolM[k];
#ifndef NDEBUG
      assert(j<n);
#endif
      // better to parallelize this loop over the columns of X
      // because no synchronization needed
      int end = MIN(ny, j-firstrow+1);
      #ifdef _OPENMP
      #pragma omp parallel for schedule(dynamic,5)
      #endif
      for (v = 0; v < end; v++) {
				Y[offset[v]+(j-(firstrow+v))] += alpha * M[k] * X[i+v*ldx];
			}
    }
  }
	delete [] offset;
	}


void doubleLexSort( int first[], int n,
		    int second[], double data[] );

void SparseStorage::symmetrize( int& info)
{
  int i, k, ku;

  int nnz = krowM[m];
  int * irowM = new int[ 2 * nnz ];
  
  info = 0;

  ku = nnz;
  for ( i = 0; i < m; i++ ) {
    // For all rows
    for( k = krowM[i]; k < krowM[i + 1]; k++ ) {
      // Loop over elements of row i
      irowM[k] = i;
      if ( i != jcolM[k] ) {
	// Not a diagonal element
	if ( ku >= len ) {
	  info++;
	} else {
	  // Add the element transpose to the scrambled matrix
	  irowM[ku] = jcolM[k];
	  jcolM[ku] = i;
	  M[ku]     = M[k];
	  ku++;
	}
      } // End not a diagonal element
    } // End loop over elements of column j
  } // End for all columns
  if( info != 0 ) {
     delete [] irowM;
     return;
  }

  nnz = ku;

  doubleLexSort( irowM, nnz, jcolM, M );
  i = 0;
  krowM[0] = 0;
  for( k = 0; k < nnz; k++ ) {
    for( ; i < irowM[k]; i++ ) {
      krowM[i + 1] = k;
    }
  }
  for( i++; i <= m; i++ ) {
    krowM[i] = nnz;
  }
  delete [] irowM;
}

void SparseStorage::getTransposePat( int row, int col,
				     int rowExtent, int colExtent,
				     int kpat[], int kcolM[], int irowM[] )
{
  int i, j, k, kin, kout;
  const int dontPermuteCols = 0, doPermuteRows = 1;
  
  kout = 0;
  for( i = row; i < row + rowExtent; i++ ) {
    kin = krowM[i];
    while( kin < krowM[i+1] && jcolM[kin] < col ) {
      kin++;
    }
	
    for( ; kin < krowM[i+1] && jcolM[kin] < col + colExtent; kin++ ) {
      irowM[kout] = i;
      kpat[kout] = kin;
      kout++;
    }
  }
  
  indexedLexSort( jcolM, kout, dontPermuteCols, irowM, 
		  doPermuteRows, kpat );

  j = 0;
  kcolM[0] = 0;
  for ( k = 0; k < kout; k++ ) {
    for( ; j + col < jcolM[ kpat[k] ]; j++ ) {
      assert( j + 1 < colExtent );
      kcolM[j + 1] = k;
    }
  }
  for( j++; j <= colExtent; j++ ) {
    kcolM[j] = kout;
  }
}

void SparseStorage::getFromPat( double data[], int ldata, int kpat[] )
{
  int k;
  for ( k = 0; k < ldata; k++ ) {
    data[k] = M[ kpat[k] ];
  }
}

void SparseStorage::randomize( double alpha, double beta, double * seed )
{
  int  i, k, NN, chosen, icurrent;
  double r;

  double scale = beta - alpha;
  double shift = alpha/scale;

  double drand( double * );

  // Knuth's algorithm for choosing length elements out of NN elts.
  NN        = m * n;
  int length = (len <= NN) ? len : NN;
  chosen    = 0;
  icurrent  = 0;
  krowM[0]  = 0;
  for ( k = 0; k < NN; k++ ) {
    r = drand( seed );
	
    if( (NN - k) * r < length - chosen ) {
      jcolM[chosen] = k % n;
      i             = k / n;

      if ( i > icurrent ) {
	for ( ; icurrent < i; icurrent++ ) {
	  krowM[icurrent + 1] = chosen;
	}
      }
      M[chosen]     = scale * (drand(seed) + shift);
      chosen++;
    }  	
  }
  for ( ; icurrent < m; icurrent++ ) {
    krowM[icurrent + 1] = length;
  }

  assert( chosen == length );

}

double SparseStorage::abmaxnorm() const
{
  double norm = 0.0;
  int nnz = this->numberOfNonZeros();
  
  int i;
  double fabsMi;
  for( i = 0; i < nnz; i++ ) {
    fabsMi = fabs( M[i] );
    if ( fabsMi > norm ) norm = fabsMi;
  }
  return norm;
}


void SparseStorage::atPutDiagonal( int idiag, OoqpVector& vvec )
{
  SimpleVector & v = dynamic_cast<SimpleVector &>(vvec);
  this->atPutDiagonal( idiag, &v[0], 1, v.n );
}

void SparseStorage::atPutDiagonal( int idiag,
				   double x[], int incx, int extent )
{
  int i;
  int info;
  for( i = idiag; i < idiag + extent; i++ ) { // Loop over elts to be put
    // Search for the diagonal elt.
    int lastk;
    lastk = krowM[i+1];
    for( int k = krowM[i]; k < lastk; k++ ) { // Loop over all elts in row
      if ( i >= jcolM[k] ) { // Found or past the diagonal
	if( i == jcolM[k] ) { // Found it, overwrite it.
	  M[k] = x[incx*(i - idiag)];
	} else { // Didn't find it, so insert it
	  this->atPutSpRow( i, &x[incx*(i - idiag)], 1, &i, info );
	}
	// Either way, bug out of the loop
	break;
      } // end if found or past the diagonal
    } // end loop over all elts
  } // end loop over elts to be put
}


/*void SparseStorage::matTransDMultMat(double* d, int** krowAtDA, int** jcolAtDA, double** AtDA)
{
  int nnzAtA = 0;
  int nnz = krowM[m];
  int ind,j;
  double val;
  
  /////////////////////////////////////////////////////
  // form the transpose 
  /////////////////////////////////////////////////////

  //build the transpose
  int* krowMt = new int[n+1];
  int* jcolMt = new int[nnz];
  double* Mt  = new double[nnz];

  //cummulative sum to find the number of elements in each row of At, ie column of A.
  int* w=new int[n];
  for(int i=0; i<n;   i++) w[i]=0;
  for(int i=0; i<nnz; i++) w[ jcolM[i] ]++;
  
  krowMt[0]=0; double sum=0.0;
  for(int i=1; i<=n; i++) {
    krowMt[i] = krowMt[i-1]+w[i-1];
    w[i-1] = krowMt[i-1];
  }

  //populate jcolMt and Mt now
  for(int i=0; i<m; i++) {
    for(int j=krowM[i]; j<krowM[i+1]; j++) {
      int ind =  w[jcolM[j]];
      jcolMt[ ind ] = i;
      Mt    [ ind ] = M[j];
      w[jcolM[j]] = ind+1; //w[jcolM[j]]++
    }
  }
  //we have the transpose
  delete[] w;

  ////////////////////////////////////////////////////
  // count the number of entries in the result AtA  //
  ////////////////////////////////////////////////////
  for(int i=0; i<n; i++) {
    for(int k=krowMt[i]; k<krowMt[i+1]; k++) {
      j = jcolMt[k]; //At[i,j] is nz

      //add the nonzero pattern of row i of A to AtA
      for(int l=krowM[j]; l<krowM[j+1]; l++)
	nnzAtA++;
      
    }
  }
  assert(nnzAtA>=0); //overflow?!?

  ////////////////////////////////////////////////////
  // alocate AtDA
  ///////////////////////////////////////////////////
  *krowAtDA = new int[n+1];
  *jcolAtDA = new int[nnzAtA];
  *AtDA     = new double[nnzAtA];

  ////////////////////////////////////////////////////
  // AtDA = A'*D*A
  ////////////////////////////////////////////////////
  double* W = new double[n];
  for(int it=0; it<n; it++) W[it] = 0.0;

  nnzAtA=0;
  for(int i=0; i<n; i++) {
    //start row i of AtDA
    krowAtDA[i]=nnzAtA;
    
    for(int pt=krowMt[i]; pt<krowMt[i+1]; pt++) { 
      k = jcolMt[pt]; //At[i,k] is non-zero
      val = Mt[pt]*d[k];

      //iterate the row k of A and scatter the values into W
      for(int p=krowM[k]; p<krowM[k+1]; p++) {
	j = jcolM[p];
	//we have A[k,j]

	jcolAtDA[nnzAtA++]=j;
	W[j] += M[p]*val;
      }
    }
    //gather the values into the i-th row AtDA
    for(int p=krowAtDA[i]; p<nnzAtA; p++) {
      j=jcolAtDA[p];
      AtDA[p]=W[j];
      W[j]=0.0;
    }
    
  }
  delete[] W;

}

*/

void SparseStorage::transpose(int* krowMt, int* jcolMt, double* Mt) const
{
  int ind, pend;//,ptend;
  /////////////////////////////////////////////////////
  // form the transpose 
  /////////////////////////////////////////////////////
  const int nnz = krowM[m];

  //cumulative sum to find the number of elements in each row of At, ie column of A.
  int* w=new int[n];
  for(int i=0; i<n;   i++) w[i]=0;
  for(int i=0; i<nnz; i++) w[ jcolM[i] ]++;

  krowMt[0]=0; //double sum=0.0;
  for(int i=1; i<=n; i++) {
    krowMt[i] = krowMt[i-1]+w[i-1];
    w[i-1] = krowMt[i-1];
  }

  //populate jcolMt and Mt now
  for(int i=0; i<m; i++) {
    pend=krowM[i+1];
    for(int j=krowM[i]; j<pend; j++) {
      ind =  w[jcolM[j]];
      jcolMt[ ind ] = i;
      Mt    [ ind ] = M[j];
      w[jcolM[j]] = ind+1; //w[jcolM[j]]++
    }
  }
  //we have the transpose
  delete[] w;
}

void SparseStorage::clear()
{
   for( int i = 0; i < len; i++ )
      M[i] = 0.0;
}

void SparseStorage::matTransDSymbMultMat(double* d,
					 int* krowMt,  int* jcolMt,  double* Mt, 
					 int** krowAtDA, int** jcolAtDA, double** AtDA)
					 
{
  int k,j, nnzAtA=0, pend,ptend;

  ////////////////////////////////////////////////////
  // count the number of entries in the result AtA  //
  ////////////////////////////////////////////////////
  char* flag=new char[n];
  for(int i=0; i<n; i++) {
    memset(flag, 0, n);

    ptend=krowMt[i+1];
    for(int pt=krowMt[i]; pt<ptend; pt++) {
      k = jcolMt[pt]; //At[i,j] is nz

      //add the nonzero pattern of row i of A to AtA
      pend=krowM[k+1];
      for(int p=krowM[k]; p<pend; p++) {
	j = jcolM[p];

	if(flag[j]==0) {
	  nnzAtA++;
	  flag[j]=1;
	}
      }
    }
  }
  assert(nnzAtA>=0); //overflow?!?

  delete[] flag;
  ////////////////////////////////////////////////////
  // alocate AtDA
  ///////////////////////////////////////////////////
  *krowAtDA = new int[n+1];
  *jcolAtDA = new int[nnzAtA];
  *AtDA     = new double[nnzAtA];
}


void SparseStorage::matTransDMultMat(double* d, 
				     int* krowMt, int* jcolMt, double* Mt,
				     int* krowAtDA, int* jcolAtDA, double* AtDA)
{
  int k,j, pend,ptend; double val;
  ////////////////////////////////////////////////////
  // AtDA = A'*D*A
  ////////////////////////////////////////////////////
  double* W = new double[n];
  char* flag=new char[n];

  for(int it=0; it<n; it++) W[it] = 0.0;

  int nnzAtA=0;
  for(int i=0; i<n; i++) {
    memset(flag, 0, n);

    //start row i of AtDA
    krowAtDA[i]=nnzAtA;
    
    ptend=krowMt[i+1];
    for(int pt=krowMt[i]; pt<ptend; pt++) { 
      k = jcolMt[pt]; //At[i,k] is non-zero
      val = Mt[pt]*d[k];

      //iterate the row k of A and scatter the values into W
      pend=krowM[k+1];
      for(int p=krowM[k]; p<pend; p++) {
	j = jcolM[p];
	//we have A[k,j]
	if(flag[j]==0) {
	  jcolAtDA[nnzAtA++]=j;
	  flag[j]=1;
	}
	
	W[j] += (M[p]*val);
      }
    }
    //gather the values into the i-th row AtDA
    for(int p=krowAtDA[i]; p<nnzAtA; p++) {
      j=jcolAtDA[p];
      AtDA[p]=W[j];
      W[j]=0.0;
    }
  }
  krowAtDA[n] =nnzAtA;
  delete[] W;
  delete[] flag;
}

void SparseStorage::matTransDinvMultMat(double* d, 
				     int* krowMt, int* jcolMt, double* Mt,
				     int* krowAtDA, int* jcolAtDA, double* AtDA)
{
  int k,j, pend,ptend; double val;
  ////////////////////////////////////////////////////
  // AtDA = A'*D*A
  ////////////////////////////////////////////////////
  double* W = new double[n];
  char* flag=new char[n];

  for(int it=0; it<n; it++) W[it] = 0.0;

  int nnzAtA=0;
  for(int i=0; i<n; i++) {
    memset(flag, 0, n);

    //start row i of AtDA
    krowAtDA[i]=nnzAtA;
    
    ptend=krowMt[i+1];
    for(int pt=krowMt[i]; pt<ptend; pt++) { 
      k = jcolMt[pt]; //At[i,k] is non-zero
      val = Mt[pt]/d[k];

      //iterate the row k of A and scatter the values into W
      pend=krowM[k+1];
      for(int p=krowM[k]; p<pend; p++) {
	j = jcolM[p];
	//we have A[k,j]
	if(flag[j]==0) {
	  jcolAtDA[nnzAtA++]=j;
	  flag[j]=1;
	}
	
	W[j] += (M[p]*val);
      }
    }
    //gather the values into the i-th row AtDA
    for(int p=krowAtDA[i]; p<nnzAtA; p++) {
      j=jcolAtDA[p];
      AtDA[p]=W[j];
      W[j]=0.0;
    }
  }
  krowAtDA[n] =nnzAtA;
  delete[] W;
  delete[] flag;
}

void SparseStorage::reduceToLower()
{
  assert(m==n); //available only for square matrices
  int newNnz=0; // the new number of nonzeros
  //get the nnz of the reduced matrix
  for(int i=0; i<n; i++) {
    int jend=krowM[i+1];
    for(int j=krowM[i]; j<jend; j++)
      if(jcolM[j]<=i) newNnz++;
  }
  
  int* newjcolM = new int[newNnz];
  double* newdM = new double[newNnz];
  int it=0;

  for(int i=0; i<n; i++) {
    int jend=krowM[i+1]; int jstart=krowM[i];
    krowM[i]=it;
    for(int j=jstart; j<jend; j++) {
      if(jcolM[j]<=i) {
	newjcolM[it] = jcolM[j];
	newdM[it] = M[j];
	it++;
      }
    }
  }
  assert(it==newNnz);
  krowM[n] = it;
  
  delete[] jcolM; delete [] M;
  jcolM = newjcolM;
  M = newdM;
  //for(int k=krowC
}

void SparseStorage::dump(const string& filename)
{
  ofstream fd(filename.c_str());
  fd << scientific; 
  fd.precision(16);
  fd << m << " " << n << " " << len << endl;

  int i;
  for (i = 0; i <= m; i++) {
    fd << krowM[i] << " ";
  }
  fd << endl;
  for (i = 0; i < len; i++) {
    fd << jcolM[i] << " ";
  }
  fd << endl;
  for (i = 0; i < len; i++) {
    fd << M[i] << " ";
  }
}

void SparseStorage::deleteEmptyRowsCols(const int* nnzRowVec, const int* nnzColVec)
{
   assert(nnzRowVec != nullptr && nnzColVec != nullptr);

   int m_new = 0;
   int n_new = 0;

   int* rowsmap = new int[m];

   for( int i = 0; i < m; i++ )
      if( nnzRowVec[i] != 0.0 )
         rowsmap[i] = m_new++;

   int* colsmap = new int[n];

   for( int i = 0; i < n; i++ )
      if( nnzColVec[i] != 0.0 )
         colsmap[i] = n_new++;

   assert(krowM[0] == 0 && "method not fully supported");

#if 0
   int nnz = 0;
   int rowcount = 0;
   int len_new = 0;

   for( int r = 0; r < m; r++ )
   {
      const int nnzRow = nnzRowVec[r];

      if( nnzRow == 0.0 )
      {
         nnz += krowM[r + 1] - krowM[r];
         continue;
      }
      krowM[rowcount] = len_new;

      for( int j = krowM[r]; j < krowM[r + 1]; j++ )
      {
         jcolM[len_new] = jcolM[colsmap[nnz]];
         M[len_new++] = M[nnz++];
      }
      rowcount++;
   }


   for( int r = rowcount; r <= m; r++ )
      krowM[r] = len_new;

   assert(nnz == len);
#endif

   m = m_new;
   n = n_new;

   delete[] colsmap;
   delete[] rowsmap;
}


void SparseStorage::getSparseTriplet_c2fortran(int*& irn, int*& jcn, double*& val) const
{
   int count = 0;
   assert(len > 0);
   assert(!irn && !jcn && !val);

   irn = new int[len];
   jcn = new int[len];
   val = new double[len];

   for( int r = 0; r < m; r++ )
   {
      for( int c = krowM[r]; c < krowM[r + 1]; c++ )
      {
         const int col = jcolM[c];
         const double value = M[c];

         irn[count] = r + 1;
         jcn[count] = col + 1;
         val[count] = value;

         count++;
      }
   }

   assert(count == len);
}


void SparseStorage::getSparseTriplet_fortran2fortran(int*& irn, int*& jcn, double*& val) const
{
   int count = 0;
   assert(len > 0);
   assert(!irn && !jcn && !val);
   assert(fortranIndexed());

   irn = new int[len];
   jcn = new int[len];
   val = new double[len];

   for( int r = 0; r < m; r++ )
   {
      for( int c = krowM[r] - 1; c < krowM[r + 1] - 1; c++ )
      {
         const int col = jcolM[c];
         const double value = M[c];

         irn[count] = r + 1;
         jcn[count] = col;
         val[count] = value;

         count++;
      }
   }

   assert(count == len);
}

void SparseStorage::deleteEmptyRows(int*& orgIndex)
{
   assert(!neverDeleteElts);
   assert(orgIndex == nullptr);

   int m_new = 0;

   // count non-empty rows
   for( int r = 0; r < m; r++ )
      if( krowM[r] != krowM[r + 1] )
         m_new++;

   int* krowM_new = new int[m_new + 1];
   orgIndex = new int[m_new + 1];

   krowM_new[0] = 0;
   m_new = 0;

   for( int r = 0; r < m; r++ )
      if( krowM[r] != krowM[r + 1] )
      {
         orgIndex[m_new] = r;
         krowM_new[++m_new] = krowM[r + 1];
      }

   assert(krowM_new[m_new] == len);
   m = m_new;

   delete[] krowM;
   krowM = krowM_new;
}


void SparseStorage::c2fortran()
{
   assert(krowM[0] == 0 && krowM[m] == len && !isFortranIndexed);

   for( int i = 0; i <= m; i++ )
      krowM[i]++;

   for( int i = 0; i < len; i++ )
      jcolM[i]++;

   isFortranIndexed = true;
}

void SparseStorage::fortran2c()
{
   assert(krowM[0] == 1 && krowM[m] == len + 1 && isFortranIndexed);

   for( int i = 0; i <= m; i++ )
      krowM[i]--;

   for( int i = 0; i < len; i++ )
      jcolM[i]--;

   isFortranIndexed = false;
}

bool SparseStorage::fortranIndexed() const
{
   return isFortranIndexed;
}

void SparseStorage::set2FortranIndexed()
{
   assert(krowM[0] == 1 && krowM[m] == len + 1);

   isFortranIndexed = true;
}

void SparseStorage::deleteZeroRowsColsSym(int*& new2orgIdx)
{
   assert(m == n);
   assert(!neverDeleteElts);
   assert(new2orgIdx == nullptr);
   assert(this->isValid());

   int* const offset = new int[m];

   for( int r = 0; r < m; r++ )
      offset[r] = 0;

   // mark rows (and columns) to be deleted
   for( int r = 0; r < m; r++ )
   {
      const int start = krowM[r];
      const int end = krowM[r + 1];

      if( start == end )
      {
         offset[r] = -1;
         continue;
      }

      int c = start;

      for( ; c < end; c++ )
         if( !PIPSisZero(M[c]) )
            break;

      // no non-zero found?
      if( c == end )
      {
         offset[r] = -1;
         continue;
      }

      offset[r] = -2;

      for( c = start; c < end; c++ )
      {
         const int col = jcolM[c];

         assert(offset[col] != 0);

         if( !PIPSisZero(M[c]) && offset[col] == -1 )
            offset[col] = -2;
      }
   }

   int rowDeletes = 0;
   int zeroEntryDeletes = 0;

   // count column offsets and entries to be deleted
   for( int r = 0; r < m; r++ )
   {
      const int start = krowM[r];
      const int end = krowM[r + 1];

      // row deleted?
      if( offset[r] == -1 )
      {
         zeroEntryDeletes += end - start;
         rowDeletes++;
         continue;
      }

      for( int c = start; c < end; c++ )
      {
         const int col = jcolM[c];
         assert(col < m);

         if( offset[col] == -1 )
         {
            assert(PIPSisZero(M[c]));
            zeroEntryDeletes++;
         }
      }

      assert(offset[r] == -2);
      offset[r] = rowDeletes;
   }

   const int m_new = m - rowDeletes;
   const int len_new = len - zeroEntryDeletes;
   assert(len_new >= 0 && m_new >= 0);

   new2orgIdx = new int[m_new];
   int* const krowM_new = new int[m_new + 1];
   int* const jcolM_new = new int[len_new];
   double* const M_new = new double[len_new];
   int m_count = 0;
   int len_count = 0;

   // fill the new arrays
   krowM_new[0] = 0;
   for( int r = 0; r < m; r++ )
   {
      if( offset[r] == -1 )
         continue;

      for( int c = krowM[r]; c < krowM[r + 1]; c++ )
      {
         const int col = jcolM[c];
         if( offset[col] == -1 )
         {
            assert(PIPSisZero(M[c]));
            continue;
         }

         assert(col - offset[col] >= 0);

         jcolM_new[len_count] = col - offset[col];

         M_new[len_count++] = M[c];
      }

      new2orgIdx[m_count] = r;
      krowM_new[++m_count] = len_count;
      assert(krowM_new[m_count] > krowM_new[m_count - 1]);
   }

   assert(m_count == m_new);
   assert(len_count == len_new);

   delete[] krowM;
   delete[] jcolM;
   delete[] M;
   delete[] offset;

   m = m_new;
   n = m_new;
   len = len_new;
   krowM = krowM_new;
   jcolM = jcolM_new;
   M = M_new;

   assert(this->isValid());
}

void SparseStorage::addNnzPerRow(int* vec) const
{
   for( int r = 0; r < m; r++ )
      vec[r] += krowM[r + 1] - krowM[r];
}

void SparseStorage::addRowSums(double* vec) const
{
   for( int r = 0; r < m; r++ )
   {
      const int end = krowM[r + 1];
      for( int c = krowM[r]; c < end; c++ )
         vec[r] += std::fabs(M[c]);
   }
}

void SparseStorage::getRowMinVec(const double* colScaleVec, double* vec) const
{
   const bool coscale = (colScaleVec != nullptr);

   if( coscale )
   {
      for( int r = 0; r < m; r++ )
      {
         double minval = vec[r];
         assert(minval >= 0.0);

         for( int i = krowM[r]; i < krowM[r + 1]; i++ )
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

         for( int i = krowM[r]; i < krowM[r + 1]; i++ )
         {
            const double absval = std::abs(M[i]);

            if( absval < minval && absval > pips_eps )
               minval = absval;
         }
         vec[r] = minval;
      }
   }
}

void SparseStorage::getRowMaxVec(const double* colScaleVec, double* vec) const
{
   const bool coscale = (colScaleVec != nullptr);

   if( coscale )
   {
      for( int r = 0; r < m; r++ )
      {
         double maxval = vec[r];
         assert(maxval >= 0.0);

         for( int i = krowM[r]; i < krowM[r + 1]; i++ )
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

         for( int i = krowM[r]; i < krowM[r + 1]; i++ )
         {
            const double absval = std::abs(M[i]);

            if( absval > maxval )
               maxval = absval;
         }
         vec[r] = maxval;
      }
   }
}

void SparseStorage::getRowMinMaxVec(bool getMin, const double* colScaleVec, double* vec) const
{
   if( getMin )
      getRowMinVec(colScaleVec, vec);
   else
      getRowMaxVec(colScaleVec, vec);
}


void SparseStorage::permuteRows(const std::vector<unsigned int>& permvec)
{
   assert(permvec.size() == size_t(m));

   if( len == 0 )
      return;

   assert(m > 0 && n > 0);

   int* jcolM_new = new int[len];
   int* krowM_new = new int[m + 1];
   double* M_new = new double[len];

   int len_new = 0;
   krowM_new[0] = 0;

   for( int r = 0; r < m; ++r )
   {
      const unsigned int r_perm = permvec[r];
      const int rowlength = krowM[r_perm + 1] - krowM[r_perm];

      assert(r_perm < static_cast<unsigned int>(m) && rowlength >= 0);

      if( rowlength > 0 )
      {
         memcpy(jcolM_new + len_new, jcolM + krowM[r_perm], rowlength * sizeof(int));
         memcpy(M_new + len_new, M + krowM[r_perm], rowlength * sizeof(double));

         len_new += rowlength;
      }

      krowM_new[r + 1] = len_new;
   }

   assert(len_new == len);

   delete[] jcolM;
   delete[] krowM;
   delete[] M;

   jcolM = jcolM_new;
   krowM = krowM_new;
   M = M_new;
}


void SparseStorage::permuteCols(const std::vector<unsigned int>& permvec)
{
   assert(int(permvec.size()) == n);

   std::vector<unsigned int> permvec_rev(n, 0);

   int* indexvec = new int[n];
   int* bufferCol = new int[n];
   double* bufferM = new double[n];

   for( int i = 0; i < n; i++ )
   {
      assert(permvec[i] < unsigned(n));

      permvec_rev[permvec[i]] = i;
   }

   for( int r = 0; r < m; ++r )
   {
      const int row_start = krowM[r];
      const int row_end = krowM[r + 1];
      const int row_length = row_end - row_start;

      if( row_length == 0 )
         continue;

      for( int c = row_start; c < row_end; c++ )
      {
         const int col = jcolM[c];
         assert(col < n);

         jcolM[c] = permvec_rev[col];
      }

      for( int i = 0; i < row_length; i++ )
         indexvec[i] = i;

      std::sort(indexvec, indexvec + row_length, index_sort(jcolM + row_start, row_length) );

      for( int i = 0; i < row_length; i++ )
      {
         assert(indexvec[i] < row_length);

         bufferCol[i] = jcolM[row_start + indexvec[i]];
         bufferM[i] = M[row_start + indexvec[i]];
      }

      memcpy(jcolM + row_start, bufferCol, row_length * sizeof(int));
      memcpy(M + row_start, bufferM, row_length * sizeof(double));
   }

   delete[] indexvec;
   delete[] bufferM;
   delete[] bufferCol;
}


void SparseStorage::sortCols()
{
   int* indexvec = new int[n];
   int* bufferCol = new int[n];
   double* bufferM = new double[n];

   for( int r = 0; r < m; ++r )
   {
      const int row_start = krowM[r];
      const int row_end = krowM[r + 1];
      const int row_length = row_end - row_start;

      if( row_length == 0 )
         continue;

      for( int i = 0; i < row_length; i++ )
         indexvec[i] = i;

      std::sort(indexvec, indexvec + row_length, index_sort(jcolM + row_start, row_length));

      for( int i = 0; i < row_length; i++ )
      {
         assert(indexvec[i] < row_length);

         bufferCol[i] = jcolM[row_start + indexvec[i]];
         bufferM[i] = M[row_start + indexvec[i]];
      }

      memcpy(jcolM + row_start, bufferCol, row_length * sizeof(int));
      memcpy(M + row_start, bufferM, row_length * sizeof(double));
   }

   delete[] indexvec;
   delete[] bufferM;
   delete[] bufferCol;
}


/*
 * computes the full sparse matrix representation from a upper triangular symmetric sparse representation
 *
 * Must be square, the storage for the full representation will be allocated within the matrix and must be released later
 */
void SparseStorage::fullMatrixFromUpperTriangular(int*& rowPtrFull, int*& colIdxFull, double*& valuesFull) const
{
   assert(n == m);

   /* cout elems per row and assert upper triangular */
   std::vector<int> nelems(n, 0);
   for(int i = 0; i < n; ++i)
   {
      // diag elem
      for(int j = krowM[i]; j < krowM[i + 1]; ++j)
      {
         assert(jcolM[j] >= i);

         if(i == jcolM[j])
            nelems[i]++;
         else
         {
            nelems[jcolM[j]]++;
            nelems[i]++;
         }
      }
   }

   // fill rowptr array
   rowPtrFull = new int[n + 1];

   rowPtrFull[0] = 0;
   for(int i = 0; i < n; ++i)
      rowPtrFull[i + 1] = rowPtrFull[i] + nelems[i];

   colIdxFull = new int[rowPtrFull[n]];
   for( int i = 0; i < rowPtrFull[n]; ++i)
      colIdxFull[i] = -1;

   valuesFull = new double[rowPtrFull[n]];

   // fill in col and value
   for(int i = 0; i < n; ++i)
   {
      int rowstart = krowM[i];
      int rowend = krowM[i+1];

      for(int k = rowstart; k < rowend; ++k)
      {
         double value = M[k];
         int col = jcolM[k];

         int kk = rowPtrFull[i];
         int colfull = colIdxFull[kk];
         while(colfull != -1)
         {
            kk++;
            colfull = colIdxFull[kk];
         }

         colIdxFull[kk] = col;
         valuesFull[kk] = value;

         if(col != i )
         {
            kk = rowPtrFull[col];
            int colfull = colIdxFull[kk];

            while(colfull != -1)
            {
               kk++;
               colfull = colIdxFull[kk];
            }

            colIdxFull[kk] = i;
            valuesFull[kk] = value;
         }
      }
    }
}


// concatenate matrices
// if "diagonal", make a block diagonal matrix:
// [ A 0 ]
// [ 0 B ]
// if not, stack the matrices:
// [ A ]
// [ B ]
// diagonal will be symmetric if the input is
/*SparseStorage::SparseStorage(const vector<SparseStorage*> &blocks, bool diagonal)
{
  assert(blocks.size() > 0);
  m = n = len = 0;
  neverDeleteElts = 0;
  for (size_t i = 0; i < blocks.size(); i++) {
    m += blocks[i]->m;
    len += blocks[i]->len;
    if (diagonal) {
      n += blocks[i]->n;
    } else {
      assert( blocks[i]->n == blocks[0]->n );
    }
  }
  if (!diagonal) {
    n = blocks[0]->n;
  }

 
  M     = new double[len];
  jcolM = new int[len];
  krowM = new int[m+1];
  
  int curnnz = 0, curn = 0, curm = 0;

  for (size_t i = 0; i < blocks.size(); i++) {
    int lnnz = blocks[i]->len;
    int ln = blocks[i]->n;
    int lm = blocks[i]->m;

    memcpy(M+curnnz, blocks[i]->M, lnnz*sizeof(double));
    memcpy(jcolM+curnnz, blocks[i]->jcolM, lnnz*sizeof(int));
    memcpy(krowM+curm, blocks[i]->krowM, lm*sizeof(int)); // skips last element

    if (diagonal) {
      for (int j = 0; j < lnnz; j++) {
        jcolM[curnnz+j] += curn;
      }
    }
    for (int j = 0; j < lm; j++) {
      krowM[curm+j] += curnnz;
    }
    curnnz += lnnz;
    curn += ln;
    curm += lm;
  }
  assert(curm == m);
  krowM[m] = len;

  SparseStorage::instances++;
}
*/
