/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

/* 2015. Modified by Nai-Yuan Chiang for NLP*/

#include <cmath>
#include <cstring>
#include <cassert>
#include <cstdio>
#include <stdlib.h>

#include "SparseStorage.h"
#include "OoqpVector.h"
#include "SimpleVector.h"

#include <fstream>
#include <climits>

using namespace std;

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

  additiveDiag = new double[n];

  SparseStorage::instances++;

}

SparseStorage::SparseStorage( int m_, int n_, int len_,
			      int * krowM_, int * jcolM_, double * M_,
			      int deleteElts, double *additiveDiag_)
{
  neverDeleteElts = (!deleteElts);

  //neverDeleteElts = 1;
  m               = m_;
  n               = n_;
  len             = len_;
  jcolM           = jcolM_;
  krowM           = krowM_;
  M               = M_;

  additiveDiag 	  = additiveDiag_;

  SparseStorage::instances++;
  
}

SparseStorage::~SparseStorage()
{
  if ( !neverDeleteElts ) {
    delete [] jcolM;
    delete [] krowM;
    delete [] M;
  }

  if(additiveDiag) 
  	delete [] additiveDiag;

  SparseStorage::instances--;
}

void SparseStorage::getSize( int& m_, int& n_ )
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
  int extent = scale.length();

  assert( extent == n );
 
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
  int extent = scale.length();

  assert( extent == m );

  int i, j, k;

  for ( i = 0; i < m; i++ ) {
    // Loop over all rows in the sparse matrix
    for( k = krowM[i]; k < krowM[i+1]; k++ ) {
      // Loop over the elements of the sparse row
      j = jcolM[k];
      M[k] = M[k] * scale[i];
    } // End loop over the elements of the sparse row
  } // End loop over all rows in the sparse matrix

}

void SparseStorage::SymmetricScale( OoqpVector& scale_in )
{
  SimpleVector & scale = dynamic_cast<SimpleVector &>(scale_in);
  int extent = scale.length();

  assert( extent == n );
  assert( extent == m );

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
  int i, j, k;

  for ( i = 0; i < m; i++ ) {
    // Loop over all rows in the sparse matrix
    for( k = krowM[i]; k < krowM[i+1]; k++ ) {
      // Loop over the elements of the sparse row
      j = jcolM[k];
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


void SparseStorage::fromAddDense( int row, int col, double * A, int lda,
				  int rowExtent, int colExtent )
{
  int i, j, k, jcurrent;

  assert( row >= 0 && row + rowExtent <= m );
  assert( col >= 0 && col + colExtent <= n );

  for ( i = row; i < row + rowExtent; i++ ) {
    jcurrent = col - 1;
    for( k = krowM[i]; k < krowM[i+1]; k++ ) {
      j = jcolM[k];
      if ( j >= col ) {
	if ( j < col + colExtent ) {
	  A[(i - row) * lda + j - col] += M[k];
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
			      double alpha, double x[], int incx )
{
  int i, j, k;
  double temp;
  for( i = 0; i < m; i++ ) {
    temp = 0;
    for( k = krowM[i]; k < krowM[i+1]; k++ ) {
      j = jcolM[k];
      temp += M[k] * x[j * incx];
#ifdef DEBUG
      assert(j<n);
#endif
    }
    y[i * incy] = beta * y[i * incy] + alpha * temp;
  }
}

void SparseStorage::transMult( double beta,  double y[], int incy,
			       double alpha, double x[], int incx )
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
#ifdef DEBUG
      assert(j<n);
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
#ifdef DEBUG
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
	Y[v +ldy*j] *= beta;  //Y[j +ldy*v] *= beta;
      }
  }
  for( i = 0; i < m; i++ ) {
    for( k = krowM[i]; k < krowM[i+1]; k++ ) {
      j = jcolM[k];
#ifdef DEBUG
      assert(j<n);
#endif
      for (int v = 0; v<ny; v++) { 
	Y[v+ldy*j] += alpha * M[k] * X[i+v*ldx];  //Y[j +ldy*v] *= beta;  
      }
    }
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
#ifdef DEBUG
      assert(j<n);
#endif
      // nonzero element at (i,j)
      // multiply with x[i] and add to y[j]
      y[j - firstrow] += alpha * M[k] * x[i];
    }
  }
}

#ifndef MAX
#define MAX(a,b) ( (a > b) ? a : b)
#endif
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
#ifdef DEBUG
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



void _doubleLexSort_SymOrder( int first[], int n,
		    int second[], double data[], int *goffOrder )
{
  int fi, se, j, k, kinc, inc;
  double dtemp, gTemp;
  const int incs[]  = {1, 5, 19, 41, 109, 209, 505,
		       929, 2161, 3905, 8929, 16001, INT_MAX};
  
  for ( k = 0; incs[k] <= n/2; k++ ) ;

  kinc = k - 1;
  // incs[kinc] is the greatest value in the sequence that is also less
  // than or equal to n/2. 
  // If n == {0,1}, kinc == -1 and so no sort will take place.

  for( ; kinc >= 0; kinc-- ) {
    // Loop over all increments
    inc = incs[kinc];

    for ( k = inc; k < n; k++ ) {
      dtemp = data[k];
	  gTemp = goffOrder[k];
      fi = first[ k ];
      se = second[ k];
      for( j = k; j >= inc; j -= inc ) {
	if ( fi < first[j - inc] ||
	     ( fi == first[j - inc] &&
	       se < second[ j - inc ]) ) {
	  data[j]         = data[j - inc];
	  goffOrder[j]	  = goffOrder[j - inc];	  
	  first[j]        = first[j - inc];
	  second[j]       = second[j - inc];
	} else {
	  break;
	}
      } 
      data[j]    = dtemp;
	  goffOrder[j] = gTemp;
      first[j]   = fi;
      second[j]  = se;
    }
  } // End loop over all increments
}


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
  if( info != 0 ) return;
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



int* SparseStorage::symmetrize_set( int& info)
{
  int i, k, ku;

  int nnz = krowM[m];
  int * irowM = new int[ 2 * nnz ];
  int *goffIDX_temp = new int[2*nnz];

  
  info = 0;

  ku = nnz;
  for ( i = 0; i < m; i++ ) {
    // For all rows
    for( k = krowM[i]; k < krowM[i + 1]; k++ ) {
      // Loop over elements of row i
      irowM[k] = i;
	  goffIDX_temp[k]=k;
      if ( i != jcolM[k] ) {
	// Not a diagonal element
	if ( ku >= len ) {
	  info++;
	} else {
	  // Add the element transpose to the scrambled matrix
	  irowM[ku] = jcolM[k];
	  jcolM[ku] = i;
	  M[ku]     = M[k];
	  goffIDX_temp[ku]=k;
	  ku++;
	}
      } // End not a diagonal element
    } // End loop over elements of column j
  } // End for all columns
  if( info != 0 ) return NULL;
  nnz = ku;

  _doubleLexSort_SymOrder( irowM, nnz, jcolM, M, goffIDX_temp );
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

  int  *goffIDX = new int[nnz];
  memcpy(&goffIDX[0],&goffIDX_temp[0],nnz*sizeof(int));
  

  delete [] goffIDX_temp;
  delete [] irowM;

  return goffIDX;
}


void SparseStorage::symmetrize_valonly( double *val_lower,int *goffIDX)
{
  int k;
  int nnz = krowM[m];
  
  for( k = 0; k < nnz; k++ ) {
    M[k] = val_lower[goffIDX[k]];
  }
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

double SparseStorage::abmaxnorm()
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

// FIXME_NY this is wrong, how to do memmove in an std::map (the index is not in ordered)
void SparseStorage::shiftRows_CorrectMap( int row, int shift, int& info, std::map<int,int> &ValIdxMap )
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
//      memmove( &ValIdxMap[ krowM[row] + shift ],
//	       &ValIdxMap[ krowM[row] ], lcopy * sizeof(int) );  // this is wrong. FIXME_NY
	  ValIdxMap.erase(krowM[row]);
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


void SparseStorage::atPutSpRow_CorrectMap( int row, double A[], int lenA,
				    int jcolA[], int& info, std::map<int,int> &ValIdxMap, int const IDXconstant )
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
    this->shiftRows_CorrectMap( row + 1, count, info, ValIdxMap);
	assert("not suppot yet" && 0);
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
		  ValIdxMap[ik] = ka + IDXconstant;
		}
		break;
      } else if ( jcolM[km] == jcolA[ka] ) {
		// The element in A will replace an element in M
		jcolM[ik] = jcolM[km];
		M[ik]     = A[ka];
		ValIdxMap[ik] = ka + IDXconstant;	
		km--; ka--; ik--;
      } else if ( jcolM[km] > jcolA[ka] ) {
		// This element is in M but not in A
		jcolM[ik] = jcolM[km];
		M[ik]     = M[km];
//		ValIdxMap[ik] = ValIdxMap[km];
		km--; ik--;
      } else {
		// The element is in A, but not in M.
		jcolM[ik] = jcolA[ka];
		assert( jcolM[ik] >= 0 && jcolM[ik] < n );
		M[ik]     = A[ka];
		ValIdxMap[ik] = ka + IDXconstant;
		ka--; ik--;
      }
    }
  }
}


void SparseStorage::copyDiagonalVal_From( int idiag, OoqpVector& vvec, 
						bool firstCall, std::map<int,int> &ValIdxMap)
{
  SimpleVector & v = dynamic_cast<SimpleVector &>(vvec);
  this->copyDiagonalVal_From( idiag, &v[0], 1, v.n, firstCall, ValIdxMap);
}


void SparseStorage::copyDiagonalVal_From( int idiag,
				   double x[], int incx, int diagLength, bool firstCall, std::map<int,int> &ValIdxMap)
{
  int i;
  int info;

  if(firstCall == true){
    for( i = idiag; i < idiag + diagLength; i++ ) { // Loop over elts to be put
      // Search for the diagonal elt.
      int lastk;
      lastk = krowM[i+1];
      for( int k = krowM[i]; k < lastk; k++ ) { // Loop over all elts in row
        if ( i >= jcolM[k] ) { // Found or past the diagonal
		  if( i == jcolM[k] ) { // Found it, overwrite it.
	  	    M[k] = x[incx*(i - idiag)];
		    ValIdxMap[k] = incx*(i - idiag);
//		    DiagIdxMap.insert( pair<int,int>(k,incx*(i - idiag)));
		  } else { // Didn't find it, so insert it
	  	    this->atPutSpRow_CorrectMap( i, &x[incx*(i - idiag)], 1, &i, info, ValIdxMap, incx*(i - idiag));
		  }
		  // Either way, bug out of the loop
		  break;
        } // end if found or past the diagonal
      } // end loop over all elts
    } // end loop over elts to be put
  }else{
    map<int,int>::iterator it;
	for( it=ValIdxMap.begin(); it!=ValIdxMap.end(); it++ ) {
	  M[it->first] = x[it->second];
	}
  }
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

void SparseStorage::transpose(int* krowMt, int* jcolMt, double* Mt)
{
  int ind, pend;//,ptend;
  /////////////////////////////////////////////////////
  // form the transpose 
  /////////////////////////////////////////////////////
  int nnz = krowM[m];

  //cummulative sum to find the number of elements in each row of At, ie column of A.
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


void SparseStorage::transpose_withOriIDX(int* krowMt, int* jcolMt, double* Mt, int* OriIDX)
{
  int ind, pend;//,ptend;
  /////////////////////////////////////////////////////
  // form the transpose 
  /////////////////////////////////////////////////////
  int nnz = krowM[m];

//  if(OriIDX==NULL)
//  	OriIDX = (int*) malloc(nnz*sizeof(int));

  //cummulative sum to find the number of elements in each row of At, ie column of A.
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
	  OriIDX[ ind ] = j;
      w[jcolM[j]] = ind+1; //w[jcolM[j]]++
    }
  }
  //we have the transpose
  delete[] w;
}

void SparseStorage::transpose_withNewIDX(int* krowMt, int* jcolMt, double* Mt, int* NewIDX)
{
  int ind, pend;//,ptend;
  /////////////////////////////////////////////////////
  // form the transpose 
  /////////////////////////////////////////////////////
  int nnz = krowM[m];

//  if(NewIDX==NULL)
//  	NewIDX = (int*) malloc(nnz*sizeof(int));

  //cummulative sum to find the number of elements in each row of At, ie column of A.
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
	  NewIDX[ j ] = ind;
      w[jcolM[j]] = ind+1; //w[jcolM[j]]++
    }
  }
  //we have the transpose
  delete[] w;
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

void SparseStorage::dump(const string& _filename)
{
  ofstream fd(_filename.c_str());
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


void SparseStorage::copyMtxFromDouble(int copyLength,double *values)
{
  if( len < copyLength ) {
    printf("copy from longer array \n");
    exit(0);
  }
  int sLen = (copyLength<len)?copyLength:len;
  
  for( int k = 0; k < sLen; k++ ) {
    M[k]     = values[k];
  }
  
}


void SparseStorage::setAdditiveDiagonal( OoqpVector& v_in )
{
  SimpleVector & Diag = dynamic_cast<SimpleVector &>(v_in);
  int extent = Diag.length();

  assert( extent == n && n == m );
 
  int i;

  for ( i = 0; i < m; i++ ) {
    additiveDiag[i] = Diag[i];
  } 

}


void SparseStorage::fromGetSpRow_WithRowStart( int row, int col,
				      double A[], int lenA, int jcolA[],
				      int& nnz,
				      int colExtent, int& info, int & rowStart )
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
  rowStart = krowM[row];
}



void SparseStorage:: fromGetDense_withMap( int row, int col, double * A,
				       int lda,
				       int rowExtent, int colExtent, int const FirstCall, std::map<int,int> &ValIdxMap )
{
	int i, j, k, jcurrent;

	assert( row >= 0 && row + rowExtent <= m );
	assert( col >= 0 && col + colExtent <= n );

  if(FirstCall==true){
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
					ValIdxMap[(i - row) * lda + j - col] = k;
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
	}// End loop over all rows in range
  }
  else{
	  map<int,int>::iterator it;
	  
	  for ( i = row; i < row + rowExtent; i++ ){
	  	jcurrent = col - 1;
		for( jcurrent++; jcurrent < col+colExtent; jcurrent++ ) {
		  A[(i - row) * lda + jcurrent - col] = 0.0;
		}
	  }
	  
	  for( it=ValIdxMap.begin(); it!=ValIdxMap.end(); it++ ) {
		A[it->first] = M[it->second];
	  }
  }
}

