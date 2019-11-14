/*
 * SparseStorageDynamic.h
 *
 *  Created on: 07.02.2018
 *      Author: bzfrehfe
 */

#ifndef PIPS_IPM_CORE_SPARSELINEARALGEBRA_SPARSESTORAGEDYNAMIC_H_
#define PIPS_IPM_CORE_SPARSELINEARALGEBRA_SPARSESTORAGEDYNAMIC_H_

#include "DoubleMatrix.h"
#include "SparseStorage.h"

typedef struct
{
   int start;
   int end;
} ROWPTRS;

struct first_is_smaller
{
    bool operator()(const std::pair<int, double>& x, const std::pair<int, double>& y) const
    {
        return x.first < y.first;
    }
};

/** A class for managing the matrix elements used by sparse matrices.
 *  @ingroup SparseLinearAlgebra
 */
class SparseStorageDynamic : public DoubleStorage {

private:
  const double spareRatio;

public:
  static int instances;

  int m;
  int n;
  int len;

  ROWPTRS * rowptr;
  int * jcolM;
  double * M;

  SparseStorageDynamic( const SparseStorage& storage, double spareRatio = 0.2 );
  SparseStorageDynamic( int m, int n, int len, double spareRatio = 0.2 );
  SparseStorageDynamic( const SparseStorageDynamic &dynamicStorage);

  ~SparseStorageDynamic();

  virtual void atPutDense( int row, int col, double * A, int lda,
            int rowExtent, int colExtent ) { assert(0 && "not implemented here"); };
  virtual void fromGetDense( int row, int col, double * A, int lda,
              int rowExtent, int colExtent ) { assert(0 && "not implemented here"); };
  virtual void atPutSpRow( int row, double A[], int lenA, int jcolA[],
                           int& info ) { assert(0 && "not implemented here"); };

  virtual void fromGetSpRow( int row, int col,
                             double A[], int lenA, int jcolA[], int& nnz,
                             int colExtent, int& info ) { assert(0 && "not implemented here"); };

  virtual void getSize( int& m, int& n );

  virtual void getDiagonal( OoqpVector& vec ) { assert(0 && "not implemented here"); };
  virtual void setToDiagonal( OoqpVector& vec ) { assert(0 && "not implemented here"); };

  virtual void atPutDiagonal( int idiag, OoqpVector& x ) { assert(0 && "not implemented here"); };
  virtual void fromGetDiagonal( int idiag, OoqpVector& x ) { assert(0 && "not implemented here"); };
  virtual void SymmetricScale ( OoqpVector& vec ) { assert(0 && "not implemented here"); };
  virtual void ColumnScale ( OoqpVector& vec ) { assert(0 && "not implemented here"); };
  virtual void RowScale ( OoqpVector& vec ) { assert(0 && "not implemented here"); };
  virtual void scalarMult( double num ) { assert(0 && "not implemented here"); };

  void addNnzPerRow(int* vec) const;

  void writeToStreamDense( std::ostream& out) const;
  void writeToStreamDenseRow( std::stringstream& out, int rowidx) const;

  void restoreOrder();

  SparseStorage* getStaticStorage(const int* rowNnz, const int* colNnz) const;
  SparseStorageDynamic* getTranspose() const;
};

typedef SmartPointer<SparseStorageDynamic>  SparseStorageDynamicHandle;

#endif /* PIPS_IPM_CORE_SPARSELINEARALGEBRA_SPARSESTORAGEDYNAMIC_H_ */
