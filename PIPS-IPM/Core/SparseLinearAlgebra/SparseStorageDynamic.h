/*
 * SparseStorageDynamic.h
 *
 *  Created on: 07.02.2018
 *      Author: bzfrehfe
 */

#ifndef PIPS_IPM_CORE_SPARSELINEARALGEBRA_SPARSESTORAGEDYNAMIC_H_
#define PIPS_IPM_CORE_SPARSELINEARALGEBRA_SPARSESTORAGEDYNAMIC_H_

#include "pipsport.h"
#include "DoubleMatrix.h"
#include "SparseStorage.h"
#include <vector>

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

  int m;      // rows
  int m_len;  // length row array
  int n;      // cols
  int len;    // length col/value array
  int len_free;

  ROWPTRS * rowptr;
  int * jcolM;
  double * M;

  /* doubles the size of rowptr */
  void extendStorageRows();

  /* compresses storage and doubles size of the col entry storage */
  void extendStorageValues();

  /* shifts all rows such that every row is again row + spareRatio length */ 
  void rebuildSpareStructure(int guaranteed_spare = 0);

public:
  static int instances;

  int getM() const { return m; };
  int getN() const { return n; };
  int getNVals() const { return len-len_free; };

  const ROWPTRS* getRowPtr() const { return rowptr; };
  const ROWPTRS getRowPtr(int i) const;

  const int* getJcolM() const { return jcolM; };
  const int getJcolM(int i) const;
  
  const double* getMat() const { return M; };
  const double& getMat(int i) const;
  void setMat(int i, double val);

  SparseStorageDynamic( const SparseStorage& storage, double spareRatio = 0.2 );
  SparseStorageDynamic( int m, int n, int len, double spareRatio = 0.2 );
  SparseStorageDynamic( const SparseStorageDynamic &dynamicStorage);

  ~SparseStorageDynamic();

  virtual void atPutDense( int row, int col, double * A, int lda,
            int rowExtent, int colExtent ) override { assert(0 && "not implemented here"); };
  virtual void fromGetDense( int row, int col, double * A, int lda,
              int rowExtent, int colExtent ) override { assert(0 && "not implemented here"); };
  virtual void atPutSpRow( int row, double A[], int lenA, int jcolA[],
                           int& info ) override { assert(0 && "not implemented here"); };

  virtual void fromGetSpRow( int row, int col,
                             double A[], int lenA, int jcolA[], int& nnz,
                             int colExtent, int& info ) override { assert(0 && "not implemented here"); };

  void getSize( int& m, int& n ) const override;

  void getDiagonal( OoqpVector& vec ) override { assert(0 && "not implemented here"); };
  void setToDiagonal( OoqpVector& vec ) override { assert(0 && "not implemented here"); };

  void atPutDiagonal( int idiag, OoqpVector& x ) override { assert(0 && "not implemented here"); };
  void fromGetDiagonal( int idiag, OoqpVector& x ) override { assert(0 && "not implemented here"); };
  void SymmetricScale ( OoqpVector& vec ) override { assert(0 && "not implemented here"); };
  void ColumnScale ( OoqpVector& vec ) override { assert(0 && "not implemented here"); };
  void RowScale ( OoqpVector& vec ) override { assert(0 && "not implemented here"); };
  void scalarMult( double num ) override { assert(0 && "not implemented here"); };

  void removeEntryAtIndex(int row, int col_idx);
  void removeEntryAtRowCol(int row, int col);

  bool addColToRow( double coeff, int col, int row );

  void clearRow( int row );
  void clearCol( int col );

  void appendRow( const SparseStorageDynamic& storage, int row );

  double rowTimesVec( const double* vec, int length, int row) const;
  void axpyWithRowAt( double alpha, double* y, int length, int row) const;
  void axpyWithRowAtPosNeg( double alpha, double * y_pos, double* y_neg, int length, int row) const;

  void scaleRow( int row, double factor );

  void addNnzPerRow(int* vec) const;

  void writeToStreamDense( std::ostream& out) const;
  void writeToStreamDenseRow( std::stringstream& out, int rowidx) const;

  void restoreOrder();

  double abmaxnorm() const override;

  bool isTransposedOf( const SparseStorageDynamic& mat_tp) const; // TODO..

  SparseStorage* getStaticStorage(const int* rowNnz, const int* colNnz) const;
  SparseStorageDynamic* getTranspose() const;

  void getRowMaxVec(const double* colScaleVec, double* vec) const;
  void getRowMinMaxVec(bool getMin, const double* colScaleVec, double* vec) const;
  void getRowMinVec(const double* colScaleVec, double* vec) const;

};

typedef SmartPointer<SparseStorageDynamic>  SparseStorageDynamicHandle;

#endif /* PIPS_IPM_CORE_SPARSELINEARALGEBRA_SPARSESTORAGEDYNAMIC_H_ */
