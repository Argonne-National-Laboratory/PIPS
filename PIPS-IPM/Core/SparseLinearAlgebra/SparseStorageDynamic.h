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

  int m;      // row
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
  void compressStorageValues();

public:
  static int instances;

  int getM() const { return m; };
  int getN() const { return n; }; // todo check for n when inserting rows cols etc
  //int getLen() const { return len; };

  const ROWPTRS* getRowPtr() const { return rowptr; };
  const ROWPTRS getRowPtr(int i) const;

  const int* getJcolM() const { return jcolM; };
  const int getJcolM(int i) const;
  
  const double* getMat() const { return M; };
  const double getMat(int i) const;



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

  virtual void getSize( int& m, int& n ) override;

  virtual void getDiagonal( OoqpVector& vec ) override { assert(0 && "not implemented here"); };
  virtual void setToDiagonal( OoqpVector& vec ) override { assert(0 && "not implemented here"); };

  virtual void atPutDiagonal( int idiag, OoqpVector& x ) override { assert(0 && "not implemented here"); };
  virtual void fromGetDiagonal( int idiag, OoqpVector& x ) override { assert(0 && "not implemented here"); };
  virtual void SymmetricScale ( OoqpVector& vec ) override { assert(0 && "not implemented here"); };
  virtual void ColumnScale ( OoqpVector& vec ) override { assert(0 && "not implemented here"); };
  virtual void RowScale ( OoqpVector& vec ) override { assert(0 && "not implemented here"); };
  virtual void scalarMult( double num ) override { assert(0 && "not implemented here"); };

  void removeEntryAtIndex(int row, int col_idx);
  void removeEntryAtRowCol(int row, int col);

  void clearRow( int row );
  void clearCol( int col );

  void appendRow( const std::vector<double>& row, const std::vector<int>& idx );

  void scaleRow( int row, double factor );

  void addNnzPerRow(double* vec) const;

  void writeToStreamDense( ostream& out) const;
  void writeToStreamDenseRow( stringstream& out, int rowidx) const;

  void restoreOrder();

  bool isTransposedOf( const SparseStorageDynamic& mat_tp) const; //todo!

  SparseStorage* getStaticStorage(double* rowNnz, double* colNnz) const;
  SparseStorageDynamic* getTranspose() const;
};

typedef SmartPointer<SparseStorageDynamic>  SparseStorageDynamicHandle;

#endif /* PIPS_IPM_CORE_SPARSELINEARALGEBRA_SPARSESTORAGEDYNAMIC_H_ */
