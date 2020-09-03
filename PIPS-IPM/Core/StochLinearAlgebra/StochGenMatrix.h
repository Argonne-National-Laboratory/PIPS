#ifndef STOCHGENMATRIX_H
#define STOCHGENMATRIX_H

#include "StochVector_fwd.h"
#include "OoqpVector_fwd.h"
#include "DoubleMatrix.h"
#include "SparseGenMatrix.h"
#include "mpi.h"
#include "pipsport.h"

#include <vector>

class StochGenMatrix : public GenMatrix {
protected:

public:
  /** Constructs a matrix having local A and B blocks having the sizes and number of nz specified by  
      A_m, A_n, A_nnz and B_m, B_n, B_nnz.
      Also sets the global sizes to 'global_m' and 'global_n'. The 'id' parameter is used 
      for output/debug purposes only.
      The matrix that will be created  has no children, just local data.*/
  StochGenMatrix(int id, 
		 long long global_m, long long global_n,
		 int A_m, int A_n, int A_nnz,
		 int B_m, int B_n, int B_nnz,
		 MPI_Comm mpiComm_);

  /** Constructs a matrix with local A, B, and Bl (linking constraints) blocks having the sizes and number of nz specified by
      A_m, A_n, A_nnz, B_m, B_n, B_nnz, and Bl_m, Bl_n, Bl_nnz. Otherwise, identical to the above constructor */
  StochGenMatrix(int id,
		 long long global_m, long long global_n,
		 int A_m, int A_n, int A_nnz,
		 int B_m, int B_n, int B_nnz,
		 int Bl_m, int Bl_n, int Bl_nnz,
		 MPI_Comm mpiComm_);

  /** Constructs a matrix with local A, B, and Bl (linking constraints) blocks set to nullptr */
  StochGenMatrix(int id,
       long long global_m, long long global_n,
       MPI_Comm mpiComm_);

  // constructor for combining scenarios
  virtual ~StochGenMatrix();

  virtual StochGenMatrix* cloneEmptyRows(bool switchToDynamicStorage = false) const;
  virtual StochGenMatrix* cloneFull(bool switchToDynamicStorage = false) const;

  virtual void AddChild(StochGenMatrix* child);

  std::vector<StochGenMatrix*> children;
  SparseGenMatrix* Amat;
  SparseGenMatrix* Bmat;
  SparseGenMatrix* Blmat;

  int id;
  long long m,n;
  MPI_Comm mpiComm;
  int iAmDistrib;
 private:
  OoqpVector* workPrimalVec;
  OoqpVector* getWorkPrimalVec(const StochVector& origin);

  /** trans mult method for children with linking constraints */
  virtual void transMult2 ( double beta,   StochVector& y,
		    double alpha,  StochVector& x,
		    OoqpVector& yvecParent, const OoqpVector& xvecl ) const;

  /** trans mult method for children; does not support linking constraints */
  virtual void transMult2 ( double beta,   StochVector& y,
		    double alpha,  StochVector& x,
		    OoqpVector& yvecParent );

  /** mult method for children; needed only for linking constraints */
  virtual void mult2 ( double beta,  OoqpVector& y,
                        double alpha, OoqpVector& x,
						   OoqpVector& yvecParent );

  /** column scale method for children */
  virtual void ColumnScale2( OoqpVector& vec, OoqpVector& parentvec );

  /** row scale method for children */
  virtual void RowScale2( OoqpVector& vec, OoqpVector* linkingvec );

  virtual void getNnzPerRow(OoqpVectorBase<int>& nnzVec, OoqpVectorBase<int>* linkParent);

  virtual void getNnzPerCol(OoqpVectorBase<int>& nnzVec, OoqpVectorBase<int>* linkParent);

  virtual void addRowSums( OoqpVector& sumVec, OoqpVector* linkParent );
  virtual void addColSums( OoqpVector& sumVec, OoqpVector* linkParent );

  /** internal method needed for handling linking constraints */
  virtual void getRowMinMaxVec( bool getMin, bool initializeVec,
        const OoqpVector* colScaleVec, const OoqpVector* colScaleParent, OoqpVector& minmaxVec, OoqpVector* linkParent);

  virtual void getColMinMaxVec( bool getMin, bool initializeVec,
        const OoqpVector* rowScaleVec, const OoqpVector* rowScaleParent, OoqpVector& minmaxVec, OoqpVector* minmaxParent );


  virtual void initTransposedChild(bool dynamic);
  virtual void initStaticStorageFromDynamic(const OoqpVectorBase<int>& rowNnzVec, const OoqpVectorBase<int>& colNnzVec,
    const OoqpVectorBase<int>* rowLinkVec, const OoqpVectorBase<int>* colParentVec);

  virtual void permuteLinkingVarsChild(const std::vector<unsigned int>& permvec);

  virtual void getLinkVarsNnzChild(std::vector<int>& vec) const;

  virtual void writeToStreamDenseChild(stringstream& out, int offset) const;
  virtual std::string writeToStreamDenseRowLink(int rowidx) const;


 public:
  virtual void updateTransposed();

  void getSize( long long& m, long long& n ) const override;
  void getSize( int& m, int& n ) const override;

  /** The actual number of structural non-zero elements in this sparse
   *  matrix. This includes so-called "accidental" zeros, elements that
   *  are treated as non-zero even though their value happens to be zero.
   */  
  virtual int numberOfNonZeros();

  int isKindOf( int matType ) const override;

  void atPutDense( int row, int col, double * A, int lda,
			   int rowExtent, int colExtent ) override;
  void fromGetDense( int row, int col, double * A, int lda,
			     int rowExtent, int colExtent ) override;
  void ColumnScale( OoqpVector& vec ) override;
  void RowScale( OoqpVector& vec ) override;
  void SymmetricScale( OoqpVector &vec) override;
  void scalarMult( double num) override;
  void fromGetSpRow( int row, int col,
			     double A[], int lenA, int jcolA[], int& nnz,
			     int colExtent, int& info ) override;
  void atPutSubmatrix( int destRow, int destCol, DoubleMatrix& M,
			       int srcRow, int srcCol,
			       int rowExtent, int colExtent ) override;
  void atPutSpRow( int col, double A[], int lenA, int jcolA[],
			   int& info ) override;

  void putSparseTriple( int irow[], int len, int jcol[], double A[], 
				int& info ) override;

  void getDiagonal( OoqpVector& vec ) override;
  void setToDiagonal( OoqpVector& vec ) override;

  /** y = beta * y + alpha * this * x */
  void mult ( double beta,  OoqpVector& y,
                      double alpha, const OoqpVector& x ) const override;

  void transMult ( double beta,   OoqpVector& y,
                           double alpha,  const OoqpVector& x ) const override;

  double abmaxnorm() const override;

  virtual void getLinkVarsNnz(std::vector<int>& vec) const;

  void writeToStream(ostream& out) const override;
  void writeToStreamDense(ostream& out) const override;
  void writeMPSformatRows(ostream& out, int rowType, OoqpVector* irhs) const override;

  /** Make the elements in this matrix symmetric. The elements of interest
   *  must be in the lower triangle, and the upper triangle must be empty.
   *  @param info zero if the operation succeeded. Otherwise, insufficient
   *  space was allocated to symmetrize the matrix.
   */
  virtual void symmetrize( int& info );

  void randomize( double alpha, double beta, double * seed ) override;

  /** initialize (dynamic) transposed matrices for A, B, Bl */
  virtual void initTransposed(bool dynamic = false);
  virtual void deleteTransposed();

  void atPutDiagonal( int idiag, OoqpVector& v ) override;
  void fromGetDiagonal( int idiag, OoqpVector& v ) override;
  void matTransDMultMat(OoqpVector& d, SymMatrix** res) override;
  void matTransDinvMultMat(OoqpVector& d, SymMatrix** res) override;

  void getNnzPerRow(OoqpVectorBase<int>& nnzVec) override
  {
     getNnzPerRow(nnzVec, nullptr);
  };

  void getNnzPerCol(OoqpVectorBase<int>& nnzVec) override
  {
     getNnzPerCol(nnzVec, nullptr);
  };

  /** fill vector with absolute minimum/maximum value of each row */
  void getRowMinMaxVec( bool getMin, bool initializeVec,
        const OoqpVector* colScaleVec, OoqpVector& minmaxVec ) override
  {
     getRowMinMaxVec(getMin, initializeVec, colScaleVec, nullptr, minmaxVec, nullptr);
  };

  /** fill vector with absolute minimum/maximum value of each column */
  void getColMinMaxVec( bool getMin, bool initializeVec,
        const OoqpVector* rowScaleVec, OoqpVector& minmaxVec ) override
  {
     getColMinMaxVec(getMin, initializeVec, rowScaleVec, nullptr, minmaxVec, nullptr);
  };

  virtual void addRowSums( OoqpVector& sumVec ) {addRowSums(sumVec, nullptr);};
  virtual void addColSums( OoqpVector& sumVec ) {addColSums(sumVec, nullptr);};

  virtual void initStaticStorageFromDynamic(const OoqpVectorBase<int>& rowNnzVec, const OoqpVectorBase<int>& colNnzVec)
  {
     initStaticStorageFromDynamic(rowNnzVec, colNnzVec, nullptr, nullptr);
  };
  virtual void freeDynamicStorage();

  /** returns Simple Vector indicating which linking rows have entries in exactly two blocks (indicated by 1.0 versus 0.0)*/
  virtual std::vector<int> get2LinkStartBlocks() const;

  virtual void updateKLinkVarsCount(std::vector<int>& linkCount) const;
  virtual void updateKLinkConsCount(std::vector<int>& linkCount) const;

  virtual void permuteLinkingVars(const std::vector<unsigned int>& permvec);
  virtual void permuteLinkingCons(const std::vector<unsigned int>& permvec);

  virtual bool isRootNodeInSync() const;

  virtual int appendRow( const StochGenMatrix& matrix_row, int child, int row, bool linking );

  /** calculate vec^T * row where row is linking or not in child child and with row index row
   *  for a linking row only the available blocks will be multiplied - currently only possible for dynamic storage! (since 
   *  this was its foremost usecase)
   */
  virtual double localRowTimesVec( const StochVector& vec, int child, int row, bool linking ) const;

  /* y += alpha * RowAt(child, row, linking) */
  virtual void axpyWithRowAt( double alpha, StochVector* y, SimpleVector* y_linking, int child, int row, bool linking) const;
  virtual void axpyWithRowAtPosNeg( double alpha, StochVector* y_pos, SimpleVector* y_link_pos, StochVector* y_neg, SimpleVector* y_link_neg, int child, int row, bool linking ) const;
};


/**
 * Dummy Class 
 */

class StochGenDummyMatrix : public StochGenMatrix {

protected:

public:

  StochGenDummyMatrix(int id)
    : StochGenMatrix(id, 0, 0, 0, 0, 0, 0, 0, 0, MPI_COMM_NULL) {};

  virtual ~StochGenDummyMatrix(){};

  void AddChild(StochGenMatrix* child) override {};

 public:
  void updateTransposed() override {};

  void getSize( int& m, int& n ) const override { m = 0; n = 0; }
  void getSize( long long& m, long long& n ) const override { m = 0; n = 0; }

  StochGenMatrix* cloneEmptyRows(bool switchToDynamicStorage = false) const override { return new StochGenDummyMatrix(id); };
  StochGenMatrix* cloneFull(bool switchToDynamicStorage = false) const  override { return new StochGenDummyMatrix(id); };


  /** The actual number of structural non-zero elements in this sparse
   *  matrix. This includes so-called "accidental" zeros, elements that
   *  are treated as non-zero even though their value happens to be zero.
   */  
  int numberOfNonZeros() override {return 0;} 

  int isKindOf( int matType ) const override ;

  void atPutDense( int row, int col, double * A, int lda,
			   int rowExtent, int colExtent ) override {};
  void fromGetDense( int row, int col, double * A, int lda,
			     int rowExtent, int colExtent ) override {};
  void ColumnScale( OoqpVector& vec ) override {};
  void RowScale( OoqpVector& vec ) override {};
  void SymmetricScale( OoqpVector &vec) override {};
  void scalarMult( double num) override {};
  void fromGetSpRow( int row, int col,
			     double A[], int lenA, int jcolA[], int& nnz,
			     int colExtent, int& info ) override {};

  void atPutSubmatrix( int destRow, int destCol, DoubleMatrix& M,
			       int srcRow, int srcCol,
			       int rowExtent, int colExtent ) override {};

  void atPutSpRow( int col, double A[], int lenA, int jcolA[],
			   int& info ) override {};

  void putSparseTriple( int irow[], int len, int jcol[], double A[], 
				int& info ) override {};

  void getDiagonal( OoqpVector& vec ) override {};
  void setToDiagonal( OoqpVector& vec ) override {};

  /** y = beta * y + alpha * this * x */
  void mult ( double beta,  OoqpVector& y,
                      double alpha, const OoqpVector& x ) const override {};

  /** mult method for children; needed only for linking constraints */
  void mult2 ( double beta,  OoqpVector& y,
                        double alpha, OoqpVector& x,
						   OoqpVector& yvecParent ) override {};

  void transMult ( double beta,   OoqpVector& y,
                           double alpha,  const OoqpVector& x ) const override {};

  void transMult2 ( double beta,   StochVector& y,
		    double alpha,  StochVector& x,
		    OoqpVector& yvecParent, const OoqpVector& xvecl ) const override {};

  void transMult2 ( double beta,   StochVector& y,
		    double alpha,  StochVector& x,
		    OoqpVector& yvecParent ) override {};

  double abmaxnorm() const override { return 0.0; };

  void permuteLinkingVarsChild(const std::vector<unsigned int>& permvec)  override {};
  void getLinkVarsNnzChild(std::vector<int>& vec) const override {};

  void getLinkVarsNnz(std::vector<int>& vec) const override {};
  void writeToStream(ostream& out) const override {};
  void writeToStreamDense(ostream& out) const override {};
  void writeMPSformatRows(ostream& out, int rowType, OoqpVector* irhs) const override {};

 private:
  void writeToStreamDenseChild(stringstream& out, int offset) const override {};
  std::string writeToStreamDenseRowLink(int rowidx) const override { return 0; };

 public:
  /** Make the elements in this matrix symmetric. The elements of interest
   *  must be in the lower triangle, and the upper triangle must be empty.
   *  @param info zero if the operation succeeded. Otherwise, insufficient
   *  space was allocated to symmetrize the matrix.
   */
  void symmetrize( int& info ) override {};

  void randomize( double alpha, double beta, double * seed ) override {};

  void atPutDiagonal( int idiag, OoqpVector& v ) override {};
  void fromGetDiagonal( int idiag, OoqpVector& v ) override {};

  void initTransposedChild(bool dynamic) override {};

  void ColumnScale2( OoqpVector& vec, OoqpVector& parentvec )override {};
  void RowScale2( OoqpVector& vec, OoqpVector* linkingvec )override {};

  void initTransposed(bool dynamic = false) override {};
  void deleteTransposed() override {};

  void getNnzPerRow(OoqpVectorBase<int>& nnzVec, OoqpVectorBase<int>* linkParent) override {};

  void getNnzPerCol(OoqpVectorBase<int>& nnzVec, OoqpVectorBase<int>* linkParent) override {};

  void getNnzPerRow(OoqpVectorBase<int>& nnzVec) override {};

  void getNnzPerCol(OoqpVectorBase<int>& nnzVec) override {};

  void getRowMinMaxVec( bool getMin, bool initializeVec,
        const OoqpVector* colScaleVec, const OoqpVector* colScaleParent, OoqpVector& minmaxVec, OoqpVector* linkParent )override {};

  void getColMinMaxVec( bool getMin, bool initializeVec,
        const OoqpVector* rowScaleVec, const OoqpVector* rowScaleParent, OoqpVector& minmaxVec, OoqpVector* minmaxParent )override {};

  void getRowMinMaxVec( bool getMin, bool initializeVec,
        const OoqpVector* colScaleVec, OoqpVector& minmaxVec )override {};

  void getColMinMaxVec( bool getMin, bool initializeVec,
        const OoqpVector* rowScaleVec, OoqpVector& minmaxVec )override {};

  void addRowSums( OoqpVector& sumVec, OoqpVector* linkParent ) override {};
  void addColSums( OoqpVector& sumVec, OoqpVector* linkParent ) override {};
  void addRowSums( OoqpVector& vec ) override {};
  void addColSums( OoqpVector& vec ) override {};

  void freeDynamicStorage() override {};
  void initStaticStorageFromDynamic(const OoqpVectorBase<int>& rowNnzVec, const OoqpVectorBase<int>& colNnzVec) {};
  void initStaticStorageFromDynamic(const OoqpVectorBase<int>& rowNnzVec, const OoqpVectorBase<int>& colNnzVec, 
    const OoqpVectorBase<int>* rowLinkVec, const OoqpVectorBase<int>* colParentVec) {};

  std::vector<int> get2LinkStartBlocks() const override { return std::vector<int>(); };

  void updateKLinkVarsCount(std::vector<int>& linkCount) const override {};
  void updateKLinkConsCount(std::vector<int>& linkCount) const override {};

  void permuteLinkingVars(const std::vector<unsigned int>& permvec) override {};
  void permuteLinkingCons(const std::vector<unsigned int>& permvec) override {};

  bool isRootNodeInSync() const override { return true; };

  int appendRow( const StochGenMatrix& matrix_row, int child, int row, bool linking ) override { assert("CANNOT APPEND ROW TO DUMMY MATRIX"); return -1; };
  double localRowTimesVec( const StochVector& vec, int child, int row, bool linking ) const override { assert("CANNOT MULTIPLY ROW WITH DUMMY MATRIX"); return -1; };

  void axpyWithRowAt( double alpha, StochVector* y, SimpleVector* y_linking, int child, int row, bool linking) const override {};
  void axpyWithRowAtPosNeg( double alpha, StochVector* y_pos, SimpleVector* y_link_pos, StochVector* y_neg, SimpleVector* y_link_neg, int child, int row, bool linking ) const override {};
};


typedef SmartPointer<StochGenMatrix> StochGenMatrixHandle;

#endif
