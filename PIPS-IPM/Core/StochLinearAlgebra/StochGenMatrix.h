#ifndef STOCHGENMATRIX_H
#define STOCHGENMATRIX_H

#include "DoubleMatrix.h"
#include "SparseGenMatrix.h"
#include <vector>
#include "mpi.h"
#include "pipsport.h"

class OoqpVector;
class StochVector;

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

  /** Constructs a matrix with local A, B, and Bl (linking constraints) blocks set to NULL */
  StochGenMatrix(int id,
       long long global_m, long long global_n,
       MPI_Comm mpiComm_);

  // constructor for combining scenarios
  virtual ~StochGenMatrix();

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
		    OoqpVector& yvecParent, OoqpVector& xvecl );

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

  virtual void getNnzPerRow(OoqpVector& nnzVec, OoqpVector* linkParent);

  virtual void getNnzPerCol(OoqpVector& nnzVec, OoqpVector* linkParent);

  virtual void addRowSums( OoqpVector& sumVec, OoqpVector* linkParent );
  virtual void addColSums( OoqpVector& sumVec, OoqpVector* linkParent );

  /** internal method needed for handling linking constraints */
  virtual void getRowMinMaxVec( bool getMin, bool initializeVec,
        const OoqpVector* colScaleVec, const OoqpVector* colScaleParent, OoqpVector& minmaxVec, OoqpVector* linkParent);

  virtual void getColMinMaxVec( bool getMin, bool initializeVec,
        const OoqpVector* rowScaleVec, const OoqpVector* rowScaleParent, OoqpVector& minmaxVec, OoqpVector* minmaxParent );


  virtual void initTransposedChild(bool dynamic);
  virtual void initStaticStorageFromDynamic(const OoqpVector& rowNnzVec, const OoqpVector& colNnzVec, const OoqpVector* rowLinkVec, const OoqpVector* colParentVec);

  virtual void permuteLinkingVarsChild(const std::vector<unsigned int>& permvec);

  virtual void getLinkVarsNnzChild(std::vector<int>& vec) const;

  virtual void writeToStreamDenseChild(stringstream& out, int offset) const;
  virtual std::string writeToStreamDenseRowLink(int rowidx) const;


 public:
  virtual void updateTransposed();

  virtual void getSize( long long& m, long long& n ) override;
  virtual void getSize( int& m, int& n ) override;

  /** The actual number of structural non-zero elements in this sparse
   *  matrix. This includes so-called "accidental" zeros, elements that
   *  are treated as non-zero even though their value happens to be zero.
   */  
  virtual int numberOfNonZeros();

  virtual int isKindOf( int matType ) const override;

  virtual void atPutDense( int row, int col, double * A, int lda,
			   int rowExtent, int colExtent ) override;
  virtual void fromGetDense( int row, int col, double * A, int lda,
			     int rowExtent, int colExtent ) override;
  virtual void ColumnScale( OoqpVector& vec ) override;
  virtual void RowScale( OoqpVector& vec ) override;
  virtual void SymmetricScale( OoqpVector &vec) override;
  virtual void scalarMult( double num) override;
  virtual void fromGetSpRow( int row, int col,
			     double A[], int lenA, int jcolA[], int& nnz,
			     int colExtent, int& info ) override;
  virtual void atPutSubmatrix( int destRow, int destCol, DoubleMatrix& M,
			       int srcRow, int srcCol,
			       int rowExtent, int colExtent ) override;
  virtual void atPutSpRow( int col, double A[], int lenA, int jcolA[],
			   int& info ) override;

  virtual void putSparseTriple( int irow[], int len, int jcol[], double A[], 
				int& info ) override;

  virtual void getDiagonal( OoqpVector& vec ) override;
  virtual void setToDiagonal( OoqpVector& vec ) override;

  /** y = beta * y + alpha * this * x */
  virtual void mult ( double beta,  OoqpVector& y,
                      double alpha, OoqpVector& x ) override;

  virtual void transMult ( double beta,   OoqpVector& y,
                           double alpha,  OoqpVector& x ) override;

  virtual double abmaxnorm() override;

  virtual void getLinkVarsNnz(std::vector<int>& vec) const;

  virtual void writeToStream(ostream& out) const override;
  virtual void writeToStreamDense(ostream& out) const override;
  virtual void writeMPSformatRows(ostream& out, int rowType, OoqpVector* irhs) const override;

  /** Make the elements in this matrix symmetric. The elements of interest
   *  must be in the lower triangle, and the upper triangle must be empty.
   *  @param info zero if the operation succeeded. Otherwise, insufficient
   *  space was allocated to symmetrize the matrix.
   */
  virtual void symmetrize( int& info );

  virtual void randomize( double alpha, double beta, double * seed ) override;

  /** initialize (dynamic) transposed matrices for A, B, Bl */
  virtual void initTransposed(bool dynamic = false);
  virtual void deleteTransposed();

  virtual void atPutDiagonal( int idiag, OoqpVector& v ) override;
  virtual void fromGetDiagonal( int idiag, OoqpVector& v ) override;
  void matTransDMultMat(OoqpVector& d, SymMatrix** res) override;
  void matTransDinvMultMat(OoqpVector& d, SymMatrix** res) override;

  virtual void getNnzPerRow(OoqpVector& nnzVec) override
  {
     getNnzPerRow(nnzVec, NULL);
  };

  virtual void getNnzPerCol(OoqpVector& nnzVec) override
  {
     getNnzPerCol(nnzVec, NULL);
  };

  /** fill vector with absolute minimum/maximum value of each row */
  virtual void getRowMinMaxVec( bool getMin, bool initializeVec,
        const OoqpVector* colScaleVec, OoqpVector& minmaxVec ) override
  {
     getRowMinMaxVec(getMin, initializeVec, colScaleVec, NULL, minmaxVec, NULL);
  };

  /** fill vector with absolute minimum/maximum value of each column */
  virtual void getColMinMaxVec( bool getMin, bool initializeVec,
        const OoqpVector* rowScaleVec, OoqpVector& minmaxVec ) override
  {
     getColMinMaxVec(getMin, initializeVec, rowScaleVec, NULL, minmaxVec, NULL);
  };

  virtual void addRowSums( OoqpVector& sumVec ) override { addRowSums(sumVec, NULL); };
  virtual void addColSums( OoqpVector& sumVec ) override { addColSums(sumVec, NULL); };

  virtual void initStaticStorageFromDynamic(const OoqpVector& rowNnzVec, const OoqpVector& colNnzVec)
  {
     initStaticStorageFromDynamic(rowNnzVec, colNnzVec, NULL, NULL);
  };
  virtual void freeDynamicStorage();

  /** returns Simple Vector indicating which linking rows have entries in exactly two blocks (indicated by 1.0 versus 0.0)*/
  virtual std::vector<int> get2LinkStartBlocks() const;

  virtual void updateKLinkVarsCount(std::vector<int>& linkCount) const;
  virtual void updateKLinkConsCount(std::vector<int>& linkCount) const;

  virtual void permuteLinkingVars(const std::vector<unsigned int>& permvec);
  virtual void permuteLinkingCons(const std::vector<unsigned int>& permvec);

  virtual bool isRootNodeInSync() const;

  virtual int appendRow( const OoqpVector& row, int child, bool linking );
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

  virtual void AddChild(StochGenMatrix* child) override {};

 public:
  virtual void updateTransposed() override {};

  virtual void getSize( int& m, int& n ) override {m=0; n=0;}
  virtual void getSize( long long& m, long long& n ) override {m=0; n=0;}

  virtual StochGenMatrix* cloneFull(bool switchToDynamicStorage = false) const  override { return new StochGenDummyMatrix(id); };


  /** The actual number of structural non-zero elements in this sparse
   *  matrix. This includes so-called "accidental" zeros, elements that
   *  are treated as non-zero even though their value happens to be zero.
   */  
  virtual int numberOfNonZeros() override {return 0;} 

  virtual int isKindOf( int matType ) const override ;

  virtual void atPutDense( int row, int col, double * A, int lda,
			   int rowExtent, int colExtent ) override {};
  virtual void fromGetDense( int row, int col, double * A, int lda,
			     int rowExtent, int colExtent ) override {};
  virtual void ColumnScale( OoqpVector& vec ) override {};
  virtual void RowScale( OoqpVector& vec ) override {};
  virtual void SymmetricScale( OoqpVector &vec) override {};
  virtual void scalarMult( double num) override {};
  virtual void fromGetSpRow( int row, int col,
			     double A[], int lenA, int jcolA[], int& nnz,
			     int colExtent, int& info ) override {};

  virtual void atPutSubmatrix( int destRow, int destCol, DoubleMatrix& M,
			       int srcRow, int srcCol,
			       int rowExtent, int colExtent ) override {};

  virtual void atPutSpRow( int col, double A[], int lenA, int jcolA[],
			   int& info ) override {};

  virtual void putSparseTriple( int irow[], int len, int jcol[], double A[], 
				int& info ) override {};

  virtual void getDiagonal( OoqpVector& vec ) override {};
  virtual void setToDiagonal( OoqpVector& vec ) override {};

  /** y = beta * y + alpha * this * x */
  virtual void mult ( double beta,  OoqpVector& y,
                      double alpha, OoqpVector& x ) override {};

  /** mult method for children; needed only for linking constraints */
  virtual void mult2 ( double beta,  OoqpVector& y,
                        double alpha, OoqpVector& x,
						   OoqpVector& yvecParent ) override {};

  virtual void transMult ( double beta,   OoqpVector& y,
                           double alpha,  OoqpVector& x ) override {};

  virtual void transMult2 ( double beta,   StochVector& y,
		    double alpha,  StochVector& x,
		    OoqpVector& yvecParent, OoqpVector& xvecl ) override {};

  virtual void transMult2 ( double beta,   StochVector& y,
		    double alpha,  StochVector& x,
		    OoqpVector& yvecParent ) override {};

  virtual double abmaxnorm() override { return 0.0; };

  virtual void permuteLinkingVarsChild(const std::vector<unsigned int>& permvec)  override {};
  virtual void getLinkVarsNnzChild(std::vector<int>& vec) const override {};

  virtual void getLinkVarsNnz(std::vector<int>& vec) const override {};
  virtual void writeToStream(ostream& out) const override {};
  virtual void writeToStreamDense(ostream& out) const override {};
  virtual void writeMPSformatRows(ostream& out, int rowType, OoqpVector* irhs) const override {};
  virtual void writeToStreamDenseChild(stringstream& out, int offset) const override {};
  virtual std::string writeToStreamDenseRowLink(int rowidx) const override { return 0; };

  /** Make the elements in this matrix symmetric. The elements of interest
   *  must be in the lower triangle, and the upper triangle must be empty.
   *  @param info zero if the operation succeeded. Otherwise, insufficient
   *  space was allocated to symmetrize the matrix.
   */
  virtual void symmetrize( int& info ) override {};

  virtual void randomize( double alpha, double beta, double * seed ) override {};

  virtual void atPutDiagonal( int idiag, OoqpVector& v ) override {};
  virtual void fromGetDiagonal( int idiag, OoqpVector& v ) override {};

  virtual void initTransposedChild(bool dynamic) override {};

  virtual void ColumnScale2( OoqpVector& vec, OoqpVector& parentvec )override {};
  virtual void RowScale2( OoqpVector& vec, OoqpVector* linkingvec )override {};

  virtual void initTransposed(bool dynamic = false) override {};
  virtual void deleteTransposed() override {};

  virtual void getNnzPerRow(OoqpVector& nnzVec, OoqpVector* linkParent) override {};

  virtual void getNnzPerCol(OoqpVector& nnzVec, OoqpVector* linkParent) override {};

  virtual void getNnzPerRow(OoqpVector& nnzVec) override {};

  virtual void getNnzPerCol(OoqpVector& nnzVec) override {};

  virtual void getRowMinMaxVec( bool getMin, bool initializeVec,
        const OoqpVector* colScaleVec, const OoqpVector* colScaleParent, OoqpVector& minmaxVec, OoqpVector* linkParent )override {};

  virtual void getColMinMaxVec( bool getMin, bool initializeVec,
        const OoqpVector* rowScaleVec, const OoqpVector* rowScaleParent, OoqpVector& minmaxVec, OoqpVector* minmaxParent )override {};

  virtual void getRowMinMaxVec( bool getMin, bool initializeVec,
        const OoqpVector* colScaleVec, OoqpVector& minmaxVec )override {};

  virtual void getColMinMaxVec( bool getMin, bool initializeVec,
        const OoqpVector* rowScaleVec, OoqpVector& minmaxVec )override {};

  virtual void addRowSums( OoqpVector& sumVec, OoqpVector* linkParent ) override {};
  virtual void addColSums( OoqpVector& sumVec, OoqpVector* linkParent ) override {};
  virtual void addRowSums( OoqpVector& vec ) override {};
  virtual void addColSums( OoqpVector& vec ) override {};

  virtual void initStaticStorageFromDynamic(const OoqpVector& rowNnzVec, const OoqpVector& colNnzVec) override {};
  virtual void initStaticStorageFromDynamic(const OoqpVector& rowNnzVec, const OoqpVector& colNnzVec, const OoqpVector* rowLinkVec, const OoqpVector* colParentVec) override {};

  virtual void freeDynamicStorage() override {};

  virtual std::vector<int> get2LinkStartBlocks() const override { return std::vector<int>(); };

  virtual void updateKLinkVarsCount(std::vector<int>& linkCount) const  override{};
  virtual void updateKLinkConsCount(std::vector<int>& linkCount) const override {};

  virtual void permuteLinkingVars(const std::vector<unsigned int>& permvec) override {};
  virtual void permuteLinkingCons(const std::vector<unsigned int>& permvec) override {};

  virtual bool isRootNodeInSync() const override { return true; };

  virtual int appendRow( const OoqpVector& row, int child, bool linking ) override { assert("CANNOT APEND ROW TO DUMMY MATRIX"); return -1; };
};


typedef SmartPointer<StochGenMatrix> StochGenMatrixHandle;

#endif
