#ifndef STOCHGENMATRIX_H
#define STOCHGENMATRIX_H

#include "StochVector_fwd.h"
#include "OoqpVector_fwd.h"
#include "DoubleMatrix.h"
#include "SparseGenMatrix.h"
#include "mpi.h"

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

  /** Constructs a matrix with local A, B, and Bl (linking constraints) blocks set to NULL */
  StochGenMatrix(int id,
       long long global_m, long long global_n,
       MPI_Comm mpiComm_);

  // constructor for combining scenarios
  //StochGenMatrix(const vector<StochGenMatrix*> &blocks); -- not needed; cpetra
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
  virtual void initStaticStorageFromDynamic(const OoqpVectorBase<int>& rowNnzVec, const OoqpVectorBase<int>& colNnzVec, const OoqpVectorBase<int>* rowLinkVec, const OoqpVectorBase<int>* colParentVec);

  virtual void permuteLinkingVarsChild(const std::vector<unsigned int>& permvec);

  virtual void getLinkVarsNnzChild(std::vector<int>& vec) const;

  virtual void writeToStreamDenseChild(stringstream& out, int offset) const;
  virtual std::string writeToStreamDenseRowLink(int rowidx) const;


 public:
  virtual void updateTransposed();

  virtual void getSize( long long& m, long long& n );
  virtual void getSize( int& m, int& n );

  /** The actual number of structural non-zero elements in this sparse
   *  matrix. This includes so-called "accidental" zeros, elements that
   *  are treated as non-zero even though their value happens to be zero.
   */  
  virtual int numberOfNonZeros();

  virtual int isKindOf( int matType ) const;

  virtual void atPutDense( int row, int col, double * A, int lda,
			   int rowExtent, int colExtent );
  virtual void fromGetDense( int row, int col, double * A, int lda,
			     int rowExtent, int colExtent );
  virtual void ColumnScale( OoqpVector& vec );
  virtual void RowScale( OoqpVector& vec );
  virtual void SymmetricScale( OoqpVector &vec);
  virtual void scalarMult( double num);
  virtual void fromGetSpRow( int row, int col,
			     double A[], int lenA, int jcolA[], int& nnz,
			     int colExtent, int& info );
  virtual void atPutSubmatrix( int destRow, int destCol, DoubleMatrix& M,
			       int srcRow, int srcCol,
			       int rowExtent, int colExtent );
  virtual void atPutSpRow( int col, double A[], int lenA, int jcolA[],
			   int& info );

  virtual void putSparseTriple( int irow[], int len, int jcol[], double A[], 
				int& info );

  virtual void getDiagonal( OoqpVector& vec );
  virtual void setToDiagonal( OoqpVector& vec );

  /** y = beta * y + alpha * this * x */
  virtual void mult ( double beta,  OoqpVector& y,
                      double alpha, OoqpVector& x );

  virtual void transMult ( double beta,   OoqpVector& y,
                           double alpha,  OoqpVector& x );

  virtual double abmaxnorm();

  virtual void getLinkVarsNnz(std::vector<int>& vec) const;

  virtual void writeToStream(ostream& out) const;
  virtual void writeToStreamDense(ostream& out) const;
  virtual void writeMPSformatRows(ostream& out, int rowType, OoqpVector* irhs) const;

  /** Make the elements in this matrix symmetric. The elements of interest
   *  must be in the lower triangle, and the upper triangle must be empty.
   *  @param info zero if the operation succeeded. Otherwise, insufficient
   *  space was allocated to symmetrize the matrix.
   */
  virtual void symmetrize( int& info );

  virtual void randomize( double alpha, double beta, double * seed );

  /** initialize (dynamic) transposed matrices for A, B, Bl */
  virtual void initTransposed(bool dynamic = false);
  virtual void deleteTransposed();

  virtual void atPutDiagonal( int idiag, OoqpVector& v );
  virtual void fromGetDiagonal( int idiag, OoqpVector& v );
  void matTransDMultMat(OoqpVector& d, SymMatrix** res);
  void matTransDinvMultMat(OoqpVector& d, SymMatrix** res);

  virtual void getNnzPerRow(OoqpVectorBase<int>& nnzVec)
  {
     getNnzPerRow(nnzVec, NULL);
  };

  virtual void getNnzPerCol(OoqpVectorBase<int>& nnzVec)
  {
     getNnzPerCol(nnzVec, NULL);
  };

  /** fill vector with absolute minimum/maximum value of each row */
  virtual void getRowMinMaxVec( bool getMin, bool initializeVec,
        const OoqpVector* colScaleVec, OoqpVector& minmaxVec )
  {
     getRowMinMaxVec(getMin, initializeVec, colScaleVec, NULL, minmaxVec, NULL);
  };

  /** fill vector with absolute minimum/maximum value of each column */
  virtual void getColMinMaxVec( bool getMin, bool initializeVec,
        const OoqpVector* rowScaleVec, OoqpVector& minmaxVec )
  {
     getColMinMaxVec(getMin, initializeVec, rowScaleVec, NULL, minmaxVec, NULL);
  };

  virtual void addRowSums( OoqpVector& sumVec ) {addRowSums(sumVec, NULL);};
  virtual void addColSums( OoqpVector& sumVec ) {addColSums(sumVec, NULL);};

  virtual void initStaticStorageFromDynamic(const OoqpVectorBase<int>& rowNnzVec, const OoqpVectorBase<int>& colNnzVec)
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

  virtual void AddChild(StochGenMatrix* child){};

 public:
  virtual void updateTransposed() {};

  virtual void getSize( int& m, int& n ){m=0; n=0;}
  virtual void getSize( long long& m, long long& n ){m=0; n=0;}

  virtual StochGenMatrix* cloneFull(bool switchToDynamicStorage = false) const { return new StochGenDummyMatrix(id); };


  /** The actual number of structural non-zero elements in this sparse
   *  matrix. This includes so-called "accidental" zeros, elements that
   *  are treated as non-zero even though their value happens to be zero.
   */  
  virtual int numberOfNonZeros(){return 0;}

  virtual int isKindOf( int matType ) const;

  virtual void atPutDense( int row, int col, double * A, int lda,
			   int rowExtent, int colExtent ){};
  virtual void fromGetDense( int row, int col, double * A, int lda,
			     int rowExtent, int colExtent ){};
  virtual void ColumnScale( OoqpVector& vec ){};
  virtual void RowScale( OoqpVector& vec ){};
  virtual void SymmetricScale( OoqpVector &vec){};
  virtual void scalarMult( double num){};
  virtual void fromGetSpRow( int row, int col,
			     double A[], int lenA, int jcolA[], int& nnz,
			     int colExtent, int& info ){};

  virtual void atPutSubmatrix( int destRow, int destCol, DoubleMatrix& M,
			       int srcRow, int srcCol,
			       int rowExtent, int colExtent ){};

  virtual void atPutSpRow( int col, double A[], int lenA, int jcolA[],
			   int& info ){};

  virtual void putSparseTriple( int irow[], int len, int jcol[], double A[], 
				int& info ){};

  virtual void getDiagonal( OoqpVector& vec ){};
  virtual void setToDiagonal( OoqpVector& vec ){};

  /** y = beta * y + alpha * this * x */
  virtual void mult ( double beta,  OoqpVector& y,
                      double alpha, OoqpVector& x ){};

  /** mult method for children; needed only for linking constraints */
  virtual void mult2 ( double beta,  OoqpVector& y,
                        double alpha, OoqpVector& x,
						   OoqpVector& yvecParent ){};

  virtual void transMult ( double beta,   OoqpVector& y,
                           double alpha,  OoqpVector& x ){};

  virtual void transMult2 ( double beta,   StochVector& y,
		    double alpha,  StochVector& x,
		    OoqpVector& yvecParent, OoqpVector& xvecl ){};

  virtual void transMult2 ( double beta,   StochVector& y,
		    double alpha,  StochVector& x,
		    OoqpVector& yvecParent ){};

  virtual double abmaxnorm(){ return 0.0; };

  virtual void permuteLinkingVarsChild(const std::vector<unsigned int>& permvec) {};
  virtual void getLinkVarsNnzChild(std::vector<int>& vec) const{};

  virtual void getLinkVarsNnz(std::vector<int>& vec) const{};
  virtual void writeToStream(ostream& out) const{};
  virtual void writeToStreamDense(ostream& out) const{};
  virtual void writeMPSformatRows(ostream& out, int rowType, OoqpVector* irhs) const{};
  virtual void writeToStreamDenseChild(stringstream& out, int offset) const{};
  virtual std::string writeToStreamDenseRowLink(int rowidx) const{return 0;};

  /** Make the elements in this matrix symmetric. The elements of interest
   *  must be in the lower triangle, and the upper triangle must be empty.
   *  @param info zero if the operation succeeded. Otherwise, insufficient
   *  space was allocated to symmetrize the matrix.
   */
  virtual void symmetrize( int& info ){};

  virtual void randomize( double alpha, double beta, double * seed ){};

  virtual void atPutDiagonal( int idiag, OoqpVector& v ){};
  virtual void fromGetDiagonal( int idiag, OoqpVector& v ){};

  virtual void initTransposedChild(bool dynamic) {};

  virtual void ColumnScale2( OoqpVector& vec, OoqpVector& parentvec ){};
  virtual void RowScale2( OoqpVector& vec, OoqpVector* linkingvec ){};

  virtual void initTransposed(bool dynamic = false) {};
  virtual void deleteTransposed() {};

  virtual void getNnzPerRow(OoqpVectorBase<int>& nnzVec, OoqpVectorBase<int>* linkParent) {};

  virtual void getNnzPerCol(OoqpVectorBase<int>& nnzVec, OoqpVectorBase<int>* linkParent) {};

  virtual void getNnzPerRow(OoqpVectorBase<int>& nnzVec) {};

  virtual void getNnzPerCol(OoqpVectorBase<int>& nnzVec) {};

  virtual void getRowMinMaxVec( bool getMin, bool initializeVec,
        const OoqpVector* colScaleVec, const OoqpVector* colScaleParent, OoqpVector& minmaxVec, OoqpVector* linkParent ){};

  virtual void getColMinMaxVec( bool getMin, bool initializeVec,
        const OoqpVector* rowScaleVec, const OoqpVector* rowScaleParent, OoqpVector& minmaxVec, OoqpVector* minmaxParent ){};

  virtual void getRowMinMaxVec( bool getMin, bool initializeVec,
        const OoqpVector* colScaleVec, OoqpVector& minmaxVec ){};

  virtual void getColMinMaxVec( bool getMin, bool initializeVec,
        const OoqpVector* rowScaleVec, OoqpVector& minmaxVec ){};

  virtual void addRowSums( OoqpVector& sumVec, OoqpVector* linkParent ) {};
  virtual void addColSums( OoqpVector& sumVec, OoqpVector* linkParent ) {};
  virtual void addRowSums( OoqpVector& vec ) {};
  virtual void addColSums( OoqpVector& vec ) {};

  virtual void initStaticStorageFromDynamic(const OoqpVectorBase<int>& rowNnzVec, const OoqpVectorBase<int>& colNnzVec) {};
  virtual void initStaticStorageFromDynamic(const OoqpVectorBase<int>& rowNnzVec, const OoqpVectorBase<int>& colNnzVec, const OoqpVectorBase<int>* rowLinkVec, const OoqpVectorBase<int>* colParentVec) {};

  virtual void freeDynamicStorage() {};

  virtual std::vector<int> get2LinkStartBlocks() const {return std::vector<int>();};

  virtual void updateKLinkVarsCount(std::vector<int>& linkCount) const {};
  virtual void updateKLinkConsCount(std::vector<int>& linkCount) const {};

  virtual void permuteLinkingVars(const std::vector<unsigned int>& permvec) {};
  virtual void permuteLinkingCons(const std::vector<unsigned int>& permvec) {};

  virtual bool isRootNodeInSync() const { return true; };
};


typedef SmartPointer<StochGenMatrix> StochGenMatrixHandle;

#endif
