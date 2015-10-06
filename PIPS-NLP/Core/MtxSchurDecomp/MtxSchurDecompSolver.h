/* PIPS-NLP                                                         	*
 * Authors: Nai-Yuan Chiang                      		*
 * (C) 2015 Argonne National Laboratory			*/

#ifndef NLPSCHURSOLVER_NLP
#define NLPSCHURSOLVER_NLP

#include "SparseSymMatrixHandle.h"
#include "DoubleLinearSolver.h"


class   SparseSymMatrix;
class   DenseSymMatrix;
class   SparseGenMatrix;
class   DenseGenMatrix;
class	SymMatrix;
class	OoqpVector;
class	SimpleVector;


enum SCHURSOLVER_INPUT_TYPE_0{
  SCHURSOLVER_INPUT_GENERIC   = 0,
  SCHURSOLVER_INPUT_SYMMETRIC = 1
};

enum SCHURSOLVER_INPUT_TYPE_1{
  SCHURSOLVER_INPUT_DENSE 	= 0,
  SCHURSOLVER_INPUT_SPARSE  = 10
};

enum SCHURSOLVER_INPUT_TYPE_2{
  SCHURSOLVER_INPUT_FULL	= 0,
  SCHURSOLVER_INPUT_UPPER  = 100,
  SCHURSOLVER_INPUT_LOWER  = 200  
};

enum SCHURSOLVER_INPUT_TYPE_3{
  SCHURSOLVER_INPUT_ROWWISE	 = 0,
  SCHURSOLVER_INPUT_COLWISE  = 1000,
};

/*	the type of matrix
** 'R': Row-wise; 'C': Column-wise;
** 'F': full matrix; 'U': upper triangular; 'L': lower triangular
** 'D':Dense; 'S': Sparse;
** 'G': general matrix; 'S': symmetric;
*/
enum SCHURSOLVER_INPUT_TYPE{
  INPUT_RFDG	= 0,
  INPUT_RFDS	= 1,
  INPUT_RFSG	= 10,
  INPUT_RFSS	= 11,
  INPUT_RUDG	= 100,
  INPUT_RUDS	= 101,
  INPUT_RUSG	= 110,
  INPUT_RUSS	= 111,
  INPUT_RLDG	= 200,
  INPUT_RLDS	= 201,
  INPUT_RLSG	= 210,
  INPUT_RLSS	= 211,			//only this has been implemented
  INPUT_CFDG	= 1000,
  INPUT_CFDS	= 1001,
  INPUT_CFSG	= 1010,
  INPUT_CFSS	= 1011,
  INPUT_CUDG	= 1100,
  INPUT_CUDS	= 1101,
  INPUT_CUSG	= 1110,
  INPUT_CUSS	= 1111,
  INPUT_CLDG	= 1200,
  INPUT_CLDS	= 1201,
  INPUT_CLSG	= 1210,
  INPUT_CLSS	= 1211   
};


class DoubleMatrix;


/**
*  fullMat = 
*    	diag1Mat 		| borderMat
*	----------------------
*     borderMat'  	| diag0Mat
*
*
*
*/


class MtxSchurDecompSolver: public DoubleLinearSolver{

  protected:

	DoubleLinearSolver *solver_schur;
	DoubleLinearSolver *solver_diag1;
//	DoubleLinearSolver *solver_border;

	int requireUpdate;
	int firstCallFlag;

	/*  the type of matrix*/
  	int mtxCase;


    /*spilit input matrix*/
	int borderDim_m, borderDim_n, diag0MatDim, diag1MatDim, fullMatDim;
	int borderMatNNz, diag0MatNNz, diag1MatNNz, fullMatNNz;	
	double *borderMat_ele, *diag0Mat_ele, *diag1Mat_ele;
	int *borderMat_rowBeg, *diag0Mat_rowBeg, *diag1Mat_rowBeg;
	int *borderMat_colIdx, *diag0Mat_colIdx, *diag1Mat_colIdx;

	int *Border_Full_eleMap;
	int *Diag1_Full_eleMap;
	int *Diag0_Full_eleMap;

	int *schurVarIDinFull;
	
	/* schur complement matrix*/
	int schurMatDim, schurMatNNz;
	double *schurMat_ele;

//	SymMatrix* kktSC;
	DenseSymMatrix* dkktSC;	

    SparseSymMatrix *diag0Mat, *diag1Mat;
	SparseGenMatrix *borderMat;
	
	SparseSymMatrix* Msys;

	// number of negative Eig Val, only used when Pardiso is called
	int localNegaEigVal;
	
    //************************************ 	Methods 		*********************************
    void newMtxSchurDecompSolver(DoubleMatrix* MatIn, const int schurDim, int *inputMatType);


    virtual void firstCall();

    // see above for the definition of Indices
    virtual void _buildDenseSchurMat_RLSS(DoubleMatrix* MatIn, int *schurVarIDX){};
    virtual void _firstCall_RLSS();

	virtual int _numericalFact_RLSS();
	

	virtual void Lsolve( OoqpVector* diag0_Rhs, OoqpVector* diag1_Rhs );
	virtual void Dsolve( OoqpVector* diag0_Rhs, OoqpVector* diag1_Rhs );
	virtual void Ltsolve( OoqpVector* diag0_Rhs, OoqpVector* diag1_Rhs );


//	void myAtPutZeros(DenseSymMatrix* mat);
//    void myAtPutZeros(DenseSymMatrix* mat, int row, int col, int rowExtent, int colExtent);


  	virtual void diagonalChanged( int idiag, int extent ) { assert(0 && "Not implemented"); }
    virtual void solve ( OoqpVector& x );


  public: 

	virtual int getSchurDim(){return schurMatDim;};
	virtual int getSchurNNz(){return schurMatNNz;};	
	
  	virtual void buildDenseSchurMat(DoubleMatrix* MatIn, int *schurVarIDX);
	virtual int matrixChanged();

	virtual void initializeKKT_Dense(DenseSymMatrix* dkktSC);
	virtual void addTermToDenseSchurCompl(SparseSymMatrix *DiagMat, SparseGenMatrix* BordMat, DenseSymMatrix& SC);
	virtual void addLnizi(SparseGenMatrix *borderMat,OoqpVector& z0_, OoqpVector& zi_);
	virtual void LniTransMult(SparseGenMatrix *borderMat, 
									SimpleVector& y, double alpha, SimpleVector& x);
	virtual void finalizeKKT(SparseSymMatrix *diag0Mat, DenseSymMatrix * kktd );
//	virtual void solveCompressed( OoqpVector& rhs_ );

	MtxSchurDecompSolver();
//	MtxSchurDecompSolver(DoubleMatrix* MatIn, const int schurDim, int *inputMatType, int* schurVarIDX);
	MtxSchurDecompSolver(DoubleMatrix* MatIn, const int schurDim, int *inputMatType, int* schurVarIDX,
								const int localNegaEigVal_in);
    ~MtxSchurDecompSolver(){};

	virtual void solve1stVarOnly(OoqpVector& rhs_, OoqpVector& sol_SC);


	private:
	  // this is equal to schurVarConIDX
//	 int* diag0VarIDinFull;
	 int* diag1VarIDinFull;

	  //if FullVarIDinDiag01[i] = 	0: 	 not in diag0 or diag1
	  //						pos: 	 in diag1
	  //						neg: 	 in diag0
	  int* FullVarIDinDiag01;
	
};






int _findIDX(const int idinFullMat, const int vlength, int *_vector);




#endif


