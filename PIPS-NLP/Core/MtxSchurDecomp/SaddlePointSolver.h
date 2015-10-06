/* PIPS-NLP                                                         	*
 * Authors: Nai-Yuan Chiang                      		*
 * (C) 2015 Argonne National Laboratory			*/

#ifndef NLPSADDLEPOINTSOLVER_NLP
#define NLPSADDLEPOINTSOLVER_NLP

#include "SparseSymMatrixHandle.h"
#include "DoubleLinearSolver.h"


class   SparseSymMatrix;
class   DenseSymMatrix;
class   SparseGenMatrix;
class   DenseGenMatrix;
class	SymMatrix;
class	OoqpVector;
class	SimpleVector;

class DoubleMatrix;


/**
*  u is the decision matrix, x is the others (state var)
*  in this solver, we assume Ax and Ax' are invertible.
*  
*  all the u var have been removed.
*
*  fullMat = 
*    	Hxx	| 		|  Ax'	|   Tx'	|			r_x
*		| Hss 	|  		|   -I	|     =   	 	r_s		(1)
*	Ax	| 	  	|(reg_d)	|		|			r_y
*	Tx	|  -I 	|  		|(reg_d)	|			r_z
*
*
*  reg_d is diagonal and ==0, we do not need it.
*
*  let C    = 	(	Ax		) ,      Hc =	( Hxx        0   )	(*)
*			( 	Tx     -I	)			(  0  	Hss )	
*			     
*
*
* (1)can be rewrited as 
*    	Hc	| 	   C'	|	=           	r_p		(2)
*	C	| 	  	| 		      	r_d
*
*  the invert mat of C is 
*  \C^{-1} =  	(	Ax^{-1}			)						
*				(      TxAx^{-1}     -I	)
*
*
* similar to the schur complement method, we solve the problem fullMat*\delta_(x,s,y,z)=r by the following sequence:
* (a) compute \delta_p = C^{-1} (r_d)
*      
* (b) compute \delta_d = C^{-T} ( r_p - Hc*delta_p)
*
*  note that input matrix fullMat is symmetric, and in a row-wise, lower-triangular form
*  
*/


class SaddlePointSolver: public DoubleLinearSolver{

  protected:

	DoubleLinearSolver *Ax_solver;

	int firstCallFlag;
	int firstSCsolve;

    /*spilit input matrix*/
	int Hxx_Dim, Hss_Dim, Ax_Dim_m, Tx_Dim_m, fullMatDim;
	int Hxx_NNz, Hss_NNz, Ax_NNz, Tx_NNz, fullMatNNz;
	double *Hxx_ele, *Hss_ele, *Ax_ele, *Tx_ele;	
	int *Hxx_rowBeg, *Hss_rowBeg, *Ax_rowBeg, *Tx_rowBeg;	
	int *Hxx_colIdx, *Hss_colIdx, *Ax_colIdx, *Tx_colIdx;	

	int *Hxx_Full_eleMap, *Hss_Full_eleMap, *Ax_Full_eleMap, *Tx_Full_eleMap;
	
	int locnx, locns, locmy;

	SparseSymMatrix *Hxx_Mat, *Hss_Mat;

	SparseGenMatrix *Ax_Mat, *Tx_Mat;

	SparseSymMatrix* Msys;

	// number of negative Eig Val, how to find it?
//	int localNegaEigVal;
	
    //************************************ 	Methods 		*********************************

    virtual void firstCall();
	virtual int _numericalFact();
	
  	virtual void diagonalChanged( int idiag, int extent ) { assert(0 && "Not implemented"); }
    virtual void solve ( OoqpVector& rhs_ );
    virtual void solve(GenMatrix& rhs_in);

  public: 

	virtual int matrixChanged(){ return this->_numericalFact();};


	SaddlePointSolver();
	SaddlePointSolver(DoubleMatrix* MatIn, 
					const int localXSize, const int localYSize, const int localZSize);
    ~SaddlePointSolver(){};


	private:

	  void SaddlePointJacSolve (OoqpVector* rhs_y_, OoqpVector* rhs_z_);
	  void SaddlePointJacTransSolve (OoqpVector* rhs_x_, OoqpVector* rhs_s_);
};








#endif




