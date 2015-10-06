/* PIPS-NLP                                                         	*
 * Authors: Nai-Yuan Chiang                      		*
 * (C) 2015 Argonne National Laboratory			*/

#ifndef NLPREDUCEDSPACESOLVER_NLP
#define NLPREDUCEDSPACESOLVER_NLP

#include "SparseSymMatrixHandle.h"
#include "DoubleLinearSolver.h"


class   SparseSymMatrix;
class   DenseSymMatrix;
class   SparseGenMatrix;
class   DenseGenMatrix;
class	SymMatrix;
class	GenMatrix;
class	OoqpVector;
class	SimpleVector;

class DoubleMatrix;


/**
*  u is the decision matrix, x is the others (state var)
*  in this solver, we assume Ax and Ax' are invertible.
*
*  fullMat = 
*    	Hxx	| Hxu	|  		|Ax'	|   Tx'	|			r_x
*	Hux	| Huu 	|   		|Au'	|   Tu'	|			r_u
*		| 	  	|  Hss  	|	|   -I	|	     =	r_s
*	Ax	| Au  	|  		|	|		|			r_y
*	Tx	| Tu  	|  -I		|	|		|			r_z
*
*  this can be reordered as
*  fullMat = 
*    	Hxx	| 		|  Hxu	|Ax'	|   Tx'	|			r_x
*		| Hss 	|   		|	|   -I	|		       r_s
*	Hux	| 	  	|  Huu  	|Au'	|   Tu'	|	   	=	r_u			(1)
*	Ax	| 	  	|  Au	|	|		|			r_y
*	Tx	|  -I 	|  Tu	|	|		|			r_z
*
*
*  let C    = 	(	Ax		) ,      Ac =	( Au ),   Huc = ( Hux  0 ) = Hcu'		(*)
*			( 	Tx     -I	)			( Tu )	
*			     
*  the invert mat of C is 
*  \C^{-1} =  	(	Ax^{-1}			)						
*				(      TxAx^{-1}     -I	)
*
*
* (1)can be rewrited as 
*    	Hc	| 	Hcu	|   C'	|            	r_p
*	Huc	| 	Huu	|   Ac'	|      =    r_u							(2)
*	C	| 	Ac  	| 		|            	r_d
*
*    	Hc	| 	C'	|   Hcu 	|            	r_p
*	C	| 	 	|   Ac 	|           	r_d							(3)
*	Huc	| 	Ac'  	|  Huu 	|      =    r_u							
*
*
*
* similar to the schur complement method, we solve the problem fullMat*\delta_(x,u,y,z)=r by the following sequence:
* (a) solve following equation to get \delta_u
*      ( Huu-Huc*C^{-1}*Ac-Ac'*C^{-T}Hcu+Ac'*C^{-T}*Hc*C^{-1}*Ac ) *\delta_u 
*   = r_u -Ac'*C^{-T}(r_p-Hc*C^{-1}r_d)-Huc*C^{-1}r_d
*
*  here we save rhs_p =  C^{-1}r_d and rhs_d = C^{-T}(r_p-Hc*C^{-1}r_d) 
*
* (b) compute \delta_p = C^{-1}(r_d-Ac*\delta_u) = rhs_p - C^{-1}Ac*\delta_u
*      
* (c) compute \delta_d 	= C^{-T} [ (r_p-Hc*C^{-1}r_d) - (Hcu-Hc*C^{-1}*Ac)*\delta_u   ]
*		or			= C^{-T} ( r_p - Hc*delta_p- Hcu*delta_u)
*		or			= rhs_d - C^{-T}  (Hcu*\delta_u -Hc*C^{-1}*Ac*\delta_u )
*
* the later terms in (b) and (c) are computed by solving
*
*    	Hc	| 	   C'	|	=           	r_p	 	=	Hcu*\delta_u
*	C	| 	  	| 		      	r_d			Ac*\delta_u
*
*
*  note that input matrix fullMat is symmetric, and in a row-wise, lower-triangular form
*  
*/


class ReducedSpaceSolver: public DoubleLinearSolver{

  protected:

	DoubleLinearSolver *Ax_solver;
	DoubleLinearSolver *SC_solver;

	int firstCallFlag;
	int firstSCsolve;

	int doBuildSc;
	// ***************    if we build schur complement
	SparseSymMatrix* kkt_diag;
	DenseSymMatrix*  kktsc;
	DoubleLinearSolver *SparseAug_solver;
	bool firstHxxUpdate,firstHssUpdate, firstAxUpdate,firstTxUpdate;
 	std::map<int,int> LocHxxMap;
 	std::map<int,int> LocHssMap;
 	std::map<int,int> LocAxMap;
	std::map<int,int> LocTxMap;
	//************************************************************

    /*spilit input matrix*/
	int Hxx_Dim, Huu_Dim, Hss_Dim, Ax_Dim_m, Tx_Dim_m, fullMatDim;
	int Hxx_NNz, Huu_NNz, Hss_NNz, Hux_NNz, Ax_NNz, Au_NNz, Tx_NNz, Tu_NNz, I_NNz, fullMatNNz;
	double *Hxx_ele, *Huu_ele, *Hss_ele, *Hux_ele, *Ax_ele, *Au_ele, *Tx_ele, *Tu_ele, *I_ele;	
	int *Hxx_rowBeg, *Huu_rowBeg, *Hss_rowBeg, *Hux_rowBeg, *Ax_rowBeg, *Au_rowBeg, 
		*Tx_rowBeg, *Tu_rowBeg, *I_rowBeg;	
	int *Hxx_colIdx, *Huu_colIdx, *Hss_colIdx, *Hux_colIdx, *Ax_colIdx, *Au_colIdx, 
		*Tx_colIdx, *Tu_colIdx, *I_colIdx;	

	int *Hxx_Full_eleMap, *Huu_Full_eleMap, *Hss_Full_eleMap, *Hux_Full_eleMap, 
		*Ax_Full_eleMap, *Au_Full_eleMap, 
		*Tx_Full_eleMap, *Tu_Full_eleMap, *I_Full_eleMap;

	int *decisionVarIDinFull;
	
	int DecisionVarDim, slackVarDim, stateVarDim, dualYDim;

	SparseSymMatrix *Hxx_Mat, *Huu_Mat, *Hss_Mat, *I_Mat;

	SparseGenMatrix *Hux_Mat, *Hxu_Mat,
					*Ax_Mat, *Au_Mat,
					*Tx_Mat, *Tu_Mat;

	SparseSymMatrix* Msys;

	// number of negative Eig Val, how to find it?
	int localNegaEigVal;
	
    //************************************ 	Methods 		*********************************


    virtual void firstCall();
	virtual int _numericalFact();
	

  	virtual void diagonalChanged( int idiag, int extent ) { assert(0 && "Not implemented"); }
    virtual void solve ( OoqpVector& rhs_ );
	virtual void solve(  GenMatrix& rhs_in);


  public: 

	virtual int getDecisionDim(){return DecisionVarDim;};

	virtual int matrixChanged(){ return this->_numericalFact();};


	ReducedSpaceSolver();
	ReducedSpaceSolver(DoubleMatrix* MatIn, const int decisionVarSize, int* decisionVarID,
								const int fullVarXSize, const int fullVarYSize, const int fullVarSSize);
    ~ReducedSpaceSolver(){};


	private:
	  // this is equal to decisionVarIDinFull
	 int* stateVarIDinFull;

	  //if FullVarIDinDiag01[i] = 	0: 	 not in diag0 or diag1
	  //						pos: 	 in diag1
	  //						neg: 	 in diag0
	  int* FullVarIDinLocal;
//	  int* FullVarIDinHuu;


	  void reducedSpaceJacSolve (OoqpVector* rhs_y_, OoqpVector* rhs_z_);
	  void reducedSpaceJacTransSolve (OoqpVector* rhs_x_, OoqpVector* rhs_s_);
	  void solveDeltaU (OoqpVector& rhs_u, OoqpVector* rhs_Full);
	  void schursolveDeltaU (OoqpVector& rhs_u, OoqpVector* rhs_Full);

	  void schursolveDeltaU_BuildSC (OoqpVector& rhs_, OoqpVector* rhs_Full);
	  void _schursolveDeltaU_BuildSC_firstCall ();

	  void addTermToDenseSchurCompl( DenseSymMatrix& SC);
};


#endif


