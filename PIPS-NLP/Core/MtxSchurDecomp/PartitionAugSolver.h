/* PIPS-NLP                                                         	*
 * Authors: Nai-Yuan Chiang                      		*
 * (C) 2015 Argonne National Laboratory			*/

#ifndef NLPPARTITIONSOLVER
#define NLPPARTITIONSOLVER

#include "SparseSymMatrixHandle.h"
#include "DoubleLinearSolver.h"

#include <map>
#include <vector>
#include <string>


class   SparseSymMatrix;
class   DenseSymMatrix;
class   SparseGenMatrix;
class   DenseGenMatrix;
class	SymMatrix;
class	GenMatrix;
class	OoqpVector;
class	SimpleVector;
class   DoubleMatrix;



/**
*  for each diagonal, after we remove all the 1st stage var
*  Aug fullMat = 
*    	H	|   A'	|							=   	r_p			(1)
*	A	|   0 	|								r_d
* = 
*	x1	|	x2	|	xp	|	y1	|	y2	|	yp
*-----------------------------------------------------------
*    	H1	| 		|  H1p	|	A1'	|   		|			p_1
*		| 	H2 	|  H2p	|		|   A2'	|		       p_2
*	Hp1	| 	Hp2 	|  Hp  	|	A1p'	|   A2p'	|	Ap'	   =	p_p			(1)
*	A1	| 	  	|  A1p	|		|		|			d_1
*		|  	A2 	|  A2p	|		|		|			d_2
*		|  	 	|  Ap	|		|		|			d_p
*
*
*  Reorder it
*	x1	|	y1	|	x2	|	y2	|	xp	|	yp
*-----------------------------------------------------------
*    	H1	| 	A1'	|  		|		|   H1p	|			p11
*	A1	| 	 	|   		|		|   A1p	|		       d_1
*		| 	  	|  H2  	|	A2'	|   H2p	|		   =	p_2			(1)
*		| 	  	|  A2	|		|   A2p	|			d_2
*	H1p'	|  	A1p'	|  H2p'	|  	A2p'	|   Hp	|	Ap'		p_p
*		|  	 	|  		|		|   Ap	|			d_p
*
*
*  can do SC on  	[	Hp	|	Ap' 	]
*				[	Ap	|		]
*
*
*
*/


class PartitionAugSolver: public DoubleLinearSolver{

  protected:

	int nb_part;

	int firstCallFlag;
	int firstSCsolve;	

	/* info of input full matrix*/
	SparseSymMatrix* M_sys;
	int M_Dim;
	int M_Nnz;
	int M_nb_Var;
//	int M_nb_slack;
	int M_nb_Con;
	
	/* info of schur complement matrix*/
	DenseSymMatrix*  dkktSC;	
	DoubleLinearSolver *schur_solver;
	int sc_Dim, sc_Nnz;
	double *sc_ele;
	int *FullIDin_Diag;

	/* info of last diag matrix*/
	SparseSymMatrix* kkt_diag_last;	
	int diag_last_Dim;
	int diag_last_Nnz;
	int* diag_last_rowBeg;	
	int* diag_last_colIdx;	
	double* diag_last_ele;	
	int* diag_last_ele_Map;		
	int diag_last_nb_Var;
	int diag_last_nb_Con;
	int *diag_last_IDinFull;

	
	/* info of diag matrix*/
	std::vector<SparseSymMatrix*> kkt_diag;	
	std::vector<DoubleLinearSolver*> diag_solver;
	std::vector<int> diag_Dim;
	std::vector<int> diag_Nnz;
	std::vector<int*> diag_rowBeg;	
	std::vector<int*> diag_colIdx;	
	std::vector<double*> diag_ele;	
	std::vector<int*> diag_ele_Map;		
	std::vector<int> diag_nb_Var;
	std::vector<int> diag_nb_Con;
	std::vector<int*> diag_IDinFull;

	/* info of bord matrix*/	
	std::vector<SparseGenMatrix*> kkt_bord;	
	std::vector<int> bord_Dim_m;
	std::vector<int> bord_Dim_n;
	std::vector<int> bord_Nnz;
	std::vector<int*> bord_rowBeg;	
	std::vector<int*> bord_colIdx;		
	std::vector<double*> bord_ele;	
	std::vector<int*> bord_ele_Map;	

	

	int *var_Part_idx;
	int *con_Part_idx;

	// number of negative Eig Val, only used when Pardiso is called
	int localNegaEigVal;



	
    //************************************ 	Methods 		*********************************
    virtual void firstCall();
	
	virtual void initializeKKT_Dense();
	virtual int _numericalFact();
	virtual void finalizeKKT();
	virtual void addTermToDenseSchurCompl();	
	virtual void addLnizi(OoqpVector& z0_, OoqpVector& zi_,const int block);
	virtual void LDLtsolve_SC( OoqpVector* lastDiag_Rhs);

  	virtual void diagonalChanged( int idiag, int extent ) { assert(0 && "Not implemented"); }
	    
	virtual void LniTransMult(SparseGenMatrix *borderMat, 
									  SimpleVector& y, double alpha, SimpleVector& x,const int block);
	  


  public: 
	PartitionAugSolver(){};
	PartitionAugSolver(DoubleMatrix* MatIn, int* var_Part_idx_in,
					int *con_Part_idx_in, const int localNegaEigVal_in,
					const int fullVarSize, const int fullConSize, const int nPart);
    ~PartitionAugSolver(){};
	
	virtual int matrixChanged(){ return this->_numericalFact();};
    virtual void solve ( OoqpVector& rhs_ );	

	virtual void solve(  GenMatrix& rhs_in);
  
};















#endif




