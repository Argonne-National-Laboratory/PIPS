/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

/* 2015. Modified by Nai-Yuan Chiang for NLP*/

#ifndef NLPGENLINSYS_NLP
#define NLPGENLINSYS_NLP

#include "LinearSystem.h"
#include "DoubleMatrixHandle.h"
#include "OoqpVector.h"

class Data;
class NlpGenVars;
class NlpGenData;
class NlpGenResiduals;
class NlpGen;
class Variables;
class Residuals;
class DoubleLinearSolver;
class LinearAlgebraPackage;
class RegularizationAlg;

#ifdef WITH_PETSC
#include "petscksp.h"
#endif


/** 
 * Linear System solvers for the general NLP formulation. This class
 * contains definitions of methods and data common to the sparse and
 * dense special cases of the general formulation. The derived classes
 * NlpGenSparseLinsys and NlpGenDenseLinsys contain the aspects that are
 * specific to the sparse and dense forms.
 *
 * @see NlpGenSparseLinsys
 * @see NlpGenDenseLinsys
 *
 * @ingroup NlpGen 
 */

class NlpGenLinsys : public LinearSystem {
protected:
  // use full Q or just diagonal part of Q;
  bool fullQ;

public:

  /** stores a critical diagonal matrix as a vector */
  OoqpVector* nomegaInv;

  /** right-hand side of the system */
  OoqpVector* rhs;

  /** copy of right-hand side of the system, to compute residual */
  OoqpVector* rhs_back;



  NlpGen * factory;


  /** dimensions of the vectors in the general NLP formulation */
  long long nx, my, mz;

  /** primal reg and dual reg **/
  double priReg,dualReg;
 
  
  /** temporary storage vectors */
  OoqpVector* dd, *dq, *temp_diagX, *temp_diagS, *temp_diagY, *temp_diagZ;
  /** diag part: X^{-1}Z */
  OoqpVector* additiveDiag;

  /** index matrices for the upper and lower bounds on x and Cx */
  OoqpVector* ixupp, *icupp, *ixlow, *iclow;


  /** dimensions of the upper and lower bound vectors */
  int nxupp, nxlow, mcupp, mclow;
  int num_slacks;

  int useRefs;

//  double linsysRes;

  /** Work vectors for iterative refinement of the XYZ linear system */
  OoqpVector *sol, *res, *resx, *ress, *resy, *resz;
  /** Work vectors for BiCGStab and iterative refinement */
  OoqpVector *sol2, *res2, *res3, *res4, *res5, 
  			 *sol2Bicg,*res2Bicg, *res3Bicg, *res4Bicg, *res5Bicg;

  // number of krylov iterations
  int KryIter;

  bool allocateSpace;
  
public:
  NlpGenLinsys(  NlpGen * factory,
		NlpGenData * data,
		LinearAlgebraPackage * la );

  NlpGenLinsys() {KryIter=0;allocateSpace=false;};

  ~NlpGenLinsys();


  /** sets up the matrix for the main linear system in "augmented
   * system" form. The actual factorization is performed by a routine
   * specific to either the sparse or dense case.
   *
   * @see NlpGenSparseLinsys::factor
   * @see NlpGenDenseLinsys::factor
   */
  virtual void factor(Data *prob, Variables *vars,RegularizationAlg *RegInfo);
  virtual void factor(Data *prob, Variables *vars);

  /** solves the system for a given set of residuals. Assembles the
   * right-hand side appropriate to the matrix factored in factor,
   * solves the system using the factorization produced there,
   * partitions the solution vector into step components, then
   * recovers the step components eliminated during the block
   * elimination that produced the augmented system form 
   * 
   * @see NlpGenSparseLinsys::solveCompressed
   * @see NlpGenDenseLinsys::solveCompressed
*/
  virtual void solve(Data *prob, Variables *vars, Residuals *res,
		     Variables *step);

  virtual void solve_NTsteps(Data *prob, Variables *vars, Residuals *resids,
		     Variables *Nstep, Variables *Tstep, Variables *NTstep);

  virtual void solve_IterRefine(Data * prob_in, Variables *vars_in,
			Residuals *res_in, Variables *step_in, Residuals *KKT_Resid_in ,Variables *KKT_sol_in);

  /** assembles a single vector object from three given vectors
   *
   * @param rhs (output) final joined vector
   * @param rhs1 (input) first part of rhs
   * @param rhs2 (input) middle part of rhs
   * @param rhs3 (input) last part of rhs
   */
  virtual void joinRHS( OoqpVector& rhs,  OoqpVector& rhs1,
			OoqpVector& rhs2, OoqpVector& rhs3 );

  /** extracts three component vectors from a given aggregated vector.
   *
   * @param vars (input) aggregated vector
   * @param vars1 (output) first part of vars
   * @param vars2 (output) middle part of vars
   * @param vars3 (output) last part of vars
   */
  virtual void separateVars( OoqpVector& vars1, OoqpVector& vars2,
			     OoqpVector& vars3, OoqpVector& vars );

  /** assemble right-hand side of augmented system and call
      solveCompressed to solve it */
  virtual void solveXYZS( OoqpVector& stepx, OoqpVector& stepy,
			  OoqpVector& stepz, OoqpVector& steps,
			  OoqpVector& ztemp, NlpGenData * data );

  /** perform the actual solve using the factors produced in factor.
   *
   * @param rhs on input contains the aggregated right-hand side of
   * the augmented system; on output contains the solution in
   * aggregated form 
   *
   * @see NlpGenSparseLinsys::solveCompressed
   * @see NlpGenDenseLinsys::solveCompressed
   */
  virtual void solveCompressed( OoqpVector& rhs ) = 0;

  virtual void solveCompressedAugXSYZ(OoqpVector& stepx, OoqpVector& steps,
						 OoqpVector& stepy,
						 OoqpVector& stepz,
						 NlpGenData* prob);

  virtual void solveCompressedAugXSYZ_PETSC(OoqpVector& stepx, OoqpVector& steps,
						 OoqpVector& stepy,
						 OoqpVector& stepz,
						 NlpGenData* prob);

  virtual void solveBiCGStab(OoqpVector& stepx, OoqpVector& steps,
					  OoqpVector& stepy, OoqpVector& stepz,
					  NlpGenData* data);


  /** places the diagonal resulting from the bounds on x into the
   * augmented system matrix ---  corresponding to pri variables X */
  virtual void putXDiagonal( OoqpVector& xdiag ) = 0;

  /** places the diagonal resulting from the bounds on Cx into the
   * augmented system matrix ---  corresponding to slack variables S */
  virtual void putSDiagonal( OoqpVector& xdiag ) = 0;

  /** places the diagonal resulting from the bounds on Cx into the
   * augmented system matrix ---  corresponding to slack variables S, z is the dual for s */
  virtual void putZDiagonal( OoqpVector& sdiag ) = 0;

  /** computes the diagonal matrices in the augmented system from the
      current set of variables */
  virtual void computeDiagonals( OoqpVector& dd, OoqpVector& omega,
				 OoqpVector& t,  OoqpVector& lambda,
				 OoqpVector& u,  OoqpVector& pi,
				 OoqpVector& v,  OoqpVector& gamma,
				 OoqpVector& w,  OoqpVector& phi );
  
virtual void solveCompressedBiCGStab(OoqpVector& stepx,
					 OoqpVector& stepy,
					 OoqpVector& stepz,
					 NlpGenData* data);
virtual void solveCompressedIterRefin(OoqpVector& stepx,
				  OoqpVector& stepy,
				  OoqpVector& stepz,
				  NlpGenData* data);

virtual void matXYZMult( double beta,  OoqpVector& res, 
			 double alpha, OoqpVector& sol, 
			 NlpGenData* data,
			 OoqpVector& solx, 
			 OoqpVector& soly, 
			 OoqpVector& solz);

virtual void matXSYZMult( double beta,  OoqpVector& res, 
			 double alpha, OoqpVector& sol, 
			 NlpGenData* data,
			 OoqpVector& solx, 
			 OoqpVector& sols,
			 OoqpVector& soly, 
			 OoqpVector& solz);


/** places the diagonal resulting from the bounds on x into the
 * augmented system matrix ---  for Regularization dual Y */
virtual void putYDualDiagonal( OoqpVector& ydiag ) = 0;


virtual void joinRHSXSYZ( OoqpVector& rhs_in,  OoqpVector& rhs1_in,
		  OoqpVector& rhs2_in, OoqpVector& rhs3_in, OoqpVector& rhs4_in );

virtual void separateVarsXSYZ( OoqpVector& x_in, OoqpVector& s_in, 
		  OoqpVector& y_in, OoqpVector& z_in, OoqpVector& vars_in);  


virtual void factorNoMatChange(Data *prob_in, Variables *vars_in,RegularizationAlg *RegInfo);


virtual void UpdateMatrices( Data * prob_in,int const updateLevel=2){ assert( "not done" && 0 );};

virtual double computeResidual(NlpGenData * data, OoqpVector & res_, OoqpVector & sol_,
			OoqpVector & solx_,OoqpVector & sols_,OoqpVector & soly_,OoqpVector & solz_);

virtual double computeResidual_FullKKT(NlpGenData * data, NlpGenResiduals *res_, NlpGenVars *sol_, NlpGenVars *var_);

virtual double eval_xWx(NlpGenData *prob, NlpGenResiduals *resid, NlpGenVars *steps);
virtual void computeQuantitiesForDualReg( 	NlpGenData *prob, NlpGenVars *vars,
															NlpGenResiduals *resid, NlpGenVars *steps,
														  	double *dualRegQuantities);



virtual void setXDiagonal( OoqpVector& xdiag )=0;
virtual void setSDiagonal( OoqpVector& sdiag )=0;  
virtual void setYDiagonal( OoqpVector& ydiag )=0;  
virtual void setZDiagonal( OoqpVector& zdiag )=0;


virtual void copyXSYZ_fromArray( OoqpVector& vec_xsyz, double* array_in, const int nb_col);
virtual void copyXSYZ_toArray( OoqpVector& vec_xsyz, double* array_in, const int nb_col);



/** 
 *
 * @
 
virtual void schur_Decomp( OoqpVector& rhs );
virtual void schur_Solve( OoqpVector& rhs );
*/

#ifdef WITH_PETSC
KSP mKsp;
Mat LinSysMat_PETSC;
PC	precond_Method;
#endif

  
};

#endif
