/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#ifndef QPGENLINSYS
#define QPGENLINSYS

#include "LinearSystem.h"
#include "DoubleMatrixHandle.h"
#include "OoqpVector.h"
#include "Observer.h"

class Data;
class QpGenData;
class QpGen;
class Variables;
class Residuals;
class DoubleLinearSolver;
class LinearAlgebraPackage;


/** 
 * Linear System solvers for the general QP formulation. This class
 * contains definitions of methods and data common to the sparse and
 * dense special cases of the general formulation. The derived classes
 * QpGenSparseLinsys and QpGenDenseLinsys contain the aspects that are
 * specific to the sparse and dense forms.
 *
 * @see QpGenSparseLinsys
 * @see QpGenDenseLinsys
 *
 * @ingroup QpGen 
 */

class QpGenLinsys : public LinearSystem, public Subject {
protected:
  /** observer pattern for convergence status of BiCGStab when calling solve */
  int bicg_conv_flag;
  int bicg_niterations;

  double bicg_resnorm;
  double bicg_relresnorm;

  int getIntValue(const std::string& s) const override;
  double getDoubleValue(const std::string& s) const override;
  bool getBoolValue(const std::string& s) const override;

  /** stores a critical diagonal matrix as a vector */
  OoqpVector* nomegaInv;

  QpGen * factory;

  /** right-hand side of the system */
  OoqpVector* rhs;

  QpGenLinsys();

  /** dimensions of the vectors in the general QP formulation */
  long long nx, my, mz;

  /** temporary storage vectors */
  OoqpVector* dd, *dq;

  /** index matrices for the upper and lower bounds on x and Cx */
  OoqpVector* ixupp, *icupp, *ixlow, *iclow;

  /** dimensions of the upper and lower bound vectors */
  long long nxupp, nxlow, mcupp, mclow;
  int useRefs;

  /** Work vectors for iterative refinement of the XYZ linear system */
  OoqpVector *sol, *res, *resx, *resy, *resz;
  /** Work vectors for BiCGStab */
  OoqpVector *sol2, *sol3, *res2, *res3, *res4, *res5;

  /// error absorbtion in linear system outer level
  const int outerSolve;
  const int innerSCSolve;

  /// parameters for the bicg solve
  const bool outer_bicg_print_statistics;

  const double outer_bicg_eps;

  const int outer_bicg_max_iter;
  const int outer_bicg_max_normr_divergences;
  const int outer_bicg_max_stagnations;

public:
  QpGenLinsys(  QpGen * factory,
		QpGenData * data,
		LinearAlgebraPackage * la );

  virtual ~QpGenLinsys();


  /** sets up the matrix for the main linear system in "augmented
   * system" form. The actual factorization is performed by a routine
   * specific to either the sparse or dense case.
   *
   * @see QpGenSparseLinsys::factor
   * @see QpGenDenseLinsys::factor
   */
  virtual void factor(Data *prob, Variables *vars);

  /** solves the system for a given set of residuals. Assembles the
   * right-hand side appropriate to the matrix factored in factor,
   * solves the system using the factorization produced there,
   * partitions the solution vector into step components, then
   * recovers the step components eliminated during the block
   * elimination that produced the augmented system form 
   * 
   * @see QpGenSparseLinsys::solveCompressed
   * @see QpGenDenseLinsys::solveCompressed
*/
  virtual void solve(Data *prob, Variables *vars, Residuals *res,
		     Variables *step);

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
			  OoqpVector& ztemp, QpGenData * data );

  /** perform the actual solve using the factors produced in factor.
   *
   * @param rhs on input contains the aggregated right-hand side of
   * the augmented system; on output contains the solution in
   * aggregated form 
   *
   * @see QpGenSparseLinsys::solveCompressed
   * @see QpGenDenseLinsys::solveCompressed
   */
  virtual void solveCompressed( OoqpVector& rhs ) = 0;

  /** places the diagonal resulting from the bounds on x into the
   * augmented system matrix */
  virtual void putXDiagonal( OoqpVector& xdiag ) = 0;

  /** places the diagonal resulting from the bounds on Cx into the
   * augmented system matrix */
  virtual void putZDiagonal( OoqpVector& zdiag ) = 0;

  /** computes the diagonal matrices in the augmented system from the
      current set of variables */
  virtual void computeDiagonals( OoqpVector& dd, OoqpVector& omega,
				 OoqpVector& t,  OoqpVector& lambda,
				 OoqpVector& u,  OoqpVector& pi,
				 OoqpVector& v,  OoqpVector& gamma,
				 OoqpVector& w,  OoqpVector& phi );
 protected:
  virtual void computeResidualXYZ(OoqpVector& sol, 
				  OoqpVector& res, 
				  OoqpVector& solx, 
				  OoqpVector& soly, 
				  OoqpVector& solz,
				  QpGenData* data);
  virtual void matXYZMult( double beta,  OoqpVector& res, 
			   double alpha, OoqpVector& sol, 
			   QpGenData* data,
			   OoqpVector& solx, 
			   OoqpVector& soly, 
			   OoqpVector& solz);

  virtual double matXYZinfnorm(
               QpGenData* data,
               OoqpVector& solx,
               OoqpVector& soly,
               OoqpVector& solz);

  virtual void solveCompressedBiCGStab(OoqpVector& stepx,
				       OoqpVector& stepy,
				       OoqpVector& stepz,
				       QpGenData* data);
  virtual void solveCompressedIterRefin(OoqpVector& stepx,
					OoqpVector& stepy,
					OoqpVector& stepz,
					QpGenData* data);

};

#endif
