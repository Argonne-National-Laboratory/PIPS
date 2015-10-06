/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

/* 2015. Modified by Nai-Yuan Chiang for NLP*/

#ifndef NLPGENSPARSELINSYS
#define NLPGENSPARSELINSYS

#include "NlpGenLinsys.h"
#include "SparseSymMatrixHandle.h"
#include <map>

class DoubleLinearSolver;

/** 
 * implements the aspects of the solvers for sparse general NLP
 * formulation that are specific to the sparse case.
 *
 * @ingroup NlpGen 
 */
class NlpGenSparseLinsys : public NlpGenLinsys {
protected:
  SparseSymMatrixHandle Mat;
  DoubleLinearSolver * solver;
public:
  NlpGenSparseLinsys(  NlpGen * factory,
		NlpGenData * data,
		LinearAlgebraPackage * la, SparseSymMatrix * Mat,
		DoubleLinearSolver * solver );

  /** perform the actual solve using the factors produced in factor.
   *
   * @param rhs on input contains the aggregated right-hand side of
   * the augmented system; on output contains the solution in
   * aggregated form
   */
  virtual void solveCompressed( OoqpVector& rhs );

  virtual void putXDiagonal( OoqpVector& xdiag );
  virtual void putSDiagonal( OoqpVector& sdiag );  
  /** places the diagonal resulting from the bounds on x into the
   * augmented system matrix ---  for Regularization dual Y */
  virtual void putYDualDiagonal( OoqpVector& ydiag );  
  virtual void putZDiagonal( OoqpVector& zdiag );

  /** calls NlpGenLinsys::factor to assemble the augmented system
   * matrix, then calls matrixChanged to factor it
   *
   * @see NlpGenLinsys::factor 
   */
  virtual void factor(Data *prob, Variables *vars,RegularizationAlg *RegInfo);
  virtual void factor(Data *prob, Variables *vars);
  
  virtual ~NlpGenSparseLinsys();

  virtual void UpdateMatrices( Data * prob_in, int const updateLevel=2);
  


  bool firstXDiagUpdate,firstSDiagUpdate,firstYDiagUpdate,firstZDiagUpdate;

  virtual void setXDiagonal( OoqpVector& xdiag );
  virtual void setSDiagonal( OoqpVector& sdiag );  
  virtual void setYDiagonal( OoqpVector& ydiag );  
  virtual void setZDiagonal( OoqpVector& zdiag );



  virtual void setAdditiveDiagonal();

//  virtual void factorNoMatChange(Data *prob_in, Variables *vars_in,RegularizationAlg *RegInfo);

  std::map<int,int> xDiagIdxMap;
  std::map<int,int> sDiagIdxMap;
  std::map<int,int> yDiagIdxMap;
  std::map<int,int> zDiagIdxMap;


  bool firstQUpdate,firstAUpdate,firstCUpdate;

  std::map<int,int> QmatIdxMap;
  std::map<int,int> AmatIdxMap;
  std::map<int,int> CmatIdxMap;



};
#endif
