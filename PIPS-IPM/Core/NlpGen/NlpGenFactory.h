#ifndef NLPGENFACT
#define NLPGENFACT

#include "ProblemFormulation.h"

class Data;
class Residuals;
class LinearSystem;
class Variables;
class LinearAlgebraPackage;
class OoqpVector;
class NlpInput;

class NlpGenFactory : public ProblemFormulation {

 public:

  virtual Data          * makeData     ( );
  virtual Residuals     * makeResiduals( Data * prob_in );
  virtual Variables     * makeVariables( Data * prob_in );
  virtual LinearSystem  * makeLinsys( Data * prob_in )
  {
      NlpGenData * prob = (NlpGenData *) prob_in;
      int n = nx + my + mz;
      
      SparseSymMatrixHandle Mat( new SparseSymMatrix( n,n + nnzQ
						      + nnzA + nnzC ) );
      
      SimpleVectorHandle v( new SimpleVector(n) );
      v->setToZero();
      Mat->setToDiagonal(*v);
      
      prob->putQIntoAt( *Mat, 0, 0 );
      prob->putAIntoAt( *Mat, nx, 0);
      prob->putCIntoAt( *Mat, nx + my, 0 );
      
      Ma57Solver * solver = new Ma57Solver( Mat );


      //QpGenSparseLinsys can be reused for NLP (this was the idea). 
      //However you need to overwrite 'matrixChanged' method to update the nlp linear system 
      //with the new Hessian, Jacobian(s) and all that stuff.
      return new QpGenSparseLinsys( this, prob, la, Mat, solver );
      
  }

};

#endif
