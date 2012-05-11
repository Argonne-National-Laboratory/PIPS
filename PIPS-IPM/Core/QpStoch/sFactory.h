/* PIPS
   Authors: Cosmin Petra
   See license and copyright information in the documentation */

#ifndef SPSTOCHFACTORY
#define SPSTOCHFACTORY

#include "QpGen.h"
#include "mpi.h"

class QpGenData;
class QpGenStochData;

class QpGenVars;
class StochInputTree;
class StochTree;
class StochSymMatrix;
class QpGenResiduals2;
class sVars;
class sLinsys;
class sLinsysRoot;
class sLinsysLeaf;
#include "StochResourcesMonitor.h"

class sFactory : public QpGen {
 protected:
  int m_blocks;
  
  int nnzQ, nnzA, nnzC;
  
 public:
  sFactory( StochInputTree* );
 protected:
  sFactory( int nx_, int my_, int mz_, int nnzQ_, int nnzA_, int nnzC_ );
  sFactory();
  
 public:

  virtual ~sFactory();

  virtual Data  * makeData();

  virtual Residuals     * makeResiduals( Data * prob_in );
  virtual Variables     * makeVariables( Data * prob_in );

  virtual LinearSystem* makeLinsys( Data * prob_in );


  virtual void joinRHS( OoqpVector& rhs_in,  OoqpVector& rhs1_in,
			OoqpVector& rhs2_in, OoqpVector& rhs3_in );

  virtual void separateVars( OoqpVector& x_in, OoqpVector& y_in,
			     OoqpVector& z_in, OoqpVector& vars_in );

  virtual sLinsysRoot* newLinsysRoot() = 0;
  virtual sLinsysRoot* newLinsysRoot(QpGenStochData* prob,
				     OoqpVector* dd,OoqpVector* dq,
				     OoqpVector* nomegaInv, OoqpVector* rhs) = 0;
  
  virtual sLinsysLeaf* newLinsysLeaf();
  virtual sLinsysLeaf* newLinsysLeaf(QpGenStochData* prob,
				     OoqpVector* dd,OoqpVector* dq,
				     OoqpVector* nomegaInv, OoqpVector* rhs);
  

  StochTree* tree;
  QpGenStochData * data;
  //  Variables

  virtual void iterateStarted();
  virtual void iterateEnded();

  QpGenResiduals2 *resid;
  vector<sVars*> registeredVars;
 
  sLinsysRoot* linsys;

  StochIterateResourcesMonitor iterTmMonitor;
  double m_tmTotal;
};

#endif
