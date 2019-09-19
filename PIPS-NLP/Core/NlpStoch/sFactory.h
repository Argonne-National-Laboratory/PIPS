/* PIPS
   Authors: Cosmin Petra
   See license and copyright information in the documentation */

/* 2015. Modified by Nai-Yuan Chiang for NLP*/

#ifndef SPSTOCHFACTORY_NLP
#define SPSTOCHFACTORY_NLP

#include "NlpGen.h"
#include "mpi.h"

class NlpGenData;
class sData;

class NlpGenVars;
class StochInputTree;
class stochasticInput;

class sTree;
class StochSymMatrix;
class sResiduals;
class sVars;
class sLinsys;
class sLinsysRoot;
class sLinsysLeaf;

class NlpInfo;

#include "StochResourcesMonitor.h"

class sFactory : public NlpGen {
 protected:
  int m_blocks;

  long long nnzQ, nnzA, nnzC;
 public:
  sFactory( stochasticInput&, MPI_Comm comm=MPI_COMM_WORLD );

  /** This is a obsolete constructor since it uses sTreeCallbacks to create
   *   data objects
   */
  sFactory( StochInputTree* ){};


 protected:
  sFactory( int nx_, int my_, int mz_, int nnzQ_, int nnzA_, int nnzC_ );
  sFactory();

 public:

  virtual ~sFactory();

  virtual Data  * makeData();
  virtual Data  * makeDataMulti();

  virtual Data	* makeData(NlpInfo *updateNlp);


  virtual Residuals     * makeResiduals( Data * prob_in );
  virtual Variables     * makeVariables( Data * prob_in );

  virtual LinearSystem* makeLinsys( Data * prob_in );


  virtual void joinRHS( OoqpVector& rhs_in,  OoqpVector& rhs1_in,
			OoqpVector& rhs2_in, OoqpVector& rhs3_in );

  virtual void separateVars( OoqpVector& x_in, OoqpVector& y_in,
			     OoqpVector& z_in, OoqpVector& vars_in );

  void joinRHSXSYZ( OoqpVector& rhs_in,  OoqpVector& rhs1_in,
			OoqpVector& rhs2_in, OoqpVector& rhs3_in, OoqpVector& rhs4_in );

  void separateVarsXSYZ( OoqpVector& x_in, OoqpVector& s_in,
			OoqpVector& y_in, OoqpVector& z_in, OoqpVector& vars_in);

  virtual sLinsysRoot* newLinsysRoot() = 0;
  virtual sLinsysRoot* newLinsysRoot(sData* prob,
				     OoqpVector* dd,OoqpVector* dq,
				     OoqpVector* nomegaInv, OoqpVector* rhs,OoqpVector* additiveDiag_) = 0;

  virtual sLinsysLeaf* newLinsysLeaf();
  virtual sLinsysLeaf* newLinsysLeaf(sData* prob,
				     OoqpVector* dd,OoqpVector* dq,
				     OoqpVector* nomegaInv, OoqpVector* rhs, OoqpVector* additiveDiag_);


  sTree * tree;
  sData * data;
  //  Variables

  virtual void iterateStarted();
  virtual void iterateEnded();


  sResiduals *resid;
  std::vector<sVars*> registeredVars;

  sLinsysRoot* linsys;

  StochIterateResourcesMonitor iterTmMonitor;
  double m_tmTotal;


  std::string datarootname;


//  int nxLOri_All, nxUOri_All,nsLOri_All,nsUOri_All;
};

#endif
