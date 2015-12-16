/* PIPS-NLP                                                         	*
 * Authors: Nai-Yuan Chiang                      		*
 * (C) 2015 Argonne National Laboratory			*/

#ifndef SPSTOCHNLPGENFACTORY_NLP
#define SPSTOCHNLPGENFACTORY_NLP

#include "NlpGen.h"
#include "mpi.h"

class NlpGenData;
class NlpGenStochData;

class NlpGenVars;
class StochInputTree;
class StochTree;
class StochSymMatrix;
class NlpGenResiduals2;
class NlpGenStochVars;
class NlpGenStochLinsys;
class NlpGenStochLinsysRoot;
class NlpGenStochLinsysLeaf;


class NlpInfo;

#include "StochResourcesMonitor.h"

class NlpGenStoch : public NlpGen {
 protected:
  int m_blocks;
  
  int nnzQ, nnzA, nnzC;
  
 public:
  NlpGenStoch( StochInputTree* );
 protected:
  NlpGenStoch( int nx_, int my_, int mz_, int nnzQ_, int nnzA_, int nnzC_ );
  NlpGenStoch();
  
 public:

  virtual ~NlpGenStoch();

  virtual Data  * makeData();
  virtual Data  * makeDataMulti(){return NULL;};
  virtual Data	* makeData(NlpInfo * updateNlp);

  virtual Residuals     * makeResiduals( Data * prob_in );
  virtual Variables     * makeVariables( Data * prob_in );

  virtual LinearSystem* makeLinsys( Data * prob_in );


  virtual void joinRHS( OoqpVector& rhs_in,  OoqpVector& rhs1_in,
			OoqpVector& rhs2_in, OoqpVector& rhs3_in );

  virtual void separateVars( OoqpVector& x_in, OoqpVector& y_in,
			     OoqpVector& z_in, OoqpVector& vars_in );
  
  virtual void joinRHSXSYZ( OoqpVector& rhs_in,  OoqpVector& rhs1_in,
			OoqpVector& rhs2_in, OoqpVector& rhs3_in, OoqpVector& rhs4_in );
  
  virtual void separateVarsXSYZ( OoqpVector& x_in, OoqpVector& s_in, 
			OoqpVector& y_in, OoqpVector& z_in, OoqpVector& vars_in);  

  virtual NlpGenStochLinsysRoot* newLinsysRoot() = 0;
  virtual NlpGenStochLinsysRoot* newLinsysRoot(NlpGenStochData* prob,
					      OoqpVector* dd,OoqpVector* dq,
					      OoqpVector* nomegaInv, OoqpVector* rhs) = 0;

  virtual NlpGenStochLinsysLeaf* newLinsysLeaf();
  virtual NlpGenStochLinsysLeaf* newLinsysLeaf(NlpGenStochData* prob,
					      OoqpVector* dd,OoqpVector* dq,
					      OoqpVector* nomegaInv, OoqpVector* rhs);




  StochTree* tree;
  NlpGenStochData * data;
  //  Variables

  virtual void iterateStarted();
  virtual void iterateEnded();


  NlpGenResiduals2 *resid;
  std::vector<NlpGenStochVars*> registeredVars;
 
  NlpGenStochLinsys* linsys;
  NlpGenStochLinsys* linsys_2;

  StochIterateResourcesMonitor iterTmMonitor;
  double m_tmTotal;

  std::string datarootname;


//  int nxLOri_All, nxUOri_All,nsLOri_All,nsUOri_All;
};

#endif

