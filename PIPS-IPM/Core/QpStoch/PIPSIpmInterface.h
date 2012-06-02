/* PIPS-IPM                                                           *
 * Author:  Cosmin G. Petra                                           *
 * (C) 2012 Argonne National Laboratory. See Copyright Notification.  */

#ifndef PIPSIPM_INTERFACE
#define PIPSIPM_INTERFACE

#include "stochasticInput.hpp"

#include "sData.h"
#include "sResiduals.h"
#include "sVars.h"
#include "StochMonitor.h"

template<class IPMSOLVER, class FORMULATION> 
class PIPSIpmInterface 
{
 public:
  PIPSIpmInterface(stochasticInput &in);
  ~PIPSIpmInterface();

  void go();

  double getObjective() const;
  //solverState getStatus();

  void setPrimalTolerance(double val);
  void setDualTolerance(double val);

  std::vector<double> getFirstStagePrimalColSolution() const;
  std::vector<double> getSecondStagePrimalColSolution(int scen) const;
  std::vector<double> getFirstStageDualColSolution() const;
  std::vector<double> getSecondStageDualColSolution(int scen) const;
  //more get methods to follow here

 protected:
  //count the number eq and ineq
  void getNum1stStgEqIneq(int& my, int &mz);
  void getNum2ndStgEqIneq(int scens, int& my, int &mz);
  
  FORMULATION* factory;
  IPMSOLVER* solver;

  PIPSIpmInterface() {};
  
};

//----------------------------------------------------------------------
// IMPLEMENTATION
//----------------------------------------------------------------------


template<class FORMULATION, class IPMSOLVER>
PIPSIpmInterface<FORMULATION, IPMSOLVER>::PIPSIpmInterface(stochasticInput &in)
{
  

  FORMULATION * factory = new FORMULATION( in );
  sData *       data    = factory->makeData();
  sVars *       vars    = factory->makeVariables( data );
  sResiduals *  resids  = factory->makeResiduals( data );

  IPMSOLVER*    solver  = new IPMSOLVER( factory, data );

  solver->addMonitor(new StochMonitor( factory ));
}

template<class FORMULATION, class IPMSOLVER>
PIPSIpmInterface<FORMULATION, IPMSOLVER>::~PIPSIpmInterface()
{ };



#endif
