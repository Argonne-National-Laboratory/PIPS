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

template<class FORMULATION, class IPMSOLVER> 
class PIPSIpmInterface 
{
 public:
  PIPSIpmInterface(stochasticInput &in);
  ~PIPSIpmInterface();

  void go();

  double getObjective() const;


  void setPrimalTolerance(double val);
  void setDualTolerance(double val);

  std::vector<double> getFirstStagePrimalColSolution() const{};
  std::vector<double> getSecondStagePrimalColSolution(int scen) const{};
  std::vector<double> getFirstStageDualColSolution() const{};
  std::vector<double> getSecondStageDualColSolution(int scen) const{};
  std::vector<double> getSecondStageDualRowSolution(int scen) const {};
  //more get methods to follow here

 protected:
 
  FORMULATION * factory;
  sData *        data;
  sVars *   vars;
  sResiduals *   resids;

  IPMSOLVER *   solver;

  PIPSIpmInterface() {};
  
};

//----------------------------------------------------------------------
// IMPLEMENTATION
//----------------------------------------------------------------------


template<class FORMULATION, class IPMSOLVER>
PIPSIpmInterface<FORMULATION, IPMSOLVER>::PIPSIpmInterface(stochasticInput &in)
{
  factory = new FORMULATION( in );
  //printf("factory created\n");

  data   = dynamic_cast<sData*>     ( factory->makeData() );
  //printf("data created\n");

  vars   = dynamic_cast<sVars*>     ( factory->makeVariables( data ) );
  //printf("variables created\n");

  resids = dynamic_cast<sResiduals*>( factory->makeResiduals( data ) );
  //printf("resids created\n");

  solver  = new IPMSOLVER( factory, data );
  //printf("solver created\n");
  solver->addMonitor(new StochMonitor( factory ));

}


template<typename FORMULATION, typename IPMSOLVER>
void PIPSIpmInterface<FORMULATION,IPMSOLVER>::go() {

  //s->monitorSelf();
  int result = solver->solve(data,vars,resids);
  
  if ( 0 == result ) {
    double objective = getObjective();
    
    cout << " " << data->nx << " variables, " << data->my  
	 << " equality constraints, " << data->mz << " inequality constraints.\n";
    
    cout << " Iterates: " << solver->iter <<",    Optimal Solution:  " 
	 << objective << endl;
  }
}

template<typename FORMULATION, typename SOLVER>
double PIPSIpmInterface<FORMULATION,SOLVER>::getObjective() const {
  return data->objectiveValue(vars);
}



template<class FORMULATION, class IPMSOLVER>
PIPSIpmInterface<FORMULATION, IPMSOLVER>::~PIPSIpmInterface()
{ 
  delete solver;
  delete resids;
  delete vars;
  delete data;
  delete factory;
};



#endif
