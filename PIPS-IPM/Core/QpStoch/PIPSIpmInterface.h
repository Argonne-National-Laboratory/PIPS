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
  PIPSIpmInterface(stochasticInput &in, MPI_Comm = MPI_COMM_WORLD);
  ~PIPSIpmInterface();

  void go();

  double getObjective() const;


  void setPrimalTolerance(double val);
  void setDualTolerance(double val);

  std::vector<double> getFirstStagePrimalColSolution() const;
  std::vector<double> getSecondStagePrimalColSolution(int scen) const;
  //std::vector<double> getFirstStageDualColSolution() const{};
  //std::vector<double> getSecondStageDualColSolution(int scen) const{};
  std::vector<double> getSecondStageDualRowSolution(int scen) const;
  //more get methods to follow here

  static bool isDistributed() { return true; }

 protected:
 
  FORMULATION * factory;
  sData *        data;
  sVars *   vars;
  sResiduals *   resids;

  IPMSOLVER *   solver;

  PIPSIpmInterface() {};
  MPI_Comm comm;
  
};

//----------------------------------------------------------------------
// IMPLEMENTATION
//----------------------------------------------------------------------


template<class FORMULATION, class IPMSOLVER>
PIPSIpmInterface<FORMULATION, IPMSOLVER>::PIPSIpmInterface(stochasticInput &in, MPI_Comm comm) : comm(comm)
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
  //solver->addMonitor(new StochMonitor( factory ));

}


template<typename FORMULATION, typename IPMSOLVER>
void PIPSIpmInterface<FORMULATION,IPMSOLVER>::go() {
  int mype;
  MPI_Comm_rank(comm,&mype);

  solver->monitorSelf();

  double tmElapsed=MPI_Wtime();
  //---------------------------------------------
  int result = solver->solve(data,vars,resids);
  //---------------------------------------------
  tmElapsed=MPI_Wtime()-tmElapsed;
  
  double objective = getObjective();
  if ( 0 == result && 0 == mype ) {
    
    cout << " " << data->nx << " variables, " << data->my  
	 << " equality constraints, " << data->mz << " inequality constraints.\n";
    
    cout << " Iterates: " << solver->iter <<",    Optimal Solution:  " 
	 << objective << endl;

    cout << "Solve time: " << tmElapsed << " seconds." << endl;

    char *var = getenv("OMP_NUM_THREADS");
    if(var != NULL) {
      int num_threads;
      sscanf( var, "%d", &num_threads );
      cout << "Num threads: " << num_threads << endl;
    }
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
}


template<class FORMULATION, class IPMSOLVER>
std::vector<double> PIPSIpmInterface<FORMULATION, IPMSOLVER>::getFirstStagePrimalColSolution() const {
	SimpleVector const &v = *dynamic_cast<SimpleVector const*>(dynamic_cast<StochVector const&>(*vars->x).vec);
	return std::vector<double>(&v[0],&v[0]+v.length());
}

template<class FORMULATION, class IPMSOLVER>
std::vector<double> PIPSIpmInterface<FORMULATION, IPMSOLVER>::getSecondStagePrimalColSolution(int scen) const {
	SimpleVector const &v = *dynamic_cast<SimpleVector const*>(dynamic_cast<StochVector const&>(*vars->x).children[scen]->vec);
	int mype;
	MPI_Comm_rank(comm,&mype);
	if (!v.length()) printf("oops, asked for scen %d on proc %d\n", scen, mype);
	assert(v.length());
	return std::vector<double>(&v[0],&v[0]+v.length());
}


/* 
TODO: This is a quick hack!
We only take the equality multipliers because these are what we need for now.
Eventually need to keep a map of inequalities and equalities
*/
template<class FORMULATION, class IPMSOLVER>
std::vector<double> PIPSIpmInterface<FORMULATION, IPMSOLVER>::getSecondStageDualRowSolution(int scen) const {
	SimpleVector const &v = *dynamic_cast<SimpleVector const*>(dynamic_cast<StochVector const&>(*vars->y).children[scen]->vec);
	assert(v.length());
	return std::vector<double>(&v[0],&v[0]+v.length());
}

#endif
