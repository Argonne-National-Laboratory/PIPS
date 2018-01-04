/* PIPS-IPM                                                           *
 * Author:  Cosmin G. Petra                                           *
 * (C) 2012 Argonne National Laboratory. See Copyright Notification.  */

#ifndef PIPSIPM_INTERFACE
#define PIPSIPM_INTERFACE

#include "stochasticInput.hpp"

#include "sData.h"
#include "sResiduals.h"
#include "sVars.h"
#include "sTree.h"
#include "StochMonitor.h"
#include "Scaler.h"
#include "ScalerFactory.h"

#include <cstdlib>

template<class FORMULATION, class IPMSOLVER> 
class PIPSIpmInterface 
{
 public:
  PIPSIpmInterface(stochasticInput &in, MPI_Comm = MPI_COMM_WORLD);
  PIPSIpmInterface(StochInputTree* in, MPI_Comm = MPI_COMM_WORLD, ScalerType scaler_type = SCALER_NONE);
  ~PIPSIpmInterface();

  void go();
  double getObjective() const;
  double getFirstStageObjective() const;


  void setPrimalTolerance(double val);
  void setDualTolerance(double val);

  std::vector<double> getFirstStagePrimalColSolution() const;
  std::vector<double> getSecondStagePrimalColSolution(int scen) const;
  //std::vector<double> getFirstStageDualColSolution() const{};
  std::vector<double> getFirstStageDualRowSolution() const;
  //std::vector<double> getSecondStageDualColSolution(int scen) const{};
  std::vector<double> getSecondStageDualRowSolution(int scen) const;
  //more get methods to follow here

  static bool isDistributed() { return true; }

 protected:
 
  FORMULATION * factory;
  sData *        data;
  sVars *        vars;
  sResiduals *   resids;

  Scaler *      scaler;
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

#ifdef TIMING
  int mype;
  MPI_Comm_rank(comm,&mype);
#endif

  factory = new FORMULATION( in, comm);
#ifdef TIMING
  if(mype==0) printf("factory created\n");
#endif

  data   = dynamic_cast<sData*>     ( factory->makeData() );
#ifdef TIMING
  if(mype==0) printf("data created\n");
#endif

  vars   = dynamic_cast<sVars*>     ( factory->makeVariables( data ) );
#ifdef TIMING
  if(mype==0) printf("variables created\n");
#endif

  resids = dynamic_cast<sResiduals*>( factory->makeResiduals( data ) );
#ifdef TIMING
  if(mype==0) printf("resids created\n");
#endif

  solver  = new IPMSOLVER( factory, data );
  solver->addMonitor(new StochMonitor( factory ));
#ifdef TIMING
  if(mype==0) printf("solver created\n");
  //solver->monitorSelf();
#endif

}

template<class FORMULATION, class IPMSOLVER>
PIPSIpmInterface<FORMULATION, IPMSOLVER>::PIPSIpmInterface(StochInputTree* in, MPI_Comm comm, ScalerType scaler_type) : comm(comm)
{

#ifdef TIMING
  int mype;
  MPI_Comm_rank(comm,&mype);
#endif

  factory = new FORMULATION( in, comm);
#ifdef TIMING
  if(mype==0) printf("factory created\n");
#endif

  data   = dynamic_cast<sData*>     ( factory->makeData() );
#ifdef TIMING
  if(mype==0) printf("data created\n");
#endif

  vars   = dynamic_cast<sVars*>     ( factory->makeVariables( data ) );
#ifdef TIMING
  if(mype==0) printf("variables created\n");
#endif

  resids = dynamic_cast<sResiduals*>( factory->makeResiduals( data ) );
#ifdef TIMING
  if(mype==0) printf("resids created\n");
#endif

  solver  = new IPMSOLVER( factory, data );
  solver->addMonitor(new StochMonitor( factory ));
#ifdef TIMING
  if(mype==0) printf("solver created\n");
  //solver->monitorSelf();
#endif

  scaler = ScalerFactory::makeScaler(data, scaler_type);
#ifdef TIMING
  if(mype==0) printf("scaler created\n");
#endif
}


template<typename FORMULATION, typename IPMSOLVER>
void PIPSIpmInterface<FORMULATION,IPMSOLVER>::go() {


  int mype;
  MPI_Comm_rank(comm,&mype);
#ifdef TIMING
  if(0 == mype) cout << "solving ..." << endl;

  if(mype==0) {
    cout << "1st stage " << data->getLocalnx() << " variables, " << data->getLocalmy() 
	 << " equality constraints, " << data->getLocalmz() << " inequality constraints." << endl;
    
    int nscens=data->children.size();
    if(nscens) {
      cout << "2nd stage " << data->children[0]->getLocalnx() << " variables, " 
	   << data->children[0]->getLocalmy() << " equality constraints, " 
	   << data->children[0]->getLocalmz() << " inequality constraints." << endl;
      
      cout << nscens << " scenarios." << endl;
      cout << "Total " << data->getLocalnx()+nscens*data->children[0]->getLocalnx() << " variables, " 
	   << data->getLocalmy()+nscens*data->children[0]->getLocalmy()  << " equality constraints, " 
	   << data->getLocalmz()+nscens*data->children[0]->getLocalmz() << " inequality constraints." << endl;
    }
  }
#endif

  double tmElapsed=MPI_Wtime();

  if( scaler )
     scaler->scale();

  //---------------------------------------------
  int result = solver->solve(data,vars,resids);
  //---------------------------------------------

  // todo: scaler->unscale(data,vars,resids);

  tmElapsed=MPI_Wtime()-tmElapsed;
#ifdef TIMING
  double objective = getObjective();
#endif

  if ( 0 == result && 0 == mype ) {
#ifdef TIMING
    //cout << " " << data->nx << " variables, " << data->my  
    // << " equality constraints, " << data->mz << " inequality constraints.\n";
    
    cout << " Iterates: " << solver->iter <<",    Optimal Solution:  " 
	 << objective << endl;

    cout << "Solve time: " << tmElapsed << " seconds." << endl;

    char *var = getenv("OMP_NUM_THREADS");
    if(var != NULL) {
      int num_threads;
      sscanf( var, "%d", &num_threads );
      cout << "Num threads: " << num_threads << endl;
    }
#endif
  }
}

template<typename FORMULATION, typename SOLVER>
double PIPSIpmInterface<FORMULATION,SOLVER>::getObjective() const {
  return data->objectiveValue(vars);
}


template<typename FORMULATION, typename SOLVER>
double PIPSIpmInterface<FORMULATION,SOLVER>::getFirstStageObjective() const {
  OoqpVector& x = *(dynamic_cast<StochVector&>(*vars->x).vec);
  OoqpVector& c = *(dynamic_cast<StochVector&>(*data->g).vec);
  return c.dotProductWith(x);
}



template<class FORMULATION, class IPMSOLVER>
PIPSIpmInterface<FORMULATION, IPMSOLVER>::~PIPSIpmInterface()
{ 
  delete solver;
  delete resids;
  delete vars;
  delete data;
  delete factory;
  delete scaler;
}


template<class FORMULATION, class IPMSOLVER>
std::vector<double> PIPSIpmInterface<FORMULATION, IPMSOLVER>::getFirstStagePrimalColSolution() const {
	SimpleVector const &v = *dynamic_cast<SimpleVector const*>(dynamic_cast<StochVector const&>(*vars->x).vec);
	return std::vector<double>(&v[0],&v[0]+v.length());
}

template<class FORMULATION, class IPMSOLVER>
std::vector<double> PIPSIpmInterface<FORMULATION, IPMSOLVER>::getSecondStagePrimalColSolution(int scen) const {
	SimpleVector const &v = *dynamic_cast<SimpleVector const*>(dynamic_cast<StochVector const&>(*vars->x).children[scen]->vec);
	//int mype;
	//MPI_Comm_rank(comm,&mype);
	//if (!v.length()) printf("oops, asked for scen %d on proc %d\n", scen, mype);
	//assert(v.length());
	if(!v.length()) 
	  return std::vector<double>(); //this vector is not on this processor
	else
	  return std::vector<double>(&v[0],&v[0]+v.length());
}

template<class FORMULATION, class IPMSOLVER>
std::vector<double> PIPSIpmInterface<FORMULATION, IPMSOLVER>::getFirstStageDualRowSolution() const 
{
  SimpleVector const &y     = *dynamic_cast<SimpleVector const*>( (dynamic_cast<StochVector const&>(*vars->y)).vec );
  SimpleVector const &z     = *dynamic_cast<SimpleVector const*>( (dynamic_cast<StochVector const&>(*vars->z)).vec );
  SimpleVector const &iclow = *dynamic_cast<SimpleVector const*>( (dynamic_cast<StochVector const&>(*vars->iclow)).vec);
  SimpleVector const &icupp = *dynamic_cast<SimpleVector const*>( (dynamic_cast<StochVector const&>(*vars->icupp)).vec);
  
  if(!y.length() && !z.length()) 
    return std::vector<double>(); //this vector is not on this processor
  else {
    std::vector<int> const &map=factory->tree->idx_EqIneq_Map;
    
    std::vector<double> multipliers(map.size());
    for(size_t i=0; i<map.size(); i++) {
      int idx=map[i];
      if(idx<0) {
	//equality
	idx=-idx-1; assert(idx>=0);
	multipliers[i]=y[idx];
      } else {
	//inequality - since, we have z-\lambda+\pi=0, where \lambda is the multiplier for low and
	//\pi is the multiplier for upp, therefore z containts the right multiplier for this row.
	assert(iclow[idx]>0 || icupp[idx]>0);
	multipliers[i]=z[idx];
      }
    }
    return multipliers;    
  }
}


template<class FORMULATION, class IPMSOLVER>
std::vector<double> PIPSIpmInterface<FORMULATION, IPMSOLVER>::getSecondStageDualRowSolution(int scen) const {
  SimpleVector const &y = *dynamic_cast<SimpleVector const*>(dynamic_cast<StochVector const&>(*vars->y).children[scen]->vec);
  SimpleVector const &z = *dynamic_cast<SimpleVector const*>(dynamic_cast<StochVector const&>(*vars->z).children[scen]->vec);
  SimpleVector const &iclow = *dynamic_cast<SimpleVector const*>(dynamic_cast<StochVector const&>(*vars->iclow).children[scen]->vec);
  SimpleVector const &icupp = *dynamic_cast<SimpleVector const*>(dynamic_cast<StochVector const&>(*vars->icupp).children[scen]->vec);
  //assert(v.length());
  if(!y.length() && !z.length()) 
    return std::vector<double>(); //this vector is not on this processor
  else {
    std::vector<int> const &map=factory->tree->children[scen]->idx_EqIneq_Map;
    
    std::vector<double> multipliers(map.size());
    for(size_t i=0; i<map.size(); i++) {
      int idx=map[i];
      if(idx<0) {
	//equality
	idx=-idx-1; assert(idx>=0);
	multipliers[i]=y[idx];
      } else {
	//inequality - since, we have z-\lambda+\pi=0, where \lambda is the multiplier for low and
	//\pi is the multiplier for upp, therefore z containts the right multiplier for this row.
	assert(iclow[idx]>0 || icupp[idx]>0);
	multipliers[i]=z[idx];
      }
    }
    return multipliers;
  }
  //return std::vector<double>(&v[0],&v[0]+v.length());
}

#endif
