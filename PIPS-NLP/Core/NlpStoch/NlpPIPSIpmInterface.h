/* PIPS-IPM                                                       *
 * Author:  Cosmin G. Petra & Nai-Yuan Chiang                     *
 * (C) 2012-2015 Argonne National Laboratory.                     */

#ifndef NLPPIPSIPM_INTERFACE
#define NLPPIPSIPM_INTERFACE

#include <cstdlib>

#include "stochasticInput.hpp"

#include "sData.h"
#include "sResiduals.h"
#include "sVars.h"
#include "StochMonitor.h"
#include "../../global_var.h"

using namespace std;

template<class FORMULATION, class IPMSOLVER,class UPDATENLP> 
class NlpPIPSIpmInterface 
{
 public:
  NlpPIPSIpmInterface(stochasticInput &in, MPI_Comm = MPI_COMM_WORLD);
  ~NlpPIPSIpmInterface();

  int go(int addslack=0);

  double getObjective() const;
  void computeProblemSize(int&, int&);
  double getFirstStageObjective() const;


  void setPrimalTolerance(double val);
  void setDualTolerance(double val);

  std::vector<double> getFirstStagePrimalColSolution() const;
  std::vector<double> getSecondStagePrimalColSolution(int scen) const;
  //std::vector<double> getFirstStageDualColSolution() const{};
  //std::vector<double> getSecondStageDualColSolution(int scen) const{};
  std::vector<double> getSecondStageDualRowSolution(int scen) const;
  //more get methods to follow here

  static bool isDistributed() { return true; }
  
  void writeSolution() const;
 protected:
 
  FORMULATION * factory;
  sData *        data;
  sVars *   vars;
  sResiduals *   resids;

  IPMSOLVER *   solver;

  UPDATENLP * updateNlpInfo;

  NlpPIPSIpmInterface() {};
  MPI_Comm comm;
  
};

//----------------------------------------------------------------------
// IMPLEMENTATION
//----------------------------------------------------------------------

template<class FORMULATION, class IPMSOLVER, class UPDATENLP>
NlpPIPSIpmInterface<FORMULATION, IPMSOLVER,UPDATENLP>::NlpPIPSIpmInterface(stochasticInput &in, MPI_Comm comm) : comm(comm)
{
  int mype; MPI_Comm_rank(comm,&mype);
  MESSAGE("enter constructor NlpPIPSIpmInterface "<<mype);
#ifdef TIMING
		double tGenTime=MPI_Wtime();
		double tGenTime2,tGenTime3;
#endif

  factory = new FORMULATION( in, comm);
  MESSAGE("factory created "<<comm);
#ifdef TIMING
  tGenTime2 = MPI_Wtime();
  MPI_Barrier(MPI_COMM_WORLD);
  
  if(0==mype){
		  cout << " " << endl;
		  cout << "FORMULATION Generation Time:   " << tGenTime2-tGenTime << endl;
		  cout << " " << endl;
  }
  if(mype==0) printf("factory created\n");
#endif

//  updateNlpInfo = new UPDATENLP();
//  data   = dynamic_cast<sData*>     ( factory->makeData(updateNlpInfo) );

  data   = dynamic_cast<sData*>	  ( factory->makeData() );
  MESSAGE("data created "<<comm);
  if(in.useInputDate == 0)
    updateNlpInfo = new UPDATENLP( data);
  else
  	updateNlpInfo = new UPDATENLP( data, in);
  MESSAGE("updateNlpInfo created");

  vars   = dynamic_cast<sVars*>     ( factory->makeVariables( data ) );

  resids = dynamic_cast<sResiduals*>( factory->makeResiduals( data ) );

  solver  = new IPMSOLVER( factory, data );
  solver->addMonitor(new StochMonitor( factory ));
}


template<typename FORMULATION, typename IPMSOLVER,typename UPDATENLP>
int NlpPIPSIpmInterface<FORMULATION,IPMSOLVER,UPDATENLP>::go(int addSlack) {


  int mype;
  MPI_Comm_rank(comm,&mype);

  if(0 == mype) cout << "solving ..." << endl;

  if(mype==0) {
    cout<< "1st stage " << data->getLocalnx() << " variables, " << data->getLocalmy()
        << " equality constraints, " << data->getLocalmz() << " inequality constraints." << endl;
  }
  
  int nscens=data->children.size();
  if(nscens) {
  	if(mype==0) {
      cout<< "2nd stage (use 1st scenario): " << data->children[0]->getLocalnx() << " variables, "
          << data->children[0]->getLocalmy() << " equality constraints, "
          << data->children[0]->getLocalmz() << " inequality constraints." << endl;
  	 
      std::cout << nscens << " scenarios." << endl;
  	}
    int sum_var=0,sum_icon=0,sum_econ=0,total_var=0,total_icon=0,total_econ=0;
    for(int j=0;j<nscens;j++){
      sum_var += data->children[j]->getLocalnx();
      sum_econ += data->children[j]->getLocalmy();
      sum_icon += data->children[j]->getLocalmz();
    }

    MPI_Allreduce(&sum_var, &total_var, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&sum_econ, &total_econ, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&sum_icon, &total_icon, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    if(mype==0) {
        std::cout << "Total " << data->getLocalnx() + total_var << " variables, "
       << data->getLocalmy()+total_econ  << " equality constraints, "
       << data->getLocalmz()+total_icon << " inequality constraints. \n" << endl;
    }
  }

  double tmElapsed=MPI_Wtime();
  //---------------------------------------------
  int result;
  if(addSlack == 0) 
  	result = solver->solve(data,vars,resids);
  else{
  	printf("Wrong! Cannot reach here! STOP"); 
  	exit(1);
  }
  //---------------------------------------------
  tmElapsed=MPI_Wtime()-tmElapsed;
#ifdef TIMING
  double objective = getObjective();
#endif

  if ( 0 == result && 0 == mype ) {
#ifdef TIMING
    //cout << " " << data->nx << " variables, " << data->my  
    // << " equality constraints, " << data->mz << " inequality constraints.\n";
    
    std::cout<< " Iterates: " << solver->iter <<",    Optimal Solution:  "
             << objective << endl;

    std::cout << "Solve time: " << tmElapsed << " seconds." << endl;

    char *var = getenv("OMP_NUM_THREADS");
    if(var != NULL) {
      int num_threads;
      sscanf( var, "%d", &num_threads );
      std::cout << "Num threads: " << num_threads << endl;
    }
#endif
    
  }

	this->writeSolution();
	return result;
}

template<typename FORMULATION, typename SOLVER, typename UPDATENLP>
double NlpPIPSIpmInterface<FORMULATION,SOLVER,UPDATENLP>::getObjective() const {
  //return data->objectiveValue(vars);
  return updateNlpInfo->ObjValue(vars);
}

template<typename  FORMULATION, typename SOLVER, typename UPDATENLP>
void NlpPIPSIpmInterface<FORMULATION,SOLVER,UPDATENLP>::computeProblemSize(int& nvar, int& ncon){
  int nscens=data->children.size();
  int sum_var=0,sum_icon=0,sum_econ=0,total_var=0,total_icon=0,total_econ=0;
  for(int j=0;j<nscens;j++){
    sum_var += data->children[j]->getLocalnx();
    sum_econ += data->children[j]->getLocalmy();
    sum_icon += data->children[j]->getLocalmz();
  }

  MPI_Allreduce(&sum_var, &total_var, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&sum_econ, &total_econ, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&sum_icon, &total_icon, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  nvar = total_var;
  ncon = total_econ + total_icon;
}

template<typename FORMULATION, typename SOLVER, typename UPDATENLP>
double NlpPIPSIpmInterface<FORMULATION,SOLVER,UPDATENLP>::getFirstStageObjective() const {
  OoqpVector& x = *(dynamic_cast<StochVector&>(*vars->x).vec);
  OoqpVector& c = *(dynamic_cast<StochVector&>(*data->grad).vec);
  return c.dotProductWith(x);
}



template<class FORMULATION, class IPMSOLVER, class UPDATENLP>
NlpPIPSIpmInterface<FORMULATION, IPMSOLVER,UPDATENLP>::~NlpPIPSIpmInterface()
{ 
  delete solver;
  delete resids;
  delete vars;
  delete data;
  delete factory;
  if(updateNlpInfo) delete updateNlpInfo;
}


template<class FORMULATION, class IPMSOLVER, class UPDATENLP>
std::vector<double> 
NlpPIPSIpmInterface<FORMULATION, IPMSOLVER,UPDATENLP>::getFirstStagePrimalColSolution() const {
	SimpleVector const &v = *dynamic_cast<SimpleVector const*>(dynamic_cast<StochVector const&>(*vars->x).vec);
	return std::vector<double>(&v[0],&v[0]+v.length());
}

template<class FORMULATION, class IPMSOLVER, class UPDATENLP>
std::vector<double> 
NlpPIPSIpmInterface<FORMULATION, IPMSOLVER,UPDATENLP>::getSecondStagePrimalColSolution(int scen) const {
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


/* 
TODO: This is a quick hack!
We only take the equality multipliers because these are what we need for now.
Eventually need to keep a map of inequalities and equalities
*/
template<class FORMULATION, class IPMSOLVER, class UPDATENLP>
std::vector<double> 
NlpPIPSIpmInterface<FORMULATION, IPMSOLVER,UPDATENLP>::getSecondStageDualRowSolution(int scen) const {
	SimpleVector const &v = *dynamic_cast<SimpleVector const*>(dynamic_cast<StochVector const&>(*vars->y).children[scen]->vec);
	//assert(v.length());
	if(!v.length())
          return std::vector<double>(); //this vector is not on this processor
	else
	return std::vector<double>(&v[0],&v[0]+v.length());
}


template<typename FORMULATION, typename SOLVER, typename UPDATENLP>
  void NlpPIPSIpmInterface<FORMULATION,SOLVER,UPDATENLP>::writeSolution() const {
  updateNlpInfo->writeSolution(vars);
}
#endif

