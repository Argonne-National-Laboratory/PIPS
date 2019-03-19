/* PIPS-IPM                                                           *
 * Author:  Cosmin G. Petra                                           *
 * (C) 2012 Argonne National Laboratory. See Copyright Notification.  */

#ifndef PIPSIPM_INTERFACE
#define PIPSIPM_INTERFACE

#include "stochasticInput.hpp"

#include "Presolver.h"
#include "sData.h"
#include "sResiduals.h"
#include "sVars.h"
#include "sTree.h"
#include "StochMonitor.h"
#include "Scaler.h"
#include <cstdlib>
#include "PreprocessFactory.h"

template<class FORMULATION, class IPMSOLVER> 
class PIPSIpmInterface 
{
 public:
  PIPSIpmInterface(stochasticInput &in, MPI_Comm = MPI_COMM_WORLD);
  PIPSIpmInterface(StochInputTree* in, MPI_Comm = MPI_COMM_WORLD,
        ScalerType scaler_type = SCALER_NONE, PresolverType presolver_type = PRESOLVER_NONE);
  ~PIPSIpmInterface();

  void go();
  double getObjective() const;
  double getFirstStageObjective() const;


  void setPrimalTolerance(double val);
  void setDualTolerance(double val);

  std::vector<double> gatherPrimalSolution() const;
  std::vector<double> gatherDualSolutionEq() const;
  std::vector<double> gatherDualSolutionIneq() const;
  std::vector<double> gatherDualSolutionVarBoundsUpp() const;
  std::vector<double> gatherDualSolutionVarBoundsLow() const;


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
  sData *        data;       // possibly presolved data
  sData *        origData;   // original data
  sVars *        vars;
  sResiduals *   resids;
  OoqpVector* dummyObj;

  Presolver*    presolver;
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

  scaler = NULL;

  solver  = new IPMSOLVER( factory, data );
  solver->addMonitor(new StochMonitor( factory ));
#ifdef TIMING
  if(mype==0) printf("solver created\n");
  //solver->monitorSelf();
#endif

}

template<class FORMULATION, class IPMSOLVER>
PIPSIpmInterface<FORMULATION, IPMSOLVER>::PIPSIpmInterface(StochInputTree* in, MPI_Comm comm, ScalerType scaler_type,
      PresolverType presolver_type) : comm(comm)
{

  int mype;
  MPI_Comm_rank(comm,&mype);

  MPI_Barrier(comm);
  const double t0 = MPI_Wtime();

  factory = new FORMULATION( in, comm);
#ifdef TIMING
  if(mype==0) printf("factory created\n");
#endif

  const PreprocessFactory& prefactory = PreprocessFactory::getInstance();

  // presolving activated?
  if( presolver_type != PRESOLVER_NONE )
  {

     origData = dynamic_cast<sData*>(factory->makeData());

/*     ofstream myfile;
     myfile.open ("PipsToMPS_original.mps");
     origData->writeMPSformat(myfile);
     myfile.close();*/

     MPI_Barrier(comm);
     const double t0_presolve = MPI_Wtime();

     presolver = prefactory.makePresolver(origData, presolver_type);

     data = dynamic_cast<sData*>(presolver->presolve());

     factory->data = data; // todo update also sTree* of factory

     MPI_Barrier(comm);
     const double t_presolve = MPI_Wtime();
     if( mype == 0 )
        std::cout << "---presolve time (in sec.): " << t_presolve - t0_presolve << std::endl;
  }
  else
  {
     data = dynamic_cast<sData*>(factory->makeData());
     origData = NULL;
     presolver = NULL;
  }
//dynamic_cast<StochVector*>
  dummyObj = (data->g->cloneFull());

/*  ofstream myfile;
  myfile.open ("PipsToMPS_prslv.mps");
  data->writeMPSformat(myfile);
  myfile.close();*/

#ifdef TIMING
  if(mype==0) printf("data created\n");
#endif

#ifdef WITH_PARDISOINDEF
  data->activateLinkStructureExploitation();
#endif

  vars   = dynamic_cast<sVars*>     ( factory->makeVariables( data ) );
#ifdef TIMING
  if(mype==0) printf("variables created\n");
#endif

  resids = dynamic_cast<sResiduals*>( factory->makeResiduals( data ) );
#ifdef TIMING
  if(mype==0) printf("resids created\n");
#endif

  scaler = prefactory.makeScaler(data, scaler_type);

#ifdef TIMING
  if(mype==0) printf("scaler created\n");
#endif

  solver  = new IPMSOLVER( factory, data );
  solver->addMonitor(new StochMonitor( factory, scaler ));
#ifdef TIMING
  if(mype==0) printf("solver created\n");
  //solver->monitorSelf();
#endif

  MPI_Barrier(comm);
  const double t1 = MPI_Wtime();

  if( mype == 0 )
     std::cout << "---reading time (in sec.): " << t1 - t0 << std::endl;
}


template<typename FORMULATION, typename IPMSOLVER>
void PIPSIpmInterface<FORMULATION,IPMSOLVER>::go() {

   int mype;
   MPI_Comm_rank(comm,&mype);

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
#ifdef TIMING
  double tmElapsed=MPI_Wtime();
#endif

  if( scaler )
  {
     MPI_Barrier(comm);
     const double t0_scaling = MPI_Wtime();

     scaler->scale();

     MPI_Barrier(comm);
     const double t_scaling = MPI_Wtime();
     if( mype == 0 )
        std::cout << "---scaling time (in sec.): " << t_scaling - t0_scaling << std::endl;
  }
  //---------------------------------------------
  const int result = solver->solve(data,vars,resids);
  //---------------------------------------------

  if( result != 0 && mype == 0 )
     std::cout << "failed to solve instance, result code: " << result << std::endl;

#ifdef TIMING
   if ( 0 != result )
      return;

   tmElapsed=MPI_Wtime()-tmElapsed;

   const double objective = getObjective();

   if( 0 == mype ) {
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
  }
#endif

}

template<typename FORMULATION, typename SOLVER>
double PIPSIpmInterface<FORMULATION,SOLVER>::getObjective() const {
  double obj = data->objectiveValue(vars);

  if( scaler )
     obj = scaler->getOrigObj(obj);

  return obj;
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
  delete origData;
  delete factory;
  delete scaler;
  delete presolver;
}


template<class FORMULATION, class IPMSOLVER>
std::vector<double> PIPSIpmInterface<FORMULATION, IPMSOLVER>::gatherPrimalSolution() const
{
  StochVector* primalOrg = NULL;

  if( scaler )
     primalOrg = dynamic_cast<StochVector*>(scaler->getOrigPrimal(*vars->x));
  else
     primalOrg = dynamic_cast<StochVector&>(*vars->x).cloneFull();

  const std::vector<unsigned int> permInv = data->getLinkVarsPermInv();

  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  for( size_t i = 0; i < permInv.size() && myrank == 0; i++ )
     std::cout << permInv[i] << std::endl;

  primalOrg->permuteVec0Entries(permInv);

  std::vector<double> primalVec = primalOrg->gatherStochVector();


  assert(dummyObj);


  const double obj = primalOrg->dotProductWith(*dummyObj);

  std::cout << "objXXX=" << obj << std::endl;

  delete primalOrg;


  return primalVec;
}

template<class FORMULATION, class IPMSOLVER>
std::vector<double> PIPSIpmInterface<FORMULATION, IPMSOLVER>::gatherDualSolutionEq() const
{
  StochVector* dualOrg = NULL;

  if( scaler )
     dualOrg = dynamic_cast<StochVector*>(scaler->getOrigDualEq(*vars->y));
  else
     dualOrg = dynamic_cast<const StochVector&>(*vars->y).cloneFull();

  const std::vector<unsigned int> permInv = data->getLinkConsEqPermInv();
  dualOrg->permuteLinkingEntries(permInv);

  std::vector<double> dualVec = dualOrg->gatherStochVector();

  delete dualOrg;

  return dualVec;
}

template<class FORMULATION, class IPMSOLVER>
std::vector<double> PIPSIpmInterface<FORMULATION, IPMSOLVER>::gatherDualSolutionIneq() const
{
  StochVector* dualOrg = NULL;

  if( scaler )
     dualOrg = dynamic_cast<StochVector*>(scaler->getOrigDualIneq(*vars->z));
  else
     dualOrg = dynamic_cast<const StochVector&>(*vars->z).cloneFull();

  const std::vector<unsigned int> permInv = data->getLinkConsIneqPermInv();

  dualOrg->permuteLinkingEntries(permInv);
  std::vector<double> dualVec = dualOrg->gatherStochVector();

  delete dualOrg;

  return dualVec;
}

template<class FORMULATION, class IPMSOLVER>
std::vector<double> PIPSIpmInterface<FORMULATION, IPMSOLVER>::gatherDualSolutionVarBoundsUpp() const
{
  StochVector* dualOrg = NULL;

  if( scaler )
     dualOrg = dynamic_cast<StochVector*>(scaler->getOrigDualVarBoundsUpp(*vars->phi));
  else
     dualOrg = dynamic_cast<const StochVector&>(*vars->phi).cloneFull();

  assert(!dualOrg->vecl);


  const std::vector<unsigned int> permInv = data->getLinkVarsPermInv();
  dualOrg->permuteVec0Entries(permInv);

  std::vector<double> dualVec = dualOrg->gatherStochVector();

  delete dualOrg;

  return dualVec;
}

template<class FORMULATION, class IPMSOLVER>
std::vector<double> PIPSIpmInterface<FORMULATION, IPMSOLVER>::gatherDualSolutionVarBoundsLow() const
{
  StochVector* dualOrg = NULL;

  if( scaler )
     dualOrg = dynamic_cast<StochVector*>(scaler->getOrigDualVarBoundsLow(*vars->gamma));
  else
     dualOrg = dynamic_cast<const StochVector&>(*vars->gamma).cloneFull();

  assert(!dualOrg->vecl);

  const std::vector<unsigned int> permInv = data->getLinkVarsPermInv();
  dualOrg->permuteVec0Entries(permInv);

  std::vector<double> dualVec = dualOrg->gatherStochVector();

  delete dualOrg;

  return dualVec;
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
