/* PIPS-IPM                                                           *
 * Author:  Cosmin G. Petra                                           *
 * (C) 2012 Argonne National Laboratory. See Copyright Notification.  */

#ifndef PIPSIPM_INTERFACE
#define PIPSIPM_INTERFACE

#include <algorithm>
#include <functional>

#include "stochasticInput.hpp"

#include "sTree.h"
#include "sData.h"
#include "sResiduals.h"
#include "sVars.h"
#include "StochMonitor.h"
#include <cstdlib>
#include <stdexcept>
#include <algorithm>


#include "PreprocessFactory.h"
#include "Scaler.h"
#include "Presolver.h"
#include "Postsolver.h"

#include "sTreeCallbacks.h"
#include "pipsport.h"

//#define PRESOLVE_POSTSOLVE_ONLY // will not call solve routine an just presolve and then postsolve the problem - for debugging presolve

template<class FORMULATION, class IPMSOLVER> 
class PIPSIpmInterface 
{
 public:
  PIPSIpmInterface(stochasticInput &in, MPI_Comm = MPI_COMM_WORLD);
  PIPSIpmInterface(StochInputTree* in, MPI_Comm = MPI_COMM_WORLD,
        ScalerType scaler_type = SCALER_NONE, PresolverType presolver_type = PRESOLVER_NONE);
  ~PIPSIpmInterface();

  void go();
  double getObjective();
  double getFirstStageObjective() const;

  void setPrimalTolerance(double val);
  void setDualTolerance(double val);

  std::vector<double> gatherPrimalSolution();
  std::vector<double> gatherDualSolutionEq();
  std::vector<double> gatherDualSolutionIneq();
  std::vector<double> gatherDualSolutionIneqUpp();
  std::vector<double> gatherDualSolutionIneqLow();
  std::vector<double> gatherDualSolutionVarBounds();
  std::vector<double> gatherDualSolutionVarBoundsUpp();
  std::vector<double> gatherDualSolutionVarBoundsLow();

  std::vector<double> getFirstStagePrimalColSolution() const;
  std::vector<double> getSecondStagePrimalColSolution(int scen) const;
  //std::vector<double> getFirstStageDualColSolution() const{};
  std::vector<double> getFirstStageDualRowSolution() const;
  //std::vector<double> getSecondStageDualColSolution(int scen) const{};
  std::vector<double> getSecondStageDualRowSolution(int scen) const;

  void postsolveComputedSolution();

  std::vector<double> gatherEqualityConsValues();
  std::vector<double> gatherInequalityConsValues();

  void getVarsUnscaledUnperm();
  void getResidsUnscaledUnperm();
  //more get methods to follow here

  static bool isDistributed() { return true; }

 protected:

  FORMULATION * factory;
  PreprocessFactory * prefactory;
  sData *        data;       // possibly presolved data
  sData *        origData;   // original data
  sVars *        vars;
  sVars *        unscaleUnpermVars;
  sVars *        postsolvedVars;

  sResiduals *   resids;
  sResiduals *   unscaleUnpermResids;
  sResiduals *   postsolvedResids;

  Presolver*    presolver;
  Postsolver* postsolver;
  Scaler *      scaler;
  IPMSOLVER *   solver;

  PIPSIpmInterface() {};
  MPI_Comm comm;
  bool ran_solver;
};

//----------------------------------------------------------------------
// IMPLEMENTATION
//----------------------------------------------------------------------


template<class FORMULATION, class IPMSOLVER>
PIPSIpmInterface<FORMULATION, IPMSOLVER>::PIPSIpmInterface(stochasticInput &in, MPI_Comm comm) :  unscaleUnpermVars(nullptr), postsolvedVars(nullptr), 
  unscaleUnpermResids(nullptr), postsolvedResids(nullptr), comm(comm), ran_solver(false)
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

  scaler = nullptr;

  solver  = new IPMSOLVER( factory, data );
  solver->addMonitor(new StochMonitor( factory ));
#ifdef TIMING
  if(mype==0) printf("solver created\n");
  //solver->monitorSelf();
#endif

}

template<class FORMULATION, class IPMSOLVER>
PIPSIpmInterface<FORMULATION, IPMSOLVER>::PIPSIpmInterface(StochInputTree* in, MPI_Comm comm, ScalerType scaler_type,
      PresolverType presolver_type) : unscaleUnpermVars(nullptr), postsolvedVars(nullptr), unscaleUnpermResids(nullptr), postsolvedResids(nullptr), comm(comm), ran_solver(false)
{
  bool postsolve = true; // todo

  int mype;
  MPI_Comm_rank(comm,&mype);

  MPI_Barrier(comm);
  const double t0 = MPI_Wtime();

  factory = new FORMULATION( in, comm);
#ifdef TIMING
  if(mype==0) printf("factory created\n");
#endif

  prefactory = new PreprocessFactory();

  // presolving activated?
  if( presolver_type != PRESOLVER_NONE )
  {

     origData = dynamic_cast<sData*>(factory->makeData());

     MPI_Barrier(comm);
     const double t0_presolve = MPI_Wtime();

     postsolver = (postsolve == true) ? prefactory->makePostsolver(origData) : nullptr;
     presolver = prefactory->makePresolver(origData, presolver_type, postsolver);

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
     origData = nullptr;
     postsolver = nullptr;
     presolver = nullptr;
  }

#if 0
  ofstream myfile;
  myfile.open ("PipsToMPS_prslv.mps");
  data->writeMPSformat(myfile);
  myfile.close();
#endif

#ifdef TIMING
  if(mype==0) printf("data created\n");
#endif

#ifdef WITH_PARDISOINDEF
  data->activateLinkStructureExploitation();
#endif

  vars   = dynamic_cast<sVars*>( factory->makeVariables( data ) );
#ifdef TIMING
  if(mype==0) printf("variables created\n");
#endif

  resids = dynamic_cast<sResiduals*>( factory->makeResiduals( data ) );
#ifdef TIMING
  if(mype==0) printf("resids created\n");
#endif

  scaler = prefactory->makeScaler(data, scaler_type);

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

#if defined(PRESOLVE_POSTSOLVE_ONLY) && !defined(NDEBUG)
  const int result = 0;
#else
  //---------------------------------------------
  const int result = solver->solve(data,vars,resids);
  //---------------------------------------------
#endif

  if( result != 0 && mype == 0 )
     std::cout << "failed to solve instance, result code: " << result << std::endl;
  
  ran_solver = true;

#ifdef TIMING
   if ( 0 != result )
      return;

   tmElapsed = MPI_Wtime()-tmElapsed;

   const double objective = getObjective();

   if( 0 == mype ) {
    //cout << " " << data->nx << " variables, " << data->my
    // << " equality constraints, " << data->mz << " inequality constraints.\n";

    cout << " Iterates: " << solver->iter <<",    Optimal Solution:  "
	 << objective << endl;

    cout << "Solve time: " << tmElapsed << " seconds." << endl;

    char *var = getenv("OMP_NUM_THREADS");
    if(var != nullptr) {
      int num_threads;
      sscanf( var, "%d", &num_threads );
      cout << "Num threads: " << num_threads << endl;
    }
  }
#endif

#if !defined(NDEBUG) && defined(PRESOLVE_POSTSOLVE_ONLY)
  postsolveComputedSolution();
#endif
}

template<typename FORMULATION, typename SOLVER>
double PIPSIpmInterface<FORMULATION,SOLVER>::getObjective() {

  if(!ran_solver)
    throw std::logic_error("Must call go() and start solution process before trying to retrieve original solution");

  if( postsolver != nullptr && postsolvedVars == nullptr)
    this->postsolveComputedSolution();

  double obj;
  if(postsolvedVars != nullptr)
    obj = origData->objectiveValue(postsolvedVars);
  else
  {
    obj = data->objectiveValue(vars);
    if( scaler )
       obj = scaler->getObjUnscaled(obj);
  }

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
  delete scaler;
  delete postsolvedResids;
  delete unscaleUnpermResids;
  delete resids;
  delete postsolvedVars;
  delete unscaleUnpermVars;
  delete vars;
  delete data;
  delete postsolver;
  delete presolver;
  delete origData;
  delete prefactory;
  delete factory;
}

template<class FORMULATION, class IPMSOLVER>
void PIPSIpmInterface<FORMULATION, IPMSOLVER>::getVarsUnscaledUnperm()
{
  assert(unscaleUnpermVars == nullptr);
  if(!ran_solver)
    throw std::logic_error("Must call go() and start solution process before trying to retrieve unscaled unpermutated solution");
  if( scaler )
  {
    sVars* unscaled_vars = dynamic_cast<sVars*>(scaler->getVariablesUnscaled(*vars));
    unscaleUnpermVars = data->getVarsUnperm(*unscaled_vars);
    delete unscaled_vars;
  }
  else
    unscaleUnpermVars = data->getVarsUnperm(*vars);

}

template<class FORMULATION, class IPMSOLVER>
void PIPSIpmInterface<FORMULATION, IPMSOLVER>::getResidsUnscaledUnperm()
{
  assert(unscaleUnpermResids == nullptr);

  if(!ran_solver)
    throw std::logic_error("Must call go() and start solution process before trying to retrieve unscaled unpermutated residuals");

  if( scaler )
  {
    sResiduals* unscaled_resids = dynamic_cast<sResiduals*>(scaler->getResidualsUnscaled(*resids));
    unscaleUnpermResids = data->getResidsUnperm(*unscaled_resids);
    delete unscaled_resids;
  }
  else
    unscaleUnpermResids = data->getResidsUnperm(*resids);
}


template<class FORMULATION, class IPMSOLVER>
std::vector<double> PIPSIpmInterface<FORMULATION, IPMSOLVER>::gatherPrimalSolution()
{
  if( unscaleUnpermVars == nullptr)
    this->getVarsUnscaledUnperm();

  if( postsolver != nullptr && postsolvedVars == nullptr)
    this->postsolveComputedSolution();

  std::vector<double> vec;
  if( postsolver == nullptr)
    vec = dynamic_cast<const StochVector&>(*unscaleUnpermVars->x).gatherStochVector();
  else
    vec = dynamic_cast<const StochVector&>(*postsolvedVars->x).gatherStochVector();

  return vec;
}

template<class FORMULATION, class IPMSOLVER>
std::vector<double> PIPSIpmInterface<FORMULATION, IPMSOLVER>::gatherDualSolutionEq()
{
  if( unscaleUnpermVars == nullptr)
    this->getVarsUnscaledUnperm();

  if( postsolver != nullptr && postsolvedVars == nullptr )
    this->postsolveComputedSolution();

  std::vector<double> vec;
  if(postsolver == nullptr)
    vec = dynamic_cast<const StochVector&>(*unscaleUnpermVars->y).gatherStochVector();
  else
    vec = dynamic_cast<const StochVector&>(*postsolvedVars->y).gatherStochVector();

  return vec;
}

template<class FORMULATION, class IPMSOLVER>
std::vector<double> PIPSIpmInterface<FORMULATION, IPMSOLVER>::gatherDualSolutionIneq()
{
  if( unscaleUnpermVars == nullptr)
    this->getVarsUnscaledUnperm();

  if( postsolver != nullptr && postsolvedVars == nullptr )
    this->postsolveComputedSolution();
  
  std::vector<double> vec;
  if( postsolver == nullptr)
    vec = dynamic_cast<const StochVector&>(*unscaleUnpermVars->z).gatherStochVector();
  else
    vec = dynamic_cast<const StochVector&>(*postsolvedVars->z).gatherStochVector();

  return vec;
}

template<class FORMULATION, class IPMSOLVER>
std::vector<double> PIPSIpmInterface<FORMULATION, IPMSOLVER>::gatherDualSolutionIneqUpp()
{
  if( unscaleUnpermVars == nullptr)
    this->getVarsUnscaledUnperm();

  if( postsolver != nullptr && postsolvedVars == nullptr )
    this->postsolveComputedSolution();

  std::vector<double> vec;
  if( postsolver == nullptr)
    vec = dynamic_cast<const StochVector&>(*unscaleUnpermVars->pi).gatherStochVector();
  else
    vec = dynamic_cast<const StochVector&>(*postsolvedVars->pi).gatherStochVector();

  return vec;
}

template<class FORMULATION, class IPMSOLVER>
std::vector<double> PIPSIpmInterface<FORMULATION, IPMSOLVER>::gatherDualSolutionIneqLow()
{
  if( unscaleUnpermVars == nullptr)
    this->getVarsUnscaledUnperm();

  if( postsolver != nullptr && postsolvedVars == nullptr )
    this->postsolveComputedSolution();
  
  std::vector<double> vec;
  if( postsolver == nullptr)
    vec = dynamic_cast<const StochVector&>(*unscaleUnpermVars->lambda).gatherStochVector();
  else
    vec = dynamic_cast<const StochVector&>(*postsolvedVars->lambda).gatherStochVector();

  return vec;
}

template<class FORMULATION, class IPMSOLVER>
std::vector<double> PIPSIpmInterface<FORMULATION, IPMSOLVER>::gatherDualSolutionVarBounds()
{
  std::vector<double> duals_varbounds_upp = gatherDualSolutionVarBoundsUpp();
  std::vector<double> duals_varbounds_low = gatherDualSolutionVarBoundsLow();

  assert(duals_varbounds_low.size() == duals_varbounds_upp.size());

  std::vector<double> duals_varbounds;
  duals_varbounds.reserve(duals_varbounds_low.size());

  std::transform(duals_varbounds_low.begin(), duals_varbounds_low.end(), duals_varbounds_upp.begin(), std::back_inserter(duals_varbounds), std::minus<double>());

  return duals_varbounds;
}

template<class FORMULATION, class IPMSOLVER>
std::vector<double> PIPSIpmInterface<FORMULATION, IPMSOLVER>::gatherDualSolutionVarBoundsUpp()
{
  if( unscaleUnpermVars == nullptr)
    this->getVarsUnscaledUnperm();

  if( postsolver != nullptr && postsolvedVars == nullptr )
    this->postsolveComputedSolution();

  std::vector<double> vec;
  if( postsolver == nullptr)
    vec = dynamic_cast<const StochVector&>(*unscaleUnpermVars->phi).gatherStochVector();
  else
    vec = dynamic_cast<const StochVector&>(*postsolvedVars->phi).gatherStochVector();

  return vec;
}

template<class FORMULATION, class IPMSOLVER>
std::vector<double> PIPSIpmInterface<FORMULATION, IPMSOLVER>::gatherDualSolutionVarBoundsLow()
{
  if( unscaleUnpermVars == nullptr)
    this->getVarsUnscaledUnperm();

  if( postsolver != nullptr && postsolvedVars == nullptr )
    this->postsolveComputedSolution();

  std::vector<double> vec;
  if( postsolver == nullptr)
    vec = dynamic_cast<const StochVector&>(*unscaleUnpermVars->gamma).gatherStochVector();
  else
    vec = dynamic_cast<const StochVector&>(*postsolvedVars->gamma).gatherStochVector();

  return vec;
}

template<class FORMULATION, class IPMSOLVER>
std::vector<double> PIPSIpmInterface<FORMULATION, IPMSOLVER>::gatherEqualityConsValues()
{
  if( unscaleUnpermResids == nullptr)
    this->getResidsUnscaledUnperm();

  if( postsolver != nullptr && postsolvedVars == nullptr )
    this->postsolveComputedSolution();

  StochVector* eq_vals = (postsolvedVars == nullptr) ? dynamic_cast<StochVector*>(unscaleUnpermResids->rA->cloneFull()) : 
    dynamic_cast<StochVector*>(postsolvedResids->rA->cloneFull());

  if( origData == nullptr || postsolvedVars == nullptr )
    eq_vals->axpy(1.0, *data->bA);
  else
    eq_vals->axpy(1.0, *origData->bA);

  std::vector<double> eq_vals_vec = eq_vals->gatherStochVector();

  delete eq_vals;

  return eq_vals_vec;
}


template<class FORMULATION, class IPMSOLVER>
std::vector<double> PIPSIpmInterface<FORMULATION, IPMSOLVER>::gatherInequalityConsValues()
{
  if( unscaleUnpermVars == nullptr)
    this->getVarsUnscaledUnperm();

  if( unscaleUnpermResids == nullptr)
    this->getResidsUnscaledUnperm();

  if( postsolver != nullptr && postsolvedVars == nullptr )
    this->postsolveComputedSolution();

  StochVector* ineq_vals = (postsolvedVars == nullptr) ? dynamic_cast<StochVector*>(unscaleUnpermResids->rC->cloneFull()) : 
    dynamic_cast<StochVector*>(postsolvedResids->rC->cloneFull());

  if( postsolvedVars == nullptr )
    ineq_vals->axpy(1.0, *unscaleUnpermVars->s);
  else
    ineq_vals->axpy(1.0, *postsolvedVars->s);

  std::vector<double> ineq_vals_vec = ineq_vals->gatherStochVector();

  delete ineq_vals;

  return ineq_vals_vec;
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
	  return std::vector<double>( &v[0], &v[0] + v.length() );
}

template<class FORMULATION, class IPMSOLVER>
std::vector<double> PIPSIpmInterface<FORMULATION, IPMSOLVER>::getFirstStageDualRowSolution() const
{
  SimpleVector const &y =
        *dynamic_cast<SimpleVector const*>((dynamic_cast<StochVector const&>(*vars->y)).vec);
  SimpleVector const &z =
        *dynamic_cast<SimpleVector const*>((dynamic_cast<StochVector const&>(*vars->z)).vec);

  if( !y.length() && !z.length() )
     return std::vector<double>(); //this vector is not on this processor
  else
  {
     std::vector<int> const &map = factory->tree->idx_EqIneq_Map;

     std::vector<double> multipliers(map.size());
     for( size_t i = 0; i < map.size(); i++ )
     {
        int idx = map[i];
        if( idx < 0 )
        {
           //equality
           idx = -idx - 1;
           assert(idx >= 0);
           multipliers[i] = y[idx];
        }
        else
        {
           //inequality - since, we have z-\lambda+\pi=0, where \lambda is the multiplier for low and
           //\pi is the multiplier for upp, therefore z containts the right multiplier for this row.
#ifndef NDEBUG
           SimpleVector const &iclow = *dynamic_cast<SimpleVector const*>((dynamic_cast<StochVector const&>(*vars->iclow)).vec);
           SimpleVector const &icupp = *dynamic_cast<SimpleVector const*>((dynamic_cast<StochVector const&>(*vars->icupp)).vec);
           assert(iclow[idx] > 0 || icupp[idx] > 0);
#endif
           multipliers[i] = z[idx];
        }
     }
     return multipliers;
  }
}


template<class FORMULATION, class IPMSOLVER>
std::vector<double> PIPSIpmInterface<FORMULATION, IPMSOLVER>::getSecondStageDualRowSolution(int scen) const 
{
  SimpleVector const &y = *dynamic_cast<SimpleVector const*>(dynamic_cast<StochVector const&>(*vars->y).children[scen]->vec);
  SimpleVector const &z = *dynamic_cast<SimpleVector const*>(dynamic_cast<StochVector const&>(*vars->z).children[scen]->vec);
  SimpleVector const &iclow = *dynamic_cast<SimpleVector const*>(dynamic_cast<StochVector const&>(*vars->iclow).children[scen]->vec);
  SimpleVector const &icupp = *dynamic_cast<SimpleVector const*>(dynamic_cast<StochVector const&>(*vars->icupp).children[scen]->vec);
  //assert(v.length());
  if( !y.length() && !z.length() )
    return std::vector<double>(); //this vector is not on this processor
  else 
  {
    std::vector<int> const &map=factory->tree->children[scen]->idx_EqIneq_Map;

    std::vector<double> multipliers(map.size());
    for(size_t i = 0; i < map.size(); i++)
    {
      int idx = map[i];
      if(idx<0) 
      {
	      //equality
	      idx =- idx - 1;
        assert( idx >= 0 );
        multipliers[i] = y[idx];
      } 
      else
      {
	      //inequality - since, we have z-\lambda+\pi=0, where \lambda is the multiplier for low and
	      //\pi is the multiplier for upp, therefore z containts the right multiplier for this row.
#ifndef NDEBUG
	      assert(iclow[idx] > 0 || icupp[idx] > 0);
	      multipliers[i] = z[idx];
#endif
      }
    }
    return multipliers;
  }
  //return std::vector<double>(&v[0],&v[0]+v.length());
}

template<class FORMULATION, class IPMSOLVER>
void PIPSIpmInterface<FORMULATION, IPMSOLVER>::postsolveComputedSolution()
{
  int my_rank;
  MPI_Comm_rank(comm,&my_rank);

  assert(origData);
  assert(data);

#if !defined(NDEBUG) && defined(PRESOLVE_POSTSOLVE_ONLY) // todo : resids for C also need recomputation.. - s variable
  resids->calcresids(data, vars);
#endif

  if( unscaleUnpermVars == nullptr)
    this->getVarsUnscaledUnperm();

  if( unscaleUnpermResids == nullptr)
    this->getResidsUnscaledUnperm();

  if( postsolvedVars != nullptr || postsolvedResids != nullptr )
    return;

  if( postsolver == nullptr )
  {
    assert( "no postsolver available" && 0 );
    return;
  }

  sTreeCallbacks& callbackTree = dynamic_cast<sTreeCallbacks&>(*origData->stochNode);
  callbackTree.switchToOriginalData();

  factory->data = origData;

  postsolvedVars = dynamic_cast<sVars*>( factory->makeVariables( origData ) );


  postsolvedResids = dynamic_cast<sResiduals*>( factory->makeResiduals( origData ) );
  postsolver->postsolve(*unscaleUnpermVars, *postsolvedVars);

  double obj_postsolved = origData->objectiveValue(postsolvedVars);
  if( my_rank == 0)
    std::cout << "Objective value after postsolve: " << obj_postsolved << std::endl;

  /* compute residuals for postprocessed solution and check for feasibility */
  postsolvedResids->calcresids(origData, postsolvedVars);
  
  const double infnorm_rA_orig = resids->rA->infnorm();
  const double infnorm_rC_orig = resids->rC->infnorm();
  const double onenorm_rA_orig = resids->rA->onenorm();
  const double onenorm_rC_orig = resids->rC->onenorm();

  const double infnorm_rA = unscaleUnpermResids->rA->infnorm();
  const double infnorm_rC = unscaleUnpermResids->rC->infnorm();
  const double onenorm_rA = unscaleUnpermResids->rA->onenorm();
  const double onenorm_rC = unscaleUnpermResids->rC->onenorm();

  const double infnorm_rA_postsolved = postsolvedResids->rA->infnorm();
  const double infnorm_rC_postsolved = postsolvedResids->rC->infnorm();
  const double onenorm_rA_postsolved = postsolvedResids->rA->onenorm();
  const double onenorm_rC_postsolved = postsolvedResids->rC->onenorm();

  if( my_rank == 0)
  {
    std::cout << "\t||.||_inf\t ||.||_1" << std::endl;
    std::cout << "Residuals of reduced problem:" << std::endl;
    std::cout << "rA:\t" << infnorm_rA_orig << "\t" << onenorm_rA_orig << std::endl;
    std::cout << "rC:\t" << infnorm_rC_orig << "\t" << onenorm_rC_orig << std::endl;
    std::cout << "Residuals after unscaling:" << std::endl;
    std::cout << "rA:\t" << infnorm_rA << "\t" << onenorm_rA << std::endl;
    std::cout << "rC:\t" << infnorm_rC << "\t" << onenorm_rC << std::endl; 
    std::cout << "Residuals after postsolve:" << std::endl;
    std::cout << "rA:\t" << infnorm_rA_postsolved << "\t" << onenorm_rA_postsolved << std::endl;
    std::cout << "rC:\t" << infnorm_rC_postsolved << "\t" << onenorm_rC_postsolved << std::endl; 
  }
}

#endif
