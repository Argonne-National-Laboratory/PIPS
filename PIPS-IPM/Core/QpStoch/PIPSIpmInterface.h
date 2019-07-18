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

#include "PreprocessFactory.h"
#include "Scaler.h"
#include "Presolver.h"
#include "Postsolver.h"

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

  std::vector<double> gatherEqualityConsValues();
  std::vector<double> gatherInequalityConsValues();

  void getUnscaledUnpermVars();
  void getUnscaledUnpermResids();
  //more get methods to follow here

  static bool isDistributed() { return true; }

 protected:
 
  FORMULATION * factory;
  sData *        data;       // possibly presolved data
  sData *        origData;   // original data
  sVars *        vars;
  sVars *        unscaleUnpermVars;
  sResiduals *   resids;
  sResiduals *   unscaleUnpermResids;

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
PIPSIpmInterface<FORMULATION, IPMSOLVER>::PIPSIpmInterface(stochasticInput &in, MPI_Comm comm) :  unscaleUnpermVars(NULL), unscaleUnpermResids(NULL), 
  comm(comm), ran_solver(false)
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
      PresolverType presolver_type) : unscaleUnpermVars(NULL), unscaleUnpermResids(NULL), comm(comm), ran_solver(false)
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

  const PreprocessFactory& prefactory = PreprocessFactory::getInstance();

  // presolving activated?
  if( presolver_type != PRESOLVER_NONE )
  {

     origData = dynamic_cast<sData*>(factory->makeData());

     MPI_Barrier(comm);
     const double t0_presolve = MPI_Wtime();

     postsolver = (postsolve == true) ? prefactory.makePostsolver(origData) : NULL;
     presolver = prefactory.makePresolver(origData, presolver_type, postsolver);

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
     postsolver = NULL;
     presolver = NULL;
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
  
  ran_solver = true;

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

   // todo postsolve an unscaled sVars object holding the solution


}

template<typename FORMULATION, typename SOLVER>
double PIPSIpmInterface<FORMULATION,SOLVER>::getObjective() const {

  if(!ran_solver)
    throw std::logic_error("Must call go() and start solution process before trying to retrieve original solution");

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
  delete unscaleUnpermResids;
  delete resids;
  delete unscaleUnpermVars;
  delete vars;
  delete data;
  delete origData;
  delete factory;
  delete scaler;
  delete postsolver;
  delete presolver;
}

template<class FORMULATION, class IPMSOLVER>
void PIPSIpmInterface<FORMULATION, IPMSOLVER>::getUnscaledUnpermVars()
{
  assert(unscaleUnpermVars == NULL);
  if(!ran_solver)
    throw std::logic_error("Must call go() and start solution process before trying to retrieve unscaled unpermutated solution");

  sVars* unscaled_vars = dynamic_cast<sVars*>(scaler->getUnscaledVariables(*vars));
  //unscaleUnpermVars = data->getUnpermVars(*unscaled_vars);

  delete unscaled_vars;
}

template<class FORMULATION, class IPMSOLVER>
void PIPSIpmInterface<FORMULATION, IPMSOLVER>::getUnscaledUnpermResids()
{
  assert(unscaleUnpermResids== NULL);

  if(!ran_solver)
    throw std::logic_error("Must call go() and start solution process before trying to retrieve unscaled unpermutated residuals");

  sResiduals* unscaled_resids = dynamic_cast<sResiduals*>(scaler->getUnscaledResiduals(*resids));
  //unscaleUnpermResids = data->getUnpermuResiduals(unscaled_resids);

  delete unscaled_resids;
}


template<class FORMULATION, class IPMSOLVER>
std::vector<double> PIPSIpmInterface<FORMULATION, IPMSOLVER>::gatherPrimalSolution()
{
  if( unscaleUnpermVars == NULL)
    this->getUnscaledUnpermVars();

  std::vector<double> vec = dynamic_cast<const StochVector&>(*unscaleUnpermVars->x).gatherStochVector();

  return vec;
}

template<class FORMULATION, class IPMSOLVER>
std::vector<double> PIPSIpmInterface<FORMULATION, IPMSOLVER>::gatherDualSolutionEq()
{
  if( unscaleUnpermVars == NULL)
    this->getUnscaledUnpermVars();

  std::vector<double> vec = dynamic_cast<const StochVector&>(*unscaleUnpermVars->y).gatherStochVector();

  return vec;
}

template<class FORMULATION, class IPMSOLVER>
std::vector<double> PIPSIpmInterface<FORMULATION, IPMSOLVER>::gatherDualSolutionIneq()
{
  if( unscaleUnpermVars == NULL)
    this->getUnscaledUnpermVars();

  std::vector<double> vec = dynamic_cast<const StochVector&>(*unscaleUnpermVars->z).gatherStochVector();

  return vec;
}

template<class FORMULATION, class IPMSOLVER>
std::vector<double> PIPSIpmInterface<FORMULATION, IPMSOLVER>::gatherDualSolutionIneqUpp()
{
  if( unscaleUnpermVars == NULL)
    this->getUnscaledUnpermVars();

  std::vector<double> vec = dynamic_cast<const StochVector&>(*unscaleUnpermVars->pi).gatherStochVector();

  return vec;
}

template<class FORMULATION, class IPMSOLVER>
std::vector<double> PIPSIpmInterface<FORMULATION, IPMSOLVER>::gatherDualSolutionIneqLow()
{
  if( unscaleUnpermVars == NULL)
    this->getUnscaledUnpermVars();

  std::vector<double> vec = dynamic_cast<const StochVector&>(*unscaleUnpermVars->lambda).gatherStochVector();

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
  if( unscaleUnpermVars == NULL)
    this->getUnscaledUnpermVars();

  std::vector<double> vec = dynamic_cast<const StochVector&>(*unscaleUnpermVars->phi).gatherStochVector();

  return vec;
}

template<class FORMULATION, class IPMSOLVER>
std::vector<double> PIPSIpmInterface<FORMULATION, IPMSOLVER>::gatherDualSolutionVarBoundsLow()
{
  if( unscaleUnpermVars == NULL)
    this->getUnscaledUnpermVars();

  std::vector<double> vec = dynamic_cast<const StochVector&>(*unscaleUnpermVars->gamma).gatherStochVector();

  return vec;
}

template<class FORMULATION, class IPMSOLVER>
std::vector<double> PIPSIpmInterface<FORMULATION, IPMSOLVER>::gatherEqualityConsValues()
{
  if( unscaleUnpermResids == NULL)
    this->getUnscaledUnpermResids();

  StochVector* eq_vals = dynamic_cast<StochVector*>(unscaleUnpermResids->rA->cloneFull());

  eq_vals->axpy(1.0, *origData->bA);

  std::vector<double> eq_vals_vec = eq_vals->gatherStochVector();

  delete eq_vals;

  return eq_vals_vec;
}


template<class FORMULATION, class IPMSOLVER>
std::vector<double> PIPSIpmInterface<FORMULATION, IPMSOLVER>::gatherInequalityConsValues()
{
  if( unscaleUnpermVars == NULL)
    this->getUnscaledUnpermVars();

  if( unscaleUnpermResids == NULL)
    this->getUnscaledUnpermResids();

  StochVector* ineq_vals = dynamic_cast<StochVector*>(unscaleUnpermResids->rC->cloneFull());;

  ineq_vals->axpy(1.0, *unscaleUnpermVars->s);

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
	  return std::vector<double>(&v[0],&v[0]+v.length());
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
