/* PIPS-IPM                                                           *
 * Author:  Cosmin G. Petra                                           *
 * (C) 2012 Argonne National Laboratory. See Copyright Notification.  */

#ifndef PIPSIPM_INTERFACE
#define PIPSIPM_INTERFACE

#include "stochasticInput.hpp"

#include "sTree.h"
#include "sData.h"
#include "sResiduals.h"
#include "sVars.h"
#include "StochMonitor.h"
#include <cstdlib>

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

  std::vector<double> gatherPrimalSolution() const;

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
  PreprocessFactory * prefactory;
  sData *        data;       // possibly presolved data
  sData *        origData;   // original data
  sVars *        vars;
  sResiduals *   resids;

  Presolver*    presolver;
  Postsolver* postsolver;
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

/*     ofstream myfile;
     myfile.open ("PipsToMPS_original.mps");
     origData->writeMPSformat(myfile);
     myfile.close();*/

     MPI_Barrier(comm);
     const double t0_presolve = MPI_Wtime();

     postsolver = (postsolve == true) ? prefactory->makePostsolver(origData) : NULL;
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
     origData = NULL;
     postsolver = NULL;
     presolver = NULL;
  }

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

   // todo postsolve an unscaled sVars object holding the solution


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
  delete scaler;
  delete resids;
  delete vars;
  delete data;
  delete postsolver;
  delete presolver;
  delete origData;
  delete prefactory;
  delete factory;
}


template<class FORMULATION, class IPMSOLVER>
std::vector<double> PIPSIpmInterface<FORMULATION, IPMSOLVER>::gatherPrimalSolution() const
{
  StochVector* primalUnscaled = NULL;

  if( scaler )
     primalUnscaled = dynamic_cast<StochVector*>(scaler->getOrigPrimal(*vars->x));

  const StochVector& primalStochVec = (scaler) ? *primalUnscaled : dynamic_cast<const StochVector&>(*vars->x);

  const SimpleVector& firstvec =  dynamic_cast<const SimpleVector&>(*primalStochVec.vec);
  const size_t nChildren = primalStochVec.children.size();

  int myrank; MPI_Comm_rank(comm, &myrank);
  int mysize; MPI_Comm_size(comm, &mysize);

  std::vector<double> primalVecLocal;

  for( size_t i = 0; i < nChildren; ++i )
  {
     const SimpleVector& vec = dynamic_cast<const SimpleVector&>(*primalStochVec.children[i]->vec);

     if( vec.length() > 0 )
        primalVecLocal.insert(primalVecLocal.end(), &vec[0], &vec[0] + vec.length());
  }

  size_t solLength = firstvec.length();

  // final vector
  std::vector<double> primalVec(0);

  if( mysize > 0 )
  {
     // get all lengths
     std::vector<int> recvcounts(mysize);
     std::vector<int> recvoffsets(mysize);

     int mylength = int(primalVecLocal.size());

     MPI_Allgather(&mylength, 1, MPI_INT, &recvcounts[0], 1, MPI_INT, comm);

     // all-gather local components
     recvoffsets[0] = 0;
     for( size_t i = 1; i < size_t(mysize); ++i )
        recvoffsets[i] = recvoffsets[i - 1] + recvcounts[i - 1];

     if( myrank == 0 )
     {
        solLength += recvoffsets[mysize - 1] + recvcounts[mysize - 1];
        primalVec = std::vector<double>(solLength);

        MPI_Gatherv(&primalVecLocal[0], mylength, MPI_DOUBLE, &primalVec[0] + firstvec.length(),
              &recvcounts[0], &recvoffsets[0], MPI_DOUBLE, 0, comm);
     }
     else
     {
        MPI_Gatherv(&primalVecLocal[0], mylength, MPI_DOUBLE, 0, &recvcounts[0], &recvoffsets[0], MPI_DOUBLE, 0, comm);
     }
  }
  else
  {
     solLength += primalVecLocal.size();

     primalVec = std::vector<double>(solLength);

     std::copy(primalVecLocal.begin(), primalVecLocal.end(), primalVec.begin() + firstvec.length());
  }

  if( myrank == 0 )
     std::copy(&firstvec[0], &firstvec[0] + firstvec.length(), &primalVec[0]);

  delete primalUnscaled;

  return primalVec;
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
