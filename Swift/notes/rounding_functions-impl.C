
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <fstream>
#include <string>
#include <vector>

// ExM data helper app
// do we need extern "C" here?
#include <data.h>

#include "rawInput.hpp"
#include "ClpRecourseSolver.hpp"


// extern "C" {
#include "rounding_functions.h"
// }

/*
g++ -c -I../../Input -I../../SharedLibraries/Cbc-2.7.6/include/coin -I../../SolverInterface -I../../PIPS-S/Basic -I../../Lagrange/RecourseSubproblemSolver -I/usr/include/mpich2 rounding_functions.cpp

*/

using namespace std;

// extern "C"
double
evaluateRecourseLP(const char *dataPath, int nScen,
                   int scen, double *candidateSolution, int CS_length)
{

  vector<double> sol(candidateSolution,candidateSolution+CS_length);

  rawInput input(string(dataPath),nScen);
  int nvar1 = input.nFirstStageVars();
  assert(nvar1 == CS_length);

  ClpRecourseSolver rsol(input, scen, sol);
  rsol.setDualObjectiveLimit(1e7);

  rsol.go();

  double obj = rsol.getObjective();
  if (rsol.getStatus() == ProvenInfeasible) {
    return COIN_DBL_MAX;
  }
  assert(rsol.getStatus() == Optimal);

  double prob = input.scenarioProbability(scen);
  const vector<double> &obj1 = input.getFirstStageObj();

  for (int k = 0; k < nvar1; k++) obj += prob*sol[k]*obj1[k];

  return obj;

}

//extern "C"
struct Data*
readConvSolution(const char *dataPath, const char *solutionPath)
{
  string problemdata = string(dataPath) + "0";
  ifstream datafd(problemdata.c_str());
  int nvar1;
  datafd >> nvar1; // read size of solution

  struct Data* data = (struct Data*) malloc(sizeof(struct Data));
  data->pointer = malloc(nvar1*sizeof(double));
  data->length = nvar1*sizeof(double);

  double *vec = reinterpret_cast<double*>(data->pointer);

  ifstream solfd(solutionPath);
  for (int i = 0; i < nvar1; i++) {
    solfd >> vec[i];
  }

  return data;
}

// extern "C"
struct Data*
roundSolution(double *convSolution, int CS_length, double cutoff)
{
  struct Data* data = (struct Data*) malloc(sizeof(struct Data));
  data->pointer = malloc(CS_length*sizeof(double));
  data->length = CS_length*sizeof(double);

  double *vec = reinterpret_cast<double*>(data->pointer);

  for (int i = 0; i < CS_length; i++) {
    if (convSolution[i] >= cutoff) {
	    vec[i] = 1.;
    } else {
      vec[i] = 0.;
    }
  }

  return data;
}


