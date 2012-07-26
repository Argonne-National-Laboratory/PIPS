
/*
 * SWIFTPIPS.CPP
 *
 * Swift/T C++ leaf functions
 */

#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <fstream>
#include <string>
#include <vector>

// ExM swig-data module:
#include <data.h>

#include "rawInput.hpp"
#include "ClpRecourseSolver.hpp"

#include "swiftpips.h"

using namespace std;

double
evaluateRecourseLP(const char *dataPath, int nScen,
                   int scen, double *candidateSolution, int CS_length)
{
  printf("evaluateRecourseLP()...\n");
  int N = CS_length/sizeof(double);
  printf("N: %i\n", N);
  vector<double> sol(candidateSolution,candidateSolution+N);

  printf("ok0\n");

  rawInput input(string(dataPath),nScen);

  printf("ok.5\n");

  int nvar1 = input.nFirstStageVars();
  printf("nvar1: %i\n", nvar1);
  assert(nvar1 == N);

  printf("ok1\n");

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

  printf("evaluateRecourseLP() done.\n");
  return obj;

}

struct Data*
readConvSolution(const char *dataPath, const char *solutionPath)
{
  printf("readConvSolution()...\n");
  printf("dataPath: %s\n", dataPath);
  string problemdata = string(dataPath) + "0";
  printf("problemdata: %s\n", problemdata.c_str());
  ifstream datafd(problemdata.c_str());
  assert(datafd.good());
  int nvar1;
  datafd >> nvar1; // read size of solution
  printf("nvar1: %i\n", nvar1);

  struct Data* data = (struct Data*) malloc(sizeof(struct Data));
  data->pointer = malloc(nvar1*sizeof(double));
  data->length = nvar1*sizeof(double);

  double *vec = reinterpret_cast<double*>(data->pointer);

  ifstream solfd(solutionPath);
  for (int i = 0; i < nvar1; i++) {
    solfd >> vec[i];
  }

  printf("readConvSolution() done.\n");
  return data;
}

struct Data*
roundSolution(double *convSolution, int CS_length, double cutoff)
{
  printf("roundSolution...\n");
  printf("convSolution: %p\n", convSolution);
  struct Data* data = (struct Data*) malloc(sizeof(struct Data));
  data->pointer = malloc(CS_length*sizeof(double));
  data->length = CS_length;

  double *vec = reinterpret_cast<double*>(data->pointer);
  int N = CS_length / sizeof(double);

  for (int i = 0; i < N; i++) {
    if (convSolution[i] >= cutoff) {
      vec[i] = 1.;
    } else {
      vec[i] = 0.;
    }
  }

  printf("roundSolution() done.\n");
  return data;
}
