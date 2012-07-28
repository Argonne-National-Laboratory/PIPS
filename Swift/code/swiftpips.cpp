
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

  rawInput input(string(dataPath),nScen,MPI_COMM_SELF);

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

/////////////////////////////////////////////////////////////////////
// Functions that work with cached "zero" (1st stage) problem data.
////////////////////////////////////////////////////////////////////

// reads "zero" problem data and returns it in a blob
struct Data*
readStg1Data(const char *dataPath)
{
  printf("readStg1Data()...\n");
  printf("dataPath: %s\n", dataPath);
  string problemdata = string(dataPath) + "0";
  printf("problemdata file: %s\n", problemdata.c_str());
  ifstream datafd(problemdata.c_str());
  assert(datafd.good());

  string data((istreambuf_iterator<char>(datafd)),istreambuf_iterator<char>()); // read into string
  datafd.close();

  size_t datalen = data.size()+1; //include trailing \0

  struct Data* data = (struct Data*) malloc(sizeof(struct Data));
  data->pointer = malloc(datalen*sizeof(char));
  data->length = datalen*sizeof(char);

  memcpy(data->pointer, data.c_str(), datalen);

  printf("read1StgData() done.\n");
  return data;
}

// cache-friendly version of readConvSolution
struct Data*
readConvSolution2(char* zerodata, int ZD_length,
		  const char *solutionPath)
		  
{
  printf("readConvSolution2()...\n");

  istringstream f1( string(zerodata) );
  f1.exceptions(ifstream::failbit | ifstream::badbit);
  int nvar1;
  f1 >> nvar1; // read size of solution
  printf("nvar1: %i\n", nvar1);

  struct Data* data = (struct Data*) malloc(sizeof(struct Data));
  data->pointer = malloc(nvar1*sizeof(double));
  data->length = nvar1*sizeof(double);

  double *vec = reinterpret_cast<double*>(data->pointer);

  ifstream solfd(solutionPath);
  for (int i = 0; i < nvar1; i++) {
    solfd >> vec[i];
  }

  printf("readConvSolution2() done.\n");
  return data;
}

// evaluateRecourseLP2 - uses cached zero data
double
evaluateRecourseLP2(const char *dataPath, int nScen,
		    double *candidateSolution, int CS_length,
		    char* zerodata, int ZD_length)
{
  printf("evaluateRecourseLP2()...\n");
  int N = CS_length/sizeof(double);
  printf("N: %i\n", N);
  vector<double> sol(candidateSolution,candidateSolution+N);

  printf("ok0\n");

  rawInput input(string(dataPath), string(zerodata), nScen);

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

  printf("evaluateRecourseLP2() done.\n");
  return obj;

}
