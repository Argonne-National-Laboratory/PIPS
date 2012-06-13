
#include <stdlib.h>

#include "functions.h"

double*
makeFeasible(const char *dataPath, int nScen, int scen,
             const double *candidateSolution)
{
  double* result = (double*) malloc(1*sizeof(double));
  return result;
}


double
evaluateSolution(const char *dataPath, int nScen,
                 const double *candidateSolution)
{
  return 0;
}

double*
readSolution(const char *dataPath, const char *solutionPath,
             int scen)
{
  double* result = (double*) malloc(1*sizeof(double));
  return result;
}

void
customReduce(double *candidateSolution, const double *nextSolution,
             int solutionLength)
{
  return;
}
