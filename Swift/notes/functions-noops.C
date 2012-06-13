
#include <stdlib.h>

#include "functions.h"

/*
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
*/

void
readSolution(const char *dataPath, const char *solutionPath,
             int scen, double** result, int* result_length)
{
  int length = 1;
  double* r = (double*) malloc(length*sizeof(double));
  *result = r;
  *result_length = length;
  return;
}

/*
void
customReduce(double *candidateSolution, const double *nextSolution,
             int solutionLength)
{
  return;
}
*/
