
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <data.h>

#include "functions.h"

struct Data*
makeFeasible(const char *dataPath, int nScen, int scen,
             const double *candidateSolution)
{
  int length = 10;
  double* result = (double*) malloc(10*sizeof(double));
  for (int i = 0; i < length; i++)
    result[i] = i;

  struct Data* data = (struct Data*) malloc(sizeof(struct Data));
  data->pointer = result;
  data->length = length*sizeof(double);

  makeFeasible( ... );

  return data;
}

double
evaluateSolution(const char *dataPath, int nScen,
                 const double *candidateSolution)
{
  return 0;
}

struct Data*
readSolution(const char *dataPath, const char *solutionPath,
             int scen)
{
  printf("readSolution: %s %s %i\n", dataPath, solutionPath, scen);

  int length = 4;
  double* result = (double*) malloc(length*sizeof(double));
  for (int i = 0; i < length; i++)
    result[i] = i;

  struct Data* data = (struct Data*) malloc(sizeof(struct Data));
  data->pointer = result;
  data->length = length*sizeof(double);
  return data;
}

void
customReduce(double *candidateSolution, const double *nextSolution,
             int solutionLength)
{
  return;
}
