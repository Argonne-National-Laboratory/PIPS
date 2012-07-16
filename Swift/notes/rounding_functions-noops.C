
#include <stdio.h>
#include <stdlib.h>

// ExM data helper app
#include <data.h>

#include "rounding_functions.h"

struct Data*
evaluateRecourseLP(const char *dataPath, int nScen,
                   int scen, const double *candidateSolution)
{
  printf("evaluateRecourseLP(dataPath=%s)\n", dataPath);

  int length = 10;

  struct Data* data = (struct Data*) malloc(sizeof(struct Data));
  data->pointer = NULL;
  data->length = length*sizeof(double);
  return data;
}

struct Data*
readConvSolution(const char *dataPath, const char *solutionPath)
{
  struct Data* data = (struct Data*) malloc(sizeof(struct Data));
  data->pointer = NULL;
  data->length = -1;
  return data;
}

struct Data* round(const double *convSolution, double cutoff)
{
  struct Data* data = (struct Data*) malloc(sizeof(struct Data));
  data->pointer = NULL;
  data->length = -1;
  return data;
}


