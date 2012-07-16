
#include <stdio.h>
#include <stdlib.h>

// ExM data helper app
#include <data.h>

#include "rounding_functions.h"

struct Data*
evaluateRecourseLP(const char *dataPath, int nScen,
                   int scen, double *candidateSolution, int CS_length)
{
  printf("evaluateRecourseLP(dataPath=%s)\n", dataPath);

  int length = 10;

  struct Data* data = (struct Data*) malloc(sizeof(struct Data));
  data->pointer = calloc(length, sizeof(double));
  data->length = length * sizeof(double);

  return data;
}

struct Data*
readConvSolution(const char *dataPath, const char *solutionPath)
{
  int length = 10;
  struct Data* data = (struct Data*) malloc(sizeof(struct Data));
  data->pointer = calloc(length, sizeof(double));
  data->length = length*sizeof(double);
  return data;
}

struct Data*
round(double *convSolution, int CS_length, double cutoff)
{
  int length = 10;
  struct Data* data = (struct Data*) malloc(sizeof(struct Data));
  data->pointer = calloc(length, sizeof(double));
  data->length = length*sizeof(double);
  return data;
}


