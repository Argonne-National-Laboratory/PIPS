
/**
   "f" of inner loop in loops-*.swift
*/
double* makeFeasible(const char *dataPath, int nScen, int scen,
                     const double *candidateSolution);

double evaluateSolution(const char *dataPath, int nScen,
                        const double *candidateSolution);

double* readSolution(const char *dataPath, const char *solutionPath,
                     int scen);

/**
   Reduce for inner loop. I'm assuming we reduce into the
   candidateSolution buffer. I can easily adapt to whatever format
   Swift needs for the reduce operation.
*/
void customReduce(double *candidateSolution, const double *nextSolution,
                  int solutionLength);
