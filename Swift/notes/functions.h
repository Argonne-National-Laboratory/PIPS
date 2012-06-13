
/*
   Lubin: The variables dataPath, nScen, and solutionPath are
   effectively compile-time (or more nicely, command-line)
   parameters.

   Wozniak: Let's keep these in the API so you can control them easily
   from Swift
*/

/**
   "f" of inner loop in loops-*.swift
*/
void makeFeasible(const char *dataPath, int nScen, int scen,
                  const double *candidateSolution,
                  double** result, int* result_length);

double evaluateSolution(const char *dataPath, int nScen,
                        const double *candidateSolution);

void readSolution(const char *dataPath, const char *solutionPath,
                  int scen, double** result, int* result_length);

/**
   Reduce for inner loop. I'm assuming we reduce into the
   candidateSolution buffer. I can easily adapt to whatever format
   Swift needs for the reduce operation.
*/
void customReduce(double *candidateSolution, const double *nextSolution,
                  int solutionLength);
