
/*
 * SWIFTPIPS.H
 *
 * Swift/T C++ leaf function headers
 */

double evaluateRecourseLP(const char *dataPath, int nScen,
				int scen, double *candidateSolution,
                                int CS_length);
struct Data* readConvSolution(const char *dataPath, const char *solutionPath);
struct Data* roundSolution(double *convSolution, int CS_length, double cutoff);
