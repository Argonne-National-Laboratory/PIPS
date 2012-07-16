// the following are for the "rounding" app:
struct Data* evaluateRecourseLP(const char *dataPath, int nScen,
				int scen, double *candidateSolution,
                                int CS_length);
struct Data* readConvSolution(const char *dataPath, const char *solutionPath);
struct Data* round(double *convSolution, int CS_length, double cutoff);
