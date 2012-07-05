// the following are for the "rounding" app:
struct Data* evaluateRecourseLP(const char *dataPath, int nScen,
				int scen, const double *candidateSolution);
struct Data* readConvSolution(const char *dataPath, const char *solutionPath);
struct Data* round(const double *convSolution, double cutoff);
