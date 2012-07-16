
#ifndef ROUNDING_FUNCTIONS_SWIFT
#define ROUNDING_FUNCTIONS_SWIFT

(blob b) evaluateRecourseLP(string dataPath, int nScen,
                            int scen, blob candidateSolution)
"rounding_functions" "0.0" "evaluateRecourseLP_turbine";

(blob b) readConvSolution(string dataPath, string solutionPath)
"rounding_functions" "0.0" "readConvSolution_turbine";

(blob b) round(blob convSolution, float cutoff)
"rounding_functions" "0.0" "round_turbine";

#endif
