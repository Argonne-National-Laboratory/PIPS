
/*
 * ROUNDING.SWIFT
 *
 * Swift/T header file
 */

#ifndef ROUNDING_FUNCTIONS_SWIFT
#define ROUNDING_FUNCTIONS_SWIFT

(float result) evaluateRecourseLP(string dataPath, int nScen,
                                   int scen, blob candidateSolution)
"swiftpips" "0.0" "evaluateRecourseLP_turbine";

(blob b) readConvSolution(string dataPath, string solutionPath)
"swiftpips" "0.0" "readConvSolution_turbine";

(blob b) roundSolution(blob convSolution, float cutoff)
"swiftpips" "0.0" "roundSolution_turbine";

#endif
