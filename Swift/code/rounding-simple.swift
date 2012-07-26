
/*
 * ROUNDING-SIMPLE.SWIFT
 *
 * Swift/T application script
 */

#include <builtins.swift>
#include <swift/stdio.swift>
#include <swift/stats.swift>

#include "rounding_functions.swift"

main {
  // problem data
  string data = "/home/wozniak/Public/PIPS.data/";
  string dataPath = data + "uc_dumps/uc_raw-4h";
  string solutionPath = data + "primalsol_conv8";
  int nScenarios = 4;
  blob s = readConvSolution(dataPath,solutionPath);

  // eventually we want to do a sweep over different values of cutoff
  // and possibly other functions in the place of "round"
  float cutoff = 0.5;
  blob r = roundSolution(s,cutoff);

  // inner loop, given "r"
  float v[];
  foreach i in [0 : nScenarios-1]
  {
    v[i] = evaluateRecourseLP(dataPath,nScenarios,i,r);
  }
  float result = sum_float(v);
  printf("result: %f\n", result);
}
