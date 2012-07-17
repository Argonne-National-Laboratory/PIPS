
#include <builtins.swift>
#include <swift/stdio.swift>

#include "rounding_functions.swift"

main {
  // problem data
  string dataPath = "file.data";
  string solutionPath = "file.solution";
  int nScenarios = 3;
  blob s = readConvSolution(dataPath,solutionPath);

  // eventually we want to do a sweep over different values of cutoff
  // and possibly other functions in the place of "round"
  float cutoff = 0.5;
  blob r = round(s,cutoff);

  // inner loop, given "r"
  float v[];
  foreach i in [0 : nScenarios-1]
  {
    // TODO: Should this return a float or a blob?
    v[i] = evaluateRecourseLP(dataPath,nScenarios,i,r);
  }
  float result = sum_float(v);
  printf("result: %f\n", result);
}
