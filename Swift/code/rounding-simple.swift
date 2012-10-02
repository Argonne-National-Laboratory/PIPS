
/*
 * ROUNDING-SIMPLE.SWIFT
 *
 * Swift/T application script
 */

#include <builtins.swift>
#include <swift/stdio.swift>
#include <swift/stats.swift>
#include <swift/unistd.swift>

#include "pips.swift"

main {
  // problem data
  string data = "/home/wozniak/PIPS-data-2/";
  // string dataPath = data + "uc_dumps/uc_raw-4h";
  string dataPath = data + "4h_dump/uc_4h";
  string solutionPath = data + "primalsol_conv8";
  int nScenarios = toint(argv("N"));
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
    printf("v[%i]=%f", i, v[i]);
 }
  float result = sum_float(v);
  printf("result: %f\n", result);

  // Output some metadata
  int procs = adlb_servers() + turbine_engines() + turbine_workers();
  printf("procs: %i", procs);
}
