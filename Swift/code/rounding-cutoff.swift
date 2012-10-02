
/*
 * ROUNDING-CUTOFF.SWIFT
 *
 * Swift/T application script
 */

#include <builtins.swift>
#include <swift/mpe.swift>
#include <swift/stdio.swift>
#include <swift/stdlib.swift>
#include <swift/stats.swift>
#include <swift/string.swift>
#include <swift/unistd.swift>

#include "pips.swift"

(float result[]) cutoffs(int count)
{
  float c = itof(count);
  result[0] = 0.0;
  foreach i in [1:count-2]
  {
    result[i] = itof(i)/c;
  }
  result[count-1] = 1.0;
}

main
{
  // problem data
  string data = "/home/wozniak/PIPS-data-2/";
  // string dataPath = data + "uc_dumps/uc_raw-4h";
  string dataPath = data + "4h_dump/uc_4h";
  string solutionPath = data + "primalsol_conv8";
  int nScenarios = toint(argv("N"));
  blob s = readConvSolution(dataPath,solutionPath);
  int nCutoffs = toint(argv("C"));

  float cutoffs[] = cutoffs(nCutoffs);

  foreach cutoff in cutoffs
  {
    blob r = roundSolution(s, cutoff);

    // inner loop, given "r"
    float v[];
    foreach i in [0 : nScenarios-1]
    {
      v[i] = evaluateRecourseLP(dataPath,nScenarios,i,r);
      // printf("v[%i]=%f", i, v[i]);
    }
    float result = sum_float(v);
    // printf("result: %f\n", result);
  }

  // Output some metadata
  int procs = adlb_servers() + turbine_engines() + turbine_workers();
  metadata(sprintf("PROCS: %i", procs));
  metadata(sprintf("ENGINES: %i", turbine_engines()));
  metadata(sprintf("WORKERS: %i", turbine_workers()));
  metadata(sprintf("SERVERS: %i", adlb_servers()));
  metadata(sprintf("JOB: %s", getenv("COBALT_JOBID")));
}
