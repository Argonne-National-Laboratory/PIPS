
#include <builtins.swift>
#include "functions.swift"

main {
  string dataPath = "problem.data";
  string solutionPath = "NOTHING";
  int scen = 3;

  blob b = readSolution(dataPath, solutionPath, scen);
}

