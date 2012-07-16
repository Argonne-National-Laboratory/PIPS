
#include <builtins.swift>
#include <swift/stdio.swift>

#include "rounding_functions.swift"

main
{
  blob result;
  blob candidateSolution;
  float t[];
  t[0] = 3;
  t[1] = 2;
  t[2] = 100;
  candidateSolution = blob_from_floats(t);

  result = evaluateRecourseLP("/path/to/data", 3, 2, candidateSolution);

  printf("OK");
}
