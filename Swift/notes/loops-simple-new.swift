
# problem data
string dataPath = "...";
string solutionPath = "...";
int nScenarios = ...;
blob s = readConvSolution(dataPath,solutionPath);

# eventually we want to do a sweep over different values of cutoff
# and possibly other functions in the place of "round"
float cutoff = 0.5;
blob r = round(s,cutoff);

# inner loop, given "r"
float v[];
foreach i in [0 : nScenarios-1]
{
  v[i] = evaluateRecourseLP(dataPath,nScenarios,i,r);
}
float result = sum(v);
printf("result: %f\n", result);
