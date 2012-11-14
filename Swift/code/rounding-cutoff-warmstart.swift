//the SWIFT syntax may be wrong, don't worry
blob s = readConvSolution(dataPath,solutionPath);

float cutoffs[] = {0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8,0.9,
                   1.0}; //the cuttoffs
blob r_sol = empty; //the rounded solution I am looking for (i.e,
feasible and giving the smallest recourse value over all cutoffs)
float rec_min = +inf;
blob st_point = empty; //the point used to warmstart

//outer loop over cutoffs
foreach cutoff in cutoffs
{
  blob r = roundSolution(s,cutoff);

  float v[];
  blob st_point_aux;

// inner loop, given "r"
  foreach i in [0 : nScenarios-1]
  {
    [v[i], st_point_aux] =
      evaluateRecourseLP(dataPath,nScenarios,i,r,st_point); //the C function
    will warm-start if st_point is nonempty
                                  if (v[i]<+inf)
                                  {
                                    st_point=st_point_aux; //warm-start from now on.
                                  }
  }

//check feasibility (=all scenarios return finite values)
  float rec = sum_float(v);
  if (rec < +inf ){

// r is feasible, check whether improves the recourse
    if(rec<rec_min) {
      rec_min = rec;
      r_sol = r;
    }
  } //end of loop over scenarios
} // end of loop over cutoffs
if (r_sol==empty)
  printf("Could not found a feasible point with the specified
cutoffs.");
else
  printf("Best upper bound is:%f", rec_min);
