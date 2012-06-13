
// Use this datum for data-aware scheduling (exm-issues #193)
@heavy blob d[];

int data_sets = 10;
foreach i in [0 : data_sets-1]
{
  // Call extension function to read data from storage
  d[i] = read_data_set(i);
}

float v[];
int trials = 1000;
int iterations = 10;
foreach j in [0 : trials-1]
{
  updateable_blob s = q(j);
  foreach i in [0 : iterations-1]
  {
    // Run extension function
    blob t = f(i, s, d[i]);

    // s = s & t
    s <merge_ones> := t;
  }
  v[j] = g(s);
}

// Using proposed timeout features (exm-issues #244)
// cutoff in seconds
float cutoff = 55*60;
// Force-close v when time expires
force v := timeout(cutoff);

float result = min_float(v);
printf("result: %f\n", result);

// sleep: Blocks on inputs, sleeps for 10 seconds, sets output
void e = sleep(10, result);
// exit: Blocks on input, nicely shuts down Swift/Turbine
exit(e);
