
// Use this datum for data-aware scheduling
@schedule blob d[];

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
  blob s = q(j);
  blob t[];
  foreach i in [0 : iterations-1]
  {
    // Run extension function
    t[i] = f(i, s, d[i]);
  }
  // Call extension function to merge results
  // m = t[0] & t[1] & t[2] ...
  blob m = merge(t);
  v[j] = m;
}

float result = min_float(v);
printf("result: %f\n", result);
