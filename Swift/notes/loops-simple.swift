
float v[];
int trials = 1000;
int iterations = 10;
foreach j in [0 : trials-1]
{
  blob s = q(j);
  blob t[];
  foreach i in [0 : iterations-1]
  {
    // Run extension function (reads data based on i)
    t[i] = f(i, s);
  }
  // Call extension function to merge results
  // m = t[0] & t[1] & t[2] ...
  blob m = merge(t);
  v[j] = m;
}

float result = min_float(v);
printf("result: %f\n", result);
