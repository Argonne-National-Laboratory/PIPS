#ifndef RAND_UTILS
#define RAND_UTILS

#include "time.h"

/* FORTRAN functions */

/** 
 * Generates normal numbers: mean=0, std=1
 *
 * 'zufalli_' must be called before to initialize the seed.
 */ 
extern "C"
void normalen_(int* n, double* numbers);

/**
 * Initializes the seed. If called with 0, the seed is initialized
 * to 1802.
 */
extern "C"
void zufalli_(int*);

/** 
 * Generates uniform numbers in [0,1]
 *
 * 'zufalli_' must be called before to initialize the seed.
 */ 
extern "C"
void zufall_(int* n, double* numbers);

/* C++ wrappers */

/** 
 * Generates normal numbers: mean=0, std=1
*/
inline void normalRand(int n, double* vals)
{
  int tm=time(NULL);
  zufalli_(&tm);
  normalen_(&n, vals);
}
inline void normalRand(double mean, double std, int n, double* vals)
{
  int tm=time(NULL);
  zufalli_(&tm);
  normalen_(&n, vals);

  for(int i=0; i<n; i++)
    vals[i] = mean+std*vals[i];
}
#endif
