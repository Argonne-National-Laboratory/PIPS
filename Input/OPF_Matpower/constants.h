#ifndef CONSTANTS_H
#define CONSTANTS_H


#include <math.h>



#define w_s (2*M_PI*freq)
#define epsilon 1E-8
#define PS_MAXLINE 1000
#define ISOLATED_BUS 4
#define REF_BUS 3
#define PV_BUS 2
#define PQ_BUS 1
#define NGEN_AT_BUS_MAX 15
#define NLOAD_AT_BUS_MAX 10
#define NPHASE 1

/* Type of variables */
#define DIFF_EQ 1 /* Differential equation */
#define ALG_EQ  0 /* Algebraic equation */


#define MAXCONNLINES 20

#endif
