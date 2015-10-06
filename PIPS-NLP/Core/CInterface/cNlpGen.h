/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */
/* 2015. Modified by Nai-Yuan Chiang for NLP */

#ifndef CNLPGEN
#define CNLPGEN

#include "OoqpMonitorData.h"

/** A data structure for holding the various bits needed to solve
 *  a general NLP. This structure is part of the CInterface
 *  to using NlpGen */
typedef struct {
  void * factory;
  void * prob;
  void * solver;
} NlpGenContext;

#ifdef __cplusplus
extern "C" {
#endif
  void NlpGenFinish ( NlpGenContext * ctx,
		     double   x[],  double gamma[],     double phi[],
		     double   y[], 
		     double   z[],  double lambda[],    double pi[],
		     double *objectiveValue,
		     int * status_code );

  void NlpGenCleanup( NlpGenContext * ctx );

  void NlpGenAddMonitor( NlpGenContext * ctx, DoItCFunc cmon,
			void * mctx );

  void NlpGenMonitorSelf( NlpGenContext * ctx );
#ifdef __cplusplus
}
#endif

#endif
