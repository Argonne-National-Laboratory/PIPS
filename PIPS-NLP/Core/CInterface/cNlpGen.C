/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */
/* 2015. Modified by Nai-Yuan Chiang for NLP */


#include "cNlpGen.h"
#include "NlpGenData.h"
#include "NlpGen.h"
#include "Solver.h"
#include "NlpGenVars.h"
#include "NlpGenResiduals.h"
#include "SimpleVector.h"
#include "OoqpMonitor.h"

#include <cstdlib>

extern "C" 
void NlpGenFinish ( NlpGenContext * ctx,
		   double   x[],  double gamma[],     double phi[],
		   double   y[], 
		   double   z[],  double lambda[],    double pi[],
		   double * objectiveValue, int * status_code )
{
  NlpGenData    * prob    = (NlpGenData *)    ctx->prob;
  NlpGen * factory = (NlpGen *) ctx->factory;
  Solver       * solver  = (Solver *)       ctx->solver;
  
  NlpGenVars      * vars    = 0;
  NlpGenResiduals * resid   = 0;
  
  try {
    vars  = (NlpGenVars * )     factory->makeVariables( prob );
    resid = (NlpGenResiduals *) factory->makeResiduals( prob );
    
    *status_code = solver->solve(prob,vars,resid);

    vars->x     ->copyIntoArray( x );
    vars->gamma ->copyIntoArray( gamma );
    vars->phi ->copyIntoArray( phi );

    vars->y->copyIntoArray( y );

    vars->z     ->copyIntoArray( z );
    vars->lambda->copyIntoArray( lambda );
    vars->pi    ->copyIntoArray( pi );
    *objectiveValue = prob->objectiveValue( vars );
  }
  catch( ... ) {
    *status_code = -1;
  }
  delete vars;
  delete resid;
}

extern "C"
void NlpGenCleanup( NlpGenContext * ctx )
{
  NlpGen * factory = (NlpGen *) ctx->factory;
  NlpGenData    * prob    = (NlpGenData *)    ctx->prob;
  Solver       * solver  = (Solver *)       ctx->solver;
  
  delete factory;       delete prob;       delete solver;
  ctx->factory = 0;     ctx->prob = 0;     ctx->solver = 0;
}

extern "C"
void NlpGenAddMonitor( NlpGenContext * ctx, DoItCFunc cmon,
		      void * mctx )
{
  COoqpMonitor * mon = new COoqpMonitor( cmon, mctx );
  ((Solver *) ctx->solver)->addMonitor( mon );
}

extern "C"
void NlpGenMonitorSelf( NlpGenContext * ctx )
{
  ((Solver *) ctx->solver)->monitorSelf();
}

