/*
 * GondzioStochSolver.C
 *
 *  Created on: Dec 7, 2017
 *      Author: Daniel Rehfeldt
 */


#include "pipsdef.h"
#include "GondzioStochSolver.h"
#include "Variables.h"
#include "Residuals.h"
#include "LinearSystem.h"
#include "Status.h"
#include "Data.h"
#include "ProblemFormulation.h"

#include "OoqpVector.h"
#include "DoubleMatrix.h"

#include "StochTree.h"
#include "QpGenStoch.h"
#include "StochResourcesMonitor.h"

#include "sVars.h"
#include "sData.h"

#include <cstring>
#include <iostream>
#include <fstream>
#include <algorithm>
using namespace std;

#include <cstdio>
#include <cassert>
#include <cmath>

#include "StochVector.h"
#include "mpi.h"
#include "QpGenVars.h"
#include "QpGenResiduals.h"

// gmu is needed by MA57!
static double gmu;

// double grnorm;
extern int gOoqpPrintLevel;

double g_iterNumber;


GondzioStochSolver::GondzioStochSolver( ProblemFormulation * opt, Data * prob, unsigned int n_linesearch_points,
      bool adaptive_linesearch )
  : GondzioSolver(opt, prob), n_linesearch_points(n_linesearch_points)
{
   assert(n_linesearch_points > 0);

   if( adaptive_linesearch )
   {
      int size;
      MPI_Comm_size(MPI_COMM_WORLD, &size);

      if( size > 1)
         this->n_linesearch_points =
               std::min(unsigned(size) + this->n_linesearch_points, max_linesearch_points);
   }

   // the two StepFactor constants set targets for increase in step
   // length for each corrector
   StepFactor0 = 0.3;
   StepFactor1 = 1.5;

#ifdef REDUCED_ACCURACY
   artol = 1.e-3;
   mutol = 1.e-5;
#else
   mutol = 1.e-6; // todo parameter
#endif

   temp_step = factory->makeVariables(prob);
}

void GondzioStochSolver::calculateAlphaWeightCandidate(Variables *iterate, Variables* predictor_step, Variables* corrector_step,
      double alpha_predictor, double& alpha_candidate, double& weight_candidate)
{
   assert(alpha_predictor > 0.0 && alpha_predictor <= 1.0);

   double alpha_best = -1.0;
   double weight_best = -1.0;
   const double weight_min = alpha_predictor * alpha_predictor;
   const double weight_intervallength = 1.0 - weight_min;

   // main loop
   for( unsigned int n = 0; n <= n_linesearch_points; n++ )
   {
      double weight_curr = weight_min + (weight_intervallength / (n_linesearch_points)) * n;

      weight_curr = min(weight_curr, 1.0);

      assert(weight_curr > 0.0 && weight_curr <= 1.0);

      temp_step->copy(predictor_step);
      temp_step->saxpy(corrector_step, weight_curr);

      const double alpha_curr = iterate->stepbound(temp_step);
      assert(alpha_curr > 0.0 && alpha_curr <= 1.0);

      if( alpha_curr > alpha_best )
      {
         alpha_best = alpha_curr;
         weight_best = weight_curr;
      }
   }

   assert(alpha_best >= 0.0 && weight_best >= 0.0);

   weight_candidate = weight_best;
   alpha_candidate = alpha_best;
}

int GondzioStochSolver::solve(Data *prob, Variables *iterate, Residuals * resid )
{
   int done;
   double mu, muaff;
   double alpha_target, alpha_enhanced, rmin, rmax;
   int status_code;
   double alpha = 1, sigma = 1;
   QpGenStoch* stochFactory = reinterpret_cast<QpGenStoch*>(factory);
   g_iterNumber = 0.0;

   gmu = 1000;
   //  grnorm = 1000;
   dnorm = prob->datanorm();
   // initialization of (x,y,z) and factorization routine.
   sys = factory->makeLinsys(prob);

   stochFactory->iterateStarted();
   this->start(factory, iterate, prob, resid, step);
   stochFactory->iterateEnded();

   iter = 0;
   NumberGondzioCorrections = 0;
   done = 0;
   mu = iterate->mu();
   gmu = mu;
   int myRank; MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

   do
   {
      iter++;
      stochFactory->iterateStarted();

      // evaluate residuals and update algorithm status:
      resid->calcresids(prob, iterate);

      //  termination test:
      status_code = this->doStatus(prob, iterate, resid, iter, mu, 0);

      if( status_code != NOT_FINISHED )
         break;

      if( gOoqpPrintLevel >= 10 )
      {
         this->doMonitor(prob, iterate, resid, alpha, sigma, iter, mu,
               status_code, 0);
      }
      // *** Predictor step ***
      resid->set_r3_xz_alpha(iterate, 0.0);

      sys->factor(prob, iterate);
      sys->solve(prob, iterate, resid, step);
      step->negate();

      alpha = iterate->stepbound(step);

      // calculate centering parameter
      muaff = iterate->mustep(step, alpha);
      sigma = pow(muaff / mu, tsig);

      if( gOoqpPrintLevel >= 10 )
      {
         this->doMonitor(prob, iterate, resid, alpha, sigma, iter, mu,
               status_code, 2);
      }

      g_iterNumber+=0.5;

      // *** Corrector step ***

      corrector_resid->clear_r1r2();

      // form right hand side of linear system:
      corrector_resid->set_r3_xz_alpha(step, -sigma * mu);
      sys->solve(prob, iterate, corrector_resid, corrector_step);
      corrector_step->negate();

      // calculate weighted predictor-corrector step
      double weight_candidate = -1.0;
      const double alpha_predictor = alpha;
      calculateAlphaWeightCandidate(iterate, step, corrector_step, alpha_predictor, alpha, weight_candidate);
      assert(weight_candidate >= 0.0 && weight_candidate <= 1.0);

      step->saxpy(corrector_step, weight_candidate);

      // prepare for Gondzio corrector loop: zero out the
      // corrector_resid structure:
      corrector_resid->clear_r1r2();

      // calculate the target box:
      rmin = sigma * mu * beta_min;
      rmax = sigma * mu * beta_max;

      NumberGondzioCorrections = 0;

      // enter the Gondzio correction loop:
      while( NumberGondzioCorrections < maximum_correctors && PIPSisLT(alpha, 1.0) )
      {

         // copy current variables into corrector_step
         corrector_step->copy(iterate);

         // calculate target steplength
         alpha_target = StepFactor1 * alpha + StepFactor0;
         if( alpha_target > 1.0 )
            alpha_target = 1.0;

         // add a step of this length to corrector_step
         corrector_step->saxpy(step, alpha_target);
         // corrector_step is now x_k + alpha_target * delta_p (a trial point)

         // place XZ into the r3 component of corrector_resids
         corrector_resid->set_r3_xz_alpha(corrector_step, 0.0);

         // do the projection operation
         corrector_resid->project_r3(rmin, rmax);

         // solve for corrector direction
         sys->solve(prob, iterate, corrector_resid, corrector_step);	// corrector_step is now delta_m

         // calculate weighted predictor-corrector step
         calculateAlphaWeightCandidate(iterate, step, corrector_step, alpha_target, alpha_enhanced, weight_candidate);

         // if the enhanced step length is actually 1, make it official
         // and stop correcting
         if( PIPSisEQ(alpha_enhanced, 1.0) )
         {
            step->saxpy(corrector_step, weight_candidate);
            alpha = alpha_enhanced;
            NumberGondzioCorrections++;

            // exit Gondzio correction loop
            break;
         }
         else if( alpha_enhanced >= (1.0 + AcceptTol) * alpha )
         {
            // if enhanced step length is significantly better than the
            // current alpha, make the enhanced step official, but maybe
            // keep correcting
            step->saxpy(corrector_step, weight_candidate);
            alpha = alpha_enhanced;
            NumberGondzioCorrections++;
         }
         else
         {
            // exit Gondzio correction loop
            break;
         }
      }

      // We've finally decided on a step direction, now calculate the
      // length using Mehrotra's heuristic.x
      alpha = finalStepLength(iterate, step);

      // alternatively, just use a crude step scaling factor.
      // alpha = 0.995 * iterate->stepbound( step );

      // actually take the step (at last!) and calculate the new mu

      iterate->saxpy(step, alpha);
      mu = iterate->mu();
      gmu = mu;

      if( 0 == myRank )
         std::cout << "final alpha: " << alpha << " mu: " << mu <<   std::endl;

      stochFactory->iterateEnded();
   }
   while( !done );

   resid->calcresids(prob, iterate);
   if( gOoqpPrintLevel >= 10 )
   {
      this->doMonitor(prob, iterate, resid, alpha, sigma, iter, mu, status_code, 1);
   }

   // print the results, if you really want to..
   // iterate->print();

   return status_code;
}


GondzioStochSolver::~GondzioStochSolver()
{
   delete temp_step;
}
