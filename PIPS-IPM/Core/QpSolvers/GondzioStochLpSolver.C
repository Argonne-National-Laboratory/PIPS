/*
 * GondzioStochLpSolver.C
 *
 *  Created on: 20.12.2017
 *      Author: Svenja Uslu
 */


#include "GondzioStochLpSolver.h"
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

static double g_iterNumber;


GondzioStochLpSolver::GondzioStochLpSolver( ProblemFormulation * opt, Data * prob, unsigned int n_linesearch_points )
  : GondzioStochSolver(opt, prob)
{


   //temp_step = factory->makeVariables(prob);
}

void GondzioStochLpSolver::calculateAlphaWeightCandidate(Variables *iterate, Variables* predictor_step, Variables* corrector_step,
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

void GondzioStochLpSolver::calculateAlphaPDWeightCandidate(Variables *iterate, Variables* predictor_step, Variables* corrector_step, double alpha_primal, double alpha_dual, double& alpha_primal_candidate, double& alpha_dual_candidate, double& weight_candidate)
{
	// todo: implement method here. think of alpha_pri * alpha_dual
}

int GondzioStochLpSolver::solve(Data *prob, Variables *iterate, Residuals * resid )
{
   int done;
   double mu, muaff;
   double alpha_target, alpha_enhanced, rmin, rmax;
   int status_code;
   double alpha = 1, sigma = 1;
   double alpha_pri = 1, alpha_dual = 1;
   double alpha_pri_target, alpha_dual_target;
   double alpha_pri_enhanced, alpha_dual_enhanced;
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

      //todo: compute primal and dual steplength (alpha_primal, alpha_dual)
      //probably also adapt calculation of centering parameter

      // instead of: alpha = iterate->stepbound(step); do:
      iterate->stepbound_primal_dual(step, alpha_pri, alpha_dual);

      // calculate centering parameter
      // instead of: muaff = iterate->mustep(step, alpha); do:
      muaff = iterate->mustep_pd(step, alpha_pri, alpha_dual);
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
      // instead of: calculateAlphaWeightCandidate(iterate, step, corrector_step, alpha, alpha, weight_candidate); do:
      calculateAlphaPDWeightCandidate(iterate, step, corrector_step, alpha_pri, alpha_dual,
    		  alpha_pri, alpha_dual, weight_candidate);

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
      // instead of: while( NumberGondzioCorrections < maximum_correctors && alpha < 1.0 ) do:
      // todo: && or || between both alphas ?
      while( NumberGondzioCorrections < maximum_correctors && alpha_pri < 1.0 && alpha_dual < 1.0 )
      {

         // copy current variables into corrector_step
         corrector_step->copy(iterate);

         // calculate target steplength
         // instead of: alpha_taget = StepFactor1 * alpha + StepFactor0; do:
         alpha_pri_target = StepFactor1 * alpha_pri + StepFactor0;
         alpha_dual_target = StepFactor1 * alpha_dual + StepFactor0;
         if( alpha_pri_target > 1.0 )
            alpha_pri_target = 1.0;
         if( alpha_dual_target > 1.0 )
            alpha_dual_target = 1.0;

         // add a step of this length to corrector_step
         // instead of: corrector_step->saxpy(step, alpha_target); do:
         corrector_step->saxpy_pd(step, alpha_pri_target, alpha_dual_target);
         // corrector_step is now x_k + alpha_target * delta_p (a trial point)

         // place XZ into the r3 component of corrector_resids
         corrector_resid->set_r3_xz_alpha(corrector_step, 0.0);

         // do the projection operation
         corrector_resid->project_r3(rmin, rmax);

         // solve for corrector direction
         sys->solve(prob, iterate, corrector_resid, corrector_step);	// corrector_step is now delta_m

         // calculate weighted predictor-corrector step
         calculateAlphaPDWeightCandidate(iterate, step, corrector_step, alpha_pri_target, alpha_dual_target,
        		 alpha_pri_enhanced, alpha_dual_enhanced, weight_candidate);

         // todo weight * corrector_step; corrector_step += step
         temp_step->copy(step);
         temp_step->saxpy(corrector_step, weight_candidate);

         // if the enhanced step length is actually 1, make it official
         // and stop correcting
         // todo: is || correct or && ?
         if( alpha_pri_enhanced == 1.0 || alpha_dual_enhanced == 1.0 )
         {
            step->copy(temp_step);
            alpha_pri = alpha_pri_enhanced;
            alpha_dual = alpha_dual_enhanced;
            NumberGondzioCorrections++;

            // exit Gondzio correction loop
            break;
         }
         else if( alpha_pri_enhanced >= (1.0 + AcceptTol) * alpha_pri || alpha_dual_enhanced >= (1.0 + AcceptTol) * alpha_dual)
         {
            // if enhanced step length is significantly better than the
            // current alpha, make the enhanced step official, but maybe
            // keep correcting
            step->copy(temp_step);
            alpha_pri = alpha_pri_enhanced;
            alpha_dual = alpha_dual_enhanced;
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
      // instead of: alpha = finalStepLength(iterate, step); do:
      finalStepLength_PD(iterate, step, alpha_pri, alpha_dual);

      // alternatively, just use a crude step scaling factor.
      // alpha = 0.995 * iterate->stepbound( step );

      // actually take the step (at last!) and calculate the new mu

      // instead of: iterate->saxpy(step, alpha); do:
      iterate->saxpy_pd(step, alpha_pri, alpha_dual);
      mu = iterate->mu();
      gmu = mu;

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


GondzioStochLpSolver::~GondzioStochLpSolver()
{
}



