/*
 * GondzioStochSolver.C
 *
 *  Created on: Dec 7, 2017
 *      Author: Daniel Rehfeldt
 */



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
double gmu;

// double grnorm;
extern int gOoqpPrintLevel;

double g_iterNumber;


GondzioStochSolver::GondzioStochSolver( ProblemFormulation * opt, Data * prob, unsigned int n_linesearch_points )
  : GondzioSolver(opt, prob), n_linesearch_points(n_linesearch_points) // todo don't call parent constructor
{
   assert(n_linesearch_points > 0);

   // the two StepFactor constants set targets for increase in step
   // length for each corrector
   StepFactor0 = 0.08; // todo change
   StepFactor1 = 1.08; // todo change
   corrector_weight = 0.0;

   temp_step = factory->makeVariables(prob);
}

double GondzioStochSolver::correctorWeight(Variables *iterate, Variables* predictor_step, Variables* corrector_step, double predictor_alpha)
{
   assert(predictor_alpha > 0.0 && predictor_alpha <= 1.0);

   double alpha_best = -1.0;
   double weight_best = 1.0;
   const double weight_min = predictor_alpha * predictor_alpha;
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

      //  std::cout << "curr alpha: " << alpha_curr << std::endl;

      if( alpha_curr > alpha_best )
      {
         alpha_best = alpha_curr;
         weight_best = weight_curr;
      }
   }

   assert(alpha_best > 0.0 && alpha_best <= 1.0);

   return weight_best;
}

int GondzioStochSolver::solve(Data *prob, Variables *iterate, Residuals * resid )
{
   int done;
   double mu, muaff;
   int StopCorrections;
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

      // line search on corrector step todo corrector_weight is member!
      corrector_weight = correctorWeight(iterate, step, corrector_step, alpha);
     // corrector_weight = 1.0;

      // calculate weighted predictor-corrector step
      step->saxpy(corrector_step, corrector_weight);

      alpha = iterate->stepbound(step);

      std::cout << "alpha " << alpha << std::endl;


      // prepare for Gondzio corrector loop: zero out the
      // corrector_resid structure:
      corrector_resid->clear_r1r2();

      // calculate the target box:
      rmin = sigma * mu * beta_min;
      rmax = sigma * mu * beta_max;

      StopCorrections = 0;
      NumberGondzioCorrections = 0;

      // enter the Gondzio correction loop:
      while( NumberGondzioCorrections < maximum_correctors && alpha < 1.0
            && !StopCorrections )
      {

         // copy current variables into corrector_step
         corrector_step->copy(iterate);

         // calculate target steplength
         alpha_target = StepFactor1 * alpha + StepFactor0;
         if( alpha_target > 1.0 )
            alpha_target = 1.0;

         // add a step of this length to corrector_step
         corrector_step->saxpy(step, alpha_target);

         // place XZ into the r3 component of corrector_resids
         corrector_resid->set_r3_xz_alpha(corrector_step, 0.0);

         // do the projection operation
         corrector_resid->project_r3(rmin, rmax);

         // solve for corrector direction
         sys->solve(prob, iterate, corrector_resid, corrector_step);


         // todo do line search on corrector_step

         // add the current step to corrector_step, and calculate the
         // step to boundary along the resulting direction
         corrector_step->saxpy(step, 1.0);
         alpha_enhanced = iterate->stepbound(corrector_step);

         // if the enhanced step length is actually 1, make it official
         // and stop correcting
         if( alpha_enhanced == 1.0 )
         {
            step->copy(corrector_step);
            alpha = alpha_enhanced;
            NumberGondzioCorrections++;
            StopCorrections = 1;
         }
         else if( alpha_enhanced >= (1.0 + AcceptTol) * alpha )
         {
            // if enhanced step length is significantly better than the
            // current alpha, make the enhanced step official, but maybe
            // keep correcting
            step->copy(corrector_step);
            alpha = alpha_enhanced;
            NumberGondzioCorrections++;
            StopCorrections = 0;
         }
         else
         {
            // otherwise quit the correction loop
            StopCorrections = 1;
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
}
