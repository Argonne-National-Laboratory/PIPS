/*
 * GondzioStochLpSolver.C
 *
 *  Created on: 20.12.2017
 *      Author: Daniel Rehfeldt, Svenja Uslu
 */

//#define PIPS_DEBUG
#include "pipsdef.h"
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

#include "sData.h"
#include "sVars.h"

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
#include "sLinsysRootAug.h"

#include <fstream>
#include <string>
#include <sstream>
#include <limits>

extern int gOoqpPrintLevel;
extern double g_iterNumber;
extern bool ipStartFound;


GondzioStochLpSolver::GondzioStochLpSolver( ProblemFormulation * opt, Data * prob)
  : GondzioStochSolver(opt, prob)
{
}

void GondzioStochLpSolver::calculateAlphaPDWeightCandidate(Variables *iterate, Variables* predictor_step,
		Variables* corrector_step, double alpha_primal, double alpha_dual,
		double& alpha_primal_candidate, double& alpha_dual_candidate,
		double& weight_primal_candidate, double& weight_dual_candidate)
{
   assert(alpha_primal > 0.0 && alpha_primal <= 1.0);
   assert(alpha_dual > 0.0 && alpha_dual <= 1.0);

   double alpha_primal_best = -1.0, alpha_dual_best = -1.0;
   double weight_primal_best = -1.0, weight_dual_best = -1.0;
   const double weight_min = alpha_primal * alpha_dual;
   const double weight_intervallength = 1.0 - weight_min;

   // main loop
   for( unsigned int n = 0; n <= n_linesearch_points; n++ )
   {
      double weight_curr = weight_min
            + (weight_intervallength / (n_linesearch_points)) * n;

      weight_curr = min(weight_curr, 1.0);

      assert(weight_curr > 0.0 && weight_curr <= 1.0);

      temp_step->copy(predictor_step);
      temp_step->saxpy(corrector_step, weight_curr);

      double alpha_primal_curr = 1.0, alpha_dual_curr = 1.0;
      iterate->stepbound_pd(temp_step, alpha_primal_curr, alpha_dual_curr);
      assert(alpha_primal_curr > 0.0 && alpha_primal_curr <= 1.0);
      assert(alpha_dual_curr > 0.0 && alpha_dual_curr <= 1.0);

      if( alpha_primal_curr > alpha_primal_best )
      {
         alpha_primal_best = alpha_primal_curr;
         weight_primal_best = weight_curr;
      }
      if( alpha_dual_curr > alpha_dual_best )
      {
         alpha_dual_best = alpha_dual_curr;
         weight_dual_best = weight_curr;
      }
   }

   assert(alpha_primal_best >= 0.0 && weight_primal_best >= 0.0);
   assert(alpha_dual_best >= 0.0 && weight_dual_best >= 0.0);

   weight_primal_candidate = weight_primal_best;
   weight_dual_candidate = weight_dual_best;

   alpha_primal_candidate = alpha_primal_best;
   alpha_dual_candidate = alpha_dual_best;
}

int GondzioStochLpSolver::solve(Data *prob, Variables *iterate, Residuals * resid )
{
   double mu;
   int status_code;
   double sigma = 1.0;
   double alpha_pri = 1.0, alpha_dual = 1.0;
   QpGenStoch* stochFactory = reinterpret_cast<QpGenStoch*>(factory);
   g_iterNumber = 0.0;

   dnorm = prob->datanorm();
   // initialization of (x,y,z) and factorization routine.
   sys = factory->makeLinsys(prob);

   setBiCGStabTol(-1);

   stochFactory->iterateStarted();
   this->start(factory, iterate, prob, resid, step);
   stochFactory->iterateEnded();

   assert(!ipStartFound);
   ipStartFound = true;
   iter = 0;
   NumberGondzioCorrections = 0;
   mu = iterate->mu();

   while( true )
   {
      iter++;
      setBiCGStabTol(iter);

      stochFactory->iterateStarted();

      // evaluate residuals and update algorithm status:
      resid->calcresids(prob, iterate);

      //  termination test:
      status_code = this->doStatus(prob, iterate, resid, iter, mu, 0);

      if( status_code != NOT_FINISHED )
         break;

      if( gOoqpPrintLevel >= 10 )
      {
         this->doMonitorPd(prob, iterate, resid, alpha_pri, alpha_dual, sigma,
               iter, mu, status_code, 0);
      }

      // *** Predictor step ***
      resid->set_r3_xz_alpha(iterate, 0.0);

      sys->factor(prob, iterate);
      sys->solve(prob, iterate, resid, step);
      step->negate();

      iterate->stepbound_pd(step, alpha_pri, alpha_dual);

      // calculate centering parameter
      const double muaff = iterate->mustep_pd(step, alpha_pri, alpha_dual);

      assert(!PIPSisZero(mu));
      sigma = pow(muaff / mu, tsig);

      if( gOoqpPrintLevel >= 10 )
      {
         this->doMonitorPd(prob, iterate, resid, alpha_pri, alpha_dual, sigma,
               iter, mu, status_code, 2);
      }

      g_iterNumber += 1.0;

      // *** Corrector step ***
      corrector_resid->clear_r1r2();

      // form right hand side of linear system:
      corrector_resid->set_r3_xz_alpha(step, -sigma * mu);

      sys->solve(prob, iterate, corrector_resid, corrector_step);
      corrector_step->negate();

      // calculate weighted predictor-corrector step
      double weight_primal_candidate, weight_dual_candidate = -1.0;

      calculateAlphaPDWeightCandidate(iterate, step, corrector_step, alpha_pri,
            alpha_dual, alpha_pri, alpha_dual, weight_primal_candidate, weight_dual_candidate);

      assert(weight_primal_candidate >= 0.0 && weight_primal_candidate <= 1.0);
      assert(weight_dual_candidate >= 0.0 && weight_dual_candidate <= 1.0);

      step->saxpy_pd(corrector_step, weight_primal_candidate,
            weight_dual_candidate);

      // prepare for Gondzio corrector loop: zero out the corrector_resid structure:
      corrector_resid->clear_r1r2();

      // calculate the target box:
      const double rmin = sigma * mu * beta_min;
      const double rmax = sigma * mu * beta_max;

      NumberGondzioCorrections = 0;

      // enter the Gondzio correction loop:
      while( NumberGondzioCorrections < maximum_correctors
            && (PIPSisLT(alpha_pri, 1.0) || PIPSisLT(alpha_dual, 1.0)) )
      {
         corrector_step->copy(iterate);

         double alpha_pri_enhanced, alpha_dual_enhanced;
         const double alpha_pri_target = std::min(1.0, StepFactor1 * alpha_pri + StepFactor0);
         const double alpha_dual_target = std::min(1.0, StepFactor1 * alpha_dual + StepFactor0);

         PIPSdebugMessage("corrector loop: %d alpha_pri: %f alpha_dual %f \n", NumberGondzioCorrections, alpha_pri, alpha_dual);

         // add a step of this length to corrector_step
         corrector_step->saxpy_pd(step, alpha_pri_target, alpha_dual_target);
         // corrector_step is now x_k + alpha_target * delta_p (a trial point)

         // place XZ into the r3 component of corrector_resids
         corrector_resid->set_r3_xz_alpha(corrector_step, 0.0);

         // do the projection operation
         corrector_resid->project_r3(rmin, rmax);

         // solve for corrector direction
         sys->solve(prob, iterate, corrector_resid, corrector_step); // corrector_step is now delta_m

         // calculate weighted predictor-corrector step
         calculateAlphaPDWeightCandidate(iterate, step, corrector_step,
               alpha_pri_target, alpha_dual_target, alpha_pri_enhanced,
               alpha_dual_enhanced, weight_primal_candidate, weight_dual_candidate);

         // if the enhanced step length is actually 1, make it official
         // and stop correcting
         if( PIPSisEQ(alpha_pri_enhanced, 1.0) && PIPSisEQ(alpha_dual_enhanced, 1.0) )
         {
            PIPSdebugMessage("both 1.0 \n");

            step->saxpy_pd(corrector_step, weight_primal_candidate, weight_dual_candidate);

            alpha_pri = alpha_pri_enhanced;
            alpha_dual = alpha_dual_enhanced;
            NumberGondzioCorrections++;

            // exit Gondzio correction loop
            break;
         }
         else if( alpha_pri_enhanced >= (1.0 + AcceptTol) * alpha_pri
               && alpha_dual_enhanced >= (1.0 + AcceptTol) * alpha_dual )
         {
            PIPSdebugMessage("both better \n");

            // if enhanced step length is significantly better than the
            // current alpha, make the enhanced step official, but maybe
            // keep correcting
            step->saxpy_pd(corrector_step, weight_primal_candidate, weight_dual_candidate);

            alpha_pri = alpha_pri_enhanced;
            alpha_dual = alpha_dual_enhanced;
            NumberGondzioCorrections++;
         }
         else if( alpha_pri_enhanced >= (1.0 + AcceptTol) * alpha_pri )
         {
            PIPSdebugMessage("primal better \n");

            step->saxpy_pd(corrector_step, weight_primal_candidate, 0.0);

            alpha_pri = alpha_pri_enhanced;
            NumberGondzioCorrections++;
         }
         else if( alpha_dual_enhanced >= (1.0 + AcceptTol) * alpha_dual )
         {
            PIPSdebugMessage("dual better \n");

            step->saxpy_pd(corrector_step, 0.0, weight_dual_candidate);

            alpha_dual = alpha_dual_enhanced;
            NumberGondzioCorrections++;
         }
         else
         {
            // exit Gondzio correction loop
            break;
         }

         temp_step->copy(iterate);
         double alpha_pri_copy = alpha_pri;
         double alpha_dual_copy = alpha_dual;
         finalStepLength_PD(temp_step, step, alpha_pri_copy, alpha_dual_copy);

         temp_step->saxpy_pd(step, alpha_pri, alpha_dual);

         if( PIPS_MPIgetRank() == 0 && temp_step->mu() < 0 )
         {
            std::cout << "SOMETHING WENT TERRIBLY WRONG IN THE CORRECTORS ..." << std::endl;
            // dump sc
            //dynamic_cast<sLinsysRootAug*>(sys)->dumpKKT(iter);
         }

      }

      // We've finally decided on a step direction, now calculate the
      // length using Mehrotra's heuristic.x
      temp_step->copy(iterate);
      finalStepLength_PD(iterate, step, alpha_pri, alpha_dual);
      // actually take the step and calculate the new mu

      iterate->saxpy_pd(step, alpha_pri, alpha_dual);
      mu = iterate->mu();

      if( mu < 0 )
      {
         double alp = alpha_pri;
         double ald = alpha_dual;
         temp_step->stepbound_pd( step, alp, ald);
         temp_step->saxpy_pd(step, alp, ald);
         double muc = temp_step->mu();
         if( PIPS_MPIgetRank() == 0 )
         {
            std::cout << "mu: " << mu << " != " << muc << std::endl;
            std::cout << "alp: " << alpha_pri << " != " << alp << std::endl;
            std::cout << "ald: " << alpha_dual << " != " << ald << std::endl;
         }
      }
      if( PIPS_MPIgetRank() == 0 && mu < 0 )
      {
         std::cout << "SOMETHING WENT TERRIBLY WRONG..." << std::endl;
         // dump sc
         dynamic_cast<sLinsysRootAug*>(sys)->dumpKKT(iter);
      }

      stochFactory->iterateEnded();
   }

   resid->calcresids(prob, iterate);
   if( gOoqpPrintLevel >= 10 )
   {
      this->doMonitorPd(prob, iterate, resid, alpha_pri, alpha_dual, sigma, iter, mu, status_code, 1);
   }

   return status_code;
}


GondzioStochLpSolver::~GondzioStochLpSolver()
{
}
