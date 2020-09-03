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
#include "sLinsysRoot.h"

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

#include <fstream>
#include <string>
#include <sstream>
#include <limits>

extern int gOoqpPrintLevel;
extern double g_iterNumber;
extern bool ipStartFound;


GondzioStochLpSolver::GondzioStochLpSolver( ProblemFormulation * opt, Data * prob, const Scaler* scaler)
  : GondzioStochSolver(opt, prob, scaler)
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
   const int my_rank = PIPS_MPIgetRank(MPI_COMM_WORLD);

   double mu, muaff;
   int status_code;
   double sigma = 1.0;
   double alpha_pri = 1.0, alpha_dual = 1.0;
   QpGenStoch* stochFactory = reinterpret_cast<QpGenStoch*>(factory);
   g_iterNumber = 0.0;

   bool small_corr_aggr = false;
   bool pure_centering_step = false;
   bool numerical_troubles = false;
   bool precond_limit = false;

   setDnorm(*prob);
   // initialization of (x,y,z) and factorization routine.
   sys = factory->makeLinsys(prob);

   // register as observer for the BiCGStab solves
   registerBiCGStabOvserver(sys);
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
      if( !pure_centering_step )
      {
         resid->set_r3_xz_alpha(iterate, 0.0);
         sys->factor(prob, iterate);
         sys->solve(prob, iterate, resid, step);
         step->negate();

         // if we are not in a refactorization - check convergence of bicgstab
         if( !bicgstab_converged && bigcstab_norm_res_rel * 1e2 * dnorm > resid->residualNorm() )
         {
            PIPSdebugMessage("Affine step computation in BiCGStab failed");
            numerical_troubles = true;

            if( !small_corr_aggr )
            {
               if( my_rank == 0 )
                  std::cout << "switching to small correctors aggressive" << std::endl;
               small_corr_aggr = true;
            }
         }
      }
      else
         step->setToZero();

      iterate->stepbound_pd(step, alpha_pri, alpha_dual);

      // calculate centering parameter
      muaff = iterate->mustep_pd(step, alpha_pri, alpha_dual);

      assert( !PIPSisZero(mu) );
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

      // if we are not in a refactorization - check convergence of bicgstab
      if( !bicgstab_converged && bigcstab_norm_res_rel * 1e2 * dnorm > resid->residualNorm() )
      {
         PIPSdebugMessage("corrector step computation in BiCGStab failed");
         numerical_troubles = true;
         if( !small_corr_aggr )
         {
            if( my_rank == 0 )
               std::cout << "switching to small correctors aggressive" << std::endl;
            small_corr_aggr = true;
         }
      }

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
      NumberSmallCorrectors = 0;

      // if small_corr_aggr only try small correctors
      bool small_corr = small_corr_aggr;

      // enter the Gondzio correction loop:
      while( NumberGondzioCorrections < maximum_correctors
            && NumberSmallCorrectors < max_additional_correctors
            && (PIPSisLT(alpha_pri, 1.0) || PIPSisLT(alpha_dual, 1.0)) )
      {
         if( dynamic_corrector_schedule )
            adjustLimitGondzioCorrectors();
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
         if( small_corr )
            corrector_resid->project_r3(rmin, std::numeric_limits<double>::infinity());
         else
            corrector_resid->project_r3(rmin, rmax);

         // solve for corrector direction
         sys->solve(prob, iterate, corrector_resid, corrector_step); // corrector_step is now delta_m

         if( !bicgstab_converged && bigcstab_norm_res_rel * 1e2 * dnorm > resid->residualNorm() / dnorm )
         {
            PIPSdebugMessage("Gondzio corrector step computation in BiCGStab failed - break corrector loop");

            // try at least one small corrector
            if( !small_corr )
            {
               if( !small_corr_aggr )
               {
                  if( my_rank == 0 )
                     std::cout << "switching to small correctors aggressive" << std::endl;
                  small_corr_aggr = true;
               }
               if( my_rank == 0 )
               {
                  std::cout << "Switching to small corrector " << std::endl;
                  std::cout << "Alpha when switching: " << alpha_pri << " " << alpha_dual << std::endl;
               }
               small_corr = true;
               continue;
            }
            else
               // exit corrector loop if small correctors have already been tried
               break;
         }

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

            if( small_corr && !small_corr_aggr )
               NumberSmallCorrectors++;

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

            if( small_corr && !small_corr_aggr )
               NumberSmallCorrectors++;

            NumberGondzioCorrections++;
         }
         else if( alpha_pri_enhanced >= (1.0 + AcceptTol) * alpha_pri )
         {
            PIPSdebugMessage("primal better \n");

            step->saxpy_pd(corrector_step, weight_primal_candidate, 0.0);

            alpha_pri = alpha_pri_enhanced;

            if( small_corr && !small_corr_aggr )
               NumberSmallCorrectors++;

            NumberGondzioCorrections++;
         }
         else if( alpha_dual_enhanced >= (1.0 + AcceptTol) * alpha_dual )
         {
            PIPSdebugMessage("dual better \n");

            step->saxpy_pd(corrector_step, 0.0, weight_dual_candidate);

            alpha_dual = alpha_dual_enhanced;

            if( !small_corr_aggr )
               NumberSmallCorrectors++;

            NumberGondzioCorrections++;
         }
         else if( additional_correctors_small_comp_pairs && !small_corr && iter >= first_iter_small_correctors)
         {
            if( alpha_pri < max_alpha_small_correctors || alpha_dual < max_alpha_small_correctors )
            {
               // try and center small pairs
               small_corr = true;
               if( my_rank == 0 )
               {
                  std::cout << "Switching to small corrector " << std::endl;
                  std::cout << "Alpha when switching: " << alpha_pri << " " << alpha_dual << std::endl;
               }
            }
            else
               // exit Gondzio correction loop
               break;
         }
         else
         {
            // exit Gondzio correction loop
            break;
         }
      }

      // We've finally decided on a step direction, now calculate the
      // length using Mehrotra's heuristic.x
      finalStepLength_PD(iterate, step, alpha_pri, alpha_dual);

      // if we encountered numerical troubles while computing the step check enter a probing round
      if( numerical_troubles )
      {
         if( !precond_limit )
            precond_limit = decreasePreconditionerImpact(sys);

         const double mu_last = iterate->mu();
         const double resids_norm_last = resid->residualNorm();

         computeProbingStep_pd(temp_step, iterate, step, alpha_pri, alpha_dual);

         resid->calcresids(prob, temp_step, false);
         const double mu_probing = temp_step->mu();
         const double resids_norm_probing = resid->residualNorm();

         const double factor = computeStepFactorProbing(resids_norm_last, resids_norm_probing,
               mu_last, mu_probing);

         alpha_pri = factor * alpha_pri;
         alpha_dual = factor * alpha_dual;

         const double alpha_max = std::max(alpha_pri, alpha_dual);

         if( restartIterateBecauseOfPoorStep( pure_centering_step, precond_limit, alpha_max ) )
            continue;
      }

      // actually take the step and calculate the new mu
      iterate->saxpy_pd(step, alpha_pri, alpha_dual);
      mu = iterate->mu();

      pure_centering_step = false;
      numerical_troubles = false;

      stochFactory->iterateEnded();
   }

   resid->calcresids(prob, iterate);
   if( gOoqpPrintLevel >= 10 )
   {
      this->doMonitorPd(prob, iterate, resid, alpha_pri, alpha_dual, sigma, iter, mu, status_code, 1);
   }

   return status_code;
}

void GondzioStochLpSolver::computeProbingStep_pd(Variables* probing_step, const Variables* iterate, const Variables* step,
        double alpha_primal, double alpha_dual) const
{
   probing_step->copy(iterate);
   probing_step->saxpy_pd(step, alpha_primal, alpha_dual);
}


GondzioStochLpSolver::~GondzioStochLpSolver()
{
}
