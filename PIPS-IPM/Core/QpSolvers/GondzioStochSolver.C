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
#include "StochOptions.h"

#include "sVars.h"
#include "sData.h"

#include <cstring>
#include <iostream>
#include <fstream>
#include <algorithm>

#include <cstdio>
#include <cassert>
#include <cmath>

#include "StochVector.h"
#include "mpi.h"
#include "QpGenVars.h"
#include "QpGenResiduals.h"
#include "QpGenLinsys.h"
#include "sLinsysRoot.h"

extern int gOoqpPrintLevel;
extern double g_iterNumber;
extern bool ipStartFound;

GondzioStochSolver::GondzioStochSolver( ProblemFormulation * opt, Data * prob, const Scaler* scaler )
  : GondzioSolver(opt, prob, scaler),
    n_linesearch_points( pips_options::getIntParameter("GONDZIO_STOCH_N_LINESEARCH")),
    dynamic_corrector_schedule( pips_options::getBoolParameter("GONDZIO_STOCH_USE_DYNAMIC_CORRECTOR_SCHEDULE") ),
    additional_correctors_small_comp_pairs( pips_options::getBoolParameter("GONDZIO_STOCH_ADDITIONAL_CORRECTORS_SMALL_VARS") ),
    max_additional_correctors( pips_options::getIntParameter("GONDZIO_STOCH_ADDITIONAL_CORRECTORS_MAX") ),
    first_iter_small_correctors( pips_options::getIntParameter("GONDZIO_STOCH_FIRST_ITER_SMALL_CORRECTORS") ),
    max_alpha_small_correctors( pips_options::getDoubleParameter("GONDZIO_STOCH_MAX_ALPHA_SMALL_CORRECTORS") ),
    NumberSmallCorrectors(0), bicgstab_skipped(false), bicgstab_converged(true), bigcstab_norm_res_rel(0.0), bicg_iterations(0),
    dynamic_bicg_tol(pips_options::getBoolParameter("OUTER_BICG_DYNAMIC_TOL"))
{
   assert(max_additional_correctors > 0);
   assert(first_iter_small_correctors >= 0);
   assert(0 < max_alpha_small_correctors && max_alpha_small_correctors < 1);
   assert(n_linesearch_points > 0);

   if( pips_options::getBoolParameter("GONDZIO_STOCH_ADAPTIVE_LINESEARCH") )
   {
      const int size = PIPS_MPIgetSize();

      if( size > 1)
         this->n_linesearch_points =
               std::min(unsigned(size) + this->n_linesearch_points, max_linesearch_points);
   }

   // the two StepFactor constants set targets for increase in step
   // length for each corrector
   StepFactor0 = 0.3;
   StepFactor1 = 1.5;

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
   const int my_rank = PIPS_MPIgetRank(MPI_COMM_WORLD);

   int done;
   double mu, muaff;
   double alpha_target, alpha_enhanced;
   int status_code;
   double alpha = 1, sigma = 1;
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
   done = 0;
   mu = iterate->mu();

   do
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
         this->doMonitor(prob, iterate, resid, alpha, sigma, iter, mu,
               status_code, 0);
      }

      // *** Predictor step ***
      if( !pure_centering_step )
      {
         resid->set_r3_xz_alpha(iterate, 0.0);

         sys->factor(prob, iterate);
         sys->solve(prob, iterate, resid, step);
         step->negate();

         if( !bicgstab_converged && bigcstab_norm_res_rel * 1e2 * dnorm >= resid->residualNorm() )
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

      alpha = iterate->stepbound(step);

      // calculate centering parameter
      muaff = iterate->mustep(step, alpha);

      assert( !PIPSisZero(mu) );
      sigma = pow(muaff / mu, tsig);

      if( gOoqpPrintLevel >= 10 )
      {
         this->doMonitor(prob, iterate, resid, alpha, sigma, iter, mu,
               status_code, 2);
      }

      g_iterNumber += 1.0;

      // *** Corrector step ***
      corrector_resid->clear_r1r2();

      // form right hand side of linear system:
      corrector_resid->set_r3_xz_alpha(step, -sigma * mu);

      sys->solve(prob, iterate, corrector_resid, corrector_step);
      corrector_step->negate();

      // if the corrector fails refactorize
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
      double weight_candidate = -1.0;
      const double alpha_predictor = alpha;

      calculateAlphaWeightCandidate(iterate, step, corrector_step, alpha_predictor, alpha, weight_candidate);

      assert(weight_candidate >= 0.0 && weight_candidate <= 1.0);

      step->saxpy(corrector_step, weight_candidate);

      // prepare for Gondzio corrector loop: zero out the
      // corrector_resid structure:
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
            && PIPSisLT(alpha, 1.0) )
      {
         if( dynamic_corrector_schedule )
            adjustLimitGondzioCorrectors();
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
         if( small_corr )
            corrector_resid->project_r3(rmin, std::numeric_limits<double>::infinity());
         else
            corrector_resid->project_r3(rmin, rmax);

         // solve for corrector direction
         sys->solve(prob, iterate, corrector_resid, corrector_step);	// corrector_step is now delta_m

         /* if a normal corrector did not converge - discard it, try a small corr one and set small correctors to aggressive */
         if( !bicgstab_converged && bigcstab_norm_res_rel * 1e2 * dnorm > resid->residualNorm() )
         {
            PIPSdebugMessage("Gondzio corrector step computation in BiCGStab failed - break corrector loop");

            // try small correctors to improve centering and numerical stability
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
                  std::cout << "Alpha when switching: " << alpha << std::endl;
               }
               small_corr = true;
               continue;
            }
            // exit corrector loop if small correctors have already been tried
            else
               break;
         }

         // calculate weighted predictor-corrector step
         calculateAlphaWeightCandidate(iterate, step, corrector_step, alpha_target, alpha_enhanced, weight_candidate);

         // if the enhanced step length is actually 1, make it official
         // and stop correcting
         if( PIPSisEQ(alpha_enhanced, 1.0) )
         {
            step->saxpy(corrector_step, weight_candidate);
            alpha = alpha_enhanced;

            if( small_corr && !small_corr_aggr )
               NumberSmallCorrectors++;

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

            if( small_corr && !small_corr_aggr )
               NumberSmallCorrectors++;

            NumberGondzioCorrections++;
         }
         /* if not done yet because correctors were not good enough - try a small corrector if enabled */
         else if( additional_correctors_small_comp_pairs && !small_corr && iter >= first_iter_small_correctors )
         {
            if( alpha < max_alpha_small_correctors )
            {
               small_corr = true;
               if( my_rank == 0 )
               {
                  std::cout << "Switching to small corrector " << std::endl;
                  std::cout << "Alpha when switching: " << alpha << std::endl;
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
      alpha = finalStepLength(iterate, step);

      // alternatively, just use a crude step scaling factor.
      // alpha = 0.995 * iterate->stepbound( step );

      // if we encountered numerical troubles while computing the step check enter a probing round
      if( numerical_troubles )
      {
         if( !precond_limit )
            precond_limit = decreasePreconditionerImpact(sys);

         const double mu_last = iterate->mu();
         const double resids_norm_last = resid->residualNorm();

         computeProbingStep(temp_step, iterate, step, alpha);

         resid->calcresids(prob, temp_step, false);
         const double mu_probing = temp_step->mu();
         const double resids_norm_probing = resid->residualNorm();

         const double factor = computeStepFactorProbing(resids_norm_last, resids_norm_probing,
               mu_last, mu_probing);

         alpha = factor * alpha;

         if( restartIterateBecauseOfPoorStep( pure_centering_step, precond_limit, alpha ) )
            continue;
      }

      // actually take the step (at last!) and calculate the new mu
      iterate->saxpy(step, alpha);
      mu = iterate->mu();

      pure_centering_step = false;
      numerical_troubles = false;

      stochFactory->iterateEnded();
   }
   while( !done );

   resid->calcresids(prob, iterate);
   if( gOoqpPrintLevel >= 10 )
   {
      this->doMonitor(prob, iterate, resid, alpha, sigma, iter, mu, status_code, 1);
   }

   return status_code;
}

void GondzioStochSolver::registerBiCGStabOvserver(LinearSystem* sys)
{
   /* every linsys handed to the GondzioStoch should be observable */
   assert( dynamic_cast<Subject*>(sys) );
   setSubject( dynamic_cast<Subject*>(sys) );
}

void GondzioStochSolver::notifyFromSubject()
{
   const Subject& subj = *getSubject();

   bicgstab_skipped = subj.getBoolValue("BICG_SKIPPED");
   if( !bicgstab_skipped )
      bicgstab_converged = subj.getBoolValue("BICG_CONVERGED");
   else
      bicgstab_converged = true;
   bigcstab_norm_res_rel = subj.getDoubleValue("BICG_RELRESNORM");
   bicg_iterations = subj.getIntValue("BICG_NITERATIONS");
   if( !bicgstab_converged )
      PIPSdebugMessage("BiGCStab had troubles converging\n");
}

void GondzioStochSolver::setBiCGStabTol(int iteration) const
{
   if( !dynamic_bicg_tol )
      return;

   assert( iteration >= -1);

   if( iteration == -1 )
      pips_options::setDoubleParameter("OUTER_BICG_TOL", 1e-10);
   else if( iteration <= 4 )
      pips_options::setDoubleParameter("OUTER_BICG_TOL", 1e-8);
   else if( iteration <= 8 )
      pips_options::setDoubleParameter("OUTER_BICG_TOL", 1e-9);
   else
      pips_options::setDoubleParameter("OUTER_BICG_TOL", 1e-10);
}

void GondzioStochSolver::adjustLimitGondzioCorrectors()
{
   assert( bicg_iterations >= 0 );
   if( dynamic_corrector_schedule )
   {
      if( bicgstab_skipped )
         maximum_correctors = 5;
      else if( bicg_iterations < 2 )
         maximum_correctors = 4;
      else if( bicg_iterations <= 15 )
         maximum_correctors = 3;
      else if( bicg_iterations < 25 )
         maximum_correctors = 2;
      else if( bicg_iterations > 35 )
         maximum_correctors = 1;
   }
}

bool GondzioStochSolver::decreasePreconditionerImpact(LinearSystem* sys) const
{
   bool success = false;
   dynamic_cast<sLinsysRoot*>(sys)->precondSC.decreaseDiagDomBound(success);
   if( !success )
   {
      if( PIPS_MPIgetRank() == 0 )
         std::cout << "Cannot increase precision in preconditioner anymore" << std::endl;
   }
   return success;
}

void GondzioStochSolver::computeProbingStep(Variables* probing_step, const Variables* iterate, const Variables* step,
      double alpha) const
{
   probing_step->copy(iterate);
   probing_step->saxpy(step, alpha);
}

double GondzioStochSolver::computeStepFactorProbing(double resids_norm_last, double resids_norm_probing,
      double mu_last, double mu_probing) const
{
   double factor = 1.0;
   const double limit_resids = std::max( artol * dnorm, resids_norm_last );

   if( resids_norm_probing > limit_resids )
   {
      const double resids_diff = resids_norm_probing - resids_norm_last;
      const double resids_max_change = limit_resids - resids_norm_last;
      assert( resids_diff > 0 ); assert( resids_max_change > 0 );
      assert( resids_max_change < resids_diff );

      factor = std::min(factor, resids_max_change / resids_diff * 0.9995 );
   }

   if( mu_probing > 10 * mu_last )
   {
      const double mu_diff = mu_probing - mu_last;
      const double mu_max_change = 10 * mu_last - mu_probing;
      assert( mu_diff > 0 ); assert( mu_max_change > 0 );
      assert( mu_max_change < mu_diff );

      factor = std::min(factor, mu_max_change / mu_diff * 0.9995 );
   }
   return factor;
}

bool GondzioStochSolver::restartIterateBecauseOfPoorStep( bool& pure_centering_step,
      bool precond_limit, double alpha_max) const
{
   const int my_rank = PIPS_MPIgetRank();

   if( !pure_centering_step && alpha_max < mutol * 1e-2 )
   {
      if( my_rank == 0 )
         std::cout << "poor step computed - trying pure centering step" << std::endl;
      pure_centering_step = true;
      return true;
   }
   else if( alpha_max < mutol * 1e-2 && pure_centering_step && !precond_limit )
   {
      if( my_rank == 0 )
         std::cout << "refactorization" << std::endl;
      return true;
   }
   return false;
}


GondzioStochSolver::~GondzioStochSolver()
{
   delete temp_step;
}
