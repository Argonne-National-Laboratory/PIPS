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
using namespace std;

#include <cstdio>
#include <cassert>
#include <cmath>

#include "StochVector.h"
#include "mpi.h"
#include "QpGenVars.h"
#include "QpGenResiduals.h"
#include "QpGenLinsys.h"

extern int gOoqpPrintLevel;
extern double g_iterNumber;
extern bool ipStartFound;


GondzioStochSolver::GondzioStochSolver( ProblemFormulation * opt, Data * prob, unsigned int n_linesearch_points,
      bool adaptive_linesearch )
  : GondzioSolver(opt, prob), n_linesearch_points(n_linesearch_points),
    additional_correctors_small_comp_pairs( pips_options::getBoolParameter("IP_GONDZIO_ADDITIONAL_CORRECTORS_SMALL_VARS") ),
    max_additional_correctors( pips_options::getIntParameter("IP_GONDZIO_ADDITIONAL_CORRECTORS_MAX") ),
    first_iter_small_correctors( pips_options::getIntParameter("IP_GONDZIO_FIRST_ITER_SMALL_CORRECTORS") ),
    max_alpha_small_correctors( pips_options::getDoubleParameter("IP_GONDZIO_MAX_ALPHA_SMALL_CORRECTORS") ),
    NumberSmallCorrectors(0), bicgstab_converged(true), bigcstab_norm_res_rel(0.0)
{
   assert(max_additional_correctors > 0);
   assert(first_iter_small_correctors >= 0);
   assert(0 < max_alpha_small_correctors && max_alpha_small_correctors < 1);
   assert(n_linesearch_points > 0);

   if( adaptive_linesearch )
   {
      const int size = PIPS_MPIgetSize(MPI_COMM_WORLD);

      if( size > 1)
         this->n_linesearch_points =
               std::min(unsigned(size) + this->n_linesearch_points, max_linesearch_points);
   }

   // the two StepFactor constants set targets for increase in step
   // length for each corrector
   StepFactor0 = 0.3;
   StepFactor1 = 1.5;

   // todo the parameters should be read in Solver.c
   if( pips_options::getBoolParameter("IP_STEPLENGTH_CONSERVATIVE") )
   {
      steplength_factor = 0.99;
      gamma_f = 0.95;
   }
   else
   {
      steplength_factor = 0.99999999;
      gamma_f = 0.99;
   }
   gamma_a = 1.0 / (1.0 - gamma_f);

   if( pips_options::getBoolParameter("IP_ACCURACY_REDUCED")  )
   {
	  artol = 1.e-3;
	  mutol = 1.e-5;
   }
   else
   {
	  artol = 1.e-4;
      mutol = 1.e-6;
   }

   if( pips_options::getBoolParameter("IP_PRINT_TIMESTAMP") )
   {
      printTimeStamp = true;
      startTime = MPI_Wtime();
   }

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
   bool do_small_correctors_aggressively = false;
   bool refactorized = false;

   dnorm = prob->datanorm();

   // initialization of (x,y,z) and factorization routine.
   sys = factory->makeLinsys(prob);

   // register as observer for the BiCGStab solves
   registerBiCGStabOvserver(sys);

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

      // if we are not in a refactorization - check convergence of bicgstab
      if( !refactorized &&
            !bicgstab_converged && bigcstab_norm_res_rel * 1e-2 > resid->residualNorm() / dnorm )
      {
         PIPSdebugMessage("Affine step computation in BiCGStab failed");
         // TODO : improve accuracy in preconditioner

         // go on, discard the affine step and try a pure centering step
         step->setToZero();
         do_small_correctors_aggressively = true;
      }

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

      // if we are not in a refactorization - check convergence of bicgstab
      if( !refactorized
            && !bicgstab_converged
            && bigcstab_norm_res_rel * 1e-2 > resid->residualNorm() / dnorm )
      {
         PIPSdebugMessage("1st corrector step computation in BiCGStab failed");
         do_small_correctors_aggressively = true;

         // TODO : use preconditioner more conservatively
         if( my_rank == 0 )
            std::cout << "refactorizing since lin solves failed" << std::endl;
         refactorized = true;
         continue;
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
      bool do_small_pairs_correction = do_small_correctors_aggressively;

      // enter the Gondzio correction loop:
      while( NumberGondzioCorrections < maximum_correctors
            && NumberSmallCorrectors < max_additional_correctors
            && PIPSisLT(alpha, 1.0) )
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
         if( do_small_pairs_correction )
            corrector_resid->project_r3(rmin, std::numeric_limits<double>::infinity());
         else
            corrector_resid->project_r3(rmin, rmax);

         // solve for corrector direction
         sys->solve(prob, iterate, corrector_resid, corrector_step);	// corrector_step is now delta_m

         if( !bicgstab_converged && bigcstab_norm_res_rel * 1e-2 > resid->residualNorm() / dnorm )
         {
            PIPSdebugMessage("Gondzio corrector step computation in BiCGStab failed - break corrector loop");

            if( !do_small_pairs_correction )
            {
               do_small_correctors_aggressively = true;
               do_small_pairs_correction = true;
               continue;
            }
            else
               // exit corrector loop if small correctors have already been tried
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

            if( do_small_pairs_correction )
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

            if( do_small_pairs_correction )
               NumberSmallCorrectors++;

            NumberGondzioCorrections++;
         }
         else if( additional_correctors_small_comp_pairs && !do_small_pairs_correction && iter >= first_iter_small_correctors )
         {
            if( alpha < max_alpha_small_correctors )
            {
               do_small_pairs_correction = true;
               if( my_rank == 0 )
               {
                  std::cout << "Small corrector " << std::endl;
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

      // actually take the step (at last!) and calculate the new mu
      iterate->saxpy(step, alpha);
      mu = iterate->mu();

      if( 0 == my_rank )
         std::cout << "final alpha: " << alpha << " mu: " << mu <<   std::endl;

      refactorized = false;
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

   bicgstab_converged = subj.getBoolValue("BICG_CONVERGED");
   bigcstab_norm_res_rel = subj.getDoubleValue("BICG_RELRESNORM");

   if( !bicgstab_converged )
      PIPSdebugMessage("BiGCStab had troubles converging\n");
}

GondzioStochSolver::~GondzioStochSolver()
{
   delete temp_step;
}
