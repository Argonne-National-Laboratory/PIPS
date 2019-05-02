/*
 * GondzioStochLpSolver.C
 *
 *  Created on: 20.12.2017
 *      Author: Svenja Uslu
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

#include <fstream>
#include <string>
#include <sstream>
#include <limits>

// gmu is needed by MA57!
static double gmu;

// double grnorm;
extern int gOoqpPrintLevel;

static double g_iterNumber;


GondzioStochLpSolver::GondzioStochLpSolver( ProblemFormulation * opt, Data * prob, unsigned int n_linesearch_points, bool adaptive_linesearch)
  : GondzioStochSolver(opt, prob, n_linesearch_points, adaptive_linesearch)
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
   int done;
   double mu, muaff;
   double rmin, rmax;
   int status_code;
   double sigma = 1;
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

   // TODO write out bounds on x as well as x

   int myRank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

   stochFactory->iterateStarted();
   this->start(factory, iterate, prob, resid, step);
   stochFactory->iterateEnded();

   {
   std::ostringstream name;
   name << myRank << "check_singleton_rows_before_stoch.txt";
   std::ofstream ofs( name.str().c_str(), std::ofstream::out);
   const StochVector& ixupp = dynamic_cast<const StochVector&>(*dynamic_cast<sData*>(prob)->ixupp);
   const StochVector& ixlow  = dynamic_cast<const StochVector&>(*dynamic_cast<sData*>(prob)->ixlow);
   const StochVector& xupp = dynamic_cast<const StochVector&>(*dynamic_cast<sData*>(prob)->bux);
   const StochVector& xlow  = dynamic_cast<const StochVector&>(*dynamic_cast<sData*>(prob)->blx);

   for( int j = 0; j < dynamic_cast<SimpleVector&>(*xupp.vec).n; ++j )
   {
      ofs << "node: " << -1 << "\tcol: " << j << "\tx = " << dynamic_cast<SimpleVector&>(*dynamic_cast<StochVector&>(*dynamic_cast<sVars*>(iterate)->x).vec)[j] << " € [" << ((dynamic_cast<SimpleVector&>(*ixlow.vec)[j] == 1.0) ? dynamic_cast<SimpleVector&>(*xlow.vec)[j] : -std::numeric_limits<double>::infinity()) << ", "
            << ((dynamic_cast<SimpleVector&>(*ixupp.vec)[j] == 1.0) ? dynamic_cast<SimpleVector&>(*xupp.vec)[j] : -std::numeric_limits<double>::infinity()) << "]" << std::endl;
   }
   for( size_t i = 0; i < xupp.children.size(); ++i )
   {
      if( xupp.children[i] )
      {
         for( int j = 0; j < dynamic_cast<SimpleVector&>(*xupp.children.at(i)->vec).n; ++j )
         {
            ofs << "node: " << i << "\tcol: " << j << "\tx = " << dynamic_cast<SimpleVector&>(*dynamic_cast<StochVector&>(*dynamic_cast<sVars*>(iterate)->x).children[i]->vec)[j] << " € [" << ( ( dynamic_cast<SimpleVector&>(*ixlow.children[i]->vec)[j] == 1.0 ) ? dynamic_cast<SimpleVector&>(*xlow.children[i]->vec)[j] :
                  -std::numeric_limits<double>::infinity() ) << ", " << ( ( dynamic_cast<SimpleVector&>(*ixupp.children[i]->vec)[j] == 1.0 ) ? dynamic_cast<SimpleVector&>(*xupp.children[i]->vec)[j] :
                        -std::numeric_limits<double>::infinity() ) << "]" << std::endl;
         }
      }
   }
   ofs.close();
   MPI_Barrier(MPI_COMM_WORLD );
   }

   iter = 0;
   NumberGondzioCorrections = 0;
   done = 0;
   mu = iterate->mu();
   gmu = mu;

   do
   {
      if(iter > 0)
      {
         std::ostringstream name;
         name << myRank << "check_singleton_rows_before_stoch_it_" << iter << ".txt";
         std::ofstream ofs( name.str().c_str(), std::ofstream::out);

         StochVector& ixupp = dynamic_cast<StochVector&>(*dynamic_cast<sVars*>(iterate)->ixupp);
         StochVector& ixlow = dynamic_cast<StochVector&>(*dynamic_cast<sVars*>(iterate)->ixlow);
         StochVector& xupp= dynamic_cast<StochVector&>(*dynamic_cast<sData*>(prob)->bux);
         StochVector& xlow = dynamic_cast<StochVector&>(*dynamic_cast<sData*>(prob)->blx);

         for( int j = 0; j < dynamic_cast<SimpleVector&>(*xupp.vec).n; ++j )
         {
            ofs << "node: " << -1 << "\tcol: " << j << "\tx = " << dynamic_cast<SimpleVector&>(*dynamic_cast<StochVector&>(*dynamic_cast<sVars*>(iterate)->x).vec)[j] << " € [" << ((dynamic_cast<SimpleVector&>(*ixlow.vec)[j] == 1.0) ? dynamic_cast<SimpleVector&>(*xlow.vec)[j] : -std::numeric_limits<double>::infinity()) << ", "
                  << ((dynamic_cast<SimpleVector&>(*ixupp.vec)[j] == 1.0) ? dynamic_cast<SimpleVector&>(*xupp.vec)[j] : -std::numeric_limits<double>::infinity()) << "]" << std::endl;
         }
         for( size_t i = 0; i < xupp.children.size(); ++i )
         {
            if( xupp.children[i] )
            {
               for( int j = 0; j < dynamic_cast<SimpleVector&>(*xupp.children.at(i)->vec).n; ++j )
               {
                  ofs << "node: " << i << "\tcol: " << j << "\tx = " << dynamic_cast<SimpleVector&>(*dynamic_cast<StochVector&>(*dynamic_cast<sVars*>(iterate)->x).children[i]->vec)[j] << " € [" << ( ( dynamic_cast<SimpleVector&>(*ixlow.children[i]->vec)[j] == 1.0 ) ? dynamic_cast<SimpleVector&>(*xlow.children[i]->vec)[j] :
                        -std::numeric_limits<double>::infinity() ) << ", " << ( ( dynamic_cast<SimpleVector&>(*ixupp.children[i]->vec)[j] == 1.0 ) ? dynamic_cast<SimpleVector&>(*xupp.children[i]->vec)[j] :
                              -std::numeric_limits<double>::infinity() ) << "]" << std::endl;
               }
            }
         }
         ofs.close();
         MPI_Barrier(MPI_COMM_WORLD);
      }

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
         this->doMonitorPd(prob, iterate, resid, alpha_pri, alpha_dual, sigma, iter, mu,
               status_code, 0);
      }
      // *** Predictor step ***
      StochVector* copy_iterate_vars = dynamic_cast<StochVector&>(*dynamic_cast<sVars*>(iterate)->x).cloneFull();
      double obj = copy_iterate_vars->dotProductWith(dynamic_cast<StochVector&>(*dynamic_cast<sData*>(prob)->g));
      double onenorm = dynamic_cast<StochVector&>(*dynamic_cast<sVars*>(iterate)->x).onenorm();
      double twonorm = dynamic_cast<StochVector&>(*dynamic_cast<sVars*>(iterate)->x).twonorm();
      double infnorm = dynamic_cast<StochVector&>(*dynamic_cast<sVars*>(iterate)->x).infnorm();
      if(myRank == 0)
         std::cout << std::endl << "GondzioStochLp\tIterate\tonenorm: " << onenorm << "\ttwonorm: " << twonorm << "\tinfnorm: " << infnorm << "\tobj: " << obj << std::endl;

      resid->set_r3_xz_alpha(iterate, 0.0);

      sys->factor(prob, iterate);
      sys->solve(prob, iterate, resid, step);
      step->negate();

      StochVector* copy_x_vars = dynamic_cast<StochVector&>(*dynamic_cast<sVars*>(step)->x).cloneFull();
      obj = copy_x_vars->dotProductWith(dynamic_cast<StochVector&>(*dynamic_cast<sData*>(prob)->g));
      onenorm = dynamic_cast<StochVector&>(*dynamic_cast<sVars*>(step)->x).onenorm();
      twonorm = dynamic_cast<StochVector&>(*dynamic_cast<sVars*>(step)->x).twonorm();
      infnorm = dynamic_cast<StochVector&>(*dynamic_cast<sVars*>(step)->x).infnorm();
      if(myRank == 0)
         std::cout << std::endl << "onenorm: " << onenorm << "\ttwonorm: " << twonorm << "\tinfnorm: " << infnorm << "\tobj: " << obj << std::endl;

      iterate->stepbound_pd(step, alpha_pri, alpha_dual);

      // calculate centering parameter
      muaff = iterate->mustep_pd(step, alpha_pri, alpha_dual);

      sigma = pow(muaff / mu, tsig);

      if( gOoqpPrintLevel >= 10 )
      {
         this->doMonitorPd(prob, iterate, resid, alpha_pri, alpha_dual, sigma, iter, mu,
               status_code, 2);
      }

      g_iterNumber+=0.5;

      // *** Corrector step ***

      corrector_resid->clear_r1r2();

      // form right hand side of linear system:
      corrector_resid->set_r3_xz_alpha(step, -sigma * mu);

      sys->solve(prob, iterate, corrector_resid, corrector_step);
      corrector_step->negate();

      onenorm = dynamic_cast<StochVector&>(*dynamic_cast<sVars*>(corrector_step)->x).onenorm();
      twonorm = dynamic_cast<StochVector&>(*dynamic_cast<sVars*>(corrector_step)->x).twonorm();
      infnorm = dynamic_cast<StochVector&>(*dynamic_cast<sVars*>(corrector_step)->x).infnorm();
      copy_x_vars = dynamic_cast<StochVector&>(*dynamic_cast<sVars*>(corrector_step)->x).cloneFull();
      obj = copy_x_vars->dotProductWith(dynamic_cast<StochVector&>(*dynamic_cast<sData*>(prob)->g));
      if(myRank == 0)
         std::cout << std::endl << "GondzioStochLP\tonenorm: " << onenorm << "\ttwonorm: " << twonorm << "\tinfnorm: " << infnorm << "\tobj: " << obj << std::endl;

      // calculate weighted predictor-corrector step
      double weight_primal_candidate, weight_dual_candidate = -1.0;

      calculateAlphaPDWeightCandidate(iterate, step, corrector_step, alpha_pri, alpha_dual,
    		  alpha_pri, alpha_dual, weight_primal_candidate, weight_dual_candidate);

      assert(weight_primal_candidate >= 0.0 && weight_primal_candidate <= 1.0);
      assert(weight_dual_candidate >= 0.0 && weight_dual_candidate <= 1.0);

      step->saxpy_pd(corrector_step, weight_primal_candidate, weight_dual_candidate);

      // prepare for Gondzio corrector loop: zero out the
      // corrector_resid structure:
      corrector_resid->clear_r1r2();

      // calculate the target box:
      rmin = sigma * mu * beta_min;
      rmax = sigma * mu * beta_max;

      NumberGondzioCorrections = 0;

      // enter the Gondzio correction loop:
      while( NumberGondzioCorrections < maximum_correctors && (PIPSisLT(alpha_pri, 1.0) || PIPSisLT(alpha_dual, 1.0)) )
      {
         PIPSdebugMessage("corrector loop: %d alpha_pri: %f alpha_dual %f \n", NumberGondzioCorrections, alpha_pri, alpha_dual);

         // copy current variables into corrector_step
         corrector_step->copy(iterate);

         // calculate target steplength
         alpha_pri_target = StepFactor1 * alpha_pri + StepFactor0;
         alpha_dual_target = StepFactor1 * alpha_dual + StepFactor0;
         if( alpha_pri_target > 1.0 )
            alpha_pri_target = 1.0;
         if( alpha_dual_target > 1.0 )
            alpha_dual_target = 1.0;

         // add a step of this length to corrector_step
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
        		 alpha_pri_enhanced, alpha_dual_enhanced, weight_primal_candidate, weight_dual_candidate);
#if 0
         // if the enhanced step length is actually 1, make it official
         // and stop correcting
         if( alpha_pri_enhanced == 1.0 || alpha_dual_enhanced == 1.0)
         {
            step->saxpy_pd(corrector_step, weight_primal_candidate, weight_dual_candidate);

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
            step->saxpy_pd(corrector_step, weight_primal_candidate, weight_dual_candidate);

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
#else
      // if the enhanced step length is actually 1, make it official
      // and stop correcting
      if( PIPSisEQ(alpha_pri_enhanced, 1.0) && PIPSisEQ(alpha_dual_enhanced, 1.0) )
      {
         PIPSdebugMessage("both one \n");

         step->saxpy_pd(corrector_step, weight_primal_candidate, weight_dual_candidate);

         alpha_pri = alpha_pri_enhanced;
         alpha_dual = alpha_dual_enhanced;
         NumberGondzioCorrections++;

         // exit Gondzio correction loop
         break;
      }
      else if( alpha_pri_enhanced >= (1.0 + AcceptTol) * alpha_pri && alpha_dual_enhanced >= (1.0 + AcceptTol) * alpha_dual)
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
   }
#endif

      // We've finally decided on a step direction, now calculate the
      // length using Mehrotra's heuristic.x
      finalStepLength_PD(iterate, step, alpha_pri, alpha_dual);

      // actually take the step (at last!) and calculate the new mu

      iterate->saxpy_pd(step, alpha_pri, alpha_dual);
      mu = iterate->mu();
      gmu = mu;

      stochFactory->iterateEnded();
   }
   while( !done );

   resid->calcresids(prob, iterate);
   if( gOoqpPrintLevel >= 10 )
   {
      this->doMonitorPd(prob, iterate, resid, alpha_pri, alpha_dual, sigma, iter, mu, status_code, 1);
   }

   // print the results, if you really want to..
   // iterate->print();

   return status_code;
}


GondzioStochLpSolver::~GondzioStochLpSolver()
{
}



