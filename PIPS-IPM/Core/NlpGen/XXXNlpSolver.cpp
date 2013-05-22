XXXNlpSolver::XXXNlpSolver(NlpGenFactory * opt, NlpDataFull * prob )
{
  factory = opt;
  step = factory->makeVariables( prob );

  //other initializations
}
XXXNlpSolver::~XXXNlpSolver()
{
  delete step;
}

// An example of the solve method
//
// 
int XXXNlpSolver::solve(Data *prob_, Variables *iterate_, Residuals * resid_ )
{
  NlpDataFull* prob = dynamic_cast<NlpDataFull*>(prob_);
  NlpVars* iterate = dynamic_cast<NlpVars*>(prob_);
  NlpResids* resid = dynamic_cast<NlpResids*>(resid_);


  //this is copy-pasted from Mehrotra's algorithm; use it for example for convex QPs
  //there is only one addition to the code, see the end of the 'while' loop
  int done;
  double mu, alpha = 1, sigma = 1, muaff;
  int status_code;

  gmu = 1000;
  //  grnorm = 1000;
  dnorm = prob->datanorm();

  // initialization of (x,y,z) and factorization routine.
  sys = factory->makeLinsys( prob );
  this->start( factory, iterate, prob, resid, step );

  iter = 0;
  done = 0;
  mu = iterate->mu();
  gmu = mu;

  do
    {
      iter ++;

      // evaluate residuals and update algorithm status:
      resid->calcresids(prob, iterate);
      
      // termination test:
      status_code = this->doStatus( prob, iterate, resid, iter, mu, 0 );
      if( status_code != NOT_FINISHED ) break;
      if( gOoqpPrintLevel >= 10 ) {
	this->doMonitor( prob, iterate, resid,
			 alpha, sigma, iter, mu, status_code, 0 );
      }

      // *** Predictor step ***

      resid->set_r3_xz_alpha(iterate, 0.0 );

      sys->factor(prob, iterate);
      sys->solve(prob, iterate, resid, step);

      step->negate();

      alpha = iterate->stepbound(step);

      // calculate centering parameter 
      muaff = iterate->mustep(step, alpha);
      sigma = pow(muaff/mu, tsig);

      assert(sigma<=1);

      // *** Corrector step ***

      // form right hand side of linear system:
      resid->add_r3_xz_alpha( step, -sigma*mu );

      sys->solve(prob, iterate, resid, step);
      step->negate();

 
      // We've finally decided on a step direction, now calculate the
      // length using Mehrotra's heuristic.
      alpha = finalStepLength(iterate, step);

      // alternatively, just use a crude step scaling factor.
      //alpha = 0.995 * iterate->stepbound( step );

      // actually take the step and calculate the new mu
      iterate->saxpy(step, alpha);
      mu = iterate->mu();
      gmu = mu;

      //
      //specific to NLP: evaluate Hessian, Jacobians and all that
      //
      prob->eval(iterate);

    } while(!done);
  
  resid->calcresids(prob,iterate);
  if( gOoqpPrintLevel >= 10 ) {
    this->doMonitor( prob, iterate, resid,
		     alpha, sigma, iter, mu, status_code, 1 );
  }

  return status_code;
}
