

int main()
{

  //create the NLP formulation factory
  NlpGenFactory * nlp  = new NlpGenFactory( );

  NlpInput* in = new HandcodedSmallExampleInput(); // or can be the Ampl input class
  NlpDataFull    * prob  = (NlpDataFull * ) nlp->makeData(in);
    
  NlpGenVars     * vars  = (NlpGenVars * )  nlp->makeVariables( prob );
  NlpGenResids   * resid = (NlpGenResids* ) nlp->makeResiduals( prob );


  XXXNlpSolver* s            = new XXXNlpSolver( nlp, prob );
  s->monitorSelf();

  int result = s->solve(prob,vars, resid);

  // bla bla bla
    
  delete s;
  delete vars;  
  delete resid;
  delete prob;
  delete in;
  delete nlp;
};

