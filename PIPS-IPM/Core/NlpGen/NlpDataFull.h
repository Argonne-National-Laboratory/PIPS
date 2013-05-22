#ifndef NLPGENFACTFULL
#define NLPGENFACTFULL

class NlpDataFull : public Data {
 public:
  SymMatrix*     H; //Hessian
  GenMatrix*     Jeq; 
  GenMatrix*     Jineq;
  OoqpVector*    grad; // gradient

  //
  // Assuming that the NLP problem is formulated as QPs are in OOQP
  //

  //blx <= x <=  bux
  OoqpVector*    bux; //upper bounds
  OoqpVector*    ixupp; //indexes  showing which variable are bounded
  OoqpVector*    blx; //lower bounds
  OoqpVector*    ixlow;

  // bl <= c(x) <= bu
  OoqpVector*    bu;
  OoqpVector*    icupp; 
  OoqpVector*    bl;
  OoqpVector*    iclow;

  long long nx, my, mz;


  //
  // this is specific to NLP 
  //
  //evaluate all data ( H,Jeq,Jineq ) at the iterate specified by 'var'
  void evalData(NlpVars* var)
  {
    input->Hessian(H,var);
    input->JacEq(H,var);
    //....
  }

  //
  // QpGenData functionality below
  //
  
 protected:
  //
  // The input class (base class)
  // Derived classes will implement functionality for dealing with AMPL, SML, or PYOMO
  // When coding a new problem by hand, a class that inherits this will be needed to be implemented

  NlpInput* input;
};
#endif
