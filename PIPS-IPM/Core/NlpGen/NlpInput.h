#ifndef NLPINPUT
#define NLPINPUT

class NlpVars;

class NlpInputFull
{
  //
  // returns the Hessian evaluated at the point specified by 'pt'
  //
  void Hessian(SymMatrix* H, NlpVars* var);
  //
  // Jacobian of the equalities
  //
  void JacEq(GenMatrix* H, NlpVars* var);
  //
  // Equality constraints function
  //
  void conEq(OoqpVector* coneq, NlpVars* var);

  // the other methods goes here: gradient, JacIneq, etc.
}

#endif
