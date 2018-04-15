/* PIPS-IPM                                                           *
 * Author:  Cosmin G. Petra                                           *
 * (C) 2012 Argonne National Laboratory. See Copyright Notification.  */

#ifndef NLPGENSTOCHVARS
#define NLPGENSTOCHVARS

#include "NlpGenVars.h"
#include "OoqpVectorHandle.h"
#include <vector>

class NlpGen;
class NlpGenData;
class LinearAlgebraPackage;
class sTree;

class sVars : public NlpGenVars {
 public:
  sVars( sTree* tree,
	 OoqpVector * ixlow_in, OoqpVector * ixupp_in,
	 OoqpVector * iclow_in, OoqpVector * icupp_in);

  /** constructor in which the data and variable pointers are set to
   *  point to the given arguments 
   *
   *  Re: global nx,my,mz: these will overwrite nx, my, and mz members
   *  in the parent NlpGenVars class
   */
  sVars( sTree* tree,
	 OoqpVector * x_in, OoqpVector * s_in,
	 OoqpVector * y_in, OoqpVector * z_in,
	 OoqpVector * v_in, OoqpVector * gamma_in,
	 OoqpVector * w_in, OoqpVector * phi_in,
	 OoqpVector * t_in, OoqpVector * lambda_in,
	 OoqpVector * u_in, OoqpVector * pi_in,
	 OoqpVector * ixlow_in, long long nxlowGlobal,
	 OoqpVector * ixupp_in, long long nxuppGlobal,
	 OoqpVector * iclow_in, long long mclowGlobal,
	 OoqpVector * icupp_in, long long mcuppGlobal,
	 long long nxGlobal, long long myGlobal, long long mzGlobal);

  virtual ~sVars();
  
  virtual void sync();

//  long long n1stSlack;
//  long long n2ndSlack;  
//  long long nSlack;

  std::vector<sVars*> children;

 protected:
  sVars() {};
  sVars(const sVars& ) {};

  void createChildren();

  void AddChild(sVars* child);

  sTree* stochNode;
  
};
#endif

