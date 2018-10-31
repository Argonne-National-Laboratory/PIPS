/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

/* 2015. Modified by Nai-Yuan Chiang for NLP*/

#ifndef NLPGENRESIDUALS
#define NLPGENRESIDUALS

#include "Residuals.h"
#include "OoqpVectorHandle.h"


class NlpGen;
class NlpGenData;
class Variables;
class NlpGenVars;
class LinearAlgebraPackage;

/** 
 * Residuals for the general NLP formulation 
 *
 * @ingroup NlpGen
 */

class NlpGenResiduals : public Residuals {
protected:
  long long nx, my, mz;

  double nxupp;
  OoqpVectorHandle ixupp;

  double nxlow;
  OoqpVectorHandle ixlow;

  double mcupp;
  OoqpVectorHandle icupp;

  double mclow;
  OoqpVectorHandle iclow;

  double perr;
  double derr;
  double comerr;
  double comerr_0;

  double Oriperr;
  double Oriderr;


  NlpGenResiduals() {  
  nrmrho = 0;
  nrmr =0;
  nrmrho0 = 0;
  psi = 0;
  mod0 = 0;
  mod = 0;
  dmod = 0;
  dWd = 0;
  thd = 0;
  dmodu = 0;
  comerr = 0.;
  };

public:
  // for barrier function
  OoqpVectorHandle rQ;	
  OoqpVectorHandle rA;
  OoqpVectorHandle rC;
  OoqpVectorHandle rz;
  OoqpVectorHandle rv;
  OoqpVectorHandle rw;
  OoqpVectorHandle rt;
  OoqpVectorHandle ru;
  OoqpVectorHandle rgamma;
  OoqpVectorHandle rphi;
  OoqpVectorHandle rlambda;
  OoqpVectorHandle rpi;

  OoqpVectorHandle rx;
  OoqpVectorHandle ry;


  OoqpVectorHandle rp_OriSys;   
  OoqpVectorHandle rd_OriSys;

  OoqpVectorHandle priWrk;	// for x
  OoqpVectorHandle dualWrk; // for y (eq cons)
  OoqpVectorHandle dualWrk_Z; // for z (ineq cons, here we do C_I(x)=s, Z is the dual for this constraint)
  OoqpVectorHandle priWrk_S; 

  OoqpVectorHandle Wd; 


  // residual norm for the full sys:  check victor's matlab
  double nrmrho;
  double nrmr;
  double nrmrho0;
  double psi;
  double mod0;
  double mod;
  double dmod;
  double dWd;
  double thd;
  double dmodu;

  double nrmc;
  double nrmrd;
  double gammaD;
  double res;
  

  NlpGenResiduals( LinearAlgebraPackage * la,
		  long long nx, long long my, long long mz,
		  OoqpVector * ixlow, OoqpVector * ixupp,
		  OoqpVector * iclow, OoqpVector * icupp );

  virtual void calcresids(Data *problem, Variables *vars);

  virtual void add_r3_xz_alpha(Variables *vars, double alpha);

  virtual void set_r3_xz_alpha(Variables *vars, double alpha);
  
  virtual void clear_r3();
  
  virtual void clear_r1r2();

  virtual void project_r3(double rmin, double rmax);

  virtual int  validNonZeroPattern();
  
  virtual ~NlpGenResiduals();
  

//  virtual void buildInitSysForY(Data *problem, Variables *vars);




  virtual double priErr();
  virtual double dualErr();

  virtual double OriPriErr();
  virtual double OriDualErr();

  virtual double comp_Err();
  virtual double comp_Err_0();

  
  virtual double getKKTRhsNorm_Primal(NlpGenData *prob, NlpGenVars *vars, 
  						const int normType=-1, const int isTrialStep=0);
  
  virtual double getKKTRhsNorm_Dual(NlpGenData *prob, NlpGenVars *vars, 
  						const int normType=-1, const int isTrialStep=0);
  
  virtual double getKKTError_Comp(NlpGenData *prob, NlpGenVars *vars, const double mu, 
  						const int normType=-1, const int isTrialStep=0);

  virtual void copyFrom(NlpGenResiduals *residual_in);

  //update lRhs for SOC: rA = AlphaStep*rA + trial_Iter->rA and rC
  virtual void updateSOCRhs(const double AlphaStep, NlpGenVars *vars_in,  NlpGenData *prob_in);

  virtual bool findSmallStep(NlpGenVars *vars, NlpGenVars *steps, const double tol_mach);

  virtual void addDampingTermToOneSidePart(const double DampingTerm);
  
  
protected:

  virtual void eval_kkt_res(NlpGenData *prob_in, NlpGenVars *vars_in);




  
};

#endif





