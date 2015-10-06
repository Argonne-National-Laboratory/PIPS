/* PIPS-IPM                                                           *
 * Author:  Cosmin G. Petra                                           *
 * (C) 2012 Argonne National Laboratory. See Copyright Notification.  */

#ifndef NLPSTOCHINTERFACE_CALLBACKS
#define NLPSTOCHINTERFACE_CALLBACKS

#include "stochasticInput.hpp"
#include "sDriver.h"

/**
 * This class is an interface to PIPS-IPM which takes the problem from 
 * stochasticInput and specify it to the solver using C-callbacks.
 * 
 * Implemented for backward compatibility. Callbacks are now obsolete,
 * PIPSIPMInterface should be used instead.
 * 
 * TO DO: implement get..Solution methods.
 */

class sInterfaceCallbacks
{
 public:
  sInterfaceCallbacks(stochasticInput &in);
  ~sInterfaceCallbacks();

  void go();

  double getObjective() const;
  //solverState getStatus();

  void setPrimalTolerance(double val);
  void setDualTolerance(double val);

  std::vector<double> getFirstStagePrimalColSolution() const;
  std::vector<double> getSecondStagePrimalColSolution(int scen) const;
  std::vector<double> getFirstStageDualColSolution() const;
  std::vector<double> getSecondStageDualColSolution(int scen) const;
  //more get methods to follow here

 protected:
  //count the number eq and ineq
  void getNum1stStgEqIneq(int& my, int &mz);
  void getNum2ndStgEqIneq(int scens, int& my, int &mz);
  
  StochInputTree* inputTree;

  void loadData();
};

#endif
