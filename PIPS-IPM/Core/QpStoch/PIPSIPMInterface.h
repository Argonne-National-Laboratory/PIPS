#ifndef QPSTOCHINTERFACE
#define QPSTOCHINTERFACE

//#include "BALPSolverInterface.hpp"

#include "stochasticInput.hpp"
#include "sDriver.h"

class PIPSIPMInterface //: public BALPSolverInterface<PIPSIPMInterface>
{
 public:
  PIPSIPMInterface(stochasticInput &in);
  ~PIPSIPMInterface();

  void loadData();
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
};

#endif
