#ifndef PYOMOINPUT_HPP
#define PYOMOINPUT_HPP

#include "stochasticInput.hpp"
//#include "sml.h"

// loads from Pyomo nl file
class PyomoInput : public stochasticInput {
public:
  PyomoInput(const std::string &root_filename, int num_scens);
  ~PyomoInput();

  virtual int nScenarios() { return nScenarios_; }
  virtual int nFirstStageVars() { return nFirstStageVars_; }
  virtual int nFirstStageCons() { return nFirstStageCons_; }
  virtual int nSecondStageVars(int scen) { return nSecondStageVars_[scen]; }
  virtual int nSecondStageCons(int scen) { return nSecondStageCons_[scen]; }
  
  virtual std::vector<double> getFirstStageColLB(); 
  virtual std::vector<double> getFirstStageColUB(); 
  virtual std::vector<double> getFirstStageObj(); 
  virtual std::vector<std::string> getFirstStageColNames(); 
  virtual std::vector<double> getFirstStageRowLB(); 
  virtual std::vector<double> getFirstStageRowUB(); 
  virtual std::vector<std::string> getFirstStageRowNames(); 
  virtual bool isFirstStageColInteger(int col) { return false; }
  virtual bool isFirstStageColBinary(int col) { return false; }
  
  virtual std::vector<double> getSecondStageColLB(int scen);
  virtual std::vector<double> getSecondStageColUB(int scen);
  virtual std::vector<double> getSecondStageObj(int scen);
  virtual std::vector<std::string> getSecondStageColNames(int scen);
  virtual std::vector<double> getSecondStageRowUB(int scen);
  virtual std::vector<double> getSecondStageRowLB(int scen);
  virtual std::vector<std::string> getSecondStageRowNames(int scen);
  virtual double scenarioProbability(int scen) { return 1.0/nScenarios_; }
  virtual bool isSecondStageColInteger(int scen, int col) { return false; }
  virtual bool isSecondStageColBinary(int scen, int col) { return false; }
  
  virtual CoinPackedMatrix getFirstStageConstraints();
  virtual CoinPackedMatrix getSecondStageConstraints(int scen);
  virtual CoinPackedMatrix getLinkingConstraints(int scen); 
  
  
  
  virtual bool scenarioDimensionsEqual() { return dimsEqual; }
  virtual bool onlyBoundsVary() { return false; } // no easy way to check
  virtual bool allProbabilitiesEqual() { return true; }
  virtual bool continuousRecourse() { return true; }
  
  
protected:
  bool dimsEqual;

  std::string nl_root_filename;  

  int nScenarios_, nFirstStageVars_, nFirstStageCons_;
  std::vector<int> nSecondStageVars_, nSecondStageCons_;
  
  // no copying
  PyomoInput(const PyomoInput&);
  PyomoInput& operator=(const PyomoInput&);
  
};

#endif
