/* PIPS
   Authors: Cosmin Petra
   See license and copyright information in the documentation */

#ifndef STOCH_RESOURCE_MON
#define STOCH_RESOURCE_MON

#include <vector>
#include <string>

//! not thread safe

enum stCommType { ctAllreduce=0, ctReduce, ctOther};

class NodeExecEntry
{
 public:
  NodeExecEntry(double localTime, double childrenTime);
  void clear();

  double tmLocal, tmChildren;

  void recLocal_start();
  void recLocal_stop();

  void recChildren_start();
  void recChildren_stop();

 protected:
  double tmOpStartLocal, tmOpStartChildren;
};

class NodeCommEntry
{
 public:
  NodeCommEntry(double time_, double size_, stCommType type_)
    : time(time_), size(size_), type(type_) {};
  double time, size;
  stCommType type;
};

using namespace std;

class StochNodeResourcesMonitor
{
 public:
  StochNodeResourcesMonitor();
  virtual ~StochNodeResourcesMonitor();

  virtual void reset();
  virtual void computeTotal();

  virtual void recFactTmLocal_start();
  virtual void recFactTmLocal_stop();
  virtual void recFactTmChildren_start();
  virtual void recFactTmChildren_stop();

  virtual void recLsolveTmLocal_start();
  virtual void recLsolveTmLocal_stop();
  virtual void recLsolveTmChildren_start();
  virtual void recLsolveTmChildren_stop();

  virtual void recDsolveTmLocal_start();
  virtual void recDsolveTmLocal_stop();
  virtual void recDsolveTmChildren_start();
  virtual void recDsolveTmChildren_stop();

  virtual void recLtsolveTmLocal_start();
  virtual void recLtsolveTmLocal_stop();
  virtual void recLtsolveTmChildren_start();
  virtual void recLtsolveTmChildren_stop();

  virtual void recSchurCom_start(double size, stCommType type);
  virtual void recSchurCom_stop(double size, stCommType type);

  virtual void recLsolveCom_start(double size, stCommType type);
  virtual void recLsolveCom_stop(double size, stCommType type);

  //record the dense matrix multiplication within Schur complement computation
  //the recorded time is part of the recFactTmChildren kept in eFact.tmChildren
  virtual void recSchurMultLocal_start();
  virtual void recSchurMultLocal_stop();
  virtual void recSchurMultChildren_start();
  virtual void recSchurMultChildren_stop();

  virtual void recReduceTmLocal_start();
  virtual void recReduceTmLocal_stop();

  virtual void recReduceScatterTmLocal_start();
  virtual void recReduceScatterTmLocal_stop();
  
  virtual void recBcastTmLocal_start();
  virtual void recBcastTmLocal_stop();
  

 public:

  NodeExecEntry eFact, eLsolve, eDsolve, eLtsolve;
  NodeExecEntry eTotal;
  NodeExecEntry eMult;
  NodeExecEntry eReduce, eReduceScatter, eBcast;
  vector<NodeCommEntry> vcSchur, vcLsolve;

 private:
  double tmOpStart;
};

class StochIterateResourcesMonitor
{
 public:
  StochIterateResourcesMonitor();
  virtual ~StochIterateResourcesMonitor();

  virtual void recIterateTm_start();
  virtual void recIterateTm_stop();
  static unsigned long getMaxMemUsage();

 public: 
  double tmIterate;
 protected:
  double tmIterateStart;

  vector<double> tmHistory;

};
#endif
