/* PIPS
   Authors: Cosmin Petra
   See license and copyright information in the documentation */

#ifndef STOCH_TREE_PRECOND
#define STOCH_TREE_PRECOND

#include "StochTree.h"
class StochTreePrecond : public StochTree {
 public:
  StochTreePrecond(StochInputTree* root);
  virtual ~StochTreePrecond();
 private:
   StochTreePrecond();

 public:
   void assignProcesses  ( );
   void assignProcesses  (MPI_Comm, vector<int>&);
 private:
   virtual void assignToPreconditioner(vector<double>& vChildsLoad, 
				       vector<vector<int> >& mChilds2Procs);
 public:
   MPI_Comm commWorkers;
};

#endif
