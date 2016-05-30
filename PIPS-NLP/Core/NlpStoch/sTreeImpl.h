/* PIPS-IPM                                                           *
 * Author:  Cosmin G. Petra                                           *
 * (C) 2012 Argonne National Laboratory. See Copyright Notification.  */

#ifndef STOCH_TREE_IMPL
#define STOCH_TREE_IMPL

#include "sTree.h"
#include "stochasticInput.hpp"

/** A full implementation of sTree that is currently used, the other implementation
 *  being sTreeCallbacks (obsolete)
 * 
 */

class sTreeImpl : public sTree
{
 public:
//  sTreeImpl(){}
  sTreeImpl(stochasticInput &in, MPI_Comm comm=MPI_COMM_WORLD);

 private: 
  sTreeImpl(int idx, stochasticInput &in);
 public:
  virtual ~sTreeImpl();

  virtual StochSymMatrix*   createQ() const;
  virtual StochVector*      createc() const;

  virtual StochVector*      createxlow()  const;
  virtual StochVector*      createixlow() const;
  virtual StochVector*      createxupp()  const;
  virtual StochVector*      createixupp() const;


  virtual StochGenMatrix*   createA() ;
  virtual StochVector*      createb() const;


  virtual StochGenMatrix*   createC() ;
  virtual StochVector*      createclow()  const;
  virtual StochVector*      createiclow() const;
  virtual StochVector*      createcupp()  const;
  virtual StochVector*      createicupp() const;

  virtual StochVector*      createCeqBody() const ;
  virtual StochVector*      createCineqBody()  const ;
  virtual StochVector*      createBarrGrad() const ;  

  virtual int nx() const;
  virtual int my() const; 
  virtual int mz() const; 
  virtual int id() const; 
  virtual int mle() const{return m_mle;};
  virtual int mli() const{return m_mli;};
  virtual void computeGlobalSizes();
  virtual void loadLocalSizes();
  virtual void  get_FistStageSize(int& nx, int& my, int& mz);

 protected:
  int m_id;
  stochasticInput &in;

  sTreeImpl* parent;

  size_t m_nx, m_my, m_mz, m_mle, m_mli;

  int compute_nFirstStageEq();
  int compute_nSecondStageEq(int scen);
};



#endif
