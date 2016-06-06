/* PIPS-NLP                                                             
 * Author: Nai-Yuan Chiang
 * (C) 2015 Argonne National Laboratory
 */

#ifndef STOCHNLPINFO
#define STOCHNLPINFO

#include "mpi.h"
#include "NlpInfo.h"

#include "OoqpVectorHandle.h"
#include "OoqpVector.h"
#include "DoubleMatrixHandle.h"

#include <vector>


class sData;
class sTree;
class SparseSymMatrix;
class SparseGenMatrix;
class NlpGenVars;
class stochasticInput;
class multiStageInputTree;


// for Jeq
//  x0        x1      x2
//  Bmat
//  Amat  Bmat
//  Amat   0  	   Bmat
//  E0     E1      E2

// for Jineq
//  x0        x1      x2
//  Dmat
//  Cmat  Dmat
//  Cmat   0  	   Dmat
// F0      F1      F2

class sInfo : public NlpInfo
{
public:

  SparseSymMatrix *Qdiag;
  SparseGenMatrix *Qborder;
  SparseGenMatrix *Amat;
  SparseGenMatrix *Bmat;
  SparseGenMatrix *Cmat;
  SparseGenMatrix *Dmat;
  SparseGenMatrix *Emat;
  SparseGenMatrix *Fmat;

  int locNx,locMy,locMz;

  std::vector<sInfo*> children;
  
  sTree* stochNode;

  sInfo* parent;

  /* MPI communicator */
  MPI_Comm mpiComm;

  sInfo();
  sInfo(sData *data_in);    
  sInfo(sData *data_in, stochasticInput &in){assert(0);}

  sInfo(int nx_in,int my_in,int mz_in,int nzH_in,int nzA_in,int nzC_in);
  sInfo( int nx_in,int my_in,int mz_in,int nzH_in,int nzA_in,int nzC_in,
				  int nxL_in,int nxU_in,int nsL_in,int nsU_in);  

  virtual ~sInfo();



  virtual void Hessian_FromSon( NlpGenVars * vars, double *tempFromParH ) = 0;
  virtual void ObjGrad_FromSon( NlpGenVars * vars, OoqpVector *grad, double *tempFromParH ) = 0;
  virtual void writeSolution( NlpGenVars * vars_){};

	virtual double ObjValue( NlpGenVars * vars)  = 0;

	virtual void ConstraintBody( NlpGenVars * vars, OoqpVector *conEq, OoqpVector *conIneq) = 0;

	virtual int ObjGrad( NlpGenVars * vars, OoqpVector *grad ) = 0;



	virtual void Hessian( NlpGenVars * vars, SymMatrix *Hess ) = 0;

	virtual void JacFull( NlpGenVars * vars, GenMatrix* JacA, GenMatrix* JacC) = 0;


	virtual void get_InitX0(OoqpVector* vX) = 0;

	virtual void createChildren( sData *data_in,stochasticInput& in) = 0;
//	void createChildren( sData *data_in);

	virtual void AddChild(sInfo* child);

	virtual void Emult(double beta,  OoqpVector& y,
			    double alpha, OoqpVector& x );
	virtual void Fmult(double beta,  OoqpVector& y,
			     double alpha, OoqpVector& x );

protected:
  void destroyChildren();

};



/** DUMMY VERSION 
 *
 */
class sInfoDummy : public sInfo {
protected:

public:
  sInfoDummy()
    : sInfo() {};

  virtual void AddChild(sInfo* child){};
  virtual void createChildren(sData* data_in, stochasticInput& in) {};

  virtual double ObjValue( NlpGenVars * vars) { return 0;};
  
  virtual void ConstraintBody( NlpGenVars * vars, OoqpVector *conEq, OoqpVector *conIneq){};

  virtual int ObjGrad( NlpGenVars * vars, OoqpVector *grad ){return 0;};

  

  virtual void Hessian( NlpGenVars * vars, SymMatrix *Hess ){};

  virtual void JacFull( NlpGenVars * vars, GenMatrix* JacA, GenMatrix* JacC){};


  virtual void get_InitX0(OoqpVector* vX){};

  virtual void Hessian_FromSon( NlpGenVars * vars, double *tempFromParH ){};
  virtual void ObjGrad_FromSon( NlpGenVars * vars, OoqpVector *grad, double *tempFromParH ){};   
  virtual void writeSolution( NlpGenVars * vars_){};

  virtual void Emult(double beta,  OoqpVector& y,
		     double alpha, OoqpVector& x ){};
  virtual void Fmult(double beta,  OoqpVector& y,
		     double alpha, OoqpVector& x ){};
};

  

#endif

