/* PIPS-NLP                                                         	*
 * Authors: Nai-Yuan Chiang                      		*
 * (C) 2015 Argonne National Laboratory			*/

#ifndef NLPINFO_CALLBACK
#define NLPINFO_CALLBACK


#include "OoqpVectorHandle.h"
#include "OoqpVector.h"
#include "DoubleMatrixHandle.h"
#include "NlpInfo.h"


typedef void * UserDataPtr;

extern "C" typedef int (*eval_f_cb)(double* vec_x, double *obj, UserDataPtr user_data);
extern "C" typedef int (*eval_g_cb)(double* vec_x, double* vec_g, UserDataPtr user_data);
extern "C" typedef int (*eval_grad_f_cb)(double* vec_x, double* vec_grad_f, UserDataPtr user_data);
extern "C" typedef int (*eval_jac_g_cb)(double* vec_x, double* vec_Jac, int* iRows, int *kCols, UserDataPtr user_data);
extern "C" typedef int (*eval_h_cb)(double* vec_x, double* vec_lambda, double* vec_Hes, int* iRows, int *kCols, UserDataPtr user_data);


class NlpGenVars;
	
class NlpInfoCallBack : public NlpInfo
{
private:
  eval_f_cb eval_f;
  eval_g_cb eval_g;
  eval_grad_f_cb eval_grad_f;
  eval_jac_g_cb eval_jac_g;
  eval_h_cb eval_h;	

  UserDataPtr usrData;

public:
  
  NlpInfoCallBack(eval_f_cb eval_f_in, eval_g_cb eval_g_in, eval_grad_f_cb eval_grad_f_in, 
  							  eval_jac_g_cb eval_jac_g_in, eval_h_cb eval_h_in, UserDataPtr user_data);    
  virtual ~NlpInfoCallBack();

  long long  nx,my,mz;
  long long nzH,nzA,nzC;
  long long nsL, nsU, nxL, nxU;
//  int *rowMap;
  
  virtual double ObjValue( NlpGenVars * vars) ;
  
  virtual void ConstraintBody( NlpGenVars * vars, OoqpVector *conEq, OoqpVector *conIneq);

  virtual int ObjGrad( NlpGenVars * vars, OoqpVector *grad );

  

  virtual void Hessian( NlpGenVars * vars, SymMatrix *Hess );

  virtual void JacFull( NlpGenVars * vars, GenMatrix* JacA, GenMatrix* JacC);

  virtual void get_InitX0(OoqpVector* vX);


  SymMatrix *Q;
  GenMatrix *A;
  GenMatrix *C;
  OoqpVector *g;
  OoqpVector *bA;

  OoqpVector *my_temp;


  int *_RowMap;
  int *_invpMap, *_xStatus, *_consStatus;
  
  int *_dA_goff, *_dC_goff;
  int *_RowMap_CLow, *_RowMap_CUp, *_RowMap_XLow, *_RowMap_XUp;

  int nzJac, nzHes;

  void setBaseInfo(int nx_in,int my_in,int mz_in,int nzH_in,int nzA_in,int nzC_in,
				  int nxL_in,int nxU_in,int nsL_in,int nsU_in);




  void _FindRowMap_AddSlack_NY(  int n_var, double *Lx, double *Ux,
														int m_con, double *Lg, double *Ug,
														int nnzJac,
			int & nx, int & nnzQ,
			int & my, int & nnzA, 
			int & mz, int & nnzC, 
			int &nnzCL, int & nnzCU,
			int &nxL, int &nxU, int &nsL, int &nsU);

  void _get_bounds( double XL[], double XU[], double GL[], double GU[],
		      double xlow[], int nx, char ixlow[],
		      double xupp[], char ixupp[],
		      double b[], int  my ,
		      double clow[], int  mz , char iclow[],
		      double cupp[], char icupp[] );

  void _get_matrices_map(
			int  nx , int nnzQ,
			int  my , int nnzA,
			int  mz , int nnzC,
			int irowQ[], int jcolQ[], double dQ[],
			int irowA[], int jcolA[], double dA[],
			int irowC[], int jcolC[], double dC[]);

private:
	int *_JacAGoffTransMap,*_JacCGoffTransMap;
	
  void doubleLexSort_ValOnly(int nzA, double *dataA, int nzC, double *dataC, double *dataJacFull);
  void doubleLexSort_Init( int first[], int n, int second[], double data[], const int ifEqCon);
  
};

#endif

