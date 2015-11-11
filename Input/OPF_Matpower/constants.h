/* PIPS-NLP                                                         	*
 * Authors: Nai-Yuan Chiang                      		*
 * (C) 2015 Argonne National Laboratory			*/

#ifndef DCOPFCONSTANTS_H
#define DCOPFCONSTANTS_H

#include <math.h>
#include <string.h>

#define w_s (2*M_PI*freq)
//#define epsilon 1E-8
#define PS_MAXLINE 1000
#define PV_BUS 2
#define PQ_BUS 1
#define NGEN_AT_BUS_MAX 15
#define NLOAD_AT_BUS_MAX 10
#define NPHASE 1

/* Type of variables */
#define DIFF_EQ 1 /* Differential equation */
#define ALG_EQ  0 /* Algebraic equation */


#define MAXCONNLINES 20

class OoqpVector;
class Variables;

/*
  Preconding information - struct that hold values for the matrix
*/
class PreCondInfo{
public:
  int n_row;
  int n_col;
  int jac_nnz;
  int hes_nnz;
  
  int Aug_nnz;

  int n_scen;

  int* jac_irow;
  int* jac_jcol;  
  int* jac_kcolbeg;
  double* jac_value;

  int* hes_irow;
  int* hes_jcol;  
  int* hes_kcolbeg;
  double* hes_value;

  int *firstVarMap;
  int *firstConMap;
  
  int **varMap;
  int **conMap;
  

  int nx_1st;
  int my_1st;
  int ms_1st;
  int n1stCon;

  OoqpVector *full_X;
  OoqpVector *full_Y;
  OoqpVector *full_step_x;
  OoqpVector *full_step_y;

  OoqpVector *grad_x;
  OoqpVector *cons_b;

  double *full_X_val;
  double *full_Y_val;
  double *full_step_x_val;
  double *full_step_y_val;

  double *grad_x_val;
  double *cons_b_val;


  Variables *curr_Iter;
  Variables *curr_Step;


  bool setGrad;

  PreCondInfo(){};
  PreCondInfo(int nRow, int nCol, int jacnnz_, int hesnnz_, int n1stVar_, int n1stEqCon_, int n1stInCon_,int n_scen_);
  ~PreCondInfo();
  
  void setAggregationVarMap( int scen,int *varMap_);
  void setAggregationConMap( int scen,int *conMap_);
  
  void setAggregationJacMatrix( int *irow_, int *jcol_, int *kcolstrart_, double *ele_);
  void setAggregationHesMatrix( int *irow_, int *jcol_, int *kcolstrart_, double *ele_);
  void setAggregation1stVarMap(int *varMap_);
  void setAggregation1stConMap(int *conMap_);


  
  void setPrecondSysRHS(double *sysRhs_){};
};

inline
PreCondInfo::PreCondInfo(int nRow, int nCol, int jacnnz_, int hesnnz_, int n1stVar_,
							int n1stEqCon_, int n1stInCon_, int n_scen_)
 : 	n_row(nRow),n_col(nCol),jac_nnz(jacnnz_), hes_nnz(hesnnz_),n_scen(n_scen_),setGrad(false),
 	full_X(NULL),full_Y(NULL),full_step_x(NULL),full_step_y(NULL),
 	grad_x(NULL),cons_b(NULL)
{ 
  jac_irow 		= new int[jac_nnz];
  jac_jcol 		= new int[jac_nnz];
  jac_kcolbeg 	= new int[nCol+1];
  jac_value 	= new double[jac_nnz];

  hes_irow 		= new int[hes_nnz];
  hes_jcol 		= new int[hes_nnz];
  hes_kcolbeg 	= new int[nCol+1];
  hes_value 	= new double[hes_nnz];

  nx_1st	= n1stVar_; 
  my_1st	= n1stEqCon_;
  ms_1st	= n1stInCon_;
  
  n1stCon	= my_1st + ms_1st;

  varMap	= new int*[n_scen];
  conMap	= new int*[n_scen];

  Aug_nnz  = jac_nnz + hes_nnz + n_col+n_row;

  grad_x_val 	= new double[n_col];
  cons_b_val 	= new double[n_row];
 
}

inline
PreCondInfo::~PreCondInfo()
{ 
  delete [] jac_irow;
  delete [] jac_jcol;
  delete [] jac_kcolbeg;
  delete [] jac_value;

  delete [] hes_irow;
  delete [] hes_jcol;
  delete [] hes_kcolbeg;
  delete [] hes_value;
  
  delete [] varMap;
  delete [] conMap;

  delete [] grad_x_val;
  delete [] cons_b_val;

}

inline void
PreCondInfo::setAggregationVarMap( int scen,int *varMap_)
{ 
  varMap[scen] = varMap_;
}

inline void
PreCondInfo::setAggregationConMap( int scen,int *conMap_)
{ 
  conMap[scen] = conMap_;
}


inline void
PreCondInfo::setAggregation1stVarMap(int *varMap_)
{ 
  firstVarMap = varMap_;
}

inline void
PreCondInfo::setAggregation1stConMap(int *conMap_)
{ 
  firstConMap = conMap_;
}

inline void
PreCondInfo::setAggregationJacMatrix( int *irow_, int *jcol_, int *kcolstrart_, double *ele_)
{ 
  memcpy( jac_irow, irow_, jac_nnz * sizeof( int ) );
  memcpy( jac_jcol, jcol_, jac_nnz * sizeof( int ) );
  memcpy( jac_kcolbeg, kcolstrart_, (n_col+1) * sizeof( int ) );
  memcpy( jac_value, ele_, jac_nnz * sizeof( double ) );
}

inline void
PreCondInfo::setAggregationHesMatrix( int *irow_, int *jcol_, int *kcolstrart_, double *ele_)
{ 
  memcpy( hes_irow, irow_, hes_nnz * sizeof( int ) );
  memcpy( hes_jcol, jcol_, hes_nnz * sizeof( int ) );
  memcpy( hes_kcolbeg, kcolstrart_, (n_col+1) * sizeof( int ) );
  memcpy( hes_value, ele_, hes_nnz * sizeof( double ) );
}

#if 0
inline void
PreCondInfo::copyVal( int *irow_, int *jcol_, int *kcolstrart_, double *ele_)
{ 
  memcpy( hes_irow, irow_, hes_nnz * sizeof( int ) );
  memcpy( hes_jcol, jcol_, hes_nnz * sizeof( int ) );
  memcpy( hes_kcolbeg, kcolstrart_, (n_col+1) * sizeof( int ) );
  memcpy( hes_value, ele_, hes_nnz * sizeof( double ) );
}
#endif
#endif
