/* PIPS-NLP                                                         	*
 * Authors: Nai-Yuan Chiang                      		*
 * (C) 2015 Argonne National Laboratory			*/

#include "NlpInfoCallBack.h"
#include "NlpGenVars.h"
#include "OoqpVector.h"
#include <cmath>
#include <climits>
#include <cassert>
#include <cstdlib>


#include "SimpleVector.h"
#include "DoubleMatrixHandle.h"
#include "DoubleMatrix.h"
#include "SparseGenMatrix.h"
#include "SparseSymMatrix.h"

//#include "pipsipmNlp_C_callbacks.h"


enum{
Fixed = -1,
NoBound = 0,
LowBound = 1, 
UppBound = 2,
BothBound = 3,
};


//extern "C" typedef int (*eval_f_cb)(double* vec_x, double *obj);
//extern "C" typedef int (*eval_g_cb)(double* vec_x, double* vec_g);
//extern "C" typedef int (*eval_grad_f_cb)(double* vec_x, double* vec_grad_f);
//extern "C" typedef int (*eval_jac_g_cb)(double* vec_x, double* vec_Jac, int* iRows, int *kCols);
//extern "C" typedef int (*eval_h_cb)(double* vec_x, double* vec_lambda, double* vec_Hes, int* iRows, int *kCols);


void NlpInfoCallBack::doubleLexSort_Init( int first[], int n, int second[], double data[], const int ifEqCon)
{
  int fi, se, j, k, kinc, inc;
  double dtemp;
  int dtemp_goff;
  const int incs[]  = {1, 5, 19, 41, 109, 209, 505,
		       929, 2161, 3905, 8929, 16001, INT_MAX};

  if(ifEqCon==1){
    _JacAGoffTransMap = (int*) malloc(n*sizeof(int));
	for(k=0;k<n;k++)
	  _JacAGoffTransMap[k] = _dA_goff[k];
  }
  else{
    _JacCGoffTransMap = (int*) malloc(n*sizeof(int));
  	for(k=0;k<n;k++)
	  _JacCGoffTransMap[k] = _dC_goff[k];
  }
	
  for ( k = 0; incs[k] <= n/2; k++ ) ;

  kinc = k - 1;
  // incs[kinc] is the greatest value in the sequence that is also less
  // than or equal to n/2. 
  // If n == {0,1}, kinc == -1 and so no sort will take place.

  for( ; kinc >= 0; kinc-- ) {
	// Loop over all increments
	inc = incs[kinc];

	for ( k = inc; k < n; k++ ) {
	  dtemp = data[k];	  	
	  if(ifEqCon==1){
		dtemp_goff = _JacAGoffTransMap[k];
	  }
	  else{
		dtemp_goff = _JacCGoffTransMap[k];
	  }

	  fi = first[ k ];
	  se = second[ k];
	  for( j = k; j >= inc; j -= inc ) {
	  	if ( fi < first[j - inc] || 
			 ( fi == first[j - inc] &&
		       se < second[ j - inc ]) ) {
		  data[j]      	= data[j - inc];
		  first[j]   	= first[j - inc];
		  second[j]  	= second[j - inc];
		  if(ifEqCon==1){
			_JacAGoffTransMap[j] = _JacAGoffTransMap[j - inc];
		  }
		  else{
		    _JacCGoffTransMap[j] = _JacCGoffTransMap[j - inc];
		  }
		} else {
		  break;
		}
	  } 
	  data[j]    	= dtemp;
	  first[j]   	= fi;
	  second[j]  	= se;
	  if(ifEqCon==1){
		_JacAGoffTransMap[j] = dtemp_goff;
	  }
      else{
		_JacCGoffTransMap[j] = dtemp_goff;
	  }	  
	}
  } // End loop over all increments
 
}




void NlpInfoCallBack::doubleLexSort_ValOnly(int nzA, double *dataA, int nzC, double *dataC, double *dataJacFull)
{
  int k;
  
  for(k=0;k<nzA;k++){
  	dataA[k] = dataJacFull[_JacAGoffTransMap[k]];	
  }
  for(k=0;k<nzC;k++){
  	dataC[k] = dataJacFull[_JacCGoffTransMap[k]];
  }
}



NlpInfoCallBack::NlpInfoCallBack(eval_f_cb eval_f_in, eval_g_cb eval_g_in, eval_grad_f_cb eval_grad_f_in, 
  							  eval_jac_g_cb eval_jac_g_in, eval_h_cb eval_h_in, UserDataPtr usr_data_in)
  							  : NlpInfo()
{
	eval_f = eval_f_in;
	eval_g = eval_g_in;
	eval_grad_f = eval_grad_f_in;
	eval_jac_g = eval_jac_g_in;
	eval_h = eval_h_in;

	usrData = usr_data_in;
	_RowMap=NULL;
	_xStatus=NULL;
	_consStatus=NULL;

}

NlpInfoCallBack::~NlpInfoCallBack()
{
}


////////////////////////////////////////////////////////////////
void NlpInfoCallBack::_FindRowMap_AddSlack_NY(  int n_var, double *Lx, double *Ux,
														int m_con, double *Lg, double *Ug,
														int  nnzJac,
			int & nx, int & nnzQ,
			int & my, int & nnzA, 
			int & mz, int & nnzC, 
			int &nnzCL, int & nnzCU,
			int &nxL, int &nxU, int &nsL, int &nsU)
{
  int i, j, k;

    printf( "CP: m_con=%d\n",  m_con);
printf( "CP: n_var=%d\n",  n_var);

  _RowMap = new int[m_con];
  _xStatus = new int[n_var];
  _consStatus = new int[m_con];

  int findEq=0, findCLow=0, findXLow=0, findCUp=0, findXUp=0;

  nx = n_var;
  my =0; mz = 0;

  printf( "CP: %d\n",  0);

  for( i = 0; i < nx; i++ ) {
  	_xStatus[i] = NoBound;

	if (Lx[i] == Ux[i]){
	  _xStatus[i] = Fixed;
	}else{
	  if( Lx[i] > -1e20 ) {
		nxL++;
		_xStatus[i] += LowBound;
	  }
	  if( Ux[i] < 1e20 ) {
		nxU++;
		_xStatus[i] += UppBound;
	  }
	}
  }
  printf( "CP: %d\n",  1);

  printf( "CP: m_con=%d\n", m_con);

  // Count the number of rows in A and C and create a map between the 
  // rows of input version of the Jacobian, and the rows of A and C.

  // Count the number of nonzeros in A, the Jacobian of the equality
  // constraints, and the number of nonzeros in C, the Jacobian of the
  // inequality constraints.
  nnzA = 0; nnzC = 0; nnzCL=0; nnzCU=0;
  for( i = 0; i < m_con; i++ ) {
    _consStatus[i] = NoBound;	
    if ( Lg[i] == Ug[i] ) {
      _RowMap[i] = - (my + 1); // Negative values indicate an equality
      _consStatus[i] = Fixed;
      my++;
    } else{
      _RowMap[i] = mz;
      mz++;
	  if( Lg[i] > -1e20 ) {
		nsL++;
		_consStatus[i] += LowBound;
	  }
	  if( Ug[i] < 1e20 ) {
		nsU++;
		_consStatus[i] += UppBound;
	  }	  
    } 
  }

  printf( "CP: %d\n",  2);



  _invpMap = new int[my+nsL+nsU+nxL+nxU];
  _RowMap_CLow = new int[m_con];
  _RowMap_CUp  = new int[m_con];
  _RowMap_XLow = new int[n_var];
  _RowMap_XUp  = new int[n_var];

  printf( "CP: %d\n",  3);
  
  for( i = 0; i < m_con; i++ ) {
    switch(_consStatus[i]){
	  case Fixed: 
	  	 _invpMap[findEq++] = i;
		 break;
	  case LowBound: 
	  	 _invpMap[my + findCLow] = i;
		 _RowMap_CLow[i] = my + findCLow++;		 
		 break;
	  case UppBound: 
	  	 _invpMap[my + nsL + findCUp] = i;
		 _RowMap_CUp[i] = my + nsL + findCUp++;			 
	  	 break;
	  case BothBound: 
	  	 _invpMap[my + findCLow] = i;
		 _RowMap_CLow[i] = my + findCLow++;
		 
		 _invpMap[my + nsL + findCUp] = i;
		 _RowMap_CUp[i] = my + nsL + findCUp++;
	  	 break;		 
	}
  }
  printf( "CP: %d\n",  4);
  

  for( i = 0; i < n_var; i++ ) {
    switch(_xStatus[i]){
	  case LowBound: 
	  	 _invpMap[my + nsL + nsU + findXLow] = i;
		 _RowMap_XLow[i] = my + nsL + nsU + findXLow++;			 
		 break;
	  case UppBound: 
	  	 _invpMap[my + nsL + nsU + nxL + findXUp] = i;
		 _RowMap_XUp[i] = my + nsL + nsU + nxL + findXUp++;			 
	  	 break;
	  case BothBound: 
		 _invpMap[my + nsL + nsU + findXLow] = i;
		 _RowMap_XLow[i] = my + nsL + nsU + findXLow++;	 
		 _invpMap[my + nsL + nsU + nxL + findXUp] = i;
		 _RowMap_XUp[i] = my + nsL + nsU + nxL + findXUp++; 
	  	 break;			 
	}	
  }
  printf( "CP: %d\n",  5);

  assert(findEq == my);
  assert(findXLow == nxL);
  assert(findCLow == nsL);
  assert(findXUp == nxU); 
  assert(findCUp == nsU);  


  double *tempX = (double*) malloc((n_var)*sizeof(double)); 
  for(i = 0; i < n_var; i++ ) tempX[i]=0.0;
  double *tempJac = (double*) malloc((nnzJac)*sizeof(double)); 
  int *tempirow = (int*) malloc((nnzJac)*sizeof(int)); 
  int *tempkCol = (int*) malloc((n_var+1)*sizeof(int)); 

  printf( "CP: %d\n",  6);

  eval_jac_g(tempX, NULL,tempirow,tempkCol,usrData);

  printf( "CP: %d\n",  7);

  printf( "CP: nnzJac = %d \n",  nnzJac);

  for(k=0; k<n_var+1;k++){
	  printf( "CP: col=%d, colstart=%d\n",  k, tempkCol[k]); 
  }

  for(k=0; k<nnzJac;k++){
	  printf( "CP: k=%d, row=%d\n",  k, tempirow[k]); 
  }


  for(k=0; k<nnzJac;k++){
  	j=tempirow[k];

	if( _RowMap[j] < 0){
	  nnzA++;
//	  assert(_consStatus[j] == Fixed);
	}else if(_consStatus[j] == LowBound){
		nnzCL++;
		nnzC++;
	}else if(_consStatus[j] == UppBound){
		nnzCU++;
		nnzC++;		
	}else if(_consStatus[j] == BothBound){
		nnzCL++;
		nnzCU++;
		nnzC++;		
	}
  }

  printf( "CP: %d\n",  8);  
  free(tempX);
  printf( "CP: %d\n",  9);  
  free(tempJac);
  printf( "CP: %d\n",  10);  
  free(tempirow);
  printf( "CP: %d\n", 11);  
  free(tempkCol);
}



void NlpInfoCallBack::setBaseInfo(int nx_in,int my_in,int mz_in,int nzH_in,int nzA_in,int nzC_in,
					int nxL_in,int nxU_in,int nsL_in,int nsU_in)
{
   	   nx=(nx_in),
	   my=(my_in),
	   mz=(mz_in),
	   nzA=(nzA_in),
	   nzC=(nzC_in),
	   nzH=(nzH_in),
	   nsL=(nsL_in),
	   nsU=(nsU_in),
	   nxL=(nxL_in),
	   nxU=(nxU_in);
}




////////////////////////////////////////////////////////////////
void NlpInfoCallBack::_get_bounds( double XL[], double XU[], double GL[], double GU[],
		      double xlow[], int nx, char ixlow[],
		      double xupp[], char ixupp[],
		      double b[], int  my ,
		      double clow[], int  mz , char iclow[],
		      double cupp[], char icupp[] )
{
  int i;
  //double (*xbnd)[2] = (double (*)[2]) LUv;
  
  // Get xlow, xupp
  for( i = 0; i < nx; i++ ) {
//  	if (XL[i] == XU[i]){
//	  printf("have fixed var, relax bound! \n");
//      xlow[i]  = XL[i];
//      ixlow[i] = 1;	  
//	  xupp[i]  = XU[i];
//      ixupp[i] = 1;
//	}else 
	if( XL[i] > -1e20 ) {
      xlow[i]  = XL[i];
      ixlow[i] = 1;
    } else {
      xlow[i] = 0; ixlow[i] = 0;
    }
    if( XU[i] < 1e20 ) {
      xupp[i]  = XU[i];
      ixupp[i] = 1;
    } else {
      xupp[i] = 0; ixupp[i] = 0;
    }
  }

  // Get b, clow, cupp
  for( int iampl = 0; iampl < my+mz; iampl++ ) {
    
    if( _RowMap[iampl]<0 ) {
      i = - ( _RowMap[iampl] + 1 ); // Recover i from the negative value
      b[i] = GL[iampl];
    } else {
      i = _RowMap[iampl];

      if( GL[iampl] > -1e20 ) {
		clow[i]  = GL[iampl];
		iclow[i] = 1;
      } else {
		clow[i] = 0; iclow[i] = 0;
      }
      if( GU[iampl] < 1e20 ) {
		cupp[i]  = GU[iampl];
		icupp[i] = 1;
      } else {
		cupp[i] = 0; icupp[i] = 0;
      }
    }
  }
}



////////////////////////////////////////////////////////////////
void NlpInfoCallBack::_get_matrices_map(
			int  nx , int nnzQ,
			int  my , int nnzA,
			int  mz , int nnzC,
			int irowQ[], int jcolQ[], double dQ[],
			int irowA[], int jcolA[], double dA[],
			int irowC[], int jcolC[], double dC[])
{			
  int i, j, k;

  // Q in upper triangular form!

  double *xeltsTemp = (double*)malloc((nx)*sizeof(double));
  double *yeltsTemp = (double*)malloc((my+mz)*sizeof(double));
  double *eltsH = (double*)malloc((nnzQ)*sizeof(double));
  int *iRowsTemp = (int*)malloc((nnzQ)*sizeof(int));  
  int *kColsTemp = (int*)malloc((nx+1)*sizeof(int)); 

  for(i = 0; i < nx; i++ ) xeltsTemp[i]=0.0;
  for(i = 0; i < my+mz; i++ ) yeltsTemp[i]=0.0;

  eval_h(xeltsTemp, yeltsTemp, NULL, iRowsTemp, kColsTemp,usrData);

 {
	int kQ = 0;
	for( j = 0; j < nx; j++ ) {
	  for( k = kColsTemp[j]; k < kColsTemp[j+1]; k++ ) {
		i = iRowsTemp[k];
		// use only the upper triangle, but transpose to lower triangle for OOQP. Note that from columwise to rowwise
		if(i>j){
		  printf("i=%d, j=%d\n", i,j);
		}
		assert(i<=j);
		irowQ[kQ] = j; jcolQ[kQ] = i; dQ[kQ] = eltsH[k];
		kQ++;
	  }
	}
	assert(kQ==nnzQ);	
  }

  free(xeltsTemp);
  free(yeltsTemp);
  free(eltsH);
  free(iRowsTemp);
  free(kColsTemp);


  int nnzJac = nnzA+nnzC;
  
  double *tempX = (double*) malloc((nx)*sizeof(double)); 
  for(i = 0; i < nx; i++ ) tempX[i]=0.0;
  double *tempJac = (double*) malloc((nnzJac)*sizeof(double)); 
  int *tempirow = (int*) malloc((nnzJac)*sizeof(int)); 
  int *tempkCol = (int*) malloc((nx+1)*sizeof(int)); 

  eval_jac_g(tempX, NULL,tempirow,tempkCol,usrData);

  _dA_goff = (int*)malloc((nnzA)*sizeof(int));
  _dC_goff = (int*)malloc((nnzC)*sizeof(int));

  
  /* find row_number of Jac*/	

  int kA = 0, kC = 0;
  for( j = 0; j < nx; j++ ) {
    for( k = tempkCol[j]; k < tempkCol[j+1]; k++ ) {
      int iampl = tempirow[k];
      if ( _RowMap[iampl] < 0 ) { // A negative value in rowMap indicates
		// and equality constraint.
		i = - ( _RowMap[iampl] + 1 );

		assert( kA < nnzA );
		irowA[kA] = i;  jcolA[kA] = j;  dA[kA] = tempJac[k];
		_dA_goff[kA]=k;

		kA++;
      } else {
     	i = _RowMap[iampl];

	  	assert( kC < nnzC );
	  	irowC[kC] = i;  jcolC[kC] = j;  dC[kC] = tempJac[k];
		_dC_goff[kC]=k;

	  	kC++;
      }
    }
  }

  // Transpose everything --- OOQP uses row-wise index, but OOPS and AMPL use col-wise index
  // NY version also save the Goff index which can be reused later by calling doubleLexSort_ValOnly
  doubleLexSort_Init( irowA, nnzA, jcolA, dA, 1);
  doubleLexSort_Init( irowC, nnzC, jcolC, dC, 0);

  free(tempX);
  free(tempJac);
  free(tempirow);
  free(tempkCol);

}





double 
NlpInfoCallBack::ObjValue( NlpGenVars * vars)
{
  SimpleVector & vec_x = dynamic_cast<SimpleVector &>(*vars->x);

  double objWrk;

  eval_f( vec_x.elements() , &objWrk,usrData);

  return objWrk;
}



int 
NlpInfoCallBack::ObjGrad( NlpGenVars * vars, OoqpVector *grad)
{
  SimpleVector & vec_x = dynamic_cast<SimpleVector &>(*vars->x);
  SimpleVector & vec_grad = dynamic_cast<SimpleVector &>(*grad);

  eval_grad_f(vec_x.elements(), vec_grad.elements(),usrData);  

  return 1;
}


void 
NlpInfoCallBack::ConstraintBody( NlpGenVars * vars, OoqpVector *conEq, OoqpVector *conIneq)
{
  int i;
  SimpleVector & vec_x = dynamic_cast<SimpleVector &>(*vars->x);
  SimpleVector & con_eq = dynamic_cast<SimpleVector &>(*conEq);
  SimpleVector & con_ineq = dynamic_cast<SimpleVector &>(*conIneq);

  double *tempCon = (double*) malloc((my+mz)*sizeof(double));  

  eval_g(vec_x.elements(), tempCon,usrData);  

  for( int iampl = 0; iampl < (my+mz); iampl++ ) {
    if( _RowMap[iampl]<0 ) {
      i = - ( _RowMap[iampl] + 1 ); 	// Recover i from the negative value: equality constraint

      con_eq[i] = tempCon[iampl];
    } else {
      i = _RowMap[iampl]; 			// Inequality constraint
	  con_ineq[i] = tempCon[iampl];
    }
  }

  free(tempCon);

}


void NlpInfoCallBack::JacFull( NlpGenVars * vars, GenMatrix* JacA, GenMatrix* JacC)
{  
  SimpleVector & vec_x = dynamic_cast<SimpleVector &>(*vars->x);
  double *tempJac = (double*) malloc((nzA+nzC)*sizeof(double)); 



  int *iRows = (int*)malloc((nzA+nzC)*sizeof(int)); 
  int *kCols = (int*)malloc((nx+1)*sizeof(int)); 

  eval_jac_g(vec_x.elements(), tempJac,iRows,kCols,usrData);
  

  SparseGenMatrix * JacASp = dynamic_cast<SparseGenMatrix *>(JacA);
  SparseGenMatrix * JacCSp = dynamic_cast<SparseGenMatrix *>(JacC);
//  JacASp->printMatrixInMatlab("matA");
//  JacCSp->printMatrixInMatlab("matC");

  doubleLexSort_ValOnly(nzA,JacASp->M(),nzC,JacCSp->M(),tempJac);

  free(tempJac);
  free(kCols);
  free(iRows);  
}


//Make sure we are using Hessian in triangular form!
void NlpInfoCallBack::Hessian( NlpGenVars * vars, SymMatrix* Hess )
{

   SimpleVector & vec_x = dynamic_cast<SimpleVector &>(*vars->x);
   SimpleVector & vec_y = dynamic_cast<SimpleVector &>(*vars->y);
   SimpleVector & vec_z = dynamic_cast<SimpleVector &>(*vars->z);

   double *dualWrk = (double*)malloc((mz+my)*sizeof(double)); 
   int *iRows = (int*)malloc((nzH)*sizeof(int)); 
   int *kCols = (int*)malloc((nx+1)*sizeof(int)); 

   int i;
  
   // assemble dual variable
   for( int iampl = 0; iampl < (mz+my); iampl++ ) {
     if( _RowMap[iampl]<0 ) {
       i = - ( _RowMap[iampl] + 1 ); // Recover i from the negative value
       dualWrk[iampl] = -vec_y[i];
     } else {
       i = _RowMap[iampl];
	   dualWrk[iampl] = -vec_z[i]; 
     }
   }

  SparseSymMatrix * HesSp = dynamic_cast<SparseSymMatrix *>(Hess);
  eval_h(vec_x.elements(), dualWrk, HesSp->M(),iRows, kCols,usrData);
  
  free(dualWrk);
  free(kCols);
  free(iRows); 
}





void
NlpInfoCallBack::get_InitX0(OoqpVector* vX)
{
  printf("Get init point \n\n\n");
  vX->print();
//  double *tempX = (double*) malloc(nx*sizeof(double));
//  ampl_get_InitX0(tempX);
//  vX->copyFromArray(tempX);
 

}


