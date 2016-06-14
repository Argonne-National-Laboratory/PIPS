/* PIPS-NLP                                                         	*
 * Authors: Nai-Yuan Chiang                      		*
 * (C) 2015 Argonne National Laboratory			*/

#include <cstdio>
#include <cassert>
#include "getAmplFunction.h"

#include <climits>

#include "asl_pfgh.h"
#include "getstub.h"


#define asl cur_ASL

// ASL options
extern int nobj;
extern int hes_obj;
extern int hes_con;
extern int hes_tri;
extern double OW[1];

static int *amplRowMap;
static int *invpRowMap;
static int *XStatus;
static int *ConStatus;
static int *JacAGoffTransMap, *JacCGoffTransMap, *dA_goff, *dC_goff, *dC_L_goff, *dC_U_goff,*dC_UL_goff;

static int * amplRowMap_CLow;
static int * amplRowMap_CUp;
static int * amplRowMap_XLow;
static int * amplRowMap_XUp;

void doubleLexSort_Init( int first[], int n, int second[], double data[], const int ifEqCon)
{
  int fi, se, j, k, kinc, inc;
  double dtemp;
  int dtemp_goff;
  const int incs[]  = {1, 5, 19, 41, 109, 209, 505,
		       929, 2161, 3905, 8929, 16001, INT_MAX};

  if(ifEqCon==1){
    JacAGoffTransMap = (int*) malloc(n*sizeof(int));
	for(k=0;k<n;k++)
	  JacAGoffTransMap[k] = dA_goff[k];
  }
  else{
    JacCGoffTransMap = (int*) malloc(n*sizeof(int));
  	for(k=0;k<n;k++)
	  JacCGoffTransMap[k] = dC_goff[k];
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
		dtemp_goff = JacAGoffTransMap[k];
	  }
	  else{
		dtemp_goff = JacCGoffTransMap[k];
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
			JacAGoffTransMap[j] = JacAGoffTransMap[j - inc];
		  }
		  else{
		    JacCGoffTransMap[j] = JacCGoffTransMap[j - inc];
		  }
		} else {
		  break;
		}
	  } 
	  data[j]    	= dtemp;
	  first[j]   	= fi;
	  second[j]  	= se;
	  if(ifEqCon==1){
		JacAGoffTransMap[j] = dtemp_goff;
	  }
      else{
		JacCGoffTransMap[j] = dtemp_goff;
	  }	  
	}
  } // End loop over all increments
 
}


void doubleLexSort_ValOnly(int nzA, double *dataA, int nzC, double *dataC, double *dataJacFull)
{
  int k;
  
  for(k=0;k<nzA;k++){
  	dataA[k] = dataJacFull[JacAGoffTransMap[k]];	
  }
  for(k=0;k<nzC;k++){
  	dataC[k] = dataJacFull[JacCGoffTransMap[k]];
  }
}


void
ampl_count_sizes_SplitSlack(  fint irow[], fint kcol[],
			int & nx, int & nnzQ,
			int & my, int & nnzA, int & mz, int & nnzC, 
			const int full_size, 
			int &nnzCL, int & nnzCU,
			int &nxL, int &nxU, int &nsL, int &nsU)
{
  int i, j, k;

  int *rowMap = new int[n_con];
  int *xStatus = new int[n_var];
  int *consStatus = new int[n_con];
  int *invpMap;

  int findEq=0, findCLow=0, findXLow=0, findCUp=0, findXUp=0;
 
  nx = n_var;
  my =0; mz = 0;

  for( i = 0; i < nx; i++ ) {
  	xStatus[i] = NoBound;

	if (LUv[i] == Uvx[i]){
	  xStatus[i] = Fixed;
	}else{
	  if( LUv[i] > -1e20 ) {
		nxL++;
		xStatus[i] += LowBound;
	  }
	  if( Uvx[i] < 1e20 ) {
		nxU++;
		xStatus[i] += UppBound;
	  }
	}
  }

  // if we use full size Hessian from AMPL
  if(1==full_size){
    // Our internal data structures need only a triangle of Q.  Count
    // the number of non-zeros in this triangle.
    nnzQ = 0;
    for( j = 0; j < n_var; j++ ) {
      for( k = kcol[j]; k < kcol[j+1]; k++ ) {
        i = irow[k];
        // count only the upper triangle
        if( i > j ) break;
        nnzQ++;
      }
    }
  }



  // Count the number of rows in A and C and create a map between the 
  // rows of ampl's version of the Jacobian, and the rows of A and C.

  // Count the number of nonzeros in A, the Jacobian of the equality
  // constraints, and the number of nonzeros in C, the Jacobian of the
  // inequality constraints.
  nnzA = 0; nnzC = 0; nnzCL=0; nnzCU=0;
  for( i = 0; i < n_con; i++ ) {
    consStatus[i] = NoBound;	
    if ( LUrhs[i] == Urhsx[i] ) {
      rowMap[i] = - (my + 1); // Negative values indicate an equality
      consStatus[i] = Fixed;
      my++;
    } else{
      rowMap[i] = mz;
      mz++;
	  if( LUrhs[i] > -1e20 ) {
		nsL++;
		consStatus[i] += LowBound;
	  }
	  if( Urhsx[i] < 1e20 ) {
		nsU++;
		consStatus[i] += UppBound;
	  }	  
    } 
  }


  invpMap = new int[my+nsL+nsU+nxL+nxU];
  int * RowMap_CLow = new int[n_con];
  int * RowMap_CUp  = new int[n_con];
  int * RowMap_XLow = new int[n_var];
  int * RowMap_XUp  = new int[n_var];

  
  for( i = 0; i < n_con; i++ ) {

    switch(consStatus[i]){
	  case Fixed: 
	  	 invpMap[findEq++] = i;
		 break;
	  case LowBound: 
	  	 invpMap[my + findCLow] = i;
		 RowMap_CLow[i] = my + findCLow++;		 
		 break;
	  case UppBound: 
	  	 invpMap[my + nsL + findCUp] = i;
		 RowMap_CUp[i] = my + nsL + findCUp++;			 
	  	 break;
	  case BothBound: 
	  	 invpMap[my + findCLow] = i;
		 RowMap_CLow[i] = my + findCLow++;
		 
		 invpMap[my + nsL + findCUp] = i;
		 RowMap_CUp[i] = my + nsL + findCUp++;
	  	 break;		 
	}
	
  }

  for( i = 0; i < n_var; i++ ) {

    switch(xStatus[i]){
	  case LowBound: 
	  	 invpMap[my + nsL + nsU + findXLow] = i;
		 RowMap_XLow[i] = my + nsL + nsU + findXLow++;			 
		 break;
	  case UppBound: 
	  	 invpMap[my + nsL + nsU + nxL + findXUp] = i;
		 RowMap_XUp[i] = my + nsL + nsU + nxL + findXUp++;			 
	  	 break;
	  case BothBound: 
		 invpMap[my + nsL + nsU + findXLow] = i;
		 RowMap_XLow[i] = my + nsL + nsU + findXLow++;	 
		 invpMap[my + nsL + nsU + nxL + findXUp] = i;
		 RowMap_XUp[i] = my + nsL + nsU + nxL + findXUp++; 
	  	 break;			 
	}
	
  }

  assert(findEq == my);
  assert(findXLow == nxL);
  assert(findCLow == nsL);
  assert(findXUp == nxU); 
  assert(findCUp == nsU);  

  cgrad *cg;

  for(j=0; j<n_con;j++){
	if( rowMap[j] < 0){
	  for(cg = Cgrad[j]; cg; cg = cg->next){
	  	nnzA++;
	  }
	  assert(consStatus[j] == Fixed);
	}else if(consStatus[j] == LowBound){
	  for(cg = Cgrad[j]; cg; cg = cg->next){
		nnzCL++;
		nnzC++;
	  }
	}else if(consStatus[j] == UppBound){
	  for(cg = Cgrad[j]; cg; cg = cg->next){
		nnzCU++;
		nnzC++;		
	  }
	}else if(consStatus[j] == BothBound){
	  for(cg = Cgrad[j]; cg; cg = cg->next){
		nnzCL++;
		nnzCU++;
		nnzC++;		
	  }
	}
  }

  amplRowMap = rowMap;
  invpRowMap = invpMap;
  XStatus = xStatus;
  ConStatus = consStatus;

  amplRowMap_CLow = RowMap_CLow;
  amplRowMap_CUp = RowMap_CUp;
  amplRowMap_XLow = RowMap_XLow;
  amplRowMap_XUp = RowMap_XUp;
  
}


void ampl_get_InitX0(double *varsX)
{
  for (int i=0; i<n_var; i++){
  	if(havex0[i])
	  varsX[i] = X0[i];
	else{
	  if (LUv[i] == Uvx[i]){
        varsX[i]  = LUv[i];
	  }else if (Uvx[i] > 1e20 && LUv[i] < -1e20){
        varsX[i]  = 0;
	  }else{	
	  	varsX[i]=0;
  	    if( LUv[i] > -1e20 ) 
          varsX[i]  += 0;//LUv[i]+1e-8;
	    if( Uvx[i] < 1e20 ) 
          varsX[i]  += 0;//Uvx[i]-1e-8;
	    if( Uvx[i] < 1e20 && LUv[i] > -1e20)
	      varsX[i] /=2;
	  }
	}
  }
}


double ampl_get_Obj(double *varsX)
{
  double objv=0;
  fint nerror = -1;

  if(varsX){
    objv = objval(nobj, varsX, &nerror);	
  }else{
    double *dxWrk = (double*) malloc (n_var*sizeof(double));
    objv = objval(nobj, dxWrk, &nerror);
    free(dxWrk);
  }
  return objv;
}



void ampl_get_Cons(double *varsX, double *consEqElt, double *consIneqElt)
{
  int i;
  fint nerror = -1;
  double *tempCon = (double*) malloc(n_con*sizeof(double));  


  if(varsX){
  	conval(varsX, tempCon, &nerror);	
  }else{
    double *dxWrk = (double*) malloc (n_var*sizeof(double));
  	conval(varsX, tempCon, &nerror);	
    free(dxWrk);
  }


  for( int iampl = 0; iampl < n_con; iampl++ ) {
    if( amplRowMap[iampl]<0 ) {
      i = - ( amplRowMap[iampl] + 1 ); 	// Recover i from the negative value: equality constraint

      consEqElt[i] = tempCon[iampl];
    } else {
      i = amplRowMap[iampl]; 			// Inequality constraint
	  consIneqElt[i] = tempCon[iampl];
    }
  }

  free(tempCon);
}



////////////////////////////////////////////////////////////////
void ampl_get_ObjGrad( double *varsX,  double Objgrad[])
{

  ograd * og;
  fint nerror = -1;

  if(varsX){
	objgrd(nobj, varsX, Objgrad, &nerror);  
  }else{
    double *dxWrk = (double*) malloc (n_var*sizeof(double));
	objgrd(nobj, dxWrk, Objgrad, &nerror);	
    free(dxWrk);
  }
 
}


////////////////////////////////////////////////////////////////
void ampl_get_bounds( 
		      double xlow[], int nx, char ixlow[],
		      double xupp[], char ixupp[],
		      double b[], int /* my */,
		      double clow[], int /* mz */, char iclow[],
		      double cupp[], char icupp[] )
{
  int i;
  
  // Get xlow, xupp
  for( i = 0; i < nx; i++ ) {
  	if (LUv[i] == Uvx[i]){
      xlow[i]  = LUv[i];
      ixlow[i] = 1;	  
	  xupp[i]  = Uvx[i];
      ixupp[i] = 1;
	}
    else if( LUv[i] > -1e20 ) {
      xlow[i]  = LUv[i];
      ixlow[i] = 1;
    } else {
      xlow[i] = 0; ixlow[i] = 0;
    }
    if( Uvx[i] < 1e20 ) {
      xupp[i]  = Uvx[i];
      ixupp[i] = 1;
    } else {
      xupp[i] = 0; ixupp[i] = 0;
    }
  }

  // Get b, clow, cupp
  for( int iampl = 0; iampl < n_con; iampl++ ) {
    
    if( amplRowMap[iampl]<0 ) {
      i = - ( amplRowMap[iampl] + 1 ); // Recover i from the negative value
      // in rowMap

      b[i] = LUrhs[iampl];
    } else {
      i = amplRowMap[iampl];

      if( LUrhs[iampl] > -1e20 ) {
		clow[i]  = LUrhs[iampl];
		iclow[i] = 1;
      } else {
		clow[i] = 0; iclow[i] = 0;
      }
      if( Urhsx[iampl] < 1e20 ) {
		cupp[i]  = Urhsx[iampl];
		icupp[i] = 1;
      } else {
		cupp[i] = 0; icupp[i] = 0;
      }
    }
  }
}

//the number of nonzeros in the sparse Hessian W of the Lagrangian (if uptri = 0) or its upper triangle (if uptri = 1)
int 
ampl_get_nnz_Hessian_Tri()
{
  int nzH = sphsetup(-1, hes_obj, hes_con, hes_tri); 
  return nzH;
}


void
ampl_get_Hessian_Tri(double *varsX, double *Helts, double *Yelts, double *Zelts)
{
  double *dualWrk = new double[n_con]; 
  int i;
  
  // Get b, clow, cupp
  for( int iampl = 0; iampl < n_con; iampl++ ) {
    
    if( amplRowMap[iampl]<0 ) {
      i = - ( amplRowMap[iampl] + 1 ); // Recover i from the negative value
      if(Yelts)
	  	dualWrk[iampl] = Yelts[i];
	  else
	  	dualWrk[iampl] = 0;
    } else {
      i = amplRowMap[iampl];
	  if(Zelts)
	  	dualWrk[iampl] = Zelts[i];
	  else
	  	dualWrk[iampl] = 0;	  
    }
  }
  
  sphes(Helts, -1, OW, dualWrk);

  delete [] dualWrk;
}


void
ampl_get_Jac(double *varsX, const int nzA, double *JacAelts, const int nzC, double *JacCelts)
{
  int i;
  fint nerror = -1;
  double *tempJac = (double*) malloc(nzc*sizeof(double)); 

  if(varsX){
	jacval(varsX, tempJac, &nerror);
  }else{
    double *dxWrk = (double*) malloc (n_var*sizeof(double));
  	jacval(varsX, tempJac, &nerror);	
    free(dxWrk);
  }

 doubleLexSort_ValOnly(nzA,JacAelts,nzC,JacCelts,tempJac);

 free(tempJac);

}


//can be called if hessian is returned in triangular form, otherwise more works need to be done as the get_matrices routine
void
ampl_GetHessianInit(double *varsX, double *Helts,double *Yelts)
{
  OW[0] = 1;

  if(varsX==NULL)
  	xunknown();
  if(Yelts)
    sphes(Helts, -1, OW, Yelts);
  else{
  	double *tempY = (double*)malloc(n_con*sizeof(double));
    sphes(Helts, -1, OW, tempY);
	delete tempY;
  }
}

////////////////////////////////////////////////////////////////
void ampl_get_matrices(  fint irow[], fint kcol[], double elts[],
			int /* nx */, int nnzQ,
			int  my , int nnzA,
			int  mz , int nnzC,
			int irowQ[], int jcolQ[], double dQ[],
			int irowA[], int jcolA[], double dA[],
			int irowC[], int jcolC[], double dC[], 
			double *dwrkX,
			const int full_size )
{			
  int i, j, k;

  // Q
  //call AMPL to get Hessian information
  double *yeltsTemp = (double*)malloc((my+mz)*sizeof(double));
  if(dwrkX)
  	ampl_GetHessianInit(dwrkX,elts,yeltsTemp);
  else
  	ampl_GetHessianInit(NULL,elts,yeltsTemp);
  free(yeltsTemp);

  if(1==full_size){
    int kQ = 0;
    for( j = 0; j < n_var; j++ ) {
      for( k = kcol[j]; k < kcol[j+1]; k++ ) {
        i = irow[k];
        // use only the lower triangle, but transpose
        if( i <= j ) {
	  	  assert( kQ < nnzQ );
	  	  irowQ[kQ] = j; jcolQ[kQ] = i; dQ[kQ] = elts[k];
	  	  kQ++;
        }
      }
    }
  }else{
	int kQ = 0;
	for( j = 0; j < n_var; j++ ) {
	  for( k = kcol[j]; k < kcol[j+1]; k++ ) {
		i = irow[k];
		// use only the upper triangle, but transpose to lower triangle for OOQP
		assert(i<=j);
		irowQ[kQ] = j; jcolQ[kQ] = i; dQ[kQ] = elts[k];
		kQ++;
	  }
	}
	assert(kQ==nnzQ);	
  }


  int *rownbsWrk = (int*)malloc((nnzA+nnzC)*sizeof(int));
  double *elWrk = (double*)malloc((nnzA+nnzC)*sizeof(double));
  dA_goff = (int*)malloc((nnzA)*sizeof(int));
  dC_goff = (int*)malloc((nnzC)*sizeof(int));
  
  /* find row_number of Jac*/	
  cgrad *cg;

  for( j=0; j < n_con; j++){
	for(cg = Cgrad[j]; cg; cg = cg->next){
	  rownbsWrk[cg->goff] = j;		  
	  elWrk[cg->goff] = cg->coef;
	}		
  }


  int kA = 0, kC = 0;
  for( j = 0; j < n_var; j++ ) {
    for( k = A_colstarts[j]; k < A_colstarts[j+1]; k++ ) {
      int iampl = rownbsWrk[k];
      if ( amplRowMap[iampl] < 0 ) { // A negative value in rowMap indicates
		// and equality constraint.
		i = - ( amplRowMap[iampl] + 1 );

		assert( kA < nnzA );
		irowA[kA] = i;  jcolA[kA] = j;  dA[kA] = elWrk[k];
		dA_goff[kA]=k;

		kA++;
      } else {
     	i = amplRowMap[iampl];

	  	assert( kC < nnzC );
	  	irowC[kC] = i;  jcolC[kC] = j;  dC[kC] = elWrk[k];
		dC_goff[kC]=k;

	  	kC++;
      }
    }
  }

  // Transpose everything --- OOQP uses row-wise index, but OOPS and AMPL use col-wise index
  // This version also save the Goff index which can be reused later by calling doubleLexSort_ValOnly
  doubleLexSort_Init( irowA, nnzA, jcolA, dA, 1);
  doubleLexSort_Init( irowC, nnzC, jcolC, dC, 0);

  free(elWrk);
  free(rownbsWrk);
  
}


void ampl_write_solution(double *varsX, double *Yelts,double *Zelts ) 
{
  double *dualWrk = (double*) malloc(n_con*sizeof(double)); 
  int i;
  
  for( int iampl = 0; iampl < n_con; iampl++ ) {  
    if( amplRowMap[iampl]<0 ) {
      i = - ( amplRowMap[iampl] + 1 ); // Recover i from the negative value
      if(Yelts)
	dualWrk[iampl] = Yelts[i];
      else
	dualWrk[iampl] = 0;
    } else {
      i = amplRowMap[iampl];
      if(Zelts)
	dualWrk[iampl] = Zelts[i];
      else
	dualWrk[iampl] = 0;  
    }
  }
  write_sol("\npipsnlp_serial:", varsX, dualWrk,NULL);
  free(dualWrk);
}



void ampl_free_mapinfo(){

	if(amplRowMap) delete [] (amplRowMap);
	if(invpRowMap) delete [] (invpRowMap);

	if(XStatus) delete [] (XStatus);
	if(ConStatus) delete [] (ConStatus);

	if(JacAGoffTransMap) free(JacAGoffTransMap);
	if(JacCGoffTransMap) free(JacCGoffTransMap);

	if(dA_goff) free(dA_goff);
	if(dC_goff) free(dC_goff);

	if(dC_L_goff) free(dC_L_goff);
	if(dC_U_goff) free(dC_U_goff);	

	if(dC_UL_goff) free(dC_UL_goff);
	if(amplRowMap_CLow) delete [] (amplRowMap_CLow);	
	if(amplRowMap_CUp) delete [] (amplRowMap_CUp);	
	if(amplRowMap_XLow) delete [] (amplRowMap_XLow);	
	if(amplRowMap_XUp) delete [] (amplRowMap_XUp);	


	
}