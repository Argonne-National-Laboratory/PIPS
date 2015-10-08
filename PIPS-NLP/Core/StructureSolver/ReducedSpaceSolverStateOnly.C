/* PIPS-NLP                                                         	*
 * Authors: Nai-Yuan Chiang                      		*
 * (C) 2015 Argonne National Laboratory			*/

#include "ReducedSpaceSolverStateOnly.h"

#include "assert.h"
#include "stdlib.h"

#include "DoubleMatrix.h"
#include "SparseSymMatrix.h"
#include "SparseGenMatrix.h"

#include "RegularizationAlg.h"
#include "NlpGenLinsys.h"

#include "UmfPackSolver.h"

#include "OoqpVector.h"
#include "SimpleVector.h"
#include "SimpleVectorHandle.h"

#include "sLinsysRoot.h"

extern int gOuterSolve;

extern double probGenTime;

#ifndef MIN
#define MIN(a,b) ((a > b) ? b : a)
#endif
#ifndef MAX
#define MAX(a,b) ((a > b) ? a : b)
#endif



ReducedSpaceSolverStateOnly::ReducedSpaceSolverStateOnly()
	: Ax_solver(NULL), 
	  firstCallFlag(1), firstSCsolve(1),fullMatDim(0),
	  Hxx_NNz(0), Hss_NNz(0), 
	  Ax_NNz(0), Tx_NNz(0), fullMatNNz(0),	
	  Hxx_ele(NULL), Hss_ele(NULL),
	  Ax_ele(NULL),Tx_ele(NULL),
	  Hxx_rowBeg(NULL),Hss_rowBeg(NULL),
	  Ax_rowBeg(NULL),Tx_rowBeg(NULL),	
	  Hxx_colIdx(NULL),Hss_colIdx(NULL),
	  Ax_colIdx(NULL),Tx_colIdx(NULL),	
	  Hxx_Full_eleMap(NULL),Hss_Full_eleMap(NULL),
	  Ax_Full_eleMap(NULL),Tx_Full_eleMap(NULL),
	  Hxx_Mat(NULL),Hss_Mat(NULL),
	  Ax_Mat(NULL),
	  Tx_Mat(NULL),
	  Msys(NULL), 
	  Hxx_Dim(0), Hss_Dim(0), Ax_Dim_m(0), Tx_Dim_m(0),	  
	  locns(0), locnx(0), locmy(0)	  
{}

ReducedSpaceSolverStateOnly::ReducedSpaceSolverStateOnly(DoubleMatrix* MatIn, 
					const int localXSize, const int localYSize, const int localZSize)
	: Ax_solver(NULL), 
	  firstCallFlag(1), firstSCsolve(1),fullMatDim(0),
	  Hxx_NNz(0), Hss_NNz(0), 
	  Ax_NNz(0), Tx_NNz(0), fullMatNNz(0),	
	  Hxx_ele(NULL), Hss_ele(NULL),
	  Ax_ele(NULL),Tx_ele(NULL),
	  Hxx_rowBeg(NULL),Hss_rowBeg(NULL),
	  Ax_rowBeg(NULL),Tx_rowBeg(NULL),	
	  Hxx_colIdx(NULL),Hss_colIdx(NULL),
	  Ax_colIdx(NULL),Tx_colIdx(NULL),	
	  Hxx_Full_eleMap(NULL),Hss_Full_eleMap(NULL),
	  Ax_Full_eleMap(NULL),Tx_Full_eleMap(NULL),
	  Hxx_Mat(NULL),Hss_Mat(NULL),
	  Ax_Mat(NULL),
	  Tx_Mat(NULL),
	  Msys(NULL),
	  Hxx_Dim(localXSize), Hss_Dim(localZSize), 
	  Ax_Dim_m(localYSize), Tx_Dim_m(localZSize),
	  locns(localZSize), 
	  locnx(localXSize), locmy(localYSize)
{
  int tempfullDim=0;
  MatIn->getSize(tempfullDim,tempfullDim);
  assert(gOuterSolve >= 3);
  assert(tempfullDim == localXSize+localYSize+localZSize+localZSize);

  fullMatDim = tempfullDim;
  fullMatNNz = ((SparseSymMatrix*)MatIn)->numberOfNonZeros();
  Msys = dynamic_cast<SparseSymMatrix*>(MatIn);

  assert(locnx==locmy);

}



void ReducedSpaceSolverStateOnly::firstCall()
{
  firstCallFlag = 0; 

  //shortcut
  int locmz = locns;

  int findNnz_x=0;

  int* krowbegFullMat 	= Msys->getStorageRef().krowM;
  int* jcolFullMat 		= Msys->getStorageRef().jcolM;
  double* eleFullMat 	= Msys->getStorageRef().M;
  

  //count nozeros in Hxx, Huu, Hux
  Hxx_NNz = krowbegFullMat[locnx];
  Hxx_rowBeg  = (int*) malloc((locnx+1)*sizeof(int));
  Hxx_colIdx 	= (int*)malloc(Hxx_NNz*sizeof(int));
  Hxx_ele 	= (double*)malloc(Hxx_NNz*sizeof(double)); 
  Hxx_Full_eleMap 	= (int*)malloc(Hxx_NNz*sizeof(int));

  Hss_NNz = krowbegFullMat[locnx+locns]-krowbegFullMat[locnx];
  Hss_rowBeg  = (int*) malloc((locns+1)*sizeof(int));
  Hss_colIdx 	= (int*)malloc(Hss_NNz*sizeof(int));
  Hss_ele 	= (double*)malloc(Hss_NNz*sizeof(double)); 
  Hss_Full_eleMap 	= (int*)malloc(Hss_NNz*sizeof(int));

  Ax_NNz = krowbegFullMat[locnx+locns+locmy]-krowbegFullMat[locnx+locns]-locmy;
  Ax_rowBeg  = (int*) malloc((locmy+1)*sizeof(int));
  Ax_colIdx 	= (int*)malloc(Ax_NNz*sizeof(int));
  Ax_ele 	= (double*)malloc(Ax_NNz*sizeof(double)); 
  Ax_Full_eleMap 	= (int*)malloc(Ax_NNz*sizeof(int));

  Tx_NNz = krowbegFullMat[locnx+locns+locmy+locns]-krowbegFullMat[locnx+locns+locmy]-locns-locns;
  Tx_rowBeg  = (int*) malloc((locns+1)*sizeof(int));
  Tx_colIdx 	= (int*)malloc(Tx_NNz*sizeof(int));
  Tx_ele 	= (double*)malloc(Tx_NNz*sizeof(double)); 
  Tx_Full_eleMap 	= (int*)malloc(Tx_NNz*sizeof(int));  

  int findNNZ=0, glbRowStart=0,locRowID=0;
  for(int j=0; j<locnx;j++){
  	Hxx_rowBeg[j] = krowbegFullMat[j];
	for( int k=krowbegFullMat[j];k<krowbegFullMat[j+1];k++){	
	  Hxx_colIdx[findNNZ] 	= jcolFullMat[k];
	  Hxx_ele[findNNZ] 		= eleFullMat[k];
	  findNNZ++;
	}
  }
  Hxx_rowBeg[locnx]=findNNZ;
  assert(findNNZ==Hxx_NNz);

  findNNZ=0;
  glbRowStart=locnx;
  for(int j=glbRowStart; j<glbRowStart+locns;j++){
  	locRowID = j-glbRowStart;
  	Hss_rowBeg[locRowID] = krowbegFullMat[j]-krowbegFullMat[glbRowStart];
	for( int k=krowbegFullMat[j];k<krowbegFullMat[j+1];k++){	
	  Hss_colIdx[findNNZ] 	= jcolFullMat[k]-glbRowStart;
	  Hss_ele[findNNZ] 		= eleFullMat[k];
	  findNNZ++;
	}
  }
  Hss_rowBeg[locns]=findNNZ;
  assert(findNNZ==Hss_NNz);

  findNNZ=0;
  glbRowStart=locnx+locns;
  for(int j=glbRowStart; j<glbRowStart+locmy;j++){
  	locRowID=j-glbRowStart;
  	Ax_rowBeg[locRowID] = krowbegFullMat[j]-krowbegFullMat[glbRowStart]-locRowID;
	for( int k=krowbegFullMat[j];k<krowbegFullMat[j+1]-1;k++){	
	  Ax_colIdx[findNNZ] 	= jcolFullMat[k];
	  Ax_ele[findNNZ] 		= eleFullMat[k];
	  findNNZ++;
	}
  }
  Ax_rowBeg[locmy]=findNNZ;
  assert(findNNZ==Ax_NNz);

  findNNZ=0;
  glbRowStart=locnx+locns+locmy;
  for(int j=glbRowStart; j<glbRowStart+locns;j++){
  	locRowID=j-glbRowStart;
  	Tx_rowBeg[locRowID] = krowbegFullMat[j]-krowbegFullMat[glbRowStart]-locRowID-locRowID;
	for( int k=krowbegFullMat[j];k<krowbegFullMat[j+1]-2;k++){	
	  Tx_colIdx[findNNZ] 	= jcolFullMat[k];
	  Tx_ele[findNNZ] 		= eleFullMat[k];
	  findNNZ++;
	}
  }
  Tx_rowBeg[locmz]=findNNZ;
  assert(findNNZ==Tx_NNz);
  assert(glbRowStart+locns==fullMatDim);
  // locmy+locns is given by dual reg, another locns is given by -I
  assert(Hxx_NNz+Hss_NNz+Tx_NNz+Ax_NNz+locns+locmy+locns==fullMatNNz);


  Hxx_Mat = new SparseSymMatrix( locnx, Hxx_NNz, 
				  Hxx_rowBeg, Hxx_colIdx, Hxx_ele);

  Hss_Mat = new SparseSymMatrix( locns, Hss_NNz, 
				  Hss_rowBeg, Hss_colIdx, Hss_ele);
  
  Ax_Mat = new SparseGenMatrix( locmy, locnx, Ax_NNz, 
				  Ax_rowBeg, Ax_colIdx, Ax_ele);
  
  Tx_Mat = new SparseGenMatrix( locns, locnx, Tx_NNz, 
				  Tx_rowBeg, Tx_colIdx, Tx_ele);

  Ax_solver = dynamic_cast<UmfPackSolver*> (new UmfPackSolver(Ax_Mat));

  Msys->printMatrixInMatlab("kktDiag");
  Ax_Mat->printMatrixInMatlab("Ax");
  Hxx_Mat->printMatrixInMatlab("Hxx");

}


int ReducedSpaceSolverStateOnly::_numericalFact ( )
{
  	
  if(1==firstCallFlag){
	firstCall();
  }
  else{
	  double* eleFullMat  = Msys->getStorageRef().M;
	  int* rowBegFullMat  = Msys->getStorageRef().krowM;

	  //update matrices with new entries
	  memcpy(&Hxx_ele[0], &eleFullMat[0], Hxx_NNz*sizeof(double));
	  memcpy(&Hss_ele[0], &eleFullMat[Hxx_NNz], Hss_NNz*sizeof(double));

	  for(int locRowID=0; locRowID<locmy;locRowID++){
		int k = rowBegFullMat[locRowID+locnx+locns];
		int rowLen = Ax_rowBeg[locRowID+1]-Ax_rowBeg[locRowID];
		memcpy(&Ax_ele[Ax_rowBeg[locRowID]], &eleFullMat[k], rowLen*sizeof(double));
	  }

	  for(int locRowID=0; locRowID<locns;locRowID++){
		int k = rowBegFullMat[locRowID+locnx+locns+locmy];
		int rowLen = Tx_rowBeg[locRowID+1]-Tx_rowBeg[locRowID];
		memcpy(&Tx_ele[Tx_rowBeg[locRowID]], &eleFullMat[k], rowLen*sizeof(double));
	  }  
  }
  
  int negEVal = 0;

  Ax_solver->matrixChanged();

  negEVal = locns+locmy;
  return negEVal;

}




void ReducedSpaceSolverStateOnly::solve ( OoqpVector& rhs_ )
{
  int locmz=locns;
  
  SimpleVector& rhs = dynamic_cast<SimpleVector&>(rhs_);
  double *x_Rhs_ele,*s_Rhs_ele, *y_Rhs_ele, *z_Rhs_ele;

  int VarIDinFull;

  x_Rhs_ele = &(rhs.elements()[0]);
  s_Rhs_ele = &(rhs.elements()[locnx]);
  y_Rhs_ele = &(rhs.elements()[locnx+locns]);
  z_Rhs_ele = &(rhs.elements()[locnx+locns+locmy]);
  
  SimpleVector  *x_Rhs =  new SimpleVector( x_Rhs_ele, locnx );
  SimpleVector  *s_Rhs =  new SimpleVector( s_Rhs_ele, locns );
  SimpleVector  *y_Rhs =  new SimpleVector( y_Rhs_ele, locmy ); 
  SimpleVector  *z_Rhs =  new SimpleVector( z_Rhs_ele, locmz );
  
  // build temp vector
  SimpleVector *x_temp=NULL, *y_temp=NULL, *s_temp=NULL, *z_temp=NULL;
  x_temp = new SimpleVector(locnx);  x_temp->copyFrom(*x_Rhs);
  y_temp = new SimpleVector(locmy);  y_temp->copyFrom(*y_Rhs);

  if(locns!=0){
    s_temp = new SimpleVector(locns);  s_temp->copyFrom(*s_Rhs);
    z_temp = new SimpleVector(locmz);  z_temp->copyFrom(*z_Rhs);
  }

  // solve \del_x  and \del_s 
  x_Rhs->copyFrom(*y_temp);							// r_d 	for x
  if(locns!=0) s_Rhs->copyFrom(*z_temp);			// r_d 	for s
  SaddlePointJacSolve(x_Rhs,s_Rhs);					// d_p =  C^{-1}r_d

  // solve \del_y  and \del_z 
  y_Rhs->copyFrom(*x_temp);							// r_p 	for y
  Hxx_Mat->mult(1,*y_Rhs,-1,*x_Rhs);				// r_p-Hc*delta_p 	for y
  if(locns!=0){ 
  	z_Rhs->copyFrom(*s_temp);						// r_p 	for z
	Hss_Mat->mult(1,*z_Rhs,-1,*s_Rhs);				// r_p-Hc*delta_p	for z
  }
  SaddlePointJacTransSolve(y_Rhs,z_Rhs);			// d_d =  C^{-1}(r_p-Hc*delta_p)



 
  delete (x_temp);
  delete (y_temp);
  if(locns!=0) delete (s_temp);
  if(locns!=0) delete (z_temp);
  delete (z_Rhs);
  delete (y_Rhs);
  delete (s_Rhs);
  delete (x_Rhs);
}


// do C^{-1} where C is the reduced Jac	(accoring to (*))
void ReducedSpaceSolverStateOnly::SaddlePointJacSolve (OoqpVector* rhs_y_, OoqpVector* rhs_z_)
{
  SimpleVector* rhs_y = dynamic_cast<SimpleVector*>(rhs_y_);
  int A_m, A_n;
  Ax_Mat->getSize(A_m,A_n);
  assert(A_m==A_n && A_m==rhs_y->length());

SimpleVector res(rhs_y->n);
res.copyFrom(*rhs_y); 
SimpleVector ress(rhs_z_->n);
ress.copyFrom(*rhs_z_);

  //do A solve
  Ax_solver->solve(*rhs_y);

  // if we have slack variable
  if(locns != 0){
  	SimpleVector* rhs_z = dynamic_cast<SimpleVector*>(rhs_z_);
    Tx_Mat->mult(-1,*rhs_z,1,*rhs_y);
  }

Ax_Mat->mult(1,res,-1,*rhs_y);
Tx_Mat->mult(1,ress,-1,*rhs_y);
ress.axpy(1.0,*rhs_z_);
double infRes = res.infnorm();
if(infRes<ress.infnorm()) infRes=ress.infnorm();
assert(infRes>=0);
if(infRes>0.1) 
  cout << " *** res from LU solver is =  " << infRes << " *** " << endl;

  
}

// do C^{-T} where C is the reduced Jac
void ReducedSpaceSolverStateOnly::SaddlePointJacTransSolve (OoqpVector* rhs_x_, OoqpVector* rhs_s_)
{
  SimpleVector* rhs_x = dynamic_cast<SimpleVector*>(rhs_x_);
  int A_m, A_n;
  Ax_Mat->getSize(A_m,A_n);
  assert(A_m==A_n && A_m==rhs_x->length());


SimpleVector res(rhs_x->n);
res.copyFrom(*rhs_x);
SimpleVector ress(rhs_s_->n);
ress.copyFrom(*rhs_s_);

  // if we have slack variable
  if(locns != 0){
  	SimpleVector* rhs_s = dynamic_cast<SimpleVector*>(rhs_s_);
	Tx_Mat->transMult(1,*rhs_x,1,*rhs_s);
	rhs_s->negate();
  }

  //do A trans solve
  Ax_solver->solveTrans(*rhs_x);

//if(rhs_x->infnorm()>1e8) 
//  cout << " *** sol from LU solver is =  " << rhs_x->infnorm() << " *** " << endl;

#if 0
FILE* outfile = fopen( "AxMat_rowise.txt", "wr" );


		int *krowM = Ax_Mat->krowM();
		int *jcolM = Ax_Mat->jcolM();
		double *M = Ax_Mat->M();
		int nnz = Ax_Mat->numberOfNonZeros();
		int n = rhs_x->n;

		fprintf(outfile,"n_dim_Ax = %d;\n\n", n); 
		fprintf(outfile,"nnz_Ax = %d;\n\n", nnz); 
		
		int findkk=0;
		fprintf(outfile,"rowId_Ax = [ ");
		for(int ii=0;ii<n;ii++)
		  for(int kk=krowM[ii];kk<krowM[ii+1];kk++)
			if(findkk%10==0){
			 fprintf(outfile,"%d, ... \n",ii);
			 findkk++;
			}
			else{
			 fprintf(outfile,"%d,", ii);
			 findkk++;
			}
		fprintf(outfile,"]; \n\n");
	
		findkk=0;
		fprintf(outfile,"colId_Ax = [ ");
		for(int kk=0;kk<nnz;kk++)
			if(findkk%10==0){
			 fprintf(outfile,"%d, ... \n", jcolM[kk]);
					 findkk++;
			}
			else{
			 fprintf(outfile,"%d,", jcolM[kk]);
					 findkk++;
			}
		fprintf(outfile,"]; \n\n");
	
		findkk=0;
		fprintf(outfile,"elts_Ax = [ ");
		for(int kk=0;kk<nnz;kk++)
			if(findkk%10==0){
			  fprintf(outfile,"%5.17g, ... \n", M[kk]);
					 findkk++;
			}
			else{
			  fprintf(outfile,"%5.17g,", M[kk]);
					 findkk++;
			}
		fprintf(outfile,"]; \n\n"); 

		fprintf(outfile,"rhs_Ax = [ ");
		for(int kk=0;kk<n;kk++)
			if(findkk%10==0){
			  fprintf(outfile,"%5.17g, ... \n", res.elements()[kk]);
					 findkk++;
			}
			else{
			  fprintf(outfile,"%5.17g,", res.elements()[kk]);
					 findkk++;
			}
		fprintf(outfile,"]; \n\n"); 

		fclose(outfile);
#endif

#if 0
FILE* outfile2 = fopen( "AxMat_colwise.txt", "wr" );

		UmfPackSolver* AXsolver = (UmfPackSolver*)Ax_solver;
		int *irowM = AXsolver->irowM;
		int *kcolbegM = AXsolver->kcolbegM;
		M = AXsolver->eleM;
		nnz = AXsolver->nnz;
		n = AXsolver->n;

		fprintf(outfile2,"n_dim_Ax = %d;\n\n", n); 
		fprintf(outfile2,"nnz_Ax = %d;\n\n", nnz); 
		
		findkk=0;
		fprintf(outfile2,"colbeg_Ax = [ ");
		for(int ii=0;ii<n+1;ii++)
		  if(findkk%10==0){
		    fprintf(outfile2,"%d, ... \n",kcolbegM[ii]);
			findkk++;
		  }
		  else{
			fprintf(outfile2,"%d,", kcolbegM[ii]);
			findkk++;
		  }
		fprintf(outfile2,"]; \n\n");
	
		findkk=0;
		fprintf(outfile2,"rowId_Ax = [ ");
		for(int kk=0;kk<nnz;kk++)
			if(findkk%10==0){
			 fprintf(outfile2,"%d, ... \n", irowM[kk]);
					 findkk++;
			}
			else{
			 fprintf(outfile2,"%d,", irowM[kk]);
					 findkk++;
			}
		fprintf(outfile2,"]; \n\n");
	
		findkk=0;
		fprintf(outfile2,"elts_Ax = [ ");
		for(int kk=0;kk<nnz;kk++)
			if(findkk%10==0){
			  fprintf(outfile2,"%5.17g, ... \n", M[kk]);
					 findkk++;
			}
			else{
			  fprintf(outfile2,"%5.17g,", M[kk]);
					 findkk++;
			}
		fprintf(outfile2,"]; \n\n"); 

		fprintf(outfile2,"rhs_Ax = [ ");
		for(int kk=0;kk<n;kk++)
			if(findkk%10==0){
			  fprintf(outfile2,"%5.17g, ... \n", res.elements()[kk]);
					 findkk++;
			}
			else{
			  fprintf(outfile2,"%5.17g,", res.elements()[kk]);
					 findkk++;
			}
		fprintf(outfile2,"]; \n\n"); 

		fclose(outfile2);
#endif



Ax_Mat->transMult(1,res,-1,*rhs_x);
Tx_Mat->transMult(1,res,-1,*rhs_s_);
ress.axpy(1.0,*rhs_s_);
double infRes = res.infnorm();
if(infRes<ress.infnorm()) infRes=ress.infnorm();
assert(infRes>=0);
if(infRes>0.1) 
  cout << " *** res from LU solver (trans) is =  " << infRes << " *** " << endl;




  
}


void ReducedSpaceSolverStateOnly::solve(GenMatrix& rhs_in)
{
  DenseGenMatrix &rhs = dynamic_cast<DenseGenMatrix&>(rhs_in);
  int N,NRHS;
  // rhs vectors are on the "rows", for continuous memory
  rhs.getSize(NRHS,N);
  assert(fullMatDim==N);
    
  for (int i = 0; i < NRHS; i++) {
    SimpleVector v(rhs[i],N);
    solve(v);
  } 

}

