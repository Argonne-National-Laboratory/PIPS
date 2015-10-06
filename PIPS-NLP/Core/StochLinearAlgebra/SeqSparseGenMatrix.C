/* PIPS-NLP                                                             
 * Author: Nai-Yuan Chiang
 * (C) 2015 Argonne National Laboratory
 */
 
#include "SeqSparseGenMatrix.h"
#include "DoubleMatrixTypes.h"
#include "StochVector.h"
#include "SimpleVector.h"
#include "sTreeMultiStage.h"

extern int gOuterSolve;

using namespace std;

SeqSparseGenMatrix::SeqSparseGenMatrix(int id_in, int rows, int cols, int nnz, 
					MPI_Comm mpiComm_, multiStageInputTree* in_)
  : id(id_in), mpiComm(mpiComm_), SparseGenMatrix(rows,cols,nnz), parent(NULL),
    iAmDistrib(0), in(in_),child(NULL), parentIsEmpty(true)
{
  if(nnz>0) 
  	isEmpty=false;
  else
  	isEmpty=true;

  if(mpiComm!=MPI_COMM_NULL) {
    int size; MPI_Comm_size(MPI_COMM_WORLD, &size);
    if(size>1) iAmDistrib=1;
  }
}

SeqSparseGenMatrix::SeqSparseGenMatrix(int id_in, int rows, int cols, int nnz, 
					sTreeMultiStage* tree_in, MPI_Comm mpiComm_, multiStageInputTree* in_)
  : id(id_in), mpiComm(mpiComm_), SparseGenMatrix(rows,cols,nnz), currTree(tree_in),parent(NULL),
    iAmDistrib(0), in(in_),child(NULL), parentIsEmpty(true)
{
  if(nnz>0) 
  	isEmpty=false;
  else
  	isEmpty=true;
   
  if(currTree && currTree->parent){
//  	int par_my = currTree->parent->my();
//	int par_nx = currTree->parent->nx();

//	int nzero = in->getCouplingConstraint(id).getNumElements();
//    parent = new SeqSparseGenMatrix(id,rows, par_nx, nzero,currTree->parent,mpiComm_, in);




//	parent->krowM()	= in->getCouplingConstraint(id).getNumElements()
//	parent->jcolM()	= in->getCouplingConstraint(id).getNumElements()
//	parent->M()		= in->getCouplingConstraint(id).getNumElements()
//    memcpy(parent->krowM(),in->getCouplingConstraint(id).getVectorStarts(),(in->getCouplingConstraint(id).getNumRows()+1)*sizeof(int));
//    memcpy(parent->jcolM(),in->getCouplingConstraint(id).getIndices(),in->getCouplingConstraint(id).getNumElements()*sizeof(int));
//    memcpy(parent->M(),in->getCouplingConstraint(id).getElements(),in->getCouplingConstraint(id).getNumElements()*sizeof(double));
	
  }

  if(mpiComm!=MPI_COMM_NULL) {
    int size; MPI_Comm_size(MPI_COMM_WORLD, &size);
    if(size>1) iAmDistrib=1;
  }
}


SeqSparseGenMatrix::~SeqSparseGenMatrix()
{
  if(parent) 
  	delete parent;

  //cout << "~~~~~~~~~~SeqSparseGenMatrix " << mStorage->refs()  << endl;
}

int SeqSparseGenMatrix::isKindOf( int type )
{
  return type == kSeqGenMatrix || type == kGenMatrix;
}

void SeqSparseGenMatrix::scalarMult( double num)
{
  SparseGenMatrix::scalarMult(num); 
  if(parent) parent->scalarMult(num);
}


/* y = beta * y + alpha * this * x */
void SeqSparseGenMatrix::mult( double beta,  OoqpVector& y_,
			   double alpha, OoqpVector& x_ )
{
  if (0.0 == alpha) {
    y_.scale( beta );
    return;
  } else {
    SparseGenMatrix::mult(beta, y_, alpha, x_);
  }
}

void SeqSparseGenMatrix::mult( double beta,  StochVector& y_,
			   double alpha, StochVector& x_ )
{
  StochVector & x = dynamic_cast<StochVector&>(x_);
  StochVector & y = dynamic_cast<StochVector&>(y_);

  //check the tree compatibility
//  if(parent){
//	assert(y.parent && x.parent);
//  }

  if (0.0 == alpha || isEmpty) {
    y.vec->scale( beta );
    return;
  } else {

	if(!isEmpty)
      SparseGenMatrix::mult(beta, *y.vec, alpha, *x.vec);

    if(parent) {
      parent->mult(1.0, y, alpha, *x.parent);
    }
  }
}


void SeqSparseGenMatrix::transMult ( double beta,   OoqpVector& y_,
				 double alpha,  OoqpVector& x_ )
{
  SparseGenMatrix::transMult(beta, y_, alpha, x_);
}


void SeqSparseGenMatrix::transMult ( double beta,   StochVector& y_,
				 double alpha,  StochVector& x_ )
{
  StochVector & x = dynamic_cast<StochVector&>(x_);
  StochVector & y = dynamic_cast<StochVector&>(y_);

#ifdef DEBUG
  //check the tree compatibility
  if(parent){
	assert(y.parent && x.parent);
  }
#endif

  SimpleVector& xvec = dynamic_cast<SimpleVector&>(*x.vec);
  SimpleVector& yvec = dynamic_cast<SimpleVector&>(*y.vec);

  if(!isEmpty)
    SparseGenMatrix::transMult(beta, *y.vec, alpha, *x.vec);

  if(parent) {
    parent->transMult(1.0, *y.parent, alpha, x);
  }

#if 1    
  //preparations for the parallel case
  int iAmSpecial = 1;
  if(iAmDistrib) {
    int rank; MPI_Comm_rank(mpiComm, &rank);
    if(rank>0) iAmSpecial = 0;
  }

  if(iAmDistrib) {
    int locn=yvec.length();
    double* buffer = new double[locn];

    MPI_Allreduce(yvec.elements(), buffer, locn, MPI_DOUBLE, MPI_SUM, mpiComm);
    yvec.copyFromArray(buffer);

    delete[] buffer;
  }
#endif
}



void SeqSparseGenMatrix::transMult ( const int setA, double beta,   StochVector& y_,
				 double alpha,  StochVector& x_, StochVector* End_Par_Pos )
{
  StochVector & x = dynamic_cast<StochVector&>(x_);
  StochVector & y = dynamic_cast<StochVector&>(y_);

  SimpleVector& xvec = dynamic_cast<SimpleVector&>(*x.vec);
  SimpleVector& yvec = dynamic_cast<SimpleVector&>(*y.vec);


  int locnx = currTree->nx();
  int locns = currTree->mz();
  int locmy = currTree->my();
  int locmz = locns;  

  int dummy, n0;
  getSize(dummy, n0);
  
  if(!child && !isEmpty){
    SimpleVector z01 (&yvec[0], n0);

    if(setA){
	  if(gOuterSolve>=3 ){
	    SimpleVector zi3 (&xvec[locnx+locns], locmy);
	    SparseGenMatrix::transMult(1.0, z01, -1.0, zi3);
	  }
	  else{
	    SimpleVector zi2 (&xvec[locnx], 	  locmy);
	    SparseGenMatrix::transMult(1.0, z01, -1.0, zi2);
	  }
    }
    else{
  	  // set C
	  if(gOuterSolve>=3 ){
	    SimpleVector zi4 (&xvec[locnx+locns+locmy], locmz);
	    SparseGenMatrix::transMult(1.0, z01, -1.0, zi4);
	  }
	  else{
	    SimpleVector zi3 (&xvec[locnx+locmy], locmz);
	    SparseGenMatrix::transMult(1.0, z01, -1.0, zi3);
	  }
    }
  }
  
  if(parent && y.parent!=End_Par_Pos) {
    parent->transMult(setA, 1.0, *y.parent, alpha, x, End_Par_Pos);
  }

#if 1    
  //preparations for the parallel case
  int iAmSpecial = 1;
  if(iAmDistrib) {
    int rank; MPI_Comm_rank(mpiComm, &rank);
    if(rank>0) iAmSpecial = 0;
  }

  if(iAmDistrib) {
    int locn=yvec.length();
    double* buffer = new double[locn];

    MPI_Allreduce(yvec.elements(), buffer, locn, MPI_DOUBLE, MPI_SUM, mpiComm);
    yvec.copyFromArray(buffer);

    delete[] buffer;
  }
#endif
}




// goal_Par +=  alpha *currAt*x.vec
void SeqSparseGenMatrix::addtransMultAtTarget ( const int setA,
				 double alpha,  StochVector& x_, SimpleVector* goal_Par )
{

  SimpleVector& xvec = dynamic_cast<SimpleVector&>(*x_.vec);
  SimpleVector& yvec = dynamic_cast<SimpleVector&>(*goal_Par);

  SeqSparseGenMatrix *tempMat = this;
  StochVector  *tempVec = x_.parent;

  while((dynamic_cast<SimpleVector*>(tempVec->vec)) != goal_Par){
	tempMat = tempMat->parent;
	tempVec = tempVec->parent;
  }


  int locnx = currTree->nx();
  int locns = currTree->mz();
  int locmy = currTree->my();
  int locmz = locns;  

  int dummy, n0;
  getSize(dummy, n0);
  
  if(!child && !isEmpty){
    SimpleVector z01 (&yvec[0], n0);

    if(setA){
	  // set A
	  if(gOuterSolve.=3 ){
	    SimpleVector zi3 (&xvec[locnx+locns], locmy);
	    tempMat->transMult(1.0, z01, -1.0, zi3);
	  }
	  else{
	    SimpleVector zi2 (&xvec[locnx], 	  locmy);
	    tempMat->transMult(1.0, z01, -1.0, zi2);
	  }
    }
    else{
  	  // set C
	  if(gOuterSolve>=3 ){
	    SimpleVector zi4 (&xvec[locnx+locns+locmy], locmz);
	    tempMat->transMult(1.0, z01, -1.0, zi4);
	  }
	  else{
	    SimpleVector zi3 (&xvec[locnx+locmy], locmz);
	    tempMat->transMult(1.0, z01, -1.0, zi3);
	  }
    }
  }

}



double SeqSparseGenMatrix::abmaxnorm()
{
  double nrm = 0.0, nrmPar=0.0;

  if(!isEmpty)
    nrm = SparseGenMatrix::abmaxnorm();
  
  if(parent) {
    nrmPar = parent->abmaxnorm();
  }

  return max(nrm,nrmPar);
}


int SeqSparseGenMatrix::numberOfNonZeros()
{
  int  nnz_total=0;

  if(!isEmpty)  
    nnz_total = SparseGenMatrix::numberOfNonZeros();
  
  if(parent) {
    nnz_total += parent->numberOfNonZeros();
  }

  return nnz_total;


}


void SeqSparseGenMatrix::AddParent(SeqSparseGenMatrix* parent_)
{
  parents.push_back(parent_);
}


bool SeqSparseGenMatrix::checkIfParentEmpty()
{
  if(parent){
	parentIsEmpty = false;
  }else{
	parentIsEmpty = true;
  }
  return parentIsEmpty;
}

