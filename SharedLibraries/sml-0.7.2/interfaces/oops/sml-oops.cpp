/* (c) 2008,2009 Jonathan Hogg and Andreas Grothey, University of Edinburgh
 *
 * This file is part of SML.
 *
 * SML is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, using version 3 of the License.
 *
 * SML is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
 * details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program. If not, see http://www.gnu.org/licenses/.
 */
/* This is the OOPS driver for the Structured Modelling Language (SML) */

#include "oops/OopsInterface.h"
#include <cmath>
#include <iostream>
#include <cassert>
#include "OOPSBlock.h"
#include "sml-oops.h"

using namespace std;

// C++ static variables defined in an *. file need to be declared here as well
FILE *printout;
double tt_start, tt_end;

static Algebra *createA(ExpandedModelInterface *em);
static Algebra *createQ(ExpandedModelInterface *em);

/* This NodeId is replaced by OOPSBlock which is a proper class carrying the 
   full information about a node in the OOPS Algebra Tree */

/** This is the identifier needed to fill in a block of the constraint 
 *  matrix */
//typedef struct NodeId_st {
//  const char *amplfile;
//  int nvar;
//  int *lvar;
//} NodeId;

#if 0
/** This is the identifier needed to fill in a block of the Hessian matrix */
typedef struct NodeIdQ_st {
  NlFile *nlfile;  //< the *.nl-file that carries the information about this
  int nrowvar;        //< number of rows that this node has
  int *lrowvar;       //< list of row indices
  int ncolvar;        //< number of columns that this node has
  int *lcolvar;       //< number of column indices
} NodeIdQ;
#endif

void FillRhsVector(Vector *vb);
void FillObjVector(Vector *vc);
void FillUpBndVector(Vector *vu);
void FillLowBndVector(Vector *vl);
void SML_OOPS_upload_sol(ExpandedModelInterface *root, Vector *vx, Vector *vy,
			 Vector *vz);

FILE *globlog = NULL;
const int prtLvl = 0;
const bool writeMPS = false;

void
SML_OOPS_driver(ExpandedModelInterface *root)
{
  Algebra *AlgAug;
  Vector *vb, *vc, *vu, *vl;

  PARALLEL_CODE(
                int    InitPar  = InitLippPar(argc, argv);
                )
  printout = stdout;
  
  HopdmOptions Opt;

  Algebra *A = createA(root);
  Algebra *Q = createQ(root);

  AlgAug = OOPSSetup(A, Q);

  // FIXME: should the stuff below be included in OOPSSetup? 
  //        that would require OOPSSetup to return vectors as well

  vb = new Vector(A->Trow, "vb");
  vc = new Vector(A->Tcol, "vc");
  vu = new Vector(A->Tcol, "vu");
  vl = new Vector(A->Tcol, "vl");

  vc->fillCallBack(FillObjVector);
  vb->fillCallBack(FillRhsVector);
  vu->fillCallBack(FillUpBndVector);
  vl->fillCallBack(FillLowBndVector);

  if (1)
  {
    FILE *mout = fopen("mat.m","w");
    // PrintMatrixMatlab(mout, A, "A");
    // PrintMatrixMatlab(mout, Q, "Q");
    vb->printMatlab(mout, "b");
    vc->printMatlab(mout, "c");
    vu->printMatlab(mout, "bu");
    vl->printMatlab(mout, "bl");
    fclose(mout);
  }

  if (writeMPS)
  {
    FILE *mps_file;
    mps_file = fopen("test.mps","w");
    Write_MpsFile(mps_file, AlgAug, vb,vc, vu, vl, 1, NULL, NULL);
    fclose(mps_file);
  }

  Vector *vx, *vy, *vz;
  PrintOptions Prt(1);

  vx = new Vector(A->Tcol, "vx");
  vy = new Vector(A->Trow, "vy");
  vz = new Vector(A->Tcol, "vz");
  PDProblem Prob(AlgAug, vb, vc, vu, vx, vy, vz);
  Prob.l = vl;
  hopdm(stdout, &Prob, &Opt, &Prt);

  /* hopdm returns the solution vector in Prob->x/y/z
     => need to recurse through the vectors and upload each bit to the 
        corresponding ModelInterface.
  */

  SML_OOPS_upload_sol(root, Prob.x, Prob.y, Prob.z);

#if 0
  FreeAlgebraAlg(A);
  FreeAlgebraAlg(Q);
#endif
}

/* ==========================================================================
Here comes the generation with all subroutines
=========================================================================== */

static Algebra* createBottom(ExpandedModelInterface *diag,
                             ExpandedModelInterface *offdiag);
static Algebra* createRhs(ExpandedModelInterface *diag,
                          ExpandedModelInterface *offdiag);
static Algebra* createBottomQ(ExpandedModelInterface *diag,
                              ExpandedModelInterface *offdiag);
static Algebra* createRhsQ(ExpandedModelInterface *diag,
                           ExpandedModelInterface *offdiag);
static void SMLCallBack(CallBackInterfaceType *cbi);
static void SMLCallBackQ(CallBackInterfaceType *cbi);


/* --------------------------------------------------------------------------
createA
--------------------------------------------------------------------------- */
/* This routine sets up the matrix A from the ExpandedModelInterface tree */
Algebra *
createA(ExpandedModelInterface *em)
{
  Algebra *Alg;

  //if (!A->localVarInfoSet) A->setLocalVarInfo();
  if (em->children.size()==0){
    // this is final node: read in *.nl file to get dimensions 

    
    OOPSBlock *obl = new OOPSBlock(em, em);

    if (em->getNLocalCons()<0){
      printf("CreateA: local block has %d constraints: Not initialised?\n",
             em->getNLocalCons());
      exit(1);
    }
    if (prtLvl>=1){
      cout << "SMLOOPS: Create leaf node: " << em->getName() << ":" << 
         em->getName() << " (" << em->getNLocalCons() << "x" << 
         em->getNLocalVars() << ")" << endl;
    }
    Alg = NewAlgebraSparse(em->getNLocalCons(), em->getNLocalVars(), 
                           (em->getName()+":"+em->getName()).c_str(),
                           (CallBackFunction)SMLCallBack, obl);

  }else{
    /* this is a complex node, set up DblBordDiag with
       - Diagonals from children (call createA again)
       - Off-diagonals with *.nl file from children and col file from parent
       - bottom from this *.nl file and col from the children              */

    if (prtLvl>=1){
      cout << "SMLOOPS: Create complex node: " << em->getName() << ":" << 
         em->getName() << " (" << em->getNLocalCons() << "x" << 
	em->getNLocalVars() << ") nchd = " <<em->children.size()<< endl;
    }

    /* every child is a diagonal block */
    Algebra **D, **B, **R;
    int nblk, i;

    nblk = (em->children).size();
    D = (Algebra **)calloc(nblk+1, sizeof(Algebra *));
    B = (Algebra **)calloc(nblk, sizeof(Algebra *));
    R = (Algebra **)calloc(nblk, sizeof(Algebra *));

    for(i=0; i<nblk; i++){
      D[i] = createA((em->children).at(i));
      B[i] = createBottom(em, (em->children).at(i));
      R[i] = createRhs(em, (em->children).at(i));
    }

    /* The final D[nblk] block is defined by local constraints/variables */
    
    // this is final node: read in *.nl file to get dimensions 
    // I suspect we can just copy in the code from the leaf node case above 

    OOPSBlock *obl = new OOPSBlock(em, em);

    D[nblk] = NewAlgebraSparse(em->getNLocalCons(), em->getNLocalVars(), 
                               (em->getName()+":"+em->getName()).c_str(),
                               (CallBackFunction)SMLCallBack, obl);

    Alg = NewAlgebraDblBordDiag(nblk, B, R, D, 
                                (em->getName()+":"+em->getName()).c_str()); 
  }

  return Alg;
}


Algebra *
createBottom(ExpandedModelInterface *diag, ExpandedModelInterface *nondiag)
{
  Algebra *Alg;
  /* This is a bottom block: 
     take the local constraints from the diag block and
     follow tree defined from the non-diag block */

  if (nondiag->children.size()==0){

    OOPSBlock *obl = new OOPSBlock(diag, nondiag);
    
    Alg = NewAlgebraSparse(diag->getNLocalCons(), nondiag->getNLocalVars(), 
                           (diag->getName()+":"+nondiag->getName()).c_str(),
                           (CallBackFunction)SMLCallBack, obl);

  }else{
    // this is going to be a BlockDense Algebra
    int nblk = nondiag->children.size();
    Algebra **B = (Algebra **)calloc(nblk+1, sizeof(Algebra *));
    
    for(int i=0;i<nblk;i++){
      B[i] = createBottom(diag, (nondiag->children).at(i));
    }

    // The right most block is made up of the diag nodes amplfile
    // and the variables from this nondiag node

    //NodeId *id = new NodeId();
    
    OOPSBlock *obl = new OOPSBlock(diag, nondiag);

    B[nblk] = NewAlgebraSparse(diag->getNLocalCons(),nondiag->getNLocalVars(), 
                           (diag->getName()+":"+nondiag->getName()).c_str(), 
                               (CallBackFunction)SMLCallBack, obl);
    Alg = NewAlgebraBlockDense(1, nblk+1, B, 
                      (diag->getName()+":"+nondiag->getName()).c_str());
  }

  return Alg;
}

Algebra *
createRhs(ExpandedModelInterface *diag, ExpandedModelInterface *nondiag)
{
  Algebra *Alg;
  /* This is a bottom block: 
     take the local variables from the diag block and
     follow tree defined from the non-diag block */

  if (nondiag->children.size()==0){

    OOPSBlock *obl = new OOPSBlock(nondiag, diag);

    Alg = NewAlgebraSparse(nondiag->getNLocalCons(), diag->getNLocalVars(), 
                           (nondiag->getName()+":"+diag->getName()).c_str(),
                           (CallBackFunction)SMLCallBack, obl);

  }else{
    // this is going to be a BlockDense Algebra
    int nblk = nondiag->children.size();
    Algebra **B = (Algebra **)calloc(nblk+1, sizeof(Algebra *));
    
    for(int i=0;i<nblk;i++){
      B[i] = createRhs(diag, (nondiag->children).at(i));
    }
    
    // The bottom node is made from this node's amplfile and the variables
    // defined in diag
    OOPSBlock *obl = new OOPSBlock(nondiag, diag);

    B[nblk] = NewAlgebraSparse(nondiag->getNLocalCons(),diag->getNLocalVars(), 
                           (nondiag->getName()+":"+diag->getName()).c_str(), 
                               (CallBackFunction)SMLCallBack, obl);

    Alg = NewAlgebraBlockDense(nblk+1, 1, B, 
                       (nondiag->getName()+":"+diag->getName()).c_str());
  }

  return Alg;
}

/* --------------------------------------------------------------------------
createQ
--------------------------------------------------------------------------- */
/* This routine sets up the matrix Q from the ExpandedModelInterface tree */
Algebra *
createQ(ExpandedModelInterface *em)
{
  Algebra *Alg;

  //if (!em->localVarInfoSet) em->setLocalVarInfo();
  if (em->children.size()==0){
    // this is final node: read in *.nl file to get dimensions 

#if 0
    NodeIdQ *id = new NodeIdQ();

    id->nlfile = em->nlfile;
    id->nrowvar = em->getNLocalVars();
    id->lrowvar = em->listOfVars;
    id->ncolvar = em->getNLocalVars();
    id->lcolvar = em->listOfVars;
#endif

    Alg = NewAlgebraSparse(em->getNLocalVars(), em->getNLocalVars(), 
                           ("Q"+em->getName()+":"+em->getName()).c_str(),
                           (CallBackFunction)SMLCallBackQ, em);
      
      
  }else{
    // Ok, should really test whether there is a DblBordDiagMatrix needed
    // or if a BlockDiagMatrix will do:
    // - Are there any cross products in the objective bewteen the 
    //   variables local to this node and the variables defined in its 
    //   children?

    /* If there are cross-products then the variables of the children
       that are refered to are part of the AMPL model at this node, but
       they are not part of the nvar/lvar list of local variables 

       AMPL will give us 
         sputinfo->hcolstarts
         sputinfo->hrownos
       for the Hessian defined in this model file. This could be either
       upper triangle only or the full matrix

       We can scan this for cross products (i.e. go through the columns
       corresponding to the local variables are see if there are any entries
       in rows corresponding to children). 
       If so then set up a DblBordDiagMatrix and let the child work out
       which of these entries are his.
       The child gets passed both ExpandedModelInterfaces (for itself and the parent)
       The objective however MUST be included in the parent part

    */
    
    /* this is a complex node, set up DblBordDiag with
       - Diagonals from children (call createQ again)
       - Off-diagonals with *.nl file from children and col file from parent
       - bottom from this *.nl file and col from the children              */


#if 0
    int nzH = (em->nlfile)->getNoHessianEntries();
    int colH = (em->nlfile)->getNoVariables(); 
    int *colbegH = (int *)malloc((colH+1)*sizeof(int));
    int *rownbsH = (int *)malloc(nzH*sizeof(int));
    (em->nlfile)->getHessianStructure(colbegH, rownbsH);

    // now scan to see of there are any off-diagonal entries

    // set marker array that indicates the local variables
    int *marker = (int*)calloc(colH, sizeof(int));
    for(int i=0;i<em->getNLocalVars();i++) marker[em->listOfVars[i]]=1;

    bool foundCross = false;
    for(int i=0;i<em->getNLocalVars();i++){
      int ix = em->listOfVars[i];
      // and scan through the Hessian structure of this row
      for(int j=colbegH[ix];j<colbegH[ix+1];j++){
        int row = rownbsH[j];
        if (marker[row]==0){
          // this is an entry in a row that does not belong to this node
          // (but presumably to a child node)
          foundCross = true;
        }
      }
    }
    free(colbegH); free(rownbsH);
    free(marker);
    
    if (foundCross){
      /* every child is a diagonal block */
      Algebra **D, **B, **R;
      int nblk, i;
      
      nblk = (em->children).size();
      D = (Algebra **)calloc(nblk+1, sizeof(Algebra *));
      B = (Algebra **)calloc(nblk, sizeof(Algebra *));
      R = (Algebra **)calloc(nblk, sizeof(Algebra *));
      
      for(i=0; i<nblk; i++){
        D[i] = createQ((em->children).at(i));
        B[i] = createBottomQ(em, (em->children).at(i));
        R[i] = createRhsQ(em, (em->children).at(i));
      }

      /* The final D[nblk] block is defined by local constraints/variables */
    
      // this is final node: read in *.nl file to get dimensions 
      // I suspect we can just copy in the code from the leaf node case above 

      NodeIdQ *id = new NodeIdQ();

      id->nlfile = em->nlfile;
      id->ncolvar = em->getNLocalVars();
      id->lcolvar = em->listOfVars;
      id->nrowvar = em->getNLocalVars();
      id->lrowvar = em->listOfVars;

      D[nblk] = NewAlgebraSparse(em->getNLocalVars(), em->getNLocalVars(), 
                            ("Q"+em->getName()+":"+em->getName()).c_str(),
                                 (CallBackFunction)SMLCallBackQ, id);

      Alg = NewAlgebraDblBordDiag(nblk, B, R, D, 
                             ("Q"+em->getName()+":"+em->getName()).c_str()); 

    }else{ // Not foundCross => setup BlockDiagMatrix
#endif
      /* every child is a diagonal block */
      Algebra **D;
      int nblk, i;
      
      nblk = (em->children).size();
      D = (Algebra **)calloc(nblk+1, sizeof(Algebra *));
      
      for(i=0; i<nblk; i++){
        D[i] = createQ((em->children).at(i));
      }

      /* The final D[nblk] block is defined by local constraints/variables */
    
      // this is final node: read in *.nl file to get dimensions 
      // I suspect we can just copy in the code from the leaf node case above 

#if 0
      NodeIdQ *id = new NodeIdQ();

      id->nlfile = em->nlfile;
      id->ncolvar = em->getNLocalVars();
      id->lcolvar = em->listOfVars;
      id->nrowvar = em->getNLocalVars();
      id->lrowvar = em->listOfVars;

      D[nblk] = NewAlgebraSparse(em->getNLocalVars(), em->getNLocalVars(), 
                              ("Q"+em->getName()+":"+em->getName()).c_str(),
                                 (CallBackFunction)SMLCallBackQ, id);
#else
      D[nblk] = NewAlgebraSparse(em->getNLocalVars(), em->getNLocalVars(), 
                              ("Q"+em->getName()+":"+em->getName()).c_str(),
                                 (CallBackFunction)SMLCallBackQ, em);
#endif

      Alg = NewAlgebraBlockDiag(nblk+1, D, 
                        ("Q"+em->getName()+":"+em->getName()).c_str()); 
#if 0      
    }
#endif
  }

  return Alg;
}

#if 0
/* ----------------------------------------------------------------------------
createBottom/RhsQ
---------------------------------------------------------------------------- */
/** These two functions will setup the Rhs/Bottom blocks of the Q matrix
 *  Since Q is symmetric, these two functions should be almost identical
 *  They are called if a parent *.nl file finds cross products in the 
 *  objective (i.e. in a column belonging to a parent variable there 
 *  is a entry in a row that does not belong to the parent - but presumably
 *  to a child.
 *  This child should look at the Hessian structure setup in the parents
 *  *.nl file (passed in diag). It should scan through the columns belonging
 *  to the parents variables and see if there are any entries in rows
 *  belonging to *this* child (names of these rows can be obtained from 
 *  'offdiag' they probably need scanning against the column name file
 *  diag->getName()+".col", to get their index numbers in the parents 
 *  numbering.
 *  Also need some treatment of complicated blocks (i.e. when the child node
 *  'offdiag' itself has children) - this should jst be a matter of 
 *  setting up the appropriate BlockDense constructor in the same manner
 *  as done for the A matrix
 *
 *  FIXME: can any of this information (blocking of the Q matrix be
 *         included as part of the ExpandedModelInterface (in order to separate
 *         frontend and backend and do as much processing in the frontend
 *         as possible?
 */

Algebra *
createBottomQ(ExpandedModelInterface *diag, ExpandedModelInterface *offdiag)
{
  printf("createBottomQ not implemented yet!\n");
  exit(1);
}

Algebra *
createRhsQ(ExpandedModelInterface *diag, ExpandedModelInterface *offdiag)
{
  printf("createBottomQ not implemented yet!\n");
  exit(1);
}
#endif

/* ---------------------------------------------------------------------------
CallBackFunction: SMLCallBack
---------------------------------------------------------------------------- */
void
SMLCallBack(CallBackInterfaceType *cbi)
{
  /* This needs to be able to fill in the local sparse nodes with the 
     information coming from the ampl file. It can be called either
     for the diagonal nodes or for the off-diagonal nodes

     It needs to know
     - name of the ampl file
     - list and number of variables to use

  */
  
  /*
   where CallBackInterface is a struct of the following form
    int nz
    int max_nz
    int *row_nbs
    int *col_beg
    int *col_len
    double *element
    void *id
  */
  
  /* id is a pointer to NodeId with components
      string amplfile;
      int nvar;
      int *lvar;
  */


  OOPSBlock *obl = (OOPSBlock*)cbi->id;

  //NodeId *id = (NodeId*)cbi->id;
  if (cbi->row_nbs==NULL){
    // only want number of nonzeros back
    cbi->nz = obl->emrow->getNzJacobianOfIntersection(obl->emcol);
    //assert(nz==cbi->nz);
  }else{
    // want to fill in matrices
    obl->emrow->getJacobianOfIntersection(obl->emcol, cbi->col_beg,
		 cbi->col_len, cbi->row_nbs, cbi->element);
  }
}

/* ---------------------------------------------------------------------------
CallBackFunction: SMLCallBack
---------------------------------------------------------------------------- */
void
SMLCallBackQ(CallBackInterfaceType *cbi)
{
  /* This needs to be able to fill in the local sparse nodes with the 
     information coming from the ampl file. It can be called either
     for the diagonal nodes or for the off-diagonal nodes

     It needs to know
     - name of the ampl file
     - list and number of variables to use

  */
  
  /*
   where CallBackInterface is a struct of the following form
    int nz
    int max_nz
    int *row_nbs
    int *col_beg
    int *col_len
    double *element
    void *id
  */
  
  /* id is a pointer to NodeId with components
      string amplfile;
      int nvar;
      int *lvar;
  */

  
  /* need to know: 
     - the amplfile (better the NlFile structure 
        (the one that has the info on the objective = the diagonal file)
     - the variable list for the diagonal part
     - the variable list for the nondiagonal part
  */
#if 0
  NodeIdQ *id = (NodeIdQ*)cbi->id;

  if (cbi->row_nbs==NULL){
    // only want number of nonzeros back
    NlFile *nlfile = id->nlfile;
    int nzH = nlfile->getNoHessianEntries();
    int colH = nlfile->getNoVariables(); 
    int *colbegH = (int *)malloc((colH+1)*sizeof(int));
    int *rownbsH = (int *)malloc(nzH*sizeof(int));
    nlfile->getHessianStructure(colbegH, rownbsH);
    
    // mark all the relevant rows
    int *marker = (int *)calloc(nlfile->getNoVariables(), sizeof(int));
    for(int i=0;i<nlfile->getNoVariables();i++) marker[i] = 0;
    for(int i=0;i<id->nrowvar;i++) marker[id->lrowvar[i]] = 1;

    // and go through all relevant columns and count the number of entries
    cbi->nz = 0;
    for(int i=0;i<id->ncolvar;i++){
      int col = id->lcolvar[i];
      // scan through this column
      for(int j=colbegH[col];j<colbegH[col+1];j++){
        int row = rownbsH[j];
        if (marker[row]==1) cbi->nz++;
      }
    }
    free(rownbsH); free(colbegH);
    free(marker);
    return;
  }else{
    // only want number of nonzeros back
    NlFile *nlfile = id->nlfile;
    int nzH = nlfile->getNoHessianEntries();
    int colH = nlfile->getNoVariables(); 
    int *colbegH = (int *)malloc((colH+1)*sizeof(int));
    int *rownbsH = (int *)malloc(nzH*sizeof(int));
    double *eltsH = (double *)malloc(nzH*sizeof(double));
    nlfile->getHessianEntries(colbegH, rownbsH, eltsH);
    
    // mark all the relevant rows
    int *marker = (int *)calloc(nlfile->getNoVariables(), sizeof(int));
    for(int i=0;i<nlfile->getNoVariables();i++) marker[i] = 0;
    for(int i=0;i<id->nrowvar;i++) marker[id->lrowvar[i]] = i+1;

    // and go through all relevant columns and copy the entries
    cbi->nz = 0;
    for(int i=0;i<id->ncolvar;i++){
      int col = id->lcolvar[i];
      cbi->col_beg[i] = cbi->nz;
      // scan through this column
      for(int j=colbegH[col];j<colbegH[col+1];j++){
        int row = rownbsH[j];
        if (marker[row]!=0) {
          cbi->element[cbi->nz] = eltsH[j];
          cbi->row_nbs[cbi->nz] = marker[row]-1; // to get the local numbering
          cbi->nz++;
        }
      }
    }
    free(colbegH); free(rownbsH); free(eltsH);
    free(marker);
    cbi->col_beg[id->ncolvar] = cbi->nz;
    for(int i=0;i<id->ncolvar;i++) 
      cbi->col_len[i] = cbi->col_beg[i+1]-cbi->col_beg[i];
    return;
  }
#else
  ExpandedModelInterface *em = (ExpandedModelInterface*)cbi->id;
  cbi->nz = 0;
  if (cbi->row_nbs){
    cbi->col_beg[em->getNLocalVars()] = cbi->nz;
    for(int i=0;i<em->getNLocalVars();i++) 
      cbi->col_len[i] = cbi->col_beg[i+1]-cbi->col_beg[i];
  }
#endif

}

/* ---------------------------------------------------------------------------
CallBackFunction: FillRhsVector
---------------------------------------------------------------------------- */
void
FillRhsVector(Vector *vb)
{
  Tree *T = vb->node;
  DenseVector *dense = GetDenseVectorFromVector(vb);

  double *checkub = new double[dense->dim];

  Algebra *A = (Algebra*)T->nodeOfAlg; // the diagonal node that spawned this tree
  OOPSBlock *obl = (OOPSBlock*)A->id; // and its id structure
  ExpandedModelInterface *emrow = obl->emrow;

  // FIXME: should the id structure include information on the ExpandedModelInterface
  //        as well? That way we could do some more sanity checks

  emrow->getRowBounds(dense->elts, checkub);

  // check that lower and upper constraint bounds are the same due to the 
  // OOPS restriction
  for(int i=0;i<dense->dim; i++){
    if (fabs(dense->elts[i]-checkub[i])>1e-6){
      cerr << "At the moment OOPS only supports equality constraints!\n";
      cerr << "Bounds for c/s " << i << " in " << emrow->getName() << ": " <<
             dense->elts[i] << " " <<  checkub[i] << endl;
      exit(1);
    }
  }
      
  delete [] checkub;
}

/* ---------------------------------------------------------------------------
CallBackFunction: FillObjVector
---------------------------------------------------------------------------- */
void
FillObjVector(Vector *vc)
{
  Tree *T = vc->node;
  DenseVector *dense = GetDenseVectorFromVector(vc);

  Algebra *A = (Algebra*)T->nodeOfAlg; // the diagonal node that spawned this tree
  OOPSBlock *obl = (OOPSBlock*)A->id;        // and its id structure
  //NodeId *id = (NodeId*)A->id;        // and its id structure

  assert(obl->nvar==T->end-T->begin);

  //double *test = new double[obl->nvar];
  //for(int i=0;i<obl->nvar;i++) test[i] = 0.;

  //obl->emrow->nlfile->getObjAMPL(obl->nvar, obl->lvar, test);
  //for (int i=0;i<obl->nvar;i++) printf("SO:%d %d\n",i,obl->lvar[i]);
  obl->emrow->getObjGradient(dense->elts);
  
  //for(int i=0;i<obl->nvar;i++){
  //  if (test[i]!=dense->elts[i]) {
  //    printf("OOPS!: %g vs %g\n",test[i], dense->elts[i]);
  //    exit(1);
  //  }
  //}

  //delete [] test;
}

/* ---------------------------------------------------------------------------
CallBackFunction: FillUpBndVector
---------------------------------------------------------------------------- */
void
FillUpBndVector(Vector *vu)
{
  Tree *T = vu->node;
  DenseVector *dense = GetDenseVectorFromVector(vu);

  Algebra *A = (Algebra *)T->nodeOfAlg; // the diagonal node that spawned this tree
  //NodeId *id = (NodeId*)A->id;        // and its id structure
  OOPSBlock *obl = (OOPSBlock*)A->id;        // and its id structure
  assert(obl->emrow==obl->emcol); // this should be a diagonal block
  //NlFile *nlf = obl->emrow->nlfile;
  ExpandedModelInterface *emrow = obl->emrow;

  assert(obl->nvar==T->end-T->begin);
  //nlf->getColUpBoundsAMPL(obl->nvar, obl->lvar, test);
  emrow->getColUpBounds(dense->elts);
}

/* ---------------------------------------------------------------------------
CallBackFunction: FillLowBndVector
---------------------------------------------------------------------------- */
void
FillLowBndVector(Vector *vl)
{
  Tree *T = vl->node;
  DenseVector *dense = GetDenseVectorFromVector(vl);

  Algebra *A = (Algebra *)T->nodeOfAlg; // the diagonal node that spawned this tree
  OOPSBlock *obl = (OOPSBlock*)A->id;        // and its id structure
  assert(obl->emrow==obl->emcol); // this should be a diagonal block
  ExpandedModelInterface *emrow = obl->emrow;


  assert(obl->nvar==T->end-T->begin);
  emrow->getColLowBounds(dense->elts);

  //for(int i=0;i<dense->dim;i++){
  //  if (fabs(dense->elts[i])>1e-6) {
  //    printf("Found lower bound !=0 (=%f) in variable %i in model %s",
  //           dense->elts[i], i, emrow->getName().c_str());
  //    printf("Currently OOPS can only cope with zero lower bounds\n");
  //    exit(1);
  //  }
  // }

}

/* ---------------------------------------------------------------------------
SML_OOPS_upload_sol
---------------------------------------------------------------------------- */
void
SML_OOPS_upload_sol(ExpandedModelInterface *root,
		    Vector *vx, Vector *vy, Vector *vz)
{
  Tree *Tx = vx->node,*Ty = vy->node;
  int nchd = root->children.size();
  
  //printf("%d: %d %d\n",nchd, Tx->nb_sons, Ty->nb_sons);
  assert((nchd==0&&Tx->nb_sons==0)||Tx->nb_sons==nchd+1);
  assert((nchd==0&&Ty->nb_sons==0)||Ty->nb_sons==nchd+1);
  if (nchd>0){
    /* The final child of the vector tree corresponds to the local variables/
       constraints of this node */
    
    /* upload local contributions */
    Vector *vxs = SubVector(vx, nchd);
    Vector *vys = SubVector(vy, nchd);
    Vector *vzs = SubVector(vz, nchd);

    DenseVector *dx = GetDenseVectorFromVector(vxs);
    DenseVector *dy = GetDenseVectorFromVector(vys);
    DenseVector *dz = GetDenseVectorFromVector(vzs);

    assert(dx->dim == root->getNLocalVars());
    assert(dy->dim == root->getNLocalCons());

    root->setPrimalSolColumns(dx->elts);
    root->setDualSolColumns(dz->elts);
    root->setDualSolRows(dy->elts);

    /* and upload a vector of zeros for the constraint slacks 
       (OOPS can only deal with equality constraints)                    */
    double *elts = new double[dy->dim];
    for(int i=0;i<dy->dim;i++) elts[i] = 0;
    root->setPrimalSolRows(elts);
    delete[] elts;

    /* recurse down the rest of the tree */

    for (int i=0;i<nchd;i++){
      ExpandedModelInterface *model = root->children[i]; 
      SML_OOPS_upload_sol(model, 
			  SubVector(vx, i), SubVector(vy, i), SubVector(vz,i));
    }
  }else{
    /* This is a root node of the model tree and the vector trees */
    DenseVector *dx = GetDenseVectorFromVector(vx);
    DenseVector *dy = GetDenseVectorFromVector(vy);
    DenseVector *dz = GetDenseVectorFromVector(vz);

    assert(dx->dim == root->getNLocalVars());
    assert(dy->dim == root->getNLocalCons());

    root->setPrimalSolColumns(dx->elts);
    root->setDualSolColumns(dz->elts);
    root->setDualSolRows(dy->elts);

    /* and upload a vector of zeros for the constraint slacks 
       (OOPS can only deal with equality constraints)                    */
    double *elts = new double[dy->dim];
    for(int i=0;i<dy->dim;i++) elts[i] = 0;
    root->setPrimalSolRows(elts);
    delete[] elts;
  }
}
