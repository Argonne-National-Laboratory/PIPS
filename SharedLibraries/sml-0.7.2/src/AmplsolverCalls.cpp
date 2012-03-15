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

#include "AmplsolverCalls.h"
#include "asl_pfgh.h"
#include <cassert>
#include <limits>
#include <string>

using namespace std;

const bool log_NL=false;

/** Constructor */
NlFile::NlFile(const string& name) :
  nlfilename(name),
  ncol(-1),
  nrow(-1),
  nzH(-1),
  nzA(-1) {
  asl_pfgh_ptr = (ASL_pfgh*) ASL_alloc(ASL_read_pfgh);
  asl_ptr = NULL;
  readNlFile();
}

/** Destructor */
NlFile::~NlFile() {

  ASL_free((ASL**) &asl_pfgh_ptr);
  if (asl_ptr) {
    ASL *asl = asl_ptr;
    free(A_vals);
    ASL_free(&asl_ptr);
  }

  map<ExpandedModel*, IndexListValue*>::iterator it;
  for (it = indexList.begin(); it != indexList.end(); ++it)
    delete (*it).second;
}

/* ===========================================================================
methods
============================================================================ */
void 
NlFile::readNlFile() {

  if (log_NL) printf("NlFile::readNlFile: (%-30s): ", nlfilename.c_str());
  ASL_pfgh *asl = asl_pfgh_ptr;

  FILE *nl = jac0dim(const_cast<char*>(nlfilename.c_str()), nlfilename.size());
  if (nl==NULL){
    printf("File not found %s\n",nlfilename.c_str());
    exit(1);
  }
  int err = pfgh_read(nl, ASL_return_read_err);
  if (err!=0){
    printf("pfgh_read returns err=%d\n",err);
    exit(1);
  }
  
  ncol = n_var;
  nrow = n_con;
  nzA = nzc;

  if(log_NL) printf("(%dx%d): %d nz\n",n_con, n_var,nzc);
}

/* ----------------------------------------------------------------------------
NlFile::getNoConstraints
---------------------------------------------------------------------------- */
int
NlFile::getNoConstraints(){
  assert(nrow >= 0);
  return nrow;
}

/* ----------------------------------------------------------------------------
NlFile::getNoVariables
---------------------------------------------------------------------------- */
int
NlFile::getNoVariables(){
  assert(ncol >= 0);
  return ncol;
}

/* ----------------------------------------------------------------------------
NlFile::getNoHessianEntries
---------------------------------------------------------------------------- */
int
NlFile::getNoHessianEntries()
{
  if (nzH>=0) return nzH;

  if(log_NL) printf("NlFile::getNzHess:  (%-30s): ",nlfilename.c_str());
  ASL_pfgh *asl = asl_pfgh_ptr;

  nzH = sphsetup(/* nobj=*/0, /*ow =*/1, /*y=*/0, /*uptri=*/0); 

  if(log_NL) printf("%d\n",nzH);
  return nzH;
}

/* ----------------------------------------------------------------------------
NlFile::getHessianStructure
---------------------------------------------------------------------------- */
void
NlFile::getHessianStructure(int *colbegH, int *rownbsH)
{
  if(log_NL) printf("NlFile::getHessStr: (%s)\n",nlfilename.c_str());
  ASL_pfgh *asl = asl_pfgh_ptr;

  nzH = sphsetup(/* nobj=*/0, /*ow =*/1, /*y=*/0, /*uptri=*/0); 
  // this should have setup the sputinfo fields

  for(int i=0;i<=ncol;i++) colbegH[i] = sputinfo->hcolstarts[i];
  for(int i=0;i<nzH;i++) rownbsH[i] = sputinfo->hrownos[i];
  // FIXME: need to copy these since sputinfo will be deallocated by the
  //        call to ASL_free. If we have te asl pointer part of the
  //        class, then we can just keep this information in memory
  //        and pass a pointer back to the caller.
}

/* ----------------------------------------------------------------------------
NlFile::getHessianEntries
---------------------------------------------------------------------------- */
void
NlFile::getHessianEntries(int *colbegH, int *rownbsH, double *eltsH)
{
  if(log_NL) printf("NlFile::getHessian: (%s)\n",nlfilename.c_str());
  ASL_pfgh *asl = asl_pfgh_ptr;

  nzH = sphsetup(/* nobj=*/0, /*ow =*/1, /*y=*/0, /*uptri=*/0); 
  // this should have setup the sputinfo fields
  sphes(eltsH, /*nobj=*/0, NULL, NULL); 

  for(int i=0;i<=ncol;i++) colbegH[i] = sputinfo->hcolstarts[i];
  for(int i=0;i<nzH;i++) rownbsH[i] = sputinfo->hrownos[i];
  // FIXME: need to copy these since sputinfo will be deallocated by the
  //        call to ASL_free. If we have te asl pointer part of the
  //        class, then we can just keep this information in memory
  //        and pass a pointer back to the caller.
}

void
NlFile::readNlFile_f() {

  ASL *asl = asl_ptr = ASL_alloc(ASL_read_f);
  FILE *nl = jac0dim(const_cast<char*>(nlfilename.c_str()), nlfilename.size());
  if (nl == NULL) {
    printf("File not found %s\n", nlfilename.c_str());
    exit(1);
  }

  // FIXME: The convenient column-wise data access is only available for the
  // LP reader ASL_read_f. For ASL_read_pfgh this results in a segmentation
  // fault due to Cgrad not being allocated (Cgrad is not allocated when
  // A_vals is set). Will f_read fall over for a QP problem? In that case we
  // need to rewrite this using the cgrad structures.

  // to say we want column-wise representation
  A_vals = (real*) Malloc(nzc * sizeof(real));

  int err = f_read(nl, ASL_return_read_err);
  if (err != 0) {
    printf("pfgh_read returns error code %d\n", err);
    exit(1);
  }
}

/* ----------------------------------------------------------------------------
getNoNonzerosAMPL
---------------------------------------------------------------------------- */
/** Return the number of nonzeros for a (vertical) slice of the
 *  constraint matrix (Jacobian) defined in this file.
 *
 *  @param nvar
 *         Number of variables (columns) in the slice.
 *  @param lvar
 *         The indices of the variables in the slice.
 *  @return Number of nonzeros in the slice of the jacobian.
 *
 *  For the part of the problem defined by the intersection of all the
 *  constraints in the *.nl file and the variables given by nvar, lvar
 *  this routine will return the nonzeros in the Jacobian.
 */
int 
NlFile::getNoNonzerosAMPL(int nvar, const int *lvar) {

  if(log_NL) printf("NlFile::getNoNz   : (%-30s): ",nlfilename.c_str());

  if (!asl_ptr)
    readNlFile_f();
  ASL *asl = asl_ptr;

  int tt_nz = 0;
  for(int i=0;i<nvar;i++){
    int col = lvar[i];
    // col==-1 indicates that this column is not present in the AMPL file
    if (col>=0)
      tt_nz += A_colstarts[col+1]-A_colstarts[col];
  }

  if (tt_nz<0){
    printf("getNoNozerosAMPL returns tt_nz = %d\n",tt_nz);
    exit(1);
  }
  if(log_NL) printf("%d\n",tt_nz);
  
  return tt_nz;
}

/* ----------------------------------------------------------------------------
fillSparseAMPL
---------------------------------------------------------------------------- */
/** Return a (vertical) slice of the constraint matrix (Jacobian)
 *  defined in this file in (columnwise) sparse matrix format (by
 *  filling in the memory locations provided).
 *
 *  @param[in] nvar
 *             Number of variables (columns) in the slice.
 *  @param[in] lvar
 *             The indices of the variables in the slice.
 *  @param[out] colbeg
 *             Pointer to column starts in rownbs, el.
 *  @param[out] collen
 *             Vector of column lengths (can be NULL).
 *  @param[out] rownbs
 *             Row indices for the sparse elements.
 *  @param[out] el
 *             The actual nonzero elements.
 *  
 *  For the part of the problem defined by the intersection of all the
 *  constraints in the *.nl file and the variables given by nvar, lvar
 *  this routine will return the Jacobian in (columnwise) sparse matrix format.
 */
void
NlFile::fillSparseAMPL(int nvar, const int *lvar,
		       int *colbeg, int *collen, int *rownbs, double *el) {

  if(log_NL) printf("NlFile::fillSparse: (%s)\n",nlfilename.c_str());

  if (!asl_ptr)
    readNlFile_f();
  ASL *asl = asl_ptr;

  int tt_nz = 0;
  for(int i=0;i<nvar;i++){
    int col = lvar[i];
    if (col==-1){
      // this column is not present in the AMPL file => empty column
      colbeg[i] = tt_nz;
      if (collen) collen[i] = 0;
    }else{
      colbeg[i] = tt_nz;
      if (collen) collen[i] = A_colstarts[col+1]-A_colstarts[col];
      for(int j=A_colstarts[col];j<A_colstarts[col+1];j++){
	rownbs[tt_nz] = A_rownos[j];
	el[tt_nz] = (double)A_vals[j];
	tt_nz++;
      }      
    }
  }
  colbeg[nvar] = tt_nz;
}

/* ----------------------------------------------------------------------------
getRowBoundsAMPL
---------------------------------------------------------------------------- */
/** Return the row bounds for the constraints defined in this *.nl file.
 *
 *  @param[out] lower
 *              The row lower bounds as a dense array.
 *  @param[out] upper
 *              The row upper bounds as a dense array.
 */
void
NlFile::getRowBoundsAMPL(double *lower, double *upper) const {

  if (log_NL) printf("NlFile::getRowBnds: (%s)\n", nlfilename.c_str());
  ASL_pfgh *asl = asl_pfgh_ptr;

  const double *elts = LUrhs;
  for (int i = 0; i < nrow; ++i) {
    *lower++ = *elts++;
    *upper++ = *elts++;
  }
}

/* ----------------------------------------------------------------------------
getObjAMPL
---------------------------------------------------------------------------- */
/** Evaluate the objective gradient (linear coefficients) for a
 *  (vertical) slice of the problem stored in the *.nl file.
 *
 *  @param[in] nvar
 *             Number of variables defining the slice.
 *  @param[in] lvar
 *             The indices of the variables defining the slice.
 *  @param[out] elts
 *             The objective gradient vector w.r.t. the variables defined in
 *             nvar/lvar.
 *
 *  @bug This only works for linear objective functions: no vector x at which
 *  the objective should be evaluated is passed in.
 *
 *  @attention This routine evaluates the second last objective
 *  function defined in the *.nl file. Standard AMPL behaviour is to
 *  evaluate the last defined objective function (unless otherwise
 *  specified). Since for the SML generated *.nl file the final
 *  objective function is the dummy objective, this routine will by
 *  default evaluate the second last objective.
 *
 *  @note SML in principle supports the definition of several
 *  objective functions. It is unclear how the user would choose
 *  them. I assume by passing some option into the 
 *  "ExpandedModel *generateSML(....)" function
 */
void 
NlFile::getObjAMPL(int nvar, int *lvar, double *elts)
{
  if(log_NL) printf("NlFile::getObj    : (%s)\n",nlfilename.c_str());
  ASL_pfgh *asl = asl_pfgh_ptr;

  if (n_obj>1){
    real *c = (real *)Malloc(n_var*sizeof(real));
    real objsign = objtype[n_obj-2] ? -1. : 1; //set to -1 for maximise
    
    for(int i=0;i<n_var;i++) c[i] = 0;
    for (ograd *og = Ograd[n_obj-2];og;og=og->next)
      c[og->varno] = objsign*og->coef;
    for(int i=0;i<nvar;i++){
      int ix = lvar[i];
      elts[i] = c[ix];
    }
    free(c);
  }
}

/* ----------------------------------------------------------------------------
getColUpBoundsAMPL
---------------------------------------------------------------------------- */
/** Return upper variable (column) bounds for a selection of the
 *  variables in the *.nl file.
 *
 *  @param[in] nvar
 *             Number of variables defining the slice.
 *  @param[in] lvar
 *             The indices of the variables defining the slice.
 *  @param[out] elts
 *             The upper bounds for the defined variables.
 */
void 
NlFile::getColUpBoundsAMPL(int nvar, int *lvar, double *elts)
{
  if (log_NL) printf("NlFile::getUpBnd  : (%s)\n", nlfilename.c_str());
  ASL_pfgh *asl = asl_pfgh_ptr;

  for(int i=0;i<nvar; i++){
    int ix = lvar[i];
    elts[i] = LUv[2*ix+1];
    if(elts[i] == negInfinity) elts[i] = - numeric_limits<double>::infinity();
    if(elts[i] == Infinity)    elts[i] =   numeric_limits<double>::infinity();
  }
}

/* ----------------------------------------------------------------------------
getColLowBoundsAMPL
---------------------------------------------------------------------------- */
/** Return lower variable (column) bounds for a selection of the
 *  variables in the *.nl file.
 *
 *  @param[in] nvar
 *             Number of variables defining the slice.
 *  @param[in] lvar
 *             The indices of the variables defining the slice.
 *  @param[out] elts
 *             The lower bounds for the defined variables.
 */
void 
NlFile::getColLowBoundsAMPL(int nvar, int *lvar, double *elts)
{
  if (log_NL) printf("NlFile::getLowBnd : (%s)\n", nlfilename.c_str());
  ASL_pfgh *asl = asl_pfgh_ptr;

  for(int i=0;i<nvar; i++){
    int ix = lvar[i];
    elts[i] = LUv[2*ix];
    if(elts[i] == negInfinity) elts[i] = - numeric_limits<double>::infinity();
    if(elts[i] == Infinity)    elts[i] =   numeric_limits<double>::infinity();
  }
}

/* -------------------------------------------------------------------------
NlFile::findIxOfLocalVarsInNlFile
-------------------------------------------------------------------------- */
/** Find the indices of the local variables of a given model in this nlfile.
 *
 *  @param[in] em
 *             The model that defines the local variables.
 *  @param[out] lvar
 *             Assumed to be allocated with em->nLocalVar elements: lvar[i]
 *             is the the index of the ith local variable in em in the nlfile.
 *
 *  @return The number of matches found.
 *
 * The given nlfile will define variables spanning different column blocks.
 * This routine scans through those variables and sees if any of them match
 * the local variables of the model 'em'. For every local variable in em
 * it returns '-1' if it is not found in the nlfile, or the position at 
 * which it is found in the nlfile otherwise. 
 *
 * @bug We should change the way this method is called, so that it can just
 * return the lvar stored on the map. Probably lvar should be returned as a 
 * C++ vector.
 */
int
NlFile::findIxOfLocalVarsInNlFile(ExpandedModel *em, int *lvar)
{
  // call the corresponding routine in the ExpandedModel class
  return em->findIxOfLocalVarsInNlFile(this, lvar);
}
