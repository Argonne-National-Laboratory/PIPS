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

#include "ExpandedModel.h"
#include "AmplModel.h"
#include "AmplsolverCalls.h"
#include "GlobalVariables.h"
#include <cassert>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <list>

using namespace std;

#if 0
  #define LogEM(X) cout << X
#else
  #define LogEM(X)
#endif

list<string> ExpandedModel::pathToNodeStack;

/** Constructor */
ExpandedModel::ExpandedModel(AmplModel *src_model) :
   nLocalVars(-1), nLocalCons(-1), localVarInfoSet(false), pvar(NULL), 
   dvar(NULL), prow(NULL), drow(NULL), src(src_model),
   nlfile(NULL) {
  // Note: Setup of model is largely done in AmplModel::createExpandedModel()
}

static unsigned
findIdxOfUncommonChar(bool isRoot, const string& line, const string& name) {
  unsigned int j = 0;
  if (!isRoot && line != "") {
    for (j = 0; j < name.size() && j < line.size(); ++j)
      if (name.at(j) != line.at(j))
        break;
    if (j < line.size() && line.at(j) == '_') ++j;
  }
  return j;
}

/* --------------------------------------------------------------------------
ExpandedModel::setLocalVarInfo
-------------------------------------------------------------------------- */
/** This routine identifies the indices of the local variables in the 
 *  *.nl file that is associated with this ExpandedModel node.
 * 
 * The routine compares the variable name stubs stored in localVarDef
 * (which are set up from the constructor using the corresponding flat
 * model (AmplModel) and the instance (defined by appropriate values
 * for the indexing variables), with the names of the variables
 * defined in the *.nl file (obtained by reading the corresponding
 * *.col file).
 *
 * The routine sets the nLocalVar, listLocalVar, nLocalCons fields of
 * the object.
 *
 * @note The method works by comparing all stubs with all names in the
 * *.col file, this uses a lot of string comparisons and is likely
 * very inefficient.
 */
void 
ExpandedModel::setLocalVarInfo()
{
  if (localVarInfoSet)
    return;

  // FIXME: pretty unelegant to have to do n*m comparisons

  // ------- read the names of constraints defined in this NlFile ------------
  // These have been written in the correct order
  const string nlrowfile = model_file + ".row";
  ifstream fin(nlrowfile.c_str());

  if (!fin) {
    cerr << "Cannot open row name file " << nlrowfile << ".\n";
    exit(1);
  }
  
  list<SymbolTable::Entry> objList = getObjList();
  string model_name = getName();
  bool isRoot = (model_name == "root");
  if (!isRoot) {
    assert(model_name.compare(0, 5, "root_") == 0);
    model_name = model_name.substr(5); // skip "root_"
  }
  string line;
  while(!fin.eof()){
    getline(fin, line);
    listOfConNames.push_back(line);
    // Trim away common characters between model name and var from varname
    unsigned int j = findIdxOfUncommonChar(isRoot, line, model_name);
    // skip empty line
    if(line.substr(j) == "") continue;
    // skip dummy row
    if(line.substr(j,5) == "dummy") {
      if (line.size() == j + 5 || line.at(j + 5) == '[') continue;
    }
    // skip any objectives
    bool isObj = false;
    for(list<SymbolTable::Entry>::const_iterator i=objList.begin(); i!=objList.end(); ++i) {
      const string id = i->id();
      if (line.substr(j, id.size()) == id) {
        if (line.size() == j + id.size() || line.at(j + id.size() == '[')) {
             isObj=true;
             break;
          }
       }
    }
    if(isObj) continue;
    listOfLocalConNames.push_back(line.substr(j));
  }

  fin.close();
  fin.clear(); // Recommended by Marco
  
  if (GlobalVariables::prtLvl >= PRINT_VERBOSE)
    cout << "Read " << listOfConNames.size() << " lines from "
         << nlrowfile << ".\n";

  // Now read vars - these actually need matching to correct order
  const string nlcolfile = model_file + ".col";
  fin.open(nlcolfile.c_str());
  if (!fin) {
    cerr << "Cannot open column name file " << nlcolfile << ".\n";
    exit(1);
  }

  list<string> colfilelist;
  while(!fin.eof()){
    getline(fin, line);
    if (line.size() > 0)
      colfilelist.push_back(line);
  }
  fin.close();
  
  if (GlobalVariables::prtLvl >= PRINT_VERBOSE)
    cout << "Read " << colfilelist.size() << " lines from "
         << nlcolfile << ".\n";

  // -------------- compare this list against the given VarDefs
  list<string>::const_iterator p, q;
  for (p = localVarDef.begin(); p != localVarDef.end(); ++p) {
    // for all variable declarations in ExpandedModel

    int len = (*p).size();
    int cnt;
    if (GlobalVariables::prtLvl >= PRINT_VERBOSE)
      cout << "Trying to match variable definition " << (*p) << ":\n";

    for (q = colfilelist.begin(), cnt = 0; q != colfilelist.end(); ++q, ++cnt) {

      // avoid trying to find impossible matches
      if ((*q).size() < len)
        continue;

      // If *q is longer than *p, it can represent the same variable as *p
      // only if it contains a '[' in its first len characters. Without this
      // check we would match 'slackTotal' with 'slack' (see test-01.mod),
      // which would produce an incorrect variable count and a repeated
      // 'slackTotal' entry in the solution file.
      if ((*q).size() > len && (*q).find('[') > len)
        continue;

      // compare the first 'len' characters from q with p
      if ((*q).compare(0, len, *p) == 0) {
        if (GlobalVariables::prtLvl >= PRINT_VERBOSE)
          cout << "  " << *q << " matches " << *p << "\n";
        listOfLocalVars.push_back(cnt);
        listOfVarNames.push_back(*q);

        // Trim away common characters between model name and var from varname
        unsigned int j = findIdxOfUncommonChar(isRoot, *p, model_name);
        listOfLocalVarNames.push_back((*q).substr(j));
      }
    }
  }
  
  // FIXME: cannot sort and uniquify just the numerical list, since these
  //        actions would need to be mirrored on the name list

  nLocalVars = listOfLocalVars.size();
  //printf("Found %d variables\n", nLocalVars);

  // open *.nl file to get number of constraints
  nLocalCons = nlfile->getNoConstraints();
  //printf("Found %d constraints\n",nLocalCons);

  localVarInfoSet = true;

  if (GlobalVariables::prtLvl >= PRINT_LOG)
    cout << "setLocalVarInfo(): " << model_file << " (" << nLocalCons <<
      "x" << nLocalVars << ")\n";
}

/* --------------------------------------------------------------------------
ExpandedModel::getLocalVarNames
-------------------------------------------------------------------------- */
/** Return the names of the variables local to this node.
 *
 *  @return Names of the local variables.
 */
const list<string>&
ExpandedModel::getLocalVarNames() const {
  assert(localVarInfoSet);
  return listOfLocalVarNames;
}

/* --------------------------------------------------------------------------
ExpandedModel::getLocalConNames
-------------------------------------------------------------------------- */
/** Return the names of the constraints local to this node.
 *
 *  @return Names of local constraints.
 */
const list<string>&
ExpandedModel::getLocalConNames() const {
  assert(localVarInfoSet);
  return listOfLocalConNames;
}

/* --------------------------------------------------------------------------
ExpandedModel::getNLocalVars
-------------------------------------------------------------------------- */
/** Return the number of variables local to this node
 *
 *  @return Number of local variables.
 */
int
ExpandedModel::getNLocalVars() const {
  assert(localVarInfoSet);
  return nLocalVars;
}

/* --------------------------------------------------------------------------
ExpandedModel::getNLocalCons
-------------------------------------------------------------------------- */
/** Return the number of constraints local to this node.
 *
 *  @return Number of local constraints.
 */
int
ExpandedModel::getNLocalCons() const {
  assert(localVarInfoSet);
  return nLocalCons;
}

/* ----------------------------------------------------------------------------
ExpandedModel::print()
---------------------------------------------------------------------------- */
/** Recursively print the contents of this instance and of its children.
 *
 *  Used only for debugging.
 */
void
ExpandedModel::print() const {

  list<string>::const_iterator p;
  list<int>::const_iterator i;

  cout << "EM: ------------------------------------------------------------\n";
  if (nlfile==NULL){
    cerr << "EM: No NlFile attached\n";
    exit(1);
  }
  cout << "EM: This is ExpandedModel: " <<  nlfile->nlfilename << "\n";
  cout << "EM: Nb local variable definitions: " << localVarDef.size() << "\n";
  for (p = localVarDef.begin(); p != localVarDef.end(); ++p)
    cout << "EM:   " << *p << "\n";
  if (!localVarInfoSet){
    cout << "EM: Further information on local variables not set\n";
  }else{
    cout << "EM: Nb local Variables: " << nLocalVars << "\n";
    for (i = listOfLocalVars.begin(); i != listOfLocalVars.end(); ++i)
      cout << *i << " ";
    cout << "\n";
    cout << "EM: Nb local Constraints: " <<  nLocalCons << "\n";
  }

  cout << "EM: Nb children: " << children.size() << "\n";
  if (children.size()>0)
    cout << "EM: now list the children:\n";
  
  for(unsigned int i=0; i<children.size(); ++i){
    ExpandedModel *em = (ExpandedModel*) children.at(i);
    em->print();
  }
}

/* ----------------------------------------------------------------------------
ExpandedModel::getNzJacobianOfIntersection
---------------------------------------------------------------------------- */
/** Return the number of nonzeros in the Jacobian of a section of the model.
 *
 *  The matrix is defined by the intersection of the local constraints in
 *  this model with the local variables of this or another model.
 *
 *  @param emcol_
 *         The model w.r.t. whose local variables the Jacobian should be
 *         evaluated. This parameter can be NULL, in which case the method
 *         works on the intersection of the local constraints with the local
 *         variables (a "diagonal" block).
 *
 *  @return The number of nonzeros in the given part of the Jacobian.
 */
int 
ExpandedModel::getNzJacobianOfIntersection(ExpandedModelInterface *emcol_)
{
  ExpandedModel *emcol = static_cast< ExpandedModel* > (emcol_);

  if (emcol==NULL) emcol = this;
  
  // need to find the indices of the variables local to emcol in the
  // nlfile given by this model.

  int nvar = emcol->getNLocalVars();
  int *lvar = new int[nvar];
  
  // find the indices of the variable local to emcol in the nlfile belonging
  // to this node.
  nlfile->findIxOfLocalVarsInNlFile(emcol, lvar);

  // and ask AMPL to evaluate the Jacobian structure
  int nz = nlfile->getNoNonzerosAMPL(nvar, lvar);

  delete [] lvar;

  return nz;
}

/* ----------------------------------------------------------------------------
ExpandedModel::getJacobianOfIntersection
---------------------------------------------------------------------------- */
/** Return the Jacobian of a section of the model in sparse matrix format.
 *
 *  The matrix is defined by the intersection of the local constraints in
 *  this model with the local variables of this or another model.
 *
 *  @param[in] emcol_
 *             The model w.r.t. whose local variables the Jacobian should be
 *             evaluated. This parameter can be NULL, in which case the method
 *             works on the intersection of the local constraints with the
 *             local variables (a "diagonal" block).
 *  @param[out] colbeg
 *             Column starts of the Jacobian.
 *  @param[out] collen
 *             Column lengths of the Jacobian (not returned if NULL on call).
 *  @param[out] rownbs
 *             Row indices of nonzeros entries.
 *  @param[out] el
 *             Values of the nonzero entries.
 *
 *  @note Parameters colbeg, collen, rownbs, el are assumes to be of
 *  appropriate dimensions before the method is called, namely
 *  colbeg[n+1], collen[n], rownbs[nz], el[nz] (n=number of columns,
 *  nz=number nonzeros). The number of nonzeros in this section of the
 *  Jacobian can be obtained from a call to getNzJacobianOfIntersection().
 */
void 
ExpandedModel::getJacobianOfIntersection(ExpandedModelInterface *emcol_, int *colbeg,
					 int *collen, int *rownbs, double *el)
{
  ExpandedModel *emcol = static_cast< ExpandedModel* > (emcol_);

  if (emcol==NULL) emcol = this;
  
  // need to find the indices of the variables local to emcol in the
  // nlfile given by this model.

  int nvar = emcol->getNLocalVars();
  int *lvar = new int[nvar];
  
  // find the indices of the variable local to emcol in the nlfile belonging
  // to this node.
  nlfile->findIxOfLocalVarsInNlFile(emcol, lvar);

  // and ask AMPL to evaluate the Jacobian structure
  nlfile->fillSparseAMPL(nvar, lvar, colbeg, collen, rownbs, el);

  delete [] lvar;
}

/* -------------------------------------------------------------------------
ExpandedModel::getRowBounds
-------------------------------------------------------------------------- */
/** Return the upper and lower bounds for the constraints in this model.
 *
 *  @param[out] lower
 *              The lower bounds on the constraints.
 *  @param[out] upper
 *              The upper bounds on the constraints.
 *
 *  The method is simply a wrapper around NlFile::getRowBoundsAMPL().
 */
void
ExpandedModel::getRowBounds(double *lower, double *upper) const {
  nlfile->getRowBoundsAMPL(lower, upper);
}

/* -------------------------------------------------------------------------
ExpandedModel::getObjGradient
-------------------------------------------------------------------------- */
/** Return the gradient of the objective defined in this model.
 * 
 *  @param[out] elts
 *              The objective gradient with respect to the local variables.
 */
void
ExpandedModel::getObjGradient(double *elts)
{
  int *lvar = new int[nLocalVars];
  nlfile->findIxOfLocalVarsInNlFile(this, lvar);
  //for (int i=0;i<nLocalVars;i++) printf("EM:%d %d\n",i,lvar[i]);

  for(int i=0;i<nLocalVars;i++) elts[i] = 0.;
  nlfile->getObjAMPL(nLocalVars, lvar, elts);
  delete [] lvar;
}

/* -------------------------------------------------------------------------
ExpandedModel::getColLowBounds
-------------------------------------------------------------------------- */
/** Return the lower bounds for the local variables defined in this model.
 * 
 *  @param[out] elts
 *              The variable lower bounds.
 */
void
ExpandedModel::getColLowBounds(double *elts)
{
  int *lvar = new int[nLocalVars];
  nlfile->findIxOfLocalVarsInNlFile(this, lvar);

  for(int i=0;i<nLocalVars;i++) elts[i] = 0.;
  nlfile->getColLowBoundsAMPL(nLocalVars, lvar, elts);
  delete [] lvar;
}

/* -------------------------------------------------------------------------
ExpandedModel::getColUpBounds
-------------------------------------------------------------------------- */
/** Return the upper bounds for the local variables defined in this model.
 *
 *  @param[out] elts
 *              The variable upper bounds.
 */
void
ExpandedModel::getColUpBounds(double *elts)
{
  int *lvar = new int[nLocalVars];
  nlfile->findIxOfLocalVarsInNlFile(this, lvar);

  for(int i=0;i<nLocalVars;i++) elts[i] = 0.;
  nlfile->getColUpBoundsAMPL(nLocalVars, lvar, elts);
  delete [] lvar;
}

/** Set up the nl file for this block */
void ExpandedModel::setupNlFile(const string& name) {
  nlfile = new NlFile(name);
  model_file = name;
}

/* -------------------------------------------------------------------------
ExpandedModel::findIxOfLocalVarsInNlFile
-------------------------------------------------------------------------- */
/** Find the indices of the local variables of this model in a given nl file.
 *
 *  @param[in] nlf
 *             The NlFile that defines local constraints and all variables
 *             that are used by these constraints.
 *  @param[out] lvar
 *             Assumed to be allocated with em->nLocalVar elements: lvar[i]
 *             is the the index of the ith local variable in em in the nlfile.
 *
 *  @return The number of matches found.
 *
 * For the full doumentation see NlFile::findIxOfLocalVarsInNlFile.
 * This method belongs logically to the NlFile class, since it calculates
 * (column) sections of the columns defined in the NlFile. However since
 * we cannot use lists in the NlFile class, the actual code is here.
 */
int
ExpandedModel::findIxOfLocalVarsInNlFile(NlFile *nlf, int *lvar) {
  ExpandedModel *em = this;
  const string& nlfilename = nlf->nlfilename;
  int nvar = em->getNLocalVars();
  int count = 0; // count number of matches

  // look up if this index set has already been calculated
  if (nlf->indexList.count(em)>0){
    // we have already calculated this list
    const IndexListValue *ilv = nlf->indexList[em];
    LogEM("<< found IndexValue " + nlfilename + ":" + em->model_file + "\n");
    assert(nvar==ilv->nvar);
    for (int i=0;i<nvar;i++){
      lvar[i] = ilv->lvar[i];
      if (lvar[i]>=0) count++;
    }
    assert(count==ilv->count);
  }else{
    LogEM("<< place IndexValue " + nlfilename + ":" + em->model_file + "\n");
    // we need to calculate it
    for(int i=0;i<nvar;i++) lvar[i] = -1;

    // ------- read the names of columns defined in this NlFile ------------
    const string nlcolfile = nlfilename + ".col";
    ifstream fin(nlcolfile.c_str());
    
    if (!fin) {
      cerr << "Cannot open column name file " << nlcolfile << ".\n";
      exit(1);
    }
    
    list<string> colfilelist;
    string line;
    getline(fin, line);
    while(!fin.eof()){
      colfilelist.push_back(line);
      getline(fin, line);
    }
  
    if (GlobalVariables::prtLvl >= PRINT_VERBOSE)
      cout << "Read " <<  colfilelist.size() << " lines from "
           << nlcolfile << ".\n";
    
    // -------------- compare this listOfVarNames against this list
    int i=0;
    list<string>::const_iterator p, q;
    for (p = em->listOfVarNames.begin(); p != em->listOfVarNames.end(); ++p) {
      // (*p) is a name of a local variable. Should see if we can find this
      // in this NlFile
      int cnt=0;
      for (q = colfilelist.begin(); q != colfilelist.end(); ++q) {
	if ((*p)==(*q)){
	  lvar[i] = cnt;
	  count++;  //increase number of matches
	  break;
	}
	cnt++;
      }
      i++;
    }
    //and place a copy on the map
    int *lvarc = new int[nvar];
    for (int j = 0; j < nvar; j++) lvarc[j] = lvar[j];
    nlf->indexList[em] = new IndexListValue(nvar,lvarc,count);
  }

  return count;
}

/* -------------------------------------------------------------------------
ExpandedModel::setPrimalSolColumns
-------------------------------------------------------------------------- */
/** Upload the local variable solutions */
void ExpandedModel::setPrimalSolColumns(const double *elts) {
   assert(localVarInfoSet);
   if(!pvar) {
     pvar = new double[nLocalVars];
   }
   for (int i = 0; i < nLocalVars; ++i)
     pvar[i] = elts[i];
}

/* -------------------------------------------------------------------------
ExpandedModel::setDualSolColumns
-------------------------------------------------------------------------- */
/** Upload the local variable duals (multipliers on bounds) */
void ExpandedModel::setDualSolColumns(const double *elts) {
   assert(localVarInfoSet);
   if(!dvar) {
     dvar = new double[nLocalVars];
   }
   for (int i = 0; i < nLocalVars; ++i)
     dvar[i] = elts[i];
}

/* -------------------------------------------------------------------------
ExpandedModel::setPrimalSolRows
-------------------------------------------------------------------------- */
/** Upload the local constraints slacks */
void ExpandedModel::setPrimalSolRows(const double *elts) {
   assert(localVarInfoSet);
   if(!prow) {
     prow = new double[nLocalCons];
   }
   for (int i = 0; i < nLocalCons; ++i)
     prow[i] = elts[i];
}

/* -------------------------------------------------------------------------
ExpandedModel::setDualSolRows
-------------------------------------------------------------------------- */
/** Upload the local constraints duals (multipliers on constraints) */
void ExpandedModel::setDualSolRows(const double *elts) {
   assert(localVarInfoSet);
   if(!drow) {
     drow = new double[nLocalCons];
   }
   for (int i = 0; i < nLocalCons; ++i)
     drow[i] = elts[i];
}

/* -------------------------------------------------------------------------
ExpandedModel::outputSolution
-------------------------------------------------------------------------- */
/** Output the solution to the supplied stream with the given indent */
void ExpandedModel::outputSolution(ostream &out, int indent) {
   assert(localVarInfoSet);

   string ind(indent, ' ');
   string ind2(ind);
   string name = getName();
   bool isRoot = (name == "root");
   if (!isRoot) {
      string pname = parent->getName();
      out << "\n" << ind << name.substr(pname.size()+1) << " {" << endl;
      ind2 += "  ";
   }

   if(pvar) {
      double *pptr = pvar;
      double *dptr = dvar;
      for(list<string>::const_iterator i=getLocalVarNames().begin(); 
          i!=getLocalVarNames().end(); ++i) {
        out << ind2 << left << setw(20) << *i
            << " Value = " << setw(15) << *(pptr++);
        if(dvar) out << " Reduced cost = " << *(dptr++);
        out << endl;
      }
      out << endl;
   }

   if(prow) {
      double *pptr = prow;
      double *dptr = drow;
      for(list<string>::const_iterator i=getLocalConNames().begin(); 
          i!=getLocalConNames().end(); ++i) {
        out << ind2 << left << setw(20) << *i
            << " Slack = " << setw(15) << *(pptr++);
        if(drow) out << " Dual = " << *(dptr++);
        out << endl;
      }
   }

   for(vector<ExpandedModelInterface*>::iterator i=children.begin(); 
         i<children.end(); ++i)
      (*i)->outputSolution(out, indent + (isRoot ? 0 : 2));

   if (!isRoot)
     out << ind << "}" << endl;
}

list<SymbolTable::Entry> ExpandedModel::getObjList() const {
   return src->getObjList();
}

/** Destructor */
ExpandedModel::~ExpandedModel() {

  delete[] pvar;
  delete[] dvar;
  delete[] prow;
  delete[] drow;

  delete nlfile;

  for (unsigned int i = 0; i < children.size(); ++i)
    delete children[i];

  if (!parent)
    delete src;
}
