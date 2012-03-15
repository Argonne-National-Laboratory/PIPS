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

#include "CompDescrParam.h"
#include "data.tab.h"
#include "GlobalVariables.h"
#include "ModelComp.h"
#include "misc.h"
#include "nodes.h"
#include <cassert>
#include <cstdlib>
#include <cstring>
#include <iostream>

using namespace std;

/* ---------------------------------------------------------------------------
CompDescrParam::CompDescrParam(ModelComp *mc, SyntaxNode *desc)
---------------------------------------------------------------------------- */
/** Parses the parameter description given in a data file.
 *
 *  This routine constructs a CompDescrParam (the actual parameter values)
 *  from a tree of SyntaxNodes originating from the data file.
 *
 *  @param mc
 *         Reference to the Component in the model file (so we can get indexing
 *         sets and dimension).
 *  @param desc
 *         The SyntaxNode tree giving the parameter value description as
 *         expressed in the data file.
 */
CompDescrParam::CompDescrParam(ModelComp *mc, SyntaxNode *desc):
  nix(0),
  indices(NULL),
  n(1),
  nread(0),
  values(NULL)
{

  /* 
     1) Match the parameter with a description in the model file and work
        out the dimension of the parameter (scalar?, vector?, array?) and
        how many parameter values to expect
     2) Parse the SyntaxNode-tree and fill in the corresponding entries of
        this object
  */


  /* ------------- 1) Work out dimension of parameter --------------------- */

  SyntaxNode *paramspec;

  /* First work out the dimension and cardinality of this parameter: 
     for this analyse the indexing expression given in the corresponding
     model component */
  SyntaxNodeIx *ix = mc->indexing;

 /* this should just be a comma separated list of sets */
  /* SyntaxNodeIx has components: 
     - int ncomp
     - SyntaxNode **sets
     - SyntaxNode *qualifier (if any ":" expression is given) */

  if (ix){
    /* if an indexing set is given, then the parameter is a vector (matrix) */
    /* work out the dimension and cardinality of the parameter */
    nsets = ix->getNComp(); // the number of indexing sets given
    indices = new Set*[nsets];
    for (int i = 0; i < nsets; ++i) {
      indices[i] = dynamic_cast<Set*>(ix->getModelComp(i)->getValue());
      if (indices[i]==NULL){
        cerr << "ERROR: Value of parameter " << mc->id
             << " given before indexing set " << ix->getModelComp(i)->id
             << " is defined\n";
        /* Of course it is possible to define the membership of the sets at 
           the same time as the parameters 
           Although a particular syntax needs to be used. Should probably have
           a separate routine for that case
        */
        exit(1);
      }
      n *= indices[i]->size(); // number of parameters is cartesian product of sets 
      nix += indices[i]->dim();
    }
  }else{ /* no indexing expression => scalar parameter */
    nsets = 0;
  }

  /* --------- 2) parse the actual parameter description ----------------- */

  /* Parameter declarations can look as follows
   - data_paramdef: PARAM param_name paramdefault_opt DEFINED paramspec_list
                    PARAM param_name paramdefault_opt value_table_list  

   - data_paramdef_alt: 
        PARAM paramdefault_opt COLON param_name_list DEFINED value_list
        PARAM paramdefault_opt COLON set_name COLON param_name_list 
            DEFINED value_list

     + paramspec: param_template_opt value value_list
                  param_template_opt value_table_list

     + value_table: (TR) COLON col_label_list DEFINED value_list
  
   NOTE: The second version is not supported yet, and will probably be
         implemented by a different function(?)

                      TOKPARAMSPECLIST
                    /                  \
   (TOKPARAMTEMPLATE)                      TOKVALUETABLELIST
            |                                      |
   ' ' (placeholder for value_list)          TOKVALUETABLE
            |                                      |
          values                 ' '(col_label_list), ' '(value_list)
                                           |                   |
                                         col_labels           values

   In any case, the tree should start with a TOKPARAMSPECLIST token

     OK, we are probably ready to analyse the actual parameter description 


   - paramspec_list: ' ' separated list of 'paramspec'

   - paramspec:  param_template_opt value value_list
     param_template_opt value_table_list
   the last two are not decided yet how they are represented on the tree
   I guess we should invent some clear tokens
  */

  values = new double[n];
  memset(values, 0, n * sizeof(double));

  // the paramdefinition is a list of PARAMSPEC's
  assert(desc->getOpCode() == TOKPARAMSPECLIST);

  // the value list at the heart of the parameter specification is a list 
  // of the form
  // object ... object entry
  // with nobj objects. The objects specify for which component of the
  // indexing set(s) the value is given. 
  // Some components of the indexing set(s) can be given by templates
  // in which case nobj is reduced
  int nobj = nix;
  for(SyntaxNode::iterator i=desc->begin(); i!=desc->end(); ++i){
    SyntaxNode *templ = NULL;
    paramspec = *i;

    // either of which can have a param_template
    if (paramspec->getOpCode() == TOKPARAMTEMPLATE) {
      SyntaxNode::iterator psi = paramspec->begin();
      templ = *psi;
      paramspec = *(++psi);
    }

    // now this can be either a "value list" or a "value table list"
    if (paramspec->getOpCode() == TOKVALUETABLELIST)
      processValueTableList(paramspec, ix);
    else {
      // this is a value list, i.e. a ' '-separated list of values
      // need to know the dimension of the parameter
      assert(paramspec->getOpCode() == ' ');
      int nval = paramspec->nchild()/(nobj+1);
      if (paramspec->nchild()%(nobj+1)!=0){
        cerr << "ERROR: Paramdef requires " << nobj << " position specifiers.\n"
             << "Length of value list is " << paramspec->nchild()
             << " which is not divisible by " << nobj+1 << ".\n";
        exit(1);
      }
      // loop through all parameter values that are given
      for(int j=0;j<nval;j++){
        // need to think how these are stored
        int pos_in_paramspec = j*(nobj+1);
        
        // get the "nobj" identifiers
        int pos_in_array = 0;
        for(int k=0;k<nobj;k++){
          pos_in_array*= indices[k]->size();
          SyntaxNode *onid = (*paramspec)[pos_in_paramspec + k];
          char *obj = (char *) onid->getValue().c_str();
          // FIXME: if we want to allow multidimensional indexing sets we
          //        need to gather indices[k]->dim entries together and 
          //        convert them into a string[]. This is then passed to
          //        findPos
          pos_in_array += indices[k]->findPos(SetElement(1,&obj));
        }          
        SyntaxNode *on = (*paramspec)[pos_in_paramspec + nobj];
        assert(on->getOpCode() == ID || on->getOpCode() == -99);
        values[pos_in_array] = ((ValueNodeBase*)on)->getFloatVal();
        nread++;
      }// end loop over j
    }

  }
  if (GlobalVariables::prtLvl >= PRINT_LOG) {
    cout << "Paramdef finished processing: " << mc->id << endl;
    cout << "  Read " << nread << " values, need " << n << endl;
  }
}

/* ---------------------------------------------------------------------------
CompDescrParam::toString
---------------------------------------------------------------------------- */
string
CompDescrParam::toString() const {
  string str = "";
  
  if (nix==0) str += "scalar";
  for(int i=0;i<nix;i++){
    if (i>0) str+="x";
    Set *st = indices[i];
    str += to_string(st->dim());
  }

  str+=": ";
  for(int i=0;i<n;i++){
    if (i>0) str += " ";
    str += to_string(values[i]);
  }

  return str;
}

static int
labelToSetElement(SyntaxNode *cl, const SyntaxNodeIx *ix,
                  Set **indices, int ixcolset) {
  IDNode *cli;
  assert(cl->getOpCode() == ID || cl->getOpCode() == LBRACKET);
  if (indices[ixcolset]->dim() > 1) {
    if (cl->getOpCode() != LBRACKET) {
      cerr << "ERROR: col_label for multidimensional set '"
           << ix->getModelComp(ixcolset)->id
           << "' must be a bracketed '(..,..)' expression.\n";
      exit(1);
    }
    if (cl->nchild() != indices[ixcolset]->dim()) {
      cerr << "ERROR: Number of entries in bracketed expression used as "
           << "col_label (" << cl->nchild() << ")\n";
      cerr << "does not match dimension (" << indices[ixcolset]->dim()
           << ") of indexing set '" << ix->getModelComp(ixcolset)->id << "'.\n";
      exit(1);
    }
    cl = cl->front(); // remove the brackets
    assert(cl->getOpCode() == COMMA);
    cli = (IDNode *) cl->front();
  }
  else {
    // this is a set of dim 1
    cli = (IDNode *) cl;
  }

  return indices[ixcolset]->findPos(SetElement(1, &cli));
}

/* ---------------------------------------------------------------------------
CompDescrParam::processValueTableList
---------------------------------------------------------------------------- */
/** Processes a part of a SyntaxNode-tree representing a value_table_list.
 *
 *  @param node
 *         The SyntaxNode of type TOKVALUETABLELIST that describes the values.
 *  @param ix
 *         The indexing expression of the ModelComp that represents the
 *         parameter (to get ix->sets_mc).
 */
void
CompDescrParam::processValueTableList(const SyntaxNode *node,
                                      const SyntaxNodeIx *ix) {
  // this is a value_table_list
  
  /* A SyntaxNode of type TOKVALUETABLELIST can have the following structure:
     - TOKVALUETABLELIST has ->nval children of type TOKVALUETABLE
     - TOKVALUETABLE has 2 children:
       + col_label_list
       + value list
     - col_label_list has opCode=' ' and is a list of objects 
     - value list has opCode=' ' and a list of objects
     - object has opCode=ID and a pointer to the object in question

                 TOKVALUETABLELIST
                          |
                   TOKVALUETABLE
                          |
      ' '(col_label_list), ' '(value_list)
              |                   |
            col_labels           values 
   
     In the simplest case TOKVALUETABLELIST gives values for a 2 dim
     parameter array. The col_label_list gives the indices for the
     first indexing set (i.e. should have a matching no of arguments: n1)
     The value_list consists of rows giving the index for the second
     indexing set and the values (i.e. should have (n1+1)*n2 arguments
  */
  
  /* In the simplest case the col_labels refer to the indexing set in pos 0
     and the row_labels to the indexing set in pos 1. For 'TR' these two
     are interchanged. When using templates these might be different
     values:  */
  
  int ixcolset = 1;
  int ixrowset = 0;
  
  /* Generalisations of this are:
     - A value table could have row/col indices reversed 
     ("TR" in SML file, how is this marked?:
     + TOKVALUETABLETR?
     + a third (dummy) argument to TOKVALUETABLE?
     - 
  */

  // loop through all value_table's
  for(SyntaxNode::iterator i=node->begin(); i!=node->end(); ++i){
    SyntaxNode *on_valtab = *i;
    assert(on_valtab->getOpCode() == TOKVALUETABLE);
    
    // get the dimensions of the value_table_list
    assert(on_valtab->nchild()==2);
    SyntaxNode::iterator ovti = on_valtab->begin();
    SyntaxNode *on_collabel = *ovti;
    ListNode *on_values = (ListNode*)*(++ovti);
    
    int ncol = on_collabel->nchild();
    int nval = on_values->nchild();
    if (nval%(ncol+1)!=0){
      cerr << "ERROR: Number of values (" << nval << ") not divisible by"
           << " col_labels+1 (" << ncol + 1 << ").\n";
      exit(1);
    }
    // dimensions seem to tally, now get a list of row and column
    // labels and look up their position in the appropriate indexing sets
    int nrow = nval/(ncol+1);
    int *rowpos = new int[nrow];
    int *colpos = new int[ncol];
    // process the col_label list first
    
    int jj=0;
    for(SyntaxNode::iterator j=on_collabel->begin();j!=on_collabel->end();++j,++jj){
      SyntaxNode *cl = *j;
      colpos[jj] = labelToSetElement(cl, ix, indices, ixcolset);
    }
    
    // also do the same loop for row_labels
    for(int j=0;j<nrow;j++){
      int pos = j*(ncol+1); // get position of next row label
      SyntaxNode *rl = (*on_values)[pos];
      rowpos[j] = labelToSetElement(rl, ix, indices, ixrowset);
    }
    // OK, we now have the lists colpos/rowpos that tell us where the
    // elements go in values/symvalues
    for (int j=0;j<ncol;j++){
      for (int k=0;k<nrow;k++){
        int poslist = (j+1)+k*(ncol+1);
        int posparam = colpos[j]+rowpos[k]*indices[ixcolset]->size();
        SyntaxNode *entry = (*on_values)[poslist];
        values[posparam] = ((ValueNodeBase *) entry)->getFloatVal();
        nread++;
      }
    }
  }
}
