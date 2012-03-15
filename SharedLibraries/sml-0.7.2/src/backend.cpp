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

#include "backend.h"
#include "AmplModel.h"
#include "GlobalVariables.h"
#include "ModelComp.h"
#include "nodes.h"
#include "sml.tab.h"
#include <cassert>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>

using namespace std;

static const int MAX_NESTED_BLOCKS = 5;

static bool prt_modwrite = false;
//produces: "Modified write (wealth), level=2, l_addIndex=2"

static void write_ampl_for_submodel(ostream& fout, AmplModel *submodel);

/* Stack of indexing expressions that are applicable to all variables that
   are printed:

   In the original AMPL-file variables are subscripted (and sets are
   indexed) *relative* to the current location in the model tree.
   
   In the AMPL-subfiles a global naming of variables is used,
   therefore indexing sets and dummy variables of the block
   definitions leading to the current node in the model tree have to
   be added to all variables.

   This is done by the addIndex stack: as the sub file writer traverses
   the tree the block-indexing expressions are added to the tree.
*/
// FIXME: this is fairly dumb at the moment: it cannot deal with 
//         multiple dimensions {i in SET1,j in SET2}
//         SET valued expressions: {i in SET1 cross SET2} 
//         or conditions:    {(i,j) in SET1:i<j}
vector<list<add_index> > l_addIndex;  /* to add to all statements */

static void
print_entry(const ModelComp *entry) {
  cout << "    " << entry->id << "\n";
  cout << "       " << *(entry->indexing) << "\n";
  cout << "       " << *(entry->attributes) << "\n";
}

static void
print_entries(const AmplModel *model, compType type) {
  list<ModelComp*>::const_iterator p;
  for (p = model->comps.begin(); p != model->comps.end(); ++p) {
    ModelComp *entry = *p;
    if (entry->type == type)
      print_entry(entry);
  }
}

static void
print_model(const AmplModel *model)
{
  ModelComp *entry;
  const AmplModel *submod;
  SyntaxNode::use_global_names=0;

  cout << "-------------------------- backend::print_model ----------------\n";
  cout << "Model: " << model->name << "\n";

  cout << "n_sets: " << model->n_sets << "\n";
  print_entries(model, TSET);

  cout << "n_cons: " << model->n_cons << "\n";
  print_entries(model, TCON);

  cout << "n_vars: " << model->n_vars << "\n";
  print_entries(model, TVAR);

  cout << "n_params: " << model->n_params << "\n";
  print_entries(model, TPARAM);

  cout << "n_obj: "<< model->n_objs << "\n";
  print_entries(model, TMAX);

  cout << "submodels: " << model->n_submodels << "\n";
  list<ModelComp*>::const_iterator p;
  for (p = model->comps.begin(); p != model->comps.end(); p++) {
    entry = *p;
    if (entry->type==TMODEL){
      print_entry(entry);
      submod = entry->other;
      print_model(submod);
    }
  }
  cout << "---END-------------------- backend::print_model ----------------\n";
}

/* ---------------------------------------------------------------------------
process_model
--------------------------------------------------------------------------- */
/* This routine will output something useful
    1) generate list of all the (sub)models defined
    2) output an AMPL file declaring each one of the submodels
   Other tasks
    3) generate the algebra tree
    4) generate the ampl script to generate a *.nl file for every
       algebra-tree node
    5) generate the corresponding structure files that partition the 
       *.nl files by columns
*/
static void fill_model_list_(AmplModel *model, list<AmplModel*>& listam);

/** Prettyprinting */
static void print_indent(ofstream& out, int k) {
  for (int j = 1; j < k; ++j)
    out << "  ";
}

int
process_model(AmplModel *model, const string& datafilename) {

  if (GlobalVariables::prtLvl >= PRINT_LOG)
    cout << "-------------- start of process_model ----------------------\n";

  /* should be called from the root model */
  assert(model->parent == NULL);

  /* 1) >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> generate list of all models */

  model->level=0; /* root is on level 0, fill_model_list also sets levels */
  list<AmplModel*> model_list;
  fill_model_list_(model, model_list);
  /* model_list[0-(n_models-1)] is now a list of all (sub)models defined
     in the ampl file */

  //printf("These are the models on the list:\n");
  if (GlobalVariables::prtLvl >= PRINT_LOG) {
    int i=0;
    for(list<AmplModel*>::iterator mli=model_list.begin();mli!=model_list.end();
        ++mli,++i){
      AmplModel *thism = *mli;
      cout << i << ": " << thism->name << " (level " << thism->level << ")\n";
    }

    cout << "----------- generate submodel files --------------\n";
  }

  /* 2) >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> output all submodel files */
  /* and */
  /* 3) >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> generate the script file */
  ofstream fscript("script.scr");
  
  /* loop over every single model defined */
  for(list<AmplModel*>::iterator mli=model_list.begin(); mli!=model_list.end();
        ++mli){
    AmplModel *this_model = *mli;
    // construct the name of the model file by concatenating the names of all
    // ancestors (and 'root' for the top-level model)
    AmplModel *tmp_model = this_model;
    string filename = tmp_model->name;
    while (tmp_model = tmp_model->parent)
      filename = tmp_model->name + "_" + filename;
    filename += ".mod";

    /* ==================== write the script file ===================== */

    fscript << "\nreset;\noption auxfiles rc;\noption presolve 0;\n";
    fscript << "model " << filename << ";\ndata ../" << datafilename << ";\n";
    
    /* the script file needs to look something like

    model submod.mod;
    data allthedata.dat;
    
    for {i in ARCS}{  
      let ARCS_SUB := {i};

      print card(COMM) > "nameofsub" &i&".crd"
      for {j in COMM}{
        let COMM_SUB := {j};
        fix all;
        unfix all;
        write("nameofsub"& i &"_"& j);
      }
    }
    */
    
    /* create a new model list that lists all models in this branch of the tree
       (i.e. all ancestors of the current model up to root */
    /* root will be entry 0, current model at entry model->level */
    AmplModel **anc_list = (AmplModel**)calloc(this_model->level+1, sizeof(AmplModel*));
    {
      AmplModel *tmp_model = this_model;
      anc_list[tmp_model->level] = tmp_model;
      while(tmp_model->parent){
        tmp_model = tmp_model->parent;
        anc_list[tmp_model->level] = tmp_model;
      }
    }

    /* go up the model_list and write a 'for{...}'/'let ...' combination
       for each model on the list */

    /* change extension of current model name to ".crd" */
    filename.erase(filename.size()-4);
    
    l_addIndex.clear();

    for(int k=1; k<=this_model->level; k++){
      AmplModel *tmp_model = anc_list[k];
      ModelComp *node = tmp_model->node; /* the node corresponding to model */
      SyntaxNode *ix = node->indexing; /* the indexing expression */
      SyntaxNode *set=NULL; /* the set that is indexed over */
      SyntaxNode *dummyVar=NULL; /* the dummy var of the indexing expression */

      // need to set SyntaxNode::default_model for component printing routines
      SyntaxNode::default_model = tmp_model; 

      /* Need to analyse the indexing expression to get the set that is
        indexed over 
        (I guess this should be done by an "indexing" object) */
      
      /* remove outside braces from indexing expression */
      list<add_index> li;
      if (ix){
        add_index ai;
        if (ix->getOpCode() == LBRACE)
          ix = ix->front();

        /* assumes that the next level is the 'IN' keyword (if present) */
        dummyVar = NULL;
        set = ix;
        if (set->getOpCode() == IN) {
          SyntaxNode::iterator ixi = set->begin();
          dummyVar = *ixi;
          set = *(++ixi);
        }
        ai.dummyVar = dummyVar;
        ai.set = set;
        li.push_back(ai);
      }

      // FIXME: how to deal with multidimensional dummy variables?
      //        need to work out the dimenension of a set?
      //        => for now simply use the original dummy variable
      //        => REQUIRE a dummy variable in block indexing
      if (set && !dummyVar){
        cerr << "ERROR: Indexing expressions on blocks need dummy variables.\n";
        cerr << "Model: " << this_model->name << " level " << k << endl;
        cerr << "Node: " << node->id << endl;
        cerr << "Indexing: " << node->indexing << "(" <<
           (void*) node->indexing << ")" << endl;
        return 1;
      }
      //print_indent(fscript, k);
      //fprintf(fscript, "for {i%d in %s }{\n", k, print_SyntaxNode(set));
      //print_indent(fscript, k + 1);
      //fprintf(fscript, "let %s_SUB := {i%d};\n", print_SyntaxNode(set), k);
      if (set){
        SyntaxNodeIDREF *set_ref = dynamic_cast<SyntaxNodeIDREF*>(set);
        if (set_ref==NULL){
          cerr << "ERROR: set reference must be SyntaxNodeIDREF\n";
          return 1;
        }
        print_indent(fscript, k);
        fscript << "for {" << dummyVar << " in " << set << " }\n";
        print_indent(fscript, k);
        fscript << "{\n";
        print_indent(fscript, k + 1);
        // the "reset data" command is needed to avoid the invalid subscript
        // error
        fscript << "reset data " << getGlobalName(set_ref->getModelComp(), set,
            SyntaxNode::default_model, NOARG) << "_SUB;\n";
        print_indent(fscript, k + 1);
        fscript << "let " << getGlobalName(set_ref->getModelComp(), set,
                                           SyntaxNode::default_model, NOARG)
                << "_SUB" << getGlobalName(set_ref->getModelComp(), set,
                                           SyntaxNode::default_model, ONLYARG)
                << " := {" << dummyVar << "};\n";
      }else{
        print_indent(fscript, k);
        fscript << "{\n";
      }
      l_addIndex.push_back(li);
    }

    /* FIXME: still need to take the "print card()" statement from the
       bit below and add it here. 
       Unclear what that really needs to do:
        - Write out all indexing sets? Just their cardinality? 
          Or the whole set? Whole set might be nice so that we can refer
          to submodels to their "proper" AMPL-name. 
        - These indexing sets might need to be subscripted as well:
          COMM might be different for each incarnation of root_MCNF
        - Should take the chance and get AMPL to write out all the information
          we need here
    */

    /* If the current model (i) has children, then write out the size
       of its childrens indexing epression to disk */
    if (this_model->n_submodels>0){
      for(list<ModelComp*>::iterator p = this_model->comps.begin();
          p!=this_model->comps.end();p++){
        ModelComp *mc = *p;
        if (mc->type != TMODEL)
          continue;

        /* found a submodel */
        AmplModel *submodel = mc->other;
        SyntaxNodeIx *ix = mc->indexing;

        // set is NULL if ix is NULL
        const SyntaxNode *set = ix->getIndexingSet();
  
        // the name of the *.crd file is the global name of the submodel
        // with all the current values of loop variables up to
        // this level attached
        
        // so buffer here should be name of current model
        string rootfilename = filename;
        filename += "_" + submodel->name;

        // for all levels from root up to here add the value of 
        // indexing variables

        filename+= "\"";
        for(int k=1;k<=this_model->level;k++){
          // get all the indexing variables and add them together (joined by &)
          SyntaxNodeIx *ixn = (anc_list[k]->node)->indexing;
          if (ixn){
            list<SyntaxNode *> dvl = ixn->getListDummyVars();
            string innerp = "";
            for(list<SyntaxNode *>::iterator q=dvl.begin();q!=dvl.end();q++){
              innerp += (*q)->print() + "&";
            }
            innerp.erase(innerp.size()-1); // delete the last '&'
            
            filename += "&\"_\"&" + innerp;
          }
        }
        filename += "&\".crd\"";
        //printf("name of file is: %s\n",buffer);

        // Need to write something like
        //  print card(indexing expression) > exact_model_name.&1.&2.txt
        // or whatever it takes to get ampl to produce filenames
        // of the form stub_1.2.crd
        // where 1, 2 are  the values of the iteration indices in the 
        // scriptfile
        if (set){
          print_indent(fscript, this_model->level + 1);
          fscript << "print card(" << set << ") > (\"" << filename << ");\n";
          filename.replace(filename.find("&\".crd\""), 7, "&\".set\"");
          print_indent(fscript, this_model->level + 1);
          fscript << "display " << set << " > (\"" << filename << ");\n";
        }else{
          fscript << "print \"0\" > (\"" << filename << ");\n";
        }

        filename = rootfilename; //delete the extension again
      }
    }

    /* write the main part of the submodel generation part in the scripts */
    print_indent(fscript, this_model->level + 1);
    fscript << "fix all;\n";
    print_indent(fscript, this_model->level + 1);
    fscript << "unfix all;\n";
    
    /* take the .mod suffix away from buffer */
    //n = strlen(buffer);buffer[n-4] = 0;
    print_indent(fscript, this_model->level + 1);
    fscript << "write (\"b" << filename << "\"";
    for(int k=1;k<=this_model->level;k++) {
      // get all the indexing variables and add them together (joined by &)
      SyntaxNodeIx *ixn = (anc_list[k]->node)->indexing;
      if (ixn){
        list<SyntaxNode *> dvl = ixn->getListDummyVars();
        string textrep = "";
        for(list<SyntaxNode *>::iterator q=dvl.begin();q!=dvl.end();q++){
          textrep += (*q)->print() + "&";
        }
        textrep.erase(textrep.size()-1); // delete the last '&'
  
        fscript << "&\"_\"&" << textrep;
      }
    }
    fscript << ");\n";
    
    /* and close all the brackets (prettyprinting) */
    for(int k=this_model->level;k>0;k--){
      print_indent(fscript, k);
      fscript << "}\n";
      //rem_from_index_stack();
      l_addIndex.pop_back();
    }
    
    /* write the submodel file */
    filename += ".mod";
    ofstream fout(filename.c_str());
    if (GlobalVariables::prtLvl >= PRINT_LOG)
      cout << "Write to model file: " << filename << endl;
    write_ampl_for_submodel(fout, *mli);
    fout.close();
    free(anc_list);
  }

  fscript.close();

  /* FIXME: while we are at it, the AMPL scripts should write out the
     cardinality of the involved indexing sets and write these
     somewhere in a form that can be retrieved by whatever routine
     generates the algebra tree. Possible ideas:
     - write a file with
        [name of block] [number of entries]
       lines that can be scanned. These lines could be generated fairly easily
       in the above script. No problem with lines generated multiple times
     - have additional scripts just for the purpose of getting at the 
       cardinality numbers
  */

  /* 3b) >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> call ampl to process script */
  {
    // call ampl to process script and analyse the output
    if (GlobalVariables::prtLvl >= PRINT_LOG)
      cout << "\nCalling AMPL to process script file... ";

    char buffer[256];
    int n_nocsobj=0; // number of model with "No constraints or objectives."
    int n_novar=0; // number of model with "No variables declared."
    int n_other=0; // other ampl error

    // "2>&1" sends stderr to stdout
    FILE *ain = popen("ampl script.scr 2>&1", "r");
    while(!feof(ain)){
      char *p;
      p = fgets(buffer, 256, ain);
      if (p){
        cout << buffer;
        if (strncmp(buffer, "No constraint",13)==0){
          n_nocsobj++;
        } else if (strncmp(buffer, "No variables",12)==0){
          n_novar++;
        } else if (strncmp(buffer, "ILOG AMPL",9)==0){
          // Do nothing, version and license string
          // eg ILOG AMPL 10.000, licensed to "university-edinburgh".
        } else if (strncmp(buffer, "AMPL Version",12)==0){
          // Do nothing version string
          // eg AMPL Version 20051214 (Linux 2.6.9-5.ELsmp)
        } else {
          n_other++;
        }
      }
    }
    pclose(ain);

    if (n_nocsobj + n_novar) {
      cout << "AMPL: Model without constraints and objectives: " <<
        n_nocsobj << "\n";
      cout << "AMPL: Model without variables                 : " <<
        n_novar << "\n";
    }
    if (n_other > 0) {
      cout << "AMPL: Unrecoverable error\n";
      return 1;
    }

    if (GlobalVariables::prtLvl >= PRINT_LOG)
      cout << "done.\n";
  }

  return 0;
}

/* ---------------------------------------------------------------------------
fill_model_list 
---------------------------------------------------------------------------- */
void
fill_model_list_(AmplModel *model, list<AmplModel*> &listam)
{
  /* fill_model_list:
     recursively creates a depth first list of all the models in the
     model tree 
     Also sets the model->level attribute of all models in the tree
     ASSUMES: that the ->level of the initial model is already set

     IN: 
       AmplModel *model:          The start node (usually root)
       list<AmplModel *> listam:  The list in which models are added
     OUT: 
       listam
       model->level
  */

  /* this works by first adding the current model to the list and then
     all its submodels. 
     It used to be that the submodels are added first, however then the
     level-calculation does not work */

  listam.push_back(model);
  if (model->parent) model->level = model->parent->level+1;

  if (model->n_submodels>0){
    for(list<ModelComp*>::iterator p=model->comps.begin(); 
          p!=model->comps.end(); ++p){
      ModelComp *comp = *p;
      if (comp->type==TMODEL)
        fill_model_list_(comp->other, listam);
    }
  }
}

/* ---------------------------------------------------------------------------
write_ampl_for_submodel
--------------------------------------------------------------------------- */
/* This is task 2: writing the ampl model for a given submodel
   
   It does the following tasks:
    - (at the moment): all set/param/var definitions above the given 
      submodel are copied verbatim
      (This could be done by using a dependency graph)
    - The declaration of the current subblock:
      - indexing 'i in ARCS' changed to 'i in ARCS_SUB'
        this indexing is added as an extra (first) parameter to all
        entities declared within this block
    - all names of var/subject to/param/set entities declared within a block
      are changed to a global name by pre-pending the name of the block to
      the entity name (to avoid local/global naming problems)
      
    - sister blocks and subblocks of the current block:
      var/set/param declarations are copied, constraints not
*/
void write_ampl_for_submodel_(ostream &fout, int thislevel, int sublevel, 
            AmplModel **list, AmplModel *submodel);
void modified_write(ostream &fout, ModelComp *comp);

void
write_ampl_for_submodel(ostream &fout, AmplModel *submodel)
{
  AmplModel *listam[MAX_NESTED_BLOCKS];
  int level;
  
  SyntaxNode::use_global_names = 1;
  l_addIndex.clear();

  if (GlobalVariables::prtLvl >= PRINT_INFO) {
    cout << "==============================================================\n";
    cout << "     ampl model for part: " << submodel->name << "\n";
    cout << "==============================================================\n";
  }

  /* need list of models at different levels from here up to root */
  {
    AmplModel *tmp;
    
    listam[0] = submodel;
    tmp = submodel;
    level = 0;
    while(tmp->parent){
      tmp = tmp->parent;
      if (++level == MAX_NESTED_BLOCKS) {
        cerr << "ERROR: Exceeded the maximum number of nested blocks "
             << "(currently set at " << MAX_NESTED_BLOCKS << ").\n";
        exit(1);
      }
      listam[level] = tmp;
    }
  }
  if (GlobalVariables::prtLvl >= PRINT_INFO) {
    cout << "-> this model is on level " <<  level << "\n";
    cout << "   Levels from top are: \n";
    for (int i = 0; i <= level; ++i)
      cout << i << ": " << listam[i]->name << "\n";
  }

  /* mark all model components that are needed by the current model */
  //ModelComp::untagAll();
  ModelComp::untagAll(AmplModel::root);
  /* and loop through all ModelComponents and mark it and its 
     dependencies as needed */
  
  for(list<ModelComp*>::iterator p = submodel->comps.begin();
      p!=submodel->comps.end();p++){
    (*p)->tagDependencies();
  }
  if (GlobalVariables::prtLvl >= PRINT_INFO) {
    cout << "processing " <<  submodel->name << "\n";
    cout << "-------> tagged now\n";
    //ModelComp::writeAllTagged();
    ModelComp::writeAllTagged(AmplModel::root);
  }

  /* now start reporting the modified model recursively */
  write_ampl_for_submodel_(fout, level, 0, listam, submodel);

  /* - indexing 'i in ARCS' changed to 'i in ARCS_SUB'
    
  This needs further clarification: The indexing expression can have several
  forms, which are changed accordingly
   - {ARCS} : just a name, no 'in' statement
      => simply replaced by ARCS_SUB 
   - {i in ARCS}: 
      => replaced by 'i in ARCS_SUB'
   - {ARCS diff i}: more complex setdefinition
      => set MCNF_IDX = ARCS diff i;
         all entities  indexed by {MCNF_IDX_SUB}
   - {i in ARCS diff i}: more complex setdefinition
      => set MCNF_IDX = ARCS diff i;
         all entities  indexed by {i in MCNF_IDX_SUB}

   it is possibly sufficient to start implementing only the first two,
   keeping in mind that the extension might be wanted later
  */
} 

/* ---------------------------------------------------------------------------
write_ampl_for_submodel_
--------------------------------------------------------------------------- */
/* This routine does part of the work in writing out an AMPL submodel file
   for the model given by 'submodel'.
   It is passed the list of ancestors of submodel (in 'list') and works
   down this list recursively. For an ancestor of the current model
   it writes out everything except the constraint definitions and submodel 
   declarations. For the current model it writes out everything.
   Global names are used throughout.
 
   This routine recursively works on the model given by list[level] and
   writes the model definition out.
   
   It ignores all submodel definitions that are not on the path to the 
   'submodel' (i.e. which are not listed on the 'list').
   
   It does special treatment for the model in question (the 'submodel')
   
   IN: 
     FILE *fout             the file where the model definition should
                            be written to
     int thislevel          level of the current node (root=0)
     int sublevel           NOT USED!
     AmplModel[] *listam     list[0] is current model, list[thislevel] is root
     AmplModel *submodel   the current model again
   OUT: none (output on data file)

   It keeps track of the stack of block indexing expressions, by putting them
   onto the l_addIndex stack
*/
void
write_ampl_for_submodel_(ostream& fout, int thislevel, int sublevel,
                         AmplModel **listam, AmplModel *submodel) {

  AmplModel *thism = listam[thislevel];
  ModelComp *comp;
  
  SyntaxNode::default_model = thism;
  // loop through all entities in this model 
  fout << "\n# start of model " << thism->global_name << "\n\n";
  for(list<ModelComp*>::iterator p = thism->comps.begin();
      p!=thism->comps.end();p++){
    comp = *p;

    if (comp->type!=TMODEL){
      // if it is not a model declaration, simply write the declaration out
      // write out *all* components of current model and everything except
      // objectives and constraints in submodels
      if (comp->type==TVAR || comp->type==TSET || comp->type==TPARAM 
          || thism==submodel){
        //SyntaxNode::default_model = comp->model;
        modified_write(fout, comp);
      }
    }

    else {
      /* this is a model declaration */
      /* check that it needs to be followed */
      if (thism!=submodel && listam[thislevel-1]==comp->other){
        list<add_index> li;
        /* ok, looks like this needs to be followed */
        /* add the indexing expression */
  
        // initialise the new entry in l_addIndex stack if not already done
        SyntaxNode *ix = comp->indexing;
        if (ix){
          add_index ai;
          
          // and place the indexing expression of this BLOCK onto the stack
          // the stack stores the dummy variable and the SET separately, 
          // so take them to bits here
          if (ix->getOpCode() == LBRACE)
            ix = ix->front(); // remove curly brackets
          if (ix->getOpCode() == IN) {
             SyntaxNode::iterator ixi = ix->begin();
            ai.dummyVar = *ixi;
            ai.set = *(++ixi);
          }else{ // no dummy variable, just a set
            ai.dummyVar = NULL;
            ai.set = ix;
          }

          /* okay we have placed the set description on the stack
             but really this should be modified:
             'i in ARCS' should read 'i in ARCS_SUB'
             ought to 'write set ARCS_SUB within ARCS'; as well
          */
          
          /* 14/03/08: what we are trying to do is to create a ModelComp
             that represents the expression
             set indset_SUB within indset;   
             and print this with modified_write. Then it should be
             automatically indexed over all subproblem indexing sets:
             set indset_SUB{ix0 in indset0_SUB} within indset[ix0];
             
             The ModelComp is of type SET, with no indexing expression
          */
          SyntaxNode *setn = ai.set;
          if (setn->getOpCode() != IDREF) {
            cerr << "At the moment can index blocks only with simple sets\n";
            cerr << "Indexing expression is " << setn->print() << "\n";
            exit(1);
          }
          SyntaxNodeIDREF *setnref = dynamic_cast<SyntaxNodeIDREF*>(setn);
          if (setnref==NULL){
            cerr << "IDREF node should be of type SyntaxNodeIDREF\n";
            exit(1);
          }
          ModelComp *setmc = setnref->getModelComp();
          {
            ModelComp newmc(setmc->id + "_SUB", TSET,
                            NULL, new SyntaxNode(WITHIN, setn));

            //newmc.model = comp->model;
            newmc.model = setmc->model;
            modified_write(fout, &newmc);
          }
          
          //fprintf(fout, "set %s_SUB within %s;\n", 
          //      print_SyntaxNode(setn), print_SyntaxNode(setn));
          /* and now modify the set declaration */

          /* newn is the new node, first copy the old one */
          SyntaxNodeIDREF *newn = setnref->clone();
          // clone the ModelComp that is referred to
          newn->ref = setmc->clone();
          // ???but associate this with the current model
          //newn->ref->model = thism;

          // set the new name by appending _SUB
          newn->getModelComp()->id = setmc->id + "_SUB";

          /* and put this on the stack */
          ai.set = newn;
          li.push_back(ai);
        }       
        l_addIndex.push_back(li);
        write_ampl_for_submodel_(fout, thislevel-1, sublevel, listam, submodel);
        SyntaxNode::default_model = thism;
        l_addIndex.pop_back();
      } /* end of (model on the current list branch) */
      else if (thislevel==0) {
        // we are in the current model and are 
        // processing a definition of a child block of the current block
        // => write out everything that is tagged
        AmplModel *childm = comp->other;
        
        //printf("\n\n-----------------------------------------------\n");
        //printf("  current model: %s\n",submodel->name);
        //printf("  these are the components needed from children:\n");
        childm->writeTaggedComponents(fout);
        SyntaxNode::default_model = thism;
        //printf("-----------------------------------------------\n\n");
        
        // For normal model components we simply call modified_write
        // on these components
        // => what does modified_write depend on?
      }
    } /* end else (end of the TMODEL branch) */
  }
}

#ifdef REDUNTANT_WAS_USED_FOR_ACL_FILES
/* --------------------------------------------------------------------------
write_columnfile_for_submodel
--------------------------------------------------------------------------- */
/* this routine just goes through the subproblem definition and prints out the
   global name of all variables defined in this subproblem
   
   Not sure at the moment what to do about indexing 
   probably just leave the stub of the variable name here and add the 
   index later (indices would have to be given as numbers, here we
   do not know the size of the set indexed over
*/
void
write_columnfile_for_submodel(ostream &fout, AmplModel *submodel)
{
  ModelComp *comp;

  for(list<ModelComp*>::iterator p = submodel->comps.begin();
      p!=submodel->comps.end(); ++p){
    comp = *p;
    if (comp->type==TVAR){
      /* print global name here: 
         either just prefix all model names up to here by a loop, or
         use the global name printing routine used elsewhere (but I
         think that also puts indexing expressions on it               */

      /* the NOARG version of getGlobalName should do the trick */
      fout << getGlobalName(comp, NULL, NULL, NOARG) << "\n";
    }
  }
}
#endif /* REDUNDANT_WAS_USED_FOR_ACL_FILES */

/* ---------------------------------------------------------------------------
modified_write
--------------------------------------------------------------------------- */
/** Writes out a component of a model.
 *
 *  Components can be modified: if this is down into a submodel, then
 *  - all declarations get new indexing expressions appended to it
 *  - all references to entities get new subscripts attached to it.
 *
 *  @param[in] fout
 *             The file to write to.
 *  @param[in] comp
 *             The component definition to write out.
 *  @pre Depends on l_addIndex: currently applicable indexing expresssions.
 *   
 *  Prints the global definition of the given ModelComponent to the given file.
 *
 *  -# get the global name of the model component
 *  -# prepend all indexing expressions on the stack to the indexing
 *      expression of this entity
 *  -# for all components that are referenced in the definition
 *    -# use their global name
 *    -# prepend the dummy variables for all indexing expressions on the stack
 *        to the argument list
 *
 *  The last part is simply done by a call to (comp->attributes)->print()
 *  (SyntaxNode::print)
 *  (with SyntaxNode::use_global_names set to true
 *   => the argument list version of ModelComp::getGlobalName is called)  
 */
void
modified_write(ostream& fout, ModelComp *comp) {

  SyntaxNode *ixsn;

  // we should check that the level of the model the component is attached to
  // tallies with the number of expressions on the indexing stack
  
  AmplModel *model = comp->model;
  int level=0;
  while (model->parent!=NULL){
    level++;
    model = model->parent;
  }
  if (prt_modwrite)
    cout << "Modified write (" << comp->id << "), level=" << level << 
      ", l_addIndex=" << l_addIndex.size() << "\n";
  assert(level <= l_addIndex.size());

  if (comp->type!=TMODEL){
    int first=1;
    /* start statement */
    fout << ModelComp::compTypes[comp->type] << " ";
    
    // find number of indexing expressions on stack
    int c_addIndex = 0;
    for (int i = 0; i < level; ++i)
      c_addIndex += l_addIndex[i].size();

    /* write name and indexing expression */
    fout << getGlobalName(comp, NULL, NULL, NOARG);
    if (c_addIndex>0 || comp->indexing)
      fout << "{";
    /* write out the additional indices */
    for (int i = 0; i < level; ++i) {
      list<add_index>& li = l_addIndex.at(i);
      for (list<add_index>::iterator p = li.begin(); p != li.end(); p++) {
        add_index ix = *p;
        if (first) {first = 0;} else {fout << ",";}
        if (ix.dummyVar)
          fout << ix.dummyVar << " in ";
        ixsn = ix.set;
        if (ixsn->getOpCode() == LBRACE)
          ixsn = ixsn->front();
        fout << ixsn;
      }
    }
    /* write out the index of this component */
    if (comp->indexing) {
      if (first) {first = 0;} else {fout << ",";}
      ixsn = comp->indexing;
      if (ixsn->getOpCode() == LBRACE)
        ixsn = ixsn->front();
      fout << ixsn;
    }
    if (c_addIndex>0 || comp->indexing)
      fout << "}";
    
    /* write out special syntax for a particular statement type */
    if (comp->type==TCON||comp->type==TMIN||comp->type==TMAX) 
      fout << ":\n   ";
    
    /* write out the rest of the expression */
    if (comp->attributes)
      fout << " " << comp->attributes;
    fout << ";\n";
    //fprintf(fout, "%s;\n",print_SyntaxNode(comp->attributes));
    //fprintf(fout, "%s;\n", print_SyntaxNodesymb(comp->attributes));
  }
}
