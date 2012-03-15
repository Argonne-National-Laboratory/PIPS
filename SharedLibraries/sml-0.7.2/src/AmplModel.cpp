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

#include "AmplModel.h"
#include "backend.h"
#include "ExpandedModel.h"
#include "GlobalVariables.h"
#include "ModelComp.h"
#include "misc.h"
#include "nodes.h"
#include "sml.tab.h"
#include <cassert>
#include <cctype>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>

using namespace std;

list<changeitem> AmplModel::changes; // initialize to empty list
AmplModel *AmplModel::root = NULL; //initialize root to NULL

extern void modified_write(ostream &fout, ModelComp *comp);

// Utility function
static bool is_int(const char *tok);

/* ---------------------------------------------------------------------------
AmplModel::AmplModel()
---------------------------------------------------------------------------- */
/** Constructor */
AmplModel::AmplModel(const string& orig_name, AmplModel *par) :
  name(orig_name),
  n_vars(0),
  n_cons(0),
  n_params(0),
  n_sets(0),
  n_objs(0),
  n_submodels(0),
  n_total(0),
  level(0),
  node(NULL),
  parent(par),
  ix(NULL)
{
  if (name != "")
   setGlobalName();
}

/** Destructor */
AmplModel::~AmplModel()
{
  for (list<ModelComp*>::iterator p = comps.begin(); p != comps.end(); p++) {
    delete *p;
  }
}

/* ---------------------------------------------------------------------------
AmplModel::setGlobalName()
---------------------------------------------------------------------------- */
void 
AmplModel::setGlobalName()
{
  if (parent==NULL){
    global_name = "";
  } else if (parent->name == "root") {
    global_name = name;
  }else{
    global_name = parent->global_name + "_" + name;
  }
}

/* ---------------------------------------------------------------------------
AmplModel::setGlobalNameRecursive()
---------------------------------------------------------------------------- */
void 
AmplModel::setGlobalNameRecursive()
{
  ModelComp *mc;
  setGlobalName();
  
  // loop through all model components
  for(list<ModelComp*>::iterator p = comps.begin();p!=comps.end();p++){
    mc = *p;
    if (mc->type==TMODEL){
      AmplModel *am = mc->other;
      am->setGlobalNameRecursive();
    }
  }
}

/** Write all tagged model components in this model and submodels to a file.
 *
 *  @param fout
 *         Name of the stream to which to write.
 */
void
AmplModel::writeTaggedComponents(ostream &fout)
{
  SyntaxNode::default_model = this;
  SyntaxNodeIx *idx = node->indexing;

  list<add_index> li;
  if (idx) {

    if (idx->getNComp() == 0) {
      cerr << "All SyntaxNodeIx should have the expression split done!\n";
      exit(1);
    }
    
    if (idx->getNComp() > 1) {
      cerr << "More than 1 indexing expression is not supported yet\n";
      cerr << "Expression: " << idx->print() << "\n";
      exit(1);
    }
    // FIXME: need to check for number of expressions and place them all
    // onto the stack
    
    //Place the indexing expression of the current model onto the addIndex stack 
    add_index ai;
    ai.dummyVar = idx->getDummyVarExpr(0);
    ai.set = idx->getSet(0);
    li.push_back(ai);
  }
  l_addIndex.push_back(li);

  // loop through all model components
  for(list<ModelComp*>::iterator p = comps.begin();p!=comps.end();p++){
    ModelComp *c = *p;
    if (c->isTagged())
      modified_write(fout, c);
    if (c->type==TMODEL) {
      AmplModel *am = c->other;
      am->writeTaggedComponents(fout);
    }
  }
  l_addIndex.pop_back();
}

/* ---------------------------------------------------------------------------
AmplModel::createExpandedModel
---------------------------------------------------------------------------- */
static string crush(const string& inst);
static void getListOfInstances(istream& file, list<string>& li);

/** Recursively create an ExpandedModel tree from the flat AmplModel.
 *
 *  This method is given a location in the Expanded tree and the corresponding
 *  node of the AmplModel tree and creates recursively the ExpandedModel 
 *  instances corresponding to it.
 *  The location in the ExpandedTree to be created is encoded in two strings:
 *  -# one gives the location in the flat AmplModel tree;
 *  -# the other gives the instance of the node.
 *
 *  @param smodelname       Name of model in flat model tree 
 *                          (this is identical to this->global_name?)
 *  @param sinstanceStub    The instance that should be created.
 *                          Concatenating of the instances up to here:
 *       INS1_INS2... where INSn is the value of the indexing variable at
 *       the n-th level.
 *
 *  The expanded tree is created by concatenating strings and reading
 *  the instance names of the nodes in the next level from the corresponding
 *  *.set file.
 *
 *  @attention It seems that the first parameter (name of flat model tree node
 *  is redundant, since this could be obtained from the current AmplModel
 *  instance.
 */
ExpandedModel*
AmplModel::createExpandedModel(const string& smodelname,
                               const string& sinstanceStub) {

  /* the route down the AmplModel tree is encoded in modelname:
         root_MOD1_MOD2 
     the instance is encoded in instancestub 
         INS1_INS2
     both together root_MOD1_MOD2_INS1_INS2.nl gives the name of the *.nl file

     for the submodel MOD3 of this node, the file
         root_MOD1_MOD2_MOD3_INS1_INS2.crd 
     gives the cardinality of this submodel. The submodel instances can be 
     obtained by either looking at the possible filename completions in the
     directory, or by looking at the set
         root_MOD1_MOD2_MOD3_INS1_INS2.set
  */ 
 
  /* - need to build this up by recursively going through the new tree
       => Can not use a routine that returns a list of ExpandedModel nodes

     - need to encode the current position on the ExpandedModel tree, this
       is given by
       + The approprite AmplModel file
       + The path up to here 
          (this is either a filename stub, encoding the AmplModel tree place 
           and the instance name at each level OR
           the AmplModel and a list of instance numbers, the names could then
           be obtained from the SETs (in the *.set files))
     - The trouble with doing anything clever/elegant is that the only 
       information on cardinality/names for the ExpandedModel tree is
       in the files on the disk (either in the filenames or in the content)

     - if we use the filename as 
  */
  if (name == "root") {
    ExpandedModel::pathToNodeStack.clear();
  }

  ExpandedModel *em = new ExpandedModel(this);
  
  // copy information on the constraints/variables?
  
  string nlfilename = smodelname;
  if (sinstanceStub != "")
    nlfilename += "_" + sinstanceStub;
  em->setupNlFile(nlfilename);

  // Open this with AMPL interface and check which variables are
  // part of this model? I guess we need to do this
  // No, just need to look at the *.col file (produced by the 'write b' 
  // command) and run a check on that.
  // What do we need to compare against? Is the stub of the variables enough?
  // Can we just run a regex match (or even simpler a match of the start of the
  // variable name? Or do we need to go into checking the correct indices.
  // => Try this on paper. root_MCNF includes Flow variables do we include
  //    all or just a subset (I guess the Flow variables are not actually local
  //    this however becomes an issue for the OOPS view
  // => write down all nodes in the tree. Vars/Constraints in each *.nl file
  //    and which ones belong to which node on both the Benders and OOPS trees

  // I assume that *all* constraints declared in the *.nl file are part of 
  // the node
  
  for(list<ModelComp*>::iterator p = comps.begin();p!=comps.end();p++){
    ModelComp *mc = *p;
    if (mc->type == TVAR){
      // need to translate this into the global name of the variable, together
      // with any applicable indexing expression
      string globname(getGlobalName(mc, NULL, NULL, NOARG));
      // and add all applicable indexing to it
      // if this is on the branch (N1,N0)->C0 then the variable indexing will
      // start with
      // ['N1','NO','CO']
      // In order to recreate this we will need to store somewhere the actual
      // path of set elements that was taken to get here
      // => Can be done by a global stack associated with ExpandedModel
      //    (similar to addIndex, but this time a C++ implementation)
      //printf("Stack of path instances is: \n");
      for(list<string>::iterator q = ExpandedModel::pathToNodeStack.begin();
          q!=ExpandedModel::pathToNodeStack.end();q++){
        // separator is '[' for first item, ',' thereafter
        if (q==ExpandedModel::pathToNodeStack.begin()){
          globname += '[';
        }else{
          globname += ',';
        }
        // the expression is (N1,N2) or N1 
        // if a simple expression (no '(') just put it to the end of 
        // globname quoted
        // otherwise, need to break this down into tokens
        // this is conceptually the same as done in breaking down the SyntaxNodeIx
        // indexing nodes for the blocks. However this time it is actual
        // instances of the dummyVariables rather than the dummy variables 
        // themselves. The expressions come from AMPL output so they have not
        // been parsed by the lexer/yacc. Need to do the parsing by hand
        if ((*q).at(0)=='('){
          char *token = strdup((*q).c_str());
          // remove first and last character
          int len = strlen(token);
          assert(token[len-1]==')');
          token[len-1]=0;
          for (char *str = token + 1; ; str = NULL) { // +1 to skip the '('
            char *tok = strtok(str, ",");
            if (tok==NULL) break;
            if (is_int(tok)){
              globname += string(tok)+",";
            }else{
              globname += "'"+string(tok)+"',";
            }
          }
          globname.resize(globname.size() - 1); // remove the last ','
          free(token);
        }else{
          if (is_int((*q).c_str())){
            globname += (*q);
          }else{
            globname += "'"+(*q)+"'";
          }
        }
      }
      //cout << globname+'\n';
      em->appendLocalVarDef(globname);
    }
    if (mc->type == TMODEL){
      int card;
      // okay this is a submodel declaration:
      // - work out the cardinality
      // - work out the set and the names of the subproblem instancestub
      // - create the ExpandedModel nodes by calling this routine
      
      string nameSetFile = smodelname + "_" + mc->id;
      if (sinstanceStub.length() > 0)
        nameSetFile += "_" + sinstanceStub;
      
      // get cardinality of this node in the ExpandedModel
      ifstream fcrd((nameSetFile+".crd").c_str());
      if (!fcrd) {
        cerr << "Cannot open file: " << nameSetFile << ".crd\n";
        exit(1);
      }

      fcrd >> card;

      // if this node is indeed repeated over an indexing set
      if (card>0){
        ifstream fset((nameSetFile+".set").c_str());
        if (!fset) {
          cerr << "Cannot open file: " << nameSetFile << ".set\n";
          exit(1);
        }

        list<string> li;
        getListOfInstances(fset, li);

        for (list<string>::iterator q =li.begin(); q != li.end(); ++q) {
          string subModelName = smodelname + "_" + mc->id;
          string subModelInst;
          if (strlen(sinstanceStub.c_str())>0) subModelInst = sinstanceStub+"_";
          subModelInst += crush(*q);
          if (GlobalVariables::prtLvl >= PRINT_LOG)
            cout << subModelName << ":" << subModelInst << endl;
          AmplModel *subampl = mc->other;
          ExpandedModel::pathToNodeStack.push_back(*q);
          ExpandedModel *subem = subampl->createExpandedModel(subModelName, subModelInst);
          ExpandedModel::pathToNodeStack.pop_back();
          subem->parent = em;
          (em->children).push_back(subem);
        }

      }else{
        // if this node is not repeated over an indexing set
        string subModelName = smodelname + "_" + mc->id;
        if (GlobalVariables::prtLvl >= PRINT_LOG)
          cout << subModelName << ":" << sinstanceStub << endl;
        AmplModel *subampl = mc->other;
        ExpandedModel *subem = subampl->createExpandedModel(subModelName, sinstanceStub);
        subem->parent = em;
        (em->children).push_back(subem);
      }
    }
  }

  em->setLocalVarInfo();
  
  return em;
}

/* ---------------------------------------------------------------------------
getListOfInstances
--------------------------------------------------------------------------- */
/* Helper routine: read the set defining statement from the given string
   extract all the set elements, crush them if mutidimensional and
   return them in a list   */
void
getListOfInstances(istream& file, list<string>& li) {

  for(string token=""; token!=":=" && !file.fail(); file >> token);
  if(file.fail()) {
    cerr << "Set definition must contain a ':='!" << endl;
    exit(1);
  }

  for(char final=' '; !file.fail() && final!=';'; ) {
     string token;
     file >> token;
     final = token.at(token.length()-1);
     if(final==';') token = token.substr(0,token.length()-1);
     li.push_back(token);
     // FIXME: multidimensional?
  }
}

/* ---------------------------------------------------------------------------
crush
---------------------------------------------------------------------------- */
/** Crushes a set element "(N0,N1)" into "N0N1". */
string crush(const string& inst) {

  if (inst[0] != '(')
    return inst;

  size_t last = inst.length() - 1, pos;
  assert(inst[last] == ')');

  // remove opening and closing brackets
  string crush = inst.substr(1, last - 1);

  // remove commas
  while ((pos = crush.find(',')) != string::npos)
    crush = crush.erase(pos, 1);

  return crush;
}

/* ---------------------------------------------------------------------------
AmplModel::print
---------------------------------------------------------------------------- */
void 
AmplModel::print() const {
  dump(cout);
}

/* ---------------------------------------------------------------------------
AmplModel::dump(const char *filename)
---------------------------------------------------------------------------- */
void 
AmplModel::dump(const char *filename) const {
  ofstream fout(filename);
  dump(fout);
}

/* ---------------------------------------------------------------------------
AmplModel::dump(ostream &fout)
---------------------------------------------------------------------------- */
void 
AmplModel::dump(ostream& fout) const {
  fout << "AM: ----------------------------------------------------------\n";
  fout << "AM: This is AmplModel (" << (void *) this << "): " << name << "\n";
  fout << "AM: global name: " << global_name << "\n";
  fout << "AM: level: " << level << "\n";
  fout << "AM: parent: " << (parent?parent->name:"NULL") << "\n";
  if (ix) {
    fout << "AM: indexing: " << ix << "\n";
    ix->dump(fout);
  }
  fout << "AM: Nb submodels  : " <<  n_submodels << "\n";
  fout << "AM: Nb sets       : " <<  n_sets << "\n";
  fout << "AM: Nb parameters : " <<  n_params << "\n";
  fout << "AM: Nb objectives : " <<  n_objs << "\n";
  fout << "AM: Nb variables  : " <<  n_vars << "\n";
  fout << "AM: Nb constraints: " <<  n_cons << "\n";
  fout << "AM: Nb objectives : " <<  n_submodels << "\n";
  fout << "AM: List components:\n";
  list<ModelComp*>::const_iterator p;
  for (p = comps.begin(); p != comps.end(); ++p)
    (*p)->dump(fout);

  if (n_submodels>0)
    fout << "AM: now list the submodels:\n";
  for (p = comps.begin(); p != comps.end(); ++p) {
    ModelComp *mc = *p;
    if (mc->type == TMODEL){
      AmplModel *am = mc->other;
      am->dump(fout);
    }
  }
}

/* ---------------------------------------------------------------------------
AmplModel::check()
---------------------------------------------------------------------------- */
/** Check instance for consistency */
void
AmplModel::check() const {

  if (name == "") {
    cerr << "AmplModel::check: AmplModel has no name given\n";
    exit(1);
  }
  if (global_name==""){
    cerr << "AmplModel has no global name set\n";
    exit(1);
  }
  if (node == NULL && name != "root") {
    cerr << "AmplModel " << name << " is not root but has no ModelComp "
       "node associated\n";
    exit(1);
  }
  
  if (n_vars<0){
    cerr << "AmplModel " << name << ": n_vars = " << n_vars << "\n";
    exit(1);
  }
  if (n_cons<0){
    cerr << "AmplModel " << name << ": n_cons = " << n_cons << "\n";
    exit(1);
  }
  if (n_params<0){
    cerr << "AmplModel " << name << ": n_params = " << n_params << "\n";
    exit(1);
  }
  if (n_sets<0){
    cerr << "AmplModel " << name << ": n_sets = "<< n_sets << "\n";
    exit(1);
  }
  if (n_objs<0){
    cerr << "AmplModel " << name << ": n_objs = " << n_objs << "\n";
    exit(1);
  }
  if (n_submodels<0){
    cerr << "AmplModel " << name << ": n_submodels = " << n_submodels << "\n";
    exit(1);
  }
  if (n_vars+n_cons+n_params+n_sets+n_submodels+n_objs!=n_total){
    cerr << "AmplModel " << name << ": n_total does not equal sum of comps\n";
    exit(1);
  }
  
  if (parent == NULL && name != "root") {
    cerr << "AmplModel " << name << " is not root but has no parent\n";
    exit(1);
  }
}

/* ---------------------------------------------------------------------------
AmplModel::addDummyObjective()
---------------------------------------------------------------------------- */
/** Add a dummy objective that uses (sums up) all variables in the model.
 *
 *  AMPL will remove variables from the model that are not used in any
 *  constraint/objective within the model. In order to prevent this,
 *  we need to add a dummy objective that uses every defined variable
 *  (by simply summing them up).
 */
void
AmplModel::addDummyObjective()
{
  vector<SyntaxNode*> list_on_sum;

  // get list of all variable declarations in the model:
  /* build up list_on_sum which is a list of all variable definitions in this
     model:
      - for variables without indexing expression a SyntaxNodeIDREF pointing
        to this ModelComp is added
      - for variables with indexing expression a SyntaxNode tree that encodes
          sum {dumN in SET} var[dumN]
        is added
  */
  list<ModelComp*>::iterator p;
  for (p = comps.begin(); p != comps.end(); ++p) {
    ModelComp *comp = *p;
    if (comp->type==TVAR){
      if (comp->indexing){
        SyntaxNodeIx *idx = comp->indexing;
        // Need to create "sum{dummy in set} var[dummy]":
        // 1) create "dummy in set"
        vector<SyntaxNode*> commaseplist;
        SyntaxNode *commasepon;
        for (int i = 0; i < idx->getNComp(); ++i) {
          SyntaxNode *dummyVar = idx->getDummyVarExpr(i);
          if (!dummyVar)
            // need to make up a dummy variable
            dummyVar = new IDNode("dum" + to_string(i));
          SyntaxNode *newon = new OpNode(IN, dummyVar, idx->getSet(i));
          commaseplist.push_back(newon);
        } // end for
        // make the commaseplist
        if (idx->getNComp() == 1)
          commasepon = commaseplist[0];
        else {
          commasepon = new ListNode(COMMA);
          for (int i = 0; i < idx->getNComp(); ++i)
            commasepon->push_back(commaseplist[i]);
        }
        SyntaxNodeIDREF *onref = new SyntaxNodeIDREF(comp);
        for (int i = 0; i < idx->getNComp(); ++i) {
          // this is the dummy variable of the i-th indexing expression
          SyntaxNode *ondum = commaseplist[i]->front();
          if (ondum->getOpCode() == LBRACKET)
            ondum = ondum->front();
          onref->push_back(ondum);
        }
        // make the sum part
        commasepon = new SyntaxNode(LBRACE, commasepon);
        list_on_sum.push_back(new SyntaxNode(SUM, commasepon, onref));
      }
      else // no indexing expression, simply add this node
        list_on_sum.push_back(new SyntaxNodeIDREF(comp));
    }
  }

  // we have now built a list of SyntaxNodes representing the components:
  // build the attribute SyntaxNode as a sum of them all
  if (list_on_sum.size()>0){
    SyntaxNode *attr = list_on_sum[0];
    for (unsigned int i = 1; i < list_on_sum.size(); ++i)
      attr = new OpNode('+', attr, list_on_sum[i]);

    ModelComp *newobj = new ModelComp("dummy", TMIN, NULL, attr);
    this->addComp(newobj);
  }

  // and recursively do this for all AmplModels below this one
  for (p = comps.begin(); p != comps.end(); ++p) {
    ModelComp *comp = *p;
    if (comp->type==TMODEL){
      // might be that we need to add the indexing expression on the stack 
      comp->other->addDummyObjective();
    }
  }
}

/* --------------------------------------------------------------------------
AmplModel::removeComp()
---------------------------------------------------------------------------- */
void
AmplModel::removeComp(const ModelComp *comp) {
  // FIXME: check of comp is indeed on the list
  
  bool found = false;
  for(list<ModelComp*>::iterator p = comps.begin();p!=comps.end();p++){
    if ((*p)->id == comp->id) {
      comps.erase(p); // this invalidates the iterator => break from loop
      found = true;
      break;
    }
  }
  if (!found){
    cerr << "ERROR: Attempting to remove component " << comp->id
         << " from model " << name << ":\nComponent is not in model.\n";
    exit(1);
  }
  n_total--;
  switch(comp->type){
    case TVAR: 
      n_vars--;
      break;
    case TCON: 
      n_cons--;
      break;
    case TPARAM: 
      n_params--;
      break;
    case TSET: 
      n_sets--;
      break;
    case TMIN: 
    case TMAX: 
      n_objs--;
      break;
    case TMODEL: 
      n_submodels--;
      break;
    default:
      cerr << "removeComp: type " << comp->type << "unrecognised\n";
      exit(1);
  }
}
/* --------------------------------------------------------------------------
AmplModel::removeComp()
---------------------------------------------------------------------------- */
void
AmplModel::addComp(ModelComp *comp)
{
  AmplModel *subm;
  //ModelComp *lastinmodel = model->last;
  switch(comp->type)
  {
    case TVAR:
      n_vars++;
      break;
    case TCON:
      n_cons++;
      break;
    case TSET:
      n_sets++;
      break;
    case TPARAM:
      n_params++;
      break;
    case TMIN:
    case TMAX:
      n_objs++;
      break;
    case TMODEL:
      n_submodels++;
      subm = comp->other;
      subm->parent = this;
      break;
    default:
      cerr << "addComp: type " << comp->type << "unrecognised\n";
      exit(1);
    }
  n_total++;
  comp->model = this;
  comps.push_back(comp);

  comp->setUpDependencies();
}

/* --------------------------------------------------------------------------
AmplModel::applyChanges()
---------------------------------------------------------------------------- */
void
AmplModel::applyChanges()
{
  for (list<changeitem>::iterator p = changes.begin();
       p != changes.end(); ++p) {
    changeitem& ch = *p;
    if (ch.action == CHANGE_REM)
      ch.model->removeComp(ch.comp);
    if (ch.action == CHANGE_ADD)
      ch.model->addComp(ch.comp);
  }
  changes.clear();
}

/* --------------------------------------------------------------------------
AmplModel::reassignDependencies()
---------------------------------------------------------------------------- */
/** Recursively recalculate dependency list and re-resolve IDREF nodes.
 *
 *  In the process of building the AmplModel tree from the StochModelTree
 *  some of the IDREF dependency nodes still point to the StochModelComp
 *  nodes from the StochModel tree (or the intermediate tree)
 *
 *  This routine goes through all components and makes sure that IDREF
 *  nodes are resolved with respect to the correct ModelComp and that
 *  the dependency lists are in order.
 *  Recursively follows down submodels.
 */
void
AmplModel::reassignDependencies()
{
  for(list<ModelComp*>::iterator p=comps.begin();p!=comps.end();p++){
    if ((*p)->type==TMODEL){
      AmplModel *submodel = (*p)->other;
      submodel->reassignDependencies();
    }else{
      (*p)->reassignDependencies();
    }
  }
}

/* --------------------------------------------------------------------------
AmplModel::findComponent(string id)
---------------------------------------------------------------------------- */
/** Find a component with name id in correct scoping order.
 *
 *  It will first search this model's SymbolTable, and if it cannot find
 *  the component it will recurse to its parent node and so on up to the root.
 */
const SymbolTable::Entry *AmplModel::findComponent(const string& id) const {
   const SymbolTable::Entry *ent = symbol_table.findSymbol(id);
   if (!ent && parent)
     ent = parent->findComponent(id);
   return ent;
}

/** Returns a list of all objective functions in context */
list<SymbolTable::Entry> AmplModel::getObjList() const {
   list<SymbolTable::Entry> result = 
      symbol_table.getListByType(SymbolTable::ST_OBJ);
   if(parent) {
      list<SymbolTable::Entry> pres = parent->getObjList();
      for(list<SymbolTable::Entry>::const_iterator i=pres.begin(); i!=pres.end(); ++i)
         result.push_back(*i);
   }
   return result;
}

/* ---------------------------------------------------------------------------
bool is_int(char *tok)
---------------------------------------------------------------------------- */
/** Checks if a cstring represents a natural number (i.e [0-9]*) or not */
bool 
is_int(const char *tok){
  const char *cp = tok;
  if (!cp) return false;
  while (*cp)
    if (!isdigit(*cp++)) return false;
  return true;
}

SyntaxNodeIDREF* AmplModel::find_var_ref_in_context(IDNode *ref) {
   for(list<ModelComp*>::iterator p=comps.begin(); p!=comps.end(); ++p){
      ModelComp *thismc = *p;
      if (ref->id() == thismc->id) {
         /* this is a match */
         if (GlobalVariables::logParseModel){
            cout << "Found Match: " << ref->id() << " refers to ";
            cout << ModelComp::nameTypes[thismc->type] << "\n";
            cout << "    " << thismc->id << "\n";
            if (thismc->indexing)
              cout << "       " << *(thismc->indexing) << "\n";
            if (thismc->attributes)
              cout << "       " << *(thismc->attributes) << "\n";
         }

         SyntaxNodeIDREF *ret;
         if(thismc->type==TMODEL) {
           ret = new SyntaxNodeIDREF(IDREFM, thismc);
         } else {
           ret = new SyntaxNodeIDREF(IDREF, thismc);
         }

         return ret;
      }
   }

   /* need also to look through parent model */
   if (parent) return parent->find_var_ref_in_context(ref);

   /* need also to look through list of local variables */

   cerr << "ERROR: Could not find ref '" << ref->id() << "' in context\n";
   exit(1);
}
