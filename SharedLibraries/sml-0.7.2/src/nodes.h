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

#ifndef NODES_H
#define NODES_H

#include <list>
#include <sstream>
#include <vector>

class ModelComp;
class AmplModel;

/** @class SyntaxNode
 *  This is a node in the operator tree.
 *
 *  All AMPL/SML expressions are broken down into a tree of operators. 
 *  The nodes of this tree are objects of type SyntaxNode.
 *
 *  Note that values[] is a list of untyped pointers. Normally these
 *  would point to further SyntaxNode structure (children of the current
 *  node). For an ID it is simply a pointer to the string containing
 *  the name.
 *
 *  There are a few special meanings of the values array depending on the
 *  type of node (i.e. the value of opCode).
 *  @bug This should probably be implemented by deriving subclasses, however
 *       an attempt for ID nodes resulted in problems with dynamic_casts
 *
 *  @bug A comma separated list is currently represented by a SyntaxNode with
 *       opCode==COMMA and values[0/1] pointing to the first and last member
 *       in a linked list of SyntaxNodes. This linked list is implemented using
 *       _indexNode objects (which is a SyntaxNode plus a pointer to next). 
 *       This is done to facilitate adding members to the list without 
 *       knowing the dimension beforehand. Better implemented by replacing
 *       the SyntaxNode::values array with a C++ list.
 */
class SyntaxNode {

 public:

  typedef std::vector<SyntaxNode*>::const_iterator iterator;

  iterator begin() const { return values.begin(); }
  iterator end  () const { return values.end();   }

  SyntaxNode* front() const { return values.front(); }
  SyntaxNode* back () const { return values.back();  }

  SyntaxNode* operator[](int i) const { return values[i]; }

  /** Clear the child list */
  virtual void clear() { opCode = 0; values.clear(); }

 protected:

  /** ID CODE of this node (a list can be found in ampl.tab.h) */
  int opCode;

  /** List of arguments.
   *
   *  For opCode==ID, there is usually one argument, this is a (char*) to the
   *  name of the entity. If there are two arguments the second argument is
   *  an (int*) stating the ancestor information as set by ancestor(1).ID in
   *  stochastic programming */
  std::vector<SyntaxNode*> values;

 public: 
  
  /* FIXME: not sure if these two should be part of the class. They are 
     global variables that affect the print() method */
  /** if use_global_names is set then the print() method will print out
   *  model component names as global names                               */
  static int use_global_names; 
  static AmplModel *default_model;

 public:
  // ------------------------ methods -----------------------------------

  /** Constructor */
  SyntaxNode(int opCode, SyntaxNode *val1=NULL, SyntaxNode *val2=NULL, SyntaxNode *val3=NULL);

  /** Copy constructor */
  SyntaxNode(const SyntaxNode &src);

  /** Destructor */
  virtual ~SyntaxNode();

  /** Retrieve the number of children */
  virtual int nchild() const { return values.size(); }

  /** Retrieve the opCode */
  int getOpCode() const { return opCode; }

  /** Recursive printing of expression */
  std::string print() const;

  /** Return the value of this node as a string */
  virtual std::string getValue() const { throw std::exception(); return "(fail)"; }

  /** Diagnostic printing */
  void dump(std::ostream& fout) const;

  // node is a dummy var -> remove (..)
  std::string printDummyVar() const;

  /** Return comma separated list of arguments for IDREF nodes */
  std::string getArgumentList() const;

  /** Find all the IDREF nodes at or below the current node */
  virtual void findIDREF(std::list<ModelComp*>& lmc) const;

  /** Find all the IDREF nodes at or below the current node */
  virtual void findIDREF(std::list<SyntaxNode*> *lnd);

  /** Find all nodes of opCode @a oc at or below the current node */
  virtual void findOpCode(int oc, std::list<SyntaxNode*> *lnd);

  /** Find the ModelComp (if it exists) referred to by this SyntaxNode.
   *
   *  If the expression given by this SyntaxNode is an immediate reference to
   *  a ModelComp then return that, otherwise return NULL.
   */
  ModelComp *findModelComp() const;

  SyntaxNode* merge(const SyntaxNode *src);

  /** Creates a deep copy of the nodes: SyntaxNodes pointed to are recreated
   *  as well. 
   *
   *  Non-SyntaxNode objects pointed to are not recreated, here just pointers
   *  are copied (->ref in the case of a SyntaxNodeIDREF object).
   *  The int/double entries pointed to by INT_VAL/FLOAT_VAL SyntaxNodes *are*
   *  recreated.
   */
  virtual SyntaxNode *deep_copy();

  /** Creates a copy of the node, reusing the pointer in the current node */
  virtual SyntaxNode *clone();

  virtual std::ostream& put(std::ostream& s) const;
  virtual SyntaxNode *push_back(SyntaxNode *newitem);
  virtual SyntaxNode *push_front(SyntaxNode *newitem);
};

/** @class StageNodeNode
 */
class StageNodeNode : public SyntaxNode {

 private:
  std::string value_;

 public:

  /** Current replacement string for the 'node' keyword */
  static std::string node;

  /** Current replacement string for the 'stage' keyword */
  static std::string stage;

  /** Constructor */
  StageNodeNode(int opCode_, const std::string& value = "") :
     SyntaxNode(opCode_), value_(value) {}

  std::ostream& put(std::ostream& s) const;
  SyntaxNode *clone() { return new StageNodeNode(opCode, value_); }
  SyntaxNode *deep_copy() { return clone(); }

  void setValue(const std::string& value) {
    value_ = value;
  }

};

/** @class SyntaxNodeIx
 *  A node on the operator tree representing an indexing expression.
 *
 *  This is a node on the operator tree that represents an indexing expression.
 *  A general indexing expression can be of the form:
 *
 *  {(i,j) in ARCS, k in SET: i>k} 
 * 
 *  which will be broken down into a list of 'dummy IN set' expressions plus
 *  an optional qualifier (the condition to the right of COLON).
 */
class SyntaxNodeIx : public SyntaxNode {

 private:

  SyntaxNodeIx(const int opCode_) :
    SyntaxNode(opCode_), qualifier(NULL), ncomp(0) {}

  //! Number of 'dummy IN set'-type expressions
  int ncomp;

  //! List of the dummyVarExpressions
  std::vector<SyntaxNode*> dummyVarExpr;

  //! List of the set expressions
  std::vector<SyntaxNode*> sets;

  //! The expresson to the right of the ':' (if present)
  SyntaxNode *qualifier;

  //! List of ModelComp for the indexing sets
  std::vector<ModelComp*> sets_mc;

 public:
  // --------------------------- methods ----------------------------------

  //! Constructor
  SyntaxNodeIx(SyntaxNode *on);

  //! Find if the indexing expression defines the given dummy variable
  SyntaxNode *hasDummyVar(const std::string& name);

  //! Return the list of all dummy variables defined by this index's expression
  std::list<SyntaxNode *> getListDummyVars() const;

  //! Retrieve the number of indexing expressions
  int getNComp() const { return ncomp; }

  //! Retrieve the dummy variable for the specified indexing expression
  SyntaxNode* getDummyVarExpr(int i) const { return dummyVarExpr[i]; }

  //! Retrieve the set for the specified indexing expression
  SyntaxNode* getSet(int i) const { return sets[i]; }

  //! Retrieve the ModelComp for the specified indexing expression
  ModelComp* getModelComp(int i) const { return sets_mc[i]; }

  //! Set up the ->sets, ->dummyVarExpr, ->ncomp, ->qualifier components
  void splitExpression();   
  
  //! Copy the node and all its subnodes into new data structures
  SyntaxNodeIx *deep_copy();    

  //! Diagnostic printing of member variables
  void printDiagnostic(std::ostream& fout) const;

  //! for nodes that are indexing expressions, get the set that is indexed over
  const SyntaxNode* getIndexingSet() const;
};

class ValueNodeBase {
  public:
   virtual ~ValueNodeBase() {}
   virtual double getFloatVal() const = 0;
};


/* ----------------------------------------------------------------------------
IDNode
---------------------------------------------------------------------------- */
/** @class IDNode
 *  A node on the tree representing a user identifier (ie variable name).
 */
class IDNode : public SyntaxNode {

 private:
   std::string name;
   long stochparent;
  
  public:
   IDNode(const std::string& name, long stochparent=0);
   std::string getValue() const { return name; }
   void findIDREF(std::list<ModelComp*> &lmc) { return; }
   void findIDREF(std::list<SyntaxNode*> *lnd) { return; }
   // We never search for ID:
   void findOpCode(int oc, std::list<SyntaxNode*> *lnd);
   std::ostream& put(std::ostream& s) const {
      return s << name;
   }
   SyntaxNode *deep_copy() { 
      return new IDNode(name, stochparent);
   }

   void setName(const std::string& id) {
     name = id;
   }

   std::string id() const {
     return name;
   }

   void setStochParent(long parent) {
     stochparent = parent;
   }

   long getStochParent() const {
     return stochparent;
   }

};


/* ----------------------------------------------------------------------------
SyntaxNodeIDREF 
---------------------------------------------------------------------------- */
/** @class SyntaxNodeIDREF
 *  A node on the tree representing a reference to a ModelComp.
 *
 *  IDREF is a SyntaxNode that represents a reference to a ModelComponent.
 */
class SyntaxNodeIDREF : public SyntaxNode {

 public:

  /** Pointer to the ModelComp referred to by this node */
  ModelComp *ref;

  /* stochrecourse was for the same purpose as stochparent, just that the
     recourse level was given as a SyntaxNode (i.e. expression to be
     eveluated by AMPL rather than an explicit  INT_VAL) */
  //?SyntaxNode *stochrecourse; //!< resourse level in stoch programming

   private:

  /** Levels above this one for which the reference is.
   *
   *  This field is only meaningful if the node represents a component
   *  in a stochastic program. In that case stochparent gives the recourse
   *  level of the component. This is the first argument in expressions
   *  such as xh(-1,i) which refers to xh[i] in the parent stage.
   */
  int stochparent;

 public:

  /** Default constructor */
  SyntaxNodeIDREF(ModelComp *r=NULL, SyntaxNode *val1=NULL);

  /** Constructor */
  SyntaxNodeIDREF(int opCode, ModelComp *r);
  
  /** Creates a shallow copy: points to the same components as the original */
  SyntaxNodeIDREF *clone();
  
  /** Creates a copy using all new datastructures (does not duplicate ref) */
  SyntaxNodeIDREF *deep_copy();

  std::ostream& put(std::ostream& s) const;

  /** Retrieve the level of the parent stage */
  int getStochParent() const { return stochparent; }

  /** Set the level of the parent stage */
  void setStochParent(int parent) { stochparent = parent; }

  ModelComp* getModelComp() const { return ref; }
};

/** @class ValueNode
 *  Represents a value.
 */
template<class T> class ValueNode : public SyntaxNode, virtual ValueNodeBase {
  public:
   const T value;

  public:
   ValueNode(const T new_value) :
      SyntaxNode(-99), value(new_value) {}
   double getFloatVal() const { return value; }
   std::string getValue() const;
   void findIDREF(std::list<ModelComp*> &lmc) { return; }
   void findIDREF(std::list<SyntaxNode*> *lnd) { return; }
   // We never search for INT_VAL or FLOAT_VAL:
   void findOpCode(int oc, std::list<SyntaxNode*> *lnd) { return; }
   std::ostream& put(std::ostream&s) const { return s << this->value; }
   SyntaxNode *deep_copy() { return new ValueNode<T>(value); }
   SyntaxNode *clone() { return deep_copy(); }
};
template<class T> std::string ValueNode<T>::getValue() const {
   std::ostringstream ost;
   ost << value;
   return ost.str();
}

/** @class ListNode
 *  Represents a comma separated list of SyntaxNodes.
 */
class ListNode: public SyntaxNode {

  public:
   ListNode(int opCode=',', SyntaxNode *val1 = NULL, SyntaxNode *val2 = NULL) :
     SyntaxNode(opCode, val1, val2) {}
   std::ostream& put(std::ostream& s) const;
};

/** @class OpNode
 *  Represents an operator.
 */
class OpNode : public SyntaxNode {
  public:
   SyntaxNode *left;
   SyntaxNode *right;

  public:
   OpNode(int opCode, SyntaxNode *op1, SyntaxNode *op2=NULL);
   std::ostream& put(std::ostream& s) const;
   OpNode *deep_copy();
   OpNode *clone();
   void findIDREF(std::list<ModelComp*>& lmc) {
      if(left) left->findIDREF(lmc);
      if(right) right->findIDREF(lmc);
   }
   void findIDREF(std::list<SyntaxNode*> *lnd) {
      if(left) left->findIDREF(lnd);
      if(right) right->findIDREF(lnd);
   }
   void findOpCode(int oc, std::list<SyntaxNode*> *lnd) {
      if(left) left->findOpCode(oc, lnd);
      if(right) right->findOpCode(oc, lnd);
   }
};

ListNode *addItemToListOrCreate(int oc, ListNode *list, SyntaxNode *newitem);
std::string print_SyntaxNodesymb(const SyntaxNode *node);

std::ostream& operator<<(std::ostream& s, const SyntaxNode &node);
std::ostream& operator<<(std::ostream& s, const SyntaxNode *node);

// Routines taken from ampl.h
SyntaxNode *findKeywordinTree(SyntaxNode *root, int oc);
SyntaxNode* find_var_ref_in_context(AmplModel *context, SyntaxNode *ref);
SyntaxNode* find_var_ref_in_indexing(const std::string& name);

#endif
