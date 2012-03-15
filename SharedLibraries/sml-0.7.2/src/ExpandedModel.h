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

#ifndef EXPANDED_MODEL_H
#define EXPANDED_MODEL_H

/* This is the Expanded version of the AMPL model 
(as opposed to the FLAT version). 

It serves to define a tree representation of the model that has a node for
EVERY submodel (rather than just every type of submodel as the FLAT model)

It would serve as an intermediate step between the AmplModel view of
the world (which seems the right structure for a Benders type solver -
complicating constraints/variables are associated with the parent node)
and the OOPS view of the world (where complicating
variables/constraints belong to off-diagonal subblocks)

This would still be in the Benders view.
*/

#include "ExpandedModelInterface.h"
#include "symtab.h"
#include <list>
#include <string>

class AmplModel;
class NlFile;

/** @class ExpandedModel
 *  A submodel (block) in the expanded model tree, which carries all the
 *  information to describe a block of the problem to the solver.
 *
 *  The ExpandedModel object describes a (sub)block of the model together
 *  with pointers to its children (depended sub-blocks). A tree of
 *  ExpandedModel objects describes the whole problem together with its
 *  structure.
 *
 * A tree of ExpandedModel objects is the means by which the problem
 * (after being parsed by SML) is communicated to a solver backend. It
 * thus forms the interface between SML and a solver.
 *
 * An ExpandedModel object has information on the variables,
 * constraints and objectives that describe this sub-block of the
 * problem as well as pointers to its children (sub-blocks).  It
 * further (through NlFile) provides routines to evaluate the
 * constraint and objectives defined in this block (and their
 * derivatives).
 *
 * Information about the sub-block is handled as a two-layer
 * structure. The first layer (represented by an ExpandedModel object)
 * provides 'static' information about the block, such as the
 * dimensions (number of local constraints and variables), list of
 * variable and constraint names and a list of children.
 *
 * The second layer (represented by an NlFile object) provides access to the 
 * 'non-static' information such as function and derivative values.
 *
 *  An ExpandedModel object roughly corresponds to a "diagonal" node in the
 *  OOPS matrix tree. The difference is that in the ExpandedModel tree
 *  complicating variables/constraints belong to the parent node, whereas
 *  in the OOPS matrix tree they would belong to the final diagonal child node.
 *
 *  The ExpandedModel tree is the natural representation for a (Benders)
 *  decomposition solver, where complicating variables/constraints do 
 *  indeed belong to the master problem.
 *
 * @details 
 * The ExpandedModel tree is created (constructed) from the (flat)
 * AmplModel tree by expanding the indexing sets associated with every
 * AmplModel object. From the AmplModel it gets passed a list of local
 * variable and constraint name stubs (in ExpandedModel::localVarDef -
 * consisting of the global entity name together with the first part
 * (identifying the particular node) of the indexing parameters). The
 * ExpandedModel object is also passed the corresponding *.nl file
 * (which provides only the local constraints, but not just the local
 * variables, but all the variables that are referred to in the local
 * constraints). The constructor then compares the (alphanumeric) list
 * of variable name stubs with the variables defined in the *.nl file
 * (from the corresponding *.col file) to generate the list of local
 * variables (ExpandedModel::listOfVarNames), and their indices within
 * the *.nl file (ExpandedModel::listOfLocalVars). This information is then
 * used by the routines in ExpandedModel and NlFile to provide an
 * interface into the problem to the solver.
 *
 * @todo The Amplsolver interface routines (currently accessed through
 * the field nlfile) should be made accessible from this object (by
 * providing some wrapper routines).
 */
class ExpandedModel : public ExpandedModelInterface {

 private:

  //! Number of local variables
  int nLocalVars;

  //! Number of local constraints
  int nLocalCons;

  /** Indicator if information on local variables (nLocalVars, listLocalVars,
   *  nLocalCons) has been obtained by comparing the localVarDef's with the
   *  *.col file.
   */
  bool localVarInfoSet;

  /* Store solutions */
  double *pvar, *dvar;
  double *prow, *drow;

  //! The flat model from which this expanded model was created
  AmplModel *src;

  /** Name of the *.nl file that describes the problem data */
  std::string model_file;

  /** The NlFile object associated with model_file. 
   *
   *  Provides interface to routines that evaluate objective and constraint 
   *  functions.
   */
  NlFile *nlfile;   

  //list of constraints
  // FIXME: can we assume that all constraints in *.nl belong to this node?
  //list of variables
  //FIXME: in which form? regex? numbers?
  // BOTH: regex at the start. Info on numbers can be generated by a 
  //       class method
 
  //! list of global names of local variables 
  std::list<std::string> listOfVarNames;

  //! list of local names of local variables 
  std::list<std::string> listOfLocalVarNames;
 
  //! list of global names of local constraints
  std::list<std::string> listOfConNames;

  //! list of local names of local constraints
  std::list<std::string> listOfLocalConNames;

  //! indices of local variables in the corresponding *.nl file
  std::list<int> listOfLocalVars;

  //! the locally applicable variable declarations
  std::list<std::string> localVarDef;
  
 public:

  /** Stack of instance names.
   *
   *  A stack of block instances that encode the current path through the
   *  ExpandedModel tree during the construction phase.
   */
  static std::list<std::string> pathToNodeStack;

  // -------------------------- methods ------------------------------------

  //! Constructor
  ExpandedModel(AmplModel *src_model);

  //! Destructor
  ~ExpandedModel();

  //! Recursively print the contents of this instance and of its children
  void print() const;
  
  //! Return the number of variables local to this node
  int getNLocalVars() const;

  //! Return the number of constraints local to this node
  int getNLocalCons() const;

  //! Return the names of the variables local to this node
  const std::list<std::string>& getLocalVarNames() const;

  //! Return the names of the constraints local to this node
  const std::list<std::string>& getLocalConNames() const;

  //! Return the number of nonzeros in the Jacobian of a section of the model
  int getNzJacobianOfIntersection(ExpandedModelInterface *emcol);

  //! Return the Jacobian of a section of the model in sparse matrix format
  void getJacobianOfIntersection(ExpandedModelInterface *emcol, int *colbeg,
				 int *collen, int *rownbs, double *el);

  //! Return the arrays of bounds for the constraints in this model
  void getRowBounds(double *lower, double *upper) const;

  //! Return the lower bounds for the local variables defined in this model
  void getColLowBounds(double *elts);

  //! Return the upper bounds for the local variables defined in this model
  void getColUpBounds(double *elts);

  //! Return the gradient of the objective defined in this model
  void getObjGradient(double *elts);

  //! Upload the local variable solutions
  void setPrimalSolColumns(const double *elts);

  //! Upload the local variable duals (multipliers on bounds)
  void setDualSolColumns(const double *elts);

  //! Upload the local constraints slacks
  void setPrimalSolRows(const double *elts);

  //! Upload the local constraints duals (multipliers on constraints)
  void setDualSolRows(const double *elts);

  //! Set up the nl file for this block
  void setupNlFile(const std::string& name);

  //! Find the indices of the local variables of this model in a given nl file
  int findIxOfLocalVarsInNlFile(NlFile *nlf, int *lvar);

  //! Return the unique name of this block
  std::string getName() const { return model_file; }

  //! Output the solution to the supplied stream with the given indent
  void outputSolution(std::ostream &out, int indent=0);

  //! Append the variable to the list of local variable declarations
  void appendLocalVarDef(const std::string& name) {
    localVarDef.push_back(name);
  }

  //! Set nLocalVar, listOfLocalVars, nLocalCons, listOfVarNames
  void setLocalVarInfo();

 private:
  std::list<SymbolTable::Entry> getObjList() const;

};

#endif
