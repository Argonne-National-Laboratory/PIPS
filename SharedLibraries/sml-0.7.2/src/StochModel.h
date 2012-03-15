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

#ifndef STOCHMODEL_H
#define STOCHMODEL_H

#include "AmplModel.h"
#include <string>
#include <vector>

class ModelComp;
class SyntaxNode;
class SyntaxNodeIDREF;

/* ------------------------------------------------------------------------ */
/** @class StochModel
 *  This class describes a stochastic model (block).
 *
 *  It will gather the information present in the SML model file for this
 *  block much in the same way that AmplModel does for ordinary blocks.
 *  The difference is that sets and parameters that define the 
 *  Scenario tree are associated with it. These are
 *  - stageset:    (ordered) set of stages
 *  - nodeset:     set of nodes
 *  - anc:         parameter array giving ancestor node for every node
 *  - prob:        parameter array giving conditional probability for each node
 *
 *  In principle the stochastic model block can also be repeated "horizontally"
 *  in the same manner as all other blocks by specifying an indexing 
 *  expression.
 *  The stochastic model block will be expanded at processing time into a 
 *  nested set of AmplModels. 
 */
class StochModel: public AmplModel{

 private:

  // FIXME: stage should be replaced by a Set (decribing the elements)

  //! The set of STAGES
  SyntaxNode *stageset;

  //! The dummy variable for the STAGES set
  IDNode *stagedummy;

  //! Explicit set of STAGES
  std::vector <std::string> stagenames;

  //! Whether stage names are symbolic or numeric
  bool is_symbolic_stages;

  //! The set of NODES
  SyntaxNode *nodeset;

  //! The dummy variable for the NODES set
  IDNode *nodedummy;

  //! The parameter array of ancestors
  SyntaxNode *anc;

  //! The parameter array of probabilities
  SyntaxNode *prob;

 public:
  // -------------------------- methods ----------------------------------
  //! Constructor 
  StochModel(SyntaxNode *onStages, SyntaxNode *onNodes, SyntaxNode *onAncs, 
	     SyntaxNode *onProb, AmplModel *parent);

  //! Expand the StochModel to a nested set of flat models
  AmplModel *expandToFlatModel();

  //! Expand the STAGES set into the actual elements to be stored in stagenames
  void expandStages();

  //! Expand the STAGES set of all StochModelComps in this model
  void expandStagesOfComp(); 

  //! Expand on AmplModel::addComp to setup stochmodel of component too
  void addComp(ModelComp *comp);

  //! Retrieve the SyntaxNode corresponding to the probability term
  SyntaxNode* getProbs() const { return prob; }

  SyntaxNodeIDREF* find_var_ref_in_context(IDNode *ref);

 private:

  //! Recursive helper function for expandToFlatModel
  void _transcribeComponents(AmplModel *current, int level);
};

#endif
