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

#ifndef MODELCOMP_H
#define MODELCOMP_H

#include "CompDescr.h"
#include <list>
#include <string>
#include <vector>

class AmplModel;
class StochModel;
class SyntaxNode;
class SyntaxNodeIx;
class SyntaxNodeIDREF;

enum { NOARG = 0,
       WITHARG = 1,
       ONLYARG = 2};

/** Possible types of a model component */
typedef enum {TVAR, TCON, TPARAM, TSET, TMIN, TMAX, TMODEL, TNOTYPE} compType;

/** @class ModelComp
 *  Object to represent a component of an AMPL/SML model/block.
 *
 *  It usually represents one line of AMPL/SML which is a definition of a
 *  variable/parameter/set/constraint/objective/block.
 *
 *  A model component is broken down into the following parts:
 *  - type: the type of the component
 *  - id: the name of the component
 *  - indexing: the indexing expression stored as a tree
 *  - attributes: a list of attributes (this includes the actual constraint
 *                definition for "subject to" components)
 */
class ModelComp{
 public:
  static const std::string nameTypes[];
  static const std::string compTypes[];

  /** Type of the component */
  compType type;

  /** Name of the component */
  std::string id;

  /** A tree of specifications (which includes :=, within, default, >=) */
  SyntaxNode *attributes;   
			 
  /**  Indexing expression */
  SyntaxNodeIx *indexing;

  /** The model this component belongs to */
  AmplModel *model;

  /** A pointer to an AmplModel structure for components of type MODEL
   *  @attention Better implemented as a subclass of ModelComp. */
  AmplModel *other;

 protected:

  /** List of all entities that this model component depends on.
   *
   *  This lists all model components used in the definition of this component.
   */
  std::list<ModelComp*> dependencies;

  /** Components can be tagged by the tagDependencies method which sets
   *  this tag for this components and everything that it depends on 
   *  (i.e. everything listed in the dependencies list).                 */
  bool tag;            //!< true if part of the current 'needed' set

  /** For sets and parameters, this points to an object that gives the
   *  values and further specific information (Set for sets)
   */
  CompDescr *value;

 public:

  /** Constructor */
  ModelComp(const std::string& id_, compType type_,
            SyntaxNode *indexing_, SyntaxNode *attrib);

  /** Default constructor */
  ModelComp(const std::string& id);

  /** Destructor */
  virtual ~ModelComp();

  /** Set up an existing model comp to specified values */
  void setTo(const std::string& id_, compType type_,
             SyntaxNodeIx *indexing_, SyntaxNode *attrib);

  /** Set up list of dependencies for this component */
  void setUpDependencies();

  /** Detailed debugging output */
  void dump(std::ostream& fout) const;

  /** Print a detailed description of this model component and all its fields */
  void print() const;

  /** Print a one line description of the object (type and name) */
  void printBrief() const;

  /** Tag this components and all its dependencies recursively */
  void tagDependencies();

  /** Recalculate dependency list and re-resolve IDREF nodes */
  void reassignDependencies();

  /** Whether the component has a tag */
  bool isTagged() const { return tag; }

  /** Set the tag to false for all models: using global_list */
  static void untagAll();  

  /** Recursively set the tag to false for all models */
  static void untagAll(AmplModel *start);  

  /** Recursively write name of all tagged components */
  static void writeAllTagged(AmplModel *start);

  /** Write definition of all tagged components to file, using global_list */
  static void modifiedWriteAllTagged(std::ostream &fout);

  /** Retrieve the value of this component (if a set or parameter) */
  CompDescr* getValue() const { return value; }

  /** Set the value for this component */
  void setValue(CompDescr *val) { value = val; }

  /** Move this model component up in the model tree */
  void moveUp(int level);

  /** Duplicate the object: shallow copy */
  virtual ModelComp *clone() const;

  /** Duplicate the object: deep copy */
  ModelComp *deep_copy() const;

  // Virtual methods implemented only for stochastic models

  virtual void setStochModel(StochModel *stoch) { throw; }
  virtual void setStageSet(SyntaxNode *stageSet) { throw; }
  virtual void setDeterministic(bool det) { throw; }
  virtual const SyntaxNode* getStageSet() const { throw; }
  virtual void addStageName(const std::string& name) { throw; }
  virtual const std::vector<std::string>& getStageNames() const { throw; }
  virtual ModelComp* transcribeToModelComp(AmplModel *current_model,
                                           const std::string& nodedummy,
                                           const std::string& stagedummy,
                                           const int level) { throw; }

 protected:

  /** Instance number of ModelComp */
  int count;

 private:

  /** Global list of all defined ModelComps */
  static std::list<ModelComp*> global_list;

  /** Number of ModelComps defined */
  static int tt_count;

  /** Find dependencies for this component */
  void findDependencies(const SyntaxNode *nd);

};

std::string
getGlobalName(const ModelComp *node, const SyntaxNode *opn,
	      const AmplModel *current_model, int witharg);
std::string
getGlobalNameNew(const ModelComp *node, const SyntaxNode *opn,
                 const AmplModel *current_model, int witharg);

#endif
