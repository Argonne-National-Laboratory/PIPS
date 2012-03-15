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
#ifndef SETNODE_H
#define SETNODE_H

#include "nodes.h"
#include "CompDescr.h"
#include "sml.tab.h"
#include <cassert>
#include <string>
#include <vector>

/** @class SetNode
 *  This class represents a set in the the syntax tree.
 */
class SetNode: public SyntaxNode, public CompDescr {
protected:
   SetNode *within;

public:
   SetNode(int opCode_, SyntaxNode *node1=NULL, SyntaxNode *node2=NULL) :
      SyntaxNode(opCode_, node1, node2), within(NULL) {}
};

/** @class SimpleSet
 *  This class represents a contiguous set defined using 1..T or similar.
 */
class SimpleSet: public SetNode {
private:
   int lower_bound_;
   SyntaxNode *lbc_;
   int upper_bound_;
   SyntaxNode *ubc_;
   int interval_;
   bool parsed_; // did we succeed at parsing, or do we need to use ampl on it?

public:
   SimpleSet(SyntaxNode *bnd1, SyntaxNode *bnd2);

   /// Retrieve the members of the set
   std::vector<std::string> members(AmplModel& context) const;
};

/** @class ListSet
 *  This class represents a set with explicitly enumerated members.
 */
class ListSet: public SetNode {
private:
   std::vector<std::string> set;
public:
   ListSet(SyntaxNode *list) :
      SetNode(LBRACE, list) {}
};

/** @class CompositeSet
 * This class represents a set composed of other sets
 */
class CompositeSet: public SetNode {
public:
   CompositeSet(int opCode_, SyntaxNode *set1, SyntaxNode *set2) :
      SetNode(opCode_, set1, set2)
   {
      assert((opCode_ == CROSS) || (opCode_ == DIFF));
   }
};

#endif /* ifndef SETNODE_H */
