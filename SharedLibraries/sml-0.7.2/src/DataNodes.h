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
#ifndef DATANODES_H
#define DATANODES_H

#include "nodes.h"
#include "data.tab.h"

/** Wraps a set template with an object list */
class SetSpec : public SyntaxNode {
  public:
   const SyntaxNode *tmpl;
   const ListNode *list;

  public:
   SetSpec(const SyntaxNode *new_tmpl, const ListNode *new_list) :
      SyntaxNode(TOKSETSPEC), tmpl(new_tmpl), list(new_list) {}
};

#endif
