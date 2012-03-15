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
#ifndef SETELH
#define SETELH

#include <string>

class IDNode;

/** @class SetElement
 *  This class describes an element of a set. It is basically an array of
 *  strings (char*) together with a size.
 */

class SetElement {

 private:
  const int n;
  std::string *val;

 public:
  /* ----------------------------- methods -------------------------------*/

  SetElement(int n=0, char **val=NULL);
  SetElement(int n, IDNode **val);
  ~SetElement();

  bool operator()(const SetElement& el1, const SetElement& el2) const;

  std::string getVal() const {
    return val[0];
  }

  std::string toString() const;

};

#endif
