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

#include "SetElement.h"
#include "nodes.h"
#include <string>
#include <cassert>

using namespace std;

SetElement::SetElement(int n_, char **new_val) :
   n(n_), val(NULL)
{
  if(n==0) return;
  assert(new_val);

  val = new string[n];
  for(int i=0; i<n; i++)
     val[i] = new_val[i];
}

SetElement::SetElement(int n_, IDNode **new_val) :
   n(n_), val(NULL)
{
  if(n==0) return;
  assert(new_val);

  val = new string[n];
  for(int i=0; i<n; i++)
    val[i] = new_val[i]->id();
}

SetElement::~SetElement() {
  // FIXME
  //if(val) delete [] val;
}

bool SetElement::operator()(const SetElement& el1, const SetElement& el2) const
{
  assert(el1.n==el2.n);
  for(int i=0;i<el1.n;i++){
    int cmp = el1.val[i].compare(el2.val[i]);
    return (cmp < 0) ? true : false;
  }
  return false;
}

string SetElement::toString() const {

  string str(val[0]);
  for(int i=1;i<n;i++){
    str += ",";
    str += val[i];
  }

  return str;
}

