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

#include "Set.h"
#include "nodes.h"
#include "data.tab.h"
#include "misc.h"
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

using namespace std;

/* ---------------------------------------------------------------------------
Set::Set(SyntaxNode *list)
---------------------------------------------------------------------------- */
/*! Constructs the Set from a list of set elements given as a tree of SyntaxNodes
 *
 *  @param list a description of the set elements as read in from the data file
 *
 * This constructor assumes that the parameter list describes the set elements
 * in the following format:
 * - The top node is of type SyntaxNode.opCode=' ', 
 *    SyntaxNode::nval=[\#items in the list]
 * - Each child describes one element of the set and is of either of the 
 *   forms
 *   + SyntaxNode.opCode=ID|INT_VAL|FLOAT_VAL, SyntaxNode.values[0] = 
 *   + (...,...,..) represented as SyntaxNode.opCode=LBRACKET with one
 *     child of type COMMA. 
 *     The COMMA SyntaxNode is a list, with number of entries equal to the 
 *     dimension of the set and each element of type opCode=ID 
 *     (carrying the actual description of the set element).
 */
Set::Set(const ListNode &list):
  elements()
{
  SyntaxNode *item;

  assert(list.getOpCode() == ' ');
  
  // have a look at the first item to get the dimension of the set
  item = list[0];
  if (item->getOpCode() == ID || item->getOpCode() ==- 99)
    this->dim_ = 1;
  else {
    // otherwise this needs to be an element of form (.., .., ..)
    assert(item->getOpCode() == LBRACKET);
    item = item->front();
    assert(item->getOpCode() == COMMA);
    this->dim_ = item->nchild();
  }
  
  // then place all the elements on the set
  for(SyntaxNode::iterator i=list.begin(); i!=list.end(); ++i){
    item = *i;
    char **array = new char*[dim_];

    // and do some checking that all elements have the same dimension
    //    this->elements.push_back(item);
    if (dim_==1) {
      array[0] = (char *) item->getValue().c_str();
      add(SetElement(1,array));
    }else{
      assert(item->getOpCode() == LBRACKET);
      item = item->front();
      assert(item->getOpCode() == COMMA);
      if (dim_==item->nchild()){
        int j = 0;
        for(SyntaxNode::iterator k=item->begin(); k!=item->end(); ++k){
          SyntaxNode *idnd = *k;
          assert(idnd->getOpCode() == ID);
          array[j++] = (char*) idnd->front();
        }
        add(SetElement(dim_, array));
        //this->elements.push_back(array);
      }else{
        cerr << "First element in set has dim=" << dim_ << " later element '"
             << *i << "' has dim=" << item->nchild() << "\n";
        exit(1);
      }
    }
  }
}

/* ---------------------------------------------------------------------------
Set::toString
---------------------------------------------------------------------------- */
string 
Set::toString() const {

  map<SetElement, int, SetElement>::const_iterator iter;
  string str="";
  for( iter = elements.begin(); iter != elements.end(); ++iter ) {
    if (iter!=elements.begin()) str += " ";
    // iter is of type 'pair*'
    SetElement element = iter->first;
    int pos = iter->second;
    str += to_string(pos) + ":" + element.getVal();
  }

  return str;
}


/* ---------------------------------------------------------------------------
Set::add(SetElement newel)
---------------------------------------------------------------------------- */
void
Set::add(const SetElement& newel) {
  int n = elements.size();
  elements.insert(pair<SetElement, int>(newel, n));
}

/* ---------------------------------------------------------------------------
Set::findPos(string *el)
---------------------------------------------------------------------------- */
int
Set::findPos(const SetElement& el) const {

  map<SetElement, int, SetElement>::const_iterator iter = elements.find(el);
  if( iter != elements.end() ) {
    return iter->second;
  }else{
    cerr << "ERROR: Element '" << el.toString() << "' is not in Set "
         << this->toString() << endl;
    exit(1);
    return -1;
  }
}

