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
#ifndef COMPDESCRPARAM_H
#define COMPDESCRPARAM_H

#include "CompDescr.h"
#include "Set.h"
#include <string>

class ModelComp;
class SyntaxNode;
class SyntaxNodeIx;

/** @class CompDescrParam
 *  This class describes a parameter: it consists of
 *   - a list of indexing sets (these might be multidimensional themselves) 
 *   - a list of parameter values (in dense format?) - 
 *         i.e. a multidimensional array
 */  
class CompDescrParam: public CompDescr{

 private:

  /** Number of indices */
  int nix;

  /** Number of indexing sets (different to nix) */
  int nsets;

  /** Pointers to the indexing sets */
  Set **indices;

  /** Total number of entries: product of the number of elements in all 
   *  indexing sets */
  int n;

  /** Number of values that are given so far */
  int nread;

  /** The array of values */
  double *values;

 public:

  // ------------------ constructors ---------------------
  
  /** Construct given data file description and model component */
  CompDescrParam(ModelComp *mc, SyntaxNode *desc); //!< Construct from data

  std::string toString() const;

 private:

  /** Service routine that processes a tree below a TOKVALUETABLELIST node */
  void processValueTableList(const SyntaxNode *node, const SyntaxNodeIx *ix);

};

#endif
