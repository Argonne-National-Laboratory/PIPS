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
#ifndef COMPDESCR
#define COMPDESCR

#include <string>

/** @class CompDescr
 *  A superclass for all Component Description classes such as Set and
 *  CompDescrParam.
 *
 *  Its main purpose is to provide a type (other than void*) that can be
 *  used in the ->value field of ModelComp.
 */
class CompDescr{
 public:
  virtual ~CompDescr() {}
  virtual std::string toString() const {return "";}
};

#endif
