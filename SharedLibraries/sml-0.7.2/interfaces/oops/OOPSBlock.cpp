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

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include "OOPSBlock.h"
#include "GlobalVariables.h"

using namespace std;

/* ----------------------------------------------------------------------------
OOPSBlock::OOPSBlock(ExpandedModelInterface*, list<string>*)
---------------------------------------------------------------------------- */
OOPSBlock::OOPSBlock(ExpandedModelInterface *rowmod, ExpandedModelInterface *colmod)
{
  /* We need to:
      - take the list of variable names from colmod (colmod->listOfVarNames)
      - compare them against the variables defined by the NlFile attached
        to rowmod (rowmod->getName()+".col")
      - colmod->listOfVarNames will give the number of columns in this block
      - need a list of indices into the NlFile for these columns
   */
  
  if (GlobalVariables::prtLvl >= PRINT_INFO) {
    cout << "-------------------------OOPS Block---------------------------\n";
    cout << "Generate OOPSBlock: col: " << colmod->getName() << 
      "/ row: " << rowmod->getName() << endl;
  }

  this->emrow = rowmod;
  this->emcol = colmod;
  this->nvar = colmod->getNLocalVars();
  this->ncon = rowmod->getNLocalCons();

#ifdef REDUNDANT
  this->nlfile = rowmod->nlfile;

  this->lvar = (int*)malloc(nvar*sizeof(int));
  
  for(int i=0;i<nvar;i++) lvar[i] = -1;

  // ------- read the names of columns defined in this NlFile ------------
  // FIXME: Should these be remembered in the NlFile rather than needing to be
  //        read in?
  ifstream fin((rowmod->getName()+".col").c_str());

  if (!fin) {
    cout << "Cannot open column name file: "+rowmod->getName()+".col";
    exit(1);
  }
  
  list<string> colfilelist;
  string line;
  getline(fin, line);
  while(!fin.eof()){
    colfilelist.push_back(line);
    getline(fin, line);
  }
  
  if (GlobalVariables::prtLvl >= PRINT_VERBOSE)
    printf("Read %d lines from file %s.col\n",colfilelist.size(),
           rowmod->getName().c_str());

  // -------------- compare this listOfVarNames against this list
  int i=0;
  for(list<string>::iterator p=colmod->listOfVarNames.begin();
      p!=colmod->listOfVarNames.end();p++){
    // (*p) is a name of a local variable. Should see if we can find this
    // in this NlFile
    int cnt=0;
    for(list<string>::iterator q=colfilelist.begin();q!=colfilelist.end();q++){
      if ((*p)==(*q)){
        lvar[i] = cnt;
        break;
      }
      cnt++;
    }
    i++;
  }

  //printf("Row model: %s\n",rowmod->getName().c_str());
  //printf("Col model: %s\n",colmod->getName().c_str());
  //printf("Nb_row = %d\n",ncon);

  if (GlobalVariables::prtLvl >= PRINT_VERBOSE) {
    cout << "The NlFile declares these variables:\n"; 
    int cnt=0;
    for(list<string>::iterator p=colfilelist.begin();p!=colfilelist.end();p++){
      cout << setw(2) << cnt << ": " << *p << endl;
      cnt++;
    }

    cout << "The column block defines " << nvar << " variables:\n";
    cnt =0;
    for(list<string>::iterator p=colmod->listOfVarNames.begin();
        p!=colmod->listOfVarNames.end();p++){
      cout << " at " << setw(2) << lvar[cnt] << ":  " << *p << endl;
      cnt++;
    }
  }
 
  if (GlobalVariables::prtLvl >= PRINT_INFO)
    cout << "OOPS Block: " << rowmod->getName() << "(rw)/" <<
       colmod->getName() << "(cl): " << ncon << "x" << nvar << "\n";

#endif  

}
