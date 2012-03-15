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
#include <cstring> // for strcmp()
#include <cstdlib> // for exit()
#include "sml.h"
#include "sml-mps.h"

using namespace std;

const string progname = "smlmps";

void writeHelp(ostream& out, const string& programname) {
   out << "Syntax:" << endl;
   out << "   " << programname
       << " [OPTIONS] modelfile datafile mpsfile\n\n";
   out << "Option summary:" << endl;
   out << " -d            Enables debug information when reading model "
      "file." << endl;
   out << " --help        Displays this help information." << endl;
   out << " modelfile     File containing SML model." << endl;
   out << " datafile      File containing SML data." << endl;
   out << " mpsfile       Filename of MPS file to write." << endl;
}

int analyseOptions(int argc, char **argv, string &modelfilename,
                   string &datafilename, string &mpsfilename, bool &debug) {
   int found = 0;

   for (int i=1;i<argc;i++){
      if(strcmp(argv[i], "-d")==0){
         debug = true;
      }
      else if(strcmp(argv[i], "--help")==0) {
         writeHelp(cout, progname);
         exit(0);
      }
      else {
         switch(found) {
         case 0:
            // first proper argument is the model file to read
            modelfilename = argv[i];
            break;
         case 1:
            // next one is data file
            datafilename = argv[i];
            break;
         default:
            mpsfilename = argv[i];
            break;
         }
         found++;
      }
   }

   if(modelfilename=="" || datafilename=="" || mpsfilename=="") {
      cerr << "ERROR: all of modelfile, datafile and mpsfile " 
         "must be supplied." << endl << endl;
      return 1;
   }

   return 0;
}

/* ----------------------------------------------------------------------------
main
---------------------------------------------------------------------------- */
int main(int argc, char **argv) {
   string modelfilename = "";
   string datafilename = "";
   string mpsfilename = "";
   bool debug = false;

   cout << "SML MPS generator, " << sml_version() << endl;
   cout << "(c) Andreas Grothey and Jonathan Hogg, "
      "University of Edinburgh 2009" << endl << endl;

   int rv = analyseOptions(argc, argv,
                           modelfilename, datafilename, mpsfilename, debug);
   if (rv)
     return rv;

   if(debug) 
      cout << "======================================================" << endl;

   ExpandedModelInterface *em = sml_generate(modelfilename, datafilename, debug);
   if (!em)
     return 1;

   SML_MPS_driver(em, mpsfilename);

   // clean up
   delete em;

   return 0;
}
